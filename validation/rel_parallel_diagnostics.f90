! Test-only observation of the actual Landau routine, in scratch builds only.
module rel_parallel_diagnostics
  use alps_var, only : iproc,sproc,ACmethod,poly_order,gamma_rel,pparbar_rel,ms,vA,kpar,kperp
  use alps_analyt, only : eval_fit,CHEB_SUCCESS
  use alps_distribution_analyt, only : distribution_analyt
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer :: trace_unit,query_unit
  double precision :: quadratic
  public :: parallel_diagnostic_open,record_parallel
contains
  subroutine parallel_diagnostic_open(b)
    double precision, intent(in) :: b
    character(80) :: name
    quadratic=b
    write(name,'(a,i0,a)') 'gamma_rank',iproc,'.csv'
    open(newunit=trace_unit,file=trim(name),status='new')
    write(trace_unit,'(a)') 'degree,species,harmonic,mode,row,omega_re,omega_im,kpar,kperp,u_re,u_im,gamma_re_bits,gamma_im_bits'
    write(name,'(a,i0,a)') 'parallel_rank',iproc,'.csv'
    open(newunit=query_unit,file=trim(name),status='new')
    write(query_unit,'(a)') 'degree,species,harmonic,mode,row,omega_re,omega_im,kpar,kperp,u_re,u_im,'// &
         'actual_re,actual_im,cheb_re,cheb_im,fd_re,fd_im,shape_re,shape_im,value_status,derivative_status'
  end subroutine

  subroutine record_parallel(om,harmonic,mode,row,u,fg,fu)
    double complex, intent(in) :: om,u,fg,fu
    integer, intent(in) :: harmonic,mode,row
    double complex :: f,du,plus,minus,fd,shape,shape_du
    double precision :: step
    integer :: status,vs,ds
    integer(int64) :: bits(2)
    if(ACmethod(sproc)/=2) return
    bits=[transfer(real(fg,kind(1.d0)),0_int64),transfer(aimag(fg),0_int64)]
    write(trace_unit,'(5(i0,","),6(es26.17e3,","),z16.16,",",z16.16)') &
         poly_order(sproc),sproc,harmonic,mode,row,real(om),aimag(om),kpar,kperp,real(u),aimag(u),bits
    ! Diagnostic-only extra calls do not occur in the production hot path.
    f=eval_fit(sproc,row,u,dfpar=du,status=status,value_status=vs,derivative_status=ds)
    if(status/=CHEB_SUCCESS) error stop 'accepted susceptibility query has failed derivative'
    step=pparbar_rel(sproc,2,2)-pparbar_rel(sproc,2,1)
    plus=eval_fit(sproc,row,u+step,status=status)
    if(status/=CHEB_SUCCESS) error stop 'diagnostic parallel plus query failed'
    minus=eval_fit(sproc,row,u-step,status=status)
    if(status/=CHEB_SUCCESS) error stop 'diagnostic parallel minus query failed'
    fd=(plus-minus)/(2.d0*step)
    shape=distribution_analyt(sproc,gamma_rel(sproc,row,1),u,1,ms(sproc),vA)
    shape_du=-2.d0*quadratic*u*shape
    write(query_unit,'(*(g0,:,","))') poly_order(sproc),sproc,harmonic,mode,row, &
         real(om),aimag(om),kpar,kperp,real(u),aimag(u),real(fu),aimag(fu),real(du),aimag(du), &
         real(fd),aimag(fd),real(shape_du),aimag(shape_du),vs,ds
  end subroutine
end module
