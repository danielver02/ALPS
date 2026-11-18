! Test-only bookkeeping. No solver, fitted-model copies or production changes.
module rel_root_diagnostics
  use alps_var
  use alps_analyt, only: rel_cheb_collection_info
  use alps_com, only: rel_cheb_distribution_info
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  integer, parameter :: dp=kind(1.d0)
  integer :: root_unit,iteration_unit,tensor_unit,query_unit,event_unit,reference_unit
  integer :: solver_id=0,last_iteration=0,initial_fits,initial_builds,initial_transfers,initial_collectives
  integer :: active_method=0,stage=0,disp_calls=0,failed_queries=0
  logical :: last_converged=.false.
  complex(dp) :: initial_frequency
  real(dp) :: shape_c(2),shape_b
  public :: root_diagnostic_open,root_context,root_solver_start,root_solver_step,root_result
  public :: root_query_record,root_disp_enter,root_tensor_record,root_query_failure
  public :: shape_c,shape_b
  public :: root_reference_close,root_accepted
contains
  logical function root_accepted(x,f)
    complex(dp),intent(in)::x,f
    root_accepted=.false.
    if(proc0) root_accepted=last_converged.and.aimag(x)<0.d0.and. &
         abs(f)/max(maxval(abs(wave))**3,tiny(1.d0))<=D_tol
  end function
  subroutine root_reference_close
    if(proc0.and.active_method==0) close(reference_unit)
  end subroutine
  subroutine root_diagnostic_open(c,b,method)
    real(dp), intent(in) :: c(2),b
    integer, intent(in) :: method
    integer :: nsr,last,calls,builds,transfers,collectives
    integer(int64) :: bytes,volume,evaluations
    shape_c=c;shape_b=b;active_method=method
    call rel_cheb_collection_info(nsr,last,initial_fits,bytes,initial_builds,evaluations)
    call rel_cheb_distribution_info(initial_transfers,initial_collectives,volume)
    if(.not.proc0) return
    open(newunit=root_unit,file='roots.csv',status='new')
    open(newunit=iteration_unit,file='iterations.csv',status='new')
    open(newunit=tensor_unit,file='tensor.csv',status='new')
    open(newunit=query_unit,file='queries.csv',status='new')
    open(newunit=event_unit,file='events.csv',status='new')
    write(root_unit,'(a)') 'method,solver,stage,kpar,kperp,initial_re,initial_im,omega_re,omega_im,iterations,solver_converged,det_re,det_im,abs_det,scaled_det,residual_ok,ez_fraction,ey_fraction,helicity,evaluations,fit_calls,build_calls,transfers,collectives'
    write(iteration_unit,'(a)') 'method,solver,stage,kpar,kperp,iteration,omega_re,omega_im,det_re,det_im,step_re,step_im,solver_converged'
    write(tensor_unit,'(a)') 'method,stage,kpar,kperp,omega_re,omega_im,species,i,j,re,im,det_re,det_im,evaluation_status'
    write(query_unit,'(a)') 'method,stage,kpar,kperp,omega_re,omega_im,species,row,u_re,u_im,wants_derivative,status,nearest_u,table_value,shape_value,relative_table_error,inside_source'
    write(event_unit,'(a)') 'method,stage,event,code,species,row,u_re,u_im'
    if(method==0) open(newunit=reference_unit,file='reference.dat',status='new')
  end subroutine

  subroutine root_context(new_stage)
    integer,intent(in)::new_stage
    stage=new_stage
  end subroutine

  subroutine root_disp_enter
    integer :: nsr,last,calls,builds,transfers,collectives
    integer(int64) :: bytes,volume,evaluations
    call rel_cheb_collection_info(nsr,last,calls,bytes,builds,evaluations)
    call rel_cheb_distribution_info(transfers,collectives,volume)
    if(calls/=initial_fits.or.builds/=initial_builds.or.transfers/=initial_transfers.or. &
         collectives/=initial_collectives) error stop 'root iteration fitted or redistributed models'
    disp_calls=disp_calls+1
  end subroutine

  subroutine root_solver_start(x)
    complex(dp),intent(in)::x
    solver_id=solver_id+1;initial_frequency=x;last_iteration=0;last_converged=.false.
    stage=2 ! recovery/scan; preflight uses stage 1.
  end subroutine

  subroutine root_solver_step(x,f,dx,j)
    complex(dp),intent(in)::x,f,dx
    integer,intent(in)::j
    last_iteration=j;last_converged=abs(dx)<D_tol.or.abs(f)==0.d0
    if(proc0) then
       write(iteration_unit,'(*(g0,:,","))') active_method,solver_id,stage,kpar,kperp,j, &
            real(x),aimag(x),real(f),aimag(f),real(dx),aimag(dx),merge(1,0,last_converged)
       flush(iteration_unit)
    endif
  end subroutine

  subroutine root_result(x,f)
    complex(dp),intent(in)::x,f
    complex(dp) :: cross(3),best(3)
    real(dp) :: norm,bestnorm,ez,ey,helicity,eta
    integer :: i,j,nsr,last,calls,builds,transfers,collectives
    integer(int64) :: bytes,volume,evaluations
    if(.not.proc0) return
    ! A null-vector polarization diagnostic, with the best-conditioned row pair.
    bestnorm=0.d0;best=(0.d0,0.d0)
    do i=1,2
       do j=i+1,3
          cross=[wave(i,2)*wave(j,3)-wave(i,3)*wave(j,2), &
               wave(i,3)*wave(j,1)-wave(i,1)*wave(j,3),wave(i,1)*wave(j,2)-wave(i,2)*wave(j,1)]
          norm=sum(abs(cross)**2)
          if(norm>bestnorm) then
             best=cross;bestnorm=norm
          endif
       enddo
    enddo
    ez=0.d0;ey=0.d0;helicity=0.d0
    if(bestnorm>0.d0) then
       ez=abs(best(3))**2/bestnorm;ey=abs(best(2))**2/bestnorm
    endif
    norm=abs(best(1))**2+abs(best(2))**2
    if(norm>0.d0) helicity=2.d0*aimag(conjg(best(1))*best(2))/norm
    eta=abs(f)/max(maxval(abs(wave))**3,tiny(1.d0))
    call rel_cheb_collection_info(nsr,last,calls,bytes,builds,evaluations)
    call rel_cheb_distribution_info(transfers,collectives,volume)
    write(root_unit,'(*(g0,:,","))') active_method,solver_id,stage,kpar,kperp, &
         real(initial_frequency),aimag(initial_frequency),real(x),aimag(x),last_iteration, &
         merge(1,0,last_converged),real(f),aimag(f),abs(f),eta,merge(1,0,eta<=D_tol), &
         ez,ey,helicity,disp_calls,calls,builds,transfers,collectives
    flush(root_unit)
    if(active_method==0) then
       write(reference_unit,'(4es26.17e3,2i4)') kpar,kperp,real(x),aimag(x),merge(1,0,last_converged),merge(1,0,eta<=D_tol)
       flush(reference_unit)
    endif
    call root_tensor_record(x,f)
  end subroutine

  subroutine root_tensor_record(x,f)
    complex(dp),intent(in)::x,f
    integer::s,i,j
    if(.not.proc0) return
    if(.not.ieee_is_finite(real(f)).or..not.ieee_is_finite(aimag(f)).or. &
         .not.all(ieee_is_finite(real(chi0))).or..not.all(ieee_is_finite(aimag(chi0)))) &
         error stop 'nonfinite fixed-frequency tensor/determinant'
    do s=1,nspec
       do i=1,3
          do j=i,3
             write(tensor_unit,'(*(g0,:,","))') active_method,stage,kpar,kperp,real(x),aimag(x), &
                  s,i,j,real(chi0(s,i,j)),aimag(chi0(s,i,j)),real(f),aimag(f),0
          enddo
       enddo
    enddo
    flush(tensor_unit)
  end subroutine

  subroutine root_query_record(om,s,row,u,want_derivative,code)
    complex(dp),intent(in)::om,u
    integer,intent(in)::s,row,code
    logical,intent(in)::want_derivative
    real(dp)::nearest,analytic,relerr,g,rmax,sourceperp,dist,bestdist
    integer::j,index
    if(.not.proc0) return
    if(stage/=1.and.stage/=3.and.code==0) return
    g=gamma_rel(s,row,0);bestdist=huge(1.d0);index=-1
    do j=0,npparbar
       if(f0_rel(s,row,j)==-1.d0) cycle
       dist=abs(pparbar_rel(s,row,j)-real(u))
       if(dist<bestdist) then
          index=j;bestdist=dist
       endif
    enddo
    if(index<0) then
       write(event_unit,'(*(g0,:,","))') active_method,stage,'no_geometric_node',code,s,row,real(u),aimag(u)
       flush(event_unit);return
    endif
    nearest=pparbar_rel(s,row,index);analytic=shape_c(s)*exp(-2.d0*g-shape_b*nearest**2)
    relerr=abs(f0_rel(s,row,index)-analytic)/max(analytic,1.d-12*maxval(f0_rel(s,:,:)))
    sourceperp=sqrt(max(0.d0,g*g-1.d0-nearest*nearest))*ms(s)/vA
    rmax=maxval(pp(s,:,:,1))
    write(query_unit,'(*(g0,:,","))') active_method,stage,kpar,kperp,real(om),aimag(om),s,row, &
         real(u),aimag(u),merge(1,0,want_derivative),code,nearest,f0_rel(s,row,index),analytic,relerr, &
         merge(1,0,sourceperp<=rmax*(1.d0+1.d-14))
  end subroutine

  subroutine root_query_failure(s,row,u,code)
    integer,intent(in)::s,row,code
    complex(dp),intent(in)::u
    failed_queries=failed_queries+1
    if(proc0) then
       write(event_unit,'(*(g0,:,","))') active_method,stage,'continuation_prerequisite',code,s,row,real(u),aimag(u)
       flush(event_unit);flush(query_unit)
    endif
  end subroutine
end module
