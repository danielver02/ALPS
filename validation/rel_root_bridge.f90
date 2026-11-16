! External test hooks keep production module dependencies and algorithms intact.
subroutine root_begin(x)
  use rel_root_diagnostics, only: root_solver_start
  implicit none
  double complex,intent(in)::x
  call root_solver_start(x)
end subroutine

subroutine root_step(x,f,dx,j)
  use rel_root_diagnostics, only: root_solver_step
  implicit none
  double complex,intent(in)::x,f,dx
  integer,intent(in)::j
  call root_solver_step(x,f,dx,j)
end subroutine

subroutine root_scan_result(x,f)
  use rel_root_diagnostics, only: root_result,root_context
  implicit none
  double complex,intent(in)::x,f
  call root_context(3)
  call root_check_queries(x)
  call root_result(x,f)
end subroutine

subroutine root_before_disp(om)
  use rel_root_diagnostics, only: root_disp_enter
  implicit none
  double complex,intent(in)::om
  call root_disp_enter
  call root_check_queries(om)
end subroutine

subroutine root_check_queries(om)
  use alps_var
  use alps_analyt, only: eval_fit,CHEB_SUCCESS
  use alps_fns, only: determine_resonances
  use alps_io, only: alps_continuation_error
  use rel_root_diagnostics, only: root_query_record,root_query_failure
  implicit none
  double complex,intent(in)::om
  double complex::pole,f,df
  double precision::step
  integer::s,nn,sgn,g,saved_sproc
  logical::rp,rm
  if(.not.proc0.or.aimag(om)>0.d0) return
  saved_sproc=sproc
  do s=1,nspec
     if(usebM(s).or..not.relativistic(s)) cycle
     sproc=s
     step=pparbar_rel(s,2,2)-pparbar_rel(s,2,1)
     do nn=0,nmax(s)
        call determine_resonances(om,nn,rp,rm)
        do sgn=1,-1,-2
           if(nn==0.and.sgn==-1) cycle
           if(sgn==1.and..not.rp) cycle
           if(sgn==-1.and..not.rm) cycle
           do g=1,ngamma-1
              pole=gamma_rel(s,g,1)*om*vA/kpar-dble(sgn*nn)*qs(s)*vA/(kpar*ms(s))
              if(real(pole)**2>gamma_rel(s,g,1)**2-1.d0) cycle
              call query(g+1,pole,.false.)
              if(g==1) then
                 call query(g,pole,.false.)
              else
                 call query(g-1,pole,.false.)
              endif
              if(ACmethod(s)==2) then
                 call query(g,pole,.true.)
              else
                 call query(g,pole+step,.false.)
                 call query(g,pole-step,.false.)
              endif
           enddo
        enddo
     enddo
  enddo
  sproc=saved_sproc
contains
  subroutine query(row,u,wants_derivative)
    integer,intent(in)::row
    double complex,intent(in)::u
    logical,intent(in)::wants_derivative
    integer::code
    if(wants_derivative) then
       f=eval_fit(s,row,u,dfpar=df,status=code)
    else
       f=eval_fit(s,row,u,status=code)
    endif
    call root_query_record(om,s,row,u,wants_derivative,code)
    if(code/=CHEB_SUCCESS) then
       call root_query_failure(s,row,u,code)
       call alps_continuation_error(s,row,code,u) ! Existing explicit fatal policy.
    endif
  end subroutine
end subroutine

double complex function root_serial_disp(om) result(d)
  ! The established Task-3D singleton assembly adapter, using production
  ! full_integrate/int_ee_rel and the same determinant algebra as disp.
  use alps_var
  use alps_fns, only: full_integrate,determine_bessel_array,determine_resonances
  use alps_fns_rel, only: int_ee_rel
  implicit none
  double complex,intent(in)::om
  integer,parameter::ti(6)=[1,2,3,1,1,2],tj(6)=[1,2,3,2,3,3]
  double complex::raw(2,3,3),terms(6),eps(3,3)
  integer::s,nn,mode
  logical::rp,rm
  if(any(usebM).or..not.all(relativistic).or.nspec/=2) error stop 'serial root fixture requires two numerical relativistic species'
  raw=(0.d0,0.d0)
  do s=1,nspec
     sproc=s;nlim=[0,nmax(s)];call determine_bessel_array
     terms=(0.d0,0.d0)
     do nn=0,nmax(s)
        call determine_resonances(om,nn,rp,rm)
        do mode=1,6
           if(nn==0.and.(mode==1.or.mode==4.or.mode==5)) cycle
           terms(mode)=terms(mode)+full_integrate(om,nn,mode,rp)
           if(nn/=0) terms(mode)=terms(mode)+full_integrate(om,-nn,mode,rm)
        enddo
     enddo
     terms(3)=terms(3)+int_ee_rel(om)
     do mode=1,6
        raw(s,ti(mode),tj(mode))=terms(mode)*ns(s)*qs(s)
     enddo
  enddo
  chi0=raw/(om*om*vA*vA)
  chi0(:,2,1)=-chi0(:,1,2);chi0(:,3,1)=chi0(:,1,3);chi0(:,3,2)=-chi0(:,2,3)
  eps=sum(raw,dim=1)
  eps(1,1)=eps(1,1)+(om*vA)**2;eps(2,2)=eps(2,2)+(om*vA)**2;eps(3,3)=eps(3,3)+(om*vA)**2
  wave=eps;wave(1,1)=eps(1,1)-kpar*kpar;wave(2,2)=eps(2,2)-kpar*kpar-kperp*kperp
  wave(3,3)=eps(3,3)-kperp*kperp;wave(1,3)=eps(1,3)+kpar*kperp
  wave(2,1)=-wave(1,2);wave(3,1)=wave(1,3);wave(3,2)=-wave(2,3)
  d=wave(1,1)*(wave(2,2)*wave(3,3)+wave(2,3)**2)+2.d0*wave(1,2)*wave(2,3)*wave(1,3)- &
       wave(1,3)**2*wave(2,2)+wave(1,2)**2*wave(3,3)
  sproc=0
end function
