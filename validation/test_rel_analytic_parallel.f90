! Manufactured-row integration test. No production fitting/transform changes.
program test_rel_analytic_parallel
  use mpi
  use alps_var
  use alps_analyt
  use alps_com, only : pass_rel_cheb_models,rel_cheb_distribution_info
  use alps_fns_rel, only : landau_integrate_rel,CBESSJ
  use alps_distribution_analyt, only : distribution_analyt,set_susceptibility_reference
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: scenario,grid,is,g,j,iw,nn,imode,ik,mode,code,file,passed=0,failed=0,total_failed
  integer :: calls,builds,nsp,ng,old_calls,old_builds,transfers,collectives,old_transfers,old_collectives
  integer(int64) :: bytes,volume,evaluations
  real(dp) :: amplitudes(2),b,u,gamma,du,dg
  complex(dp) :: omega,got,want,fd,old_fd,f,d,coarse(2,3,-1:1,2,2),analytic_gamma
  logical :: legacy_only
  character(32) :: option

  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0; call get_command_argument(1,option); legacy_only=option=='legacy'
  nspec=2; nspec_rel=2; ngamma=6; pi=4.d0*atan(1.d0); vA=.6d0; kpar=.4d0
  writeOut=.false.; profile_startup=.false.
  allocate(ms(2),qs(2),relativistic(2),usebM(2),ACmethod(2),poly_order(2))
  ms=[1.2d0,2.d0]; qs=[1.d0,-1.d0]; relativistic=.true.; usebM=.false.; ACmethod=2
  allocate(n_fits(2),fit_type(2,1),perp_correction(2,1),param_fit(2,0:ngamma,3,1))
  n_fits=1; perp_correction=2.d0
  amplitudes=[.7d0,.23d0]
  if(proc0) then
     open(newunit=unit_error,file='continuation-errors.log',status='new')
     open(newunit=file,file='legacy.dat',status='new')
  endif

  do scenario=1,2
     b=0.d0; poly_order=0; fit_type=4
     if(scenario==2) then
        b=.18d0; poly_order=2; fit_type=5
     endif
     call set_susceptibility_reference(amplitudes,b)
     do grid=1,3
        npparbar=6*2**(grid-1)
        if(allocated(gamma_rel)) deallocate(gamma_rel,pparbar_rel,f0_rel)
        allocate(gamma_rel(2,0:ngamma,0:npparbar),pparbar_rel(2,0:ngamma,0:npparbar),f0_rel(2,0:ngamma,0:npparbar))
        do is=1,2
           do g=0,ngamma
              gamma=1.d0+2.d0*g/ngamma
              param_fit(is,g,1,1)=amplitudes(is); param_fit(is,g,2,1)=b; param_fit(is,g,3,1)=0.d0
              do j=0,npparbar
                 u=-1.2d0+2.4d0*j/npparbar
                 gamma_rel(is,g,j)=gamma; pparbar_rel(is,g,j)=u
                 f0_rel(is,g,j)=real(shape(is,g,cmplx(u,0.d0,dp)),dp)
                 if(u*u>gamma*gamma-1.d0) f0_rel(is,g,j)=-1.d0
              enddo
           enddo
        enddo
        if(option=='fatal_derivative') then
           ! A finite value but unrepresentable physical-coordinate derivative.
           ! This must take the established status/communicator failure path.
           poly_order=1
           do is=1,2
              do g=0,ngamma
                 do j=0,npparbar
                    pparbar_rel(is,g,j)=1.d-310*(-1.d0+2.d0*j/npparbar)
                    f0_rel(is,g,j)=10.d0**(-1.d0+2.d0*j/npparbar)
                 enddo
              enddo
           enddo
        endif
        ACmethod=2; call initialize_rel_cheb(code); call check(code==CHEB_SUCCESS,'root initialization')
        call pass_rel_cheb_models(code); call check(code==CHEB_SUCCESS,'coefficient distribution')
        call rel_cheb_collection_info(nsp,ng,old_calls,bytes,old_builds,evaluations)
        call rel_cheb_distribution_info(old_transfers,old_collectives,volume)
        if(option=='fatal_derivative') then
           f=eval_fit(1,1,(0.d0,0.d0),status=code)
           call check(code==CHEB_SUCCESS.and.ieee_is_finite(real(f,dp)),'value alone is representable')
           sproc=1; kperp=.2d0
           got=landau_integrate_rel((0.d0,0.d0),0,3)
           error stop 'unrepresentable analytic derivative silently accepted'
        endif
        do is=1,2
           sproc=is
           do iw=1,3
              omega=(.18d0,-.025d0)
              if(iw==2) omega=(.31d0,-.08d0)
              if(iw==3) omega=(.18d0,0.d0)
              do nn=-1,1
                 do ik=1,2
                    kperp=0.d0
                    if(ik==2) kperp=.2d0
                    do imode=1,2
                       mode=1
                       if(imode==2) mode=3
                       want=reference(omega,nn,mode,.true.,.false.)
                       fd=reference(omega,nn,mode,.false.,.false.)
                       ACmethod=0; got=landau_integrate_rel(omega,nn,mode)
                       call close(got,fd,2.d-12,'analytic-function mode retains parallel FD')
                       if(proc0) write(file,'(5i4,2es26.17e3)') scenario,npparbar,is,iw,nn,real(got,dp),aimag(got)
                       ACmethod=1; got=landau_integrate_rel(omega,nn,mode)
                       call close(got,fd,2.d-12,'functional mode retains parallel FD')
                       if(proc0) write(file,'(5i4,2es26.17e3)') scenario,npparbar,is,iw,nn,real(got,dp),aimag(got)
                       if(legacy_only) cycle
                       ACmethod=2; got=landau_integrate_rel(omega,nn,mode)
                       call close(got,want,3.d-11,'Chebyshev analytic parallel derivative; Gamma FD retained')
                       if(grid==1) then
                          coarse(is,iw,nn,ik,imode)=got
                       else
                          call close(got,coarse(is,iw,nn,ik,imode),3.d-11,'analytic result independent of parallel step')
                       endif
                       if(scenario==2.and.grid==1.and.iw==1.and.nn==0.and.mode==3) then
                          call check(abs(want-fd)>1.d-6*max(abs(want),1.d-10),'parallel FD differs detectably from analytic derivative')
                          analytic_gamma=reference(omega,nn,mode,.true.,.true.)
                          call check(abs(want-analytic_gamma)>1.d-4*abs(want),'Gamma FD distinguishable from analytic Gamma derivative')
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
        call rel_cheb_collection_info(nsp,ng,calls,bytes,builds,evaluations)
        call rel_cheb_distribution_info(transfers,collectives,volume)
        call check(calls==old_calls.and.builds==old_builds,'integration does not refit')
        call check(transfers==old_transfers.and.collectives==old_collectives,'integration does not redistribute coefficients')
     enddo
  enddo
  if(.not.legacy_only) call adapter_checks
  if(proc0) then
     close(file); close(unit_error)
  endif
  call mpi_allreduce(failed,total_failed,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
  write(*,'(a,i0,a,i0,a,i0)') 'Analytic parallel checks rank=',iproc,' passed=',passed,' failed=',failed
  if(total_failed/=0) call mpi_abort(MPI_COMM_WORLD,98,ierror)
  if(proc0) print '(a)', 'PASS: analytic parallel derivative with unchanged Gamma stencil and integration'
  call mpi_finalize(ierror)
contains
  subroutine adapter_checks
    ! Independent analytic checks through the SAME adapter used by Landau,
    ! including drift and independently reported reconstruction conditions.
    real(dp) :: cg(2,0:1,0:1),cu(2,0:1,0:1),cf(2,0:1,0:1)
    real(dp) :: gg,x,gd,beta,theta,started,finished
    complex(dp) :: z,want,dwant,v,d
    integer :: kind,is,row,j,k,code,vs,ds
    gd=1.d0/sqrt(1.d0-.3d0**2); beta=.3d0; theta=.7d0
    ACmethod=2
    do kind=1,3
       poly_order=0
       if(kind==2) poly_order=2
       if(kind==3) poly_order=1
       do is=1,2
          do row=0,ngamma
             gg=gamma_rel(is,row,0)
             do j=0,npparbar
                x=pparbar_rel(is,row,j)
                f0_rel(is,row,j)=amplitudes(is)*exp(-2.d0*gg)
                if(kind==2) f0_rel(is,row,j)=amplitudes(is)*exp(-2.d0*gg-.18d0*x*x)
                if(kind==3) f0_rel(is,row,j)=amplitudes(is)*exp(-gd*(gg-beta*x)/theta)
                if(x*x>gg*gg-1.d0) f0_rel(is,row,j)=-1.d0
             enddo
          enddo
       enddo
       call initialize_rel_cheb(code); call check(code==CHEB_SUCCESS,'analytic adapter row initialization')
       call pass_rel_cheb_models(code); call check(code==CHEB_SUCCESS,'analytic adapter row distribution')
       do is=1,2
          do row=1,ngamma
             gg=gamma_rel(is,row,0)
             do k=1,3
                z=(.2d0,0.d0)
                if(k==2) z=cmplx(.95d0*sqrt(gg*gg-1.d0),.02d0,dp)
                if(k==3) z=(1.3d0,.08d0)
                want=amplitudes(is)*exp(-2.d0*gg); dwant=(0.d0,0.d0)
                if(kind==2) then
                   want=want*exp(-.18d0*z*z); dwant=-.36d0*z*want
                elseif(kind==3) then
                   want=amplitudes(is)*exp(-gd*(gg-beta*z)/theta); dwant=gd*beta*want/theta
                endif
                v=eval_fit(is,row,z,dfpar=d,status=code,value_status=vs,derivative_status=ds)
                call check(code==CHEB_SUCCESS.and.vs==CHEB_SUCCESS.and.ds==CHEB_SUCCESS,'analytic adapter statuses')
                call close(v,want,1.d-11,'analytic adapter physical value')
                call check(abs(d-dwant)<=1.d-11*max(abs(want),abs(dwant)),'isotropic/anisotropic/drifting physical-u derivative')
             enddo
          enddo
       enddo
    enddo
    cg=1.5d0; cu(:,:,0)=-1.d-300; cu(:,:,1)=1.d-300
    cf(:,0,0)=1.d-301; cf(:,0,1)=1.d-299; cf(:,1,:)=0.d0
    if(proc0) call build_rel_cheb_rows(cg,cu,cf,[1,2],[1,1],code)
    call pass_rel_cheb_models(code)
    call check(code==CHEB_SUCCESS,'separate underflow/derivative model initialization')
    v=eval_fit(1,0,(-50.d-300,0.d0),dfpar=d,status=code,value_status=vs,derivative_status=ds)
    call check(code==CHEB_LOG_UNDERFLOW.and.vs==CHEB_LOG_UNDERFLOW.and.ds==CHEB_SUCCESS,'value underflow distinct from representable derivative')
    call check(v==(0.d0,0.d0).and.abs(d)>0.d0,'underflow does not replace representable derivative with zero')
    call check(abs(d/(log(10.d0)*1.d-50)-1.d0)<1.d-11,'independent scaled log-chain derivative')
    v=eval_fit(1,0,(0.d0,0.d0),status=code,value_status=vs,derivative_status=ds)
    call check(code==CHEB_SUCCESS.and.ds==CHEB_NOT_REQUESTED,'value-only adapter does not request derivative')
    v=eval_fit(1,1,(0.d0,0.d0),dfpar=d,status=code,value_status=vs,derivative_status=ds)
    call check(code==CHEB_NO_POSITIVE_SAMPLES.and.vs==code.and.ds==code,'unavailable row statuses preserved')
    call check(ieee_is_nan(real(v,dp)).and.ieee_is_nan(real(d,dp)),'unavailable row never returns zero')
    call rel_cheb_collection_info(nsp,ng,old_calls,bytes,old_builds,evaluations)
    call rel_cheb_distribution_info(old_transfers,old_collectives,volume)
    call cpu_time(started)
    do j=1,10000
       v=eval_fit(1,0,(0.d0,0.d0),dfpar=d,status=code)
       if(code/=CHEB_SUCCESS) error stop 'local derivative timing evaluation failed'
    enddo
    call cpu_time(finished)
    if(proc0) write(*,'(a,es12.4)') 'Local value+derivative seconds/call: ',(finished-started)/10000.d0
    call rel_cheb_collection_info(nsp,ng,calls,bytes,builds,evaluations)
    call rel_cheb_distribution_info(transfers,collectives,volume)
    call check(calls==old_calls.and.builds==old_builds,'repeated analytic adapter evaluations do not initialize or fit')
    call check(transfers==old_transfers.and.collectives==old_collectives,'repeated analytic adapter evaluations do not communicate')
    cu(:,:,0)=-1.d300; cu(:,:,1)=1.d300
    if(proc0) call build_rel_cheb_rows(cg,cu,cf,[1,2],[1,1],code)
    call pass_rel_cheb_models(code)
    v=eval_fit(1,0,(0.d0,0.d0),dfpar=d,status=code,value_status=vs,derivative_status=ds)
    call check(code==CHEB_LOG_UNDERFLOW.and.vs==CHEB_SUCCESS.and.ds==CHEB_LOG_UNDERFLOW,'derivative underflow distinct from representable value')
    call check(abs(v)>0.d0.and.d==(0.d0,0.d0),'derivative underflow keeps physical value')
    cu(:,:,0)=-1.d-20; cu(:,:,1)=1.d-20; cf(:,0,0)=1.d299; cf(:,0,1)=1.d301
    if(proc0) call build_rel_cheb_rows(cg,cu,cf,[1,2],[1,1],code)
    call pass_rel_cheb_models(code)
    v=eval_fit(1,0,(0.d0,0.d0),dfpar=d,status=code,value_status=vs,derivative_status=ds)
    call check(code==CHEB_LOG_DERIVATIVE_FAILURE.and.vs==CHEB_SUCCESS.and.ds==CHEB_LOG_OVERFLOW,'hard derivative failure distinct from value success')
    call check(ieee_is_nan(real(v,dp)).and.ieee_is_nan(real(d,dp)),'hard derivative failure publishes only NaNs')
  end subroutine adapter_checks

  function shape(is,g,z) result(value)
    integer, intent(in) :: is,g
    complex(dp), intent(in) :: z
    complex(dp) :: value
    value=distribution_analyt(is,gamma_rel(is,g,0),z,1,ms(is),vA)
  end function
  function reference(om,harmonic,component,analytic_u,analytic_g) result(value)
    complex(dp), intent(in) :: om
    integer, intent(in) :: harmonic,component
    logical, intent(in) :: analytic_u,analytic_g
    complex(dp) :: value,pole,phi,fg,fu,jb,zb,tensor
    real(dp) :: dg,du,weight,gg
    integer :: row
    dg=gamma_rel(sproc,2,2)-gamma_rel(sproc,1,2)
    du=pparbar_rel(sproc,2,2)-pparbar_rel(sproc,2,1)
    value=(0.d0,0.d0)
    do row=1,ngamma-1
       gg=gamma_rel(sproc,row,1)
       pole=gg*om*vA/kpar-dble(harmonic)*qs(sproc)*vA/(kpar*ms(sproc))
       if(real(pole,dp)**2>gg*gg-1.d0) cycle
       phi=shape(sproc,row,pole)
       ! Closed-form results of the unchanged Gamma FD on exp(-2 Gamma).
       fg=-sinh(2.d0*dg)*phi/dg
       if(row==1) fg=(exp(-2.d0*dg)-1.d0)*phi/dg
       if(analytic_g) fg=-2.d0*phi
       fu=-2.d0*b*pole*phi
       if(.not.analytic_u) fu=-phi*exp(-b*du*du)*sinh(2.d0*b*pole*du)/du
       zb=(kperp*ms(sproc)/(vA*qs(sproc)))*sqrt(gg*gg-1.d0-pole*pole)
       call CBESSJ(zb,abs(harmonic),jb)
       if(component==3) then
          tensor=jb*jb*pole*pole ! T_zz, the unchanged Bessel helper.
       elseif(kperp==0.d0) then
          tensor=(0.d0,0.d0)
          if(abs(harmonic)==1) tensor=.25d0*(gg*gg-1.d0-pole*pole)
       else
          tensor=dble(harmonic*harmonic)*jb*jb/(kperp*ms(sproc)/(vA*qs(sproc)))**2
       endif
       weight=1.d0
       if(row==ngamma-1) weight=.5d0
       value=value-weight*(om*fg+(kpar/vA)*fu)*tensor
    enddo
    value=value*(0.d0,1.d0)*dg*pi*2.d0*pi*(qs(sproc)*vA/(kpar*ms(sproc)))*(ms(sproc)/vA)**3
  end function
  subroutine check(ok,label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    if(ok) then
       passed=passed+1
    else
       failed=failed+1; print '(a,i0,2a)', 'FAIL rank=',iproc,': ',label
    endif
  end subroutine
  subroutine close(actual,expected,tol,label)
    complex(dp), intent(in) :: actual,expected
    real(dp), intent(in) :: tol
    character(*), intent(in) :: label
    call check(ieee_is_finite(real(actual,dp)).and.ieee_is_finite(aimag(actual)),trim(label)//' finite')
    call check(abs(actual-expected)<=tol*max(abs(expected),1.d-10),label)
  end subroutine
end program
