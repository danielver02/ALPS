! Fixed-frequency validation only. Production integration routines are called
! unchanged; no root finder, derivative substitution, or support modification.
program test_rel_susceptibility
  use mpi
  use alps_var
  use alps_io, only : init_param,alps_error_init,read_f0
  use alps_analyt
  use alps_com, only : pass_instructions,pass_distribution,pass_rel_cheb_models,rel_cheb_distribution_info
  use alps_fns, only : derivative_f0,determine_nmax,split_processes,determine_bessel_array, &
       disp,full_integrate,determine_resonances
  use alps_fns_rel, only : int_ee_rel
  use alps_distribution_analyt, only : set_susceptibility_reference
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  integer, parameter :: dp=kind(1.d0),points=5
  integer, parameter :: ti(6)=[1,2,3,1,1,2],tj(6)=[1,2,3,2,3,3]
  character(16), parameter :: methods(0:2)=[character(16)::'analytic','functional','cheb']
  character(16) :: case_kind,study
  character(128) :: label
  integer, allocatable :: degrees(:)
  integer :: status,id,p,method,s,i,j,g,nsr,last,calls,builds,transfers,collectives, &
       initial_calls,initial_builds,initial_transfers,initial_collectives, &
       out,summary,observations,rows_out,queries_out,failed_queries_out,failures,uflow,missing,other,requested, &
       nred,nsample,nrank,nmin,nmax_eff,extrap,total_failures
  integer(int64) :: bytes,volume,evaluations
  real(dp) :: c(2),b,dg,du,norm,lower,upper,slice,weight,phi,peak,scale,magnitude,nan
  real(dp) :: re_err,im_err,scaled_re,scaled_im,eps_re,eps_im,absolute,log_res, &
       fit_abs,fit_rel,interp_rel,value,max_x,refscale
  real(dp), allocatable :: centers(:,:),scales(:,:),left(:,:),right(:,:)
  complex(dp) :: omega,determinant,reference_d(points),functional_d(points),z,fval,derivative,delta
  complex(dp) :: reference_chi(2,3,3,points),functional_chi(2,3,3,points),actual(2,3,3)
  type(cheb_model) :: model
  type(rel_cheb_row_info) :: info

  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0; sproc=0
  call get_command_argument(2,case_kind); call get_command_argument(3,study)
  degrees=[2,8,16,30,40]
  if(study=='refine') degrees=[30,40]
  if(study=='smoke') degrees=[2,30]
  if(proc0) then
     call alps_error_init
     call init_param ! Valid, supported AC_method=1 input; guard is unchanged.
  endif
  call pass_instructions
  if(any(usebM)) error stop 'bimaxwellian proxies excluded from this numerical validation'
  if(nspec/=2.or..not.all(relativistic)) error stop 'fixture requires two numerical relativistic species'
  b=.15d0
  if(case_kind=='iso') b=0.d0
  if(proc0) then
     allocate(f0(nspec,0:nperp,0:npar)); call read_f0
     call derivative_f0
     call determine_param_fit ! Existing functional fitting of this exact completed table.
     c=0.d0
     ! Independent adjacent-panel normalization of the KNOWN ANALYTIC SHAPE.
     ! This is an analytic-continuation control, not a complex table interpolant.
     do s=1,nspec
        dg=gamma_rel(s,1,0)-gamma_rel(s,0,0); du=pparbar_rel(s,0,1)-pparbar_rel(s,0,0)
        norm=0.d0
        do g=0,ngamma
           slice=0.d0
           do j=0,npparbar-1
              if(f0_rel(s,g,j)==-1.d0.or.f0_rel(s,g,j+1)==-1.d0) cycle
              slice=slice+.5d0*gamma_rel(s,g,j)*du*( &
                   exp(-2.d0*gamma_rel(s,g,j)-b*pparbar_rel(s,g,j)**2)+ &
                   exp(-2.d0*gamma_rel(s,g,j+1)-b*pparbar_rel(s,g,j+1)**2))
           enddo
           weight=1.d0
           if(g==0.or.g==ngamma) weight=.5d0
           norm=norm+weight*slice*dg
        enddo
        c(s)=1.d0/(norm*2.d0*pi*(ms(s)/vA)**3)
     enddo
     deallocate(f0)
  endif
  call mpi_bcast(c,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
  call set_susceptibility_reference(c,b)
  ! Legacy coefficient buffers still belong to pass_distribution, although this
  ! harness selects method 2 internally only after reading supported method 1.
  poly_kind=1; poly_order=degrees(1); ACmethod=2
  if(proc0) call initialize_rel_cheb(status)
  call pass_distribution ! Existing data broadcasts and one coefficient transfer.
  if(proc0) then
     allocate(centers(nspec,0:ngamma),scales(nspec,0:ngamma),left(nspec,0:ngamma),right(nspec,0:ngamma))
     open(newunit=out,file='tensor.csv',status='new')
     open(newunit=summary,file='summary.csv',status='new')
     open(newunit=rows_out,file='rows.csv',status='new')
     open(newunit=observations,file='observations.dat',status='new')
     if(study=='queries') then
        open(newunit=queries_out,file='queries.csv',status='new')
        open(newunit=failed_queries_out,file='failed_queries.csv',status='new')
        write(failed_queries_out,'(a)') 'degree,point,species,row,u_re,u_im,status'
        write(queries_out,'(a)') 'degree,point,species,row,u_re,u_im,abs_x,fit_re,fit_im,analytic_re,analytic_im,relative_error'
     endif
     write(out,'(a)') 'case,ngamma,npparbar,degree,point,method,species,i,j,ref_re,ref_im,re,im,delta_re,delta_im,abs_delta,normalized_re,normalized_im,chi_scale'
     write(summary,'(a)') 'case,ngamma,npparbar,degree,point,method,status,max_abs_re,max_abs_im,max_normalized_re,max_normalized_im,tensor_scaled_re,tensor_scaled_im,det_re,det_im,underflow,unavailable,other_failure,extrapolation_queries,max_abs_x'
     write(rows_out,'(a)') 'case,ngamma,npparbar,degree,species,min_degree,max_degree,reduced,sample_limited,rank_limited,max_log_L2,max_fit_abs,max_fit_rel,max_interpolation_rel'
  endif
  nan=ieee_value(0.d0,ieee_quiet_nan)
  do id=1,size(degrees)
     requested=degrees(id)
     if(id>1) then
        ACmethod=2; poly_order=requested
        if(proc0) call initialize_rel_cheb(status)
        call pass_rel_cheb_models(status)
        if(status/=CHEB_SUCCESS) error stop 'coefficient initialization failed'
     endif
     if(proc0) then
        do s=1,nspec
           nred=0; nsample=0; nrank=0; nmin=huge(0); nmax_eff=0
           log_res=0.d0; fit_abs=0.d0; fit_rel=0.d0; interp_rel=0.d0; peak=maxval(f0_rel(s,:,:))
           do g=0,ngamma
              call get_rel_cheb_row(s,g,model,info,status) ! Setup diagnostics only, never in disp.
              if(status/=CHEB_SUCCESS.and.status/=CHEB_REDUCED) error stop 'unexpected unavailable fitting row'
              centers(s,g)=model%center; scales(s,g)=model%scale
              left(s,g)=model%retained_min; right(s,g)=model%retained_max
              if(status==CHEB_REDUCED) nred=nred+1
              if(model%sample_limited) nsample=nsample+1
              if(model%rank_reduced) nrank=nrank+1
              nmin=min(nmin,model%effective_degree); nmax_eff=max(nmax_eff,model%effective_degree)
              log_res=max(log_res,model%residual_norm)
              do j=0,npparbar
                 if(f0_rel(s,g,j)==-1.d0) cycle
                 call cheb_evaluate(model,pparbar_rel(s,g,j),value,status)
                 if(status/=CHEB_SUCCESS) error stop 'real retained-node reconstruction failed'
                 phi=c(s)*exp(-2.d0*gamma_rel(s,g,j)-b*pparbar_rel(s,g,j)**2)
                 fit_abs=max(fit_abs,abs(value-f0_rel(s,g,j)))
                 fit_rel=max(fit_rel,abs(value-f0_rel(s,g,j))/max(f0_rel(s,g,j),1.d-12*peak))
                 interp_rel=max(interp_rel,abs(phi-f0_rel(s,g,j))/max(phi,1.d-12*peak))
              enddo
           enddo
           write(rows_out,'(*(g0,:,","))') trim(case_kind),ngamma,npparbar,requested,s,nmin,nmax_eff,nred,nsample,nrank, &
                log_res,fit_abs,fit_rel,interp_rel
        enddo
     endif
     call rel_cheb_collection_info(nsr,last,initial_calls,bytes,initial_builds,evaluations)
     call rel_cheb_distribution_info(initial_transfers,initial_collectives,volume)
     do p=1,points
        call wave_point(p,omega)
        if(id==1) then
           ACmethod=0
           call response(omega,reference_d(p))
           if(proc0) reference_chi(:,:,:,p)=chi0
           ACmethod=1
           call response(omega,functional_d(p))
           if(proc0) functional_chi(:,:,:,p)=chi0
        endif
        do method=0,2
           ACmethod=method
           failures=0; uflow=0; missing=0; other=0; extrap=0; max_x=0.d0
           if(method==2) then
              if(proc0) call preflight(omega,failures,uflow,missing,other,extrap,max_x)
              call mpi_bcast(failures,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
           endif
           status=CHEB_SUCCESS
           if(failures/=0) status=CHEB_EVALUATION_FAILURE
           if(status==CHEB_SUCCESS) then
              select case(method)
              case(0)
                 determinant=reference_d(p)
                 if(proc0) actual=reference_chi(:,:,:,p)
              case(1)
                 determinant=functional_d(p)
                 if(proc0) actual=functional_chi(:,:,:,p)
              case(2)
                 call response(omega,determinant)
                 if(proc0) actual=chi0
              end select
           else
              determinant=cmplx(nan,nan,dp); actual=determinant
           endif
           if(proc0) then
              re_err=0.d0; im_err=0.d0; eps_re=0.d0; eps_im=0.d0; scaled_re=0.d0; scaled_im=0.d0
              do s=1,nspec
                 magnitude=maxval(abs(reference_chi(s,:,:,p)))
                 scale=1.d-6*magnitude
                 do i=1,6
                    z=reference_chi(s,ti(i),tj(i),p)
                    fval=actual(s,ti(i),tj(i)); delta=fval-z
                    if(status==CHEB_SUCCESS) then
                       if(.not.ieee_is_finite(real(fval,dp)).or..not.ieee_is_finite(aimag(fval))) &
                            error stop 'nonfinite susceptibility output'
                       re_err=max(re_err,abs(real(delta,dp))); im_err=max(im_err,abs(aimag(delta)))
                       eps_re=max(eps_re,abs(real(delta,dp))/max(abs(real(z,dp)),scale))
                       eps_im=max(eps_im,abs(aimag(delta))/max(abs(aimag(z)),scale))
                       scaled_re=max(scaled_re,abs(real(delta,dp))/magnitude)
                       scaled_im=max(scaled_im,abs(aimag(delta))/magnitude)
                    endif
                    write(out,'(*(g0,:,","))') trim(case_kind),ngamma,npparbar,requested,p,trim(methods(method)),s,ti(i),tj(i), &
                         real(z,dp),aimag(z),real(fval,dp),aimag(fval),real(delta,dp),aimag(delta),abs(delta), &
                         abs(real(delta,dp))/max(abs(real(z,dp)),scale),abs(aimag(delta))/max(abs(aimag(z)),scale),scale
                    if(status==CHEB_SUCCESS) then
                       write(label,'(a,a,i0,a,i0,a,i0,a,i0,a,i0)') trim(methods(method)),'_d',requested,'_p',p,'_s',s, &
                            '_i',ti(i),'_j',tj(i)
                       write(observations,'(a,4es26.17e3)') trim(label)//' ',real(fval,dp),aimag(fval),1.d-13,1.d-12
                    endif
                 enddo
              enddo
              if(status/=CHEB_SUCCESS) then
                 re_err=nan; im_err=nan; eps_re=nan; eps_im=nan; scaled_re=nan; scaled_im=nan
              endif
              write(summary,'(*(g0,:,","))') trim(case_kind),ngamma,npparbar,requested,p,trim(methods(method)),status, &
                   re_err,im_err,eps_re,eps_im,scaled_re,scaled_im,real(determinant,dp),aimag(determinant), &
                   uflow,missing,other,extrap,max_x
              write(label,'(a,a,i0,a,i0)') trim(methods(method)),'_status_d',requested,'_p',p
              write(observations,'(a,4es26.17e3)') trim(label)//' ',dble(status),0.d0,1.d-13,1.d-12
           endif
           call rel_cheb_collection_info(nsr,last,calls,bytes,builds,evaluations)
           call rel_cheb_distribution_info(transfers,collectives,volume)
           if(calls/=initial_calls.or.builds/=initial_builds.or.transfers/=initial_transfers.or. &
                collectives/=initial_collectives) error stop 'response evaluation refit or retransmitted coefficients'
        enddo
     enddo
     if(proc0) write(*,'(a,a,3i5)') 'Completed fixed frequencies: ',trim(case_kind),ngamma,npparbar,requested
  enddo
  if(proc0) then
     close(out); close(summary); close(rows_out); close(observations)
     if(study=='queries') then
        close(queries_out); close(failed_queries_out)
     endif
     write(*,'(a)') 'PASS: finite-frequency observations; no roots, refitting, or coefficient communication during response'
  endif
  call underflow_contract
  call mpi_finalize(ierror)
contains
  subroutine underflow_contract
    real(dp) :: cg(2,0:1,0:1),cu(2,0:1,0:1),cf(2,0:1,0:1)
    complex(dp) :: v
    integer :: code
    cg=1.5d0; cu(:,:,0)=-1.d0; cu(:,:,1)=1.d0
    cf(:,0,0)=1.d-300; cf(:,0,1)=1.d-301; cf(:,1,:)=0.d0
    ACmethod=2
    if(proc0) call build_rel_cheb_rows(cg,cu,cf,[1,2],[1,1],code)
    call pass_rel_cheb_models(code)
    if(code/=CHEB_SUCCESS) error stop 'underflow contract initialization'
    v=eval_fit(1,0,(1000.d0,0.d0),status=code)
    if(code/=CHEB_LOG_UNDERFLOW.or.v/=(0.d0,0.d0)) error stop 'valid underflow was not preserved'
    v=eval_fit(1,1,(0.d0,0.d0),status=code)
    if(code/=CHEB_NO_POSITIVE_SAMPLES.or..not.ieee_is_nan(real(v,dp))) error stop 'unavailable row confused with underflow'
    v=eval_fit(1,0,cmplx(nan,0.d0,dp),status=code)
    if(code/=CHEB_INVALID_COORDINATES.or..not.ieee_is_nan(real(v,dp))) error stop 'invalid evaluation confused with underflow'
    v=eval_fit(1,0,(-2000.d0,0.d0),status=code)
    if(code/=CHEB_LOG_OVERFLOW.or..not.ieee_is_nan(real(v,dp))) error stop 'overflow confused with underflow'
    v=eval_fit(1,0,(0.d0,0.d0),status=code)
    if(code/=CHEB_SUCCESS) error stop 'underflow invalidated a valid model'
    if(proc0) write(*,'(a)') 'PASS: underflow/unavailable/invalid/overflow statuses remain distinct'
  end subroutine underflow_contract
  subroutine wave_point(point,om)
    integer, intent(in) :: point
    complex(dp), intent(out) :: om
    kpar=.1d0; kperp=.001d0
    select case(point)
    case(1); om=(10.d0,-.001d0)
    case(2); om=(.062713d0,.0001d0)
    case(3); om=(.062713d0,-.0001d0)
    case(4); om=(.062713d0,-.005d0); kperp=.02d0
    case(5); om=(.12d0,-.001d0); kpar=.2d0; kperp=.05d0
    end select
    call determine_nmax
    sproc=0
    if(nproc>1) then
       call split_processes
       if(sproc/=0) call determine_bessel_array
    endif
  end subroutine
  subroutine response(om,d)
    complex(dp), intent(in) :: om
    complex(dp), intent(out) :: d
    complex(dp) :: raw(2,3,3),terms(6),eps(3,3)
    integer :: is,nn,mode
    logical :: rp,rm
    if(nproc>1) then
       d=disp(om) ! Actual production susceptibility assembly and MPI reduction.
       return
    endif
    ! disp reserves rank 0 for coordination, so a singleton communicator cannot
    ! use it for kinetic work. Assemble its six components from the SAME public
    ! production full_integrate/int_ee_rel calls; never duplicate quadrature.
    raw=(0.d0,0.d0)
    do is=1,nspec
       sproc=is; nlim=[0,nmax(is)]; call determine_bessel_array
       terms=(0.d0,0.d0)
       do nn=0,nmax(is)
          call determine_resonances(om,nn,rp,rm)
          do mode=1,6
             if(nn==0.and.(mode==1.or.mode==4.or.mode==5)) cycle
             terms(mode)=terms(mode)+full_integrate(om,nn,mode,rp)
             if(nn/=0) terms(mode)=terms(mode)+full_integrate(om,-nn,mode,rm)
          enddo
       enddo
       terms(3)=terms(3)+int_ee_rel(om)
       do mode=1,6
          raw(is,ti(mode),tj(mode))=terms(mode)*ns(is)*qs(is)
       enddo
    enddo
    chi0=raw/(om*om*vA*vA)
    chi0(:,2,1)=-chi0(:,1,2); chi0(:,3,1)=chi0(:,1,3); chi0(:,3,2)=-chi0(:,2,3)
    eps=sum(raw,dim=1)
    eps(1,1)=eps(1,1)+(om*vA)**2
    eps(2,2)=eps(2,2)+(om*vA)**2
    eps(3,3)=eps(3,3)+(om*vA)**2
    wave=eps; wave(1,1)=eps(1,1)-kpar*kpar; wave(2,2)=eps(2,2)-kpar*kpar-kperp*kperp
    wave(3,3)=eps(3,3)-kperp*kperp; wave(1,3)=eps(1,3)+kpar*kperp
    wave(2,1)=-wave(1,2); wave(3,1)=wave(1,3); wave(3,2)=-wave(2,3)
    d=wave(1,1)*(wave(2,2)*wave(3,3)+wave(2,3)**2)+2.d0*wave(1,2)*wave(2,3)*wave(1,3)- &
         wave(1,3)**2*wave(2,2)+wave(1,2)**2*wave(3,3)
    sproc=0
  end subroutine
  subroutine preflight(om,nfail,nu,nmissing,nother,nextrap,xmax)
    complex(dp), intent(in) :: om
    integer, intent(out) :: nfail,nu,nmissing,nother,nextrap
    real(dp), intent(out) :: xmax
    integer :: is,nn,sign_n,ig
    complex(dp) :: pole
    real(dp) :: step
    logical :: rp,rm
    nfail=0; nu=0; nmissing=0; nother=0; nextrap=0; xmax=0.d0
    if(aimag(om)>0.d0) return
    do is=1,nspec
       sproc=is; step=pparbar_rel(is,0,1)-pparbar_rel(is,0,0)
       do nn=0,nmax(is)
          call determine_resonances(om,nn,rp,rm)
          do sign_n=1,-1,-2
             if(nn==0.and.sign_n==-1) cycle
             if(sign_n==1.and..not.rp) cycle
             if(sign_n==-1.and..not.rm) cycle
             do ig=1,ngamma-1
                pole=(gamma_rel(is,ig,0)*om-sign_n*nn*qs(is)/ms(is))*vA/kpar
                if(real(pole,dp)**2>gamma_rel(is,ig,0)**2-1.d0) cycle
                call query(is,ig+1,pole,nfail,nu,nmissing,nother,nextrap,xmax)
                if(ig==1) then
                   call query(is,ig,pole,nfail,nu,nmissing,nother,nextrap,xmax)
                else
                   call query(is,ig-1,pole,nfail,nu,nmissing,nother,nextrap,xmax)
                endif
                call query(is,ig,pole+step,nfail,nu,nmissing,nother,nextrap,xmax)
                call query(is,ig,pole-step,nfail,nu,nmissing,nother,nextrap,xmax)
             enddo
          enddo
       enddo
    enddo
    sproc=0
  end subroutine preflight
  subroutine query(is,ig,z,nfail,nu,nmissing,nother,nextrap,xmax)
      integer, intent(in) :: is,ig
      integer, intent(inout) :: nfail,nu,nmissing,nother,nextrap
      real(dp), intent(inout) :: xmax
      complex(dp), intent(in) :: z
      complex(dp) :: f,want
      integer :: code
      f=eval_fit(is,ig,z,status=code)
      xmax=max(xmax,abs((z-centers(is,ig))/scales(is,ig)))
      if(real(z,dp)<left(is,ig).or.real(z,dp)>right(is,ig)) nextrap=nextrap+1
      if(code==CHEB_SUCCESS) then
         if(study=='queries') then
            want=c(is)*exp(-2.d0*gamma_rel(is,ig,0)-b*z*z)
            write(queries_out,'(*(g0,:,","))') requested,p,is,ig,real(z,dp),aimag(z), &
                 abs((z-centers(is,ig))/scales(is,ig)),real(f,dp),aimag(f),real(want,dp),aimag(want), &
                 abs(f-want)/max(abs(want),tiny(1.d0))
         endif
         return
      endif
      if(study=='queries') write(failed_queries_out,'(*(g0,:,","))') requested,p,is,ig,real(z,dp),aimag(z),code
      nfail=nfail+1
      select case(code)
      case(CHEB_LOG_UNDERFLOW); nu=nu+1
      case(CHEB_NO_SAMPLES,CHEB_NO_POSITIVE_SAMPLES,CHEB_INVALID_MODEL); nmissing=nmissing+1
      case default; nother=nother+1
      end select
  end subroutine query
end program
