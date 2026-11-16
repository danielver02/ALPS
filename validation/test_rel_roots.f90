! Controlled root validation, using the production rtsec and om_scan routines.
program test_rel_roots
  use mpi
  use alps_var
  use alps_io, only: init_param,alps_error_init,read_f0
  use alps_analyt
  use alps_com, only: pass_instructions,pass_distribution
  use alps_fns, only: derivative_f0,determine_nmax,split_processes,determine_bessel_array, &
       rtsec,disp,om_scan
  use alps_distribution_analyt, only: set_susceptibility_reference
  use rel_root_diagnostics, only: root_diagnostic_open,root_context,root_tensor_record,root_reference_close,root_accepted
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer,parameter::dp=kind(1.d0)
  character(32)::case_kind,arg
  integer::method,code,s,g,j,unit,flag,nref,i,converged(6),residual_ok(6)
  real(dp)::c(2),b,norm,slice,dg,du,w,omega_re,omega_im,reference_k(6),reference_perp(6)
  complex(dp)::root,initial,d,refs(6)
  logical::has_seed,accepted
  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0;sproc=0
  call get_command_argument(2,case_kind);call get_command_argument(3,arg);read(arg,*)method
  if(method<0.or.method>2) error stop 'test method must be 0,1,2'
  if(proc0) then
     call alps_error_init;call init_param ! Normal, supported method-1 input.
  endif
  call pass_instructions
  if(any(usebM).or..not.all(relativistic).or.nspec/=2) error stop 'root fixture requires two numerical relativistic species'
  if(secant_method/=1) error stop 'this harness observes the existing rtsec pathway'
  b=.15d0;if(case_kind=='iso') b=0.d0
  ACmethod=method;poly_order=2
  c=0.d0
  if(proc0) then
     allocate(f0(nspec,0:nperp,0:npar));call read_f0;call derivative_f0
     do s=1,nspec
        dg=gamma_rel(s,1,0)-gamma_rel(s,0,0);du=pparbar_rel(s,0,1)-pparbar_rel(s,0,0)
        norm=0.d0
        do g=0,ngamma
           slice=0.d0
           do j=0,npparbar-1
              if(f0_rel(s,g,j)==-1.d0.or.f0_rel(s,g,j+1)==-1.d0) cycle
              slice=slice+.5d0*gamma_rel(s,g,j)*du*(exp(-2.d0*gamma_rel(s,g,j)-b*pparbar_rel(s,g,j)**2)+ &
                   exp(-2.d0*gamma_rel(s,g,j+1)-b*pparbar_rel(s,g,j+1)**2))
           enddo
           w=1.d0;if(g==0.or.g==ngamma) w=.5d0
           norm=norm+w*slice*dg
        enddo
        c(s)=1.d0/(norm*2.d0*pi*(ms(s)/vA)**3)
     enddo
     call qualify_table
     open(newunit=unit,file='completed.bin',access='stream',form='unformatted',status='new')
     write(unit) gamma_rel,pparbar_rel,f0_rel
     close(unit)
     call initialize_rel_cheb(code)
     if(code/=CHEB_SUCCESS) error stop 'root model initialization failed'
     call determine_param_fit
     deallocate(f0)
  endif
  call mpi_bcast(c,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
  call set_susceptibility_reference(c,b)
  ! Existing pass_distribution retains a legacy polynomial buffer even when
  ! relativistic rows are stored separately. It is not used for evaluation.
  call pass_distribution
  call root_diagnostic_open(c,b,method)
  initial=wroots(2) ! Existing fixture's ordinary-mode initial guess; not invented.
  if(method==0.and.case_kind=='aniso') then
     ! Follow the same mode into anisotropy using the recovered isotropic root,
     ! rather than inventing a new frequency or changing the solver.
     inquire(file='seed.dat',exist=has_seed)
     if(has_seed) then
        open(newunit=unit,file='seed.dat',status='old')
        read(unit,*) reference_k(1),reference_perp(1),omega_re,omega_im,converged(1),residual_ok(1)
        close(unit)
        if(converged(1)/=1.or.residual_ok(1)/=1) error stop 'unqualified isotropic branch seed'
        initial=cmplx(omega_re,omega_im,dp)
     endif
  endif
  nroots=1
  call setup_wave(.1d0,.001d0)
  if(method/=0) then
     call read_references
     ! Preflight EVERY reference point before allowing any Chebyshev iteration.
     do i=1,nref
        call setup_wave(reference_k(i),reference_perp(i))
        call root_context(1);d=disp(refs(i));call root_tensor_record(refs(i),d)
     enddo
     if(converged(1)/=1.or.residual_ok(1)/=1) error stop 'unqualified initial reference root: preflight recorded, recovery withheld'
     ! Put BOTH native rtsec seeds below the near-degenerate transverse pair.
     ! An upper-side seed reaches its y-polarized partner instead. Keep rtsec's
     ! own separation, stopping rules and branch behavior entirely unchanged.
     initial=refs(1)*(1.d0-2.d0*D_prec)
     call setup_wave(.1d0,.001d0)
  else
     call root_context(1);d=disp(initial);call root_tensor_record(initial,d)
  endif
  flag=0;root=rtsec(disp,initial,1,flag)
  d=disp(root);call root_scan_result(root,d)
  if(.not.ieee_is_finite(real(root)).or..not.ieee_is_finite(aimag(root))) error stop 'nonfinite root'
  accepted=root_accepted(root,d)
  call mpi_bcast(accepted,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  wroots(1)=root
  ! Existing scan predictor/jump/branch-selection logic remains in om_scan.
  if(allocated(scan)) deallocate(scan)
  allocate(scan(1))
  scan(1)%type_s=4;scan(1)%n_out=5;scan(1)%n_res=1;scan(1)%diff=.002d0
  scan(1)%log_scan=.false.;scan(1)%eigen_s=.false.;scan(1)%heat_s=.false.
  if(accepted) call om_scan(1)
  if(method==0) then
     call root_reference_close
     call read_references
     do i=1,nref
        call setup_wave(reference_k(i),reference_perp(i))
        call root_context(1);d=disp(refs(i));call root_tensor_record(refs(i),d)
     enddo
  endif
  if(proc0) print '(a)', 'COMPLETE: controlled root attempts recorded; convergence/branch acceptance checked separately'
  call mpi_finalize(ierror)
contains
  subroutine setup_wave(par,perp)
    real(dp),intent(in)::par,perp
    kpar=par;kperp=perp;sproc=0
    call determine_nmax
    if(nproc>1) then
       call split_processes
       if(sproc/=0) call determine_bessel_array
    endif
  end subroutine
  subroutine read_references
    integer::ios
    nref=0
    if(proc0) then
       open(newunit=unit,file='reference.dat',status='old')
       do i=1,6
          read(unit,*,iostat=ios) reference_k(i),reference_perp(i),omega_re,omega_im,converged(i),residual_ok(i)
          if(ios/=0) exit
          refs(i)=cmplx(omega_re,omega_im,dp);nref=i
       enddo
       close(unit)
    endif
    call mpi_bcast(nref,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    if(nref<1) error stop 'no reference attempts available'
    call mpi_bcast(refs,6,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
    call mpi_bcast(reference_k,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    call mpi_bcast(reference_perp,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    call mpi_bcast(converged,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call mpi_bcast(residual_ok,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  end subroutine
  subroutine qualify_table
    real(dp)::phi,err,peak,l1,base,maxerr,insideerr,coreerr,sourceperp,rmax,du,dg,weight
    integer::s,g,j,source_count,valid
    open(newunit=unit,file='table_accuracy.csv',status='new')
    write(unit,'(a)') 'species,C,b,valid_samples,inside_source,max_relative,max_inside_relative,max_core_relative,weighted_L1_relative'
    do s=1,nspec
       maxerr=0.d0;insideerr=0.d0;coreerr=0.d0;l1=0.d0;base=0.d0;valid=0;source_count=0
       peak=maxval(f0_rel(s,:,:));rmax=maxval(pp(s,:,:,1))
       do g=0,ngamma
          do j=0,npparbar
             if(f0_rel(s,g,j)==-1.d0) cycle
             phi=c(s)*exp(-2.d0*gamma_rel(s,g,j)-b*pparbar_rel(s,g,j)**2)
             err=abs(f0_rel(s,g,j)-phi)/max(phi,1.d-12*peak)
             valid=valid+1;maxerr=max(maxerr,err)
             if(phi>1.d-6*peak) coreerr=max(coreerr,err)
             sourceperp=sqrt(max(0.d0,gamma_rel(s,g,j)**2-1.d0-pparbar_rel(s,g,j)**2))*ms(s)/vA
             if(sourceperp<=rmax*(1.d0+1.d-14)) then
                insideerr=max(insideerr,err);source_count=source_count+1
             endif
             weight=gamma_rel(s,g,j)
             l1=l1+weight*abs(f0_rel(s,g,j)-phi);base=base+weight*phi
          enddo
       enddo
       ! Discrete weighted L1 on supplied samples, not a new physical quadrature.
       write(unit,'(*(g0,:,","))') s,c(s),b,valid,source_count,maxerr,insideerr,coreerr,l1/base
    enddo
    close(unit)
  end subroutine
end program
