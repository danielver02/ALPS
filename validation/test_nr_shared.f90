program test_nr_shared
 use mpi
 use alps_var
 use alps_io, only:init_param,alps_error_init
 use alps_analyt
 use alps_com, only:pass_instructions,pass_distribution
 use alps_fns, only:landau_integrate,determine_bessel_array
 use, intrinsic::ieee_arithmetic
 use legacy_glls_reference
 implicit none
 integer::s,r,j,k,code,n,rep,unit,comparison_unit,checks=0,failures=0,total
 double precision::c(0:4),nan,coefficient_checksum
 double complex::z,v,d,old,od,want
 character(16)::option
 type(cheb_model)::model
 call mpi_init(ierror)
 call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror);call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
 proc0=iproc==0
 if(proc0) then
  call alps_error_init;call init_param
 endif
 call get_command_argument(2,option)
 call pass_instructions
 call ensure_historical_storage
 call check(allocated(historical_coefficients),"validation-owned coefficient storage")
 if(proc0)allocate(f0(nspec,0:nperp,0:npar))
 if(proc0)then
 open(newunit=comparison_unit,file="comparisons.csv",status="new")
 write(comparison_unit,'(a)') "phase,representation,species,query,legacy_re,legacy_im,shared_re,shared_im,delta_re,delta_im,derivative_abs_error,legacy_d_re,legacy_d_im,shared_d_re,shared_d_im,derivative_delta_re,derivative_delta_im"
 endif
 ! Test validation-only historical transfer in all-degree-one and mixed cases.
 do k=1,2
  poly_order=1;if(k==2)poly_order(2)=4
  if(allocated(historical_coefficients))deallocate(historical_coefficients)
  if(proc0.or.k==2)allocate(historical_coefficients(nspec,0:nperp,0:maxval(poly_order)))
  logfit=.true.
  if(proc0) then
   historical_coefficients=0.d0
   do s=1,nspec
    do r=0,nperp
     do j=0,npar
      pp(s,r,j,2)=-2.d0+4.d0*j/npar;pp(s,r,j,1)=.1d0*r
     enddo
     historical_coefficients(s,r,0)=-2.d0;historical_coefficients(s,r,1)=.3d0
     if(poly_order(s)>1)historical_coefficients(s,r,2)=.02d0
     do j=0,npar
      z=cmplx(pp(s,r,j,2)/2.d0,0.d0,kind(1.d0))
      f0(s,r,j)=10.d0**(-2.d0+.3d0*real(z))
      if(poly_order(s)>1)f0(s,r,j)=f0(s,r,j)*10.d0**(.02d0*(2.d0*real(z)**2-1.d0))
     enddo
    enddo
   enddo
   call legacy_determine_param_fit
  else
   ! Worker fixture storage is validation-owned; production has no coefficient buffer.
   if(.not.allocated(historical_coefficients))allocate(historical_coefficients(nspec,0:nperp,0:maxval(poly_order)))
   historical_coefficients=0.d0
  endif
  call pass_distribution
  call initialize_historical_nr(code);call check(code==CHEB_SUCCESS,'local import')
  call check(maxval(abs(historical_coefficients(:,:,0)+2.d0))<1.d-12,"broadcast constant coefficient")
  call check(maxval(abs(historical_coefficients(:,:,1)-.3d0))<1.d-12,"broadcast linear coefficient")
  coefficient_checksum=sum(historical_coefficients)
  do rep=0,1
   logfit=rep==1
   call initialize_historical_nr(code);call check(code==CHEB_SUCCESS,'representation import')
   do s=1,nspec
    do j=0,6
     z=cmplx(-1.5d0+.5d0*j,.05d0,kind(1.d0))
     call select_nr_test_evaluation(.false.);old=legacy_forward_evaluate(s,1,z,poly_order(s),historical_coefficients(s,1,0:poly_order(s)),od)
     call select_nr_test_evaluation(.true.);v=eval_fit(s,1,z,dfpar=d,status=code)
     call check(code==CHEB_SUCCESS,'shared status')
     call check(abs(v-old)<2.d-12*max(1.d0,abs(old)),'same coefficient value')
     call check(abs(d-od)<2.d-12*max(1.d0,abs(od)),'same coefficient derivative')
     if(proc0)write(comparison_unit,'(*(g0,:,","))')k,rep,s,j,real(old),aimag(old),real(v),aimag(v),real(v-old),aimag(v-old),abs(d-od),real(od),aimag(od),real(d),aimag(d),real(d-od),aimag(d-od)
    enddo
   enddo
  enddo
  call check(sum(historical_coefficients)==coefficient_checksum,'coefficients immutable')
 enddo
 if(proc0)open(newunit=unit,file='policies.csv',status='new')
 if(proc0)write(unit,'(a)')'query,legacy_re,legacy_im,shared_re,shared_im,status'
 do j=1,3
  z=(2.01d0,.01d0);if(j==2)z=(0.d0,2.01d0)
  if(j==3)then
   historical_log_max=-10.d0;z=(0.d0,0.d0)
  endif
  call select_nr_test_evaluation(.false.);old=eval_fit(1,1,z)
  call select_nr_test_evaluation(.true.);v=eval_fit(1,1,z,status=code)
  call check(old==(0.d0,0.d0).and.abs(v)>0.and.code==CHEB_SUCCESS,'intentional suppression removed')
  if(proc0)write(unit,'(*(g0,:,","))')j,real(old),aimag(old),real(v),aimag(v),code
 enddo
 if(proc0)then
 close(unit);close(comparison_unit)
 endif
 nan=ieee_value(0.d0,ieee_quiet_nan)
 do n=0,4
  c=0.d0;c(0)=.2d0;if(n>0)c(n)=.3d0
  call cheb_from_coefficients(c(:n),n,1.d0,2.d0,CHEB_DIRECT,model,code)
  call check(code==CHEB_SUCCESS.and.model%coefficients_imported.and.model%rank==-1,'import metadata')
  call check(ieee_is_nan(model%residual_norm).and.model%fit_attempts==0,'no fabricated fit diagnostics')
  call cheb_evaluate(model,(1.d0,0.d0),v,code,d)
  if(n==0)call check(d==(0.d0,0.d0),'constant derivative')
  if(n==1)call check(abs(d-(.15d0,0.d0))<1.d-15,'linear physical derivative')
  if(n==2)call check(abs(v-(-.1d0,0.d0))<1.d-15,'full c0 convention')
  z=(.2d0,.1d0);want=(0.d0,0.d0)
  if(n==1)want=(.15d0,0.d0)
  if(n==2)want=.6d0*z
  if(n==3)want=.15d0*(12.d0*z*z-3.d0)
  if(n==4)want=.15d0*(32.d0*z**3-16.d0*z)
  call cheb_evaluate(model,1.d0+2.d0*z,v,code,d)
  call check(code==CHEB_SUCCESS.and.abs(d-want)<1.d-13,'independent complex imported derivative')
 enddo
 call cheb_from_coefficients([nan],0,0.d0,1.d0,CHEB_DIRECT,model,code)
 call check(code==CHEB_NONFINITE_COEFFICIENTS,'nonfinite coefficients')
 call cheb_from_coefficients([1.d0],0,0.d0,0.d0,CHEB_DIRECT,model,code)
 call check(code==CHEB_INVALID_NORMALIZATION,'invalid scale')
 call cheb_from_coefficients([1.d0],-1,0.d0,1.d0,CHEB_DIRECT,model,code)
 call check(code==CHEB_INVALID_DEGREE,'invalid degree')
 call cheb_from_coefficients([1.d0],0,0.d0,1.d0,99,model,code)
 call check(code==CHEB_INVALID_REPRESENTATION,'invalid representation')
 call cheb_from_coefficients([400.d0],0,0.d0,1.d0,CHEB_LOG10,model,code)
 call cheb_evaluate(model,(0.d0,0.d0),v,code)
 call check(code==CHEB_LOG_OVERFLOW.and.ieee_is_nan(real(v)),'explicit overflow')
 call cheb_from_coefficients([-400.d0],0,0.d0,1.d0,CHEB_LOG10,model,code)
 call cheb_evaluate(model,(0.d0,0.d0),v,code)
 call check(code==CHEB_LOG_UNDERFLOW.and.v==(0.d0,0.d0),'explicit underflow')
 ! Landau legacy gate and corrected shared path, exposed only in scratch.
 sproc=1;kpar=.1d0;kperp=0.d0;nlim=[0,0]
 call determine_bessel_array
 logfit=.true.;historical_log_max=18.d0
 call initialize_historical_nr(code)
 if(proc0)then
 open(newunit=unit,file='stencil_policy.csv',status='new')
 write(unit,'(a)')'selector,legacy_re,legacy_im,shared_re,shared_im'
 endif
 do j=0,1
 poly_derivative_type=j
 call select_nr_test_evaluation(.false.);old=landau_integrate((0.d0,-.202d0),0,3)
 call select_nr_test_evaluation(.true.);v=landau_integrate((0.d0,-.202d0),0,3)
 call check(old==(0.d0,0.d0).and.abs(v)>1.d-6,'actual Landau disk gate difference')
 if(proc0)write(unit,'(*(g0,:,","))')j,real(old),aimag(old),real(v),aimag(v)
 enddo
 if(proc0)close(unit)
 ! A valid all-zero direct distribution integrates to zero without a reset.
 historical_coefficients=0.d0;logfit=.false.
 call initialize_historical_nr(code);call select_nr_test_evaluation(.true.)
 v=landau_integrate((0.d0,-.01d0),0,3)
 call check(v==(0.d0,0.d0),'valid zero integrand gives zero contour')
 historical_coefficients(:,:,0)=2.d0;historical_coefficients(:,:,1)=.3d0
 historical_coefficients(:,3,:)=0.d0
 call initialize_historical_nr(code)
 v=landau_integrate((0.d0,-.01d0),0,3)
 want=(0.d0,0.d0)
 do r=1,nperp
 if(r==3)cycle
 coefficient_checksum=1.d0
 if(r>=nperp-1)coefficient_checksum=.5d0
 want=want+coefficient_checksum*pp(sproc,r,1,1)*.15d0*.01d0
 enddo
 want=want*cmplx(0.d0,1.d0,kind(1.d0))*.1d0*2*pi*pi
 call check(abs(v-want)<1.d-10.and.abs(v)>0,'later genuine zero preserves accumulated contour')
 call direct_glls_checks
 call release_nr_cheb
 call eval_nr_cheb(1,1,(0.d0,0.d0),v,code)
 call check(code==CHEB_INVALID_MODEL,'released collection')
 if(option=='fatal')then
 call select_nr_test_evaluation(.true.)
 v=eval_fit(1,1,(0.d0,0.d0))
 error stop 'missing NR model silently accepted'
 endif
 call mpi_allreduce(failures,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
 write(*,'(a,i0,a,i0,a,i0)')'NR shared rank=',iproc,' passed=',checks,' failed=',failures
 if(total/=0)call mpi_abort(MPI_COMM_WORLD,99,ierror)
 call mpi_finalize(ierror)
contains
 subroutine direct_glls_checks
 integer::s,r,j,k,code,file
 double precision::u
 double complex::z,v,d,f,df,legacy,ld
 logfit=.false.;poly_order=4
 if(proc0)then
 do s=1,nspec
 do r=0,nperp
 do j=0,npar
 u=pp(s,r,j,2)
 f0(s,r,j)=2.d0+.2d0*pp(s,r,j,1)+.3d0*u-.2d0*u*u+.05d0*u**4
 enddo
 enddo
 enddo
 call legacy_determine_param_fit
 endif
 call pass_distribution
 call initialize_historical_nr(code);call check(code==CHEB_SUCCESS,'direct GLLS models')
 if(proc0)then
 open(newunit=file,file='direct_glls.csv',status='new')
 write(file,'(a)')'species,query,re,im,expected_re,expected_im,derivative_re,derivative_im,expected_derivative_re,expected_derivative_im,status'
 endif
 do s=1,nspec
 do k=1,3
 z=(-1.2d0,.05d0)
 if(k==2)z=(.4d0,-.1d0)
 if(k==3)z=(2.1d0,.3d0)
 f=2.d0+.2d0*pp(s,1,0,1)+.3d0*z-.2d0*z*z+.05d0*z**4
 df=.3d0-.4d0*z+.2d0*z**3
 call select_nr_test_evaluation(.true.)
 v=eval_fit(s,1,z,dfpar=d,status=code)
 call check(code==CHEB_SUCCESS.and.abs(v-f)<1.d-11,'independent direct GLLS value')
 call check(abs(d-df)<1.d-11,'independent direct GLLS derivative')
 call select_nr_test_evaluation(.false.)
 legacy=eval_fit(s,1,z,dfpar=ld)
 if(k<3)then
 call check(abs(v-legacy)<1.d-11.and.abs(d-ld)<1.d-11,'direct fitted inside legacy agreement')
 else
 call check(legacy==(0.d0,0.d0).and.abs(v)>0.d0,'direct fitted extrapolation policy')
 endif
 if(proc0)write(file,'(*(g0,:,","))')s,k,real(v),aimag(v),real(f),aimag(f),real(d),aimag(d),real(df),aimag(df),code
 enddo
 enddo
 if(proc0)close(file)
 end subroutine
 subroutine check(ok,label)
 logical,intent(in)::ok
 character(*),intent(in)::label
 if(ok)then
 checks=checks+1
 else
 failures=failures+1;print *,'FAIL ',label
 endif
 end subroutine
end program
