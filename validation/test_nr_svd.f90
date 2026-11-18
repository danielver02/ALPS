program test_nr_svd
 use mpi
 use alps_var
 use alps_io, only:init_param,alps_error_init
 use alps_analyt
 use alps_com, only:pass_instructions,pass_nr_cheb_models,nr_cheb_distribution_info
 use, intrinsic::ieee_arithmetic
 use, intrinsic::iso_fortran_env,only:int64
 use legacy_glls_reference, only: legacy_determine_param_fit,reference_coefficients
 implicit none
 type(cheb_model)::model
 type(nr_cheb_row_info)::info
 integer::phase,policy,s,r,j,code,unit,query_unit,checks=0,failures=0,total,expected_code,vs,ds
 integer::collection_species,rows,builds,fits,evals,setup,oldbuilds,oldfits,tx,collect,oldtx,oldcollect,at,n
 integer,allocatable::meta(:,:)
 real(kind(1.d0)),allocatable::mapping(:,:),coef(:,:),legacy(:,:,:)
 integer(int64),allocatable::bits(:),small(:),large(:)
 double precision::u,b,nan,oldsum
 double complex::z,v,d,want,wantd,old,oldd
 call mpi_init(ierror)
 call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror);call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
 proc0=iproc==0
 if(proc0)then
 call alps_error_init;call init_param
 endif
 call pass_instructions
 ! An inactive/missing model must never reach historical coefficient data.
 call release_nr_cheb
 z=(0.d0,0.d0)
 v=eval_fit(1,0,z,dfpar=d,status=code,value_status=vs,derivative_status=ds)
 call check(code==CHEB_INVALID_MODEL.and.vs==code.and.ds==code,'missing model statuses')
 call check(.not.ieee_is_finite(real(v)).and..not.ieee_is_finite(real(d)),'missing model NaN outputs')
 v=eval_fit(1,0,z,status=code,derivative_status=ds)
 call check(code==CHEB_INVALID_MODEL.and.ds==CHEB_NOT_REQUESTED,'missing value-only model')
 poly_kind(1)=2
 v=eval_fit(1,0,z,status=code)
 call check(code==CHEB_INVALID_REPRESENTATION,'unsupported kind explicit status')
 poly_kind(1)=1
 if(proc0)then
 allocate(f0(nspec,0:nperp,0:npar))
 open(newunit=unit,file='rows.csv',status='new')
 write(unit,'(a)')'phase,policy,species,row,status,state,requested,effective,original,retained,rank,sample_limited,rank_reduced,lower_zeros,upper_zeros,original_min,original_max,retained_min,retained_max,log_L2,svd_tolerance'
 open(newunit=query_unit,file='queries.csv',status='new')
 write(query_unit,'(a)')'phase,species,row,status,value_re,value_im,derivative_re,derivative_im,expected_re,expected_im,expected_derivative_re,expected_derivative_im'
 endif
 nan=ieee_value(0.d0,ieee_quiet_nan)
 do phase=0,8
 policy=NR_FIT_NATIVE
 if(phase==0.or.phase==2)policy=NR_FIT_LEGACY_COMPAT
 poly_order=[4,1];logfit=.true.
 if(phase>=4.and.phase<=6)logfit=.false.
 if(phase==5)poly_order=0
 if(phase==6)poly_order=1
 if(phase==7)poly_order=0
 if(phase==8)then
 logfit=.false.;useBM(2)=.true.
 endif
 if(proc0)then
 do s=1,nspec
 b=.02d0;if(s==2)b=0.d0
 do r=0,nperp
 do j=0,npar
 u=1.d0+3.d0*j/npar
 if(r==7.and.phase>=1.and.phase<=3)u=2.d0
 if(r==0.and.phase==3)then
 if(j==1)u=2.2d0
 if(j==2)u=1.3d0
 endif
 pp(s,r,j,2)=u;pp(s,r,j,1)=.1d0*r
 f0(s,r,j)=10.d0**(-2.d0+.3d0*u-b*u*u)
 if(.not.logfit(s))then
 f0(s,r,j)=2.d0+.3d0*u
 if(poly_order(s)>=4)f0(s,r,j)=f0(s,r,j)+.04d0*u*u+.01d0*u**4
 if(phase==5)f0(s,r,j)=2.d0
 endif
 if(phase==7)f0(s,r,j)=1.d-310
 enddo
 if(phase>=1.and.phase<=3)then
 if(r==0)f0(s,r,0)=0.d0
 if(r==1)then
 f0(s,r,npar)=0.d0
 endif
 if(r==2)then
 f0(s,r,0:1)=0.d0;f0(s,r,npar-1:npar)=0.d0
 endif
 if(r==3)f0(s,r,3)=0.d0
 if(r==4)f0(s,r,:)=0.d0
 if(r==5)f0(s,r,2)=-.1d0
 if(r==6)f0(s,r,2)=nan
 if(r==6.and.phase==3)then
 f0(s,r,2)=1.d0;pp(s,r,2,2)=nan
 endif
 if(r==8)then
 f0(s,r,:)=0.d0;u=pp(s,r,3,2)
 f0(s,r,3)=10.d0**(-2.d0+.3d0*u-b*u*u)
 endif
 endif
 if(phase==4.and.r==4)f0(s,r,:)=0.d0
 enddo
 enddo
 if(phase==0)then
 call legacy_determine_param_fit
 allocate(legacy(nspec,0:nperp,0:4));legacy=reference_coefficients
 endif
 oldsum=sum(f0,mask=ieee_is_finite(f0))
 endif
 call initialize_nr_cheb_svd(policy,code)
 if(proc0)call check(sum(f0,mask=ieee_is_finite(f0))==oldsum,'builder preserves VDF values')
 call check(code==CHEB_SUCCESS,'root setup / worker no-op')
 call nr_cheb_collection_info(collection_species,rows,builds,fits,evals,setup)
 if(.not.proc0)call check(fits==0.and.builds==0,'workers never fit')
 if(proc0.and.phase==8)call check(fits==nperp+1,'bimaxwellian species excluded from fitting')
 call pass_nr_cheb_models(code);call check(code==CHEB_SUCCESS,'single initialization transfer')
 call nr_cheb_distribution_info(tx,collect);call check(tx==1.and.collect==6,'fixed collectives')
 call set_nr_cheb_evaluation(.true.)
 call nr_cheb_collection_info(collection_species,rows,oldbuilds,oldfits,evals,setup)
 oldtx=tx;oldcollect=collect
 call nr_cheb_export(meta,mapping,coef,code)
 n=size(meta)+size(mapping)+2*size(coef)
 allocate(bits(n),small(n),large(n));at=0
 bits(1:size(meta))=int(reshape(meta,[size(meta)]),int64);at=size(meta)
 bits(at+1:at+size(mapping))=transfer(mapping,[0_int64],size(mapping));at=at+size(mapping)
 bits(at+1:at+size(coef))=transfer(coef,[0_int64],size(coef));at=at+size(coef)
 do r=0,nperp
 do s=1,nspec
 call get_nr_cheb_row(s,r,model,info,code)
 expected_code=CHEB_SUCCESS
 if(phase==8.and.s==2)expected_code=CHEB_INVALID_MODEL
 if(phase==1.or.phase==3)then
 if(r==3)expected_code=CHEB_INTERIOR_ZERO
 if(r==4)expected_code=CHEB_NO_POSITIVE_SAMPLES
 if(r==5)expected_code=CHEB_NEGATIVE_SAMPLES
 if(r==6)expected_code=CHEB_INVALID_VALUES
 if(phase==3.and.r==6)expected_code=CHEB_INVALID_COORDINATES
 endif
 if(phase==2)then
 if(r==5)expected_code=CHEB_NEGATIVE_SAMPLES
 if(r==6)expected_code=CHEB_INVALID_VALUES
 if(r==7)expected_code=CHEB_INVALID_COORDINATES
 endif
 call check(code==expected_code.or.(expected_code==CHEB_SUCCESS.and.code==CHEB_REDUCED),'row classification')
 if(proc0)write(unit,'(*(g0,:,","))')phase,info%policy,s,r,code,info%state,model%requested_degree,model%effective_degree, &
 model%original_sample_count,model%retained_sample_count,model%rank,model%sample_limited,model%rank_reduced, &
 model%lower_tail_zeros,model%upper_tail_zeros,info%original_min,info%original_max,model%retained_min,model%retained_max, &
 model%residual_norm,model%svd_tolerance
 bits(at+1:at+size(coef,1))=0_int64
 if(allocated(model%derivative_coefficients))then
 bits(at+1:at+size(model%derivative_coefficients))=transfer(model%derivative_coefficients,[0_int64],size(model%derivative_coefficients))
 endif
 at=at+size(coef,1)
 z=(2.4d0,.05d0);want=cmplx(nan,nan,kind(1.d0));wantd=want
 if(code==CHEB_SUCCESS.or.code==CHEB_REDUCED)then
 call check(model%rank==model%effective_degree+1,'accepted rank')
 if(phase==0.and.proc0)call check(maxval(abs(model%coefficients-legacy(s,r,0:poly_order(s))))<1.d-11,'matched GLLS/SVD coefficients')
 if((phase==1.or.phase==3).and.r==2.and.s==1)call check(model%effective_degree==2.and.model%sample_limited,'sparse degree reduction')
 if((phase==1.or.phase==3).and.r==8)call check(model%effective_degree==0.and.model%retained_sample_count==1,'singleton')
 if(phase==1.and.r==7)call check(model%effective_degree==0.and.model%rank_reduced,'rank refit')
 if(phase==1.and.r==1)call check(model%retained_min==1.d0.and.model%retained_max==3.5d0,'trimmed actual interval')
 do j=1,10
 v=eval_fit(s,r,z,dfpar=d,status=code)
 call check(code==CHEB_SUCCESS,'local successful evaluation')
 enddo
 b=.02d0;if(s==2)b=0.d0
 if(logfit(s))then
 want=10.d0**(-2.d0+.3d0*z-b*z*z);wantd=log(10.d0)*want*(.3d0-2.d0*b*z)
 if((phase==1.or.phase==3).and.(r==7.or.r==8))then
 u=2.d0;if(r==8)u=2.5d0
 want=cmplx(10.d0**(-2.d0+.3d0*u-b*u*u),0.d0,kind(1.d0));wantd=(0.d0,0.d0)
 endif
 else
 want=2.d0+.3d0*z;wantd=(.3d0,0.d0)
 if(poly_order(s)>=4)then
 want=want+.04d0*z*z+.01d0*z**4;wantd=wantd+.08d0*z+.04d0*z**3
 endif
 if(phase==5)then
 want=(2.d0,0.d0);wantd=(0.d0,0.d0)
 endif
 if(phase==4.and.r==4)then
 want=(0.d0,0.d0);wantd=(0.d0,0.d0)
 endif
 endif
 if(phase==7)then
 want=cmplx(1.d-310,0.d0,kind(1.d0));wantd=(0.d0,0.d0)
 call check(abs(v/1.d-310-1.d0)<1.d-10,'positive subnormal retained/reconstructed')
 elseif(phase/=2)then
 call check(abs(v-want)<1.d-10*max(abs(want),1.d-12),'independent complex reconstruction')
 call check(abs(d-wantd)<1.d-10*max(abs(want),1.d-12),'independent physical derivative')
 endif
 else
 call eval_nr_cheb(s,r,z,v,code,d)
 call check(code==expected_code.and.ieee_is_nan(real(v)).and.ieee_is_nan(real(d)),'unavailable is not zero')
 endif
 if(phase==2)then
 want=cmplx(nan,nan,kind(1.d0));wantd=want ! no independent floor-fit oracle
 endif
 if(proc0)write(query_unit,'(*(g0,:,","))')phase,s,r,code,real(v),aimag(v),real(d),aimag(d),real(want),aimag(want),real(wantd),aimag(wantd)
 enddo
 enddo
 call mpi_allreduce(bits,small,n,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierror)
 call mpi_allreduce(bits,large,n,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierror)
 call check(all(small==large),'bit-identical metadata/maps/value and derivative coefficients')
 deallocate(bits,small,large,meta,mapping,coef)
 call nr_cheb_collection_info(collection_species,rows,builds,fits,evals,setup)
 call nr_cheb_distribution_info(tx,collect)
 call check(builds==oldbuilds.and.fits==oldfits.and.tx==oldtx.and.collect==oldcollect,'no fitting or transfer during queries')
 enddo
 ! Mixed numerical methods: shared activation for species 1 must not route
 ! a functional species through its unavailable/unselected NR model slot.
 useBM(2)=.false.;ACmethod(2)=1
 fit_type(2,1)=1;perp_correction(2,1)=.75d0
 param_fit(2,:,1,1)=1.2d0;param_fit(2,:,2,1)=.8d0;param_fit(2,:,3,1)=.1d0
 z=(.2d0,-.1d0)
 want=1.2d0*exp(-.75d0*pp(2,1,1,1)**2)*exp(-.8d0*(z-.1d0)**2)
 call check(nr_cheb_evaluation_active(1).and..not.nr_cheb_evaluation_active(2),'mixed numerical eligibility')
 call nr_cheb_collection_info(collection_species,rows,builds,fits,evals,setup);n=evals
 v=eval_fit(2,1,z,status=code)
 call check(code==CHEB_SUCCESS.and.abs(v-want)<1.d-13,'functional complex evaluation with active NR species')
 call nr_cheb_collection_info(collection_species,rows,builds,fits,evals,setup)
 call check(evals==n,'functional species bypasses NR model evaluation')
 ACmethod(2)=2;useBM(2)=.true.
 call initialize_nr_cheb_svd(99,code);call pass_nr_cheb_models(code)
 call check(code==CHEB_INVALID_VALUES,'failed root setup reaches all ranks')
 if(proc0)then
 close(unit);close(query_unit)
 endif
 call mpi_allreduce(failures,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
 write(*,'(a,i0,a,i0,a,i0)')'NR SVD rank=',iproc,' passed=',checks,' failed=',failures
 if(total/=0)call mpi_abort(MPI_COMM_WORLD,99,ierror)
 call mpi_finalize(ierror)
contains
 subroutine check(ok,label)
 logical,intent(in)::ok
 character(*),intent(in)::label
 if(ok)then
 checks=checks+1
 else
 failures=failures+1;print *,'FAIL rank',iproc,': ',label
 endif
 end subroutine
end program
