program test_nr_zero_return
 use mpi
 use alps_var
 use alps_io,only:init_param,alps_error_init
 use alps_com,only:pass_instructions,pass_nr_cheb_models,nr_cheb_distribution_info
 use alps_analyt
 use alps_fns,only:landau_integrate,landau_integrate_before,landau_integrate_candidate,determine_bessel_array
 use nr_zero_diagnostics,only:context,b5_begin,b5_counts,b5_recording
 use, intrinsic::ieee_arithmetic
 use legacy_glls_reference
 implicit none
 integer::s,r,j,c,selector,code,unit,reference_unit,samples,builds,fits,evals,setup,nb,nf,tx,nc,checks=0,failed=0,total
 integer::previous_tx,previous_nc,rows
 double precision::u,step,h,weight,nan,a
 double complex::omega,z,v,d,old,new,actual,want,acc,local,fu,fp,fm,gp,gm,fg,expected_old,previous
 character(20)::stage,option
 logical::shared,stopped
 call mpi_init(ierror)
 call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror);call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
 proc0=iproc==0
 if(proc0)then
 call alps_error_init;call init_param
 endif
 call get_command_argument(2,stage);call get_command_argument(3,option)
 call pass_instructions
 call ensure_historical_storage
 logfit=.false.;poly_order=2;sproc=1;kpar=.125d0;kperp=0.d0;nlim=0;ms=1.d0;qs=1.d0
 do s=1,nspec
 do r=0,nperp
 do j=0,npar
 pp(s,r,j,1)=.1d0*r;pp(s,r,j,2)=-2.d0+4.d0*j/npar
 enddo
 enddo
 enddo
 call determine_bessel_array
 call check(.not.nr_cheb_evaluation_active(1),'legacy remains default')
 call select_nr_test_evaluation(.true.)
 call check(nr_cheb_evaluation_active(1),'shared method two selected')
 ACmethod(1)=1;call check(.not.nr_cheb_evaluation_active(1),'functional path retains legacy policy');ACmethod(1)=2
 relativistic(1)=.true.;call check(.not.nr_cheb_evaluation_active(1),'relativistic path not selected');relativistic(1)=.false.
 call select_nr_test_evaluation(.false.)
 if(proc0)then
 allocate(f0(nspec,0:nperp,0:npar))
 open(newunit=unit,file='contours.csv',status='new')
 open(newunit=reference_unit,file='stencil_reference.csv',status='new')
 write(reference_unit,'(a)')'case,selector,row,status,F_re,F_im,kernel_analytic_Fu_re,kernel_analytic_Fu_im,Fplus_re,Fplus_im,Fminus_re,Fminus_im,Fperpplus_re,Fperpplus_im,Fperpminus_re,Fperpminus_im,selected_reference_Fu_re,selected_reference_Fu_im,reference_Fperp_re,reference_Fperp_im,local_re,local_im,accumulator_re,accumulator_im'
 write(unit,'(a)')'case,selector,legacy_re,legacy_im,old_shared_re,old_shared_im,candidate_re,candidate_im,actual_shared_re,actual_shared_im,reference_re,reference_im'
 endif
 nan=ieee_value(0.d0,ieee_quiet_nan)
 omega=cmplx(0.d0,-.015625d0,kind(1.d0));z=omega/kpar;step=4.d0/npar
 do c=0,8
 kpar=.125d0;omega=cmplx(0.d0,-.015625d0,kind(1.d0));z=omega/kpar
 logfit=c>=7
 historical_coefficients=0.d0
 do s=1,nspec
 do r=0,nperp
 historical_coefficients(s,r,0)=2.d0;historical_coefficients(s,r,1)=.3d0
 if((c==1.and.r==3).or.(c==2.and.r==2).or.(c==6.and.(r==0.or.r==nperp)))historical_coefficients(s,r,:)=0.d0
 if(c==3.and.r==3)historical_coefficients(s,r,:)=historical_coefficients(s,r,:)*1.d-14
 if(c==4)then
 historical_coefficients(s,r,:)=[2.265625d0,-2.d0,2.d0]
 endif
 if(c==5.and.r==2)historical_coefficients(s,r,:)=[2.015625d0,0.d0,2.d0]
 if(c>=7)historical_coefficients(s,r,:)=[-2.d0,.3d0,0.d0]
 enddo
 enddo
 if(c<8)then
 call initialize_historical_nr(code)
 else
 if(proc0)then
 do s=1,nspec
 do r=0,nperp
 do j=0,npar
 u=pp(s,r,j,2);f0(s,r,j)=0.d0
 if(abs(u)<=.5d0)f0(s,r,j)=(2.d0+.01d0*r)*exp(-u*u)
 enddo
 enddo
 enddo
 endif
 call initialize_nr_cheb_svd(NR_FIT_NATIVE,code)
 call check(code==CHEB_SUCCESS,'native trimmed row setup')
 call pass_nr_cheb_models(code)
 endif
 call check(code==CHEB_SUCCESS,'models available')
 call select_nr_test_evaluation(.true.)
 if(c==4)then
 call eval_nr_cheb(1,1,z+step,v,code,d)
 call check(code==CHEB_SUCCESS.and.v==(0.d0,0.d0).and.abs(d)>0,'valid zero parallel neighbor and finite derivative')
 endif
 if(c==5)then
 call eval_nr_cheb(1,2,z,v,code,d)
 call check(code==CHEB_SUCCESS.and.v==(0.d0,0.d0).and.d==2*z,'valid zero current row and finite derivative')
 endif
 call nr_cheb_collection_info(samples,rows,nb,nf,evals,setup)
 call nr_cheb_distribution_info(previous_tx,previous_nc)
 do selector=0,1
 poly_derivative_type=selector
 kpar=.125d0;omega=cmplx(0.d0,-.015625d0,kind(1.d0));z=omega/kpar
 if(c==1)then
 kpar=.1d0;omega=(0.d0,-.01d0);z=omega/kpar
 endif
 if(c==7)then
 omega=cmplx(0.d0,-.2525d0,kind(1.d0));z=omega/kpar
 endif
 call context(c,selector)
 ! independent oracle: n=0, kperp=0 => J0=1, Tzz=z**2.
 acc=(0.d0,0.d0);expected_old=acc;stopped=.false.
 do j=1,nperp+1
 r=j;if(j==nperp)r=0
 if(j==nperp+1)r=nperp
 h=1.d0;if(r==0.or.r>=nperp-1)h=.5d0
 if(selector==1)then
 fu=derivative(r,z,c)
 else
 fu=(shape(r,z+step,c)-shape(r,z-step,c))/(2*step)
 endif
 if(r==0)then
 fg=(shape(1,z,c)-shape(0,z,c))/.1d0
 elseif(r==nperp)then
 fg=(shape(r,z,c)-shape(r-1,z,c))/.1d0
 else
 fp=shape(r,z+step,c);fm=shape(r,z-step,c)
 gp=shape(r+1,z,c);gm=shape(r-1,z,c)
 if(.not.stopped.and.(abs(fp)==0.or.abs(fm)==0.or.abs(gp)==0.or.abs(gm)==0))then
 stopped=.true.;expected_old=(0.d0,0.d0)
 endif
 fg=(gp-gm)/.2d0
 endif
 call eval_nr_cheb(1,r,z,v,code,d)
 call check(code==CHEB_SUCCESS,'finite current value and analytic derivative')
 call check(abs(d-derivative(r,z,c))<2.d-12*max(1.d0,abs(d)),'independent analytic derivative at possible zero')
 local=-h*z*z*(1.d0/abs(kpar))*((.1d0*r*fu-z*fg)*kpar+omega*fg)
 acc=acc+local
 if(proc0)write(reference_unit,'(*(g0,:,","))')c,selector,r,code,real(v),aimag(v),real(d),aimag(d), &
 real(shape(r,z+step,c)),aimag(shape(r,z+step,c)),real(shape(r,z-step,c)),aimag(shape(r,z-step,c)), &
 real(shape(min(nperp,r+1),z,c)),aimag(shape(min(nperp,r+1),z,c)), &
 real(shape(max(0,r-1),z,c)),aimag(shape(max(0,r-1),z,c)),real(fu),aimag(fu),real(fg),aimag(fg),real(local),aimag(local),real(acc),aimag(acc)
 if(.not.stopped)expected_old=expected_old+local
 enddo
 want=acc*cmplx(0.d0,1.d0,kind(1.d0))*.1d0*2*pi*pi
 expected_old=expected_old*cmplx(0.d0,1.d0,kind(1.d0))*.1d0*2*pi*pi
 old=landau_integrate_before(omega,0,3)
 new=landau_integrate_candidate(omega,0,3)
 actual=landau_integrate(omega,0,3)
 call check(abs(new-old)<1.d-12.or.c/=0,'positive control no change')
 call check(abs(new-want)<2.d-12*max(1.d0,abs(want)),'candidate matches independent numerator/quadrature')
 call check(abs(old-expected_old)<2.d-12*max(1.d0,abs(expected_old)),'old gate reference')
 if(c==7)call check(abs(old-new)<1.d-12.and.abs(new)>0,'suppression removed by evaluator, independent of gate')
 if(stage=='pre')then
 call check(actual==old,'pre-change production equals old gate')
 else
 call check(actual==new,'corrected production equals candidate')
 endif
 if(c/=8)then
 call select_nr_test_evaluation(.false.)
 previous=landau_integrate(omega,0,3)
 v=landau_integrate_before(omega,0,3)
 call check(previous==v,'legacy retains old gate')
 if(c==7)call check(previous==(0.d0,0.d0).and.abs(new)>0,'legacy disk suppression distinguished')
 call select_nr_test_evaluation(.true.)
 else
 previous=cmplx(nan,nan,kind(1.d0)) ! no legacy fit for native-only control
 endif
 if(proc0)write(unit,'(*(g0,:,","))')c,selector,real(previous),aimag(previous),real(old),aimag(old),real(new),aimag(new),real(actual),aimag(actual),real(want),aimag(want)
 omega=cmplx(0.d0,-.015625d0,kind(1.d0));z=omega/kpar
 enddo
 call nr_cheb_collection_info(samples,rows,builds,fits,evals,setup)
 call nr_cheb_distribution_info(tx,nc)
 call check(builds==nb.and.fits==nf.and.tx==previous_tx.and.nc==previous_nc,'no fitting/communication during integration')
 enddo
 call method_zero_policy_checks
 call evaluation_dependency_checks
 if(proc0)then
 close(unit);close(reference_unit)
 endif
 ! Process-isolated failure checks exercise status-less eval_fit via the integrator.
 if(option/='')then
 logfit=.true.;historical_coefficients=0.d0;historical_coefficients(:,:,0)=-2.d0
 if(option=='underflow')historical_coefficients(:,2,0)=-400.d0
 if(option=='overflow')historical_coefficients(:,2,0)=400.d0
 if(option=='derivative')then
 deallocate(historical_coefficients);allocate(historical_coefficients(nspec,0:nperp,0:3));poly_order=3
 historical_coefficients=0.d0;historical_coefficients(:,:,0)=2.d0
 a=scale(1.d0,1019);historical_coefficients(:,1,1)=5.5d0*a;historical_coefficients(:,1,3)=2.d0*a
 omega=(0.d0,0.d0);poly_derivative_type=1
 endif
 call initialize_historical_nr(code)
 call check(code==CHEB_SUCCESS,'failure fixture imports finite coefficients')
 if(option=='unavailable'.or.option=='invalid')then
 if(proc0)then
 f0=1.d0
 if(option=='unavailable')f0(:,2,:)=0.d0
 if(option=='invalid')f0(:,2,3)=-1.d0
 endif
 call initialize_nr_cheb_svd(NR_FIT_NATIVE,code);call pass_nr_cheb_models(code)
 endif
 if(option=='missing')call release_nr_cheb
 call select_nr_test_evaluation(.true.)
 if(option=="inactive")call set_nr_cheb_evaluation(.false.)
 call context(99,0)
 v=landau_integrate(omega,0,3)
 call check(.false.,'failed evaluation must abort before use')
 endif
 call mpi_allreduce(failed,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
 write(*,'(a,i0,a,i0,a,i0)')'Zero return rank=',iproc,' passed=',checks,' failed=',failed
 if(total/=0)call mpi_abort(MPI_COMM_WORLD,99,ierror)
 call mpi_finalize(ierror)
contains
 subroutine evaluation_dependency_checks
 integer::degree,selected,row,idx,category,edge,request,file,n0(3,2,2),new_counts(3,2,2)
 integer::saved_builds,saved_fits,saved_tx,saved_collect,old_evaluations
 double complex::got,prior,oracle,point,du,dp
 double precision::half
 character(80)::filename
 write(filename,'(a,i0,a)')'query_counts_rank',iproc,'.csv'
 open(newunit=file,file=trim(filename),status='new')
 write(file,'(a)')'degree,selector,variant,kind,endpoint,derivative,count,contour_re,contour_im'
 point=cmplx(0.d0,-.125d0,kind(1.d0));omega=point*kpar
 logfit=.false.
 do degree=0,4
 if(degree==2.or.degree==3)cycle
 poly_order=degree
 deallocate(historical_coefficients);allocate(historical_coefficients(nspec,0:nperp,0:degree))
 historical_coefficients=0.d0
 do row=0,nperp
 historical_coefficients(:,row,0)=2.d0+.2d0*row
 if(degree>=1)historical_coefficients(:,row,1)=.6d0
 if(degree==4)then
 historical_coefficients(:,row,0)=1.92d0+.2d0*row
 historical_coefficients(:,row,2)=-.04d0;historical_coefficients(:,row,4)=.04d0
 endif
 enddo
 call initialize_historical_nr(code);call check(code==CHEB_SUCCESS,'B5 known-coefficient model import')
 call select_nr_test_evaluation(.true.)
 do row=0,nperp
 call eval_nr_cheb(1,row,point,v,code,d)
 du=(0.d0,0.d0)
 if(degree>=1)du=.3d0
 if(degree==4)du=du-.2d0*point+.08d0*point**3
 call check(code==CHEB_SUCCESS.and.abs(v-b5_polynomial(row,point,degree))<2.d-13,'B5 independent physical distribution at complex resonance')
 call check(abs(d-du)<2.d-13,'B5 independent physical parallel derivative at fixed row')
 enddo
 do selected=0,1
 poly_derivative_type=selected
 call nr_cheb_collection_info(samples,rows,saved_builds,saved_fits,old_evaluations,setup)
 call nr_cheb_distribution_info(saved_tx,saved_collect)
 call b5_begin(degree,selected,'original',.false.)
 prior=landau_integrate_candidate(omega,0,3);n0=b5_counts
 call b5_begin(degree,selected,'optimized',.false.)
 got=landau_integrate(omega,0,3);new_counts=b5_counts
 b5_recording=.false.
 call check(got==prior,'B5 retained contour arithmetic byte-identical')
 oracle=(0.d0,0.d0)
 do idx=1,nperp+1
 row=idx
 if(idx==nperp)row=0
 if(idx==nperp+1)row=nperp
 half=1.d0
 if(row==0.or.row==nperp.or.row==nperp-1)half=.5d0
 du=(0.d0,0.d0)
 if(degree>=1)du=.3d0
 if(degree==4)du=du-.2d0*point+.08d0*point**3
 if(selected==0)du=(b5_polynomial(row,point+step,degree)-b5_polynomial(row,point-step,degree))/(2*step)
 dp=2.d0
 oracle=oracle-half*point**2/abs(kpar)*((.1d0*row*du-point*dp)*kpar+omega*dp)
 enddo
 oracle=oracle*cmplx(0.d0,1.d0,kind(1.d0))*.1d0*2*pi*pi
 call check(abs(got-oracle)<2.d-12,'B5 independent polynomial/discrete contour oracle')
 call check(sum(new_counts(2,:,:))==2*nperp,'B5 required perpendicular evaluations retained')
 if(selected==0)then
 call check(all(n0==new_counts),'B5 FD exact call sequence counts unchanged')
 call check(sum(new_counts(3,:,:))==2*(nperp+1),'B5 FD shifted parallel evaluations retained')
 call check(sum(new_counts(:,:,2))==0,'B5 FD does not request analytic derivative')
 else
 call check(sum(new_counts(3,:,:))==0,'B5 analytic no shifted parallel queries')
 call check(sum(n0(3,:,:))==2*(nperp-1),'B5 original redundant interior parallel query count')
 call check(sum(new_counts(:,:,2))==nperp+1,'B5 analytic derivative requested at every row')
 call check(new_counts(1,2,1)==2.and.new_counts(1,2,2)==2,'B5 endpoint center values and derivatives retained')
 call check(sum(n0)-sum(new_counts)==2*(nperp-1),'B5 exact removed query count')
 endif
 call nr_cheb_collection_info(samples,rows,builds,fits,evals,setup);call nr_cheb_distribution_info(tx,nc)
 call check(builds==saved_builds.and.fits==saved_fits.and.tx==saved_tx.and.nc==saved_collect,'B5 queries never refit or transfer models')
 call check(evals-old_evaluations==sum(n0)+sum(new_counts),'B5 query trace matches kernel evaluation counter')
 do request=1,2
 do edge=1,2
 do category=1,3
 write(file,'(*(g0,:,","))')degree,selected,'original',category,edge-1,request==2,n0(category,edge,request),real(prior),aimag(prior)
 write(file,'(*(g0,:,","))')degree,selected,'optimized',category,edge-1,request==2,new_counts(category,edge,request),real(got),aimag(got)
 enddo
 enddo
 enddo
 enddo
 enddo
 ! All required points of this exact quartic succeed. The deliberately injected
 ! status-18 query is row 1 at u=0.5-0.125i, needed only by FD interior calls.
 poly_derivative_type=1
 call b5_begin(4,1,'unused_analytic',.true.)
 got=landau_integrate(omega,0,3);b5_recording=.false.
 call check(got==prior,'B5 unused-query failure does not affect analytic contour')
 if(option=='unused_fd'.or.option=='unused_old')then
 selected=0
 if(option=='unused_old')selected=1
 poly_derivative_type=selected
 call b5_begin(4,selected,trim(option),.true.)
 if(option=='unused_fd')got=landau_integrate(omega,0,3)
 if(option=='unused_old')got=landau_integrate_candidate(omega,0,3)
 call check(.false.,'B5 required or original-unused shifted failure must abort')
 endif
 close(file)
 poly_order=2;poly_derivative_type=1
 deallocate(historical_coefficients);allocate(historical_coefficients(nspec,0:nperp,0:2));historical_coefficients=0.d0
 end subroutine
 function b5_polynomial(row,u,degree)result(value)
 integer,intent(in)::row,degree
 double complex,intent(in)::u
 double complex::value
 value=2.d0+.2d0*row
 if(degree>=1)value=value+.3d0*u
 if(degree==4)value=value-.1d0*u*u+.02d0*u**4
 end function
 subroutine method_zero_policy_checks
 integer::method,zero_case,row,idx,first_zero,file,stat,old_counts(3,2,2)
 double precision::delta,half,parallel_step,coordinate
 double complex::point,plus,minus,upper,lower,du,dp,term,sum_all,sum_gate,prior,oracle,old_value,candidate,actual_value
 logical::gate
 character(80)::filename
 write(filename,'(a,i0,a)')'method_policy_rank',iproc,'.csv'
 open(newunit=file,file=trim(filename),status='new')
 write(file,'(a)')'method,zero_case,first_zero,prior_re,prior_im,oracle_re,oracle_im,before_re,before_im,candidate_re,candidate_im,actual_re,actual_im'
 poly_derivative_type=0;kpar=.125d0;omega=(0.d0,-.015625d0);point=omega/kpar
 parallel_step=abs(pp(1,2,2,2)-pp(1,2,1,2))
 do method=0,1
 do zero_case=0,1
 ACmethod=method
 delta=.1d0
 do s=1,nspec
 fit_type(s,1)=1;perp_correction(s,1)=0.d0
 do row=0,nperp
 pp(s,row,:,1)=delta*row
 param_fit(s,row,1,1)=2.d0+.2d0*row
 if(method==1.and.zero_case==1.and.row==3)param_fit(s,row,1,1)=0.d0
 param_fit(s,row,2,1)=.7d0;param_fit(s,row,3,1)=.1d0
 enddo
 enddo
 ! Test-only endpoint induces a floating-point zero after resolved inner rows.
 ! The oracle retains the integrator's existing fixed-delta stencil.
 if(method==0.and.zero_case==1)pp(:,nperp,:,1)=32.d0
 sum_all=(0.d0,0.d0);sum_gate=sum_all;prior=sum_all;gate=.false.;first_zero=-1
 do idx=1,nperp+1
 row=idx
 if(idx==nperp)row=0
 if(idx==nperp+1)row=nperp
 half=1.d0
 if(row==0.or.row==nperp.or.row==nperp-1)half=.5d0
 plus=method_shape(row,point+parallel_step,method,zero_case,delta)
 minus=method_shape(row,point-parallel_step,method,zero_case,delta)
 du=(plus-minus)/(2*parallel_step)
 if(row==0)then
 upper=method_shape(1,point,method,zero_case,delta)
 lower=method_shape(0,point,method,zero_case,delta)
 dp=(upper-lower)/delta
 elseif(row==nperp)then
 upper=method_shape(row,point,method,zero_case,delta)
 lower=method_shape(row-1,point,method,zero_case,delta)
 dp=(upper-lower)/delta
 else
 upper=method_shape(row+1,point,method,zero_case,delta)
 lower=method_shape(row-1,point,method,zero_case,delta)
 dp=(upper-lower)/(2*delta)
 if(.not.gate.and.(abs(plus)==0.or.abs(minus)==0.or.abs(upper)==0.or.abs(lower)==0))then
 first_zero=row;prior=sum_gate;sum_gate=(0.d0,0.d0);gate=.true.
 endif
 endif
 coordinate=pp(1,row,0,1)
 term=-half*point*point/abs(kpar)*((coordinate*du-point*dp)*kpar+omega*dp)
 sum_all=sum_all+term
 if(.not.gate)sum_gate=sum_gate+term
 enddo
 oracle=sum_all*cmplx(0.d0,1.d0,kind(1.d0))*delta*2*pi*pi
 sum_gate=sum_gate*cmplx(0.d0,1.d0,kind(1.d0))*delta*2*pi*pi
 call context(100+method,0)
 call b5_begin(-1,0,"method"//achar(48+method)//"_old"//achar(48+zero_case),.false.)
 old_value=landau_integrate_before(omega,0,3);old_counts=b5_counts;b5_recording=.false.
 candidate=landau_integrate_candidate(omega,0,3)
 call b5_begin(-1,0,"method"//achar(48+method)//"_now"//achar(48+zero_case),.false.)
 actual_value=landau_integrate(omega,0,3);b5_recording=.false.
 call check(all(b5_counts==old_counts),"B5 method-0/1 evaluation counts unchanged")
 call check(abs(candidate-oracle)<=1.d-10*max(abs(oracle),tiny(1.d0)),'method-0/1 independent discrete numerator oracle')
 call check(actual_value==old_value,'method-0/1 actual policy unchanged')
 if(zero_case==1)then
 call check(first_zero>1.and.abs(prior)>0,'method-0/1 zero after nonzero accumulation')
 call check(actual_value==(0.d0,0.d0).and.abs(candidate)>0,'method-0/1 zero retains whole-contour return')
 else
 call check(first_zero==-1.and.actual_value==candidate,'method-0/1 positive control identical')
 endif
 write(file,'(*(g0,:,","))')method,zero_case,first_zero,real(prior),aimag(prior),real(oracle),aimag(oracle), &
 real(old_value),aimag(old_value),real(candidate),aimag(candidate),real(actual_value),aimag(actual_value)
 enddo
 enddo
 close(file)
 ACmethod=2;poly_derivative_type=1
 do s=1,nspec
 do row=0,nperp
 pp(s,row,:,1)=.1d0*row
 enddo
 enddo
 end subroutine
 function method_shape(row,u,method,zero_case,delta)result(value)
 integer,intent(in)::row,method,zero_case
 double precision,intent(in)::delta
 double complex,intent(in)::u
 double complex::value
 double precision::amplitude
 if(method==0)then
 ! Independent expression for the repository's species-1 analytic Maxwellian.
 value=pi**(-1.5d0)*exp(-(pp(1,row,0,1)**2+u*u))
 else
 amplitude=2.d0+.2d0*row
 if(zero_case==1.and.row==3)amplitude=0.d0
 value=amplitude*exp(-.7d0*(u-.1d0)**2)
 endif
 end function

 function shape(row,u,which)result(f)
 integer,intent(in)::row,which
 double complex,intent(in)::u
 double complex::f
 f=2.d0+.15d0*u
 if((which==1.and.row==3).or.(which==2.and.row==2).or.(which==6.and.(row==0.or.row==nperp)))f=0.d0
 if(which==3.and.row==3)f=f*1.d-14
 if(which==4)f=(u-.5d0)**2+.015625d0
 if(which==5.and.row==2)f=u*u+.015625d0
 if(which==7)f=10.d0**(-2.d0+.15d0*u)
 if(which==8)f=(2.d0+.01d0*row)*exp(-u*u)
 end function
 function derivative(row,u,which)result(f)
 integer,intent(in)::row,which
 double complex,intent(in)::u
 double complex::f
 f=.15d0
 if((which==1.and.row==3).or.(which==2.and.row==2).or.(which==6.and.(row==0.or.row==nperp)))f=0.d0
 if(which==3.and.row==3)f=f*1.d-14
 if(which==4)f=2*(u-.5d0)
 if(which==5.and.row==2)f=2*u
 if(which==7)f=log(10.d0)*.15d0*shape(row,u,which)
 if(which==8)f=-2*u*shape(row,u,which)
 end function
 subroutine check(ok,label)
 logical,intent(in)::ok
 character(*),intent(in)::label
 if(ok)then
 checks=checks+1
 else
 failed=failed+1;print *,'FAIL rank=',iproc,trim(label),c,selector
 endif
 end subroutine
end program
