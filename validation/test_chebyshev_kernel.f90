! Independent direct-polynomial tests: no production fitting/evaluation wrapper.
program test_chebyshev_kernel
  use alps_analyt, only : cheb_model,cheb_fit,cheb_evaluate,CHEB_SUCCESS,CHEB_REDUCED, &
       CHEB_NO_SAMPLES,CHEB_INVALID_COORDINATES,CHEB_INVALID_NORMALIZATION,CHEB_INVALID_DEGREE, &
       CHEB_INVALID_VALUES,CHEB_INVALID_TOLERANCE,CHEB_INVALID_MODEL,CHEB_EVALUATION_FAILURE, &
       CHEB_NONFINITE_COEFFICIENTS,CHEB_DIAG_VALID,CHEB_DIAG_RANGE,CHEB_DIAG_FAILED
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer, parameter :: degrees(4)=[3,6,12,30]
  type(cheb_model) :: model,reference,constant_model,linear_model
  integer :: status,i,j,n,passed=0,failed=0,degree
  real(dp) :: u(9),f(9),value,derivative,step,pi,nan,inf,expected,rnorm
  real(dp) :: nodes(41),samples(41),cluster(6),target(6),loose_residual
  real(dp) :: hmax,a,b,gap,multiply_a,multiply_b,nearby,smallest
  logical :: halting_before(3),halting_after(3)
  complex(dp) :: z,result,dresult,want,dwant,plus,minus

  pi=4.d0*atan(1.d0)
  call ieee_get_halting_mode(ieee_invalid,halting_before(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,halting_before(2))
  call ieee_get_halting_mode(ieee_overflow,halting_before(3))
  nan=ieee_value(0.d0,ieee_quiet_nan); inf=ieee_value(0.d0,ieee_positive_inf)
  u=[-2.d0,-1.3d0,-0.7d0,-0.1d0,0.d0,0.4d0,1.1d0,1.7d0,3.d0]
  f=3.25d0
  call cheb_fit(u,f,0,model,status)
  call check(status==CHEB_SUCCESS.and.model%rank==1,'constant fit and rank')
  constant_model=model
  call cheb_evaluate(model,cmplx(20.d0,3.d0,dp),result,status,dresult)
  call close(result,cmplx(3.25d0,0.d0,dp),3.d-15,'constant complex extrapolation')
  call check(status==CHEB_SUCCESS.and.dresult==(0.d0,0.d0),'degree-zero derivative is exactly zero')
  call cheb_evaluate(model,huge(1.d0),value,status,derivative)
  call check(status==CHEB_SUCCESS.and.derivative==0.d0,'constant avoids unnecessary coordinate overflow')
  call cheb_fit(u,0.d0*u,0,model,status)
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==0.d0.and.derivative==0.d0,'legitimate zero polynomial succeeds')
  call check(model%relative_residual==0.d0,'zero-data relative residual is defined')

  f=-0.7d0+2.d0*u
  call cheb_fit(u,f,1,model,status)
  call check(status==CHEB_SUCCESS.and.model%rank==2,'linear fit and rank')
  linear_model=model
  call check(abs(model%center-0.5d0)<1.d-15.and.abs(model%scale-2.5d0)<1.d-15,'actual-coordinate mapping')
  do i=1,5
     z=cmplx(-5.d0+2.d0*i,0.2d0*i,dp)
     call cheb_evaluate(model,z,result,status,dresult)
     call close(result,-0.7d0+2.d0*z,2.d-14,'linear complex value')
     call close(dresult,cmplx(2.d0,0.d0,dp),2.d-14,'linear physical derivative')
     call check(status==CHEB_SUCCESS,'linear complex evaluation status')
  enddo
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call close(cmplx(value,0.d0,dp),cmplx(1.3d0,0.d0,dp),2.d-14,'real generic value')
  call check(abs(derivative-2.d0)<2.d-14,'real generic derivative')

  ! Reverse a nonuniform array; expected function and derivative use powers,
  ! independently of the Chebyshev recurrence used by the implementation.
  u=u(9:1:-1); f=2.d0-3.d0*u+0.5d0*u*u
  call cheb_fit(u,f,2,model,status)
  call check(status==CHEB_SUCCESS.and.model%effective_degree==2,'unsorted nonuniform quadratic')
  z=cmplx(4.d0,-0.7d0,dp)
  call cheb_evaluate(model,z,result,status,dresult)
  call close(result,2.d0-3.d0*z+0.5d0*z*z,4.d-14,'quadratic outside interval')
  call close(dresult,-3.d0+z,4.d-14,'quadratic derivative and constant coefficient convention')
  call check(status==CHEB_SUCCESS,'quadratic extrapolation status')
  rnorm=0.d0
  do i=1,9
     call cheb_evaluate(model,u(i),value,status)
     rnorm=rnorm+(value-f(i))**2
  enddo
  call check(model%residual_norm<2.d-13.and.sqrt(rnorm)<2.d-13,'exact-fit residual diagnostics')
  call check(model%svd_tolerance==9.d0*epsilon(1.d0),'documented dimension-scaled default RCOND')

  ! A least-squares residual which is not zero, with an independent mean/RMS.
  call cheb_fit([0.d0,1.d0,2.d0],[1.d0,2.d0,6.d0],0,model,status)
  call cheb_evaluate(model,0.d0,value,status)
  expected=3.d0
  call check(abs(value-expected)<2.d-15,'constant least-squares mean')
  call check(abs(model%residual_norm-sqrt(14.d0))<1.d-14,'residual L2 norm')
  call check(abs(model%rms_residual-sqrt(14.d0/3.d0))<1.d-14,'RMS residual')
  call check(abs(model%relative_residual-sqrt(14.d0/41.d0))<1.d-14,'relative residual')
  call check(model%residual_norm_status==CHEB_DIAG_VALID.and.model%rms_residual_status==CHEB_DIAG_VALID.and. &
       model%relative_residual_status==CHEB_DIAG_VALID,'ordinary diagnostics all available')

  ! Mixed cubic: P'= (c1+3*c3)T0 + 4*c2*T1 + 6*c3*T2, then divide by h.
  nodes(1:6)=[-1.d0,-.7d0,-.2d0,.3d0,.8d0,1.d0]
  samples(1:6)=2.d0+3.d0*nodes(1:6)+5.d0*(2.d0*nodes(1:6)**2-1.d0)+ &
       7.d0*(4.d0*nodes(1:6)**3-3.d0*nodes(1:6))
  call cheb_fit(2.d0+3.d0*nodes(1:6),samples(1:6),3,model,status)
  call check(status==CHEB_SUCCESS,'independent mixed cubic fit')
  call check(maxval(abs(model%derivative_coefficients-[24.d0,20.d0,42.d0]))<1.d-12,'mixed cubic derivative coefficients')
  z=cmplx(.2d0,.1d0,dp)
  call cheb_evaluate(model,2.d0+3.d0*z,result,status,dresult)
  call close(dresult,(24.d0+20.d0*z+42.d0*(2.d0*z*z-1.d0))/3.d0,1.d-13,'mixed cubic complex physical derivative')

  ! Known Chebyshev modes are generated/evaluated with independent trigonometric
  ! identities, including complex points and analytic endpoint derivatives.
  do n=1,4
     degree=degrees(n)
     do i=1,41
        nodes(i)=3.d0+2.d0*cos(pi*dble(i-1)/40.d0)
        samples(i)=cos(dble(degree)*acos(max(-1.d0,min(1.d0,(nodes(i)-3.d0)/2.d0))))
     enddo
     call cheb_fit(nodes,samples,degree,model,status)
     call check(status==CHEB_SUCCESS.and.model%rank==degree+1,'higher-order full-rank fit')
     call check(abs(model%coefficients(degree)-1.d0)<2.d-12,'known highest Chebyshev coefficient')
     do j=0,degree-1
        call check(abs(model%coefficients(j))<2.d-12,'known zero Chebyshev coefficient')
     enddo
     do j=1,5
        select case(j)
        case(1); z=cmplx(0.2d0,0.d0,dp)
        case(2); z=cmplx(0.3d0,-0.2d0,dp)
        case(3); z=cmplx(1.1d0,0.05d0,dp)
        case(4); z=cmplx(1.d0,0.d0,dp)
        case(5); z=cmplx(-1.d0,0.d0,dp)
        end select
        want=cos(dble(degree)*acos(z))
        if(j==4) then
           dwant=cmplx(dble(degree*degree)/2.d0,0.d0,dp)
        elseif(j==5) then
           dwant=cmplx(dble((-1)**(degree-1)*degree*degree)/2.d0,0.d0,dp)
        else
           dwant=dble(degree)*sin(dble(degree)*acos(z))/sqrt(1.d0-z*z)/2.d0
        endif
        call cheb_evaluate(model,3.d0+2.d0*z,result,status,dresult)
        call close(result,want,3.d-11,'higher-degree independent complex value')
        call close(dresult,dwant,3.d-11,'higher-degree independent derivative')
        call check(status==CHEB_SUCCESS,'higher-degree evaluation status')
     enddo
  enddo
  z=cmplx(3.4d0,0.2d0,dp); step=1.d-6
  call cheb_evaluate(model,z,result,status,dresult)
  call cheb_evaluate(model,z+step,plus,status)
  call cheb_evaluate(model,z-step,minus,status)
  call close((plus-minus)/(2.d0*step),dresult,2.d-8,'finite-difference derivative cross-check')

  call cheb_fit([-1.d0,0.d0,1.d0],[1.d0,0.d0,1.d0],20,model,status)
  call check(status==CHEB_REDUCED.and.model%effective_degree==2.and.model%sample_limited,'sample-limited reduction')
  call check(model%requested_degree==20.and.model%sample_count==3,'requested degree and sample metadata')
  call cheb_fit([7.d0],[2.d0],20,model,status)
  call check(status==CHEB_REDUCED.and.model%effective_degree==0.and.model%scale>0.d0,'singleton constant mapping')
  call cheb_evaluate(model,cmplx(-7.d0,9.d0,dp),result,status,dresult)
  call close(result,cmplx(2.d0,0.d0,dp),1.d-15,'singleton extrapolation')
  call check(dresult==(0.d0,0.d0),'singleton derivative')

  call cheb_fit([2.d0,2.d0,2.d0,2.d0],[1.d0,2.d0,3.d0,6.d0],3,model,status)
  call check(status==CHEB_REDUCED.and.model%rank_reduced.and.model%degenerate_interval,'degenerate rank reduction')
  call check(model%effective_degree==0.and.model%rank==1.and.model%initial_rank==1,'degenerate numerical rank')
  call check(model%fit_attempts==2,'rank-deficient model is refitted')
  call cheb_evaluate(model,0.d0,value,status)
  call check(abs(value-3.d0)<3.d-15,'reduced constant equals independently computed mean')

  ! Two distinct nodes with repetitions: a degree-one refit must give these
  ! group means, not coefficients truncated from a cubic minimum-norm solution.
  call cheb_fit([-1.d0,-1.d0,1.d0,1.d0],[1.d0,3.d0,6.d0,8.d0],3,model,status)
  call check(status==CHEB_REDUCED.and.model%effective_degree==1.and.model%fit_attempts==2,'repeated nodes refit linear')
  call cheb_evaluate(model,-1.d0,value,status)
  call check(abs(value-2.d0)<5.d-15,'refitted left group mean')
  call cheb_evaluate(model,1.d0,value,status)
  call check(abs(value-7.d0)<5.d-15,'refitted right group mean')
  call check(abs(model%residual_norm-2.d0)<5.d-15,'repeated-node residual')
  call cheb_evaluate(model,cmplx(.3d0,.2d0,dp),result,status,dresult)
  call close(dresult,cmplx(2.5d0,0.d0,dp),3.d-15,'derivative after rank-limited refit')

  cluster=[-1.d0,-1.d0+1.d-6,-1.d0+2.d-6,1.d0-2.d-6,1.d0-1.d-6,1.d0]
  target=cluster**2
  call cheb_fit(cluster,target,2,model,status,rcond=1.d-4)
  call check(status==CHEB_REDUCED.and.model%rank_reduced,'explicit loose rank tolerance reduces degree')
  loose_residual=model%residual_norm
  call cheb_fit(cluster,target,2,reference,status,rcond=1.d-12)
  call check(status==CHEB_SUCCESS.and.reference%rank==3,'tighter rank tolerance retains full rank')
  call check(reference%residual_norm<loose_residual/100.d0,'rank tolerance residual diagnostics')
  call check(model%svd_tolerance==1.d-4.and.reference%svd_tolerance==1.d-12,'selected rank thresholds recorded')

  do i=1,9
     u(i)=1000.d0+dble(i-1)*1.d-10
  enddo
  f=1.d0+3.d0*(u-u(1))/(u(9)-u(1))
  call cheb_fit(u,f,1,model,status)
  call check(status==CHEB_SUCCESS.and.model%scale>0.d0,'narrow nondegenerate interval')
  do i=1,9
     call cheb_evaluate(model,u(i),value,status)
     call check(abs(value-f(i))<2.d-14,'narrow interval reconstruction')
  enddo
  call cheb_fit([-huge(1.d0),huge(1.d0)],[1.d0,1.d0],0,model,status)
  call check(status==CHEB_SUCCESS.and.model%scale==huge(1.d0),'large endpoints avoid mapping overflow')
  call cheb_fit([-1.d0,1.d0],[huge(1.d0)/2.d0,huge(1.d0)/2.d0],0,model,status)
  call check(status==CHEB_SUCCESS,'large values avoid RHS/norm overflow')

  ! Failed refits must invalidate an earlier good model.
  call cheb_fit([real(dp)::],[real(dp)::],0,model,status)
  call check(status==CHEB_NO_SAMPLES.and..not.allocated(model%coefficients),'empty input invalidates output')
  call cheb_fit([1.d0],[1.d0],-1,model,status)
  call check(status==CHEB_INVALID_DEGREE,'invalid degree')
  call cheb_fit([1.d0,2.d0],[1.d0],1,model,status)
  call check(status==CHEB_INVALID_VALUES,'mismatched sample lengths')
  call cheb_fit([nan],[1.d0],0,model,status)
  call check(status==CHEB_INVALID_COORDINATES,'NaN sample coordinate')
  call cheb_fit([inf],[1.d0],0,model,status)
  call check(status==CHEB_INVALID_COORDINATES,'infinite sample coordinate')
  call cheb_fit([1.d0],[nan],0,model,status)
  call check(status==CHEB_INVALID_VALUES,'NaN sample value')
  call cheb_fit([1.d0],[inf],0,model,status)
  call check(status==CHEB_INVALID_VALUES,'infinite sample value')
  call cheb_fit([1.d0],[1.d0],0,model,status,rcond=-1.d0)
  call check(status==CHEB_INVALID_TOLERANCE,'negative RCOND rejected')
  call cheb_fit([1.d0],[1.d0],0,model,status,rcond=1.d0)
  call check(status==CHEB_INVALID_TOLERANCE,'unit RCOND rejected')
  call cheb_fit([1.d0],[1.d0],0,model,status,rcond=nan)
  call check(status==CHEB_INVALID_TOLERANCE,'nonfinite RCOND rejected')
  call cheb_fit([0.d0,ieee_next_after(0.d0,1.d0)],[1.d0,2.d0],1,model,status)
  call check(status==CHEB_INVALID_NORMALIZATION,'unrepresentable half-width rejected')
  call cheb_fit([-1.d0,0.d0,1.d0],[1.d308,-1.d308,1.d308],2,model,status)
  call check(status==CHEB_NONFINITE_COEFFICIENTS,'unrepresentable derivative coefficients rejected')

  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_INVALID_MODEL.and.ieee_is_nan(value),'failed model cannot evaluate as zero')
  model=linear_model; model%scale=0.d0
  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_INVALID_MODEL.and.ieee_is_nan(value),'zero-scale model rejected')
  model=linear_model; model%coefficients(0)=nan
  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_INVALID_MODEL.and.ieee_is_nan(value),'nonfinite model coefficient rejected')
  model=linear_model; deallocate(model%derivative_coefficients)
  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_INVALID_MODEL,'incomplete model rejected')
  call cheb_evaluate(linear_model,nan,value,status)
  call check(status==CHEB_INVALID_COORDINATES.and.ieee_is_nan(value),'nonfinite real argument')
  call cheb_evaluate(linear_model,cmplx(0.d0,inf,dp),result,status)
  call check(status==CHEB_INVALID_COORDINATES.and.ieee_is_nan(real(result,dp)),'nonfinite complex argument')
  call cheb_evaluate(linear_model,huge(1.d0),value,status)
  call check(status==CHEB_EVALUATION_FAILURE.and.ieee_is_nan(value),'overflow returns explicit failure under traps')
  model=linear_model; model%center=0.d0; model%scale=tiny(1.d0)/16.d0
  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_SUCCESS,'value-only call skips derivative evaluation')
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call check(status==CHEB_EVALUATION_FAILURE.and.ieee_is_nan(derivative),'physical derivative overflow is explicit')

  ! R1: rounded H-b and H/b guards previously passed these overflowing cases.
  ! Models here specify exact coefficients; every operation is exercised through
  ! the public evaluator, including the real wrapper and complex component paths.
  hmax=huge(1.d0); a=nearest(hmax,-1.d0); gap=spacing(hmax); b=1.5d0*gap
  call linear_polynomial(model,a,b)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call failed_real(value,derivative,status,'rounded-bound final addition')
  call cheb_evaluate(model,cmplx(1.d0,1.d0,dp),result,status,dresult)
  call failed_complex(result,dresult,status,'complex rounded-bound final addition')
  call linear_polynomial(model,a,gap)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax.and.derivative==gap,'nearby exact H addition retained')
  call linear_polynomial(model,-a,-gap)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==-hmax.and.derivative==-gap,'negative exact boundary retained')
  call linear_polynomial(model,-a,-b)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call failed_real(value,derivative,status,'negative rounded-bound addition')
  call linear_polynomial(model,hmax,-hmax/2.d0)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax/2.d0,'opposite-sign near-limit addition')

  ! Exact binary64 fixture from the review: both factors finite, rounded H/b
  ! equals the first factor, but their product overflows. TRANSFER only supplies
  ! reproducible test data; it is not signaling-NaN classification machinery.
  multiply_a=transfer(int(z'7FE388D6D20C39B7',int64),1.d0)
  multiply_b=transfer(int(z'3FFA35C2A3BC0F30',int64),1.d0)
  call linear_polynomial(model,0.d0,multiply_b)
  call cheb_evaluate(model,multiply_a,value,status,derivative)
  call failed_real(value,derivative,status,'rounded-bound Clenshaw multiplication')
  call cheb_evaluate(model,cmplx(0.d0,multiply_a,dp),result,status,dresult)
  call failed_complex(result,dresult,status,'imaginary multiplication boundary')
  nearby=nearest(multiply_a,-1.d0)
  call cheb_evaluate(model,nearby,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(value),'nearby representable multiplication')
  call cheb_evaluate(model,cmplx(0.d0,nearby,dp),result,status,dresult)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(aimag(result)),'nearby imaginary multiplication')
  call linear_polynomial(model,0.d0,hmax)
  call cheb_evaluate(model,1.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax.and.derivative==hmax,'unit product at H')
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==0.d0.and.derivative==hmax,'zero product at H')
  call linear_polynomial(model,0.d0,hmax/2.d0)
  call cheb_evaluate(model,2.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax,'non-unit exact product at H')

  call linear_polynomial(model,0.d0,1.d0)
  model%center=-b
  call cheb_evaluate(model,a,value,status,derivative)
  call failed_real(value,derivative,status,'coordinate subtraction overflow')
  model%center=-gap
  call cheb_evaluate(model,a,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax,'nearby exact coordinate subtraction')
  model%center=0.d0; model%scale=.5d0
  call cheb_evaluate(model,hmax/2.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.value==hmax.and.derivative==2.d0,'exact division at H')
  call cheb_evaluate(model,nearest(hmax/2.d0,1.d0),value,status,derivative)
  call failed_real(value,derivative,status,'coordinate division overflow')
  call cheb_evaluate(model,cmplx(0.d0,hmax/2.d0,dp),result,status,dresult)
  call check(status==CHEB_SUCCESS.and.aimag(result)==hmax,'imaginary coordinate division at H')
  call cheb_evaluate(model,cmplx(0.d0,nearest(hmax/2.d0,1.d0),dp),result,status,dresult)
  call failed_complex(result,dresult,status,'imaginary coordinate division overflow')
  model%scale=1.d0
  call cheb_evaluate(model,cmplx(hmax,hmax,dp),result,status,dresult)
  call check(status==CHEB_SUCCESS.and.real(result,dp)==hmax.and.aimag(result)==hmax,'finite complex components, large modulus')

  ! A finite polynomial may have an unrepresentable intermediate: reject safely,
  ! rather than overflowing while accumulating a complex recurrence coefficient.
  call linear_polynomial(model,0.d0,1.d0)
  deallocate(model%coefficients,model%derivative_coefficients)
  model%requested_degree=2; model%effective_degree=2; model%rank=3; model%sample_count=3
  allocate(model%coefficients(0:2),model%derivative_coefficients(0:1))
  model%coefficients=[0.d0,a,b]; model%derivative_coefficients=[a,4.d0*b]
  call cheb_evaluate(model,cmplx(.5d0,.5d0,dp),result,status,dresult)
  call failed_complex(result,dresult,status,'complex recurrence intermediate addition')
  model%coefficients(2)=gap; model%derivative_coefficients(1)=4.d0*gap
  call cheb_evaluate(model,cmplx(.5d0,.5d0,dp),result,status)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(real(result,dp)).and.ieee_is_finite(aimag(result)), &
       'nearby finite complex recurrence, value only')
  call cheb_evaluate(model,cmplx(.5d0,.5d0,dp),result,status,dresult)
  call failed_complex(result,dresult,status,'finite value but overflowing complex derivative')
  call linear_polynomial(model,0.d0,hmax/2.d0)
  model%scale=.5d0
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.derivative==hmax,'physical derivative exact boundary')
  model%coefficients(1)=nearest(hmax/2.d0,1.d0); model%derivative_coefficients(0)=model%coefficients(1)
  call cheb_evaluate(model,0.d0,value,status)
  call check(status==CHEB_SUCCESS.and.value==0.d0,'value-only call at derivative overflow boundary')
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call failed_real(value,derivative,status,'physical derivative boundary overflow')

  ! R3: valid full-rank fit, unrepresentable L2 norm, finite RMS and relative norm.
  call cheb_fit([-1.d0,1.d0],[-hmax,hmax],0,model,status)
  call check(status==CHEB_SUCCESS.and.model%fit_status==CHEB_SUCCESS,'diagnostic range does not invalidate fit')
  call check(model%rank==1.and.model%lapack_info==0.and.model%effective_degree==0,'diagnostic range preserves rank and degree')
  call check(model%residual_norm_status==CHEB_DIAG_RANGE.and.ieee_is_nan(model%residual_norm),'unrepresentable L2 diagnostic flagged')
  call check(model%rms_residual_status==CHEB_DIAG_VALID,'extreme RMS diagnostic available')
  if(model%rms_residual_status==CHEB_DIAG_VALID) call check(abs(model%rms_residual/hmax-1.d0)<3.d-15,'extreme RMS value')
  call check(model%relative_residual_status==CHEB_DIAG_VALID,'extreme relative diagnostic available')
  if(model%relative_residual_status==CHEB_DIAG_VALID) call check(abs(model%relative_residual-1.d0)<3.d-15,'extreme relative value')
  call cheb_evaluate(model,0.d0,value,status,derivative)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(value).and.derivative==0.d0,'fit remains evaluable after diagnostic overflow')
  call cheb_fit([0.d0,0.d0],[-hmax,hmax],5,model,status)
  call check(status==CHEB_REDUCED.and.model%fit_status==CHEB_REDUCED.and.model%rank_reduced, &
       'diagnostic range preserves adaptive-reduction success')
  call check(model%residual_norm_status==CHEB_DIAG_RANGE,'reduced fit retains diagnostic availability')
  smallest=ieee_next_after(0.d0,1.d0)
  call cheb_fit([0.d0,1.d0,2.d0,3.d0,4.d0],[smallest,0.d0,0.d0,0.d0,0.d0],0,model,status)
  call check(status==CHEB_SUCCESS,'fit survives diagnostic underflow')
  call check(model%rms_residual_status==CHEB_DIAG_RANGE.and.ieee_is_nan(model%rms_residual), &
       'unrepresentable positive RMS is not reported as zero')
  call cheb_fit([0.d0,1.d0],[hmax,smallest],0,model,status)
  call check(status==CHEB_SUCCESS,'fit survives unavailable diagnostic computation')
  call check(model%residual_norm_status==CHEB_DIAG_FAILED.and.ieee_is_nan(model%residual_norm), &
       'unrepresentable residual rescaling marks computation unavailable')

  call ieee_get_halting_mode(ieee_invalid,halting_after(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,halting_after(2))
  call ieee_get_halting_mode(ieee_overflow,halting_after(3))
  call check(all(halting_before.eqv.halting_after),'caller halting modes unchanged')
  write(*,'(a,3l2)') 'IEEE halting modes (invalid,zero,overflow):',halting_after

  write(*,'(a,i0,a,i0)') 'Chebyshev kernel: passed=',passed,', failed=',failed
  if(failed/=0) error stop 1
contains
  subroutine linear_polynomial(m,c0,c1)
    type(cheb_model), intent(out) :: m
    real(dp), intent(in) :: c0,c1
    m%requested_degree=1; m%effective_degree=1; m%rank=2; m%sample_count=2; m%fit_status=CHEB_SUCCESS
    allocate(m%coefficients(0:1),m%derivative_coefficients(0:0))
    m%coefficients=[c0,c1]; m%derivative_coefficients=c1
  end subroutine
  subroutine failed_real(v,dv,code,label)
    real(dp), intent(in) :: v,dv
    integer, intent(in) :: code
    character(*), intent(in) :: label
    call check(code==CHEB_EVALUATION_FAILURE,label//' status')
    call check(ieee_is_nan(v).and.ieee_is_nan(dv),label//' atomic NaN outputs')
  end subroutine
  subroutine failed_complex(v,dv,code,label)
    complex(dp), intent(in) :: v,dv
    integer, intent(in) :: code
    character(*), intent(in) :: label
    call check(code==CHEB_EVALUATION_FAILURE,label//' status')
    call check(ieee_is_nan(real(v,dp)).and.ieee_is_nan(aimag(v)).and. &
         ieee_is_nan(real(dv,dp)).and.ieee_is_nan(aimag(dv)),label//' atomic NaN outputs')
  end subroutine
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if(condition) then
       passed=passed+1
    else
       failed=failed+1
       write(*,'(a)') 'FAIL: '//label
    endif
  end subroutine
  subroutine close(actual,expected,tolerance,label)
    complex(dp), intent(in) :: actual,expected
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label
    logical :: finite
    finite=ieee_is_finite(real(actual,dp)).and.ieee_is_finite(aimag(actual))
    if(.not.finite) then
       call check(.false.,label)
       return
    endif
    call check(abs(actual-expected)<=tolerance*max(1.d0,abs(expected)),label)
  end subroutine
end program test_chebyshev_kernel
