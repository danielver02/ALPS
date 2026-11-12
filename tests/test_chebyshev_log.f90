! Independent analytic references for the generalized LOG10 kernel.
program test_chebyshev_log
  use alps_analyt
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: dp=kind(1.d0)
  type(cheb_model) :: model,saved,direct
  integer :: status,vs,ds,i,passed=0,failed=0
  real(dp), parameter :: ln10=log(10.d0)
  real(dp) :: u(9),f(9),v,d,p,pu,step,err1,err2,hmax,smallest,bound,nan,inf
  complex(dp) :: z,value,derivative,want,dwant,raw,draw,plus,minus
  logical :: before(3),after(3)

  call ieee_get_halting_mode(ieee_invalid,before(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,before(2))
  call ieee_get_halting_mode(ieee_overflow,before(3))
  nan=ieee_value(0.d0,ieee_quiet_nan); inf=ieee_value(0.d0,ieee_positive_inf)
  hmax=huge(1.d0); smallest=nearest(0.d0,1.d0)
  u=[-2.d0,-1.1d0,-.2d0,.3d0,1.1d0,1.8d0,2.1d0,3.d0,4.d0]
  f=10.d0**(.3d0-.7d0*u)
  call cheb_fit(u,f,1,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_SUCCESS.and.model%representation==CHEB_LOG10,'linear logarithmic fit')
  call check(model%original_sample_count==9.and.model%retained_sample_count==9.and.model%sample_count==9, &
       'unfiltered sample counts')
  call check(model%center==1.d0.and.model%scale==3.d0,'non-unit physical mapping')
  call check(model%residual_norm_status==CHEB_DIAG_VALID.and.model%residual_norm<1.d-13,'log-space residual')
  saved=model
  do i=1,6
     select case(i)
     case(1); z=cmplx(-1.d0,0.d0,dp)
     case(2); z=cmplx(.2d0,1.d-8,dp)
     case(3); z=cmplx(1.2d0,.7d0,dp)
     case(4); z=cmplx(-4.d0,-1.1d0,dp)
     case(5); z=cmplx(5.d0,2.d0,dp)
     case(6); z=cmplx(0.d0,-3.d0,dp)
     end select
     want=exp(ln10*(.3d0-.7d0*z)); dwant=-.7d0*ln10*want
     call cheb_evaluate(model,z,value,status,derivative,vs,ds)
     call check(status==CHEB_SUCCESS.and.vs==CHEB_SUCCESS.and.ds==CHEB_SUCCESS,'complex linear status')
     call close(value,want,3.d-14,'complex exponential analytic value')
     call close(derivative,dwant,4.d-14,'complex exponential analytic derivative')
     call cheb_evaluate_polynomial(model,z,raw,status,draw)
     call close(raw,.3d0-.7d0*z,2.d-14,'raw log polynomial')
     call close(draw,cmplx(-.7d0,0.d0,dp),2.d-14,'raw physical log derivative')
     call check(status==CHEB_SUCCESS,'raw log evaluation status')
  enddo
  call cheb_evaluate(model,.5d0,v,status,d,vs,ds)
  call relative(v,10.d0**(-.05d0),3.d-14,'real logarithmic reconstruction')
  call relative(d,-.7d0*ln10*10.d0**(-.05d0),4.d-14,'real logarithmic derivative')
  call cheb_evaluate_polynomial(model,.5d0,p,status,pu)
  call check(abs(p+.05d0)<2.d-15.and.abs(pu+.7d0)<2.d-15,'raw real generic')
  call cheb_evaluate(model,.5d0,v,status,value_status=vs,derivative_status=ds)
  call check(status==CHEB_SUCCESS.and.ds==CHEB_NOT_REQUESTED,'value-only status metadata')

  ! log10(A exp(-alpha*(u-u0)^2)) is quadratic. No production reference fitter.
  f=2.7d0*exp(-.4d0*(u-.6d0)**2)
  call cheb_fit(u(9:1:-1),f(9:1:-1),2,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_SUCCESS.and.model%rank==3,'unsorted Gaussian log fit')
  do i=1,5
     z=cmplx(-2.d0+1.4d0*i,.15d0*(i-2),dp)
     want=2.7d0*exp(-.4d0*(z-.6d0)**2); dwant=-.8d0*(z-.6d0)*want
     call cheb_evaluate(model,z,value,status,derivative)
     call check(status==CHEB_SUCCESS,'Gaussian continuation status')
     call close(value,want,5.d-14,'Gaussian analytic complex value')
     call close(derivative,dwant,7.d-14,'Gaussian analytic complex derivative')
  enddo
  z=cmplx(1.3d0,.2d0,dp)
  call cheb_evaluate(model,z,value,status,derivative)
  step=1.d-3
  call cheb_evaluate(model,z+step,plus,status)
  call cheb_evaluate(model,z-step,minus,status)
  err1=abs((plus-minus)/(2.d0*step)-derivative)
  step=step/2.d0
  call cheb_evaluate(model,z+step,plus,status)
  call cheb_evaluate(model,z-step,minus,status)
  err2=abs((plus-minus)/(2.d0*step)-derivative)
  call check(err2<.3d0*err1.and.err2<1.d-6,'centered finite difference converges quadratically')

  ! Unsorted zeros on both tails; selection uses coordinates, never indices.
  u=[4.d0,-4.d0,1.d0,-1.d0,0.d0,-3.d0,3.d0,2.d0,-2.d0]
  f=[0.d0,0.d0,10.d0,.1d0,1.d0,0.d0,0.d0,0.d0,0.d0]
  call cheb_fit(u,f,9,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_REDUCED.and.model%sample_limited.and.model%effective_degree==2,'degree after filtering')
  call check(model%original_sample_count==9.and.model%retained_sample_count==3.and.model%sample_count==3, &
       'filtered sample diagnostics')
  call check(model%lower_tail_zeros==3.and.model%upper_tail_zeros==3,'physical lower/upper zero counts')
  call check(model%retained_min==-1.d0.and.model%retained_max==1.d0,'positive coordinate interval')
  call check(model%rank==3.and.model%requested_degree==9,'accepted retained rank')
  call cheb_evaluate(model,3.d0,v,status,d)
  call relative(v,1000.d0,3.d-13,'retained interval is not physical support')
  call relative(d,ln10*1000.d0,3.d-13,'outside-retained-interval derivative')

  do i=1,2
     if(i==1) then
        call cheb_fit([-2.d0,-1.d0,0.d0,1.d0],[0.d0,0.d0,1.d0,10.d0],1,model,status,representation=CHEB_LOG10)
        call check(model%lower_tail_zeros==2.and.model%upper_tail_zeros==0,'lower tail only')
     else
        call cheb_fit([-1.d0,0.d0,1.d0,2.d0],[.1d0,1.d0,0.d0,0.d0],1,model,status,representation=CHEB_LOG10)
        call check(model%lower_tail_zeros==0.and.model%upper_tail_zeros==2,'upper tail only')
     endif
     call check(status==CHEB_SUCCESS,'one-sided tail success')
  enddo
  call cheb_fit([2.d0,-1.d0,0.d0,1.d0],[0.d0,0.d0,7.d0,0.d0],8,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_REDUCED.and.model%effective_degree==0.and.model%retained_sample_count==1,'singleton positive')
  call check(model%lower_tail_zeros==1.and.model%upper_tail_zeros==2.and.model%scale==1.d0,'singleton tail metadata')
  call cheb_evaluate(model,cmplx(hmax,hmax,dp),value,status,derivative)
  call close(value,cmplx(7.d0,0.d0,dp),2.d-15,'singleton ignores irrelevant huge coordinates')
  call check(status==CHEB_SUCCESS.and.derivative==(0.d0,0.d0),'constant logarithm has exact zero derivative')
  call cheb_fit([1.d0,0.d0,1.d0,0.d0],[10.d0,1.d0,10.d0,1.d0],3,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_REDUCED.and.model%rank_reduced.and.model%effective_degree==1,'repeated nodes rank reduction')
  call check(model%rank==2.and.model%initial_rank==2.and.model%fit_attempts==2,'accepted reduced rank diagnostics')
  call cheb_evaluate(model,.5d0,v,status,d)
  call relative(d,ln10*sqrt(10.d0),3.d-14,'rank-reduced log derivative')
  call cheb_fit([2.d0,2.d0],[1.d0,100.d0],3,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_REDUCED.and.model%degenerate_interval.and.model%effective_degree==0,'degenerate log coordinates')
  call cheb_evaluate(model,4.d0,v,status,d)
  call relative(v,10.d0,2.d-15,'degenerate row mean in log space')
  call check(d==0.d0,'degenerate row derivative')
  call cheb_fit([1.d0,nearest(1.d0,1.d0)],[1.d0,10.d0],1,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_SUCCESS.and..not.model%degenerate_interval,'narrow nondegenerate retained interval')
  call cheb_evaluate_polynomial(model,1.d0,p,status,pu)
  call relative(pu,1.d0/spacing(1.d0),3.d-14,'narrow interval physical derivative scaling')

  call cheb_fit([0.d0,1.d0,2.d0],[1.d0,0.d0,10.d0],2,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INTERIOR_ZERO,'interior zero')
  call cheb_fit([0.d0,0.d0],[1.d0,0.d0],1,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INTERIOR_ZERO,'coincident zero and positive are inconsistent')
  call cheb_fit([0.d0,1.d0],[-1.d0,1.d0],1,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_NEGATIVE_SAMPLES,'negative sample')
  call cheb_fit([0.d0,1.d0],[0.d0,-0.d0],1,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_NO_POSITIVE_SAMPLES,'all-zero row')
  call cheb_fit([real(dp)::],[real(dp)::],1,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_NO_POSITIVE_SAMPLES,'empty log row')
  call cheb_fit([nan],[1.d0],0,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_COORDINATES,'NaN coordinate')
  call cheb_fit([inf],[1.d0],0,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_COORDINATES,'infinite coordinate')
  call cheb_fit([0.d0],[nan],0,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_VALUES,'NaN sample')
  call cheb_fit([0.d0],[inf],0,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_VALUES,'infinite sample')
  call cheb_fit([0.d0],[1.d0],-1,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_DEGREE,'negative degree')
  call cheb_fit([0.d0],[1.d0,2.d0],0,model,status,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_VALUES,'mismatched shapes')
  call cheb_fit([0.d0],[1.d0],0,model,status,representation=99)
  call bad_fit(CHEB_INVALID_REPRESENTATION,'invalid fit representation')
  call cheb_fit([0.d0],[1.d0],0,model,status,rcond=1.d0,representation=CHEB_LOG10)
  call bad_fit(CHEB_INVALID_TOLERANCE,'log fit retains tolerance validation')

  call cheb_fit([0.d0],[smallest],0,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_SUCCESS.and.model%retained_sample_count==1,'positive least subnormal is retained')
  call check(abs(model%coefficients(0)-log10(smallest))<1.d-12,'subnormal was not floored')
  call cheb_evaluate(model,0.d0,v,status,d)
  call check((status==CHEB_SUCCESS.or.status==CHEB_LOG_UNDERFLOW).and.d==0.d0,'least subnormal reconstruction classified')
  call check(v==smallest.or.(status==CHEB_LOG_UNDERFLOW.and.v==0.d0),'least subnormal is not replaced by positive floor')
  call cheb_fit([0.d0,1.d0],[1.d-310,1.d-309],1,model,status,representation=CHEB_LOG10)
  call check(status==CHEB_SUCCESS,'two positive subnormal samples')
  call cheb_evaluate(model,.5d0,v,status,d)
  call check(status==CHEB_SUCCESS,'representable subnormal success')
  call relative(v,10.d0**(-309.5d0),3.d-12,'subnormal interpolation')
  call relative(d,ln10*10.d0**(-309.5d0),3.d-12,'subnormal derivative')
  call log_linear(nearest(log10(smallest),1.d0),0.d0)
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_SUCCESS.and.v==smallest,'just above least-subnormal exponent')
  call log_linear(nearest(log10(smallest),-1.d0),0.d0)
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_LOG_UNDERFLOW.and.v==0.d0,'just below least-subnormal exponent')

  ! Construct public models with exactly controlled log coefficients to isolate
  ! reconstruction boundaries from SVD roundoff. All evaluations remain public.
  bound=log10(hmax)-4.d0*spacing(log10(hmax))
  call log_linear(bound,0.d0)
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(v).and.d==0.d0,'near-largest reconstructed value')
  call relative(v,10.d0**bound,3.d-13,'near-largest analytic reconstruction')
  call log_linear(log10(hmax)+4.d0*spacing(log10(hmax)),0.d0)
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_OVERFLOW.and.vs==CHEB_LOG_OVERFLOW,'exponential overflow status')
  call check(ieee_is_nan(v).and.ieee_is_nan(d),'overflow atomic NaNs')
  call cheb_evaluate_polynomial(model,0.d0,p,status,pu)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(p),'overflow leaves raw model evaluable')

  do i=1,3
     select case(i)
     case(1); p=log10(tiny(1.d0))+.1d0
     case(2); p=log10(tiny(1.d0))-.1d0
     case(3); p=-320.d0
     end select
     call log_linear(p,0.d0)
     call cheb_evaluate(model,0.d0,v,status,d)
     call check(status==CHEB_SUCCESS.and.v>0.d0.and.d==0.d0,'normal/subnormal representable value')
     call relative(v,10.d0**p,3.d-12,'normal/subnormal analytic reference')
     if(i==1) then
        call check(ieee_class(v)==ieee_positive_normal,'normal classification')
     else
        call check(ieee_class(v)==ieee_positive_denormal,'subnormal classification')
     endif
  enddo
  call log_linear(-400.d0,0.d0)
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_UNDERFLOW.and.vs==CHEB_LOG_UNDERFLOW.and.ds==CHEB_SUCCESS,'underflow differs from zero derivative')
  call check(v==0.d0.and.d==0.d0.and.model%fit_status==CHEB_SUCCESS,'underflow does not invalidate fit')
  call log_linear(-hmax,0.d0)
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_UNDERFLOW.and.v==0.d0.and.d==0.d0.and.ds==CHEB_SUCCESS,'extreme finite negative exponent')
  call log_linear(hmax,0.d0)
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_LOG_OVERFLOW.and.ieee_is_nan(v).and.ieee_is_nan(d),'extreme finite positive exponent')

  call log_linear(0.d0,1.d0)
  z=cmplx(0.d0,1.d6,dp)
  call cheb_evaluate(model,z,value,status,derivative)
  want=exp(cmplx(0.d0,ln10*1.d6,dp))
  call check(status==CHEB_SUCCESS,'large but resolvable complex phase')
  call close(value,want,3.d-14,'large phase reference')
  call close(derivative,ln10*want,3.d-14,'large phase derivative')
  call cheb_evaluate(model,cmplx(0.d0,hmax,dp),value,status,derivative)
  call check(status==CHEB_PHASE_FAILURE,'unreliable extreme phase explicitly rejected')
  call complex_nan(value,derivative,'phase failure atomic NaNs')
  call cheb_evaluate(model,cmplx(0.d0,1.d9,dp),value,status)
  call check(status==CHEB_PHASE_FAILURE,'phase precision threshold')

  call log_linear(308.d0,2.d0)
  call cheb_evaluate(model,0.d0,v,status)
  call check(status==CHEB_SUCCESS.and.ieee_is_finite(v),'value-only avoids unrequested derivative overflow')
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_DERIVATIVE_FAILURE.and.vs==CHEB_SUCCESS.and.ds==CHEB_LOG_OVERFLOW,'derivative overflow independently classified')
  call check(ieee_is_nan(v).and.ieee_is_nan(d),'derivative overflow atomic NaNs')
  call log_linear(308.d0,.1d0)
  call cheb_evaluate(model,cmplx(0.d0,2.d0,dp),value,status,derivative)
  want=1.d308*exp(cmplx(0.d0,.2d0*ln10,dp))
  call check(status==CHEB_SUCCESS,'complex value and derivative near upper range')
  call component_relative(value,want,3.d-13,'near-range complex value')
  call component_relative(derivative,.1d0*ln10*want,3.d-13,'near-range complex derivative')
  call log_linear(-400.d0,1.d308)
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_UNDERFLOW.and.vs==CHEB_LOG_UNDERFLOW.and.ds==CHEB_SUCCESS,'underflow value with finite derivative')
  call check(v==0.d0.and.d>0.d0,'derivative not lost with underflow value')
  call relative(d,ln10*1.d-92,3.d-13,'compensated derivative analytic reference')
  ! ln(10)*f alone would overflow here, but the final derivative is small.
  call log_linear(308.d0,1.d-300)
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_SUCCESS,'scaled derivative avoids unnecessary intermediate overflow')
  call relative(d,ln10*1.d8,3.d-13,'scaled derivative analytic reference')
  call log_linear(-300.d0,1.d-100)
  call cheb_evaluate(model,0.d0,v,status,d,vs,ds)
  call check(status==CHEB_LOG_UNDERFLOW.and.vs==CHEB_SUCCESS.and.ds==CHEB_LOG_UNDERFLOW,'derivative-only underflow')
  call check(v>0.d0.and.d==0.d0,'finite value retained with underflowing derivative')
  call log_linear(-400.d0,1.d308)
  z=cmplx(0.d0,1.d-308,dp)
  call cheb_evaluate(model,z,value,status,derivative,vs,ds)
  want=ln10*1.d-92*exp(cmplx(0.d0,ln10,dp))
  call check(status==CHEB_LOG_UNDERFLOW.and.ds==CHEB_SUCCESS,'complex derivative survives value underflow')
  call component_relative(derivative,want,5.d-13,'complex compensated derivative')

  model=saved; model%representation=999
  call cheb_evaluate(model,0.d0,v,status,d)
  call check(status==CHEB_INVALID_REPRESENTATION.and.ieee_is_nan(v).and.ieee_is_nan(d),'invalid model representation')
  call cheb_evaluate_polynomial(model,0.d0,p,status)
  call check(status==CHEB_INVALID_REPRESENTATION,'raw evaluator validates representation')
  model=saved
  call cheb_evaluate(model,nan,v,status,d)
  call check(status==CHEB_INVALID_COORDINATES.and.ieee_is_nan(v).and.ieee_is_nan(d),'NaN logarithmic coordinate')
  call cheb_evaluate(model,cmplx(0.d0,inf,dp),value,status,derivative)
  call check(status==CHEB_INVALID_COORDINATES,'infinite logarithmic coordinate')
  call complex_nan(value,derivative,'invalid coordinate atomic NaNs')
  ! Existing optional positional RCOND and default direct behavior still work.
  call cheb_fit([-1.d0,0.d0,1.d0],[2.d0,3.d0,4.d0],1,direct,status,1.d-13)
  call check(status==CHEB_SUCCESS.and.direct%representation==CHEB_DIRECT,'source-compatible direct fit')
  call cheb_evaluate(direct,cmplx(.4d0,.2d0,dp),value,status,derivative)
  call cheb_evaluate_polynomial(direct,cmplx(.4d0,.2d0,dp),raw,status,draw)
  call check(value==raw.and.derivative==draw,'direct physical/raw evaluations identical')

  call ieee_get_halting_mode(ieee_invalid,after(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,after(2))
  call ieee_get_halting_mode(ieee_overflow,after(3))
  call check(all(before.eqv.after),'logarithmic kernel preserves caller halting modes')
  write(*,'(a,i0,a,i0)') 'Chebyshev logarithmic kernel: passed=',passed,', failed=',failed
  if(failed/=0) error stop 1
contains
  subroutine log_linear(c0,c1)
    real(dp), intent(in) :: c0,c1
    call cheb_fit([-1.d0,1.d0],[1.d0,10.d0],1,model,status,representation=CHEB_LOG10)
    if(status/=CHEB_SUCCESS) error stop 'test model setup failed'
    model%coefficients=[c0,c1]; model%derivative_coefficients=c1
  end subroutine
  subroutine bad_fit(expected,label)
    integer, intent(in) :: expected
    character(*), intent(in) :: label
    call check(status==expected,label//' fit status')
    call check(model%fit_status==status.and.model%effective_degree==-1.and. &
         .not.allocated(model%coefficients).and..not.allocated(model%derivative_coefficients),label//' invalidated model')
  end subroutine
  subroutine complex_nan(a,b,label)
    complex(dp), intent(in) :: a,b
    character(*), intent(in) :: label
    call check(ieee_is_nan(real(a,dp)).and.ieee_is_nan(aimag(a)).and. &
         ieee_is_nan(real(b,dp)).and.ieee_is_nan(aimag(b)),label)
  end subroutine
  subroutine component_relative(actual,expected,tolerance,label)
    complex(dp), intent(in) :: actual,expected
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label
    call relative(real(actual,dp),real(expected,dp),tolerance,label//' real')
    call relative(aimag(actual),aimag(expected),tolerance,label//' imaginary')
  end subroutine
  subroutine relative(actual,expected,tolerance,label)
    real(dp), intent(in) :: actual,expected,tolerance
    character(*), intent(in) :: label
    if(.not.ieee_is_finite(actual)) then
       call check(.false.,label)
    elseif(expected/=0.d0) then
       call check(abs(actual/expected-1.d0)<=tolerance,label)
    else
       call check(actual==0.d0,label)
    endif
  end subroutine
  subroutine close(actual,expected,tolerance,label)
    complex(dp), intent(in) :: actual,expected
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label
    if(.not.ieee_is_finite(real(actual,dp)).or..not.ieee_is_finite(aimag(actual))) then
       call check(.false.,label)
    else
       call check(abs(actual-expected)<=tolerance*max(1.d0,abs(expected)),label)
    endif
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
end program test_chebyshev_log
