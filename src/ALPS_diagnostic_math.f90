! Numerical experiments only. None of these routines supplies solver values.
module alps_diagnostic_math
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: fitted_sample, sample_log_chebyshev, chebyshev_pair, compare_brackets
  public :: complex_magnitude, complex_log_magnitude, finite_complex, unavailable
  public :: derivative_difference, difference_metrics, guarded_product, safe_ratio
  public :: sample_limited, sample_domain, sample_small, sample_overflow, sample_invalid, sample_underflow
  integer, parameter :: sample_limited=1, sample_domain=2, sample_small=4
  integer, parameter :: sample_overflow=8, sample_invalid=16, sample_underflow=64
  type :: fitted_sample
     double complex :: polynomial, derivative_z, value, analytic
     double precision :: log_function, log_derivative
     integer :: flags=0
     logical :: analytic_available=.false.
  end type
contains
  double precision function unavailable()
    unavailable=ieee_value(0.d0,ieee_quiet_nan)
  end function

  logical function finite_complex(z)
    double complex, intent(in) :: z
    finite_complex=ieee_is_finite(real(z)).and.ieee_is_finite(aimag(z))
  end function

  double precision function complex_log_magnitude(z) result(value)
    double complex, intent(in) :: z
    double precision :: a,b
    value=unavailable()
    if (ieee_is_nan(real(z)).or.ieee_is_nan(aimag(z))) return
    if (.not.finite_complex(z)) then
       value=ieee_value(0.d0,ieee_positive_inf); return
    endif
    a=max(abs(real(z)),abs(aimag(z))); b=min(abs(real(z)),abs(aimag(z)))
    if (a==0.d0) then
       value=ieee_value(0.d0,ieee_negative_inf)
    else
       value=log10(a)+0.5d0*log10(1.d0+(b/a)**2)
    endif
  end function

  double precision function complex_magnitude(z) result(value)
    double complex, intent(in) :: z
    double precision :: a,b,factor
    value=complex_log_magnitude(z)
    if (.not.finite_complex(z)) return
    a=max(abs(real(z)),abs(aimag(z))); b=min(abs(real(z)),abs(aimag(z)))
    value=0.d0
    if (a==0.d0) return
    factor=sqrt(1.d0+(b/a)**2)
    if (a>huge(a)/factor) then
       value=ieee_value(0.d0,ieee_positive_inf)
    else
       value=a*factor
    endif
  end function

  subroutine chebyshev_pair(a,z,polynomial,derivative_z,ok)
    double precision, intent(in) :: a(0:)
    double complex, intent(in) :: z
    double complex, intent(out) :: polynomial,derivative_z
    logical, intent(out) :: ok
    double complex :: t0,t1,t2,b1,b2,b0,d1,d2,d0
    double precision :: cap
    integer :: j
    ok=.false.
    polynomial=cmplx(unavailable(),unavailable(),kind(1.d0)); derivative_z=polynomial
    ! A conservative bound leaves ample headroom for three-factor recurrences.
    cap=huge(1.d0)**0.125d0/100.d0
    if (.not.finite_complex(z).or.any(.not.ieee_is_finite(a))) return
    if (complex_magnitude(z)>cap.or.maxval(abs(a))>cap) return
    ! Match the production polynomial's forward accumulation. Compute the
    ! derivative independently with differentiated Clenshaw (d/dz).
    t0=(1.d0,0.d0); t1=z; polynomial=a(0)*t0
    if (size(a)>1) polynomial=polynomial+a(1)*t1
    do j=2,size(a)-1
       if (max(complex_magnitude(t0),complex_magnitude(t1),complex_magnitude(polynomial))>cap) return
       t2=2.d0*z*t1-t0; polynomial=polynomial+a(j)*t2
       t0=t1; t1=t2
    enddo
    b1=(0.d0,0.d0); b2=b1; d1=b1; d2=b1
    do j=size(a)-1,1,-1
       if (max(complex_magnitude(b1),complex_magnitude(b2),complex_magnitude(d1), &
            complex_magnitude(d2))>cap) return
       b0=a(j)+2.d0*z*b1-b2
       d0=2.d0*b1+2.d0*z*d1-d2
       b2=b1; b1=b0; d2=d1; d1=d0
    enddo
    if (max(complex_magnitude(b1),complex_magnitude(d1),complex_magnitude(d2))>cap) return
    derivative_z=b1+z*d1-d2
    ok=finite_complex(polynomial).and.finite_complex(derivative_z)
  end subroutine

  double precision function safe_ratio(numerator,denominator) result(value)
    double precision, intent(in) :: numerator,denominator
    value=unavailable()
    if (.not.ieee_is_finite(numerator).or..not.ieee_is_finite(denominator)) return
    if (denominator<=0.d0.or.numerator<0.d0) return
    if (numerator==0.d0) then
       value=0.d0
    elseif (log10(numerator)-log10(denominator)>log10(huge(1.d0))-1.d-12) then
       value=ieee_value(0.d0,ieee_positive_inf)
    else
       value=numerator/denominator
    endif
  end function

  subroutine sample_log_chebyshev(a,z,norm2,limit,sample)
    double precision, intent(in) :: a(0:),norm2,limit
    double complex, intent(in) :: z
    type(fitted_sample), intent(out) :: sample
    double precision :: phase,magnitude,ceiling
    logical :: ok,cutoff
    sample%flags=0; sample%analytic_available=.false.
    sample%value=cmplx(unavailable(),unavailable(),kind(1.d0)); sample%analytic=sample%value
    sample%log_function=unavailable(); sample%log_derivative=unavailable()
    if (complex_magnitude(z)>1.d0) sample%flags=sample_domain
    call chebyshev_pair(a,z,sample%polynomial,sample%derivative_z,ok)
    if (.not.ok.or..not.ieee_is_finite(norm2).or.norm2==0.d0) then
       sample%flags=ior(sample%flags,sample_invalid); return
    endif
    sample%log_function=real(sample%polynomial)
    sample%log_derivative=real(sample%polynomial)+log10(log(10.d0))+ &
         complex_log_magnitude(sample%derivative_z)-log10(abs(norm2))
    ! Production applies the domain test before forming or limiting P.
    if (iand(sample%flags,sample_domain)==0.and.real(sample%polynomial)>limit) &
         sample%flags=ior(sample%flags,sample_limited)
    if (real(sample%polynomial)<log10(tiny(1.d0))) sample%flags=ior(sample%flags,sample_small)
    ceiling=log10(huge(1.d0))-1.d-12
    if (real(sample%polynomial)>ceiling.or.sample%log_derivative>ceiling) &
         sample%flags=ior(sample%flags,sample_overflow)
    cutoff=iand(sample%flags,sample_limited+sample_domain)/=0
    if (cutoff) then
       ! Do not evaluate an exponential the production evaluator suppresses.
       sample%value=(0.d0,0.d0); return
    endif
    if (real(sample%polynomial)>ceiling.or.abs(aimag(sample%polynomial))>huge(1.d0)/log(10.d0)) return
    sample%value=10.d0**sample%polynomial
    if (.not.finite_complex(sample%value)) sample%flags=ior(sample%flags,sample_invalid)
    if (sample%value==(0.d0,0.d0)) sample%flags=ior(sample%flags,sample_underflow)
    if (sample%log_derivative>ceiling) return
    if (sample%derivative_z==(0.d0,0.d0)) then
       sample%analytic=(0.d0,0.d0)
    else
       phase=aimag(sample%polynomial)*log(10.d0)+ &
            atan2(aimag(sample%derivative_z),real(sample%derivative_z))
       if (norm2<0.d0) phase=phase+acos(-1.d0)
       magnitude=10.d0**sample%log_derivative
       sample%analytic=cmplx(magnitude*cos(phase),magnitude*sin(phase),kind(1.d0))
    endif
    sample%analytic_available=finite_complex(sample%analytic)
  end subroutine

  double complex function guarded_product(a,b) result(value)
    double complex, intent(in) :: a,b
    value=cmplx(unavailable(),unavailable(),kind(1.d0))
    if (.not.finite_complex(a).or..not.finite_complex(b)) return
    if (complex_log_magnitude(a)+complex_log_magnitude(b)>log10(huge(1.d0))-2.d0) return
    value=a*b
  end function

  double complex function derivative_difference(plus,minus,h) result(value)
    double complex, intent(in) :: plus,minus
    double precision, intent(in) :: h
    value=cmplx(unavailable(),unavailable(),kind(1.d0))
    if (.not.finite_complex(plus).or..not.finite_complex(minus).or.h<=0.d0) return
    if (max(complex_magnitude(plus),complex_magnitude(minus))>huge(1.d0)/4.d0) return
    if (complex_log_magnitude(plus-minus)-log10(2.d0*h)>log10(huge(1.d0))-2.d0) return
    value=(plus-minus)/(2.d0*h)
  end function

  subroutine difference_metrics(analytic,fd,function_scale,grid_step,step,metrics,resolved)
    double complex, intent(in) :: analytic,fd
    double precision, intent(in) :: function_scale,grid_step,step
    double precision, intent(out) :: metrics(5)
    logical, intent(out) :: resolved
    double precision :: amplitude_a,amplitude_fd,scale,floor
    metrics=unavailable(); resolved=.false.
    if (.not.finite_complex(analytic).or..not.finite_complex(fd)) return
    amplitude_a=complex_magnitude(analytic); amplitude_fd=complex_magnitude(fd)
    if (max(amplitude_a,amplitude_fd)>huge(1.d0)/4.d0) return
    metrics(1)=complex_magnitude(fd-analytic)
    floor=max(tiny(1.d0),sqrt(epsilon(1.d0))*safe_ratio(function_scale,grid_step))
    scale=max(amplitude_a,floor)
    metrics(2)=safe_ratio(metrics(1),scale)
    if (amplitude_a>safe_ratio(epsilon(1.d0)*function_scale,grid_step)) &
         metrics(3)=safe_ratio(amplitude_fd,amplitude_a)
    resolved=min(amplitude_a,amplitude_fd)>safe_ratio(epsilon(1.d0)*function_scale,step)
    if (resolved) metrics(4)=atan2(sin(atan2(aimag(fd),real(fd))-atan2(aimag(analytic),real(analytic))), &
         cos(atan2(aimag(fd),real(fd))-atan2(aimag(analytic),real(analytic))))
    metrics(5)=safe_ratio(epsilon(1.d0)*function_scale,step)
  end subroutine

  subroutine compare_brackets(om,k,mass,charge,n,pperp,pres,dfpar,dfperp,analytic,b,terms,ok)
    double complex, intent(in) :: om,pres,dfpar,dfperp,analytic
    double precision, intent(in) :: k,mass,charge,pperp
    integer, intent(in) :: n
    double complex, intent(out) :: b(4),terms(6)
    logical, intent(out) :: ok
    double precision :: cap
    b=cmplx(unavailable(),unavailable(),kind(1.d0)); terms=b(1); ok=.false.
    cap=sqrt(huge(1.d0))/100.d0
    if (.not.finite_complex(dfpar).or..not.finite_complex(dfperp).or..not.finite_complex(pres)) return
    if (.not.finite_complex(om).or.mass==0.d0) return
    if (abs(mass)<1.d0/(cap**0.25d0)) return
    if (max(complex_magnitude(dfpar),complex_magnitude(dfperp),complex_magnitude(pres), &
         complex_magnitude(om),abs(k),abs(pperp),abs(dble(n)*charge),1.d0/abs(mass))>cap**0.25d0) return
    ! Keep the source's multiply-then-divide grouping for the original bracket.
    b(1)=((pperp*dfpar-pres*dfperp)*k/(mass)+om*dfperp)
    b(2)=(k*pperp/mass)*dfpar+(dble(n)*charge/mass)*dfperp
    terms(1)=(pperp*dfpar)*k/mass
    terms(2)=(-pres*dfperp)*k/mass
    terms(3)=om*dfperp
    terms(4)=terms(2)+terms(3)
    terms(5)=(om-k*pres/mass)*dfperp
    terms(6)=(dble(n)*charge/mass)*dfperp
    if (finite_complex(analytic)) then
       if (complex_magnitude(analytic)<=cap**0.25d0) then
          b(3)=((pperp*analytic-pres*dfperp)*k/(mass)+om*dfperp)
          b(4)=(k*pperp/mass)*analytic+(dble(n)*charge/mass)*dfperp
       endif
    endif
    ok=finite_complex(b(1)).and.finite_complex(b(2))
  end subroutine
end module alps_diagnostic_math
