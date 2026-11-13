! Test-only analytic CONTINUATION control. This is not an independent dielectric
! solver or a complex interpolator for arbitrary transformed tables.
module alps_distribution_analyt
  implicit none
  private
  double precision :: amplitude(2)=1.d0,quadratic=0.d0
  public :: distribution_analyt,set_susceptibility_reference
contains
  subroutine set_susceptibility_reference(c,b)
    double precision, intent(in) :: c(2),b
    amplitude=c; quadratic=b
  end subroutine
  double complex function distribution_analyt(is,gamma,u,call_ID,mass_gd,vA_gd) result(f)
    integer, intent(in) :: is,call_ID
    double precision, intent(in) :: gamma,mass_gd,vA_gd
    double complex, intent(in) :: u
    if(call_ID/=1) error stop 'analytic control accepts transformed continuation coordinates only'
    f=amplitude(is)*exp(-2.d0*gamma-quadratic*u*u)
  end function
end module
