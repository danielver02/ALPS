! Test-only phase-space densities. C is explicit; no normalization is implicit.
module relativistic_references
  use iso_fortran_env, only : dp=>real64
  implicit none
  private
  public :: dp,juttner,anisotropic,drifting,to_relativistic,from_relativistic
contains
  pure subroutine juttner(g,u,theta,c,f,fg,fu)
    real(dp), intent(in) :: g,theta,c
    complex(dp), intent(in) :: u
    complex(dp), intent(out) :: f,fg,fu
    f=cmplx(c*exp(-g/theta),0._dp,dp)
    fg=-f/theta; fu=cmplx(0._dp,0._dp,dp)
  end subroutine
  pure subroutine anisotropic(g,u,a,b,c,f,fg,fu)
    real(dp), intent(in) :: g,a,b,c
    complex(dp), intent(in) :: u
    complex(dp), intent(out) :: f,fg,fu
    f=c*exp(-a*g-b*u*u)
    fg=-a*f; fu=-2*b*u*f
  end subroutine
  pure subroutine drifting(g,u,theta,beta,c,f,fg,fu)
    real(dp), intent(in) :: g,theta,beta,c
    complex(dp), intent(in) :: u
    complex(dp), intent(out) :: f,fg,fu
    real(dp) :: gd
    gd=1/sqrt(1-beta*beta)
    f=c*exp(-gd*(g-beta*u)/theta)
    fg=-gd*f/theta; fu=gd*beta*f/theta
  end subroutine
  pure subroutine to_relativistic(pperp,ppar,va,mass,g,u)
    real(dp), intent(in) :: pperp,ppar,va,mass
    real(dp), intent(out) :: g,u
    u=va*ppar/mass
    g=sqrt(1+(va*pperp/mass)**2+u*u)
  end subroutine
  pure subroutine from_relativistic(g,u,va,mass,pperp,ppar)
    real(dp), intent(in) :: g,u,va,mass
    real(dp), intent(out) :: pperp,ppar
    pperp=mass/va*sqrt(g*g-1-u*u)
    ppar=mass/va*u
  end subroutine
end module
