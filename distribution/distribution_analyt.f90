!=============================================================================!
!=============================================================================!
!!*ALPS                                                                  *!!
!!                                                                           !!
!!Kristopher Klein, Daniel Verscharen                                        !!
!!
!!Space Science Center, University of New Hampshire                          !!
!!                                                                           !!
!!*INPUT/OUTPUT FUNCTIONS                                                   *!!
!=============================================================================!
!=============================================================================!
!
module alps_distribution_analyt
  implicit none

  public :: distribution_analyt

contains

double complex function distribution_analyt(is,pperp,ppar_C)
  implicit none
  integer :: is
  double precision :: pperp
  double complex :: ppar_C
  double complex :: f0


  ! This function defines the f0 table for the creation of distributions and for the analytic continuation.
  ! Ensure that the distribution is normalised. generate_distribution can do this automatically, but the fit
  ! routine of ALPS will not do this.
  ! Define the function for however many species you would like to do this.
  select case(is)
    case(1)
      f0=123.d0
    case(2)
      f0=321.d0
    end select


distribution_analyt=f0

end function

end module
