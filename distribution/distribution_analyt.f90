!===============================================================================
!I                                                                             I
!I                              A  L  P  S                                     I
!I                     Arbitrary Linear Plasma Solver                          I
!I                                                                             I
!I                              Version 1.0                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)                                I
!I                                                                             I
!===============================================================================

module alps_distribution_analyt
  implicit none

  public :: distribution_analyt

contains

double complex function distribution_analyt(is,pperp,ppar)
  implicit none
  integer :: is
  double precision :: pperp
  double complex :: ppar
  double complex :: f0


  ! This function defines the f0 table for the creation of distributions and for the analytic continuation.
  ! Ensure that the distribution is normalised. generate_distribution can do this automatically, but the fit
  ! routine of ALPS will not do this.
  ! Define the function for however many species you would like to do this.

  ! Remember that ppar is a complex variable. For most functions, this should not make a difference.
  ! In the following, define the complex f0 function of pperp and ppar for any number of species as needed.
  ! Only those species that are set to have either distributionS=0 in the generate_distribution input file or
  ! those species that have ff=0 in the ALPS input file will be treated in this way.
  ! Any other cases for species will be ignored.

  select case(is)

    case(1) ! Species 1
      f0=1.d0


    case(2) ! Species 2
      f0=1.d0


    end select


distribution_analyt=f0

end function

end module
