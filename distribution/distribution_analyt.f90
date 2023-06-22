!Copyright (c) 2023, Kristopher G. Klein and Daniel Verscharen
!All rights reserved.
!
!This source code is licensed under the BSD-style license found in the
!LICENSE file in the root directory of this source tree.
!
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
  !! This module evaluates a pre-defined function for f0. This is an alternative
  !! to the use of an f0 table or a bi-Maxwellian approximation.
  !! This function can be used to define f0 for the integration or for the analytic
  !! continuation (or both).
  implicit none

  public :: distribution_analyt

contains

double complex function distribution_analyt(is,pperp,ppar)
!! This function returns the pre-defined function as f0.
  implicit none

  integer, intent(in) :: is
  !! Index of species.

  double precision :: pperp
  !! Perpendicular momentum.

  double complex :: ppar
  !! Parallel momentum.

  double complex :: f0
  !! Value of the distribution function.

  double precision :: pi
  !! Pi.

  double precision :: beta
  !! Plasma beta.

  double precision :: ms
  !! Mass of species.

  pi=atan(1.d0)*4.d0


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

      ! The example below illustrates how to set up a Maxwellian with beta = 1:
      beta=1.d0
      ms=1.d0

      f0=(pi**(-1.5d0) /((ms * beta )**(3.d0/2.d0) )) * exp( -( ppar**2/( beta * ms)&
        + (pperp**2)/( beta * ms ) ) )

    case(2) ! Species 2

      ! The example below illustrates how to set up a Maxwellian with beta = 1:
      beta=1.d0
      ms=1.d0/1836.d0

      f0=(pi**(-1.5d0) /((ms * beta )**(3.d0/2.d0) )) * exp( -( ppar**2/( beta * ms)&
        + (pperp**2)/( beta * ms ) ) )

    end select


distribution_analyt=f0

end function

end module
