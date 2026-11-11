!Copyright (c) 2026, Kristopher G. Klein and Daniel Verscharen
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
!I                              Version 1.2                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (daniel.verscharen@rub.de)                              I
!I                                                                             I
!===============================================================================

module alps_distribution_analyt
  !! This module evaluates a pre-defined function for f0. This is an alternative
  !! to the use of an f0 table or a bi-Maxwellian approximation.
  !! This function can be used to define f0 for the integration or for the analytic
  !! continuation (or both).
  !! It works in the non-relativistic case as well as the relativistic case.
  implicit none

  public :: distribution_analyt

contains

double complex function distribution_analyt(is,pperp,ppar,call_ID,mass_gd,vA_gd)
!! This function returns the pre-defined function as f0.
  implicit none

  integer, intent(in) :: is
  !! Index of species.

  integer, intent(in) :: call_ID
  !! ID to describe where this function is called from. 1 = analytic continuation, 2 = generate_distribution.

  double precision, intent(in) :: mass_gd, vA_gd
  !! Mass and Alfven speed when calling from generate_distribution.

  double precision :: pperp, gamma_rel
  !! Perpendicular momentum and Lorentz factor.

  double complex :: ppar, pparbar
  !! Parallel momentum (non-relativistic and relativistic).

  double complex :: f0
  !! Value of the distribution function.

  double precision :: pi
  !! Pi.

  double precision :: beta
  !! Plasma beta.

  double precision :: ms
  !! Mass of species.

  ! Define parameters:
  pi=atan(1.d0)*4.d0

  if (call_ID.EQ.1) then
    gamma_rel=pperp
    pparbar=ppar
  else if (call_ID.EQ.2) then
        gamma_rel=sqrt(1.d0+(pperp**2+real(ppar)**2)*vA_gd*vA_gd/(mass_gd*mass_gd))
		    pparbar = ppar*vA_gd/mass_gd
  end if
 

  ! This function defines the f0 table for the creation of distributions and for the analytic continuation.
  ! Ensure that the distribution is normalised. generate_distribution can do this automatically, but the fit
  ! routine of ALPS will not do this.
  ! Define the function for however many species you would like to do this.

  ! Remember that ppar and pparbar are a complex variables. For most functions, this should not make a difference.
  ! In the following, define the complex f0 function of (pperp, ppar) or (gamma_rel,pparbar) for any number of species as needed.
  ! Only those species that are set to have either distributionS=0 in the generate_distribution input file or
  ! those species that have ff=0 in the ALPS input file will be treated in this way.
  ! Any other cases for species will be ignored.

  ! For relativistic cases: When using distribution_analyt for the analytic continuation, the function must be defined in terms of (gamma_rel, pparbar). When using distribution_analyt to generate a new distribution via generate_distribution, the function can be defined either in terms of (pperp, ppar) or in terms of (gamma_rel, pparbar). 


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
