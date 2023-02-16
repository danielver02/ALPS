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

module alps_var
  use mpi
  implicit none


  private

  !I/O Variables
  character(500) :: runname   !Root of input file name
  integer :: option           !System Option: Chooses Style of Run
                              !See README file for list of options
  integer :: nroots           !Number of Dispersion Solutions
                              !under consideration
  integer :: nroots_max       !Number of Dispersion Solutions
                              !found in frequency map scan
  logical :: use_map          !Choice of
    !T: searching for roots over a map in complex frequency space
    !OR
    !F: input (nroots) guesses for solutions
  logical :: writeOut =.true.   !Write or Spress output to screen

  integer :: unit_error

  !MPI variables
  integer :: nproc           !Number of Processors (1 to nproc)
  integer :: iproc           !Number of local processor (0 to nproc-1)
  logical :: proc0           !T if iproc=0
  integer :: ierror          !Integer error flag

  !Plasma Parameters
  !Wave vector normalized by inverse reference inertial length
  double precision :: kperp  !Perpendicular Wave number;
  double precision :: kpar   !Parallel Wave number;

  double precision :: vA     !reference Alfven Velocity /c
  double precision :: Bessel_zero=1.d-45	! Use Bessel functions until
  					! the maximum is less than this value
  !The proton typically serves as the reference species.

  integer :: nspec           !Number of plasma species
  integer :: nspec_rel		 !Number of relativistic plasma species

  !Solutions for dispersion surface
  double complex, dimension(:), allocatable :: wroots !(1:numroots)
  !maximum # of roots
  integer :: numroots = 100
  !Limits on (real, imaginary) map search
  double precision :: omi,omf, gami,gamf
  !T: log spacing for complex frequency map search
  !F: linear spacing for complex frequency map search
  logical :: loggridw,loggridg
  !number of points in (real, imaginary) frequency space grid
  integer :: ni=128, nr=128

  !Number of Momentum Space Grid Points
  integer :: nperp, npar
  integer :: ngamma=100, npparbar=200

  ! How many points in ipar should be used to calculate the principal value integral?
  integer :: positions_principal=5

  ! Threshold for analytical principal-value integration
  double precision :: Tlim=0.01d0

  ! The species number on which this process is working
  integer :: sproc

  ! Maximum number of iterations in secant method
  integer :: numiter=50

  ! Value of "zero" for secant method
  double precision :: D_threshold=1.d-5

  ! size of bounding region for secant method
  double precision :: D_prec=1.d-5

  ! size of allowable difference between roots
  double precision :: D_gap =1.d-5

  !pi
  double precision :: pi

  !Name of input files for distributions
  character (75) :: arrayName

  !Background Distribution Function Array
  !f0(1:nspec,0:nperp,0:npar)
  double precision, dimension(:,:,:), allocatable :: f0
  double precision, dimension(:,:,:), allocatable :: f0_rel

  !Perpendicular and parallel Derivatives of f0
  !df0(1:nspec,0:nperp,0:npar,1:2)
  !with index 1-> dvperp
  !     index 2-> dvpar
  double precision, dimension(:,:,:,:), allocatable :: df0
  double precision, dimension(:,:,:,:), allocatable :: df0_rel

  !Momentum Space Array for f0
  !pp(1:nspec,0:nperp,0:npar,1:2)
  !with index 1->vperp
  !     index 2->vpar
  double precision, dimension(:,:,:,:), allocatable :: pp

  double precision, dimension(:,:,:), allocatable :: gamma_rel,pparbar_rel

  !number of n values to sum over
  !nmax(1:nspec)
  integer, dimension(:), allocatable :: nmax

  !Lower and Upper limits for n values for
  !iproc to sum over
  integer :: nlim(2)

  !Ratios of species density
  double precision, dimension(:), allocatable :: ns
  !charge
  double precision, dimension(:), allocatable :: qs
  !and mass
  double precision, dimension(:), allocatable :: ms
  !relative to the reference value.
  logical, dimension(:), allocatable :: relativistic !use relativistic treatment; set for each species

  !Wave Equation Tensor
  !wave(1:3, 1:3)
  double complex, dimension(:,:), allocatable :: wave
  double complex, dimension(:,:,:), allocatable :: chi0

  ! Array of Bessel functions:
  ! bessel_array(nlim(1):nlim(2)+1,0:nperp)
  double precision, dimension(:,:), allocatable :: bessel_array

  !Fit Parameters for Hybrid-Analytical Continuation Method
  integer, dimension(:), allocatable :: n_fits !number of fitted functions, per species
  integer, dimension(:,:), allocatable :: fit_type !type of analytic function to be fit
  !fit_type(1:nspec,1:maxval(nfits))
  integer :: maxsteps_fit=500 !maximum number of fitting iterations
  double precision :: lambda_initial_fit=1.d0 !Inital Levenberg-Marquardt damping parameter
  double precision :: lambdafac_fit=1.d1 !Adjustment factor for Levenberg-Marquardt damping parameter
  double precision :: epsilon_fit=1.d-8 !Convergence for fit
  !Fit output
  !param_fit(1:nspec,nperp,4,maxval(n_fits))
  double precision, dimension(:,:,:,:), allocatable :: param_fit
  double precision, dimension(:,:), allocatable :: perp_correction
  logical :: fit_check=.true. !T-> output fitted functions to ASCII file for each species

  logical :: determine_minima=.true. ! after map search, determine minima and refine?
  integer :: n_resonance_interval=100 ! how many steps should be used to integrate around the resonance

  integer :: scan_option=1 !select case for scans.
  !1: consecutive scans along input paths in wavevector space
  !2: double scans of two selected parameters

  integer :: n_scan=0 !number of wavevector scans.
  !must be set to 2 for scan_option=2
  !must be 1 or larger for scan_option=1
  !n_scan=0 turns off wavevector scans

  logical :: use_secant=.true. !use secant or rtsec methods

  logical, dimension(:), allocatable :: logfit ! Use logarithmic fitting?


  logical, dimension(:), allocatable :: usebM ! Use bi-Maxwellian calculation from NHDS?


  integer, dimension(:), allocatable ::  bMnmaxs ! Maximum number of n for bi-Maxwellian calculation from NHDS

  double precision, dimension(:), allocatable ::  bMBessel_zeros ! Besser-zero for bi-Maxwellian calculation from NHDS

  double precision, dimension(:), allocatable ::  bMbetas ! Species beta for bi-Maxwellian calculation from NHDS

  double precision, dimension(:), allocatable ::  bMalphas ! Species temperature anisotropy for bi-Maxwellian calculation from NHDS

  double precision, dimension(:), allocatable ::  bMpdrifts ! Species drift momentum for bi-Maxwellian calculation from NHDS


  public :: scanner
     type :: scanner
        double precision :: range_i       !initial value
        double precision :: range_f       !final value
        logical :: log_scan   !T-> log, F-> linear scan
        logical :: heat_s     !T-> heating calc; F-> no heating
        logical :: eigen_s    !T-> eigen calc;   F-> no eigen
        integer :: type_s     !Type of parameter scan
        integer :: n_out     !Number of steps
        integer :: n_res      !scan resolution
        double precision :: diff, diff2 !step size for scanned parameter
     end type scanner
     !-=-=-=-=-=-=-=-=-
     !Defines nature of parameter scans:
     !     Type: 0 k_0-> k_1
     !           1 theta_0 -> theta_1
     !           2 k_fixed angle
     !     Type: 3 kperp
     !           4 kpar
     !-=-=-=-=-=-=-=-=-
  type (scanner), dimension (:), allocatable :: scan

  double precision :: kperp_last, kpar_last
  double precision :: kperp_0, kpar_0

  public :: nproc, iproc, proc0, ierror
  public :: runname, option, writeOut
  public :: kperp, kpar, nspec, use_map, wroots
  public :: loggridw, loggridg, omi, omf, gami, gamf
  public :: arrayName, nperp, npar, f0, pp, df0, bessel_array
  public :: nmax, nlim, wave, numiter, D_threshold, D_prec, D_gap, chi0
  public :: ns, qs, ms, vA, pi, Bessel_zero,sproc
  public :: ni, nr, positions_principal, Tlim
  public :: n_fits, maxsteps_fit,  lambda_initial_fit, lambdafac_fit, epsilon_fit
  public :: param_fit, fit_check, fit_type, perp_correction
  public :: nroots, nroots_max, numroots
  public :: determine_minima, n_resonance_interval
  public :: unit_error, scan_option, n_scan, scan
  public :: kperp_last, kpar_last, kperp_0, kpar_0
  public :: use_secant, relativistic, logfit, usebM
  public :: f0_rel,df0_rel,nspec_rel,gamma_rel,pparbar_rel,ngamma,npparbar
  public :: bMnmaxs,bMBessel_zeros,bMbetas,bMalphas,bMpdrifts

 end module alps_var
