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

module alps_var
  !!Contains all global variables.
  use mpi
  implicit none


  private

  !I/O Variables:
  character(500) :: runname
  !!Root of input file name.
  
  character(500) :: foldername
  !!Directory of input file name.
  
  integer :: nroots
  !!Number of dispersion solutions under consideration.

  integer :: nroots_max
  !!Number of dispersion solutions found in frequency map scan.
  
  logical :: use_map
  !!Choice of:
  !! -T: searching for roots over a map in complex frequency space;
  !! -F: input (nroots) guesses for solutions.
  
  logical :: writeOut =.true.
  !! Write or suppress output to screen.

  integer :: unit_error
  !! Output unit for error file.
  
  !MPI variables:
  integer :: nproc
  !! Total number of processors.
  
  integer :: iproc
  !! Number of local processor.
  
  logical :: proc0
  !! T if iproc=0.
  
  integer :: ierror
  !! Integer error flag.

  !Plasma Parameters:
  double precision :: kperp
  !! Perpendicular wavenumber, normalized by inverse reference inertial length, \(k_\perp d_p\).
  
  double precision :: kpar
  !! Parallel wavenumber, normalized by inverse reference inertial length, \(k_\parallel d_p\).

  double precision :: vA
  !! Alfven Velocity, normalized to speed of light, \(v_{Ap}/c\).
  
  double precision :: Bessel_zero=1.d-45
  !! Calculated Bessel functions until the maximum is less than this value.

  integer :: nspec
  !!Number of plasma components.
  
  integer :: nspec_rel
  !!Number of relativistic plasma components.


  double complex, dimension(:), allocatable :: wroots !(1:numroots)
  !!Dispersion solutions, ranging from 1 to numroots.

  integer :: numroots = 100
  !! Maximum number of solutions.
  
  double precision :: omi
  !!Smallest \(\omega_{\textrm{r}}/\Omega_p\) value for complex map search.
  
  double precision :: omf
  !!Largest \(\omega_{\textrm{r}}/\Omega_p\) value for complex map search.

  double precision :: gami
  !!Smallest \(\gamma/\Omega_p\) value for complex map search.

  double precision :: gamf  
  !!Largest \(\gamma/\Omega_p\) value for complex map search.
  
  logical :: loggridw
  !!Linear (F) or Log (T) spacing for \(\omega_{\textrm{r}}/\Omega_p\) map search.

  logical :: loggridg
  !!Linear (F) or Log (T) spacing for \(\gamma/\Omega_p\) map search.

  integer :: nr=128
  !!Number of points in \(\omega_{\textrm{r}}/\Omega_p\) grid.

  integer :: ni=128
  !!Number of points in \(\gamma/\Omega_p\) grid.

  integer :: nperp
  !!Number of perpendicular momentum space grid points, \(N_\perp\).

  integer :: npar
  !!Number of parallel momentum space grid points, \(N_\parallel\).
  
  integer :: ngamma=100
  !!Number of grid points in relativitic \(\Gamma=\sqrt{1+\frac{p_\perp^2+p_\parallel^2}{m_j^2c^2}}\),
  !!\(N_\Gamma\) (Eqn. 3.14).

  integer :: npparbar=200
  !!Number of grid points in dimensionless paralell momentum
  !!\(\bar{p}_\parallel = p_\parallel/m_j c\), \(N_{\bar{p}_\parallel}\).

  integer :: positions_principal=5
  !! How many parallel momentum steps distant from the resonant momentum
  !! are included in the numerical calculation of Eqn 3.5, \(M_I\).

  double precision :: Tlim=0.01d0
  !! Threshold for analytical principal-value integration for
  !! evaluating Eqn 3.6 and 3.7, \(t_{\textrm{lim}}\).

  integer :: sproc
  !! The species number on which this process is working.

  integer :: numiter=50
  !! Maximum number of iterations in secant method.

  double precision :: D_threshold=1.d-5
  !! Minimum threshold for secant method.

  double precision :: D_prec=1.d-5
  !! Size of bounding region for secant method.

  double precision :: D_gap =1.d-5
  !! Size of allowable difference between roots.

  double precision :: pi
  !! The ratio of a circle's circumference to its diameter.

  character (75) :: arrayName
  !! Name of input files for distributions.

  double precision, dimension(:,:,:), allocatable :: f0
  !! Background distribution function array \(f_{0j}\);
  !! (1:nspec,0:nperp,0:npar).

  double precision, dimension(:,:,:), allocatable :: f0_rel
  !! Relativistic background distribution function array \(f_{0j}\);
  !! (1:nspec,0:ngamma,0:npparbar)
    
  double precision, dimension(:,:,:,:), allocatable :: df0
  !!Perpendicular and parallel derivatives of \(f_{0j}\);
  !!(1:nspec,0:nperp,0:npar,1:2),
  !!with \(\partial_{p_\perp} f_{0j} \) index 1
  !!with \(\partial_{p_\parallel} f_{0j} \) index 2

  double precision, dimension(:,:,:,:), allocatable :: df0_rel
  !!Derivatives of \(f_{0j}\);
  !!(1:nspec,0:nperp,0:npar,1:2),
  !!with \(\partial_{\Gamma} f_{0j} \) index 1
  !!and \(\partial_{\bar{p}_\parallel} f_{0j} \) index 2.
  
  double precision, dimension(:,:,:,:), allocatable :: pp
  !!Momentum Space Array for \(f_{0j}\);
  !!(1:nspec,0:nperp,0:npar,1:2)
  !!with \(p_\parallel/m_p v_A\) index 1
  !!and \(p_\perp/m_p v_A\) index 2.

  double precision, dimension(:,:,:), allocatable :: gamma_rel
  !!Relativistic momentum space array of \(\Gamma\);
  !!(1:nspec,0:ngamma,0:npparbar).

  double precision, dimension(:,:,:), allocatable :: pparbar_rel
  !!Relativistic momentum space array of \(\bar{p}_\parallel\);
  !!(1:nspec,0:ngamma,0:npparbar).
  
  integer, dimension(:), allocatable :: nmax
  !!number of n values to sum over, (1:nspec).

  integer :: nlim(2)
  !!Lower and Upper limits for n values for iproc to sum over.

  
  double precision, dimension(:), allocatable :: ns
  !!Ratio of species density to reference \(n_j/n_p\), (1:nspec).

  double precision, dimension(:), allocatable :: qs
  !!Ratio of species charge to reference \(q_j/q_p\), (1:nspec).

  double precision, dimension(:), allocatable :: ms
  !!Ratio of species mass to reference \(m_j/m_p\), (1:nspec).
    
  logical, dimension(:), allocatable :: relativistic
  !!Use relativistic treatment; (1:nspec).

  double complex, dimension(:,:), allocatable :: wave
  !!Wave Equation Tensor (1:3,1:3).
  
  double complex, dimension(:,:,:), allocatable :: chi0
  !!Susceptibility Tensor (1:nspec,1:3,1:3).

  double precision, dimension(:,:), allocatable :: bessel_array
  !! Array of Bessel functions; (nlim(1):nlim(2)+1,0:nperp).

  !Fit Parameters for Hybrid-Analytical Continuation Method:
  integer, dimension(:), allocatable :: n_fits
  !!Number of fitted functions, (1:nspec)

  integer, dimension(:,:), allocatable :: fit_type
  !! Type of analytic function to be fit, (1:nspec,1:maxval(nfits));
  !! 1) Maxwellian,
  !! 2) Kappa,
  !! 3) Juettner with \(p_\perp,p_\parallel\),
  !! 4) Juettner with \(\Gamma,\bar{p}_\parallel\), constant \(\bar{p}_\parallel\),
  !! 5) Juettner with \(p_\perp,p_\parallel\); variable \(\bar{p}_\parallel\),
  !! 6) Bi-Moyal distribution.
  
  integer :: maxsteps_fit=500
  !!Maximum number of fitting iterations.
  
  double precision :: lambda_initial_fit=1.d0
  !!Inital Levenberg-Marquardt damping parameter.
  
  double precision :: lambdafac_fit=1.d1
  !!Adjustment factor for Levenberg-Marquardt damping parameter.
  
  double precision :: epsilon_fit=1.d-8
  !!Convergence for Levenberg-Marquardt fit.

  !Fit output:
  double precision, dimension(:,:,:,:), allocatable :: param_fit
  !!Fit parameters, (1:nspec,0:nperp,4,maxval(n_fits)).
  
  double precision, dimension(:,:), allocatable :: perp_correction
  !!This parameter, \(y\) in Eqn. B1, compensates for the strong
  !! \(p_\perp\) dependence of \(u_1\), making the fit more reliable.
  
  logical :: fit_check=.true.
  !!If true, output fitted functions to ASCII file for each species.

  logical :: determine_minima=.true.
  !!If true, after map search, determine minima and refine solutions.

  integer :: n_resonance_interval=100
  !! Mow many steps should be used to integrate around the resonance,
  !!\(M_P\), used for integrating near poles (see section 3.1).

  integer :: scan_option=1
  !!!Select case for scans;
  !!1) consecutive scans along input paths in wavevector space
  !!2) double scans of two selected parameters

  integer :: n_scan=0 !number of wavevector scans.
  !must be set to 2 for scan_option=2
  !must be 1 or larger for scan_option=1
  !n_scan=0 turns off wavevector scans

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
  public :: runname, foldername, writeOut
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
  public :: relativistic, logfit, usebM
  public :: f0_rel,df0_rel,nspec_rel,gamma_rel,pparbar_rel,ngamma,npparbar
  public :: bMnmaxs,bMBessel_zeros,bMbetas,bMalphas,bMpdrifts

 end module alps_var
