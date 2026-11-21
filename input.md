# ALPS Input

This is a reference for the key input parameters used by ALPS.

## Namelists in execution input files.

The following namelists and associated input parameters are read in by ALPS from the input file.

### *&system*  
General system parameters.

**`profile_startup`**
Optional logical in `&system`, default `.false.`. 
When enabled, rank 0 reports elapsed wall time in seconds using `MPI_Wtime`.
This includes 

- VDF input, 
- original-grid processing, 
- relativistic preparation (if applicable), 
- spline matrix/solve/evaluation, 
- postprocessing,
- diagnostic output, 
- analytic-continuation fitting, and total 
- VDF startup before the first dispersion calculation. 

Species timings are summed. 
Detailed reporting is independent of `writeOut`.

**`kperp`**  
Initial perpendicular wavevector $k_{\perp} d_{ref}$.

**`kpar`**  
Initial parallel wavevector $k_{\parallel} d_{ref}$.

**`nspec`**   
Number of plasma species.

**`nroots`**   
Number of dispersion solutions to find and follow.

**`use_map`**   
Choice of:  

- True: Searching for roots over a map in complex frequency space (see &maps_1 namelist).  
- False: Input `nroots` guesses for solutions (see &guess_1 namelist).

**`writeOut`**  
Write or suppress output to screen.

**`nperp`**  
Perpendicular momentum space resolution, $N_{\perp}$.
The input file must have $N_{\perp}+1$ values spanning parallel momentum space.

**`npar`**  
Parallel momentum space resolution, $N_{\parallel}$.
The input file must have $N_{\parallel}+1$ values spanning parallel momentum space.

**`vA`**  
Alfven velocity using the mass density of the reference species, 
normalized to the speed of light, $v_{A,ref}/c$.

**`arrayName`**  
Name of input array, located in 'distribution' folder.

**`Bessel_zero`**  
Maximum amplitude of Bessel function to determine `nmax`.

**`secant_method`**  
Selection for root finding method.

- 0: *secant method from NHDS.* Uses full secant steps and declares convergence when the **magnitude of the dispersion determinant** falls below `D_threshold`.
- 1: *rtsec method from PLUME.* Applies half of each calculated secant step and stops when the calculated **step magnitude** falls below `D_tol`, or the determinant is exactly zero.
- 2: *Secant method with oscillation control (default):* Adds adaptive damping, caps steps at 10% of the current frequency magnitude, and switches to a regularized finite-difference Newton step after repeated detections of revisited frequencies; convergence uses `D_threshold`.

All three use `D_prec` to construct a second initial frequency and `numiter` to limit iterations.

**`numiter`**  
Maximum number of iterations in secant method.

**`D_threshold`**  
Minimum threshold for secant methods 0 and 2.

**`D_prec`**  
Size of bounding region for secant method.

**`D_gap`**  
Size of allowable difference between roots.

**`D_tol`**  
Tolerance for secant method = 1, rtsec.

**`positions_principal`**  
Number of parallel momentum steps distant from the resonant momentum
included in the numerical calculation of Verscharen et al 2018, Eqn 3.5, $M_{I}$.

**`n_resonance_interval`**  
How many steps should be used to integrate around the resonance,
$M_{P}$, used for integrating near poles (see section 3.1 of Verscharen et al 2018).

**`Tlim`**  
Threshold for analytical principal-value integration, $t_{\mathrm{lim}}$.

**`determine_minima`**  
If true, after map search, determine minima and refine solutions.

**`scan_option`**  
Select operation (default 1):

- 1: Consecutive scans along input paths in wavevector space,  
- 2: Double scan over wavevector plane.
- 3: Nyquist calculation at the input `kperp,kpar`; bypass root solving and scans.
- 4: Nyquist calculations along the existing `scan_input_N` wavevector paths.
- 5: Nyquist wavevector-plane scan using exactly two different `scan_input_N` directions.

For Nyquist operation details, see [nyquist.md](nyquist.md).

**`n_scan`**  
Number of wavevector scans.  
0 turns off wavevector scans.  
Must be 1 or larger for `scan_option`=1.  
Must be set to 2 for `scan_option`=2.  
For `scan_option=3`, `n_scan`, `nroots`, and `use_map` are ignored; no guesses,
maps, or scan instructions are required. Options 1 and 2 retain their existing
root-search/map behavior, including when `n_scan=0`.
Option 4 requires `n_scan>=1`; it ignores `nroots,use_map`, reads species and
`scan_input_1` through `scan_input_n_scan`, and bypasses guesses, maps and root tracking.
Option 5 requires `n_scan=2`; it likewise bypasses guesses/maps/root tracking.
The two scan types must differ and both must be in 1–4 (type 0 is unsupported).

### *&nyquist*

Namelist for Nyquist operations (`scan_option=3`, `4`, or `5`).
The parameters and defaults are listed below. For contour methods, growth-rate
bracketing, scan geometry, output formats, and limitations, see
[Nyquist Stability Analysis](nyquist.md).

Option 3 also accepts nonrelativistic numerical VDFs (`use_bM=F`, `relat=F`)
and mixed analytic/numerical species for initial C1 validation. Use the normal
`distribution/<arrayName>.<species>.array` input and continuation settings;
`bM_*` parameters apply only to analytic species. Numerical species require a
finite, nonnegative, uniform rectangular grid (`nperp,npar>=4`; strictly positive
samples for methods 1/2 when `log_fit=T`). Density/mass must be positive and charge nonzero for
every species. Relativistic and cold analytic Nyquist responses remain unsupported.
Options 4/5 and maximum-growth searches remain restricted to warm analytic
bi-Maxwellians. For numerical option 3 set `nyq_find_max_growth=F`.
The numerical winding describes the computed boundary trajectory; physical root
counting requires momentum-grid, support and integration-continuity validation.

```fortran
&nyquist
  nyq_geometry='rectangle'
  nyq_omega_min=-1, nyq_omega_max=1
  nyq_gamma_min=0.01, nyq_gamma_max=0.2
  nyq_sampling%initial_intervals=8
  nyq_sampling%max_evaluations=200000
  nyq_save_samples=F
/
```

| Parameter | Meaning / default |
|---|---|
| `nyq_geometry` | `'rectangle'` (default) or `'horizontal'`, lowercase. |
| `nyq_omega_min`, `nyq_omega_max` | Required finite ordered real-frequency bounds, in reference gyrofrequency units. |
| `nyq_gamma_min` | Finite strictly positive imaginary frequency; default `1.d-4`. |
| `nyq_gamma_max` | Required only for rectangles; finite and greater than `nyq_gamma_min`. |
| `nyq_sampling` | The A1 `nyq_options` object; set individual components with `%` syntax below. |
| `nyq_save_samples` | Default `F`; opt into ordered sample storage and ASCII output. |
| `nyq_normalize_rectangle` | Default `F`; diagnostic normalized rectangle for comparison with horizontal phases. Applies only to rectangles. |
| `nyq_find_max_growth` | Default `F`; enable bounded maximum-growth estimation with cached closed strips. Rectangles only. |
| `nyq_growth_tolerance` | Default `1.d-4`; finite positive absolute bracket-width tolerance for `'linear'`, dimensionless relative bracket-width tolerance for `'log'`. Both use strict inequality. |
| `nyq_growth_max_iterations` | Default 32; positive integer limit on bisection steps. |
| `nyq_growth_spacing` | Default `'linear'`; `'log'` enables geometric trials and nearest-original-boundary strips. Used/validated only when maximum-growth estimation is enabled. |

The `nyq_sampling` components use A1 defaults directly; there is no second set
of sampling controls or implementation:

| Component | Default |
|---|---:|
| `initial_intervals` (per segment) | 8 |
| `max_depth` (per initial interval) | 24 |
| `max_evaluations` (logical requests over all passes of each contour) | 200000 |
| `max_passes` | 6 (must be at least 3) |
| `phase_limit` (radians) | 0.7 |
| `curvature_limit` | 0.15 |
| `phase_tolerance` (absolute radians across passes) | `1e-8` |
| `integer_tolerance` (closed winding units) | `1e-8` |
| `log_magnitude_limit` | 1.5 |
| `local_ratio` | 0.05 |

For `scan_option=3`–`5`, warm analytic bi-Maxwellians are the validated
Nyquist backend in this source version; nonrelativistic numerical-VDF and cold
cases require separate validation. A rectangle requires `nyq_gamma_max`;
`nyq_find_max_growth=T` requires a rectangle. The `solution/` directory must exist.

### *&fitting_params*

Optional global Levenberg–Marquardt (LM) settings shared by all numerical species using `AC_method=1`, 
for both nonrelativistic and relativistic cases. 
Omission of the full list or individual keys retains the defaults below.

```fortran
&fitting_params
    maxsteps_fit       = 500
    lambda_initial_fit = 1.d0
    lambdafac_fit      = 1.d1
    epsilon_fit        = 1.d-8
/
```

**`maxsteps_fit`** (integer): 
maximum LM fitting iterations.

**`lambda_initial_fit`** (double precision): 
initial LM damping parameter.

**`lambdafac_fit`** (double precision): 
damping multiplier when a trial increases the sum of squared residuals, and divisor otherwise.

**`epsilon_fit`** (double precision): 
convergence threshold for the absolute change in the sum of squared residuals, not the residual itself.

LM parameter files are written during fitting; fit-result diagnostics are
controlled by `writeOut`.

Methods 0 and 2 do not require this name list. 
Method 2 uses `&poly_spec_j`; mixed method-1
and method-2 species use their respective configurations independently.

### *&relativistic_params*

Settings shared by all relativistic species. 
ALPS reads this group when any `&spec_*` block sets `relat=T`, 
after reading all species and before sizing fit arrays. 

```fortran
&relativistic_params
  ngamma = 500
  npparbar = 500
  rel_cache = 'auto'
  rel_cache_dir = 'distribution/'
/
```

The cache is a reusable snapshot of the completed relativistic interpolation. 
ALPS still reads the original ASCII VDF, but with rel_cache='auto' it first checks whether it has already calculated the transformed grids, normalized distribution, and derivatives for those inputs. 
These arrays are saved in a binary file, which preserves every floating-point bit without rounding values into ASCII text. 
The long filename is an automatically generated fingerprint of the VDF contents, interpolation settings, and computational build. 
Matching inputs can reuse the same file across runs or species; changing relevant inputs produces a different fingerprint. 
Before loading, ALPS checks compatibility and a checksum to detect damage. 
If the cache is missing or invalid, ALPS simply calculates the interpolation again and saves a fresh copy. 
The original VDF files remain the authoritative inputs, and fitting and dispersion calculations continue normally.

An absent or empty group uses `ngamma=100`, `npparbar=200`, `rel_cache='off'`, and
`rel_cache_dir='distribution/'`. 
Partial groups retain defaults for omitted fields.
When no species is relativistic, the entire group is ignored, including validation.

**`ngamma`**
Relativistic momentum space resolution, $N_{\Gamma}$.
Requires `ngamma >= 2` when a relativistic species is present.

**`npparbar`**
Relativistic parallel momentum space resolution, $N_{\bar{p}_{\parallel}}$.
Requires `npparbar >= 4` to keep the existing three-point resonance stencil in the interior.

**`rel_cache`**
Optional `&relativistic_params` character parameter, default `'off'`. Modes are case-insensitive:

- `'off'`: construct the relativistic representation normally, without cache identity or file access.
- `'auto'`: load a fully validated entry, or construct and save one on a miss.
- `'refresh'`: construct every relativistic species and atomically replace its entry.

**`rel_cache_dir`**
Optional `&relativistic_params` character parameter, default `'distribution/'`.
The directory must already exist.
Only rank 0 accesses entries.
Missing/unwritable storage, incompatible entries, and integrity failures fall back to normal calculation. Invalid modes are input errors.
Existing relativistic diagnostic files are still written, including on cache hits and when `writeOut=.false.`.
Analytic-continuation fitting is always performed normally.

Entries are identified by numerical VDF contents, coordinates, construction settings, and strict executable/library compatibility, independently of filenames and species
indices.

### *&guess_m*  
Initial guess of complex frequency for $m$th solution.  
Only used when `use_map`=.false.  
Need to have number of name lists equal to `nroots`.

**`g_om`**  
Guess for real solution $\omega_{r}/\Omega_{ref} $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_{ref} $.


### *&maps_1*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{r}/\Omega_{ref}$ map search.
Spacing automatically calculated between `omi` and `omf`.  

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{ref}$ map search.
Spacing automatically calculated between `gami` and `gamf`  

**`omi`**  
Smallest $\omega_{r}/\Omega_{ref}$ value for complex map search.

**`omf`**  
Largest $\omega_{r}/\Omega_{ref}$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_{ref}$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_{ref}$ value for complex map search.

**`ni`**  
Number of $\gamma/\Omega_{p}$ points in frequency grid.

**`nr`**  
Number of $\omega_{r}/\Omega_{ref}$ points in frequency grid.


### *&spec_j*  
Species parameters list for distribution $f_{j}$.
The first species is set as the reference.

**`nn`**  
Relative density $n_{j}/n_{ref}$.

**`qq`**  
Relative charge $q_{j}/q_{ref}$.

**`mm`**  
Relative mass $m_{j}/m_{ref}$.

**`ff`**  
Number of fitted functions for analytical continuation calculation, 
when using `AC_method`=1.

**`relat`**  
Treat $f_{j}$ as non-relativistic or relativistic. 
Analytic bimaxwellian species (`use_bM=T`) require `relat=F`.

**`log_fit`**  
Use linear or $\log_{10}$ fitting routine. 
Relativistic numerical Chebyshev
continuation (`relat=T`, `AC_method=2`) supports only `log_fit=T`.

**`use_bM`**  
Use actual numerical integration (F) or 
bi-Maxwellian/cold-plasma proxy via NHDS routines (T),
with parameters read in from &bM_spec_j namelist.

**`AC_method`**  
Choose the method for the evaluation of the analytic continuation:

- 0: Use the function that is defined analytically in distribution/distribution_analyt.f90
- 1: Use the fit routine as defined in the &ffit_j_k namelist.
- 2: Use a polynomial basis representation as defined in the &poly_spec_j namelist. 
*The limits of this analytic continuation for 
$|\gamma|\gtsrim |\omega_{r}|$ are being investigated*.


Relativistic numerical method 2 uses native SVD/Clenshaw models at fixed Lorentz factor, with `kind=1`, positive requested `order`, and `log_fit=T`. 
Original source samples must be finite and strictly positive on a finite, increasing
rectangular momentum grid.
Zeros are rejected before the log-TPS transformation,
including when caching is enabled. 
No positive floor is substituted.

These models represent the **completed transformed distribution** on the existing Gamma/parallel-momentum grid. Its global natural-log thin-plate spline extends
beyond the supplied source momentum rectangle at some cone-admissible nodes.
Those values are TPS extensions, not source observations. ALPS reports source and
completed bounds, interior/exterior counts, affected Gamma rows and fresh/cache
origin once per selected species. Counts are recomputed from the actual geometry;
cache hits verify identity against the current source arrays and coordinates.
Fit status, degree, rank and residual do not certify source support or accurate
complex continuation. Numerical range or unavailable-row failures remain explicit.
Original-grid fit files and fitted-moment totals omit these fixed-Gamma models;
original VDF moments are unchanged.


### *&ffit_j_k*
Initial Fit Values for species $j$, function $k$.

**`fit_type_in`**  
Kind of fit function:

- 1: Maxwellian,  

$$F_M(\hat{p}\_{\parallel})=u_{1}\mathrm{exp}[-y{\hat{p}}^2\_{\perp}-u_{2}(\hat{p}\_{\parallel}-u_{3})^2]$$

- 2: Kappa,  

$$F_{\kappa}(\hat{p}\_{\parallel})=u_{1}[1+u_2({\hat{p}}\_{\parallel}-u_{3})^2+yu_{5} {\hat{p}}^2\_{\perp}]^{u_{4}}.$$

- 3: Juettner with $p_{\perp},p_{\parallel}$,  

$$F_{J}(\hat{p}\_{\perp},\hat{p}\_{\parallel})=
u_{1}\mathrm{exp}\left[-u_{2}\sqrt{1+\frac{\hat{p}^2\_{\perp}+(\hat{p}^2\_{\parallel}-u_3)^2 v_A^2}{m_{j}^2 c^2}}\right].$$

- 4: Juettner with variable $\Gamma$, constant $\bar{p}_{\parallel}$,  

$$F_{J}(\Gamma)= u_{1} \mathrm{exp}[-y \Gamma].$$

- 5: Juettner with $p_{\perp},p_{\parallel}$; variable $\bar{p}_{\parallel}$,  

$$F_{\kappa}(\hat{p}\_{\perp},\hat{p}\_{\parallel})=u_{1}\mathrm{exp}[-y \hat{p}\_{\perp}]\mathrm{exp}[-u_{2}*(\hat{p}\_{\parallel}-u_{3})^2].$$

- 6: Bi-Moyal distribution

$$F_{bMo}(\hat{p}\_{\perp},\hat{p}\_{\parallel})= u_{1} \mathrm{exp}[0.5 (y u_4 \hat{p}^2\_{\perp} + u_{2} (\hat{p}\_{\parallel} -u_{3})^2 -\mathrm{exp}(y u_{4} \hat{p}^2\_{\perp} + u_{2} (\hat{p}\_{\parallel}-u_{3})^2) )].$$

**`fit_1`-`fit_5`**  
Fit parameters, $u_{1}$ - $u_{5}$, defined in the above equations for each of the types of fit functions.
Not all parameters will be used for all functions.  
Suggested values for parameters generated by generate_distribution.

**`perpcorr`**  
This parameter, $y$ in Eqn. B1, compensates for the strong
$p_{\perp}$ dependence of $u_1$, making the fit more reliable.


### *&bM_spec_j*
Bi-Maxwellian/cold-plasma parameters; for species j.
Only used if `use_bM=T`.

**`bM_nmaxs`**  
Maximum number of resonances to consider.

**`bM_Bessel`**  
Precision threshold for $I_n$.

**`bM_betas`**  
$\beta_{\parallel,j}$ of bi-Maxwellian distribution $f_{j}$. If this variable is set to 0.d0, then the code will treat the given species with the susceptibility from cold-plasma theory.

**`bM_alphas`**  
$T_{\perp,j}/T_{\parallel,j}$ of bi-Maxwellian distribution $f_{j}$.

**`bM_pdrifts`**  
Relative drift of bi-Maxwellian distribution $f_{j}$ or the cold plasma species in units of $m_{ref} v_{A,ref}$.


### *&poly_spec_j*
Input for the polynomial representation of the input distribution for the analytical continuation.
Only used if `AC_method=2`.

**`kind`**  
Type of the basis polynomial:

- 1: Chebychev

**`order`**  
Maximum order of the basis polynomial.

Native LOG10 continuation reconstructs the physical distribution using the shared
Chebyshev kernel. Overflow, underflow and invalid evaluations have explicit
numerical statuses; there is no adjustable logarithmic suppression threshold.
The former `log_max` namelist setting is unsupported and must be removed from inputs.

**`derivative_type`**  
Controls the **nonrelativistic** parallel derivative in the Landau integral
when using polynomial continuation. 
Relativistic AC method 2 always uses the shared
analytic derivative with respect to `u=P_parallel/(m_s*c)` at fixed Gamma; its
Gamma finite-difference stencil remains unchanged.

- 0: Finite Difference Derivative
- 1: Analytic Derivative of Polynomial Basis Function (default)

Omitting `derivative_type` in `&poly_spec_j` selects `1` for that species.
Set `derivative_type=0` explicitly to use the finite-difference derivative.

### *&scan_input_l*
Inputs for scanning parameter space for $l$th scan.  

**`scan_type`**  
Type of parameter scan:

- 0: Current value of $\textbf{k}$ to $k\_{\perp}$=`swi` and $k\_{\parallel}$ =`swf`.   
- 1: $\theta_0 \rightarrow \theta_1$ at fixed $|k|$ from current value of $\theta=\mathrm{atan}(k\_{\perp}/k\_{\parallel})$ to `swf`.  
- 2: Wavevector scan at fixed angle $\theta_{k,B}$ to $|k|$ =`swf`.  
- 3: $k\_{\perp}$ scan with constant $k\_{\parallel}$ to $k\_{\perp}$=`swf`.  
- 4: $k\_{\parallel}$ scan with constant $k\_{\perp}$ to $k\_{\parallel}$=`swf`.  

**`swi`**  
Target perpendicular component for `scan_type=0`; unused for other scan types.

**`swf`**  
Scan variable to define end of scan through wavevector space.

**`swlog`**  
Use $\log_{10}$ (T) or linear (F) spacing.

**`ns`**  
Number of output steps after the initial wavevector. An ordinary segment file
has `ns+1` rows; option 4 writes shared segment boundaries only once.

**`nres`**  
Root-tracking steps per output interval. Option 4 evaluates only output points,
using the same coordinate indices without the intermediate tracking work.

**`heating`**  
Calculates heating rates if true.

**`eigen`**  
Calculates eigenfunctions if true.     

## Namelists in distribution input files.

The following namelists and associated input parameters are read in by the routine `generate_distribution` from input files. This generates `nspec` arrays of distributions that are read in by ALPS

### *&system*
General system parameters.

**`nspec`**

Number of plasma species.

**`beta`**

Reference plasma beta.

**`vA`**

Reference Alfven velocity, normalized to the speed of light, $v_{A,ref}/c$.

**`nperp`**  

Perpendicular momentum space resolution, $N_{\perp}$.
The input file will have $N_{\perp}+1$ values spanning parallel momentum space.
Set the same values of `nperp` in both the distribution input file and ALPS execution input file.

**`npar`**  

Parallel momentum space resolution, $N_{\parallel}$.
The input file will have $N_{\parallel}+1$ values spanning parallel momentum space.
Set the same values of `npar` in both the distribution input file and ALPS execution input file.

**`maxP`**

Maximum value of reference Alfven momentum, $m_{ref} v_{A,ref}$.

**`writeName`**

Name of output arrays.

### *&spec_j*

**`ms_read`**

Mass ratio $m_j/m_{ref}$.

**`taus`**

Temperature ratio $T_{\parallel,j}/T_{\parallel,ref}$.

**`alphs`**

Temperature anisotropy $T_{\perp,j}/T_{\parallel,j}$.

**`ps`**

Normalized drift momentum $m_j v_{drift}/m_{ref} v_{A,ref}$.

**`kappas`**

$\kappa$ index.
Only used if `distribution`=2.

**`distributions`**

Type of distribution:

- 0: Use distribution from distribution_analyt.
- 1: Bi-Maxwellian distribution.
- 2: Bi-Kappa distribution.
- 3: Anisotropic Juettner distribution.
- 4: Bi-Moyal distribution.

**`autoscaleS`**

If `autoscaleS`=.true., define maximum momentum for each array automatically using global value `maxP`.

**`maxPperpS`**

If `autoscaleS`=.false., force maximum perpendicular momentum for species array.

**`maxPparS`**

If `autoscaleS`=.false., force maximum parallel momentum for species array.




