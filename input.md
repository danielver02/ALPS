# ALPS Input

This is a reference for the key input parameters used by ALPS.

## Namelists in input files.

The following namelists and associated input parameters are read in by ALPS from the input file.

### *&system*  
General system parameters.

**`kperp`**  
Initial perpendicular wavevector \\(k_{\perp} d_{p}\\).

**`kpar`**  
Initial parallel wavevector $k_{\parallel} d_{p}$.

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

**`ngamma`**  
Relativistic momentum space resolution, $N_{\Gamma}$.

**`npparbar`**  
Relativistic parallel momentum space resolution, $N_{\bar{p}_{\parallel}}$.

**`vA`**  
Reference Alfven velocity, normalized to the speed of light, $v_{A}/c$.

**`arrayName`**  
Name of input array, located in 'distribution' folder.

**`Bessel_zero`**  
Maximum amplitude of Bessel function to determine `nmax`.

**`numiter`**  
Maximum number of iterations in secant method.

**`D_threshold`**  
Minimum threshold for secant method.

**`D_prec`**  
Size of bounding region for secant method.

**`D_gap`**  
Size of allowable difference between roots.

**`positions_principal`**  
Number of parallel momentum steps distant from the resonant momentum
included in the numerical calculation of Eqn 3.5, $M_{I}$.

**`n_resonance_interval`**  
How many steps should be used to integrate around the resonance,
$M_{P}$, used for integrating near poles (see section 3.1).

**`Tlim`**  
Threshold for analytical principal-value integration, $t_{\mathrm{lim}}$.

**`maxsteps_fit=500`**  
Maximum number of fitting iterations.

**`lambda_initial_fit`**  
Inital Levenberg-Marquardt damping parameter.

**`lambdafac_fit`**  
Adjustment factor for Levenberg-Marquardt damping parameter.

**`epsilon_fit`**  
Convergence for Levenberg-Marquardt fit.

**`fit_check`**  
If true, output fitted functions for each species to file in distribution directory.

**`determine_minima`**  
If true, after map search, determine minima and refine solutions.

**`scan_option`**  
Select case for wavevector scans:

- 1: Consecutive scans along input paths in wavevector space,  
- 2: Double scan over wavevector plane.

**`n_scan`**  
Number of wavevector scans.  
0 turns off wavevector scans.  
Must be 1 or larger for `scan_option`=1.  
Must be set to 2 for `scan_option`=2.  


### *&guess_m*  
Initial guess of complex frequency for $m$th solution.  
Only used when `use_map`=.false.  
Need to have number of name lists equal to `nroots`.

**`g_om`**  
Guess for real solution $\omega_{r}/\Omega_{p} $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_{p} $.


### *&maps_1*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{r}/\Omega_{p}$ map search.
Spacing automatically calculated between `omi` and `omf`.  

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{p}$ map search.
Spacing automatically calculated between `gami` and `gamf`  

**`omi`**  
Smallest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`omf`**  
Largest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_{p}$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_{p}$ value for complex map search.

**`ni`**  
Number of $\gamma/\Omega_{p}$ points in frequency grid.

**`nr`**  
Number of $\omega_{r}/\Omega_{p}$ points in frequency grid.


### *&spec_j*  
Species parameters list for distribution $f_{j}$.

**`nn`**  
Relative density $n_{j}/n_{p}$.

**`qq`**  
Relative charge $q_{j}/q_{p}$.

**`mm`**  
Relative mass $m_{j}/m_{p}$.

**`ff`**  
Number of fitted functions for analytical continuation calculation.

**`relat`**  
Treat $f_{j}$ as non-relativistic or relativistic.

**`log_fit`**  
Use linear or $\log_{10}$ fitting routine.

**`use_bM`**  
Use actual numerical integration (F) or bi-Maxwellian/cold-plasma proxy via NHDS routines,
with parameters read in from &bM_spec_j namelist.

**`AC_method`**  
Choose the method for the evaluation of the analytic continuation:

- 0: Use the function that is defined analytically in distribution/distribution_analyt.f90
- 1: Use the fit routine as defined in the &ffit_j_k namelist.
- 2: Use a polynomial basis representation as defined in the &poly_spec_j namelist. This method should only be used if $|\gamma|\ll |\omega_{r}|$.


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

$$F_{\kappa}(\hat{p}\_{\perp},\hat{p}\_{\parallel})=
u_1\mathrm{exp}[-y \hat{p}\_{\perp}]
\mathrm{exp}[-u_{2}*(\hat{p}\_{\parallel}+u_{3})^2]
.$$

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
Relative drift of bi-Maxwellian distribution $f_{j}$ or the cold plasma species in units of $m_{p} v_{A,p}$.


### *&poly_spec_j*
Input for the polynomial representation of the input distribution for the analytical continuation.
Only used if `AC_method=2`.

**`kind`**  
Type of the basis polynomial:

- 1: Chebychev

**`order`**  
Maximum order of the basis polynomial.

**`poly_log_max`**  
When using logfit for the polynomial representation, set all output values to zero if the log(fit_function_poly) is greater than this variable.


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
Scan variable to define end of scan through wavevector space (only for `scan_type=1`).

**`swf`**  
Scan variable to define end of scan through wavevector space.

**`swlog`**  
Use $\log_{10}$ (T) or linear (F) spacing.

**`ns`**  
Number of output scan values.

**`nres`**  
Resolution between output scan values.

**`heating`**  
Calculates heating rates if true.

**`eigen`**  
Calculates eigenfunctions if true.     
