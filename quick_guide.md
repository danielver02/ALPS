# Quick Guide

lorem $\sqrt{3x-1}+(1+x)^2$ ipsum

$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$

This is a quick guide that gives a reference for the key parameters used by ALPS.

## Namelists in input files.

The following namelists and associated input parameters are read in by ALPS from the input file.

### *&system*  
General system parameters.

**`kperp`**  
Initial perpendicular wavevector $k_\perp d_p$.

**`kpar`**  
Initial parallel wavevector $k_\parallel d_p$.

**`nspec`**   
Number of plasma species.

**`nroots`**   
Number of dispersion solutions to find and follow.

**`use_map`**   
Choice of:  
    -True: Searching for roots over a map in complex frequency space [see &maps_1 namelist].  
    -False: Input `nroots` guesses for solutions [see &guess_* namelist].

**`writeOut`**  
Write or suppress output to screen.

**`nperp`**  
Perpendicular momentum space resolution, $N_\perp$.

**`npar`**  
Parallel momentum space resolution, $N_\parallel$.

**`ngamma`**  
Relativistic momentum space resolution, $N_\Gamma$.

**`npparbar`**  
Relativistic parallel momentum space resolution, $N_{\bar{p}_\parallel}$.

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
included in the numerical calculation of Eqn 3.5, $M_I$.

**`n_resonance_interval`**  
How many steps should be used to integrate around the resonance,
$M_P$, used for integrating near poles (see section 3.1).

**`Tlim`**  
Threshold for analytical principal-value integration, $t_\textrm{lim}$.

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
Select case for scans;  
  -1 Consecutive scans along input paths in wavevector space,  
  -2 Double scan over wavevector plane.
  
**`n_scan`**  
Number of wavevector scans.  
Must be set to 2 for scan_option=2;  
Must be 1 or larger for scan_option=1.  
0 turns off wavevector scans.


### *&guess_m*  
Initial guess of complex frequency for $m^{\textrm{th}}$ solution.  
Only used when `use_map`=.false.  
Need to have number of name lists equal to `nroots`.

**`g_om`**  
Guess for real solution $\omega_{\textrm{r}/\Omega_p $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_p $.


### *&maps_1*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{\textrm{r}}/\Omega_p$ map search.

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_p$ map search.

**`omi`**  
Smallest $\omega_{\textrm{r}}/\Omega_p$ value for complex map search.

**`omf`**  
Largest $\omega_{\textrm{r}}/\Omega_p$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_p$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_p$ value for complex map search.

**`ni`**  
Number of $\gamma/\Omega_p$ points in frequency grid.

**`nr`**  
Number of $\omega_{\textrm{r}}/\Omega_p$ points in frequency grid.


### *&spec_j*  
Species parameters list for species j.

**`nn`**  
Relative density for $f_j$, $n_j/n_p$.

**`qq`**  
Relative charge for $f_j$, $q_j/q_p$.

**`mm`**  
Relative mass for $f_j$, $m_j/m_p$.

**`ff`**  
Number of fitted functions for $f_j$.

**`relat`**  
Treat species as non-relativistic or relativistic.

**`log_fit`**  
Use linear or $\log_{10}$ fitting routine.

**`use_bM`**  
Use actual numerical integration or bi-Maxwellian proxy via NHDS,
with parameters read in from &bM_spec_j namelist.


### *&ffit_j_k*
Initial Fit Values for species $j$, function $k$.

**`fit_type_in`**  
Kind of fit function:

0. Analytic function (KGK: Need to add details on this functionality).  

1. Maxwellian,  
$$F_M(\hat{p}_\parallel)=
u_1 \mathrm{exp}[-y \hat{p}^2_\perp-u_2(\hat{p}_\parallel-u_3)^2 ]. $$

2. Kappa,  
$$F_{\kappa}(\hat{p}_\parallel)= u_1[1+u_2(\hat{p}_\parallel-u_3)^2+y \hat{p}_\perp^2]^{u_4}.$$

3. Juettner with $p_\perp,p_\parallel$,  
$$F_{J}(\hat{p}_\perp,\hat{p}_\parallel)=
u_1 \mathrm{exp}[-u_2].$$

4. Juettner with variable $\Gamma$, constant $\bar{p}_\parallel$,  
$$F_{J}(\Gamma)= u_1 \mathrm{exp}[-y \Gamma].$$

5. Juettner with $p_\perp,p_\parallel$; variable $\bar{p}_\parallel$,  
$$F_{\kappa}()= .$$

6. Bi-Moyal distribution  
$$F_{\kappa}()= .$$

**`fit_1`-`fit_5`**  
Fit parameters, $u_1$-$u_5$, defined in the above equations for each of the types of fit functions.
Not all parameters will be used for all functions.  
Suggested values for parameters generated by generate_distribution.

**`perpcorr`**  
This parameter, $y$ in Eqn. B1, compensates for the strong
$p_\perp$ dependence of $u_1$, making the fit more reliable.


### *&bM_spec_j*
Bi-Maxwellian parameters; for species j.
Only used if use_bM=T.

**`bM_nmaxs`**
Maximum number of resonances to consider.

**`bM_Bessel`**
Precision threshold for $I_n$.

**`bM_betas`**
$\beta_{\parallel,j}$ of biMaxwellian distribution $f_j$.

**`bM_alphas`**
$T_{\perp,j}/T_{\parallel,j} $ of biMaxwellian distribution $f_j$.

**`bM_pdrifts`**
Relative drift of biMaxwellian distribution $f_j$,
in units of $m_p v_{A,p}$.


### *&scan_input_l*
Inputs for scanning parameter space for l$^{\textrm{th}$ scan.  

**`scan_type`**
Type of parameter scan;  
0: Current value of $\textbf{k}$ to
 $k_\perp$=range$_\textrm{i}$ and $k_\parallel $=range$_\textrm{f}$.   
1: $\theta_0 \rightarrow \theta_1$ at fixed $|k|$
 from current value of $\theta=\mathrm{atan}(k_\perp/k_\parallel)$
 to range$_\textrm{f}$.  
2: Wavevector scan at fixed angle $\theta_{k,B}$ to $|k|$=range$_\textrm{f}$.  
3: $k_\perp$ scan with constant $k_\parallel$.  
4: $k_\parallel$ scan with constant $k_\perp$.  

**`swi`**
Initial scan value.

**`swf`**  
Final scan value.

**`swlog`**
Use $\mathrm{log}_{10}$ (T) or linear (F) spacing.

**`ns`**
Number of output scan values.

**`nres`**
Resolution between output scan values.

**`heating`**
Calculates heating rates if true.

**`eigen`**  =.false.
Calculates eigenfunctions if true.     




