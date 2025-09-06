# ALPS Output

ALPS writes output solutions to the `/solution` directory.  
All output file names start with the name of the input file used in running ALPS, e.g. 
```
mpirun -np 4 ./src/ALPS filename.in
```
will produce output files all starting with the string *filename*.

## *filename*.map

Value of the dispersion tensor $\mathcal{D}(\omega_{\textrm{r}},\gamma)$ on a defined complex frequency grid.  
Solutions to the dispersion relation satisfy $|\mathcal{D}|  =0$.
This file is generated from the *map_search* subroutine in ALPS_fns.f90, and invoked when `use_map` =.true. .  

The data is ordered in columns as  
1. $\omega_r$  
2. $\gamma$   
3. $\log_{10} |\mathcal{D}|$  
4. Re $[|\mathcal{D}|]$  
5. Im $[|\mathcal{D}|]$  

The *&maps_1* namelist in *filename*.in determines the structure of *filename*.map.  
The range of $\omega_{\textrm{r}}/\Omega_{ref}$ is from `omi` to `omi` with `nr` steps. Logorithmic or linear spacing is selected with `loggridw`.
The range of $\gamma_{\textrm{r}}/\Omega_{ref}$ is from `gami` to `gami` with `ni` steps. Logorithmic or linear spacing is selected with `loggridg`.

## *filename*.roots

Identified solutions to the dispersion relation $|\mathcal{D}|  =0$, calculated using *refine_guess* in ALPS_fns.f90 when `determine_minima` is set to true.  

The data is ordered as  
1. Solution number
2. $\omega_r/\Omega_{ref}$  
3. $\gamma/Omega_{ref}$   
4. $\log_{10} |\mathcal{D}|$  
5. Re $[|\mathcal{D}|]$  
6. Im $[|\mathcal{D}|]$  

The routine uses  either the coarse dispersion tensor map generated from the *map_search* subroutine (in the case of `use_map` = .true.)  
or from the input guesses (for `use_map` = .false.).  
Only the first `nroots` solutions will be identified and written to file.

## *filename*.scan_*scan_type_l*.root_m

The complex frequencies associated with solution *m* calculated from `om_scan`  
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_{ref}$
2. $k_\parallel d_{ref}$  
3. $\omega_{\textrm{r}}/\Omega_{ref}$   
4. $\gamma/\Omega_{ref}$   

See the *&scan_input* namelist description in the Quick Guide for details on determining the kind of wavevector scan.  
This same data structure is preserved for the output from `om_double_scan`.

## *filename*.eigen_*scan_type_l*.root_m

The eigenfunctions associated with solution *m* calculated from `om_scan` when `eigen` is set to .true.
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_{ref}$
2. $k_\parallel d_{ref}$  
3. $\omega_{\textrm{r}}/\Omega_{ref}$   
4. $\gamma/\Omega_{ref}$   
5. Re $[E_x]$ 
6. Im $[E_x]$ 
7. Re $[E_y]$ 
8. Im $[E_y]$ 
9. Re $[E_z]$ 
10. Im $[E_z]$ 
11. Re $[B_x]$ 
12. Im $[B_x]$ 
13. Re $[B_y]$ 
14. Im $[B_y]$ 
15. Re $[B_z]$ 
16. Im $[B_z]$  
17. [+6(is-1)] Re $[\delta U_{x,is}]$   
18. [+6(is-1)] Im $[\delta U_{x,is}]$   
19. [+6(is-1)] Re $[\delta U_{y,is}]$   
20. [+6(is-1)] Im $[\delta U_{y,is}]$   
21. [+6(is-1)] Re $[\delta U_{z,is}]$   
22. [+6(is-1)] Im $[\delta U_{z,is}]$   
17. [+6(`nspec`)+2(is-1)] Re $[\delta n_{is}]$   
18. [+6(`nspec`)+2(is-1)] Im $[\delta n_{is}]$ 

This same data structure is preserved for the output from `om_double_scan`.

## *filename*.heat_*scan_type_l*.root_m

The heating rates associated with solution *m* calculated from `om_scan` when `heating` is set to .true.
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_{ref}$
2. $k_\parallel d_{ref}$  
3. $\omega_{\textrm{r}}/\Omega_{ref}$   
4. $\gamma/\Omega_{ref}$   
5. [+(is-1)] $\gamma_{is}/\omega_{\textrm{r}}$

This same data structure is preserved for the output from `om_double_scan`.

## *filename*.heat_mech*scan_type_l*.root_m

The heating rates associated with Landau damping, Transit Time damping, and the $n=\pm 1$ cyclotron resonance are computed following [Huang et al 2024 JPP](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174) and expressed in Appendix C of [Klein & Verscharen 2025 PoP](https://doi.org/10.1063/5.0286477). The values associated with solution *m* are calculated in the routine `om_scan` when `heating` is set to .true. in the *&scan_input_l* namelist.

The Transit time damping terms given by the sum of $\gamma_{is}^{yy}/\omega_{\textrm{r}}$ and $\gamma_{is}^{yz}/\omega_{\textrm{r}}$.
The Landau damping terms by the sum of $\gamma_{is}^{zy}/\omega_{\textrm{r}}$ and $\gamma_{is}^{zz}/\omega_{\textrm{r}}$.
Unlike Huang et al 2024 JPP, we separate the n=-1 and n=+1 terms.
These calculations are accurate in weak damping limit $\gamma < \omega_{\textrm{r}}$.


The data is ordered as  
1. $k_\perp d_{ref}$
2. $k_\parallel d_{ref}$  
3. $\omega_{\textrm{r}}/\Omega_{ref}$   
4. $\gamma/\Omega_{ref}$   
5. [+4(is-1)] $\gamma_{is}^{TTD}/\omega_{\textrm{r}}$
6. [+4(is-1)] $\gamma_{is}^{LD}/\omega_{\textrm{r}}$
7. [+4(is-1)] $\gamma_{is}^{n=+1}/\omega_{\textrm{r}}$
8. [+4(is-1)] $\gamma_{is}^{n=-1}/\omega_{\textrm{r}}$
