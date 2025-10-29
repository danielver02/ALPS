#ALPS Best Practices

##Momentum Space Resolution
If the input distribution is symmetric and resonances occur close to ppar=0, it is sometimes useful to work with an odd number of `npar` steps (e.g., `npar=401` instead of `npar=400`). This centres the distribution not right on one grid point, so that the code can evaluate the resonance integral more safely.

###Convergence Tests
In general, there is no clear rule that prescribes the optimal resolution of the f0-table. Therefore, it is good practice to run convergence tests with different numbers of `nperp` and `npar`. When the ALPS results do not depend on a further increase in the numbers `nperp` and `npar` anymore, then a numerically suitable resolution has been found. This point is particularly important regarding `npar`.

##Zero-Frequency Solutions
When `ALPS` calculates solutions with very small Re(omega) - including solutions for non-propagating modes such as mirror modes of oblique firehose solutions - the code sometimes struggles with the resonance handling of the n=0 resonance near ppar=0. The reason is that small changes in Re(omega) lead to fine structure in the integrand near the resonance.

To overcome this issue for non-relativistic particle species, we recommend performing these calculations in a different reference frame, drifting with respect to the plasma frame along the background magnetic field. Due to Galilean invariance, the physics of the solutions does not depend on the reference frame. 

The easiest approach is to move into a reference frame that moves with 1 Alfvén speed of the reference species with respect to the plasma frame. This is achieved by shifting the ppar coordinates by the appropriate value. For the reference species (in most cases, protons), the drift momentum would be 1, and all other species would be shifted by a drift momentum corresponding to their mass ratio.

The initial guess for the real frequency is then Doppler-shifted by `kpar` multiplied with the drift speed of the reference frame. If 1 Alfvén speed of the reference species is used, then the Doppler shift is simply `kpar`. The solutions can be easily transferred back into the plasma frame by subtracting `kpar` times the drift speed of the reference frame from the real part of the frequency of the solution.

##Relativistic Caveats
For relativistic calculations with `AC_method=1`, it is generally *not* recommended to use `fit_type_in=3`. This fit type is only useful when Jüttner-like f0-tables are used, but without a relativistic treatment of the dispersion relation. For all other cases, it is recommended to use `fit_type_in=4` or `fit_type_in=5`, depending on whether the distribution depends on pparbar or not (i.e., when anisotropies are to be accounted for).

When the error message `ERROR: Resonance integration covers entire subluminal cone.` appears, this means that a resonance occurs at a low gamma value, where not sufficient pparbar points are defined. This error can be resolved in multiple ways. For once, a reduction in `ngamma` lowers the resolution of the gamma-grid, leading to fewer gamma-steps with low pparbar. In addition, an increase in `npparbar` increases the number of points in pparbar-space and thus the number of steps over which the code can treat the resonances appropriately. A combination of both steps (lowering `ngamma`, increasing `npparbar`) often resolves this error. If not, it may help to reduce `positions_principal`.

The autoscaling function in generate_distribution does not work very well for highly relativistic input distributions. When the f0-table is to be created with `generate_distribution`, it is recommended to turn `autoscaleS=F` and to choose manually the values for `maxPperpS` and `maxPparS`. When running `generate_distribution` on the input file, the integration should lead to a value close to unity. This result indicates that the momentum space is well resolved.
