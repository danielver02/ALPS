# Tutorial

This is a tutorial for ALPS. It will guide you through the setting up of some basic input files, the running of the code, and the basic output. For more details, we refer to the [ALPS Input](input.md) page, the [ALPS Output](output.md) page, and the [ALPS Documentation](http://alps.space).

## Authors

Kristopher Klein   (kgklein@arizona.edu)  
Daniel Verscharen  (d.verscharen@ucl.ac.uk)

## Contents

1. Before getting started
2. Installing ALPS
3. Setting up input distributions
4. Running ALPS on f0-tables
5. Interpolation of input distributions
6. Using analytical expressions for the background distribution
7. Running bi-Maxwellian species


## 1. Before getting started

Before starting with the steps described in this tutorial, we recommend that you familiarise yourself with the code paper.

[Verscharen, D., Klein, K. G., Chandran, B. D. G., Stevens, M. L., Salem, C. S.,
and Bale, S. D.: ALPS: the Arbitrary Linear Plasma Solver, J. Plasma Phys. 84,
905840403, 2018, doi: 10.1017/S0022377818000739](http://doi.org/10.1017/S0022377818000739)

You don't need to go through all details, but it is certainly helpful to know what ALPS does and doesn't calculate. The code paper also explains the numerical techniques used in the code, and the [ALPS Documentation](http://alps.space) often refers explicitly to equations and sections in the code paper. We also recommend checking the [Readme](README.md) file.

## 2. Installing ALPS

This tutorial assumes that you have a working copy of ALPS on your computer, including all the required dependencies. You can find the installation guide [here](INSTALL.md). Make sure you have a version of ALPS that compiled completely without error messages after typing

>     ./configure  
>     make

## 3. Setting up input distributions

### 3.1 The format of f0-tables

ALPS calculates the linear Vlasov-Maxwell dispersion relation for a plasma with arbitrary background distribution functions $f_{0j}$. These input distributions are called **f0-tables**.

All f0-tables for ALPS are stored in the folder `./distribution` within the ALPS directory structure. The format of f0-tables is a simple ASCII table with three columns. The columns can be separated with tabs or spaces, but don't add empty lines in the file. The three columns have the following meanings:

1. Perpendicular momentum $p_{\perp}$ normalised to $m_{\mathrm p}v_{\mathrm A}$
2. Parallel momentum $p_{\parallel}$ normalised to $m_{\mathrm p}v_{\mathrm A}$
3. Value of the background distribution function $f_{0j}$ at momentum $(p_{\perp},p_{\parallel})$ normalised to $(m_{\mathrm p}v_{\mathrm A})^{-3}$

In these definitions, $m_{\mathrm p}$ is the proton mass and $v_{\mathrm A}$ is the proton Alfvén speed. The distribution function $f_{0j}$ must be normalised to one if integrated over cylindrical momentum space.

### 3.2 Generating f0-tables

In most cases, you may want to use an f0-table generated from observations or numerical simulations. If you can create such a table directly in the ALPS format, that would be ideal. Simply store it in an ASCII file with the following format:

>     dist_name.N.array

where *dist_name* is the name that you want to give the distribution, and *N* is an integer number indicating the plasma species (e.g., 1 for protons, 2 for electrons).

In some cases, you may need to interpolate distributions from an irregular momentum grid to a regular grid in the ALPS coordinates. ALPS comes with an interpolation routine, which is described [further below](#5.-interpolation-of-input-distributions).

In other cases, you may want to create an f0-table based on a pre-defined function (e.g., Maxwellian, bi-Maxwellian, or $\kappa$-distribution). For this tutorial, let's create a simple Maxwellian distribution and run this through ALPS.

For this exercise, use your terminal to go into the ./distribution folder:

>     cd ./distribution

In the standard ALPS package, you'll find a number of input files to generate distributions. Let's open the file `test_ICW_dist.in` which includes all parameters needed to create Maxwellian f0-tables for the protons and for the electrons. The various parameters are explained in comments in this input file. Feel free to modify some of these parameters if you want to experiment with it a bit. In general, it's a good idea to use one of the test files as the basis to create new input files for your purposes.

When a quantity has the index `ref` added to it, this refers to the reference species. For example, `m_j/m_ref` stands for the mass of the current species *j* in units of the mass of the reference species, which in the example case is the protons.

Let's generate the distributions for this example. Simply execute the command

>     ./generate_distribution test_ICW_dist.in

The code now generates the following output:

```
Species  1
 Integration:       9.9958E-01
 Norm:              5.9862E-02
 pperp_max:         1.0392E+01
 ppar_max:          6.0000E+00

 Fit type:          1
 ideal fit_1:       5.9862E-02
 ideal fit_2:       1.0000E+00
 ideal fit_3:       0.0000E+00
 ideal perpcorr:    3.3333E-01
============================
Species  2
 Integration:       9.9958E-01
 Norm:              1.4128E+04
 pperp_max:         1.4003E-01
 ppar_max:          1.4003E-01

 Fit type:          1
 ideal fit_1:       1.4128E+04
 ideal fit_2:       1.8360E+03
 ideal fit_3:       0.0000E+00
 ideal perpcorr:    1.8360E+03
============================
```

Species 1 corresponds to the protons, and species 2 corresponds to the electrons in this case. The output gives us some information about the integration and normalisation of $f_{0j}$ as well as some suggestions for the fit parameters (*ideal fit* and *ideal perpcorr*) that will become important later for the hybrid-analytic continuation.

In addition to this output, the code has also created two files: `test_ICW.1.array` and `test_ICW.2.array`. These are the the two f0-tables, which now include the pre-defined Maxwellian distributions for the protons and electrons. Let's open the proton file `test_ICW.1.array` and have a look:

```
0.0000000000000000       -6.0000000000000000        1.3885214326235464E-017
0.0000000000000000       -5.9500000000000002        2.5237337794355358E-017
0.0000000000000000       -5.9000000000000004        4.5641827072255333E-017
0.0000000000000000       -5.8499999999999996        8.2131741067563143E-017
0.0000000000000000       -5.7999999999999998        1.4705763083579771E-016
0.0000000000000000       -5.7500000000000000        2.6199477385677698E-016
0.0000000000000000       -5.7000000000000002        4.6443636702139106E-016
0.0000000000000000       -5.6500000000000004        8.1919697073946868E-016
...                      ...                        ...
```

As discussed above, the first column is the normalised perpendicular momentum, the second column is the normalised parallel momentum, and the third column is the value of the distribution function. The file contains 29,161 lines, which is just the combination of all perpendicular and parallel momentum steps. These have been defined by the lines

```
nperp=120  
npar=240
```

in the file `test_ICW_dist.in`. Actually, the number 29,161 corresponds to 121\*241, which is greater than 120\*240 due to the inner and outer boundaries of the integration space. The variables `nperp` and `npar` refer to the integration domain only, which is one step smaller in each dimension than the f0-table. This is an important point to consider when setting up f0-tables:

> If you have an f0-table of X-by-Y entries in momentum space, set `nperp=X-1` and `npar=Y-1` in ALPS.  
> Also ensure that the perpendicular momentum grid begins at a momentum of zero.

Now that we have two f0-tables in the correct format, let's run ALPS on these!




## 4. Running ALPS on f0-tables

Like for the generation of the input distribution above, we will run one of the test input files for this tutorial to demonstrate the running of ALPS. You can find this test input in the folder `./tests` within the ALPS directory structure. Assuming that your terminal is still in `./distribution`, let's change directory to `./tests`:

>     cd ../tests

The input file that we want to use is called `test_ICW.in`. The name of the input file doesn't need to be the same as the name of the f0-table files (*dist_name*), but it helps to use a consistent nomenclature. Let's open `test_ICW.in`. It includes all the parameters that ALPS needs to run. Like in the case of the file `test_ICW_dist.in` above, the input file is commented to help you understand the meaning of the parameters. For more details and a reference to all of the parameters, please have a look at our [ALPS Input](input.md) page.

We focus on a few key entries here for now. First, it's important to give the code the correct `nperp` and `npar` values. They are defined in the same way as above. So, here is a reminder:

> If you have an f0-table of X-by-Y entries in momentum space, set `nperp=X-1` and `npar=Y-1` in ALPS.

In the example file, we see that `nperp=120` and `npar=240`, which are the correct entries for our example.

Further down, you find the line

```
arrayName='test_ICW'
```

This line defines the name of the f0-tables that ALPS will use. Since we tell ALPS to work with two species (from the line `nspec=2`), the code will look for two files

>     ./distribution/test_ICW.1.array  
>     ./distribution/test_ICW.2.array

which we had created earlier.

If you scroll further down, you'll find the following block:

```
!Initial guess of complex frequency for first solution
!Only used when use_map=.false.
!Need to have # of name lists equal to nroots
&guess_1
g_om=1.5508d-01   !real frequency
g_gam=-1.3429d-05 !imaginary frequency
/
```

This block defines the *initial guess* for the frequency that the code will use to search for solutions. By picking the appropriate guess, you can select which plasma mode the code should follow. Remember that all frequencies are normalised to the gyro-frequency of the protons.

Further down, you can find definitions for the species. The first species is defined through the following:

```
!Species parameters list; species 1
&spec_1
nn=1.d0  !relative density; n_j/n_ref
qq=1.d0  !relative charge;  q_j/q_ref
mm=1.d0  !relative mass;    m_j/m_ref
ff=1     !# of fits to be used.
relat=F  !relativistic (T) or non-relativistic (F) species
log_fit=T!log (T) or linear (F) fit
use_bM=F !numerically integrate array (F) or use bi-Maxwellian (T)
/

!Initial Fit Values; species 1, function 1
&ffit_1_1
fit_type_in=1   !Kind of fit function (see documentation)
fit_1=5.986D-2  !Suggested values for parameters,
fit_2=1.d0      !e.g. generated by generate_distribution
fit_3=0.d0
fit_4=0.d0
fit_5=0.d0
perpcorr=3.33D-1   !renormalization factor
/
```

You can recognise many of the parameters from our `test_ICW_dist.in` file, like the charge and the mass of the particle species The block `&ffit_1_1` includes the fit parameters for the hybrid-analytic continuation. As you can see, this file includes the suggestions for the ideal fit that `generate_distribution` has given us earlier.

The key point about the fit parameters is that, for all damped solutions, you want to find a good representation of the distribution function so that the Landau-contour integral is precise. The procedure is described in the code paper in Section 3.2. If you have any closed mathematical expression to use, you can also input this for the analytic continuation as described [here](#6.-using-analytical-expressions-for-the-background-distribution).

Once we are happy with the parameters in `test_ICW.in`, let's run ALPS. Depending on your MPI configuration, you can run the code directly with the following commands (it's best to go back into the main folder of the ALPS directory structure):

>     cd ../
>     mpirun -np 4 ./src/ALPS tests/test_ICW.in

This command will run 4 instances of ALPS via MPI. Depending on your computer and the MPI setup, you may need to add the option `--oversubscribe` to ensure that 4 instances of MPI can actually be executed. The ALPS binary (`./src/ALPS`) expectes only one command-line argument: the name of the input file to run on. The number of processes (in this case 4), must always be an even number greater than 2, even if your computer doesn't have four or more cores. In those cases, the `--oversubscribe` option is particularly important. In general, the code scales quite well, especially for calculations of the dispersion relation with many plasma species and high perpendicular wavenumbers (compared to the species gyro-radii).

If everything was successful so far, you'll now start to get the ALPS standard output:

```
===========================================================
I                                                         I
I                       A  L  P  S                        I
I              Arbitrary Linear Plasma Solver             I
I                                                         I
I                       Version X.X                       I
I                                                         I
I  Kristopher Klein   (kgklein@arizone.edu)               I
I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)            I
I                                                         I
===========================================================
All processes are up and running.
Reading from Input File: test_ICW
GUESS ROUTINE:
Intial Guess   1 :     1.5508E-01   -1.3429E-05
SPECIES PARAMETERS:
Species   1 :
ns/nREF =  1.0000E+00 | qs/qREF =  1.0000E+00 | ms/mREF =  1.0000E+00
Number of fitted functions =    1
Relativistic effects = F
Species   2 :
ns/nREF =  1.0000E+00 | qs/qREF = -1.0000E+00 | ms/mREF =  5.4466E-04
Number of fitted functions =    1
Relativistic effects = F
-=-=-=-=-=-=-=-=-
```

This output summarises a lot of the information that is useful. We recommend that you read this output carefully, or even that you pipe the output into a log file for later reference by using

>     mpirun -np 4 ./src/ALPS tests/test_ICW.in > test_ICW.output &

Further down in the ALPS standard output, you'll find lines similar to this:

```
kperp:     1.0000E-03 kpar:     1.0275E-01
 Converged after iteration    7
 D(    1.4337E-01   -8.9137E-06)=    -2.0138E-18   -1.2447E-18
kperp:     1.0000E-03 kpar:     1.0557E-01
 Converged after iteration    8
 D(    1.4729E-01   -8.6524E-06)=    -2.9456E-19   -1.6708E-19
kperp:     1.0000E-03 kpar:     1.0846E-01
 Converged after iteration    8
 D(    1.5132E-01   -8.3974E-06)=     7.8943E-19    5.5227E-19
```

These lines tell us that the code is actually now finding solutions, as it scans through the wavevector space. The scan options had been defined in `test_ICW.in`, and details on how this works are given in our [ALPS Input](input.md) page.

When the code has finished, it has produced a number of output files, which you can find in the folder `./solution` in the ALPS directory structure. Let's look at the file `test_ICW.scan_kpara_1.root_1`. As the name suggests, this file contains the scan result for a scan along the parallel wavenumber $k_{\parallel}$ for root number 1 (we only scan along one root in this example).

The file has the following format:

```
1.0000E-03    1.0000E-01    1.4424E-01   -1.0449E-05
1.0000E-03    1.0275E-01    1.4337E-01   -8.9137E-06
1.0000E-03    1.0557E-01    1.4729E-01   -8.6524E-06
1.0000E-03    1.0846E-01    1.5132E-01   -8.3974E-06
1.0000E-03    1.1144E-01    1.5546E-01   -8.1484E-06
1.0000E-03    1.1450E-01    1.5971E-01   -7.9055E-06
...           ...           ...          ...
```

The format of the columns is as follows:

1. Wavenumber $k_{\perp}$ perpendicular to the background magnetic field in units of $\Omega_{\mathrm p}/v_{\mathrm A}$
2. Wavenumber $k_{\parallel}$ parallel to the background magnetic field in units of $\Omega_{\mathrm p}/v_{\mathrm A}$
3. Real part of the wave frequency $\omega$ in units of $\Omega_{\mathrm p}$
4. Imaginary part of the wave frequency $\omega$ in units of $\Omega_{\mathrm p}$

In these definitions, $\Omega_{\mathrm p}$ is the proton gyro-frequency.

Congratulations! You have found the solutions to the Maxwellian example case for the ion-cyclotron wave. Now you can also experiment with the other test cases in the folders `./distribution` and `./tests`. The code also includes a full test suite which you can launch with

>     cd ./tests  
>     ./run_test_suite.sh

The given test cases cover a range of typical applications, so there should be good starting points for your purposes.

The code also creates additional output files for the heating contributions and for the eigenfunctions of the solutions. The general output format of ALPS is described on our [ALPS Output](output.md) page.


## 5. Interpolation of input distributions

In many cases, you may want to use an f0-table based on a data file that is not in the same format as the ALPS format (e.g., from spacecraft observations) or isn't equally spaced in momentum space as required by ALPS. For those cases, ALPS provides an interpolation routine, which you can find in the folder `./interpolation`. The ALPS code also comes with a test case to illustrate the the use of the interpolation routine. We'll go through this example here.

Let's have a look at the file `test_interp_coarse.array`. This file includes a table of the distribution function, but not in the format as needed by ALPS:

```
3.48994945E-04   9.99390800E-03  0.999899983      
3.67403193E-03   9.30061750E-03  0.999899983     
6.57521421E-03   7.53435818E-03  0.999899983    
8.71784426E-03   4.89889691E-03  0.999899983    
9.85473860E-03   1.69827254E-03  0.999899983    
9.85473767E-03  -1.69827335E-03  0.999899983    
...              ...             ...  
```

The overall format is the same as in the f0-table files above in terms of the meaning of the columns:

1. Perpendicular momentum $p_{\perp}$ normalised to $m_{\mathrm p}v_{\mathrm A}$
2. Parallel momentum $p_{\parallel}$ normalised to $m_{\mathrm p}v_{\mathrm A}$
3. Value of the background distribution function $f_{0j}$ at momentum $(p_{\perp},p_{\parallel})$ normalised to $(m_{\mathrm p}v_{\mathrm A})^{-3}$

Now let's look at the input file `test_interp.in` for the interpolation routine. This file includes a number of comments to help you with the setting up of the interpolation parameters. The file asks the interpolation routine to take the irregular table from `test_interp_coarse.array` and to interpolate it onto a grid in the ALPS format with `nperp=50` and `npar=100`. If needed, the code can also re-normalise and scale the distribution function depending on your needs.

Let's run the interpolation routine on the input file:

>     ./interpolation test_interp.in

If successful, we will get the following output:

```
Number of points in the coarse grid:          700

Properties of the fine grid:
Number of grid points in nperp:                50
Number of grid points in npar:                100
Maximum Pperp:                         2.9663E+00
Minimum Pperp:                         0.0000E+00
Maximum Ppar:                          3.0082E+00
Minimum Ppar:                         -3.0082E+00

Writing output to file test_interp.in.array   
```                                                                      

At the same time, the folder now includes a new file `test_interp.in.array`, which includes the interpolated f0-table in the correct format for ALPS:

```
0.0000000000000000       -3.0081663100000000        2.5699885852486149E-005
0.0000000000000000       -2.9480029837999999        3.4671386965111912E-005
0.0000000000000000       -2.8878396575999998        4.7473668668616121E-005
0.0000000000000000       -2.8276763314000002        6.5558031421796410E-005
0.0000000000000000       -2.7675130052000001        9.0732738996456109E-005
0.0000000000000000       -2.7073496790000000        1.2528393376824723E-004
0.0000000000000000       -2.6471863527999999        1.7212176003060674E-004
0.0000000000000000       -2.5870230265999998        2.3495069283188794E-004
0.0000000000000000       -2.5268597004000002        3.1846300082066743E-004
...                       ...                       ...
```

You can copy this file over into the `./distribution` folder and give it the appropriate naming with the convention given above. Then ALPS can use this interpolated file for the calculation.


## 6. Using analytical expressions for the background distribution

In some instances, you may have an analytical expression for the distribution function (which is not covered by `generate_distribution`) that you want to run through ALPS. For example, you may want to analyse the instability of a model distribution that was derived analytical from a model. ALPS can also do this.

For such a calculation, we need to define the background distribution in the code directly. Let's open the file `./solution/distribution_analyt.f90`. This is a Fortran file that allows you to feed an analytical expression for the background distribution into ALPS.

If you have the standard version of ALPS from the repository, this f90 file already contains some variable declarations and an example case. Let's scroll down to the following part:

```
select case(is)

  case(1) ! Species 1

    ! The example below illustrates how to set up a Maxwellian with beta = 1:
    beta=1.d0
    ms=1.d0

    f0=(pi**(-1.5d0) /((ms * beta )**(3.d0/2.d0) )) * exp( -( ppar**2/( beta * ms)&
      + (pperp**2)/( beta * ms ) ) )

  case(2) ! Species 2
     ...
```

The different *cases* here refer to the species. This allows you to use different equations for, say, electrons and protons. Whatever is defined as the variable `f0` as a function of the variables `pperp` and `ppar` inside these *cases* defines the f0-table. If you require additional variables or functions to define your model distribution, please feel free to define them locally in this Fortran function.

**Important**: You must ensure that the defined background distribution function is normalised to one. The code doesn't automatically normalise the distribution. However, the ALPS standard output will list the integration of the distribution, so that you can double check.

The example above defines a simple Maxwellian for the protons and for the electrons as a simple example.

Once you have defined the background distribution in `distribution_analyt.f90`, it's important to re-compile the code, so that it includes the defined function:

>     cd ../  
>     make

For the next step, we need to create the f0-table based on the defined function. This is done with `generate_distribution` as for the other analytical cases. We have an example input file `./distribution/test_analytical_dist.in`, which does exactly that. For both species, this file defines the distribution type as zero in the lines:

```
distributions=0  !Type of distribution (see documentation)
```

If you set `distributions=0`, the function `generate_distribution` will use the function defined in `distribution_analyt.f90` to create an f0-table. Let's run

>     cd ./distribution    
>     ./generate_distribution test_analytical_dist.in

The standard output explains a couple of details, which we had already encountered above. Most importantly, the programme now creates two new files:

>     test_analytical.1.array   
>     test_analytical.2.array

These files include the f0-tables according to the analytical function.


Now it's important to tell ALPS to use the defined distribution functions. This happens in the ALPS input file. In the folder `./tests`, there is an example file called `test_analytical.in`, which we will use in this example. As expected, it commands ALPS to use the f0-tables that we have just created:

```
arrayName='test_analytical'
```

It also makes sense to use the pre-defined background distribution function from `distribution_analyt.f90` also for the analytical continuation. You can achieve this by setting

```
ff=0
```

in the block of the corresponding species. If you set `ff=0`, the code will look up the function in `distribution_analyt.f90` to evaluate the Landau-contour integral.

If you now run this test case through

>     mpirun -np 4 ./src/ALPS ./tests/test_analytical.in

the code will fully calculate the dispersion relation based on the distribution function defined in `distribution_analyt.f90`.

A particular strength is that you can use that distribution for the creation of the f0-table, the Landau-contour integral (analytic continuation), or for both. Sometimes, you may only want to use it for the f0-table and rely on ALPS's internal fitting. If that is the case, simply set `ff` to the corresponding number of fits and define the fits as usual. At other times, you may only want to use the function from `distribution_analyt.f90` as the fit function, but use a different f0-table. In that is the case, point ALPS to your usual f0-table files, but set `ff=0`.


## 7. Running bi-Maxwellian species

In some instances, it may be useful to assume a Maxwellian or bi-Maxwellian distribution for one or more species in the system. For example, you may want to calculate the proton susceptibilities based on an f0-table, but you're happy to work with a simple bi-Maxwellian electron population.

In that case, ALPS includes an implementation of the [NHDS code](https://github.com/danielver02/NHDS) to accelerate the calculation. ALPS then does not integrate over the bi-Maxwellian distribution explicitly but uses the known analytical approximations.

The relevant option for a bi-Maxwellian calculation is the logical flag `use_bM` in the ALPS input file under the section for the corresponding species. If you set
```
use_bM=T
```
for any species, its contributions to the dielectric tensor will be calculated with the NHDS routines. Also the Landau-contour integration will be done through NHDS.

Once `use_bM` is set to true, ALPS will look for the parameters used to define the bi-Maxwellian properties of the corresponding species. These are given in a species-dependent block in the ALPS input file, which may look like this:

```
!Bi-Maxwellian parameters; for species 1
!Only used if use_bM=T
&bM_spec_1
bM_nmaxs=500           !Maximum number of resonaces to consider
bM_Bessel_zeros=1.d-50 !Amplitude limit for Bessel function
bM_betas=1.d0          !Plasma beta for species 1
bM_alphas=1.d0         !Tperp/Tpar for species 1
bM_pdrifts=0.d0        !Momentum drift, norm. to m_ref v_A,ref
/
```

The ALPS code suite includes a test case called `./tests/test_bimax.in` that uses the NHDS routines for a quick calculation of the dispersion relation for a bi-Maxwellian plasma.
