# A  L  P  S  Arbitrary Linear Plasma Solver

## Authors

Kristopher Klein   (kgklein@arizona.edu)
Daniel Verscharen  (d.verscharen@ucl.ac.uk)

This is the ALPS code: the Arbitrary Linear Plasma Solver.

## WARNING

IN ITS CURRENT FORM, THE CODE IS STILL UNDER DEVELOPMENT. WE RECOMMEND THAT YOU
USE IT VERY CAREFULLY OR ONLY IN DIRECT COMMUNICATION WITH THE CODE DEVELOPMENT
TEAM. THE CODE WILL BE MADE MORE USEABLE AND SUSTAINABLE IN THE FUTURE.

## CONTENTS

1. What is ALPS?
2. Acknowledgements
3. Installing the ALPS Code
4. Running the ALPS Code
5. Documentation for Input/Output Data
6. List of Error Codes
7. License

## WHAT IS ALPS?

ALPS is a parallelised numerical code that solves the Vlasov-Maxwell dispersion
relation in hot (even relativistic) magnetised plasma. ALPS allows for any
number of particle species with arbitrary gyrotropic equilibrium distribution
functions supporting waves with any direction of propagation with respect to
the background magnetic field.

If you use the code for a science publication, please provide the code website
on github.com/danielver02/ALPS and cite the code paper:

Verscharen, D., Klein, K. G., Chandran, B. D. G., Stevens, M. L., Salem, C. S.,
and Bale, S. D.: ALPS: the Arbitrary Linear Plasma Solver, J. Plasma Phys. 84,
905840403, 2018, doi: 10.1017/S0022377818000739

## ACKNOWLEDGEMENTS

The development of the ALPS code was supported by NASA Grant NNX16AG81G. We will
present document more details about the numerics on the website
http://www.alps.space. The code developers appreciate support from the UK Science
and Technology Facilities Council (STFC) Ernest Rutherford Fellowship ST/P003826/1,
STFC Consolidated Grants ST/S000240/1 and ST/W001004/1, and the Open Source
Software Sustainability Funding programme from UCL's Advanced Research Computing
Centre and UCL's eResearch Domain. We appreciate software engineering support by
David Stansby (UCL).

##  INSTALLING THE ALPS CODE

For advice on the installation of the code, please check the `INSTALL` file in the
main directory.

##  RUNNING THE ALPS CODE

ALPS works with input files that specify the plasma and numerical parameters for
the calculation. We recommend that you start by checking out the provided test
cases as a guidance for the creation of input files. These test cases are listed
in the scripts `run_test.sh` and `run_test_suite.sh`. All associated input files have
a name starting with `test_`.

You can execute the ALPS code through the following command:

```
mpirun -np <NP> ./src/ALPS <input_file.in>
```

where `<NP>` is the number of processors you want to use. This number must be greater
than or equal to 4, and it must be an even number. `<input_file.in>` is the input file
that includes all parameters for your run.

On some systems, depending on the MPI configuration, the oversubscribe flag is
required. In this case, the above command must be replaced with

```
mpirun -np <NP> --oversubscribe ./src/ALPS <input_file.in>
```

## DOCUMENTATION FOR INPUT/OUTPUT DATA

### Input Velocity Distribution Function
  - folder: `./distribution`
  - (optionally created by subroutine generate_distribution
      `distribution/generate_distribution.f90`)
  - file: `distribution/<arrayName>.<is>.array`
  - data: `pp(is,iperp,ipar,1:2)`,`f0(is,iperp,ipar)`

### Fit Parameters for Hybrid Analytic Continuation
  - folder: `./distribution`
  - created by subroutine `determine_fit_parameters` (`ALPS_analyt.f90`)
  - file: `distribution/<runname>.fit_parameters.<is>.out`
  - data: `iperp`, `params`

### Result from Fit Routine for Hybrid Analytic Continuation
  - folder: `./distribution`
  - created by `subroutine output_fit` (`ALPS_analyt.f90`)
  - file: `distribution/<runname>.fit_parameters.<is>.out`
  - data: `pp(is,iperp,ipar,1:2)`,`fit_result(is,iperp,ipar)`,
      `abs(fit_result(is,iperp,ipar,1:2)- f0(is,iperp,ipar))/f0(is,iperp,ipar)`

### Map (grid in omega,gamma) of Dispersion Solutions
  - created by `subroutine map_search` (`ALPS_fns.f90`)
  - file: `solution/<runname>.map`
  - data: `ir`,`ii`,`om(ir,ii)`,`log10(val(ir,ii))`,`cal(ir,ii)`

### Listing of Roots
  - created by `subroutine refine_guess` (`ALPS_fns.f90`)
  - file: `solution/<runname>.roots`
  - data: `iw`,`wroots(iw)`,`log10(abs(tmpDisp))`,`tmpDisp`

### Scans of the Dispersion Solutions through k-space
  - created by `subroutine om_scan` (`ALPS_fns.f90`)
  - file: `solution/<runname>.scan_<option>.root_<rootnumber>`
  - data: `kperp`,`kpar`,`om(ir,ii)`

### Evaluation of the Eigen Vectors
  - created by `subroutine om_scan` (`ALPS_fns.f90`)
  - file: `solution/<runname>.eigen_<option>.root_<rootnumber>`
  - data: `kperp`,`kpar`,`om(ir,ii)`,`ef`,`bf`,`Us`,`ds`

### Evaluation of the Heating Rates
  - created by `subroutine om_scan` (`ALPS_fns.f90`)
  - file: `solution/<runname>.heat_<option>.root_<rootnumber>`
  - data: `kperp`,`kpar`,`om(ir,ii)`,`Ps`

## LIST OF ERROR CODES

```
error_id=0
   - if mod(nproc,2) .ne. 0
error_id=1
   - Does requested input namelist exist?
error_id=2
   - if ((1+npar).LT.n_params)
```

### ERRORS TO ADD:
- if the size of the input distribution function array is not equal to
the allocated distribution function array in ALPS

## LICENSE

BSD 2-Clause License

Copyright (c) 2023, Kristopher G. Klein and Daniel Verscharen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
