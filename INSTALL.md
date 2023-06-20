# A  L  P  S  Arbitrary Linear Plasma Solver

These are the install instructions for the ALPS code: the Arbitrary Linear
Plasma Solver.

## Authors

Kristopher Klein   (kgklein@arizona.edu)
Daniel Verscharen  (d.verscharen@ucl.ac.uk)

## CONTENTS

1. Requirements and Dependencies
2. Setting up Computer Environments
3. Getting the ALPS Code
4. Installing the ALPS Code
5. Execution of Test Runs

## REQUIREMENTS AND DEPENDENCIES

ALPS has the following requirements:

- A UNIX, Linux, or macOS operating system with a working shell
- GNU make, in some cases it is useful to have the autoconf/automake tools
- Fortran 90 compiler (e.g., gfortran) - we recommend using the latest version
  of the compiler to avoid any surprises in the evaluation.
- MPI (e.g., openmpi) - likewise, this should also be the latest version
- BLAS and LAPACK - these two libraries are used for the polyharmonic spline
  interpolation in ALPS. They are directly linked during the compilation

## SETTING UP COMPUTER ENVIRONMENTS

For Ubuntu and macOS users, the following instructions have proven to be useful.
On both systems, we recommend deactivating potential anaconda installations
that could interfere with the ALPS installation:

```
conda deactivate
```

On Ubuntu, the following installation routines obtain the necessary software
packages for ALPS:

```
sudo apt-get install -y libopenmpi-dev
sudo apt-get install -y libopenblas-dev libblas-dev liblapack-dev
```

If an older compiler or MPI version is still installed, it may be necessary to
deinstall this before using apt-get.


On macOS, homebrew is a good way to install the necessary packages:

```
brew install gcc openmpi
```

## GETTING THE ALPS CODE

We recommend pulling the latest version of ALPS from GitHub. For this, go to
the directory where you want to install ALPS and execute the following command:

```
git clone https://github.com/danielver02/ALPS
```

Alternatively, you can also go to the website https://github.com/danielver02/ALPS
directly and download the source code from there. The advantage of using the git
command is that you can now contribute to the development of the ALPS code. If
you make any changes to the code, GitHub will run automatic tests (via workflows)
to ensure that the changes do not break the code.

## INSTALLING THE ALPS CODE

If all requirements are available, the code can be compiled with the following
commands:

```
./configure
make
sudo make install           (this option is only required if you want to
                             make the ALPS executable available to all users)
```

If either of these steps fails, we recommend starting with

```
autoreconf -i -f
```

before the execution of `./configure`.

## EXECUTION OF TEST RUNS

ALPS comes with a selection of test runs that cycle through various test
problems. To execute a small set of tests, execute the following shell script:

```
./run_test.sh
```

This script will test the interpolation routine, the routine to generate pre-
described distribution functions, and a simply fast dispersion relation.

To execute a more complete set of test problems, execute the following shell
script:

```
./run_test_suite.sh
```

This script will test the interpolation routine and then a number of ALPS test
cases, including the generation of the relevant distribution functions. The test
routine outputs time stamps during each steps to compare the speed of the ALPS
code and to facilitate scaling tests. For each test, the script will explain
whether errors occurred or not in any of the tests. The code output itself will
be piped into the .out and .error files in the home directory.

---

For advice on running ALPS, please see the README file.
