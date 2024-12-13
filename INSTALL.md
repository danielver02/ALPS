# ALPS: The Arbitrary Linear Plasma Solver

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
- MPI (e.g., Open MPI or MPICH) - likewise, this should also be the latest version
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
brew install gcc open-mpi
```

There have been reports about issues with Open MPI when used with certain compilers on macOS. In that case, it is worth trying MPICH instead of Open MPI:

```
brew unlink open-mpi
brew install mpich
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

The `Makefile.in` in the repository has been generated with `automake-1.15` from `Makefile.am`. If you have a different version of `automake`, `make` may fail. In that case, start again with

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

## Building the documentation

ALPS uses [Ford](https://forddocs.readthedocs.io/en/latest/) to build its documentation. The documentation is automatically built and deployed to [github.io](https://danielver02.github.io/ALPS/) by the [doc workflow](https://github.com/danielver02/ALPS/blob/master/.github/workflows/doc.yml). To build the documentation locally, follow the [Build documentation](https://github.com/danielver02/ALPS/blob/07a4f8dc996ff76729edeedf5c2a0dc1a5c3028b/.github/workflows/doc.yml#L25-L32) step in the workflow, summarized here:
1. Install `ford` by e.g. `pip install ford`. See [Ford documentation](https://forddocs.readthedocs.io/en/latest/) for details
2. Create a `docs` directory by `mkdir docs`
3. Add a line `title: Readme` to the top of [README.md](./README.md) and copy it to `docs/index.md`
4. Add a line `title: Install` to the top of [INSTALL.md](./INSTALL.md) and copy it to `docs/INSTALL.md`
5. Run `ford ford_project.md`
6. Open `docs/index.html` in a browser to view the documentation

### Adding static pages to documentation

The [README.md](./README.md) and [INSTALL.md](./INSTALL.md) files are added to the [Ford documentation](https://danielver02.github.io/ALPS/) as static pages. You can add more static pages to the documentation by
1. Add the content in a markdown file to the repository.
2. Add a `title: ` line to the beginning of the file and copy it to `docs/` in the [doc workflow](https://github.com/danielver02/ALPS/blob/master/.github/workflows/doc.yml). See steps 3-4 in the previous section, or the [Build documentation](https://github.com/danielver02/ALPS/blob/07a4f8dc996ff76729edeedf5c2a0dc1a5c3028b/.github/workflows/doc.yml#L25-L32) step in the workflow.
3. Add the name of the markdown file as a new line under `ordered_subpage` in [ford_project.md](./ford_project.md)
