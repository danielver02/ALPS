#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(ALPS, 1.0, d.verscharen@ucl.ac.uk)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR(src/)

#AC_CONFIG_MACRO_DIRS([m4])

: ${FCFLAGS=""}

# Checks for programs.
AC_PROG_CC([gcc])
AC_PROG_FC([mpifort mpif90 f90 fort gfortran])
AC_LANG([Fortran])

# Checks for libraries.


AC_SEARCH_LIBS(dgesv,blas)
AC_SEARCH_LIBS(dgesv,openblas)
AC_SEARCH_LIBS(dgesv,lapack)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
#AC_SEARCH_LIBS([h5open_f],[hdf5_fortran])

#AC_ARG_WITH([hdf5],[AS_HELP_STRING([−-with−hdf5=DIR],[enable path for HDF5 library])],[WITH_HDF5="yes" HDF5PATH="$withval"], [WITH_HDF5="no"])
#AM_CONDITIONAL(WITH_HDF5, test x$WITH_HDF5 = xyes)


AC_OUTPUT(Makefile src/Makefile interpolation/Makefile distribution/Makefile)
