# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.72])
AC_INIT([ALPS], [1.0], [d.verscharen@ucl.ac.uk])
AC_CONFIG_SRCDIR([src/])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])

# Optional: macro directory if using custom macros
# AC_CONFIG_MACRO_DIR([m4])

# Set default Fortran flags if not already set
: ${FCFLAGS=""}

# Checks for compilers
AC_PROG_CC
AC_PROG_FC([mpifort mpif90 f90 fort gfortran])
AC_LANG([Fortran])

# -----------------------
# Optional: BLAS/LAPACK support
# -----------------------

AC_ARG_WITH([blas],
    [AS_HELP_STRING([--with-blas=LIBS], [Specify BLAS libraries to link, e.g., -lopenblas])],
    [BLAS_LIBS="$withval"],
    [BLAS_LIBS="-lblas"]
)

AC_ARG_WITH([lapack],
    [AS_HELP_STRING([--with-lapack=LIBS], [Specify LAPACK libraries to link, e.g., -llapack])],
    [LAPACK_LIBS="$withval"],
    [LAPACK_LIBS="-llapack"]
)

# Basic linking test to ensure these libraries work
AC_MSG_CHECKING([for usable BLAS/LAPACK])
SAVE_LIBS="$LIBS"
LIBS="$LIBS $BLAS_LIBS $LAPACK_LIBS"
AC_LINK_IFELSE([AC_LANG_CALL([Fortran], [dgesv])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_FAILURE([Could not link with specified BLAS/LAPACK libraries])]
)
LIBS="$SAVE_LIBS"

AC_SUBST([BLAS_LIBS])
AC_SUBST([LAPACK_LIBS])

# -----------------------
# Optional: HDF5 support (commented out)
# -----------------------
# AC_ARG_WITH([hdf5],
#     [AS_HELP_STRING([--with-hdf5=DIR], [Specify path to HDF5 installation])],
#     [WITH_HDF5=yes; HDF5PATH="$withval"],
#     [WITH_HDF5=no])
# AM_CONDITIONAL([WITH_HDF5], [test "x$WITH_HDF5" = "xyes"])

# Example: Check for h5open_f in hdf5_fortran
# if test "x$WITH_HDF5" = "xyes"; then
#     AC_CHECK_LIB([hdf5_fortran], [h5open_f], [], [
#         AC_MSG_ERROR([HDF5 library not found or unusable])
#     ])
# fi

# -----------------------
# Output files
# -----------------------
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    interpolation/Makefile
    distribution/Makefile
])
AC_OUTPUT
