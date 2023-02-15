AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=zero,overflow -Wunused -funroll-loops -g -fbacktrace

ALPS_LDADD = -lm -llapack -lblas

bin_PROGRAMS = ALPS

ALPS_SOURCES = ALPS_var.f90 ALPS_io.f90 ALPS_com.f90 ALPS_analyt.f90 ALPS_fns_rel.f90 ALPS_fns.f90
ALPS_SOURCES += ALPS.f90
