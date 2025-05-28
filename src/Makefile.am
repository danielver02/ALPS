AUTOMAKE_OPTIONS = subdir-objects

# Fortran Compiler and Flags
AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=zero,overflow -Wunused -funroll-loops -g -fbacktrace 

# Linker Libraries
ALPS_LDADD = -lm $(LAPACK_LIBS) $(BLAS_LIBS)

# Output Program
bin_PROGRAMS = ALPS

# Sources (used only for dependency tracking and rule generation)
ALPS_SOURCES = \
    ALPS_var.f90 \
    ALPS_io.f90 \
    ALPS_com.f90 \
    ALPS_analyt.f90 \
    ALPS_NHDS.f90 \
    ALPS_fns_rel.f90 \
    ALPS_fns.f90 \
    ALPS.f90 \
    ../distribution/distribution_analyt.f90


# Dependency Tracking
AM_CPPFLAGS = -MMD -MP
DEPDIR = .deps
-include $(wildcard $(DEPDIR)/*.d)
