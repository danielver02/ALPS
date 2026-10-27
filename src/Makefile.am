AUTOMAKE_OPTIONS = subdir-objects

# Fortran Compiler and Flags
#AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=zero,overflow -Wunused -funroll-loops -g -fbacktrace
#With newer MPIs, ALPS crashes during MPI initialization with -ffpe-trap=zero. 
AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=overflow -Wunused -funroll-loops -g -fbacktrace 

# Linker Libraries
ALPS_LDADD = -lm $(LAPACK_LIBS) $(BLAS_LIBS)

# Output Program
bin_PROGRAMS = ALPS

# Sources (used only for dependency tracking and rule generation)
ALPS_SOURCES = \
    ALPS_rel_cache.f90 \
    ALPS_var.f90 \
    ALPS_diagnostic_math.f90 \
    ALPS_diagnostics.f90 \
    ../distribution/distribution_analyt.f90 \
    ALPS_io.f90 \
    ALPS_com.f90 \
    ALPS_analyt.f90 \
    ALPS_NHDS.f90 \
    ALPS_fns_rel.f90 \
    ALPS_fns.f90 \
    ALPS.f90 
    


# Dependency Tracking
AM_CPPFLAGS = -MMD -MP
DEPDIR = .deps
-include $(wildcard $(DEPDIR)/*.d)

# Explicit module ordering, including clean parallel builds.
ALPS_diagnostics.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_diagnostic_math.$(OBJEXT)
ALPS_io.$(OBJEXT) ALPS_com.$(OBJEXT) ALPS_analyt.$(OBJEXT): ALPS_diagnostics.$(OBJEXT)
ALPS_analyt.$(OBJEXT): ALPS_io.$(OBJEXT) ../distribution/distribution_analyt.$(OBJEXT)
ALPS_NHDS.$(OBJEXT): ALPS_var.$(OBJEXT)
ALPS_fns_rel.$(OBJEXT): ALPS_rel_cache.$(OBJEXT) ALPS_analyt.$(OBJEXT)
ALPS_fns.$(OBJEXT): ALPS_io.$(OBJEXT) ALPS_fns_rel.$(OBJEXT) ALPS_NHDS.$(OBJEXT)
ALPS.$(OBJEXT): ALPS_com.$(OBJEXT) ALPS_fns.$(OBJEXT)
