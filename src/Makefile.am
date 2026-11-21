AUTOMAKE_OPTIONS = subdir-objects

# Fortran Compiler and Flags
#AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=zero,overflow -Wunused -funroll-loops -g -fbacktrace
#With newer MPIs, ALPS crashes during MPI initialization with -ffpe-trap=zero. 
AM_FCFLAGS = -O4 -fbounds-check -ffpe-trap=overflow -Wunused -funroll-loops -g -fbacktrace 

# Linker Libraries
ALPS_LDADD = distribution_analyt.$(OBJEXT) -lm $(LAPACK_LIBS) $(BLAS_LIBS)

# Output Program
bin_PROGRAMS = ALPS

# Sources (used only for dependency tracking and rule generation)
ALPS_SOURCES = \
    ALPS_rel_cache.f90 \
    ALPS_var.f90 \
    ALPS_io.f90 \
    ALPS_analyt.f90 \
    ALPS_com.f90 \
    ALPS_NHDS.f90 \
    ALPS_fns_rel.f90 \
    ALPS_fns.f90 \
    ALPS_nyquist.f90 \
    ALPS_nyquist_growth.f90 \
    ALPS_nyquist_driver.f90 \
    ALPS.f90 
    


# Dependency Tracking
AM_CPPFLAGS = -MMD -MP
DEPDIR = .deps
-include $(wildcard $(DEPDIR)/*.d)

# Each build directory owns its module interfaces and objects. In particular,
# do not share distribution/distribution_analyt.o with generate_distribution:
# its .mod is emitted in the compiler's working directory.
CLEANFILES = *.mod *.smod

distribution_analyt.$(OBJEXT): $(srcdir)/../distribution/distribution_analyt.f90
	$(AM_V_FC)$(FCCOMPILE) -c -o $@ $<

# Direct local USE dependencies (MPI and intrinsic modules are compiler inputs).
# Object prerequisites order clean/parallel builds and propagate provider edits.
ALPS_io.$(OBJEXT): ALPS_var.$(OBJEXT)
ALPS_analyt.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_io.$(OBJEXT) distribution_analyt.$(OBJEXT)
ALPS_com.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_analyt.$(OBJEXT)
ALPS_NHDS.$(OBJEXT): ALPS_var.$(OBJEXT)
ALPS_fns_rel.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_io.$(OBJEXT) ALPS_rel_cache.$(OBJEXT) ALPS_analyt.$(OBJEXT)
ALPS_fns.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_io.$(OBJEXT) ALPS_analyt.$(OBJEXT) ALPS_fns_rel.$(OBJEXT) ALPS_NHDS.$(OBJEXT)
ALPS_nyquist_growth.$(OBJEXT): ALPS_nyquist.$(OBJEXT)
ALPS_nyquist_driver.$(OBJEXT): ALPS_nyquist.$(OBJEXT) ALPS_nyquist_growth.$(OBJEXT) ALPS_var.$(OBJEXT) ALPS_io.$(OBJEXT) ALPS_fns.$(OBJEXT)
ALPS.$(OBJEXT): ALPS_var.$(OBJEXT) ALPS_io.$(OBJEXT) ALPS_analyt.$(OBJEXT) ALPS_com.$(OBJEXT) ALPS_fns.$(OBJEXT) ALPS_nyquist_driver.$(OBJEXT)
