#################################################################### OVERVIEW
#  Makefile for ALPS:
#  (name of ALPS to be determined)
#
#   VERSION 0.0
#
#   Kristopher Klein and Daniel Verscharen
#
#   Version notes:
#
#  LAST UPDATE:  2016/05/20
###############################################################################
# SET VARIABLES FOR SYSTEM, and PROFILING
# SYSTEM options:ubik,trillian,mercer
ALPS_SYSTEM=wub
PROFILE= false
###############################################################################
###############################################################################
# An example execution:
# mpirun -np 16 ./ALPS.e sample.in
#

PACK = Makefile \
	src/*.f90 \
	*.in \
	README* \
	distribution/*.in \
	distribution/*.f90 \
	interpolation/*.f90

#SET FLAGS FOR Wub (Ubuntu 14.10 Laptop)
#MPI compiled with gfortran
ifeq ($(ALPS_SYSTEM),wub)
	COMP= mpifort
	STDCOMP:= gfortran
	STDCOMPOPTS := -fcheck=all -fbacktrace -fbounds-check -ffpe-trap=zero,overflow -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -g -fbacktrace -lm -I./include/
	COMPOPTS := -O3 -g -fbounds-check -ffpe-trap=zero,overflow -ffast-math -Wunused -funroll-loops -g -fbacktrace -lm -I./include/
	#FLAGS= -DDOUBLE
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

ifeq ($(ALPS_SYSTEM),Mac10)
COMP= mpifort-openmpi-gcc10
STDCOMP:= gfortran
STDCOMPOPTS := -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -fbounds-check -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -g -fbacktrace -lm -I./include/
COMPOPTS := -fbounds-check -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch -ffast-math -Wunused -funroll-loops -g -fbacktrace -lm -I./include/
#FLAGS= -DDOUBLE
ifeq ($(PROFILE),true)
	FLAGS += -g
endif
endif


#SET FLAGS FOR Trillian
#MPI compiled with ifort (I assume?)
ifeq ($(ALPS_SYSTEM),trillian)
	#COMP= ifort
	COMP= ftn
	STDCOMP:= gfortran
	STDCOMPOPTS := -fcheck=all -fbacktrace -fbounds-check -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -g -fbacktrace -lm -I./include/
	FLAGS= -O3 -DDOUBLE
	# This are good options for cray-ftn:
	# COMPOPTS := -e m -lm -I./include/
	# This are good options for GFortran:
	 COMPOPTS := -ffast-math -funroll-loops -lm -I./include/
	# These are good options for PG Fortran:
	# COMPOPTS :=  -fast -Mipa=fast -Mfprelaxed -lm -I./include/
	# These are good options for Intel Fortran:
	#COMPOPTS := -unroll-aggressive -opt-prefetch -fast -I./include/
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#Other systems to be added...





SRCS := $(shell ls src/*.f90)
OBJS:=  $(SRCS:src/%.f90=obj/%.o)

# Let's use different COMPOPTS depending on the system.
#COMPOPTS := -ffpe-trap=precision -fcheck=all -fbacktrace -fbounds-check -Wconversion -Wsurprising -ffpe-trap=zero -ffpe-trap=overflow -ffpe-trap=underflow -ffixed-line-length-none -ffree-line-length-none -Wunused -lm -I./include/
#COMPOPTS := -fcheck=all -fbacktrace -fbounds-check -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -lm -I./include/

VPATH= src:include:/usr/include/

###############################################################################
all: ALPS generate_distribution interpolation tidyup

ALPS: $(OBJS) | solution
	$(COMP) $(COMPOPTS) $(FLAGS) -o ALPS.e  $(LIBS) $(OBJS)

###############################################################################

###############################################################################

tidyup:	| include
	if ls *.mod 1> /dev/null 2>&1; then mv *.mod include/; fi


clean:
	rm -f obj/*.o
	rm -f include/*.mod
	rm -f ALPS*.e
	rm -f distribution/generate_distribution.e
	rm -f interpolation/interpolation.e


interpolation:interpolation/interpolation.e

generate_distribution:distribution/generate_distribution.e

tar:
	mkdir pack_ALPS_`date +'%y%m%d'`
	rsync -R $(PACK) pack_ALPS_`date +'%y%m%d'`
	tar -cvf  pack_ALPS_`date +'%y%m%d'`.tar pack_ALPS_`date +'%y%m%d'`
	rm -r pack_ALPS_`date +'%y%m%d'`
#	tar -cvf  pack_ALPS_`date +'%y%m%d'`.tar $(PACK)

map:
#	mpirun -np 8 ./ALPS.e map.in
	mpirun -np 8 ./ALPS.e map_max_C.in

guess:
	mpirun -np 4 ./ALPS.e guess.in

run:
	mpirun -np 8 ./ALPS.e guess_scan.in
#	mpirun -np 8 ./ALPS.e guess_double_scan.in
#	mpirun -np 8 ./ALPS.e guess_max.in
#	mpirun -np 4 ./ALPS.e map_scan.in

#########Rules
obj/%.o:src/%.f90
	$(COMP) $(COMPOPTS) $(FLAGS) -c  $< -o $@

interpolation/interpolation.e:interpolation/interpolation.f90
	$(STDCOMP) $(STDCOMPOPTS)  ./interpolation/interpolation.f90 -o ./interpolation/interpolation.e

distribution/generate_distribution.e:distribution/generate_distribution.f90
	$(STDCOMP) $(STDCOMPOPTS)  ./distribution/generate_distribution.f90 -o ./distribution/generate_distribution.e

#########Dependencies
obj/ALPS_var.o:
obj/ALPS_io.o:     obj/ALPS_var.o
obj/ALPS_fns_rel.o:    obj/ALPS_var.o obj/ALPS_analyt.o
obj/ALPS_fns.o:    obj/ALPS_var.o obj/ALPS_analyt.o obj/ALPS_fns_rel.o
obj/ALPS_com.o:    obj/ALPS_var.o
obj/ALPS_analyt.o: obj/ALPS_var.o obj/ALPS_io.o
obj/ALPS.o:	obj/ALPS_var.o obj/ALPS_io.o obj/ALPS_fns_rel.o obj/ALPS_fns.o obj/ALPS_com.o obj/ALPS_analyt.o

$(OBJS): | obj

obj: | include
	mkdir obj

include:
	mkdir include

solution:
	mkdir solution
