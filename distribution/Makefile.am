AM_FCFLAGS = -O4 -fcheck=all -fbacktrace -fbounds-check -ffpe-trap=zero,overflow -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -g -fbacktrace -I../src

generate_distribution_LDADD = -lm -llapack -lblas

bin_PROGRAMS = generate_distribution

generate_distribution_SOURCES = distribution_analyt.f90 generate_distribution.f90
