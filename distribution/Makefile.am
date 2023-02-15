AM_FCFLAGS = -O4 -fcheck=all -fbacktrace -fbounds-check -ffpe-trap=zero,overflow -Wconversion -Wsurprising -ffixed-line-length-none -ffree-line-length-none -Wunused -g -fbacktrace

generate_distribution_LDADD = -lm -llapack -lblas

bin_PROGRAMS = generate_distribution

generate_distribution_SOURCES = generate_distribution.f90
