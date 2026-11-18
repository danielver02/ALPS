#!/usr/bin/env bash
# Serial kernel tests; mpifort builds the legacy dependencies of ALPS_analyt.
# No MPI initialization/launcher or existing binaries used. The transformed-row
# fixture writes distributions/cache entries only in this run's scratch directory.
# Overrides: MPIFC, FFLAGS, ALPS_TEST_LIBS, TMPDIR (must be outside the repo).
# --matrix selects O0/O2 with traps off/on (overrides FFLAGS for these runs).
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case "${1:-}" in
    --matrix)
        for optimization in -O0 -O2; do
            for traps in off on; do
                flags="$optimization -g -fcheck=all -fbacktrace -ffree-line-length-none"
                [[ $traps != on ]] || flags+=" -ffpe-trap=invalid,zero,overflow"
                printf '\nConfiguration: %s, traps %s\n' "$optimization" "$traps"
                FFLAGS="$flags" bash "$repo/tests/run_chebyshev_kernel.sh"
            done
        done
        exit 0;;
    '') ;;
    *) echo 'Usage: bash tests/run_chebyshev_kernel.sh [--matrix]' >&2; exit 2;;
esac
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the source tree' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-chebyshev.XXXXXX")
printf 'Chebyshev kernel artifacts: %s\n' "$scratch"
mkdir "$scratch/source"
cp "$repo/src/ALPS_var.f90" "$repo/src/ALPS_io.f90" "$repo/src/ALPS_analyt.f90" "$scratch/source/"
cp "$repo/distribution/distribution_analyt.f90" "$repo/tests/test_chebyshev_kernel.f90" \
   "$repo/tests/test_chebyshev_log.f90" "$repo/tests/test_chebyshev_lapack_failure.f90" "$scratch/source/"
cp "$repo/tests/relativistic_references.f90" "$repo/tests/test_rel_chebyshev_rows.f90" \
   "$repo/tests/test_rel_chebyshev_transformed.f90" "$repo/src/ALPS_rel_cache.f90" \
   "$repo/src/ALPS_fns_rel.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -fbacktrace -ffree-line-length-none -ffpe-trap=invalid,zero,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
fc=${MPIFC:-mpifort}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
cd "$scratch"
"$fc" "${flags[@]}" -c source/ALPS_var.f90 source/ALPS_io.f90 source/distribution_analyt.f90 source/ALPS_analyt.f90
objects=(ALPS_var.o ALPS_io.o distribution_analyt.o ALPS_analyt.o)
# The real-library executable verifies DGELSD availability and actual SVD fits.
"$fc" "${flags[@]}" -o kernel_test "${objects[@]}" source/test_chebyshev_kernel.f90 "${libs[@]}"
./kernel_test
"$fc" "${flags[@]}" -o log_test "${objects[@]}" source/test_chebyshev_log.f90 "${libs[@]}"
./log_test
"$fc" "${flags[@]}" -o lapack_failure_test "${objects[@]}" source/test_chebyshev_lapack_failure.f90 "${libs[@]}"
./lapack_failure_test
"$fc" "${flags[@]}" -c source/relativistic_references.f90
"$fc" "${flags[@]}" -o row_test "${objects[@]}" relativistic_references.o source/test_rel_chebyshev_rows.f90 "${libs[@]}"
./row_test
"$fc" "${flags[@]}" -c source/ALPS_rel_cache.f90 source/ALPS_fns_rel.f90
"$fc" "${flags[@]}" -o transformed_row_test "${objects[@]}" ALPS_rel_cache.o ALPS_fns_rel.o \
    relativistic_references.o source/test_rel_chebyshev_transformed.f90 "${libs[@]}"
mkdir distribution cache
./transformed_row_test
printf 'PASS: serial kernel and LAPACK fault tests; generated files remain in %s\n' "$scratch"
