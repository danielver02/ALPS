#!/usr/bin/env bash
# Task 3C: fresh, isolated initialization/transfer/continuation tests.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the source tree' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-chebyshev-mpi.XXXXXX")
printf 'Relativistic Chebyshev MPI artifacts: %s\n' "$scratch"
mkdir "$scratch/source"
cp "$repo/src/ALPS_var.f90" "$repo/src/ALPS_io.f90" "$repo/src/ALPS_analyt.f90" \
   "$repo/src/ALPS_com.f90" "$repo/src/ALPS_rel_cache.f90" "$repo/src/ALPS_fns_rel.f90" \
   "$repo/distribution/distribution_analyt.f90" "$repo/tests/test_rel_chebyshev_mpi.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -fbacktrace -ffree-line-length-none -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}; mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
cd "$scratch"
"$fc" "${flags[@]}" -c source/ALPS_var.f90 source/ALPS_io.f90 source/distribution_analyt.f90 \
    source/ALPS_analyt.f90 source/ALPS_com.f90 source/ALPS_rel_cache.f90 source/ALPS_fns_rel.f90
"$fc" "${flags[@]}" -o mpi_rows ALPS_var.o ALPS_io.o distribution_analyt.o ALPS_analyt.o \
    ALPS_com.o ALPS_rel_cache.o ALPS_fns_rel.o source/test_rel_chebyshev_mpi.f90 "${libs[@]}"
for np in 1 2 4; do
    mkdir -p "np$np/distribution"
    (
        cd "np$np"
        timeout 180 "$mpi" "${mpi_flags[@]}" -np "$np" ../mpi_rows >stdout.log 2>stderr.log
        cat stdout.log
        rg -q 'PASS: MPI construction, distribution, continuation and lifecycle' stdout.log
    )
done
mkdir -p fatal/distribution
(
    cd fatal
    set +e
    timeout 60 "$mpi" "${mpi_flags[@]}" -np 2 ../mpi_rows fatal >stdout.log 2>stderr.log
    status=$?
    set -e
    [[ $status -ne 0 && $status -ne 124 ]]
    rg -q 'ERROR: relativistic Chebyshev continuation:' stdout.log
    rg -q 'species=4 igamma=3 status=15' stdout.log
    echo 'PASS: unavailable continuation without status aborts the communicator'
)
printf 'PASS: 1/2/4-rank row integration; all artifacts remain in %s\n' "$scratch"
