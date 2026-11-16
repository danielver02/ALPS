#!/usr/bin/env bash
# Actual Landau routine, known row functions, inherited MPI initialization only.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the repository' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-analytic-parallel.XXXXXX")
printf 'Analytic parallel derivative artifacts: %s\n' "$scratch"
mkdir "$scratch/source" "$scratch/build" "$scratch/before"
cp "$repo"/src/ALPS_{var,io,analyt,com,rel_cache,fns_rel}.f90 "$scratch/source/"
cp "$repo/validation/rel_susceptibility_reference.f90" "$repo/validation/test_rel_analytic_parallel.f90" "$scratch/source/"
if [[ -n ${ALPS_PARALLEL_BEFORE_SOURCE:-} ]]; then
    cp "$ALPS_PARALLEL_BEFORE_SOURCE" "$scratch/source/ALPS_fns_rel_before.f90"
else
    git -C "$repo" show HEAD:src/ALPS_fns_rel.f90 > "$scratch/source/ALPS_fns_rel_before.f90"
fi
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -fbacktrace -ffree-line-length-none -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}; mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
(cd "$scratch/source" && sha256sum *.f90) > "$scratch/source.sha256"
cd "$scratch/build"
"$fc" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/rel_susceptibility_reference.f90 \
    ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 ../source/ALPS_rel_cache.f90 ../source/ALPS_fns_rel.f90
objects=(ALPS_var.o ALPS_io.o rel_susceptibility_reference.o ALPS_analyt.o ALPS_com.o ALPS_rel_cache.o)
"$fc" "${flags[@]}" -o derivative_test "${objects[@]}" ALPS_fns_rel.o ../source/test_rel_analytic_parallel.f90 "${libs[@]}"
for np in 1 2 4; do
    mkdir "$scratch/np$np"
    (
        cd "$scratch/np$np"
        timeout 120 "$mpi" "${mpi_flags[@]}" -np "$np" ../build/derivative_test >stdout.log 2>stderr.log
        cat stdout.log
        rg -q 'PASS: analytic parallel derivative with unchanged Gamma stencil and integration' stdout.log
    )
done
# Compare supported non-Chebyshev modes with the pre-change routine, byte for byte.
cd "$scratch/before"
"$fc" "${flags[@]}" -I../build -c ../source/ALPS_fns_rel_before.f90 -o ALPS_fns_rel.o
baseline_objects=()
for object in "${objects[@]}"; do baseline_objects+=("../build/$object"); done
"$fc" "${flags[@]}" -I../build -o before_test "${baseline_objects[@]}" ALPS_fns_rel.o ../source/test_rel_analytic_parallel.f90 "${libs[@]}"
mkdir run
(cd run && timeout 120 "$mpi" "${mpi_flags[@]}" -np 1 ../before_test legacy >stdout.log 2>stderr.log)
cmp run/legacy.dat ../np1/legacy.dat
printf 'PASS: %s supported-method Landau records byte-identical before/after\n' "$(wc -l < run/legacy.dat)"
mkdir "$scratch/fatal"
(
    cd "$scratch/fatal"
    set +e
    timeout 60 "$mpi" "${mpi_flags[@]}" -np 2 ../build/derivative_test fatal_derivative >stdout.log 2>stderr.log
    status=$?
    set -e
    [[ $status -ne 0 && $status -ne 124 ]]
    rg -q 'ERROR: relativistic Chebyshev continuation:' stdout.log
    rg -q 'igamma=1 status=10' stdout.log
    if rg -q 'SIGFPE|Floating point exception' stderr.log; then exit 1; fi
    echo 'PASS: finite value with unrepresentable analytic derivative fails through explicit continuation status'
)
printf 'PASS: analytic parallel integration tests; artifacts remain in %s\n' "$scratch"
