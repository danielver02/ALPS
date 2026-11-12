#!/usr/bin/env bash
# Task 1: the established generate_distribution -> MPI ALPS workflow, isolated.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
mode=${1:---check}
case "$mode" in --check|--capture|--convergence) ;; *) echo 'Usage: bash tests/run_numerical_baselines.sh [--check|--capture|--convergence]' >&2; exit 2;; esac
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the source tree' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-numerical.XXXXXX")
printf 'Task 1 artifacts: %s\n' "$scratch"
mkdir "$scratch/build" "$scratch/source" "$scratch/seed"
mkdir "$scratch/source/src" "$scratch/source/distribution" "$scratch/source/tests"
# Copy only consumed sources/fixtures. No prebuilt .mod/.o, cache, or result enters the build.
cp "$repo"/src/*.f90 "$scratch/source/src/"
cp "$repo/distribution/distribution_analyt.f90" "$repo/distribution/generate_distribution.f90" "$scratch/source/distribution/"
cp "$repo"/tests/*.f90 "$scratch/source/tests/"
cp "$repo/tests/run_numerical_baselines.sh" "$scratch/source/tests/"
cp "$repo/tests/test_relativistic.in" "$repo/tests/test_chebyshev.in" "$scratch/source/tests/"
cp "$repo/distribution/test_relativistic_dist.in" "$repo/distribution/test_bimax_dist.in" "$scratch/source/distribution/"
cp "$scratch/source/distribution/"*.in "$scratch/seed/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=bounds -fbacktrace -ffree-line-length-none -ffpe-trap=overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--lopenblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}
mpi=${MPIEXEC:-mpirun}
# Reproducible default BLAS reduction order; the chosen values are recorded.
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
{
    date -u '+UTC %Y-%m-%dT%H:%M:%SZ'
    git -C "$repo" rev-parse HEAD
    git -C "$repo" status --short
    "$fc" --version
    "$fc" --showme || true
    "$mpi" --version
    pkg-config --modversion openblas || true
    printf 'Compiler: %s\nFlags: %s\nLibraries: %s\nMPI flags: %s\nRanks: %s\n' "$fc" "${flags[*]}" "${libs[*]}" "${mpi_flags[*]}" "${ALPS_TEST_NP:-4}"
    printf 'OPENBLAS_NUM_THREADS=%s\nOMP_NUM_THREADS=%s\n' "$OPENBLAS_NUM_THREADS" "$OMP_NUM_THREADS"
    env | sort | rg '^(MKL_|BLIS_|GOTO_|OPENBLAS_|OMP_|VECLIB_)' || true
} > "$scratch/provenance.txt"
git -C "$repo" diff --binary > "$scratch/source.patch"
(cd "$scratch/source" && find . -type f -print0 | sort -z | xargs -0 sha256sum) > "$scratch/source.sha256"
cd "$scratch/build"
sources=(ALPS_rel_cache ALPS_var ALPS_io ALPS_analyt ALPS_com ALPS_NHDS ALPS_fns_rel ALPS_fns)
"$fc" "${flags[@]}" -c ../source/src/ALPS_rel_cache.f90 ../source/src/ALPS_var.f90 ../source/distribution/distribution_analyt.f90
for source in "${sources[@]:2}"; do
    "$fc" "${flags[@]}" -c "../source/src/$source.f90"
done
objects=(ALPS_rel_cache.o ALPS_var.o distribution_analyt.o ALPS_io.o ALPS_com.o ALPS_analyt.o ALPS_NHDS.o ALPS_fns_rel.o ALPS_fns.o)
"$fc" "${flags[@]}" -c ../source/tests/relativistic_references.f90
"$fc" "${flags[@]}" -o ALPS "${objects[@]}" ../source/src/ALPS.f90 "${libs[@]}"
"$fc" "${flags[@]}" -o numerical_baselines "${objects[@]}" relativistic_references.o ../source/tests/test_numerical_baselines.f90 "${libs[@]}"
"$fc" "${flags[@]}" -o coordinate_contract relativistic_references.o ../source/tests/test_coordinate_contract.f90
"$fc" "${flags[@]}" -o make_case relativistic_references.o ../source/tests/make_baseline_case.f90
"$fc" "${flags[@]}" -o compare_baselines ../source/tests/compare_baselines.f90
"$fc" "${flags[@]}" -o compare_convergence ../source/tests/compare_convergence.f90
"$fc" "${flags[@]}" -o comparison_fixtures ../source/tests/test_comparison_fixtures.f90
"$fc" "${flags[@]}" -o generate_distribution distribution_analyt.o ../source/distribution/generate_distribution.f90 "${libs[@]}"
ldd ALPS >> "$scratch/provenance.txt"
sha256sum ALPS numerical_baselines >> "$scratch/provenance.txt"
mkdir "$scratch/coordinates"
(
    cd "$scratch/coordinates"
    ../build/coordinate_contract
    ../build/comparison_fixtures
    ../build/compare_baselines reader_reference.dat reader_reference.dat
    for bad in changed nonfinite empty; do
        if ../build/compare_baselines reader_reference.dat "reader_$bad.dat" >"reader_$bad.log" 2>&1; then
            echo "FAIL: numerical comparator accepted $bad data" >&2
            exit 1
        fi
    done
    echo 'PASS: numerical reader rejects changed, nonfinite, and missing records'
)
(cd "$scratch/seed" && ../build/generate_distribution test_relativistic_dist.in >relativistic.log 2>&1)
(cd "$scratch/seed" && ../build/generate_distribution test_bimax_dist.in >bimax.log 2>&1)
run_case() {
    local name=$1 kind=$2 ng=$3 nu=$4 nr=$5
    mkdir -p "$scratch/$name/distribution" "$scratch/$name/solution"
    cp "$scratch/seed/"*.array "$scratch/$name/distribution/"
    (
        cd "$scratch/$name"
        ../build/make_case "$scratch/source" "$kind" "$ng" "$nu" "$nr"
        timeout 300 "$mpi" "${mpi_flags[@]}" -np "${ALPS_TEST_NP:-4}" ../build/numerical_baselines case.in "$kind" >stdout.log 2>stderr.log
        rg -q 'PASS: finite numerical observations' stdout.log
        # A Fortran STOP in legacy fitting can exit successfully; completion is mandatory.
        ! rg -qi 'Fortran runtime error|out of bounds|Fit matrix.*singular|Fit matrix inversion failed' stdout.log stderr.log
        printf '%s: ' "$name"
        if [[ $mode == --check ]]; then
            ../build/compare_baselines "$repo/tests/data/relativistic_baselines/$kind.dat" observations.dat
        else
            echo 'captured; observations and root status retained'
        fi
        cat root_status.dat
    )
}
for kind in rel4 rel5 nr_fd nr_analytic; do
    run_case "$kind" "$kind" 100 200 100
done
if [[ $mode == --convergence ]]; then
    for kind in rel4 rel5; do
        run_case "${kind}_grid2" "$kind" 200 400 100
        run_case "${kind}_grid3" "$kind" 500 500 100
        run_case "${kind}_res2" "$kind" 100 200 200
        run_case "${kind}_res4" "$kind" 100 200 400
    done
    {
        for kind in rel4 rel5; do
            for refined in grid2 grid3 res2 res4; do
                printf '%s -> %s_%s\n' "$kind" "$kind" "$refined"
                "$scratch/build/compare_convergence" "$scratch/$kind/observations.dat" "$scratch/${kind}_$refined/observations.dat"
            done
            printf '%s_grid2 -> %s_grid3\n' "$kind" "$kind"
            "$scratch/build/compare_convergence" "$scratch/${kind}_grid2/observations.dat" "$scratch/${kind}_grid3/observations.dat"
        done
    } > "$scratch/convergence.txt"
    cat "$scratch/convergence.txt"
fi
printf 'Complete. All generated files are under %s\n' "$scratch"
