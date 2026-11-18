#!/usr/bin/env bash
# Task 3D scientific observations, not a baseline capture or a root search.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case "${1:-}" in '') mode=all;; --smoke) mode=smoke;; --queries) mode=queries;; --case) [[ $# == 4 ]] || exit 2; mode=case; case_kind=$2; case_ng=$3; case_nu=$4;; *) echo 'Usage: bash validation/run_rel_susceptibility.sh [--smoke|--queries|--case {iso,aniso,coarse} ngamma npparbar]' >&2; exit 2;; esac
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the source tree' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-susceptibility.XXXXXX")
printf 'Fixed-frequency susceptibility artifacts: %s\n' "$scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/*.f90 "$scratch/source/"
cp "$repo/validation/rel_susceptibility_reference.f90" "$repo/validation/make_rel_susceptibility_case.f90" \
   "$repo/validation/test_rel_susceptibility.f90" "$repo/validation/compare_baselines.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=bounds -fbacktrace -ffree-line-length-none -ffpe-trap=overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}; mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
{
    git -C "$repo" rev-parse HEAD
    git -C "$repo" status --short
    "$fc" --version
    printf 'Flags: %s\nLibraries: %s\n' "${flags[*]}" "${libs[*]}"
} > "$scratch/provenance.txt"
(cd "$scratch/source" && sha256sum *.f90) > "$scratch/source.sha256"
cd "$scratch/build"
"$fc" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/rel_susceptibility_reference.f90 \
    ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 ../source/ALPS_rel_cache.f90 ../source/ALPS_NHDS.f90 \
    ../source/ALPS_fns_rel.f90 ../source/ALPS_fns.f90
objects=(ALPS_var.o ALPS_io.o rel_susceptibility_reference.o ALPS_analyt.o ALPS_com.o ALPS_rel_cache.o ALPS_NHDS.o ALPS_fns_rel.o ALPS_fns.o)
"$fc" "${flags[@]}" -o susceptibility "${objects[@]}" ../source/test_rel_susceptibility.f90 "${libs[@]}"
"$fc" "${flags[@]}" -o make_case ../source/make_rel_susceptibility_case.f90
"$fc" "${flags[@]}" -o compare_mpi ../source/compare_baselines.f90
ldd susceptibility >> "$scratch/provenance.txt"
run_case() {
    local kind=$1 ng=$2 nu=$3 study=$4 np=$5 name
    name="${kind}_g${ng}_u${nu}_np${np}"
    mkdir -p "$scratch/$name/distribution" "$scratch/$name/solution"
    (
        cd "$scratch/$name"
        ../build/make_case "$repo" "$kind" "$ng" "$nu"
        timeout 900 "$mpi" "${mpi_flags[@]}" -np "$np" ../build/susceptibility case.in "$kind" "$study" >stdout.log 2>stderr.log
        rg -q 'PASS: finite-frequency observations' stdout.log
        rg -q 'PASS: underflow/unavailable/invalid/overflow' stdout.log
        printf 'Completed %s\n' "$name"
    )
}
if [[ $mode == case ]]; then
    run_case "$case_kind" "$case_ng" "$case_nu" queries 4
elif [[ $mode == queries ]]; then
    run_case aniso 50 100 queries 4
elif [[ $mode == smoke ]]; then
    run_case iso 50 100 smoke 1
    run_case iso 50 100 smoke 4
    "$scratch/build/compare_mpi" "$scratch/iso_g50_u100_np1/observations.dat" "$scratch/iso_g50_u100_np4/observations.dat"
else
    for kind in iso aniso coarse; do
        run_case "$kind" 50 100 all 4
        for grid in '100 100' '200 100' '50 200' '50 400' '100 200' '200 400'; do
            read -r ng nu <<< "$grid"
            run_case "$kind" "$ng" "$nu" refine 4
        done
        run_case "$kind" 50 100 all 1
        "$scratch/build/compare_mpi" "$scratch/${kind}_g50_u100_np1/observations.dat" \
            "$scratch/${kind}_g50_u100_np4/observations.dat" > "$scratch/${kind}_mpi_comparison.txt"
        cat "$scratch/${kind}_mpi_comparison.txt"
    done
fi
for file in tensor summary rows; do
    first=true
    for path in "$scratch"/*np4/"$file.csv"; do
        if $first; then cat "$path"; first=false; else tail -n +2 "$path"; fi
    done > "$scratch/$file.csv"
done
printf 'Complete: all complex components, errors and row diagnostics in %s/{tensor,summary,rows}.csv\n' "$scratch"
