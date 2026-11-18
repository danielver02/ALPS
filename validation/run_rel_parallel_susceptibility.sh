#!/usr/bin/env bash
# Paired Task-3D harnesses: identical production mathematics and completed VDF;
# scratch instrumentation changes only the method-2 parallel derivative.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the repository' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-parallel-susceptibility.XXXXXX")
printf 'Paired susceptibility artifacts: %s\n' "$scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/*.f90 "$scratch/source/"
cp "$repo/validation/rel_susceptibility_reference.f90" "$repo/validation/make_rel_susceptibility_case.f90" \
   "$repo/validation/rel_parallel_diagnostics.f90" "$repo/validation/compare_baselines.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -fbacktrace -ffree-line-length-none -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}; mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
cd "$scratch/build"
"$fc" "${flags[@]}" -c ../source/ALPS_rel_cache.f90 ../source/ALPS_var.f90 ../source/ALPS_io.f90 \
    ../source/rel_susceptibility_reference.f90 ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 \
    ../source/ALPS_NHDS.f90 ../source/rel_parallel_diagnostics.f90
common=(ALPS_rel_cache.o ALPS_var.o ALPS_io.o rel_susceptibility_reference.o ALPS_analyt.o ALPS_com.o ALPS_NHDS.o rel_parallel_diagnostics.o)
"$fc" "${flags[@]}" -o make_case ../source/make_rel_susceptibility_case.f90
"$fc" "${flags[@]}" -o compare_mpi ../source/compare_baselines.f90
for derivative in fd analytic; do
    awk -v target=landau -v derivative="$derivative" -f "$repo/validation/instrument_rel_analytic_parallel.awk" \
        "$scratch/source/ALPS_fns_rel.f90" > "$scratch/source/landau_$derivative.f90"
    awk -v target=driver -v derivative="$derivative" -f "$repo/validation/instrument_rel_analytic_parallel.awk" \
        "$repo/validation/test_rel_susceptibility.f90" > "$scratch/source/driver_$derivative.f90"
    "$fc" "${flags[@]}" -c "$scratch/source/landau_$derivative.f90" -o ALPS_fns_rel.o
    "$fc" "${flags[@]}" -c ../source/ALPS_fns.f90
    "$fc" "${flags[@]}" -o "susceptibility_$derivative" "${common[@]}" ALPS_fns_rel.o ALPS_fns.o \
        "$scratch/source/driver_$derivative.f90" "${libs[@]}"
    for kind in iso aniso; do
        for np in 1 4; do
            name="${kind}_${derivative}_np${np}"
            mkdir -p "$scratch/$name/distribution" "$scratch/$name/solution"
            (
                cd "$scratch/$name"
                ../build/make_case "$repo" "$kind" 50 200
                study=all; [[ $kind != iso ]] || study=smoke
                timeout 900 "$mpi" "${mpi_flags[@]}" -np "$np" "../build/susceptibility_$derivative" case.in "$kind" "$study" >stdout.log 2>stderr.log
                rg -q 'PASS: finite-frequency observations' stdout.log
                rg -q 'PASS: underflow/unavailable/invalid/overflow' stdout.log
            )
            echo "PASS: $name"
        done
        ./compare_mpi "$scratch/${kind}_${derivative}_np1/observations.dat" \
            "$scratch/${kind}_${derivative}_np4/observations.dat" > "$scratch/${kind}_${derivative}_mpi.txt"
    done
done
for kind in iso aniso; do
    for np in 1 4; do
        for ((rank=0;rank<np;rank++)); do
            cmp "$scratch/${kind}_fd_np${np}/gamma_rank${rank}.csv" "$scratch/${kind}_analytic_np${np}/gamma_rank${rank}.csv"
        done
        for derivative in fd analytic; do
            rg '^(analytic|functional)_' "$scratch/${kind}_${derivative}_np${np}/observations.dat" \
                > "$scratch/${kind}_${derivative}_np${np}/supported.dat"
        done
        cmp "$scratch/${kind}_fd_np${np}/supported.dat" "$scratch/${kind}_analytic_np${np}/supported.dat"
    done
    awk -f "$repo/validation/summarize_rel_parallel.awk" "$scratch/${kind}_fd_np1/tensor.csv" \
        "$scratch/${kind}_analytic_np1/tensor.csv" > "$scratch/${kind}_paired_tensor.csv"
done
(cd "$scratch/source" && sha256sum *.f90) > "$scratch/source.sha256"
echo 'PASS: Gamma derivatives bit-identical for paired runs on every rank'
echo 'PASS: analytic-shape/functional control records byte-identical between paired runs'
printf 'PASS: paired fixed-frequency results and actual-query diagnostics in %s\n' "$scratch"
