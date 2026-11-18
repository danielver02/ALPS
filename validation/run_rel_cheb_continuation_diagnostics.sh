#!/usr/bin/env bash
# Fresh snapshots; existing staged/unstaged tests and all production code stay unchanged.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case "${1:-}" in '') control=false;; --with-isotropic) control=true;; *) echo 'Usage: bash validation/run_rel_cheb_continuation_diagnostics.sh [--with-isotropic]' >&2; exit 2;; esac
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the repository' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-cheb-continuation.XXXXXX")
printf 'Continuation diagnostic artifacts: %s\n' "$scratch"
mkdir -p "$scratch/source" "$scratch/build" "$scratch/fixture/tests"
cp "$repo"/src/*.f90 "$scratch/source/"
cp "$repo/validation/rel_susceptibility_reference.f90" "$repo/validation/make_rel_susceptibility_case.f90" \
   "$repo/validation/rel_cheb_continuation_diagnostics.f90" "$scratch/source/"
cp "$repo/validation/test_rel_susceptibility.f90" "$scratch/source/original_driver.f90"
cp "$repo/tests/test_relativistic.in" "$scratch/fixture/tests/"
awk -f "$repo/validation/instrument_rel_cheb_continuation.awk" "$scratch/source/original_driver.f90" > "$scratch/source/diagnostic_driver.f90"
diff -u "$scratch/source/original_driver.f90" "$scratch/source/diagnostic_driver.f90" > "$scratch/instrumentation.patch" || [[ $? == 1 ]]
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=bounds -fbacktrace -ffree-line-length-none -ffpe-trap=overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort}; mpi=${MPIEXEC:-mpirun}; np=${ALPS_TEST_NP:-4}; nu=${ALPS_DIAG_NPPAR:-200}
case "$nu" in 100|200) ;; *) echo 'ALPS_DIAG_NPPAR must be 100 or 200' >&2; exit 2;; esac
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
{
    git -C "$repo" rev-parse HEAD
    git -C "$repo" status --short
    "$fc" --version
    "$mpi" --version
    printf 'Flags: %s\nLibraries: %s\nMPI ranks: %s\nParallel intervals: %s\nBLAS/OMP threads: %s/%s\n' "${flags[*]}" "${libs[*]}" "$np" "$nu" "$OPENBLAS_NUM_THREADS" "$OMP_NUM_THREADS"
} > "$scratch/provenance.txt"
(cd "$scratch" && sha256sum source/*.f90 fixture/tests/*.in) > "$scratch/source.sha256"
cd "$scratch/build"
"$fc" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/rel_susceptibility_reference.f90 \
    ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 ../source/ALPS_rel_cache.f90 ../source/ALPS_NHDS.f90 \
    ../source/ALPS_fns_rel.f90 ../source/ALPS_fns.f90 ../source/rel_cheb_continuation_diagnostics.f90
objects=(ALPS_var.o ALPS_io.o rel_susceptibility_reference.o ALPS_analyt.o ALPS_com.o ALPS_rel_cache.o ALPS_NHDS.o ALPS_fns_rel.o ALPS_fns.o rel_cheb_continuation_diagnostics.o)
"$fc" "${flags[@]}" -o diagnostics "${objects[@]}" ../source/diagnostic_driver.f90 "${libs[@]}"
"$fc" "${flags[@]}" -o make_case ../source/make_rel_susceptibility_case.f90
ldd diagnostics >> "$scratch/provenance.txt"
cases=(aniso); if $control; then cases+=(iso); fi
for case_kind in "${cases[@]}"; do
    mkdir -p "$scratch/$case_kind/distribution" "$scratch/$case_kind/solution"
    (
        cd "$scratch/$case_kind"
        ../build/make_case "$scratch/fixture" "$case_kind" 50 "$nu"
        timeout 900 "$mpi" "${mpi_flags[@]}" -np "$np" ../build/diagnostics case.in "$case_kind" queries > stdout.log 2>stderr.log
        rg -q 'PASS: coefficient/query/truncation observations' stdout.log
        rg -q 'PASS: finite-frequency observations' stdout.log
        printf 'PASS: %s diagnostics and selected fixed-frequency susceptibilities\n' "$case_kind"
    )
done
printf 'CSV results and instrumented source: %s\n' "$scratch"
printf 'Plots: gnuplot -e "data_dir=\x27%s/aniso\x27; out_dir=\x27%s\x27" %s/validation/plot_rel_cheb_continuation.gp\n' "$scratch" "$scratch" "$repo"
