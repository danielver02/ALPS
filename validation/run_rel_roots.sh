#!/usr/bin/env bash
# Task 3E: controlled attempts, never a baseline capture or production activation.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case "${1:-}" in '') mode=all;; --smoke) mode=smoke;; *) echo 'Usage: bash validation/run_rel_roots.sh [--smoke]' >&2;exit 2;; esac
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the repository' >&2;exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-rel-roots.XXXXXX")
printf 'Controlled relativistic root artifacts: %s\n' "$scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/*.f90 "$scratch/source/"
cp "$repo/validation/rel_susceptibility_reference.f90" "$repo/validation/make_rel_susceptibility_case.f90" \
   "$repo/validation/rel_root_diagnostics.f90" "$repo/validation/rel_root_bridge.f90" "$repo/validation/test_rel_roots.f90" \
   "$repo/validation/compare_baselines.f90" "$scratch/source/"
awk -f "$repo/validation/instrument_rel_root.awk" "$scratch/source/ALPS_fns.f90" > "$scratch/source/ALPS_fns_observed.f90"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -fbacktrace -ffree-line-length-none -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
fc=${MPIFC:-mpifort};mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
{
 git -C "$repo" rev-parse HEAD
 git -C "$repo" status --short
 "$fc" --version
 printf 'Flags: %s\nLibraries: %s\nMPI flags: %s\n' "${flags[*]}" "${libs[*]}" "${mpi_flags[*]}"
} > "$scratch/provenance.txt"
(cd "$scratch/source" && sha256sum *.f90) > "$scratch/source.sha256"
cd "$scratch/build"
"$fc" "${flags[@]}" -c ../source/ALPS_rel_cache.f90 ../source/ALPS_var.f90 ../source/ALPS_io.f90 \
 ../source/rel_susceptibility_reference.f90 ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 \
 ../source/ALPS_NHDS.f90 ../source/ALPS_fns_rel.f90 ../source/rel_root_diagnostics.f90 ../source/ALPS_fns_observed.f90
objects=(ALPS_rel_cache.o ALPS_var.o ALPS_io.o rel_susceptibility_reference.o ALPS_analyt.o ALPS_com.o ALPS_NHDS.o ALPS_fns_rel.o rel_root_diagnostics.o ALPS_fns_observed.o)
"$fc" "${flags[@]}" -o roots "${objects[@]}" ../source/rel_root_bridge.f90 ../source/test_rel_roots.f90 "${libs[@]}"
"$fc" "${flags[@]}" -o make_case ../source/make_rel_susceptibility_case.f90
"$fc" "${flags[@]}" -o compare_mpi ../source/compare_baselines.f90
printf 'case,np,method,exit_status,completed\n' > "$scratch/process_status.csv"
kinds=(iso aniso);ranks=(1 4);methods=(0 1 2)
if [[ $mode == smoke ]]; then kinds=(iso);ranks=(4);methods=(0);fi
for kind in "${kinds[@]}"; do
 for np in "${ranks[@]}"; do
  for method in "${methods[@]}"; do
   name="${kind}_np${np}_method${method}"
   mkdir -p "$scratch/$name/distribution" "$scratch/$name/solution"
   (
    cd "$scratch/$name"
    ../build/make_case "$repo" "$kind" "${ALPS_ROOT_NG:-50}" "${ALPS_ROOT_NU:-200}"
    # Established Task-1 rtsec settings; no solver algorithm/tolerance changes.
    awk '/^secant_method=/{ $0="secant_method=1" } /^numiter[ =]/{ $0="numiter=150" } /^D_tol=/{ $0="D_tol=1.d-11" } {print}' case.in > root.in
    if [[ $kind == aniso && $method == 0 && -f ../iso_np${np}_method0/reference.dat ]]; then
     head -n 1 "../iso_np${np}_method0/reference.dat" > seed.dat
    fi
    if [[ $method != 0 ]]; then
     if [[ ! -f ../${kind}_np${np}_method0/reference.dat ]]; then
      printf '%s,%s,%s,99,0\n' "$kind" "$np" "$method" >> "$scratch/process_status.csv"
      echo 'MISSING REFERENCE: method not attempted' > stdout.log
      exit 0
     fi
     cp "../${kind}_np${np}_method0/reference.dat" reference.dat
    fi
    set +e
    timeout "${ALPS_ROOT_TIMEOUT:-900}" "$mpi" "${mpi_flags[@]}" -np "$np" ../build/roots root.in "$kind" "$method" >stdout.log 2>stderr.log
    status=$?
    set -e
    completed=0
    if [[ $status == 0 ]] && rg -q '^COMPLETE: controlled root attempts' stdout.log; then completed=1;fi
    printf '%s,%s,%s,%s,%s\n' "$kind" "$np" "$method" "$status" "$completed" >> "$scratch/process_status.csv"
    printf 'Recorded %s: exit=%s completed=%s (acceptance requires residual and branch checks)\n' "$name" "$status" "$completed"
   )
  done
 done
done
cat "$scratch/process_status.csv"
# Summarize and compare without treating process completion as root acceptance.
bash "$repo/validation/summarize_rel_root_run.sh" "$scratch"
printf 'Complete: root/iteration/tensor/query CSVs and failure logs remain in %s\n' "$scratch"
