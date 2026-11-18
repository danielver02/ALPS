#!/usr/bin/env bash
# Negative control at the SAME anisotropic parameters: the original generic
# seed reaches a longitudinal candidate, not the selected transverse branch.
# Reuse the built driver; preserve its rejected-reference preflight evidence.
set -euo pipefail
[[ $# == 1 ]] || { echo 'Usage: bash validation/run_rel_root_failure_check.sh ROOT_RUN_DIRECTORY' >&2;exit 2; }
scratch=$(realpath "$1")
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case "$scratch/" in "$repo/"*) echo 'Root run must be outside the repository' >&2;exit 2;; esac
[[ -x $scratch/build/roots && -f $scratch/aniso_np4_method0/root.in ]]
mkdir "$scratch/rejected_initial_seed"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
mpi=${MPIEXEC:-mpirun}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
printf 'method,exit_status,interpretation\n' > "$scratch/rejected_initial_seed/process_status.csv"
for method in 0 1 2; do
 dir="$scratch/rejected_initial_seed/method$method"
 mkdir -p "$dir/distribution" "$dir/solution"
 cp "$scratch/aniso_np4_method0/root.in" "$dir/"
 cp "$scratch/aniso_np4_method0/distribution/"*.array "$dir/distribution/"
 if [[ $method != 0 ]];then cp "$scratch/rejected_initial_seed/method0/reference.dat" "$dir/";fi
 (
  cd "$dir"
  set +e
  timeout 240 "$mpi" "${mpi_flags[@]}" -np 4 "$scratch/build/roots" root.in aniso "$method" >stdout.log 2>stderr.log
  status=$?
  set -e
  if [[ $method == 0 ]];then
   [[ $status == 0 ]]
   awk -F, 'NR==2 {if($11!=1 || $16!=0 || $17<.99) exit 1;found=1} END {if(!found || NR!=2) exit 1}' roots.csv
   printf '0,%s,rejected_longitudinal_candidate_no_scan\n' "$status" >> ../process_status.csv
  else
   [[ $status -ne 0 && $status -ne 124 ]]
   rg -q 'unqualified initial reference root: preflight recorded, recovery withheld' stderr.log
   [[ $(wc -l < tensor.csv) == 13 && $(wc -l < roots.csv) == 1 ]]
   if rg -q 'SIGFPE|Floating point exception' stderr.log;then exit 1;fi
   printf '%s,%s,preflight_recorded_recovery_withheld\n' "$method" "$status" >> ../process_status.csv
  fi
 )
done
cat "$scratch/rejected_initial_seed/process_status.csv"
echo 'PASS: small-step termination is not accepted as reference convergence; rejected-reference tensors retained'
