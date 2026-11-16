#!/usr/bin/env bash
# Analyze existing attempts without rerunning the solver or touching references.
set -euo pipefail
[[ $# == 1 ]] || { echo 'Usage: bash validation/summarize_rel_root_run.sh ROOT_RUN_DIRECTORY' >&2;exit 2; }
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(realpath "$1")
case "$scratch/" in "$repo/"*) echo 'Root run must be outside the repository' >&2;exit 2;; esac
first=true
for kind in iso aniso; do
 for np in 1 4; do
  for method in 0 1 2; do
   actual="$scratch/${kind}_np${np}_method${method}"
   reference="$scratch/${kind}_np${np}_method0"
   [[ -f $actual/roots.csv && -f $reference/roots.csv ]] || continue
   awk -v testcase="$kind" -v np="$np" -f "$repo/validation/summarize_rel_roots.awk" \
     "$reference/roots.csv" "$actual/roots.csv" > "$actual/comparison.csv"
   if [[ -s $actual/comparison.csv ]];then
    if $first;then cp "$actual/comparison.csv" "$scratch/root_comparisons.csv";first=false
    else tail -n +2 "$actual/comparison.csv" >> "$scratch/root_comparisons.csv";fi
   fi
   awk -F, 'NR>1 {printf "omega_%d %.17g %.17g 1e-10 1e-10\n",$2,$8,$9;printf "residual_%d %.17g 0 1e-11 0\n",$2,$15}' \
     "$actual/roots.csv" > "$actual/mpi.dat"
   awk -F, 'NR>1 {print $2,$11,$16,($18>.5),($17>.5)}' "$actual/roots.csv" > "$actual/categorical.dat"
  done
 done
 for method in 0 1 2; do
  one="$scratch/${kind}_np1_method${method}"
  four="$scratch/${kind}_np4_method${method}"
  if [[ -s $one/mpi.dat && -s $four/mpi.dat ]];then
   cmp "$one/categorical.dat" "$four/categorical.dat"
   "$scratch/build/compare_mpi" "$one/mpi.dat" "$four/mpi.dat" > "$scratch/${kind}_method${method}_mpi.txt"
   cat "$scratch/${kind}_method${method}_mpi.txt"
  fi
 done
done
echo 'COMPLETE: root comparison tables and MPI checks; inspect damped_candidate and numerical statuses'
