#!/usr/bin/env bash
# Reuse the isolated Task 1 NR workflow. All hooks and selectors are test-only.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
root=$(realpath "${TMPDIR:-/tmp}")
case "$root/" in "$repo/"*)exit 2;;esac
scratch=$(mktemp -d "$root/alps-nr-svd-physics.XXXXXX")
echo "NR SVD physics artifacts: $scratch"
shape=${1:-maxwell}
case "$shape" in maxwell|kappa);; *)echo 'Expected maxwell or kappa' >&2;exit 2;;esac
for backend in legacy shared svd_compat svd_native;do
 export ALPS_NR_SVD_BACKEND=$backend
 awk -v repo="$repo" -v shape="$shape" '
 /^repo=/{print "repo="repo;next}
 /^for kind in rel4 rel5 nr_fd nr_analytic; do/{print "for kind in nr_fd nr_analytic; do";next}
 /^cd "\$scratch\/build"/ {
 if(shape=="kappa") print "sed -i '\''s/distributions=1/distributions=2/'\'' \"$scratch/seed/test_bimax_dist.in\""
 print "for pair in driver:validation/test_numerical_baselines.f90 fns:src/ALPS_fns.f90;do"
 print "  awk -v mode=\"${pair%%:*}\" -v shape=\""shape"\" -f \"$repo/validation/nr_svd_instrument.awk\" \"$scratch/source/${pair#*:}\" > \"$scratch/instrumented.tmp\""
 print "  mv \"$scratch/instrumented.tmp\" \"$scratch/source/${pair#*:}\""
 print "done"
 }
 /^objects=/{print "\"$fc\" \"${flags[@]}\" -c \"$repo/validation/nr_svd_physics.f90\"";sub(/\)$/," nr_svd_physics.o)")}
 /        cat root_status.dat/{next}
 /build\/compare_baselines.*\$repo\/tests\/data/ {if(shape=="kappa"){print "            echo '\''Kappa: comparison only; no historical reference asserted'\''";next}}
 {print}' "$repo/validation/run_numerical_baselines.sh" > "$scratch/$backend.sh"
 ALPS_TEST_LIBS="${ALPS_TEST_LIBS:--llapack -lblas}" bash "$scratch/$backend.sh" --check > "$scratch/$backend.log" 2>&1
 cat "$scratch/$backend.log"
 baseline_root=$(awk '/Task 1 artifacts:/{print $4;exit}' "$scratch/$backend.log")
 if [[ $backend == legacy ]];then legacy_root=$baseline_root;fi
 for selector in nr_fd nr_analytic;do
  awk -v shape="$shape" -v backend="$backend" -v selector="$selector" -f "$repo/validation/compare_nr_svd_physics.awk" \
   "$legacy_root/$selector/observations.dat" "$baseline_root/$selector/observations.dat" > "$scratch/${backend}_${selector}_comparison.csv"
 done
done
