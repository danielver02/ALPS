#!/usr/bin/env bash
# Existing NR baseline workflow, scratch-only backend/zero-return observations.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
root=$(realpath "${TMPDIR:-/tmp}")
case "$root/" in "$repo/"*)exit 2;;esac
scratch=$(mktemp -d "$root/alps-nr-shared-physics.XXXXXX")
echo "NR paired physics artifacts: $scratch"
for backend in legacy shared;do
 export ALPS_NR_SHARED_BACKEND=$backend ALPS_NR_SHARED_REPO=$repo
 awk -v repo="$repo" '
 /^repo=/{print "repo="repo;next}
 /^for kind in rel4 rel5 nr_fd nr_analytic; do/{print "for kind in nr_fd nr_analytic; do";next}
 /^cd "\$scratch\/build"/ {
 print "awk '\''{print} /^  call pass_distribution$/ {print \"  block\";print \"    use alps_analyt, only: CHEB_SUCCESS\";print \"    use legacy_glls_reference, only: initialize_historical_nr,select_nr_test_evaluation\";print \"    integer :: shared_status\";print \"    character(16) :: backend\";print \"    call get_environment_variable(\\\"ALPS_NR_SHARED_BACKEND\\\",backend)\";print \"    if(trim(backend)==\\\"shared\\\") then\";print \"      call initialize_historical_nr(shared_status)\";print \"      if(shared_status/=CHEB_SUCCESS) error stop \\\"NR shared initialization\\\"\";print \"      call select_nr_test_evaluation(.true.)\";print \"    endif\";print \"  end block\"}'\'' \"$scratch/source/validation/test_numerical_baselines.f90\" > \"$scratch/source/validation/driver.tmp\""
 print "mv \"$scratch/source/validation/driver.tmp\" \"$scratch/source/validation/test_numerical_baselines.f90\""
 print "awk '\''/landau_integrate =cmplx/ {print \"     write(*,*) \\\"NR_ZERO_RETURN\\\",sproc,iperp,nn,mode,om,landau_integrate\"} {print}'\'' \"$scratch/source/src/ALPS_fns.f90\" > \"$scratch/source/src/observed.tmp\""
 print "mv \"$scratch/source/src/observed.tmp\" \"$scratch/source/src/ALPS_fns.f90\""
 }
 {print}' "$repo/validation/run_numerical_baselines.sh" > "$scratch/$backend.sh"
 ALPS_TEST_LIBS="${ALPS_TEST_LIBS:--llapack -lblas}" bash "$scratch/$backend.sh" --check > "$scratch/$backend.log" 2>&1
 cat "$scratch/$backend.log"
done
