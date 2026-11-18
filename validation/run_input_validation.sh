#!/usr/bin/env bash
# Task 0 alone: compile the current input reader and run in a fresh directory.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
temp_root=$(realpath "${TMPDIR:-/tmp}")
case "$temp_root/" in "$repo/"*) echo 'TMPDIR must be outside the source tree' >&2; exit 2;; esac
scratch=$(mktemp -d "$temp_root/alps-input-validation.XXXXXX")
printf 'Task 0 artifacts: %s\n' "$scratch"
cd "$scratch"
read -r -a flags <<< "${FFLAGS:--O0 -g -fcheck=all -fbacktrace -ffree-line-length-none}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
"${MPIFC:-mpifort}" "${flags[@]}" -c "$repo/src/ALPS_var.f90" "$repo/src/ALPS_io.f90"
"${MPIFC:-mpifort}" "${flags[@]}" -o input_validation ALPS_var.o ALPS_io.o "$repo/validation/test_input_validation.f90"
for mode in invalid_first invalid_second valid_mixed valid_rel4 valid_rel5 obsolete_logmax supported_rel invalid_kind invalid_bm invalid_grid; do
    mkdir "$mode"
    (
        cd "$mode"
        set +e
        timeout 60 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np 4 ../input_validation validation.in "$mode" >stdout.log 2>stderr.log
        status=$?
        set -e
        if [[ $mode == obsolete_logmax ]]; then
            [[ $status -ne 0 && $status -ne 124 ]]
            rg -q 'ERROR: log_max is no longer supported in &poly_spec.' validation.log
            ! rg -q 'Attempting to read f0' stdout.log
        elif [[ $mode == invalid_* ]]; then
            [[ $status -ne 0 && $status -ne 124 ]]
            case "$mode" in
                invalid_first) rg -q 'requires log_fit=T' validation.log;;
                invalid_second) rg -q 'positive requested polynomial degree' validation.log;;
                invalid_kind) rg -q 'ERROR:' validation.log;;
                invalid_bm) rg -q 'use_bM=T with relat=T is unsupported' validation.log;;
                invalid_grid) rg -q 'ERROR:' validation.log;;
            esac
            ! rg -qi 'out of bounds|Fortran runtime error|Attempting to read f0' stdout.log stderr.log
        else
            [[ $status -eq 0 ]]
            rg -q 'PASS: supported species selections preserved' stdout.log
        fi
        printf 'PASS: %s\n' "$mode"
    )
done
