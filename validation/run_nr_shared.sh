#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
root=$(realpath "${TMPDIR:-/tmp}")
case "$root/" in "$repo/"*) exit 2;;esac
scratch=$(mktemp -d "$root/alps-nr-shared.XXXXXX")
echo "NR shared artifacts: $scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/ALPS_{rel_cache,var,io,analyt,com,NHDS,fns_rel,fns}.f90 "$scratch/source/"
cp "$repo/distribution/distribution_analyt.f90" "$repo/validation/test_nr_shared.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Expose the unchanged private routine only in the scratch test build.
sed -i 's/private :: integrate_res, landau_integrate, resU/private :: integrate_res, resU/' "$scratch/source/ALPS_fns.f90"
awk '/landau_integrate =cmplx/ {print "     write(*,*) \"NR_ZERO_RETURN\",sproc,iperp,nn,mode,om,landau_integrate"} {print}' "$scratch/source/ALPS_fns.f90" > "$scratch/source/observed"
mv "$scratch/source/observed" "$scratch/source/ALPS_fns.f90"
cp "$repo/validation/legacy_glls_reference.f90" "$scratch/source/"
cp "$repo/validation/legacy_glls_reference.f90" "$scratch/source/"
python3 "$repo/validation/instrument_nr_history.py" "$scratch/source"
cd "$scratch/build"
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/distribution_analyt.f90 ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 ../source/ALPS_rel_cache.f90 ../source/ALPS_NHDS.f90 ../source/ALPS_fns_rel.f90 ../source/ALPS_fns.f90
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/legacy_glls_reference.f90
"${MPIFC:-mpifort}" "${flags[@]}" -o shared *.o ../source/test_nr_shared.f90 "${libs[@]}"
for np in 1 2 4;do
 mkdir -p "$scratch/np$np/distribution"
 awk '/^order=/{ $0="order=1" } /^n_scan=/{ $0="n_scan=0" } /^fit_check=/{ $0="fit_check=F" } /^writeOut=/{ $0="writeOut=F" } {print}' "$repo/tests/test_chebyshev.in" > "$scratch/np$np/test.in"
 (cd "$scratch/np$np";timeout 120 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np "$np" ../build/shared test.in >stdout.log 2>stderr.log;tail -5 stdout.log)
done

for np in 2 4;do
 cmp "$scratch/np1/comparisons.csv" "$scratch/np$np/comparisons.csv"
 cmp "$scratch/np1/direct_glls.csv" "$scratch/np$np/direct_glls.csv"
 cmp "$scratch/np1/stencil_policy.csv" "$scratch/np$np/stencil_policy.csv"
done
mkdir -p "$scratch/fatal/distribution"
cp "$scratch/np1/test.in" "$scratch/fatal/"
(cd "$scratch/fatal"
 set +e
 timeout 120 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np 2 ../build/shared test.in fatal >stdout.log 2>stderr.log
 code=$?
 set -e
 [[ $code -ne 0 && $code -ne 124 ]]
 rg -q 'ERROR: nonrelativistic Chebyshev continuation:' stdout.log
 rg -q 'iperp=1 status=9' stdout.log
 ! rg -q 'SIGFPE|Floating point exception' stderr.log
)
echo 'PASS: shared NR values/derivatives, exact cross-rank tables, low-degree transfer and explicit unavailable-model failure'
