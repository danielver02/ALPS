#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
root=$(realpath "${TMPDIR:-/tmp}")
case "$root/" in "$repo/"*)exit 2;;esac
scratch=$(mktemp -d "$root/alps-nr-svd.XXXXXX")
echo "NR SVD artifacts: $scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/ALPS_{var,io,analyt,com}.f90 "$scratch/source/"
cp "$repo/distribution/distribution_analyt.f90" "$repo/validation/test_nr_svd.f90" "$scratch/source/"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cp "$repo/validation/legacy_glls_reference.f90" "$scratch/source/"
cd "$scratch/build"
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/distribution_analyt.f90 ../source/ALPS_analyt.f90 ../source/ALPS_com.f90
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/legacy_glls_reference.f90
"${MPIFC:-mpifort}" "${flags[@]}" -o svd *.o ../source/test_nr_svd.f90 "${libs[@]}"
for np in 1 2 4;do
 mkdir -p "$scratch/np$np/distribution"
 awk '/^order=/{ $0="order=4" } /^nperp=/{ $0="nperp=8" } /^npar=/{ $0="npar=6" } /^n_scan=/{ $0="n_scan=0" } /^fit_check=/{ $0="fit_check=F" } /^writeOut=/{ $0="writeOut=F" } {print}' "$repo/tests/test_chebyshev.in" > "$scratch/np$np/test.in"
 (cd "$scratch/np$np";timeout 120 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np "$np" ../build/svd test.in >stdout.log 2>stderr.log;tail -5 stdout.log)
done
