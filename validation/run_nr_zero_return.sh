#!/usr/bin/env bash
# Before/candidate variants are generated only in scratch; no second production algorithm.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
stage=${1:-post}
case "$stage" in pre|post);;*)echo 'Expected pre or post' >&2;exit 2;;esac
root=$(realpath "${TMPDIR:-/tmp}")
case "$root/" in "$repo/"*)exit 2;;esac
scratch=$(mktemp -d "$root/alps-nr-zero.XXXXXX")
echo "NR zero-return artifacts: $scratch"
mkdir "$scratch/source" "$scratch/build"
cp "$repo"/src/ALPS_{rel_cache,var,io,analyt,com,NHDS,fns_rel,fns}.f90 "$scratch/source/"
cp "$repo/distribution/distribution_analyt.f90" "$repo/validation/test_nr_zero_return.f90" "$repo/validation/nr_zero_trace.f90" "$scratch/source/"
for variant in before candidate;do
 awk -v variant="$variant" -f "$repo/validation/nr_zero_variants.awk" "$scratch/source/ALPS_fns.f90" > "$scratch/$variant.f90"
done
awk -v before="$scratch/before.f90" -v candidate="$scratch/candidate.f90" '
 /^end module/{while((getline l<before)>0)print l;while((getline l<candidate)>0)print l}
 {print}' "$scratch/source/ALPS_fns.f90" > "$scratch/functions.tmp"
mv "$scratch/functions.tmp" "$scratch/source/ALPS_fns.f90"
sed -i 's/private :: integrate_res, landau_integrate, resU/private :: integrate_res, resU/' "$scratch/source/ALPS_fns.f90"
read -r -a flags <<< "${FFLAGS:--O2 -g -fcheck=all -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow}"
read -r -a libs <<< "${ALPS_TEST_LIBS:--llapack -lblas}"
read -r -a mpi_flags <<< "${MPIEXEC_FLAGS:---oversubscribe}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cp "$repo/validation/legacy_glls_reference.f90" "$scratch/source/"
python3 "$repo/validation/instrument_nr_history.py" "$scratch/source"
python3 "$repo/validation/instrument_nr_queries.py" "$scratch/source"
cd "$scratch/build"
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/ALPS_var.f90 ../source/ALPS_io.f90 ../source/distribution_analyt.f90 ../source/ALPS_analyt.f90 ../source/ALPS_com.f90 ../source/ALPS_rel_cache.f90 ../source/ALPS_NHDS.f90 ../source/ALPS_fns_rel.f90 ../source/nr_zero_trace.f90 ../source/ALPS_fns.f90
"${MPIFC:-mpifort}" "${flags[@]}" -c ../source/legacy_glls_reference.f90
"${MPIFC:-mpifort}" "${flags[@]}" -o zero *.o ../source/test_nr_zero_return.f90 "${libs[@]}"
for np in 1 2 4;do
 mkdir -p "$scratch/np$np/distribution"
 awk '/^order=/{ $0="order=2" } /^nperp=/{ $0="nperp=4" } /^npar=/{ $0="npar=8" } /^n_scan=/{ $0="n_scan=0" } /^fit_check=/{ $0="fit_check=F" } /^writeOut=/{ $0="writeOut=F" } {print}' "$repo/tests/test_chebyshev.in" > "$scratch/np$np/test.in"
 (cd "$scratch/np$np";timeout 120 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np "$np" ../build/zero test.in "$stage" >stdout.log 2>stderr.log;tail -4 stdout.log)
 python3 - "$scratch/np$np/query_trace_rank0.csv" <<'PYTRACE'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1])))
def selected(variant,degree=None,selector=None,remove_shifts=False):
 return [[v for k,v in r.items()if k!='variant']for r in rows if r['variant']==variant
  and(degree is None or int(r['degree'])==degree)and(selector is None or int(r['selector'])==selector)
  and(not remove_shifts or r['kind']!='3')]
for method in [0,1]:
 for zero in [0,1]:assert selected(f'method{method}_old{zero}')==selected(f'method{method}_now{zero}')
for degree in [0,1,4]:
 assert selected('original',degree,0)==selected('optimized',degree,0)
 assert selected('original',degree,1,True)==selected('optimized',degree,1)
print('PASS: exact method-0/1 and FD query order; analytic trace removes only interior shifts')
PYTRACE
 for ((rank=1;rank<np;rank++));do
 cmp "$scratch/np$np/integrands_rank0.csv" "$scratch/np$np/integrands_rank$rank.csv"
 cmp "$scratch/np$np/method_policy_rank0.csv" "$scratch/np$np/method_policy_rank$rank.csv"
 cmp "$scratch/np$np/query_trace_rank0.csv" "$scratch/np$np/query_trace_rank$rank.csv"
 cmp "$scratch/np$np/query_counts_rank0.csv" "$scratch/np$np/query_counts_rank$rank.csv"
 done
done
for np in 2 4;do
 cmp "$scratch/np1/contours.csv" "$scratch/np$np/contours.csv"
 cmp "$scratch/np1/method_policy_rank0.csv" "$scratch/np$np/method_policy_rank0.csv"
 cmp "$scratch/np1/query_counts_rank0.csv" "$scratch/np$np/query_counts_rank0.csv"
 cmp "$scratch/np1/query_trace_rank0.csv" "$scratch/np$np/query_trace_rank0.csv"
done
if [[ $stage == post ]];then
 for item in missing:9 inactive:9 unavailable:15 invalid:16 underflow:19 overflow:18 derivative:21 unused_fd:18 unused_old:18;do
 option=${item%:*};status=${item#*:}
 mkdir -p "$scratch/$option/distribution";cp "$scratch/np1/test.in" "$scratch/$option/"
 (cd "$scratch/$option"
 set +e
 timeout 120 "${MPIEXEC:-mpirun}" "${mpi_flags[@]}" -np 2 ../build/zero test.in "$stage" "$option" >stdout.log 2>stderr.log
 code=$?
 set -e
 [[ $code -ne 0 && $code -ne 124 ]]
 rg -q "ERROR: nonrelativistic Chebyshev continuation:" stdout.log
 rg -q "status=$status([^0-9]|$)" stdout.log
 ! rg -q 'SIGFPE|Floating point exception|Fortran runtime error' stderr.log
 )
 echo "PASS: $option explicit continuation status $status"
 done
fi
echo 'PASS: manufactured zero contours, independent stencils, legacy control and MPI equality'
