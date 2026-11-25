# SPAN-i 250320-210001 — single-point attempt and prepared scans

**The single-point calculation was attempted once and failed during input validation. No winding was calculated. Neither scan was executed.** Current production also explicitly forbids numerical species in Nyquist scans, so the scan files are blocked templates, not ready scientific runs. Production sources and both supplied scientific files are unchanged.

## Original configuration

The original `inputs/SPANi-250320-210001-v08a_map_k05.in` specifies `(kperp,kpar)=(0.002926159066102,0.1463079533051)` in `k*d_p` units and `vA/c=1e-4`.

| Species | Representation | Retained parameters |
|---|---|---|
| Proton | Actual supplied numerical VDF; nonrelativistic | n=q=m=1, `AC_method=2,ff=2,log_fit=T`, Chebyshev kind 1/order 30 |
| Electron | Analytic bi-Maxwellian; nonrelativistic | n=1,q=-1,m=0.00054466, beta_parallel=0.2241188310194, anisotropy=0.6567986982607, drift momentum=-0.00001215584585925 |

The electrons are anisotropic/drifting; no species was replaced. All new species/continuation groups match the original exactly. Numerical controls remain `Bessel_zero=1e-40,positions_principal=3,n_resonance_interval=100,Tlim=1e-4`. Normal centered VDF derivatives and the unspecified method-2 `derivative_type` default (one for parallel continuation differentiation) are preserved.

The original map requests 12 roots, ±0.2194619299576 real frequency, gamma [-0.004,0.0025], `ni=nr=128`. Root settings are `secant_method=2,numiter=50,D_threshold=1e-30,D_prec=1e-5,D_gap=1e-5,D_tol=1e-7`. These root controls remain present but unused: Nyquist disables roots/maps.

Read-only inspection of `distribution/SPANi-250320-210001-v08a.1.array` found **80,802 finite, strictly positive samples**, 201×402 nodes (`nperp=200,npar=401` intervals). Perpendicular support is [0,2.361192], parallel support [-2.018287,1.081979]; repeated axes are rectangular. VDF values span `7.137809e-14`–`2.535044`. Values are read directly in ALPS's reference-proton momentum/density normalization; the file supplies no independent SI/calibration metadata. The production-style full-node rectangle integral is 1.00008586713131, with mean parallel momentum -0.0223182266264271. No renormalization was applied.

## Generated configurations and contour

All files below are under `inputs/`, prefixed `SPANi-250320-210001-v08a_`:

| Suffix | Requested calculation | Execution |
|---|---|---|
| `nyq_point.in` | Option 3, original wavevector | One failed startup attempt |
| `nyq_parallel31.in` | Option 4, kperp exactly 0, kpar 0.01–10 | **Not executed** |
| `nyq_plane31x31.in` | Option 5, both components 0.01–10 | **Not executed** |

Every file retains `arrayName='SPANi-250320-210001-v08a'`, independently of its output basename. Separate run directories under `validation/spani_250320_210001/results/` contain copied inputs, separate `solution/` directories, and links to the original measured array. No input/VDF alias or production I/O change is needed.

The explicit fixed rectangle is **omega/Omega_p ∈ [-0.2194619299576,0.2194619299576], gamma/Omega_p ∈ [1e-4,0.0025]**. Real/upper bounds come from the existing map; the strictly positive lower boundary is the supported default. Growth maximization is disabled and samples requested. This narrow initial rectangle is not a global root inventory; its adequacy across three decades in k requires review before future scans. No automatic frequency rescaling is implemented.

A1 defaults are unchanged: 8 seed intervals/edge; depth 24; 200,000 evaluations; up to 6 passes; phase/curvature/log-magnitude/local-ratio limits 0.7/0.15/1.5/0.05; phase and integer tolerances 1e-8. Acceptance requires three agreeing global phases, including an independent seeded pass.

Static review of `scan_read` and `scan_wavevector` confirms `ns=30,nres=1,swlog=T,swf=10` gives `10**(-2+j/10)`, j=0…30. The parallel file uses type 4; the plane uses type 3 outer/type 4 inner, giving 961 points with kpar fastest and all corners included. Static coordinate tables are retained; **they are not scan execution outputs**. Tensor code has explicit exact-kperp-zero limits, so no small nonzero substitute was made. However, `src/ALPS_nyquist_driver.f90:77–78` rejects numerical species for options 4/5. That backend restriction remains intact.

Supported scan output conventions would be `.nyquist.scan` diagnostics, `.nyquist.summary`, and four-column `.nyquist.dat`: `kperp kpar Wn gamma_max`. Wn denotes the original rectangle winding; failed counts remain NaN, accepted zero remains zero, and growth is NaN here because maximization is disabled. Plane diagnostic indices are one-based. No such scientific files were generated.

## Actual single-point result

The unchanged production code was built using `validation/build_nyquist.sh` in external scratch, with `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`, LAPACK/BLAS. Exactly one `mpirun --oversubscribe -np 4` job ran at the original wavevector, with single-threaded BLAS/OpenMP.

It successfully read the supplied proton table, then reported:

```text
ERROR: Numerical Nyquist requires a uniform rectangular momentum grid.
MPI_ABORT ... rank 0 ... errorcode 31
```

The validator (`ALPS_nyquist_driver.f90:107–131`) extrapolates each axis from its first stored increment. Rounded coordinates fail its existing tolerance:

| Grid diagnostic | Perpendicular | Parallel |
|---|---:|---:|
| First stored increment | 0.01180596 | 0.007732 |
| Maximum affine-grid residual | 5.20e-7 | 2.66e-4 |
| Permitted residual | 2.361192e-8 | 2.018287e-8 |

80,684 nodes violate a comparison; the first is `(0,2)`. Decimal rounding is consistent with the discrepancy, but the exact original higher-precision grid is unavailable and was not reconstructed.

| Requested diagnostic | Outcome |
|---|---|
| Measured VDF load | Completed; finite/positive independent audit passes |
| Numerical derivatives / continuation initialization | **Not reached** |
| Harmonic cutoff / worker assignments / Bessel arrays | **Not reached** |
| Accepted winding / A1 status / passes | **Unavailable: A1 never entered**, not accepted zero |
| Phase / minimum abs(D) / refinements / depth | No contour samples or refinement performed |
| Physical determinant calls, ranks 0/1/2/3 | 0/0/0/0 inferred from startup control flow; no instrumented counters |
| MPI exit | 31; workers terminated by MPI_ABORT |
| Launcher elapsed time | 0.53 seconds, including launch/abort, not Nyquist evaluation time |
| Failed frequencies / conditioning | No frequency evaluated; no conditioning diagnostic |

Logs are retained in `validation/spani_250320_210001/results/point/`: `stdout.txt`, `stderr.txt`, input-specific `.log`, `resources.txt`, and `mpi_exit_code.txt`. No `.nyquist.summary` was produced. Build provenance, VDF audit, static checks and preservation hashes are in the parent results directory. This failure establishes neither stability nor instability. No tolerances, contour limits or distribution data were adjusted after it.

## Launch commands for later review

**These scan commands were not executed and will fail with the current backend.** Resolve numerical-scan support, the grid contract and contour coverage in a separate authorized task first. From the repository root:

```bash
REPO=$PWD
BUILD=$(mktemp -d /tmp/alps-spani-reviewed-build.XXXXXXXX)
bash validation/build_nyquist.sh "$BUILD"

# Parallel scan, later only:
(cd "$REPO/validation/spani_250320_210001/results/parallel31" && \
 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 mpirun --oversubscribe -np 4 "$BUILD/ALPS" \
 SPANi-250320-210001-v08a_nyq_parallel31.in > stdout.txt 2> stderr.txt)

# Plane scan, later only:
(cd "$REPO/validation/spani_250320_210001/results/plane31x31" && \
 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 mpirun --oversubscribe -np 4 "$BUILD/ALPS" \
 SPANi-250320-210001-v08a_nyq_plane31x31.in > stdout.txt 2> stderr.txt)
```

Use fresh directories after any future attempt. The single-point reproduction helper `validation/spani_250320_210001/run_single.sh` builds/runs only the point in fresh scratch; it was syntax-checked, **not executed for a second physical attempt**. Details and the direct executed command are retained in [validation notes](validation/spani_250320_210001/README.md).

Both prepared scan directories have no run logs or scientific outputs. Protected source/original-input/VDF hashes pass unchanged. No additional wavevectors, ordinary roots, maps, growth searches, coherent integration or C3 continuation were run. **Stopped for review of the failure and unsupported scans.**

The repository ignores `.in`/`.log` files; `results/prepared_inputs_and_logs.tar.gz` additionally retains all three generated inputs and actual run/build logs so they are not lost from a normal Git-based handoff. The measured array is not bundled or modified. `validation/spani_250320_210001/SHA256SUMS` verifies the retained evidence from the repository root.
