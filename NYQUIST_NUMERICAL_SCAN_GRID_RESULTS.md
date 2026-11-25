# Rounded momentum grids and numerical-VDF Nyquist scans

**Both production changes are implemented and independently validated.** Optimized/debug synthetic comparisons pass. The one authorized SPAN-i point returns **accepted winding 1, A1 status 0, three agreeing passes**, with the unchanged VDF, original wavevector and prescribed rectangle. Neither full SPAN-i scan was executed.

## Changes and original defect

The old Nyquist validator extrapolated from the first printed increment. For the measured SPAN-i parallel axis this used 0.007732 instead of the full-sequence spacing near 0.007731336, accumulating a 2.66e-4 error. Simply accepting that error while retaining the old increments would leave derivatives, quadrature and fits inconsistent.

Production changes are confined to:

- `src/ALPS_io.f90`: decimal-coordinate precision parsing, affine-axis feasibility and in-memory canonicalization helpers.
- `src/ALPS_nyquist_driver.f90`: use that validator for options 3/4/5; remove the numerical-species scan prohibition; retain all other species/continuation/growth restrictions; add numerical worker/call diagnostics and avoid presenting unused analytic parameters in numerical scan headers.
- `src/ALPS.f90`: invoke numerical validation/canonicalization for all three Nyquist modes, immediately after VDF reading and before derivatives or fits.
- `input.md`, `nyquist.md`: supported modes, canonicalization and physical limitations.

`ALPS_fns.f90`, `ALPS_fns_rel.f90`, Bessel functions, `int_ee`, derivative formulas, continuation models, determinant assembly, harmonic stopping/assignment, A1/A3 and scan geometry are unchanged. Production C2b remains the susceptibility; experimental C2d is not linked.

## Precision model and affine feasibility

A decimal token with exponent e and d digits after the decimal point has quantum `10**(e-d)` and half-rounding radius `r=0.5*10**(e-d)`. Trailing digits, including zeros, count. E/D exponents and signs are handled; plain integer literals have exact integer semantics. Unsupported spellings/exponent ranges are rejected explicitly. This is a nearest-rounded serialization model, not a guarantee for arbitrarily truncated or undocumented-format data.

Each coordinate's interval includes its own decimal radius plus an IEEE allowance `32*epsilon*axis_scale`; coordinate reconstruction is verified with the corresponding arithmetic slack. This bound is independent of filenames and is not a user-adjustable physical tolerance. Every original coordinate is checked, repeated-coordinate intervals are intersected at each axis index, observed monotonicity is checked along every row/column, and endpoint and local-increment constraints are retained. A local increment's decimal allowance is the sum of its endpoint radii. Nonrectangular or nonmonotonic arrays are rejected.

For axis index i, let the common repeated-coordinate interval be `[L_i,U_i]`. A uniform model is `a+i*b`. Every pair i<j implies

```text
b >= (L_j-U_i)/(j-i)
b <= (U_j-L_i)/(j-i).
```

Intersect all these inequalities with positive spacing. A nonnegative perpendicular origin additionally constrains `b <= U_i/i` for i>0 and requires `U_0>=0`. Choose the slope closest to the full-axis least-squares estimate within that feasible interval, then choose an intercept in `[max(L_i-i*b),min(U_i-i*b)]`. A written zero perpendicular origin is preserved whenever feasible. No least-squares residual multiplier is used to expand the decimal bounds.

Already uniform high-precision axes, whose full-sequence residual is within floating-point roundoff, retain their original coordinate bits. Rounded accepted axes are replaced by the inferred affine sequence **in memory only**. There is no interpolation, resampling, smoothing or normalization of f0. The inferred grid need not uniquely recover the unknown pre-serialization grid; it is a demonstrably feasible member of the recorded rounding intervals.

### Actual SPAN-i coordinate model

Read-only tests of the actual 201×402 coordinate arrays accept:

| Axis | Inferred spacing | Maximum coordinate residual | Largest decimal half-radius |
|---|---:|---:|---:|
| Perpendicular | 0.011805960135134907 | 4.98648617286e-7 | 5e-7 |
| Parallel | 0.007731335947368724 | 4.97736805771e-7 | 5e-7 |

Every node uses its own radius, often much smaller than the maximum shown. IEEE allowances are about 1.68e-14 and 1.43e-14 respectively. Maximum local increment residuals are about 9.60135e-7 and 6.64053e-7, checked against summed endpoint radii. Canonical support is `[0,2.3611920270269815]` perpendicular and `[-2.0182867767369101,1.0819789381579481]` parallel. All changes lie inside the input precision model. The measured VDF file and all its values remain untouched.

## Initialization and wavevector lifecycle

| Stage | Coordinate/state handling |
|---|---|
| Ordinary `read_f0` | Loads original values and coordinates. The new validator rereads coordinate tokens only to obtain their precision. |
| Nyquist validation | Rank zero checks values/geometry and replaces only pp axes where canonicalization is needed. Runs once for modes 3,4,5. |
| `derivative_f0` | Uses canonical coordinates and the unchanged original value at every array index. |
| Native Chebyshev/other fits | Initialized afterward, with the same canonical nodes. No continuation setup is skipped. |
| `pass_distribution` / native model transfer | Broadcasts canonical pp, derivative arrays and continuation models. Synthetic startup observers compare rank records. |
| Initial nmax/split/Bessel | Uses canonical pp and the original wavevector through the existing routines. |
| Scan point | Existing scan coordinates are broadcast; every rank updates kperp/kpar. |
| Changed kperp | Existing `update_scan_wavevector` calls unchanged nmax, split and Bessel routines, including halos and allocations. |
| Pure kpar change | Bessel/harmonic state is retained. Resonance momenta and wave-matrix k terms are calculated from current global k in each collective disp. |
| New contour | Existing `run_nyquist` constructs fresh cache/counters. Numerical growth maximization remains forbidden, so its exact-frequency growth cache is inactive; no cross-wavevector determinant reuse occurs. |

The equilibrium, numerical derivatives and fits are independent of wavevector and reused throughout. Workers remain in the existing species/harmonic decomposition, cooperating on one wavevector at a time. No new MPI scheduling or communicator is introduced. Exact kperp zero uses existing tensor limits and zero-argument Bessel values.

## Focused synthetic regressions

The reference is the ordinary ALPS-generated bi-Maxwellian proton table, 32×64 intervals/four-width support, beta_parallel=3, anisotropy=2, with analytic isotropic electrons and the existing C1 contour `[-1,1]×[0.01,0.2]`. Numerical derivatives and production C2b are retained. Seven unique independent option-3 coordinates are compared with line/plane points; all available complete complex determinant sample files are compared byte for byte, not merely by winding.

Optimized and debug results agree:

| Scan | kperp | kpar | Accepted Wn | Physical calls per rank |
|---|---:|---:|---:|---:|
| Line | 0 | .3 | 3 | 628 |
| Line | .001 | .3 | 3 | 628 |
| Line | .3 | .3 | 2 | 626 |
| Parallel | 0 | .3 | 3 | 628 |
| Parallel | 0 | .35 | 3 | 612 |
| Plane | .2 | .25 | 2 | 604 |
| Plane | .2 | .3 | 2 | 608 |
| Plane | .3 | .25 | 4 | 616 |
| Plane | .3 | .3 | 2 | 626 |

These are the computed **coarse production C2b inventories**, not imposed analytic counts. For example, the corresponding analytic response gives two on the near-parallel contour and three at (.3,.25); those known finite-grid/integration discrepancies are not repaired by this task. The acceptance here is unchanged single-point science and correct scan/state reuse.

Every point in both builds has A1 status zero, three agreeing passes, aligned detailed/four-column compact rows, and equal physical calls on all ranks. Scan point samples match the independent option-3 files exactly. All seven numerical option-3 sample files match the pre-change executable exactly. Analytic-only line, parallel and plane compact files and all diagnostic scientific columns (excluding timings) are byte-identical to pre-change output.

The kperp-changing line exercises numerical proton nmax 1,15,43, with proton rank ranges respectively `(0:0,1:1)`, `(0:7,8:15)`, `(0:21,22:43)`. Rank three retains the analytic electron assignment. Pure kpar changes avoid harmonic/Bessel reconstruction. The independently observed canonical coordinates/derivatives agree across ranks.

Grid mathematical tests cover exact/byte-preserved axes, limited-decimal uniform axes, perturbations exceeding precision, genuinely nonuniform axes, nonmonotonic sequences, inconsistent repeated-coordinate intervals, D-exponent precision and a feasible nonnegative-origin corner case. The latter constrains origin and slope jointly rather than rejecting an unfortunate unconstrained fit candidate.

## Rounded initialization and recoverable failure

Both builds additionally initialize the same synthetic f0 from (a) limited-decimal coordinates and (b) the corresponding canonical coordinates printed at full precision. Complete native determinant samples are byte-identical. Startup observations show identical pp and derivative-array checksums on all four ranks. A final high-precision control is also byte-identical to the pre-change option-3 baseline. No distribution value is changed between these serialized representations.

A validation-only callback injects evaluator status 77 after an actual collective determinant at the interior kperp=.001 point of the three-point line. Both builds finish collectively with the required nonzero aggregate status, preserve all three rows, write `NaN NaN` at the failed point, and then converge at kperp=.3. The last point's complete determinant samples match its independent pre-change single point exactly. The callback is confined to copied validation sources; production susceptibility and failure handling are unchanged.

There are **54 synthetic four-rank jobs** across optimized/debug, totaling **56,902 collective determinant calls**. Expected fault-job exits are separate from regression failures. The core matrix and additional phases were executed sequentially; the final shell entry combines those same phases. A final nonnegative-origin feasibility correction affects only the rounded-grid branch: it was tested in both builds with mathematical, rounded/canonical and high-precision physical controls before the measured point. The earlier exact-grid comparisons return before that branch. No claim is made that the final convenience wrapper was rerun from fresh scratch after assembling the phases.

An independent AWK audit of all 80,802 SPAN-i records, using the inferred affine coefficients but not calling the production fitter, finds **zero** violations of the individual decimal intervals plus IEEE slack. Minimum half-decimal radii are 5e-9 perpendicular and 5e-11 parallel; maximum normalized residuals are 0.999999665 and 0.999999718. Thus the actual dataset does not pass merely by comparison with the largest global decimal radius.

## One SPAN-i physical point

Exactly one measured-VDF MPI job was run, after all synthetic gates passed, using the final production executable. Source/input paths remain those in `NYQUIST_SPANI_250320_210001_RESULTS.md`; its earlier failed attempt and outputs remain unchanged. This run used a fresh external directory and the exact prepared `inputs/SPANi-250320-210001-v08a_nyq_point.in`, copied without changes.

- Wavevector: `(kperp,kpar)=(0.002926159066102,0.1463079533051)`.
- Rectangle: real frequency ±0.2194619299576, gamma `[1e-4,0.0025]`, in reference gyrofrequency units; raw ALPS determinant.
- Numerical protons: n=q=m=1, 200×401 grid intervals, AC_method=2, Chebyshev kind 1/order 30, default parallel continuation derivative type one; original normal tabulated derivatives.
- Analytic electrons: n=1,q=-1,m=0.00054466,beta_parallel=0.2241188310194,anisotropy=0.6567986982607,drift momentum=-0.00001215584585925.
- Numerical controls: Bessel_zero=1e-40,Tlim=1e-4,positions_principal=3,n_resonance_interval=100. Growth maximization remains disabled.

Canonicalization completes before ordinary derivatives and native continuation preparation. The resulting production density integral is **1.0000000005559402**, with proton current **-0.022318227589042488** in the existing normalization. Net charge is about 5.56e-10 and net current 3.78e-11. This change from the former first-increment quadrature diagnostic is caused by consistent coordinates, **not VDF renormalization**.

Actual worker assignments are rank one/protons n=0…6, rank two/protons n=7…13, and rank three/analytic electrons n=0…1; rank zero coordinates. Proton nmax is 13, analytic assignment nmax is one. Existing allocation rules give the required one-order halos; no assignments or stopping criteria are changed.

| Diagnostic | Actual measured-point result |
|---|---|
| Launcher exit / A1 status / evaluator status | 0 / 0 / 0 |
| Accepted winding | **1** |
| Global passes | **3**, all phases 6.2831853071795862 |
| Total phase | **2*pi** |
| Logical requests / physical evaluations | 1,030 / 1,030 |
| Actual per-rank physical calls | rank 0=1,030; rank 1=1,030; rank 2=1,030; rank 3=1,030 |
| Exact-frequency cache hits | 0; growth search/cache disabled |
| Refinements / deepest level | 291 / 9 |
| Largest accepted phase increment | 0.457503683619 |
| Minimum final-pass sampled raw abs(D) | **1.1013477830750554e-8** |
| Location of that boundary sample | omega=0.000107141048244064, gamma=0.0001 |
| Retained final-pass samples | 415, including closure |
| Shared startup time | 9.353178 s |
| Native Nyquist wall time | **473.095719 s** |
| Launcher elapsed time | 482.95 s, includes startup/launch/finalization |

No failed determinant frequency or numerical warning was reported; stderr is empty. The legacy diagnostic line “Standard error of the estimate” is not a runtime error and is not a native SVD residual certificate: native rows are explicitly excluded from those legacy fit totals. No independent wave-matrix condition number is available from this production output, and none is inferred from the determinant minimum.

The result is an **accepted winding of the computed numerical boundary trajectory**. Its physical root-count interpretation still requires the usual numerical-VDF representation/continuity evidence. The location of the smallest sampled boundary determinant is not a root estimate. No root frequency, maximum growth, or global stability classification is claimed; no automatic follow-up was performed because the winding is nonzero.

## Reproduction, preservation and limits

`bash validation/run_nyquist_numerical_scan_grid.sh` reproduces the synthetic suite in fresh scratch using the archived pre-change source and the ordinary distribution generator. Set `GRID_SCAN_SCRATCH` to choose a fresh external directory. It performs **no measured SPAN-i determinant calculation or scan**; the measured file is read only by coordinate tests. The separate guarded `run_spani_point.sh STUDY_DIRECTORY` reproduces only the original measured point after the synthetic gates and source hashes pass. That command was executed once for this study. Both scripts retain nonzero outcomes rather than silently treating failures as zero windings.

Compiler/toolchain: GNU Fortran 11.4/Open MPI 4.1.2, LAPACK/BLAS, four ranks, one BLAS/OpenMP thread. Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags are `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`. Initial and finalized source/executable provenance is retained separately. Builds and runs are isolated outside the checkout; original/prepared SPAN-i inputs, the measured array, old scientific outputs and historical validation artifacts pass preservation checks. All unrelated production Fortran sources compare byte for byte with the pre-task snapshot.

New validation code is under `validation/numerical_scan_grid/`; compact native evidence and the independent pre-change source are under `tests/data/nyquist_numerical_scan_grid/`. Archived input/log files avoid the repository's `.in`/`.log` ignore rules. Final source patches, all complex sample trajectories, detailed/compact tables, worker observations, fault outcomes, physical stdout/stderr, summary and samples are retained. No existing fixtures or expected scientific values were changed.

Limits: the precision model assumes ordinary decimal rounding, not arbitrary truncation; plain integers have exact semantics. Extremely small/large unsupported serialization exponents are rejected explicitly. A table compatible with multiple affine grids does not uniquely identify the pre-serialization coordinates. Already uniform high-precision input is preserved to machine precision, rather than forcing a gratuitous reconstruction. This canonicalization is scoped to Nyquist modes; ordinary root/map input behavior is unchanged. Focused scan comparisons cover numerical protons plus analytic electrons; other existing option-3 numerical/analytic combinations use the same generic species loop and unchanged backend, without a new species-convergence matrix in this task.

Both full SPAN-i scan inputs and their old prepared directories remain unexecuted; this task made no full-scan runtime estimate. The historical template comments describe the former backend restriction, now removed, but frequency coverage across their broad k ranges still requires review. **Stop after this one successful point.** No susceptibility, Bessel, int_ee, continuation, Nyquist tolerance, MPI decomposition, measured f0 value, experimental C2d or C3 change is made.
