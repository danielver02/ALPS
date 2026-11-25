# C2f.4a — Oblique coherent susceptibility validation

**Completed on 2026-11-24: qualified / partial physical validation.** The
frozen entry executed cleanly from start to finish in optimized and debug
builds. All 96 A1 contours converged, and tensor/root errors decrease under
momentum refinement. The analytic inventory is two roots in the original
rectangle and **three** in the lower rectangle. The 32×64 both-numerical model
misses the weak third mode; 64×128 and 96×192 recover it. That scientific
comparison remains FAIL in both builds. Arithmetic, Bessel/harmonic, cache,
tensor-assembly, MPI and regression gates pass. Production remains unchanged.

## 1. Source provenance and reused implementation

The equilibrium has beta_parallel,p=beta_parallel,e=3, proton anisotropy 2,
isotropic electrons, unit equal species densities, zero drift, electron mass
`5.44617021487e-4`, vA/c=`1e-4`, and `(kperp,kpar)=(.3,.3)`. The isotropic
control changes only proton anisotropy to one and uses its corresponding
four-width momentum support. All new physical calculations retain positive
kpar; signed real frequency represents propagation direction.

The production baseline is Git revision
`99a087271928846cb8dbc49507d5a3a0f4d68d84` plus the preserved prior working-tree
changes. The launcher records that state and independently reapplies the two
qualified C2f.1 patches to their historical snapshot before comparing current
production byte for byte. Production hashes remain:

| Source | SHA-256 |
|---|---|
| `src/ALPS_fns.f90` | `0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488` |
| `src/ALPS_fns_rel.f90` | `c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6` |

The susceptibility comes from the existing C2f.3/C2f.2/C2d builders and
`helpers_cached.inc`/`model.f90`/`work.f90`; C2f.3a `build_wide.sh` supplies the
unchanged restricted extended Cauchy arithmetic. On this GNU compiler it uses
18-decimal-digit real(10) Cauchy intermediates and weights, returning through
the existing double interface. It is not wholesale real128 susceptibility.
Default production remains C2b. Bessel values, corrected int_ee, derivatives,
fits, Landau continuation, harmonic stopping/assignments, collective disp,
determinant assembly, root solvers, A1/A3 and scan machinery are unchanged.

New files under `validation/c2f4a/` prepare inputs/coordinates, qualify observed
results and reduce diagnostics. There is no new susceptibility formulation.
The additional copied-source observations are all-perpendicular-node stored
Bessel snapshots and compact tensor logging. `BENCH` keeps the original
IO-disabled observer for fair timing. An analytic pilot comparison verifies
compact and original logging produce byte-identical determinants; numerical
one-frequency comparisons test the same property in the frozen entry.

## 2. Independently established analytic inventory

The ordinary analytic solver gives the oscillatory pair approximately

```
+/-0.428821435684 + 0.028412697885 i.
```

This differs from the supplied `.4280853535+.0283366840i` seed, which is the
historical finite-grid production numerical-proton result. That seed is not
used as a prescribed coherent or analytic answer.

Independent analytic A1 calculations reveal an additional zero below gamma=.01:

| Original real bounds | Gamma bounds | Analytic accepted winding |
|---|---|---:|
| [-1,1] | [.01,.2] | 2 |
| [-1,1] | [.001,.2] | 3 |

Both have status zero and three independent global passes. The extra mode is
independently bracketed on the positive imaginary axis and located by the
unchanged ordinary solver near `omega=i*.0094243124`. The difference in the
accepted counts independently confines one additional zero to the intervening
strip. It is described here as an aperiodic analytic mode; no polarization-based
mode identification is claimed. Isotropic analytic counts are measured
separately rather than inferred from the unstable case.

This inventory is the important oblique reference. Requiring winding two on
both rectangles would be incorrect. The final driver remeasures the analytic
reference before judging any numerical representation.

## 3. Bessel and harmonic audit

The primary numerical limits are independently confirmed as **43 for protons**
and **19 for electrons**. The audit reads actual metadata, momentum files and
stored arrays/halos. The unchanged first-below-1e-45 rule is replayed with the
independent C2e real128 series/Fourier reference, extracted verbatim. Stored
values are also compared directly with the qualified production BESSJ API.
No array is renormalized and no stopping rule is modified.

All signed arguments and actual local/halo orders are retained. Reference
checks cover Jn, independently differentiated Jn', n*Jn/x with exact-zero limits,
recurrences and the finite sums. Output separates individual function error,
independently summed omitted J² tail, finite derivative endpoint term,
function-sum error, and double accumulation error. The finite sum of Jn*Jn'
is compared with `-J_N*J_(N+1)`, not incorrectly forced to zero.

For the 64×128 both-numerical calculation, maximum absolute production value and
derivative errors are about `1.78e-16`; nJn/x error is `2.11e-16`. Electron
errors are below `6.2e-17`. Stored-array differences from production are exactly
zero. Proton rank ranges are 0..21 and 22..43; the electron owner has 0..19,
with the original neighboring halos. Actual harmonic returns plus one
complementary moment reconstruct the six species components within `2.5e-16`
absolute at this resolution, and `4.45e-16` over the complete audit. The final evidence includes all grids, isotropic controls
and the diagnostic +8-order configuration.

The sole enlarged-harmonic physical configuration uses both numerical species
at 64×128, with eight added orders through the existing diagnostic builder.
It changes neither the stopping criterion nor the assignment algorithm. Its
matched tensor/wave/determinant comparison distinguishes finite-tail effects
from floating-point regrouping caused by refreshed harmonic ranges.

Both builds' +8 runs observe limits 51/27 and exactly identical wave matrices
and determinants at all 377 query coordinates; maximum species-component
difference is `7.75e-121`. At the original 43/19 limits, the independently
summed maximum omitted J² tails at primary-grid arguments are approximately
`1.34e-94` (proton) and `8.30e-102` (electron). Maximum finite sumJ² rounding
residuals are about `1.0e-15` and `2.2e-16`, respectively. Individual-value
errors and summation roundoff dominate these tiny truncation tails. There is
no evidence that harmonic truncation causes the weak-mode grid discrepancy.

At the `1.5+.02i` primary higher-harmonic probe, actual proton harmonic
six-component norms are about `.07094, .52181` for n=-2,+2, and
`.005746, .01679` for n=-3,+3; multiple harmonics appreciably contribute.
These are complex tensor contribution norms, not positive probability weights;
the separate complementary moment and cancellation remain part of assembly.

## 4. Numerical inputs and matched tensor convergence

Normal ALPS generation/differentiation prepares 32×64, 64×128 and 96×192 grids.
At each resolution, all mixed representations use identical generated numerical
arrays; only analytic versus numerical response selection changes. Proton
unstable support is `pperp_max=4*sqrt(6)`, `|ppar_max|=4*sqrt(3)`; electron
support is `4*sqrt(3*m_e)` in each direction. Isotropic proton support is
`4*sqrt(3)` in each direction. No finite-grid density renormalization or
analytic-derivative substitution is introduced.

The primary 64×128 actual proton and electron derivative dumps are each
byte-identical to the retained C2f.3a `actual_derivatives.tar.gz` arrays. The
comparison is retained in `diagnostics/derivative_preservation.txt`; the oblique
study changes the wavevector, not these supplied numerical derivatives.

The physical query deck is an exact-coordinate union of the two retained
ordered 257-point contours, sixteen representative probes, positive-axis
samples and two static-limit observations: 377 unique frequencies. Every
representation evaluates these coordinates collectively. Projection to each
ordered contour uses exact complex-frequency equality and actual evaluator
success; there is no interpolation or approximate reuse. Independently adapted
A1 samples are never compared pointwise.

Probes include all three analytic roots, both real-frequency signs, points near
the lower boundary, and `+/-1.5+.02i` to exercise higher harmonics. Actual return
streams identify the harmonic/component contributions at the latter frequency.
All six species components, full wave matrices and determinants are retained.
Absolute errors remain when relative component errors are excluded near zero.
Independent numerical-proton and numerical-electron tensors are postassembled
and compared with the native both-numerical result using the existing component
and determinant-perturbation bounds. Observed convergence orders are measured,
not fixed to two.

Optimized matched-determinant errors (median / maximum relative error over
the identical ordered 257 frequencies, dimensionless) are:

| Numerical species | Grid | Original contour | Lower contour |
|---|---|---:|---:|
| Protons | 32×64 | .014572 / 11.964 | .014572 / .85271 |
| Protons | 64×128 | .0036589 / 3.0236 | .0036174 / .21657 |
| Protons | 96×192 | .0016313 / 1.3430 | .0016072 / .094674 |
| Electrons | 32×64 | .0098031 / 8.8588 | .0098031 / .62361 |
| Electrons | 64×128 | .0023422 / 2.2256 | .0023422 / .15664 |
| Electrons | 96×192 | .0010688 / .98882 | .0010688 / .069592 |
| Both | 32×64 | .011855 / 20.464 | .011899 / 1.4512 |
| Both | 64×128 | .0030474 / 5.2269 | .0030596 / .37167 |
| Both | 96×192 | .0014122 / 2.3274 | .0014185 / .16396 |

The original contour passes close to the weak analytic root: minimum matched
analytic |D| is `4.0289e-7`, versus `6.0877e-8` on the lower contour. At the
worst relative-error point on the original contour, both-numerical absolute
errors decrease `8.2447e-6 -> 2.1059e-6 -> 9.3767e-7`. These are distinct from
the maximum absolute errors over the whole curve (`.054681 -> .013971 ->
.0063951`). Relative errors above one are retained; no uniform boundary-error
certificate is claimed. Both-numerical maximum cumulative phase error on the
original contour decreases `.38411 -> .10032 -> .044642` radians. On the
coarse lower contour the total phase difference is `2*pi`, consistent with
the independently accepted inventory mismatch.

Species errors are assessed before determinant assembly. Using
`sqrt(sum_{six independent components}|delta chi|^2)` divided by the same
analytic six-component norm (not an off-diagonal-weighted 3×3 Frobenius norm),
the original-contour results are:

| Species | Grid | Median relative norm | Maximum relative norm | Maximum absolute norm |
|---|---|---:|---:|---:|
| Proton | 32×64 | .00862227 | .0391086 | .0260856 |
| Proton | 64×128 | .00214933 | .00977619 | .00648080 |
| Proton | 96×192 | .000956492 | .00434149 | .00288189 |
| Electron | 32×64 | .0176341 | .529230 | .135598 |
| Electron | 64×128 | .00437689 | .132042 | .0338330 |
| Electron | 96×192 | .00199085 | .0586026 | .0152862 |

Measured median-norm orders are `2.00418, 1.99681` for protons and
`2.01039, 1.94290` for electrons. All six individual complex components,
absolute/relative errors and near-zero exclusions are retained in the raw
matched tensors and component statistics. Maximum errors are not replaced
by the much smaller medians. These finite-grid errors remain relevant to
weak modes even when the A1 count agrees.

## 5. Growing modes and coarse aperiodic-mode discrepancy

At 32×64 the both-numerical oscillatory pair is approximately
`+/-.423592734181+.027569785471i`. The ordinary solver's returned pair is
independently checked for the existing raw residual, drift-free reflection,
actual wave-matrix agreement and positive imaginary-frequency trajectory.

The aperiodic mode is much more sensitive to resolution. Frozen-run coarse
results give:

| Representation | 32×64 aperiodic growth |
|---|---:|
| Analytic / analytic | .0094243124 |
| Numerical protons / analytic electrons | .0022853185 |
| Analytic protons / numerical electrons | .0042257066 |
| Both numerical | No positive-axis bracket in the sampled range |

The coarse both-numerical lower rectangle accepts winding **2**, not the analytic
**3**. This is retained as a failed scientific inventory comparison. The
oscillatory pair remains present; this does not mean the equilibrium is stable.
The missing branch is not assigned gamma=0 or an unconverged midpoint.

An initial unbracketed pilot third-root attempt returned gamma about `1.63e-6`
because raw |D| becomes small near the origin. Its D/omega² remains finite and
it is not accepted as a root. The final workflow attempts an aperiodic ordinary
solve only after an actual finite positive-axis sign bracket, then requires the
return to remain in that bracket. No solver or stopping tolerance changes.
A missing axis bracket alone is not a global zero-exclusion result; the
independent A1 counts establish the prescribed-box inventory.

The frozen ordinary-root results are identical in optimized and debug builds:

| Numerical species | Grid | Positive oscillatory real frequency | Growth | Aperiodic growth |
|---|---|---:|---:|---:|
| Protons | 32×64 | .424732524720 | .027386198756 | .002285318523 |
| Protons | 64×128 | .427802708667 | .028141504028 | .007664710099 |
| Protons | 96×192 | .428368738349 | .028292157440 | .008647135350 |
| Electrons | 32×64 | .427664076380 | .028580294737 | .004225706621 |
| Electrons | 64×128 | .428532785960 | .028454895666 | .008135081018 |
| Electrons | 96×192 | .428693215280 | .028431488617 | .008852881699 |
| Both | 32×64 | .423592734181 | .027569785471 | No resolved bracket |
| Both | 64×128 | .427514757894 | .028184635695 | .006369754970 |
| Both | 96×192 | .428240656582 | .028311123957 | .008075085805 |

The stronger oscillatory roots and the resolved weak branch approach the
analytic values under refinement. Every returned root in the frozen run has
the unchanged ordinary raw residual below `1e-13` (maximum `8.99e-14`). The
independently solved oscillatory partners have maximum reflection error
`7.88e-14`. Full complex signed frequencies and residuals are archived.

Using the complex oscillatory-root error against the independently solved
analytic partner, observed optimized orders are `1.99965, 2.00034` (protons),
`2.00316, 2.00119` (electrons), and `1.99741, 1.99971` (both), for the two
successive refinement ratios. These are measurements, not acceptance
tolerances. Aperiodic roots retain larger relative errors at the same grids.
At 96×192 their growth-rate errors are approximately 8.25% (protons),
6.06% (electrons), and 14.3% (both). Recovering the lower-box inventory at
64×128 and 96×192 is not a claim that this weak growth rate is precisely
resolved. No contour or frequency cutoff is moved to hide that sensitivity.
No coherent calculation below the real axis is used to infer where the coarse
missing mode went.

A read-only cofactor diagnostic at the analytic aperiodic frequency shows that
coarse physical tensor errors, rather than cache disagreement, dominate. At
32×64 the proton/electron yy errors are about `-.01446` / `-.00708` in the
actual ALPS tensor normalization. Their linear cofactor effects account for
about 92% of the net linear determinant perturbation. The full linear prediction
is about `-7.56e-6`, compared with the measured `-7.35e-6`; nonlinear terms
explain the difference. Individual xx/zz/xz effects are larger but largely
cancel. This diagnostic is not a compensating correction.

The evidence supports a finite-grid transverse-response bias in a weak mode.
It does not uniquely separate finite differencing, derivative-PWL representation
and perpendicular quadrature. In particular longitudinal static cancellation
does not constrain yy when ky=0. Preserving that identity is essential but
insufficient to guarantee the physical accuracy of this branch. No derivative,
quadrature or coherent-formulation repair is attempted.

At that same fixed analytic aperiodic coordinate the proton yy error decreases
`-.0144577 -> -.00362931 -> -.00160914`, while the electron yy error decreases
`-.00707966 -> -.00176622 -> -.000783773`. The exact both-numerical determinant
error decreases `7.3547e-6 -> 1.8785e-6 -> 8.3625e-7`; the independently
computed linear cofactor predictions are `-7.5585e-6, -1.8912e-6, -8.3876e-7`.
The reduction supports the finite-grid interpretation without assigning the
coarse missing zero a damped frequency or changing the susceptibility.

## 6. Independent contours and isotropic controls

The unchanged original rectangles are evaluated independently with A1 for all
four representations and all three grid sizes in optimized/debug builds.
Results include accepted signed winding, status, three-pass convergence, total
phase, sampled minimum |D|, adaptive work and per-rank calls. Scientific
count mismatches are retained while the remaining valid cases continue.
The unstable A1 results below are identical in both builds:

| Numerical species | Grid | Original winding / calls | Lower winding / calls |
|---|---|---|---|
| None | All three | 2 / 630 | 3 / 624 |
| Protons | 32×64 | 2 / 606 | 3 / 620 |
| Protons | 64×128 | 2 / 618 | 3 / 622 |
| Protons | 96×192 | 2 / 622 | 3 / 622 |
| Electrons | 32×64 | 2 / 614 | 3 / 624 |
| Electrons | 64×128 | 2 / 618 | 3 / 622 |
| Electrons | 96×192 | 2 / 628 | 3 / 624 |
| Both | 32×64 | 2 / 608 | **2 / 638** |
| Both | 64×128 | 2 / 614 | 3 / 622 |
| Both | 96×192 | 2 / 618 | 3 / 622 |

Every listed result has A1 status zero and three converged global passes.
Total phase is the corresponding accepted integer times `2*pi`. The bold
entry is an accepted finite-grid count but fails the analytic-inventory
comparison. All isotropic configurations at all three resolutions in both builds
accept zero on both rectangles. Final evidence retains the individual minima,
phase, adaptive refinements/depth and rank counters, rather than substituting
the fixed-query phase diagnostics for A1 results.

Across all 96 contours, status is zero and passes=3; refinements range from
53 to 95, maximum depth is 9, and maximum absolute phase residual from the
accepted integer times `2*pi` is `1.78e-15`. Individual sampled minima |D|
range from `9.14e-9` to `2.54e-5`; their coordinates and all work counts are
in `results/contours.dat`. Original/lower adaptive sample files are
byte-identical across builds (48 file pairs). No stability claim is made
outside these domains. A fixed 257-point phase sum remains a diagnostic
trajectory, not a substitute for independent A1 acceptance.

## 7. Precision, cache, tensor identities and analyticity

At 64×128, the predetermined 49-point cache deck contains a uniform contour
subsample and all sixteen probes, including the analytic roots. The unchanged
criterion is `|Dc-Dd| <= 2e-12*max(1e-6,|Dc|)`. No points are removed after
observing results. The frozen optimized **and debug** comparisons all pass:

| Numerical species | Maximum relative D difference, either build |
|---|---:|
| Protons | 4.46747e-16 |
| Electrons | 9.42511e-15 |
| Both | 6.21225e-15 |

All 49 coordinates are checked at the original criterion. A full 257-point
direct comparison is not repeated because these representative gates pass;
the physical analytic/coherent trajectories each retain all 257 points.
The old C2f.3a arithmetic alternatives study is not repeated.

Existing C2f.3a primitives, moments, static/full-tensor identities and exact
cache-invalidation fixtures are reused. Representative independent real/imaginary
frequency perturbations test numerical CR consistency. At the analytic
oscillatory-root coordinate the both-numerical determinant CR residual
falls approximately `2.56e-5 -> 2.57e-7 -> 8.40e-10` with tenfold step
reductions. Weak-component absolute diagnostics remain available.

All six physical CR gates pass. Their four centres are the analytic positive
oscillatory root, `.4+.001i`, `1.5+.02i`, and the proton usable-support-edge
probe `2.04598501644+.01i`. Reused optimized/debug mathematical tests also
pass: independent signed affine Cauchy error `4.71e-16`, weight-reuse and
frequency/species/kpar/grid invalidation error `2.81e-15`, alternating
mass-scaled-species error `3.56e-15`, and exact static longitudinal residual
`9.89e-17`. No inherited tolerances are changed. Numerical tensor/CR records
retain both species separately and the assembled response.

Fixed derivative data/domains, frequency-independent Bessel coefficients,
Cauchy kernels and polynomial frequency dependence retain the established
upper-half-plane holomorphy. Numerical CR evidence and finite-grid convergence
are separate from that mathematical property. No Tlim or support-entry switch,
new continuation or damped coherent response is introduced.

## 8. MPI and performance

All actual determinant requests use the original collective disp and harmonic
assignments. No ranks are assigned frequencies, contour pieces or wavevectors.
Physical calls and worker/rank workload counters are compared explicitly.
Primary benchmarks use identical 49 frequencies for C2b production, cached
coherent and direct coherent at 64×128. Original IO-disabled BENCH binaries
exclude the compact observer's tensor writes; constructor/initialization timing
is separate. No optimization is performed.

Cauchy-weight reuse is distinct from the Nyquist exact-frequency determinant
cache. The unchanged single-point driver enables that determinant cache with
`nyq_find_max_growth`; this study keeps growth search disabled, so the A1
summaries correctly report cache_enabled=F, zero determinant-cache hits and
logical requests equal to physical calls. Reused Cauchy weights operate within
each actual determinant evaluation and do not remove or approximate MPI
frequency requests.

The completed study has **237 MPI jobs, 69,823 collective disp calls and 948
rank records, with zero per-job rank-call mismatches**. These totals include
study probes, roots, audits and benchmarks; the separately run compatibility
suite is not included. Numerical-proton groups are 0..21 and 22..43 on two
workers; electron-only groups are 0..9 and 10..19 on two workers. With both
species numerical, two workers retain the proton groups and the remaining
worker owns electron 0..19. Halos retain the existing neighbouring entries.

Measured optimized wall seconds per collective determinant, using the maximum
rank's accumulated MPI_Wtime interval divided by 49 requests:

| Numerical species | Production C2b | Cached coherent | Direct coherent |
|---|---:|---:|---:|
| Protons | .0801202 | .117809 | 1.31735 |
| Electrons | .0299972 | .0531951 | .595960 |
| Both | .0794083 | .116744 | 1.30585 |

Production/cached values average two repetitions; direct has one repetition
per representation, so these are representative measurements rather than a
general scaling claim. Weight reuse gives approximately **11.2× acceleration
over direct coherent evaluation**, compared with the earlier nearly parallel
6.6× reference. It remains **1.47×, 1.77× and 1.47× slower than C2b** for these
three configurations. No claim of universal accuracy or speed superiority
over C2b follows from this experiment.

Instrumented totals over one 49-frequency cached run are:

| Numerical species | Kernel-cell operations | Derivative-node reads | Weight sets | Sum of rank weight-cache bytes |
|---|---:|---:|---:|---:|
| Protons | 1,074,276 | 406,946,862 | 4,263 | 709,864 |
| Electrons | 481,572 | 181,126,638 | 1,911 | 319,336 |
| Both | 1,555,848 | 588,073,500 | 6,174 | 869,532 |

Direct kernel-cell totals are 201,871,278 / 89,850,222 / 291,721,500; derivative
reads are unchanged and no reusable weight sets are stored. These are logical
instrumented operation counts, not hardware-instruction counts. Cache bytes
measure weight storage, not whole-process RSS. Per-rank counters retain the
unequal harmonic workload. In the both-numerical cached run, signed-harmonic
weight-set counts per worker are 2,107 / 2,156 / 1,911; all four ranks execute
49 collective disp calls.

Bessel arrays and harmonic limits are constructed once per numerical worker
per initialization, never per frequency. Across the primary cached benchmarks,
maximum worker array-construction times are about .423/.0912/.424 ms and
harmonic-limit times .610/.124/.715 ms for proton/electron/both cases. Maximum
worker weight-construction totals over 49 requests are .312/.143/.314 s.
These isolate the existing counters; full VDF reading/fitting startup is not
newly profiled. Larger harmonic limits are a physical workload difference,
not an opportunity to change decomposition. No optimization was added.

## 9. Regression and reproducibility

Run `bash validation/run_nyquist_c2f4a.sh` with an optional fresh external
`C2F4A_SCRATCH`. The frozen final run is
`/tmp/alps-c2f4a-acceptance-final-20261124`. It builds current qualified
production sources, uses ordinary VDF generation, and runs four-rank jobs
sequentially. No pre-existing /tmp artifacts are required by that entry.

GNU Fortran 11.4.0 / Open MPI 4.1.2 with existing system BLAS/LAPACK/OpenBLAS
are used; OMP_NUM_THREADS and OPENBLAS_NUM_THREADS are both one.
Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`;
debug uses `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
A limited nearly parallel numerical-proton control, ordinary production oblique
root control, and existing A1/A2/legacy compatibility gate protect construction.
The old electron support/resolution and low-frequency matrices and twelve
compatibility suites are not wholesale repeated.

The final entry was frozen before execution and completed with exit=0 and
the final COMPLETE line. Executed entry/validation files still match their
frozen copies. Expected PROBE exit 1/status 6 is
not an A1 failure. Execution, communication, input and preservation failures
abort; scientific discrepancies remain explicit FAIL and yield PARTIAL while
measurements continue. Both scientific FAIL records are the same coarse-grid
inventory discrepancy; missing coarse aperiodic axis brackets are explicit NA,
not gamma=0. No numerical evaluator or MPI failure occurred.

Both builds pass the reused C2f.3a mathematical tests, all six primary cache
and CR gates, native species assembly, harmonic audits and +8 control, nearly
parallel control and existing `run_nyquist.sh` A1/A2/legacy compatibility gate.
The nearly parallel 257-point numerical-proton reference comparison has maximum
relative D difference `2.523e-14`. Ordinary C2b oblique control roots remain
`+/-.42808535353519+.02833668396181i`, consistent with the preserved C2f.1
reference at the existing `1e-10` control bound. The coherent method is not
activated in ordinary production. Unrelated inherited harness issues are not
claimed to be repaired, and the twelve older compatibility suites are not
claimed to have been rerun.

Cross-build postchecks find byte-identical tensors for all twelve matched
configurations, twelve oscillatory root files, nine aperiodic root files, and
all 48 original/lower adaptive-sample file pairs. Preservation checks pass for
121 original source/report files, 998 recent historical files and 402 frozen
source files. No production hashes or expected historical results were
updated. Historical C2d–C2f.3a reports and fixtures remain unchanged.

Compact evidence, inputs, actual derivative arrays, generated source variants,
mathematical logs, per-case samples/tensors/counters, provenance and successful
preservation checks are retained under `tests/data/nyquist_c2f4a/`. See its
README and SHA256SUMS. Absolute paths in provenance identify the executed
snapshot; reproduction does not depend on that old scratch directory.

## 10. Limitations and disposition

No production-default adoption, non-Maxwellian calculation, new oblique family,
wavevector scan, maximum-growth search or C3 damped continuation is performed.
**Disposition: partial physical-validation success, with all execution and
numerical-consistency gates passed.** The oblique six-component tensors and
oscillatory roots converge toward the analytic reference, the weak analytic
mode is recovered at both finer grids, and no unexplained cache, harmonic,
tensor-assembly or MPI inconsistency is found. The coarse both-numerical
lower-box inventory failure remains, and even 96×192 underestimates this weak
growth by about 14.3%. Therefore the tested coarse grid cannot be recommended
for identifying this weak instability, and a winding match at a finer grid
does not establish precise growth-rate accuracy.

Finite support is fixed, as requested. This study does not independently
separate finite-difference derivative error, derivative-PWL representation and
perpendicular quadrature error. It does not certify roots outside the tested
rectangles or arbitrarily near zero, establish a damped causal continuation,
or generalize the reference's minimum usable grid to other equilibria.

For C2f.4b, retain the validated arithmetic, source construction and independent
analytic/tensor controls, and explicitly qualify momentum resolution when weak
modes determine the inventory. No kernel correction is justified by these
data alone. This is a recommendation for further validation, not production
adoption. The present task stops here.
