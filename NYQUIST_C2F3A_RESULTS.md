# C2f.3a — Cancellation-aware Cauchy arithmetic

This is a focused arithmetic experiment in validation copies. Production ALPS,
its C2b default, the physical distributions, derivatives, susceptibility
identities, MPI assignments, and scientific acceptance tolerances are unchanged.
The original C2f.3 report and archived failed comparisons are preserved.

**Completed 2026-11-24: success for the focused arithmetic study.** Restricted
extended Cauchy arithmetic improves independently measured accuracy and passes
all four original strict comparisons without changing their threshold. The
complete final frozen entry finished with exit status zero. Production adoption
is neither performed nor implied.

## 1. Baseline and reproduction

The source baseline is Git revision `99a087271928846cb8dbc49507d5a3a0f4d68d84`
with the pre-existing working-tree changes recorded by the launcher. Qualified
C2f.1 patches are independently reapplied to the historical snapshot and compared
byte for byte with the current production files before building candidates.

| Production file | SHA-256 |
|---|---|
| `src/ALPS_fns.f90` | `0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488` |
| `src/ALPS_fns_rel.f90` | `c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6` |

The experiment reuses C2f.3's source builder, retained 64×128 inputs, normal VDF
generator/derivatives, 257 ordered frequencies, actual tensor observer, and
C2d cached/uncached helpers. It does not reconstruct the physical study matrix.
Only numerical electrons/analytic protons (E), both numerical species (B), and
one numerical-proton control are generated. Original momentum support and all
A1 parameters are retained.

The unchanged `validation/c2d/compare.f90` gate is

```
|D_cached-D_direct| <= 2e-12 * max(1e-6, |D_cached|).
```

All matched coordinates have sufficient determinant magnitude for the relative
part to apply. Each of the four original comparisons has **15 failing rows**:

| Representation | Builds | Maximum relative difference | Maximum absolute difference |
|---|---|---:|---:|
| E | optimized/debug | 9.40378384050979e-12 | 9.19561877216816e-12 |
| B | optimized/debug | 9.17168989873391e-12 | 9.21052450995151e-12 |

The failing coordinate sets are identical. On the bottom edge, gamma=.01, the
real frequencies are `-.625,-.59375,-.5625,-.53125,-.4375,.34375,.40625,.4375,
.53125,.5625,.59375,.625,.75,.8125`; the remaining point is `-.5625+.2i`.
Every failed coordinate, all six component differences for both species, and
cofactor-predicted determinant changes are retained in `original_*.frequencies`,
`.components`, and `.determinants`.

At `-.53125+.01i`, the original E determinants are
`-.009570272584982713-.003346477743117304i` and
`-.009570272584930033-.003346477743196769i`.
The absolute difference is `9.53401921535831e-14`, amplified to the maximum
relative discrepancy by the small determinant. Linear tensor cofactors predict
the measured determinant changes with a maximum absolute residual about
`4.7e-16` over the failed coordinates. The observed error is dominated by electron
xx, yy and xy; zz changes are much smaller. This identifies tensor roundoff
amplification rather than a different zero inventory.

## 2. Arithmetic localization and independent reference

At electron perpendicular index 16 and the retained worst-relative frequency,
for n=−1 the moment has magnitude `3541.610980795522` and the Cauchy contribution
divided by kpar has magnitude `3540.586514537703`. Their difference has magnitude
about `1.025`. For n=+1 the corresponding Cauchy magnitude is
`3542.636038589406`. Harmonic contributions subsequently cancel further.
The original assignment and reduction order is retained; no harmonics are
paired across ranks.

The two local differences have magnitudes `1.0246476861` and `1.0252394326`.
This loses roughly 3.5 decimal digits before the subsequent harmonic cancellation.
It amplifies ordinary rounding in the Cauchy dot/cell reductions; there is no
evidence of stale species weights or a mismatch in the fixed-cell formulas.

The independent diagnostic reads the actual initialized double-precision grid
and df0 arrays. It uses actual production int_T/BESSJ coefficients, independent
real128 complex-log cell primitives, exact PWL polynomial moments, and real128
parallel/perpendicular/harmonic accumulation. Independent real128 direct-cell
and nodal-weight expressions agree at approximately 1e-30 on the retained local
checks. Neither double implementation is designated the reference.

Two reference variants distinguish arithmetic effects:

1. Freeze the actually constructed double numerator combinations and resonance
   momentum, then integrate/accumulate in real128.
2. Rebuild those combinations and the resonance in real128 from the identical
   double derivative, grid, frequency and physical inputs.

No analytic VDF derivatives are substituted. Bessel coefficients, pi and the
perpendicular spacing remain those of the actual double source; the reference
does not remove finite-grid, derivative or finite-support physical errors.
The four diagnostic double paths reproduce their actual MPI tensors within the
retained `2e-14` verification tolerance.

`reference_*.stages` separates direct cell accumulation and primitive error,
cached dot accumulation, product and weight-construction errors, moment and
numerator construction, perpendicular accumulation, and harmonic/reduction
accumulation. Unweighted errors have different units/scales and cannot be
compared as if they were final tensor errors. In the electron failed-coordinate
study, the largest original Cauchy discrepancies are of order 2e-10 before
physical weighting. Plain accumulation accounts for most of that error;
primitive/weight/product rounding remains after compensation. Moment and
numerator-combination rounding are common additional sources. The final
component/reference and conditional determinant comparisons below assess their
combined effect.

| Local error source, before tensor weighting | Maximum absolute discrepancy on the 15 electron frequencies |
|---|---:|
| Direct cell accumulation alone | 1.784e-10 |
| Rounded direct primitives, excluding accumulation | 2.102e-11 |
| Cached dot accumulation alone | 2.354e-10 |
| Rounded cached products | 5.035e-12 |
| Rounded weight coefficients | 1.827e-11 |
| Double polynomial moments | 1.627e-12 |
| Double numerator construction, nodal maximum | 1.122e-10 |

Perpendicular-sum discrepancies peak at `2.20e-13` in their pre-normalized
units; subsequent harmonic/group accumulation discrepancies are below
`4.51e-16` in tensor units. These maxima occur at different intermediate
quantities; they must not be added as a scalar determinant error budget.

The determinant reference is explicitly conditional: it replaces the measured
numerical tensor(s) in an actual wave matrix with the real128 tensor reference,
then evaluates a generic 3×3 determinant in real128. Analytic-species/vacuum
inputs and inherited wave-matrix construction error remain double precision.
It measures propagation of the studied arithmetic errors; it is not a second
production determinant or an independently quad-precision physical ALPS backend.

## 3. Alternatives and selected experimental change

**Compensated double reductions.** Componentwise Neumaier accumulation is applied
to cached weight products and uncached unchanged cell primitives. Moments,
perpendicular/harmonic sums and MPI reductions are initially left unchanged.
At i=16,n=−1 it reduces local errors from `1.34e-12` cached / `1.84e-12` direct
to `1.94e-14` in each path. It improves full tensors, but the complete E/B
trajectories still fail: maximum relative differences are respectively
`2.87064e-12` / `2.80973e-12`. These failures are retained. Compensating the local
moment in this slice provides no further improvement; it is not adopted.

**Restricted extended Cauchy arithmetic.** The second candidate mechanically
instantiates the unchanged C2c kernel with `selected_real_kind(18)` intermediates,
stores its reusable weights in that kind, and applies those weights in that
kind before returning through the existing complex(8) integral interface.
Uncached cell evaluation and accumulation use the same extended kind but retain
their independent direct-cell evaluation order.

On GNU Fortran 11.4 this is x87 real(10): an 18-decimal-digit type with a 64-bit
significand and 16-byte real storage. The independent real128 reference has
33 decimal digits. This is **not** wholesale real128 susceptibility evaluation.
Numerator combinations, PWL moments, polynomial-power recurrences, outer sums,
Bessel values, MPI reduction and determinant assembly remain double precision.

The only arithmetic transformations are in new validation files:
`arithmetic.f90`, `transform.awk`, `wide.f90`, `wide.awk`, and their builders.
The fixed-cell logarithm formula, far-field series, branch choices, support and
C2d tensor identities are unchanged. Exact cache keys and invalidation are
unchanged. There is no rounding-for-agreement, discarded contribution, empirical
correction, full-determinant interpolation, new input selector, or production
source change.

The extended path improves independently measured tensor accuracy as well as
agreement. On the 15 failed E frequencies, the largest frozen-reference tensor
errors reproduced in both final build configurations are approximately:

| Arithmetic | Cached | Direct |
|---|---:|---:|
| Original double | 3.13e-13 | 2.10e-13 |
| Compensated double | 1.28e-13 | 1.21e-13 |
| Extended Cauchy only | 3.62e-14 | 3.59e-14 |

These are maximum absolute component errors, not a claim of improvement at
every individual sample. Weak-component relative errors and absolute zz errors
are also retained. Common double numerator/moment errors remain, so bit identity
and arbitrarily high susceptibility accuracy are not asserted.

Rebuilding the numerator and resonance from the same derivative inputs in
real128 exposes the remaining common double arithmetic: maximum electron
component errors are `3.90e-13` / `3.08e-13` for original cached/direct, versus
`1.68e-13` for either extended path. The proton reference errors are around
`1.2e-15` for the selected method. Both references therefore support improvement,
but cached/direct agreement is appreciably better than absolute accuracy of the
whole double calculation.

For example, the conditional determinant error maxima on these 15 frequencies
are:

| Arithmetic | E: max relative, frozen / rebuilt reference | B: max relative, frozen / rebuilt reference |
|---|---:|---:|
| Original cached | 2.042e-11 / 2.048e-11 | 1.999e-11 / 1.999e-11 |
| Original direct | 1.103e-11 / 1.697e-11 | 1.084e-11 / 1.656e-11 |
| Extended cached | 1.639e-12 / 1.462e-11 | 1.667e-12 / 1.428e-11 |
| Extended direct | 1.639e-12 / 1.462e-11 | 1.667e-12 / 1.428e-11 |

Against the rebuilt reference, maximum absolute E determinant errors decrease
from `3.996e-12` / `7.222e-13` to `6.151e-13` / `6.170e-13`; corresponding B
values decrease from `4.013e-12` / `7.290e-13` to `6.221e-13` / `6.242e-13`.
The compensated candidate does not uniformly improve determinant errors even
when its tensor maxima improve, and is not selected. No universal `2e-12`
absolute-accuracy claim follows from passing the pairwise gate. Further
reorganization of common numerator/moment arithmetic is outside this task.

## 4. Focused acceptance and physical controls

The completed frozen acceptance run is
`/tmp/alps-c2f3a-acceptance-final-20261123`; the portable entry is
`bash validation/run_nyquist_c2f3a.sh`. It finished from start to finish with
`entry_status.txt: exit=0 line=1` and a final COMPLETE line. Its frozen driver was
not edited while running. Compact portable evidence is archived under
`tests/data/nyquist_c2f3a/`.

The completed gates are the original four strict comparisons, both 257-frequency
trajectories, original/new static and full-tensor fixtures, Cauchy primitives,
cache invalidation, independently perturbed CR checks, signed ordinary growing
roots, both standard A1 rectangles, a numerical-proton control, one low-frequency
raw/scaled exclusion box, and copied-source nonfinite failure/shutdown.
No full resolution/support matrix or unrelated analytic suite is repeated.

| Representation | Optimized maximum relative cached/direct D difference | Debug maximum |
|---|---:|---:|
| Numerical electrons, analytic protons | 1.23376230900643e-13 | 1.23376230900643e-13 |
| Both numerical species | 1.24242220541721e-13 | 1.24242220541721e-13 |

Maximum absolute pairwise differences on the complete trajectories are
`3.03478e-13` / `3.04262e-13`. Maximum cumulative phase differences are
`2.75335e-14` / `2.93099e-14` radians, with total phase differences zero and
`1.77636e-15`. Comparing old cached with new cached instead gives maximum
cumulative phase changes about `2.1e-11` radians and zero total-loop change.
No actual physical root or accepted count changes.

| Representation | New positive-frequency root | Maximum change of either partner from C2f.3 |
|---|---|---:|
| E | .396754648208647 + .0859766939815021i | 1.8433e-13 |
| B | .394870463511772 + .0858748871269201i | 1.2773e-13 |

Both negative-real-frequency partners were solved independently at positive
kpar. Maximum raw root residual is about `3.8e-14`; the existing residual and
reflection-symmetry gates pass. Optimized and debug roots agree in the retained
output. These are the established finite-grid roots, not new analytic references.

All eight fresh standard A1 calculations have status zero, winding two and three
accepted global passes:

| Representation | gamma=.01 rectangle: calls / refinements / depth | gamma=.001 rectangle |
|---|---:|---:|
| E | 628 / 90 / 6 | 704 / 128 / 9 |
| B | 630 / 91 / 6 | 702 / 127 / 9 |

Both builds reproduce those work counts. The representative both-numerical
low-frequency box is `[-5e-5,5e-5] × [5e-6,2e-5]`: raw D accepts zero with
468 calls; D/omega² accepts zero with 448 calls. Each has status zero and three
passes in both builds. The origin is excluded, so the scaling introduces no
enclosed zero or pole. No exclusion claim below the established conditioning
limit is added. The limited numerical-proton 257-point control differs from its
C2f.2/C2f.3 reference by at most `2.523e-14` relative D.

At `-.53125+.01i`, independently perturbed determinant CR residuals for step
`1e-3,1e-4,1e-5` are `7.30e-5,7.30e-7,1.50e-8` (E), and
`6.37e-5,6.38e-7,1.44e-8` (B). Both builds pass the unchanged strict CR checker.
All six component and summed-species absolute/norm-relative diagnostics are
retained; weak component relative residuals can be roundoff limited.

There are 56 four-rank jobs, 224 rank metric records and 16,024 actual collective
determinant calls. Per-job call counts agree across all four ranks. Actual grid,
derivative, Bessel-array and harmonic-assignment diagnostics compare byte for
byte between each original and extended cached calculation. E retains numerical
electron rank ranges 0..5 and 6..10; B retains proton ranges 0..7 and 8..15 and
electron range 0..10. No cross-rank harmonic pairing is introduced.

The two deliberate nonfinite tests return unresolved status 2, not an accepted
zero, after three collective calls and coordinated shutdown. Fixed-frequency
PROBE jobs intentionally exit through the retained adapter with summary status
6 after their requested samples; they are not failed adaptive A1 contours.
The compact `adaptive.dat` separates the twelve valid adaptive contours and two
deliberate failures from those diagnostic exits.

The mathematical preflight preserves exact zero in the coupled parallel static
fixture; the full tensor static contraction is approximately `9.89e-17`.
The historical naive static residuals `−6.41814e-4,−1.61609e-4,−7.19970e-5`
are reproduced as historical comparisons, not imposed on the coherent model.
Extended signed-affine Cauchy tests have maximum absolute error `4.71e-16` against
an independent real128 primitive. Original cache/invalidation and alternating
species/different-grid tests pass at their unchanged `2e-13` tolerance.

## 5. Performance and limitations

Cauchy kernel operation counts are unchanged: two basis evaluations per cached
cell/weight set, one direct evaluation per uncached cell. Extended complex
weights occupy 32 bytes instead of 16; exact double cache keys are not enlarged.
No additional explicit allocation or persistent equilibrium cache is introduced.
Uncached conversion temporaries scale with the parallel grid length. The final
IO-disabled 257-frequency benchmarks separate the incremental cached cost from
the speedup relative to uncached coherent integration and from C2b production.

The measured optimized collective determinant times are maxima over rank metrics
(initialization excluded). Cached rows have two sequential repeats; direct rows
have one IO-disabled repeat, with the earlier diagnostic trajectory providing
an additional consistency observation rather than an IO-matched repeat.

| Case / arithmetic | Seconds per 257 frequencies | Mean seconds per determinant |
|---|---:|---:|
| E original cached | 7.582–7.595 | .029528 |
| E extended cached | 7.505–7.660 | .029504 |
| E original direct | 50.054 | .194764 |
| E extended direct | 80.853 | .314601 |
| B original cached | 14.582–14.824 | .057210 |
| B extended cached | 14.394–14.398 | .056015 |
| B original direct | 96.073 | .373825 |
| B extended direct | 157.848 | .614195 |

The cached change is indistinguishable from ordinary timing variation for E
and about −2% in the measured B pair; no optimization speedup is claimed.
Extended uncached evaluation costs about 62–64% more. Reuse is about 10.7× / 11.0×
faster than the **extended** direct method, whose increased cost enlarges that
ratio. Against the original direct implementation, the new cached path retains
about 6.6× / 6.7× benefit, consistent with the prior approximately 6.6× reference.
This study does not rebenchmark C2b or imply that coherent integration is faster
than production C2b.

Per trajectory, E retains 1,360,044 cached kernel-cell operations and 5,397 weight
sets; B retains 3,367,728 and 13,364. Direct counts remain 250,928,118 and
624,260,196. Reported live cache footprints summed over active ranks increase
from 87,544 to 172,888 bytes (E) and 172,036 to 340,692 bytes (B).
These are array-counter footprints, not process peak-RSS measurements. No new
explicit allocation is added, but wider weight storage and conversion/expression
temporaries are real costs. All numerator-read counters remain unchanged.

The selected extended type and performance are compiler-dependent. Portability
to a compiler whose first type with 18 digits is real128 has not been qualified.
The represented susceptibility remains mathematically holomorphic in the upper
half-plane: frequency-independent rounding does not introduce a new integration
switch. Finite-difference CR tests assess numerical evaluation separately.
This study makes no damped-continuation, production-adoption, electron-resolution,
non-Maxwellian or broader-wavevector claim.

## 6. Reproduction, evidence and disposition

Run:

```sh
C2F3A_SCRATCH=/tmp/alps-c2f3a-new bash validation/run_nyquist_c2f3a.sh
```

Use a fresh scratch directory. The entry freezes sources/fixtures, regenerates
only the necessary 64×128 VDFs, executes MPI jobs sequentially, retains the
original failed gates, and checks production/historical-source preservation.
GNU Fortran 11.4.0, Open MPI 4.1.2 and the existing system BLAS/LAPACK are used.
Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`;
debug uses `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
Toolchain, linked-library paths, intrinsic Bessel symbols, Git revision,
pre-existing changes and source SHA-256 identities are archived.

Both configurations pass original and modified mathematical fixtures and all
focused physical gates. Existing original-source/report/archive checksums and
the final frozen-source manifest verify successfully. The four old strict
failures and the two compensated failures remain explicitly recorded as such;
they are not silently reclassified as passes. No inherited unrelated harness
failure has been repaired or claimed fixed.

The preliminary orchestration attempt was stopped after its active MPI job
finished because the new timing reducer selected the wrong metrics column.
That reporting bug was corrected and smoke-tested before freezing the final
entry. The interrupted attempt is not counted as acceptance. The complete final
entry was then executed cleanly without editing its running driver.

**Recommendation:** retain the restricted extended Cauchy arithmetic in the
experimental coherent validation path. It improves measured accuracy against
two independent real128 arithmetic references, restores all four strict
comparisons and preserves the required mathematical/physical/MPI controls with
practical weight reuse. Keep production C2b and all defaults unchanged. A future
production review must address compiler precision support and the separate
coherent susceptibility/causal-continuation requirements; this task does not
authorize that adoption or further C2f/C3 studies.
