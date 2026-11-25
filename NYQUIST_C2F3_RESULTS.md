# C2f.3 — Numerical electrons and simultaneously numerical species

**Scientific disposition: partial success.** Optimized and debug physical results
converge toward the analytic bi-Maxwellian reference. The coherent method removes
the resolved C2b electron artifact on the tested contours. The existing strict
cached/uncached determinant comparison fails for electrons and both-numerical
species; its tolerance has not been changed. Independent real128 diagnostics
identify amplified floating-point cancellation. Production remains unchanged.

**Execution status: complete.** The final frozen entry executed cleanly from
start to finish at `/tmp/alps-c2f3-acceptance-final-20261123`, exited zero, and
completed optimized/debug physics, mathematical checks, all twelve analytic
compatibility suites, collection, reductions and source-preservation checks.
It ran on November 23, 2026, approximately 15:15–18:04 America/Phoenix
(2 h 48 min, including builds and diagnostics). Executor completion does not
override the four retained strict numerical-comparison failures.

## 1. Baseline and experimental construction

Git revision: `99a087271928846cb8dbc49507d5a3a0f4d68d84`, with pre-existing working-tree
changes preserved. Qualified production SHA-256 values:

| Source | SHA-256 |
|---|---|
| `src/ALPS_fns.f90` | `0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488` |
| `src/ALPS_fns_rel.f90` | `c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6` |

The entry independently reapplies the C2f.1 patches to retained before sources,
then requires exact agreement with current production. No production source,
default, input selector, derivative, continuation, root solver, harmonic cutoff,
MPI assignment, Nyquist kernel or growth-search algorithm is changed.

`validation/c2f3/build.sh` delegates to C2f.2's builder, reusing unchanged C2d
`helpers_cached.inc`, `model.f90`, `work.f90` and C2c's Cauchy kernel. The uncached
comparison uses the preserved C2d uncached helper. Only observational startup
instrumentation is added. Current variants link corrected production BESSJ
**directly**: symbol inspection identifies the shared module routine and its
intrinsic `jn/j0/j1` dependencies; no accurate-array substitution is present.

The coherent positive-imaginary-frequency path combines the n=0 longitudinal
resonant and nonresonant terms once. Its `coherent_ee` returns the complementary
`(1-J0²)E` moment. It bypasses the ordinary `int_ee` body; the corrected production
corner is neither reversed nor added a second time. Fresh proton-only matched
trajectories reproduce all three archived C2f.2 optimized trajectories exactly.

All changes are validation files under `validation/c2f3/`, its shell entry, this
report and new `tests/data/nyquist_c2f3/` artifacts. The historical C1/C2d/C2e/C2f.2
reports and scientific fixtures remain unchanged.

Validation file responsibilities:

| Files under `validation/c2f3/` | Role |
|---|---|
| `snapshot.sh`, `driver.sh` and top-level `run_nyquist_c2f3.sh` | Frozen execution, provenance, build configurations and suite orchestration |
| `build.sh`, `audit.inc`, `fault.sh` | Reuse C2f.2 builders, observe actual initialization, inject one post-collective failure |
| `cases.sh`, `prepare.sh`, `run_case.sh`, `physical.sh`, `gate_A.sh` | Normal VDF generation, sequential jobs, staged physical comparisons and actual call checks |
| `analyze.f90`, `qualify_roots.f90`, `cr_absolute.f90` | Matched tensors, independent matrix determinant, ordinary-root qualification and derivative diagnostics |
| `math.sh`, `tools.sh`, `test_species_weights.f90`, `test_electron_roundoff.f90`, `support.f90` | Existing mathematical gates, species-key tests, actual-data real128 arithmetic and exact support moments |
| `decompose.sh`, `collect.sh`, `review.sh`, `archive.sh` | Read-only return-piece attribution, compact reductions and lossless evidence packaging |

## 2. Historical C1 reconstruction

The 64×128 C1 `electron_only.in` and generator input compare byte for byte with
retained fixtures. Normal generation also produces byte-identical electron data
under the original and new case names. Historical physics sources are retained
under `tests/data/nyquist_c2f3/historical_c1/`, with hashes and provenance. The
historical numerical source identity agrees with the independently recorded C2b
before-source hash. It is linked in a separate build with current common startup
and observational adapters; this is not presented as an archived full C1 binary.

That reconstruction returns **W=3, status 0, three passes, 586 disp calls**, exactly
reproducing the C1 diagnostic. Corrected production C2b still returns three at
64×128. Corrected production Bessel and corner changes alone therefore do not
resolve the original electron discrepancy.

## 3. Electron normalization, support and derivative audit

The unchanged equilibrium has proton beta 3 and anisotropy 2, electron beta 3 and
anisotropy 1, equal densities, no drift, `me=0.000544617021487`, `qe=-1`,
`vA/c=1e-4`, and `(kperp,kpar)=(0.001,0.3)`. All physical tests retain positive
kpar and independently solve signed-real-frequency partners.
Frequencies are in reference-proton gyrofrequency units, wavevectors in inverse
proton inertial-length units, and momenta in reference-proton Alfvén momentum
units throughout the comparison.

Electron S4 support is `pperp ∈ [0,0.16168369438931063]`,
`ppar ∈ [-0.16168369438931063,0.16168369438931063]` in reference-proton momentum
units. At 64×128 both spacings are approximately `0.002526307724832979`.
Grid labels denote ALPS interval counts: full arrays include endpoints, while
normal centered numerical derivatives use interior indices. The audit reads
actual broadcast pp and df0, checks generated-file coordinates exactly, and never
substitutes analytic derivatives. Normal method-2/log continuation initialization,
kind 1/order 8 fits, differentiation and harmonic setup are retained.

Exact finite-support integration independently gives S4 density
`.9999998720475691`, conditional parallel beta `2.999998476211816`, and
perpendicular beta `2.999994598311006`. Actual generated-grid trapezoidal moments:

| Grid | Density | Parallel beta | Perpendicular beta |
|---|---:|---:|---:|
| 32×64 | .997391610930386 | 2.999998362597994 | 3.007851753671465 |
| 64×128 | .999348573248434 | 2.999998447534884 | 3.001950426961898 |
| 96×192 | .999710468790161 | 2.999998463443704 | 3.000863163118764 |
| 128×256 | .999837095085248 | 2.999998469025323 | 3.000483029867505 |

Measured drifts are consistent with roundoff. Full rectangle-rule and trapezoidal
moments are retained; these full-grid moments are not mislabeled as the exact
moments of the integrator's interior derivative representation.

For electron grids 32,64,96 the relative L2 derivative errors (perpendicular,
parallel) are `(0.00897912,0.01003096)`, `(0.00225264,0.00251803)`, and
`(0.00100182,0.00111997)`. The isotropic identity residual
`S=u*dfpar-p*dfperp` has relative L2 norms `.00712797,.00178272,.000792377`.
It is a real finite-difference representation error and decreases approximately
quadratically; it has not been silently removed or renormalized.

Stored Bessel values, signed arguments, halos, source mass/charge and assignments
are audited on workers. Stored values agree exactly with corrected production
BESSJ. Primary harmonic limits remain proton 15 and electron 10. The electron
finite-sum J² defect is about `1.13e-16`, with predicted static longitudinal
Bessel error about `-4.33e-18` on the representative 64×128 audit.

## 4. Electron tensor and determinant convergence

All six returned ALPS species-tensor components are compared at identical complex
frequencies. These are source-normalized ALPS tensor contributions, not a newly
introduced conventional dielectric normalization. Absolute complex errors are
retained. Relative component errors are omitted when the reference is below
`sqrt(epsilon)` times the species tensor norm; no scientific acceptance tolerance
is changed by that reporting convention.

On the original 257-point matched contour:

| Electron grid | Median relative D error | Maximum relative D error | Maximum absolute D error | Maximum relative electron zz error |
|---|---:|---:|---:|---:|
| 32×64 | .00455931 | .448094 | .0346509 | .0183984 |
| 64×128 | .00107673 | .113343 | .00821549 | .00459520 |
| 96×192 | .000573430 | .0504093 | .00439844 | .00206107 |
| 128×256 | .000287658 | .0283140 | .00218479 | .00115764 |

Maximum absolute xx errors decrease approximately `.001412,.0003524,.0001566,
.00008807`; near-zero mixed components are assessed with absolute errors. Maximum
cumulative phase errors decrease approximately `.4671,.1079,.04683,.02601` radians;
total-loop phase errors are at roundoff. Complete lower-contour and component
statistics are retained. Finite sampled agreement is not a proof of a uniform
boundary error bound.

At fixed spacing, S3/48×96, S4/64×128 and S6/96×192 give original-contour median
relative D errors `.00587191,.00107673,.00106978`. Extending four to six widths
changes the result much less than grid refinement. S3 has a larger median but a
smaller worst-case relative error, demonstrating metric-dependent cancellation;
truncating the distribution is not recommended as an accuracy correction.
Refining a fixed-support grid retains a finite-support error floor; the observed
trend does not imply that resolution alone exactly recovers the infinite-support
analytic distribution.

## 5. Electron windings and the extra C2b zero

The original rectangle is `[-1,1] × [0.01,0.2]`; the lower rectangle changes only
the bottom to `.001`. Sampling controls and three-pass acceptance are unchanged.

| Electron grid | C2b W, original/lower | Coherent W, original/lower | Coherent calls, original/lower |
|---|---|---|---|
| 32×64 | 2 / 2 | 2 / 2 | 634 / 720 |
| 64×128 | 3 / 3 | 2 / 2 | 628 / 704 |
| 96×192 | 3 / 3 | 2 / 2 | 630 / 706 |
| 128×256 | not run | 2 / 2 | 630 / 694 |

Both builds reproduce these coherent counts with status zero and three passes.
The analytic reference independently gives two on both rectangles. C2b's coarse
winding match does not establish accuracy: its original-contour median relative
D error is about `.1559` at 32×64, versus `.004559` for coherent.

The ordinary C2b solver resolves extra nearly imaginary candidates:

| Grid | Real frequency | Growth | Absolute D residual |
|---|---:|---:|---:|
| 64×128 | 1.34320e-13 | .0360215233715567 | 1.41202e-14 |
| 96×192 | -1.79791e-13 | .0127767817142450 | 9.54e-15 |

Unchanged local boxes with real bounds ±.005 and imaginary bounds
`gamma_candidate ± .003` independently accept W=1 for C2b and W=0 for coherent,
with 448 calls and three passes each. Coherent D at the candidate frequencies is
approximately `-.00174522` and `-.000239360`, respectively. The finite-frequency
matrix qualification and local contours avoid mistaking tiny near-origin D for
a resolved zero. No zero was bracketed on the sampled C2b 32×64 imaginary axis
from .005 through .2; this is not an exclusion outside that interval.

Actual electron n=0 zz plus the once-counted nonresonant contribution at 64×128
is `+.01602502379964` for C2b versus `-7.36889305e-6` for coherent at i.001;
at i.01 it is `+.01545363883633` versus `-.00073675023347`. The external electron
charge factor is included in both pieces. The electron Tlim switch height is
about 5.51, so C2b uses its small-imaginary approximation throughout these boxes.
The observed discrepancy is consistent with its previously diagnosed finite-core
and cancellation limitations. C2f.3 makes no repair to that approximation.
A converged C2b boundary winding is not automatically a physical root count.

## 6. Electron isotropic control

With isotropic analytic protons, analytic W=0 is independently established on
both standard rectangles. Coherent numerical electrons return zero at all three
grids and both heights, in both builds. C2b gives 0/0 at 32×64 and 1/1 at both
64×128 and 96×192. The electron tensors in unstable-proton and isotropic-proton
64×128 cases match exactly at every common frequency; no different electron
input explains the isotropic extra winding. These are finite-domain statements.

## 7. Both-species assembly and static structure

All four analytic/numerical representation combinations use identical physical
parameters and matched coordinates. Species tensors are added before the wave
matrix determinant is taken. No species determinants are added.
The generated numerical proton and electron files compare byte for byte between
their separate-species and both-species cases at every primary grid size.

Native both-numerical tensors agree exactly with the separately measured species
tensors at all three grids. Independent generic 3×3 determinants differ from the
native determinant by at most `8.75e-15,6.66e-15,5.77e-15`, respectively. The
comparison accounts for matrix perturbation/roundoff near cancellation rather
than demanding relative accuracy for an arbitrarily small D.

Actual-source coupled parallel identity, static zero, PWL moment, six-component
full-tensor identity and static Ward tests pass optimized and debug. Electron
small-positive-frequency contractions retain the coupled structure and no old
Bessel static floor appears. Historical audit fields describing a hypothetical
old corner defect are not interpreted as an extra correction to coherent E.

## 8. Independent proton/electron refinement and its limitation

ALPS has one global `nperp,npar` pair. Changing that architecture is outside scope.
Native both-species adaptive contours and ordinary roots therefore use equal
grid dimensions. Independent refinement uses actual species tensors at precisely
matched frequencies: first hold proton 96×192 and refine electrons, then hold
electron 96×192 and refine protons, followed by native diagonal runs. The complete
3×3 tensor-combination matrix is retained and diagonal assembly is checked.

Original-contour median relative D errors at fixed proton 96 are
`.00545898,.00227809,.00182018` for electron 32,64,96. At fixed electron 96 they
are `.0130674,.00357912,.00182018` for proton 32,64,96. Species tensor errors are
retained independently; cancellation of total D error alone is not accepted as
accuracy evidence.

Off-diagonal outputs are matched-frequency diagnostics, **not** independent A1
counts or ordinary root solves on unequal grids. No determinant interpolation,
new susceptibility or hidden unequal-grid input implementation is introduced.

## 9. Physical growing roots

The unchanged ordinary solver independently resolves both real-frequency partners
at positive kpar, with `D_threshold=1e-13`. The analytic positive root is
`.39666286088749 + .08598637345950 i`.

| Numerical representation | Grid | Positive real frequency | Growth | Complex root error |
|---|---|---:|---:|---:|
| Electrons | 32×64 | .397030221485648 | .0859472898085411 | 3.69434e-4 |
| Electrons | 64×128 | .396754648208787 | .0859766939813827 | 9.22963e-5 |
| Electrons | 96×192 | .396703465239665 | .0859821781452990 | 4.08205e-5 |
| Electrons | 128×256 | .396685545682750 | .0859840974509767 | 2.27987e-5 |
| Both | 32×64 | .389520320307861 | .0855823918788390 | 7.15396e-3 |
| Both | 64×128 | .394870463511675 | .0858748871270028 | 1.79586e-3 |
| Both | 96×192 | .395865727275250 | .0859359093121872 | 7.98729e-4 |

The observed error decreases approximately quadratically, without imposing that
order as a universal VDF guarantee. Negative partners, complex residuals, tensor
errors at identical frequencies, and matrix-scaled qualification are retained.
Only strictly upper-half-plane coherent evaluations are interpreted physically;
no coherent damped response is claimed.

For comparison, the corrected-production C2b electron-only positive roots are
`.397280640631510+.0857337643591784i`,
`.396817116718169+.0859236159099872i`, and
`.396731217995583+.0859586115368088i` at the three grids. Those physical-pair
solutions also improve with refinement, despite the independently resolved
additional C2b zero. Ordinary convergence to the known pair therefore does not
validate the complete computed inventory.

## 10. Both-numerical contours

| Grid | Coherent W, original/lower | Calls, original/lower | Isotropic W, original/lower | C2b W, original/lower |
|---|---|---|---|---|
| 32×64 | 2 / 2 | 624 / 694 | 0 / 0 | 2 / 2 |
| 64×128 | 2 / 2 | 630 / 702 | 0 / 0 | 3 / 3 |
| 96×192 | 2 / 2 | 630 / 702 | 0 / 0 | 3 / 3 |

Both compiler configurations reproduce status zero and three passes for all these
physical contours. Original-contour coherent median relative D errors are
`.0153522,.00393766,.00182018`; maximum relative errors are
`.200911,.0486920,.0215573`, and maximum absolute errors are approximately
`.118019,.0302749,.0140396`. All component errors, sampled minima and phases are
retained. The two independently solved roots and accepted two-count inventory
agree within these user-prescribed domains; no global stability claim follows.

## 11. Low frequency, conditioning and analyticity

At 64×128, both electron-only and both-numerical coherent cases independently
accept W=0, status zero and three passes on all existing narrow, broad and deep
C2e windows, for raw D and diagnostic D/omega². This gives 24 accepted exclusions
across two builds. The boxes exclude the origin, so scaling has no enclosed pole
or zero and preserves winding. The deepest lower bound is 1e-8.

At i1e-8 the conservative one-ulp input-matrix D bound is about 4.7% of |D| for
both representations; at i1e-9 it exceeds 200%. Counts arbitrarily closer to the
origin are not certified. Scaled reciprocal conditioning, absolute D, diagnostic
D/omega² and raw-versus-quad matrix determinants are retained.

For fixed derivative data/support, coherent numerator frequency dependence is
polynomial and Cauchy weights are analytic in Im(a)>0. Positive kpar maps positive
Im(omega) there for either charge. No Tlim, moving window or resonance-entry
integration switch is introduced. All four coherent determinant CR studies
(electron/both × optimized/debug) pass the existing strict C2c checker.
Near i.001, electron-only D residuals fall `4.09e-5 → 4.09e-7 → 4.19e-9` with
shrinking perturbation. Weak component relative CR residuals can be ill-conditioned:
actual xz derivative discrepancies correspond to differences near machine epsilon
in tensor values. Absolute residuals are retained rather than hiding these floors.
Mathematical holomorphy of the represented model, floating-point evaluation and
convergence to the physical VDF remain distinct claims.

## 12. Weight reuse, strict reproducibility exception and performance

Existing frequency/species/grid/kpar invalidation tests and a new alternating
mass-scaled species test pass. The latter uses different grid lengths and signed
charges, with maximum absolute kernel/weight discrepancy `4.45e-15`. No approximate
frequency matching, conjugated application or determinant interpolation is used.

However, the **unchanged** C2f.2 `2e-12` relative-D comparison fails in both builds:

| Case | Median relative cached/uncached D difference | Maximum relative | Maximum absolute |
|---|---:|---:|---:|
| Electron-only | 5.20e-13 | 9.40378e-12 | 9.19562e-12 |
| Both | 5.21e-13 | 9.17169e-12 | 9.21052e-12 |

Fifteen electron-only matched rows exceed the retained criterion. At the worst
relative point, `-.53125+.01i`, absolute D difference is `9.53e-14`. Electron tensor
differences are below `4.8e-13` absolute, dominated by transverse components;
zz differs by at most about `3.6e-15`. Total-loop phases agree and maximum
cumulative phase differences are below `1e-11` radians.

A read-only real128 calculation on actual initialization derivatives diagnoses
cancellation: at perpendicular index 16, terms of magnitude about 3542 subtract
to about 1.025 for each n=±1, followed by opposite-harmonic cancellation. Cached
and uncached double results each differ from independent real128 by roughly
`1.4e-13–1.8e-12` before tensor weighting, a few ulps of the large operands.
The data identify arithmetic amplification without evidence of stale species
weights. Neither evaluation order is designated exact. The strict gates remain
**FAIL**; no tolerance or integrator is modified to obtain a pass.

Two sequential IO-disabled 257-frequency optimized repeats measured:

| Case | Implementation | Seconds/trajectory, range | Mean seconds/disp |
|---|---|---:|---:|
| Electrons | Cached coherent | 7.545–7.693 | .02965 |
| Electrons | Uncached coherent | 49.920–50.135 | .19466 |
| Electrons | Production C2b | 4.155–4.267 | .01638 |
| Both | Cached coherent | 14.688–14.724 | .05722 |
| Both | Uncached coherent | 97.640–97.677 | .37999 |
| Both | Production C2b | 9.926–10.412 | .03957 |

Reuse accelerates coherent integration about 6.57× and 6.64×, respectively, but
coherent is about 1.81× and 1.45× slower than production C2b. These are different
comparisons and neither establishes general non-Maxwellian performance.

Across workers, cached electron trajectories use 1,360,044 Cauchy-cell kernel
calls, 5,397 weight sets and 87,544 cache bytes (43,772 maximum/rank); both cases
use 3,367,728 calls, 13,364 sets and 172,036 bytes (64,132 maximum/rank).
Uncached kernel counts are 250,928,118 and 624,260,196. Existing counters record
505,839,222 and 1,258,429,284 resonant-helper derivative-node reads, respectively;
they are instrumentation counters, not a hardware measure of every memory read.
Bessel construction/cutoff timing is separate, occurs once per worker startup,
and is not charged as per-frequency recomputation.

## 13. Regression, MPI and reproduction

Run `bash validation/run_nyquist_c2f3.sh`; optionally set `C2F3_SCRATCH` to a fresh
directory outside the checkout. It freezes sources, validation programs and
fixtures before compilation and executes four-rank jobs sequentially. No old
/tmp directory is a dependency. GNU Fortran 11.4.0 and Open MPI 4.1.2 are used
with system BLAS/LAPACK; exact library paths, flags and linked symbols are saved.
Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug
uses `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.

Completed scientific execution contains 454 study MPI jobs and 165 adaptive
contours: 161 accepted, two deliberately nonfinite failures and two evaluation-
budget failures. Injected cases retain no accepted integer, use 3 or 10 physical
calls on every rank, and finalize collectively. No unexpected evaluator failure
was observed. Rank call counts agree for all study jobs; initial assignments,
array halos, cache work and species construction are recorded. No communicator,
harmonic distribution or work-scheduling change is made.
The existing fixed-frequency `PROBE` adapter deliberately exits 1 after
coordinated finalization; the wrapper requires that observer exit, complete
data and equal four-rank counters. Those expected probe exits are not numerical
evaluator failures or substitutes for an accepted contour.

Optimized/debug Cauchy primitives, PWL convergence/moments, actual-source static
and tensor identities, production Bessel API and cache invalidation pass.
All twelve analytic compatibility suite executions pass: A1/A2, A3 growth, log
growth, mirror, line and plane, including their established legacy scientific
comparisons. Fresh numerical-proton production controls reproduce
`.395411725171728+.086220209204739i` and
`.405337296518808-.001709481527499i`, confirming the default remains C2b.
Inherited unrelated C2f.1 harness defects are not claimed repaired or fully rerun.

Coverage distinction: fresh isotropic adaptive grids cover fully analytic,
analytic-proton/numerical-electron, and both-numerical cases. The numerical-
proton/analytic-electron isotropic Nyquist controls remain the archived C2f.2
controls; they were not repeated as an additional grid in this entry. Fresh
proton-only unstable contours, matched tensors and growing/damped ordinary
production controls were repeated as described above.

Execution attempts are distinguished explicitly:

- The first attempt was stopped for a sign error in optional return-piece
  postprocessing; actual D/tensors were unaffected. It is marked ABORTED.
- A second stopped at the strict electron cache gate; that failure evidence is
  preserved under `tests/data/nyquist_c2f3/cache_diagnostic/`.
- `/tmp/alps-c2f3-complete-20261123` completed optimized/debug physics and all twelve
  compatibility suites, then failed packaging on an absent optional matched file.
  Corrected collection and read-only reductions were rerun successfully against
  all completed results. Scientific tolerances and sources were not changed.
- `/tmp/alps-c2f3-acceptance-final-20261123` is the fresh end-to-end replay with
  optional-file packaging fixed and preflighted. **It completed from start to
  finish with executor exit zero.** Both production/source-history preservation
  checks and the frozen source-manifest check passed. Its scientific acceptance
  summary remains PARTIAL because all four strict cache comparisons fail.

The frozen driver remained unedited throughout execution. Executor completion and
scientific disposition are separate: `acceptance.txt` reports PARTIAL while any
retained strict cache comparison fails. Original-source and historical-report
hash checks are repeated at completion. Compact tables, normal generator inputs,
matched coordinates, tensor examples, counter data and reproduction programs
are retained; large quadrature dumps remain in isolated scratch.
The final compact archive is under `tests/data/nyquist_c2f3/results/`; its
`SHA256SUMS` covers the packaged evidence. The executed entry and all its frozen
C2f.3 scientific scripts match the delivered versions. Ordered adaptive samples
and complete matched complex tensor trajectories are retained losslessly.

## 14. Disposition and next work

The principal physical objectives succeed for this bi-Maxwellian reference:
coherent electron and both-species tensors/roots converge, both standard unstable
boxes count two, isotropic controls count zero, and tested resolved low-frequency
windows have no extra count. Historical C1's W=3 is reproduced and separated from
the coherent result. No new susceptibility formulation was needed in this study.

Overall disposition remains **partial success** because the strict old
cached/uncached relative-D gate is exceeded and native unequal-grid ordinary
root studies cannot be performed through ALPS's global dimension interface.
Both limitations are exposed rather than hidden by tolerance changes or
interpolated determinants. Before production adoption, review cancellation-aware
arithmetic for far electron harmonics and the separate consistent damped-
continuation requirements. Any such correction is a new task, not implemented
here. Broader physical validation requires review of these results first.

Production remains C2b; no selector, coherent damped backend, numerical-VDF scan,
maximum-growth search, non-Maxwellian case or new oblique study is introduced.
