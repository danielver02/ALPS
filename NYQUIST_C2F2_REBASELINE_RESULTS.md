# C2f.2 — Corrected-production coherent susceptibility rebaseline

Result: the corrected-production coherent susceptibility reproduces C2e.
Optimized/debug scientific, compatibility and preservation gates passed.
Production source and defaults remain unchanged.

## Scope and source baseline

This is a compatibility revalidation of the experimental C2d **upper-half-plane**
response, not a production integration-method change. All new machinery is under
`validation/c2f2/` with entry point `validation/run_nyquist_c2f2.sh`.
Production source, existing C2d/C2e programs, historical reports and scientific
reference tables are protected by before/after SHA-256 checks.

The identified Git revision is `99a087271928846cb8dbc49507d5a3a0f4d68d84`.
The workspace includes substantial pre-existing modifications; the revision
alone does not identify the tested science. The scratch provenance directory
records the working-tree state, prior diff, source hashes and linked objects.
Relevant current production identities are:

| Source | SHA-256 |
|---|---|
| `src/ALPS_fns.f90` | `0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488` |
| `src/ALPS_fns_rel.f90` | `c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6` |
| C2c `cauchy.f90` | `9d06e2a07622c67b8259e95294097a12c7dc9a046d6c16504b9bbe4d147d6eba` |
| C2d `model.f90` | `1a3c7823f0ccc5e29a159c53a6edf90a92fd5dc81473b5ef8e5c8e46c3c920e9` |
| C2d `work.f90` | `818b6071eafd4332ecae2d1017878cb527cb5eabdccc91c5724fe2ce4b2faca8` |
| C2d `helpers_cached.inc` | `1b78911f018b395433a2a8b96d7d0bf2441ad39d16448865089a39e8e9d13f63` |

`provenance.sh` compares the existing C2e manifest without updating it. Only
the two explicit, independently qualified C2f.1 corrections may differ.
The original numerical-source hash is
`905416da5e1a2dfdb03974bb49e167595b1f11ab70e7c1b77ca1a3c1f0058435`;
the original real-Bessel module hash is
`8a8b7e3260681afb0cefe59915101d57a9d99be1f40120214ffe5d75a85f79ee`.
Both match the exact retained before files in `tests/data/bessel_intee/before/`.
Applying the separate qualified C2f.1 patches to those files reproduces current
production byte for byte. Any additional baseline difference is an error.

## Source construction and actual Bessel linkage

The historical C2d `variant.awk` inserts `integrate_coherent` into
`full_integrate` only for a nonrelativistic species and strictly positive
imaginary frequency. Both the nonresonant and resonant pieces then use the
same PWL **derivative** basis, exact polynomial moments and C2c Cauchy kernel.
Cached helpers call `c2d_work:prepare` for species/harmonic weights.

The historical C2e `variants.awk` can additionally replace stored-array BESSJ
calls with `c2e_bessel`, select intrinsic values only for numerical
nonrelativistic species, change the original corner factor, and optionally
add eight diagnostic harmonics. It never changes shared BESSJ itself.
Those original builders remain unchanged; blindly applying their corner
substitution to current production would fail its original-source match guard.
That is an obsolete construction assumption, not a production defect.

The new builder makes the variants explicit:

| Variant | Real Bessel object | Stored-array substitution | Susceptibility |
|---|---|---|---|
| current | Current corrected production | None; no local Bessel wrapper | Cached coherent C2d |
| historical | Exact pre-C2f.1 module | Original C2e intrinsic-array substitution | Same cached coherent C2d |
| uncached | Current corrected production | None | Original uncached coherent C2d |
| extra | Current corrected production | None | Cached coherent C2d, diagnostic +8 orders |
| production | Current corrected production | None | Default C2b, observation only |

Copied counters, constructor snapshots and timing are observational. The
current source has no `c2e_bessel` call or definition and no local BESSJ
implementation. Its numerical object has an unresolved reference to the
shared module's BESSJ, and that linked production object refers to `jn`,
`j0`, `j1`. The public production API is also executed against the independent
616-row real128 C2f.1 references, and actual stored arrays/halos are compared
with the reconstructed historical accurate constructor. Compilation alone
is not used as proof of which implementation supplies the values.

## Why the corrected corner is not counted twice

For positive imaginary frequency the copied `int_ee` returns `coherent_ee()`
before reaching its ordinary production body. That body retains the C2f.1
correction, and remains authoritative in the production control and outside
the coherent upper-half-plane dispatch.

For each perpendicular node, let `E=q/m*(u*M1(B)-M2(A))`. The `n=0,zz`
resonant contribution combines `J0^2*E` algebraically before numerical
integration, leaving the explicit-frequency expression
`-q*omega*u*J0^2/kpar * (M0(B)+a*C0(B))`. `coherent_ee()` contributes only
`(1-J0^2)*E`. Thus E is included exactly once, including its finite-harmonic
complement; no independent old/new corner offset is added. The underlying
six-mode tensor and determinant assembly in `disp` are unchanged.

There are two source-level calls to `int_ee` on the n=0 owner: one updates
`schi(s,3,3)`, the other updates the separate diagnostic `schi_low(s,3,3,0)`.
Both are intercepted by the same early return. This repeats computation for
bookkeeping; `schi_low` is reduced separately and is not added to `chi` a
second time. The physical tensor therefore contains one complementary E term,
not two. This distinction also explains repeated local moment reads in work
counters without implying doubled susceptibility.

The strongest compatibility comparison uses historical original-corner sources
and current corrected-corner sources with identical coherent helpers. At the
same frequencies, all six complex species components and the assembled wave
matrix are compared, in addition to D. The unreachable original versus
corrected corner cannot change that upper-half-plane result. No original
production factor is restored, and the coherent moment receives no second fix.

## Mathematical gates

The unchanged C2c tests exercise exact fixed-cell primitives, constant and
complex-linear data, curved polynomials/Gaussians, signed imaginary
resonances, endpoint/support behavior and stable complex-log evaluation.
The C2d tests exercise exact polynomial moments, second-order PWL convergence,
source-level tensor identities, the static parallel identity and weight reuse.
Original tolerances are retained.

The actual-source parallel fixtures reproduce:

| Parallel intervals | Historical naive C2c residual | Coupled residual |
|---:|---:|---:|
| 64 | -6.418143596134129e-4 | 0 |
| 128 | -1.616090865725894e-4 | 0 |
| 192 | -7.199702955351661e-5 | 0 |

The maximum pointwise continuum-identity residual is `2.844e-14`. These
nonzero naive residuals remain historical observations, not corrected-model
expectations. The copied corrected-corner fixture checks the common nodal
continuum rule with the original `2e-13` absolute tolerance; its measured
rounding difference is `1.42e-14`. It does not rewrite C2d's historical
corner-attribution fixture.

The full-tensor scalar tests retain all six source rearrangements and the
static Ward limit. Complex-frequency perturbations independently test tensor
and determinant derivatives; the physical groups include internal momentum
locations, both sides of the support, near-axis points and the growing-root
region. Frequency-independent Bessel coefficients do not alter the
holomorphy argument: fixed-domain Cauchy transforms, polynomial frequency
factors and finite harmonic sums have no new frequency-dependent selector.
Numerical CR residuals are interpreted with truncation and roundoff, not as
proof that a floating-point function is exactly analytic.

## Physical fixtures and measured results

The normal ALPS generator prepares numerical protons on 32×64, 64×128 and
96×192 grids, with four-width support, beta_parallel,p=3 and anisotropy 2.
Electrons are analytic isotropic bi-Maxwellians with beta_parallel,e=3 and
mass `5.44617021487e-4`. Densities are one, drifts zero, vA/c=`1e-4`,
`(kperp,kpar)=(.001,.3)`. Normal numerical differentiation and method-2
continuation initialization are retained. No numerical electrons, non-Maxwellian
study, new oblique case or growth-rate search is introduced.

### Static tensor and physical growing roots

At omega=i*1e-8, the measured real raw determinants are:

| Grid | Corrected-production coherent D | Archived accurate-Bessel D |
|---|---:|---:|
| 32×64 | -1.519965021994278e-16 | -1.51996502199e-16 |
| 64×128 | -1.542586658124704e-16 | -1.54258665812e-16 |
| 96×192 | -1.542051306832260e-16 | -1.54205130683e-16 |

The residual imaginary parts are approximately `-3.06e-23`, `9.75e-23`,
and `1.18e-24`, respectively. The analytic reference is about `-1.54618e-16`.
The current and reconstructed historical six-component tensors and wave
matrices agree exactly at these probes and all 257 original-box frequencies.
Lower-box determinants are separately compared with the archived C2e curves.
The old positive `O(1e-10)` determinant floor is absent. The Bessel-defect
moment prediction is now at rounding scale; at 64×128 its z contraction
prediction is `-7.73e-19`, versus the historical `+8.48843e-10`.
The observer's retained `original/intrinsic` column labels now mean **shared
production call / direct intrinsic call** in current runs: they are equal,
not two evaluations of the historical inaccurate function.

Independent ordinary growing roots at positive kpar are:

| Grid | Positive real frequency | Growth | Complex error from analytic |
|---|---:|---:|---:|
| 32×64 | .38916307329839805 | .08561150418134517 | .00750915048 |
| 64×128 | .39477930797795668 | .08588395280663391 | .00188633548 |
| 96×192 | .39582524753381620 | .08593998388469128 | .000838896968 |

The independently solved negative partners reproduce the C2e reflection
symmetry. Maximum raw determinant residual is below `3.8e-14`; roots and
their residuals are retained at full precision. Agreement with the archived
roots is value-for-value in the optimized run, not merely within the existing
`1e-10` root comparison tolerance. The decreasing grid error is retained;
matching historical results does not remove finite-grid representation error.

### Matched trajectories and standard contours

For both sets of 257 ordered frequencies on each grid, current optimized
determinants match the archived C2e accurate-Bessel curves exactly. All three
original-box comparisons against the independently reconstructed historical
build also have zero determinant and full-tensor difference. Thus the
comparison does not rely solely on archived rounded root coordinates.

| Grid | Median relative D error | Maximum relative, original | Maximum relative, lower | Maximum absolute, original |
|---|---:|---:|---:|---:|
| 32×64 | .0127795 | .274638 | .528795 | .0940136 |
| 64×128 | .00322331 | .0675131 | .129543 | .0235663 |
| 96×192 | .00143553 | .0298145 | .0586804 | .0104674 |

These errors compare with the unchanged analytic determinants at identical
frequencies. They are not differences between unrelated adaptive grids.
The maximum original-contour cumulative phase errors are `.215799`,
`.0580640`, and `.0259886` radians; total-loop phase differences are at
rounding scale. Absolute errors and minimum reference/candidate magnitudes
are retained in `errors.dat`, including the lower contour's small analytic D.

Independent A1 evaluations use unchanged defaults and give:

| Grid | Gamma bounds | Accepted W | Calls on each rank | Passes |
|---|---|---:|---:|---:|
| 32×64 | [.01,.2] | 2 | 616 | 3 |
| 32×64 | [.001,.2] | 2 | 702 | 3 |
| 64×128 | [.01,.2] | 2 | 630 | 3 |
| 64×128 | [.001,.2] | 2 | 690 | 3 |
| 96×192 | [.01,.2] | 2 | 630 | 3 |
| 96×192 | [.001,.2] | 2 | 686 | 3 |

All are status zero with total phase `4*pi`, matching the archived counts
and work exactly. The isotropic numerical-proton control has accepted winding
zero on both prescribed rectangles at all three grids. Neither count is a
global stability classification outside the specified frequency box.

### Analyticity diagnostics

The unchanged full-tensor mathematical fixture gives static longitudinal norm
`1.00614e-16` and maximum nodal identity error `5.11788e-16`.
Its independent zz CR residual decreases from `8.56e-7` at step `1e-4`
to `8.56e-9` at `1e-5` and `4.04e-11` at `1e-6`.
Actual physical tensor and determinant CR groups cover internal nodes,
half cells, support endpoints and exterior points, at gamma .001 and .01.
All existing determinant decrease tests pass in the optimized three-grid
matrix. For example at omega=.4+.001i on 64×128, the full proton-tensor
residual decreases from order `1e-6` to order `1e-10` before increasing at
the smallest perturbation from roundoff. Near-origin and growing-root
groups are retained separately rather than imposing monotonic improvement
into the floating-point noise regime.

At the physical growing-root coordinate on 64×128, the proton-tensor CR
residuals are `6.5233e-4`, `6.5252e-6`, `6.5253e-8` for successive tenfold
step reductions. At omega=i*1e-5 the determinant residuals are
`4.2347e-7`, `4.2227e-9`, `9.0012e-10` for steps `1e-6,1e-7,1e-8`.
The tensor difference at this near-origin point is already roundoff dominated
at the smallest steps. No sampling or derivative tolerance was relaxed.

## Low-frequency exclusion and conditioning

Each window below is tested independently with raw D and diagnostic
D/omega², on **all three** momentum grids. In the optimized matrix every
one of the 30 tests accepts winding zero, status zero and three passes. The
same complete matrix also passes in debug, for 60 accepted exclusions:

| Real bounds | Gamma bounds | Raw W | D/omega² W |
|---|---|---:|---:|
| ±5e-5 | [5e-6,2e-5] | 0 | 0 |
| ±5e-5 | [1e-6,1e-5] | 0 | 0 |
| ±5e-5 | [2e-5,5e-5] | 0 | 0 |
| ±2e-4 | [1e-7,1e-3] | 0 | 0 |
| ±2e-4 | [1e-8,1e-3] | 0 | 0 |

These reproduce the archived C2e accurate-Bessel results and exclude the old
`~1.38e-5` artifact in the tested resolved windows. The origin is outside
each rectangle, so dividing by omega² has zero net closed-contour winding.
It is a validation-only callback, not a production determinant change.
The broad scaled control also accepts zero with initial intervals doubled
to 16 (896 calls), and with the diagnostic +8 harmonic limit. The enlarged
limit gives zero difference on all 257 matched original-contour determinants.

The existing C2e real128 postprocessing acts on the **actual saved wave
matrix**. At 64×128 and gamma=1e-8, D/omega² is about `1.54259` and the
scaled reciprocal condition measure about `.0154`; the estimated input-ulp
perturbation ratio is about `.047`. At gamma=1e-9 that ratio exceeds one.
Nearer-origin raw determinant smallness is therefore not accepted as a root
or a zero-exclusion certificate. The original uncertainty diagnostics and
probes down to 1e-10 are retained with this explicit limitation. No ordinary
solver "zero" inferred solely from an absolute raw-D threshold is promoted
to a physical or resolved numerical mode.

## Weight reuse and default-production separation

The unchanged cache tests pass the original `2e-13` absolute weight/kernel
tolerance, with measured maximum error `1.77636e-15`. They independently
invalidate on frequency, species, kpar and parallel-grid changes, and show
that repeating a species/harmonic set does not rebuild it. The application is
`sum(weights*G)`, without complex conjugation. No determinant interpolation,
approximate-frequency reuse or new MPI request is introduced.

On the optimized matched 64×128 trajectory, cached versus uncached determinant
maximum relative difference is `1.23204e-14`, passing the unchanged C2d
comparison tolerance. The n=0 owner records:

| Implementation | Kernel-cell calls | Derivative-node reads | Weight sets | Cache bytes |
|---|---:|---:|---:|---:|
| Cached | 971460 | 357788718 | 3855 | 64132 |
| Uncached | 177485742 | 357788718 | 0 | 0 |

The cached work exactly reproduces C2d/C2e. All four ranks still participate
in every one of the 257 determinant calls. The weights are local to the
unchanged harmonic assignment; no frequency or wavevector worker groups exist.
Primary nmax remains 15 and current/historical stored array/halo records agree.

Ordinary default-production controls, built without a coherent dispatch,
recover the corrected C2f.1 positive roots:

| Control | Real frequency | Gamma |
|---|---:|---:|
| Numerical-proton growing | .39541172517172751 | .08622020920473941 |
| Isotropic numerical-proton damped | .40533729651880795 | -.00170948152749903 |

Their signed partners are independently solved and residuals meet `1e-13`.
These values differ from the coherent growing roots as expected, confirming
that the experiment is not the production default. The damped control uses
the unchanged C2b integrator and fitted Landau continuation, not a new C2d
causal continuation. The existing mismatch between grid and fitted boundary
residues remains outside this rebaseline's scope.

## Representative performance

Three sequential four-rank 257-frequency repeats disable diagnostic I/O and
use identical inputs and frequencies. Rank-0 collective determinant time is:

| Implementation | Seconds per trajectory, range | Mean seconds per call |
|---|---:|---:|
| Cached coherent | 11.1007–11.4689 | .04390585 |
| Uncached coherent | 75.3694–75.7411 | .29417507 |

The measured ratio is approximately **6.70×**, consistent with preservation
of the historical acceleration. This is a within-build comparison, not a
claim that the rebaseline optimizes C2d. Weight construction totals about
`.361–.363 s` on the n=0 worker per cached trajectory. The one-time Bessel
constructor takes roughly `28.6–28.9 microseconds` on that worker, and cutoff
determination about `80.7–81.6 microseconds`. There is one constructor and
one cutoff determination, not one per frequency. Full per-rank work and call
counts are retained; timing does not imply performance for other VDFs.

## Reproduction and retained artifacts

Run from a checkout containing the preserved C2e/C2f.1 reference fixtures:

```sh
bash validation/run_nyquist_c2f2.sh
```

The script creates a fresh isolated scratch study by default; optionally set
`C2F2_SCRATCH` to a fresh directory outside the repository. It requires no
previous `/tmp` data. Distributions, numerical derivatives and fitted startup
representations are prepared through normal ALPS paths. Every MPI case uses
four ranks, with BLAS/OpenMP threads set to one, and jobs execute sequentially.
The measured study is `/tmp/alps-c2f2-current/`.

GNU Fortran 11.4.0 and Open MPI 4.1.2 are the identified toolchain, with the
existing BLAS/LAPACK libraries. Optimized flags are
`-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags are
`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
Mathematical helpers additionally retain their existing full checks.
`provenance/toolchain.txt` records actual linked libraries and versions;
common-build manifests and object/module/executable hashes identify linkage.

New files are the shell entry point, `validation/c2f2/{build,physical,
provenance,collect,bessel}.sh`, the two Fortran utilities
`compare_tensors.f90` and `bessel_identities.f90`, their README, this report,
and compact measured data under `tests/data/nyquist_c2f2/`. No old builder,
scientific fixture, production file or documentation of a historical result
is rewritten. The corrected-corner mathematical expectation is generated only
in a copied fixture and checked at its original tolerance.

Full raw harmonic/tensor streams, generated VDF arrays, ordinary root
observations, inputs and source variants remain in scratch. The compact
export includes matched curves and analytic error/phase statistics, all
adaptive statuses/work, per-rank calls, static tensors/contractions,
conditioning and CR diagnostics, Bessel identities, root data and mathematical
logs. A representative complete historical/current six-component tensor
comparison is retained at 64×128. Its `gate.txt` is written only after
checking archived roots, exact standard-contour work, accepted low-frequency
counts, isotropic controls, production controls and cache-work signatures.
Fixed-frequency PROBE's deliberate shutdown status is **not** an accepted
adaptive winding.

## Remaining limitations

This task revalidates an existing experimental upper-half-plane model. It
does not make C2d the production default, change C2b's Tlim treatment, validate
causal damped continuation, or extend numerical-species scans/growth searches.
The default production control retains C2f.1's corrected real Bessel and
corner implementations and the existing fitted Landau contribution.

Finite-grid determinant and root errors remain measurable and decrease under
the tested refinement; exact reproduction of C2e is not a new proof of
physical accuracy for arbitrary VDFs. Near-origin conditioning below roughly
1e-8 remains unresolved in double precision. Finite contours do not establish
global stability outside their boundaries.

The three inherited harness issues recorded in `ALPS_BESSEL_INTEE_RESULTS.md`
(debug manufactured metadata, historical normalization, archived mixed-map
comparison) are separate from this rebaseline. They are not repaired or
relabeled here, and the full C2f.1 relativistic matrix is not redundantly
rerun when no production source changes. Existing analytic Nyquist
compatibility suites provide the additional validation-construction gate.

During harness development, editing the entry script while its first matrix
was executing caused Bash to resume at an obsolete file offset after the
physical block (`al: command not found`). All physical jobs had completed.
The remaining unchanged compatibility/collection stages were resumed from an
immutable tail script. The finalized entry point is syntax-checked and every
stage is executed, but a second full fresh replay of the final convenience
script is not claimed. This operational interruption changed no scientific
code, expected value, tolerance or completed result.

## Final validation disposition

| Gate | Optimized | Debug |
|---|---|---|
| Independent linked production Bessel API, 616 rows | PASS | PASS |
| Loaded-node Bessel identities, all three grids | PASS | PASS |
| C2c primitives and direct signed-b tests | PASS | PASS |
| C2d moments, PWL convergence, static/tensor identities, cache invalidation | PASS | PASS |
| Three-grid growing partners and archived root residuals | PASS | PASS |
| Both 257-point matched contours, all three grids | Exact archived agreement | Exact archived agreement |
| Reconstructed historical/current full tensor and wave comparisons | Exact agreement | Exact agreement |
| Six standard unstable contours | W=2; archived calls | W=2; archived calls |
| Six isotropic contours | W=0 | W=0 |
| Five low-frequency windows × raw/scaled × three grids | All 30 W=0 | All 30 W=0 |
| Enlarged-harmonic and stricter-sampling controls | PASS | PASS |
| Cached/uncached equivalence | 1.23204e-14 max relative | 1.23204e-14 max relative |
| Default C2b growing/damped compatibility | PASS | PASS |
| Existing analytic A1/A2, growth, log, mirror, line, plane suites | All six PASS | All six PASS |

All 88 independent adaptive runs are status zero with three convergence
passes: 12 have winding two, 76 have winding zero. Standard unstable calls
remain exactly 616/702, 630/690 and 630/686 by grid and rectangle in **both**
builds. Existing analytic suites retain their legacy scientific-output
comparisons, aligned scan outputs and MPI fault-recovery checks.

The 146 numerical validation jobs record 55,932 collective determinant calls
in total, with identical physical-call counts on all four ranks in every
job. This total includes diagnostic and timing jobs; it excludes the two
fresh analytic-reference trajectories and the separate analytic compatibility
suites. No new evaluator failure, numerical mismatch, harmonic-allocation
change or MPI synchronization failure was observed. Runtime array/halo
comparisons retain primary nmax=15 and identical historical/current assignments.
All protected production, experimental and historical-reference hashes pass.

The collector initially selected the 32×64 reference root for every grid;
its missing grid filter was corrected and all reference checks rerun, without
changing tolerances or reference values. The final export also explicitly
separates the six fixed-frequency lower-trajectory diagnostics from accepted
adaptive runs. These were validation bookkeeping errors, not differences in
the susceptibility. Their development logs remain in scratch.

**Decision: successful rebaseline.** The production Bessel implementation
supplies the accurate coefficients directly, the coherent E moment is included
exactly once, archived upper-half-plane results are reproduced, and reusable
weights remain active. This provides the requested baseline for separately
scoped broader physical validation. It authorizes neither production-default
adoption of C2d nor causal damped continuation. No numerical-electron,
both-species numerical, non-Maxwellian or new oblique Nyquist study was begun.
