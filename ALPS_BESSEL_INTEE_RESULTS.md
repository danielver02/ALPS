# C2f.1 — Real Bessel and nonresonant corner corrections

## Scope and review boundaries

This work treats two independent defects. `validation/bessel_intee/int_ee.patch`
is a one-factor correction to the nonrelativistic quadrature. The separate
`bessel.patch` contains the accepted shared real-Bessel replacement. Neither patch changes
the continuum susceptibility formulas, resonance integration, continuation,
MPI decomposition, harmonic stopping criterion, root solvers or Nyquist methods.

The exact pre-change `ALPS_fns.f90` and `ALPS_fns_rel.f90` are retained in
`tests/data/bessel_intee/before/`, with the prior working-tree patch and complete
source/executable provenance in the isolated study `/tmp/alps-c2f1/`.
The C2b changes and unrelated pre-existing working-tree modifications were not
reverted. C2d/C2e mathematical implementations remain unchanged.

Both corrections are now applied to production. The shared replacement was
validated in isolated source trees before adoption; the corner-only fix was
independently accepted first. Original failed harness gates remain explicitly
identified below rather than being relabeled as successful regressions.

## 1. Shared real-Bessel call-site audit

The original public real functions are `BESSJ`, `BESSJ0`, and `BESSJ1` in
`ALPS_fns_rel.f90`. There are two production source files with their callers:

| Caller | Argument and orders | Use and preserved behavior |
|---|---|---|
| `ALPS_fns:determine_nmax` | `kperp*pp(s,iperp,1,1)/qs(s)`, positive integer orders | Maximum over the original perpendicular nodes; unchanged first-below-`Bessel_zero` stopping rule and process-count padding |
| `ALPS_fns:determine_bessel_array` | Same signed real argument; local `nlim(1)-1:nlim(2)+1` | Same allocation, complete halo values, and special negative-one halo handling |
| `ALPS_fns:int_T` | Stored values and neighboring orders | Same derivative recurrence, tensor factors, and exact parallel special cases; no new frequency-loop Bessel calculation |
| `ALPS_fns_rel:int_T_rel` | `kperp*ms(s)/(vA*qs(s))*sqrt(Gamma**2-1-pbar_parallel**2)` | Signed charge argument, positive-order calls with explicit negative-harmonic parity, same derivative/tensor formulas and relativistic quadrature geometry |
| `BESSJ0`, `BESSJ1` | Public real scalar functions | Historically also seeds for BESSJ; public names/types retained |

Relativistic real tensor evaluation already computes these values at its real
quadrature nodes; it does not universally use the nonrelativistic stored array.
The candidate changes the accuracy of those existing calls, not their location
or scheduling. In particular it adds no frequency-dependent constructor.

`int_T_res_rel` uses the **separate complex** `CBESSJ` path, with its existing
complex argument and parity conventions. CBESSJ and every complex-continuation
formula are unchanged. Analytic bi-Maxwellian response uses its independent
NHDS implementation; it is not routed through the replaced real routines.

## 2. Independent pre-change reproduction

The original functions were executed before editing production code. The C2e
real128 reference and the new broader signed-argument audit reproduce:

* Small nonzero J0 tends to `57568490574/57568490411`, with intercept error
  `2.831409995351919e-9`; the small-argument squared-sum defect is approximately
  `5.66282e-9`, while the exact-zero special return is one.
* J1's small-argument slope error is `7.600609031044314e-11`. Its default-real
  asymptotic phase literal gives a promoted error of
  `4.827995985223765e-7` relative to `3*pi/4`. Truncated J0 phase/amplitude and
  rational coefficients are additional approximation errors.
* The `X>N` recurrence selector sends every negative argument to Miller's
  downward recurrence. Its finite starting order can be inadequate there:
  `BESSJ(7,-32)` is about `137.245`, versus the independent `0.0762103487`.
* The old single `1e-10` rescaling is not a general intermediate-overflow
  safeguard. An independent floating-point-trap run at high order and tiny
  argument overflows inside the recurrence, before normalization. The broad
  audit records explicitly skipped unsafe legacy evaluations as NaN rather
  than treating them as successful results. The intrinsic is tested there too.

`legacy_bessel.f90` preserves the original three functions with renamed
procedures for direct comparison. Its mathematics, coefficients, starting
order and scaling are unchanged. Existing C2e reference tables are retained
separately; expected values are not inferred from the corrected production code.

## 3. Replacement strategies and accuracy

Option 1 uses Fortran 2008 `BESSEL_JN`, `BESSEL_J0`, and `BESSEL_J1`, behind
the existing ALPS interfaces. GNU Fortran 11.4.0 supports these interfaces in
the tested environment. No library or input parameter is added.

Option 2 was substantively tested: accurate intrinsic low-order seeds,
`abs(x)` recurrence selection, explicit argument parity, a larger Miller
starting order (`IACC=200`), and stronger `1e100`/`1e-100` scaling. This avoids
preserving the inaccurate rational approximations merely for compatibility.
It still requires intrinsic seeds, so it offers no compiler-portability
advantage over Option 1 and maintains more numerical code.

The broader audit uses independent real128 series or periodic Fourier
integration, with the latter repeated at double node count. It covers signed
arguments through 1024, orders below/near/above the argument through 1088,
actual proton/electron argument scales, original regime switches near eight,
zero and arguments down to `1e-300`. There are 616 reference rows, including
points near genuine J0/J1 zeros. Parity of
both signed argument and signed order, independent derivatives, order recurrence
and public J0/J1 consistency are checked. Relative errors are reported only
when the reference magnitude exceeds `1e-12`; absolute errors remain primary
near zeros or underflow.

| Method | Maximum absolute error | Maximum relative error on the stated subset |
|---|---:|---:|
| Original (finite, safe calls) | 1.37168858931e2 | 2.44378513088e3 |
| Intrinsic | 4.65204229623e-16 | 1.04826605765e-14 |
| Repaired recurrence | 5.80100851196e-16 | 8.74561370486e-15 |

Optimized and debug audits give the same maxima. A separate executable linked
against the **actual candidate ALPS module**, rather than a reimplementation,
checks all reference rows through BESSJ and the public low-order interfaces.
It passes in both builds (maximum absolute error against rounded reference
values `4.64906e-16`). No Bessel arrays are renormalized.

For the same optimized scalar timing workload, the original, intrinsic and
repaired-recurrence methods take `.1750`, `.3009` and `.2264` seconds
respectively. The repaired recurrence is faster than the intrinsic in this
bounded workload but offers no accuracy or portability advantage and retains
more algorithmic complexity. Physical costs are measured separately below;
the intrinsic was selected for numerical correctness and simplicity rather
than scalar benchmark speed.

C2e's retained independent identity tables separate value error, omitted
harmonic tail and summation roundoff. The derivative identity is checked against
an independently differentiated reference; merely constructing a derivative
from neighboring values would not validate it. In particular the finite
`sum(Jn*Jn')` equals `-J_N*J_(N+1)`, not zero. The symmetric weighted sum
`sum(n*Jn**2)` cancels algebraically under parity and cannot by itself establish
value accuracy. Exact-zero tensor limits avoid evaluating `n*Jn/x` at zero.
In the actual int_T test extended to `kperp=1e-12`, the original values leave
a `5.66282e-9` error floor, while corrected values give the expected `1.4e-12`
approach-to-limit error. At exactly zero both retain the existing special
limit within `5.6e-17`. No tensor limiting formula is changed.

The intrinsic is the preferred candidate: it is simple, independently accurate
in the tested domain and removes both low-order approximation defects and
incorrect signed recurrence selection. This is a bounded compiler/domain
validation, not a universal guarantee for arbitrary orders or arguments.

## 4. Harmonic stopping and MPI protection

The stopping algorithm, `Bessel_zero`, process padding, `split_processes`,
communicators and collective determinant assembly are not edited.
The signed-support replay includes zero, reference proton and electron scales,
relativistic scales, and moderate/large arguments. All 20 tested signed support
cases retain the original cutoff, including 15 at the primary reference,
43 at `kperp=.3`, 108 for support argument 32, and 245 for 128.
The independent C2e tests already establish that the primary omitted tail is
negligible compared with the old value-approximation error.
An independent real128 series additionally checks the cutoff and predecessor
coefficients. At support argument 128, `|J244|=1.063422722107e-45` and
`|J245|=2.997597420562e-46` qualify the unchanged cutoff; intrinsic relative
error in J245 is `1.35e-15`. The tail after 245 is approximately `1.53388e-92`;
a series/geometric bound beyond another 128 orders is `6.27e-243`.

Ordinary-run observers record actual species limits, rank assignments,
allocated bounds and all stored local/halo orders at selected perpendicular
nodes. They introduce no MPI operations. Physical call-count checks compare
all ranks, including ranks without a numerical contribution.
Observed changes in array **values** are expected; allocation/assignment
changes require independent explanation, not an adjusted stopping criterion.

## 5. Nonrelativistic int_ee correction

Only the first usable perpendicular/parallel corner is changed:

```fortran
! original
pp(sproc,1,1,1)*df0(sproc,1,1,1)
! corrected
pp(sproc,1,1,2)*df0(sproc,1,1,1)
```

The continuum numerator is `p_parallel*(p_perp*df_parallel-p_parallel*df_perp)`.
All other corners already use the second momentum factor. The special corner
weight and final trapezoid normalization yield the original-minus-corrected
integral difference

```
delta_E = pi*q/m * delta_p_perp*delta_p_parallel
                  * p_parallel*(p_parallel-p_perp)*df_perp.
```

No endpoint weight, domain, resonance formula, derivative, relativistic integral
or continuation term changes. The preserved original int_ee is directly
callable only in validation. Constant, polynomial and Gaussian derivative
fixtures, all three grid shapes and growing/real/damped arguments pass in
optimized and debug builds.

The independent C2e physical predictions retained for comparison are:

| Proton grid | Original-minus-corrected integral |
|---|---:|
| 32×64 | -5.05089351866e-9 |
| 64×128 | -4.00403595412e-10 |
| 96×192 | -1.01397106636e-10 |

In the assembled tensor `delta_chi_zz=-ns*qs*delta_E`. The corresponding
determinant change is
`delta_D=delta_chi_zz*(Wxx*Wyy+Wxy**2)`. Matched-frequency checks retain the
density and signed charge explicitly, including numerical-electron controls.
Unrelated tensor entries and the corner-only Landau streams must remain
unchanged. The coherent C2d moment already has the intended expression and
receives no additional correction.

## 6. Physical and compatibility validation

Four isolated configurations start from the exact preserved source snapshot:
I original, II Bessel-only, III corner-only, IV combined. Ordinary existing
root paths are used for analytic and numerical growing/damped partners, the
three proton resolutions, numerical electrons, oblique propagation and the
exact/approaching parallel limit. Existing matched-frequency diagnostics
supplement those ordinary runs; no new Nyquist physics is introduced.

The pre-change complete newer relativistic susceptibility matrix passed:
three equilibrium kinds, seven resolution combinations each and corresponding
one-versus-four-rank controls. The same candidate optimized matrix also passed.
Across the 21 four-rank comparisons the scientific changes are quantified
against the frozen same-version baseline, not a different historical source.

The historical `run_numerical_baselines.sh --convergence` command failed
**before physics** in the unchanged tree: `instrument_nr_history.py` assumes a
single-line `legacy_zero_return` assignment, which C2b has split across lines.
This failure is retained as a failed gate; it is not called a full pass and its
scientific fixtures/tolerances are not rewritten. The repository's validation
README additionally documents an inherited rel4 normalization failure in that
historical runner. The available newer relativistic suites are run separately.

One original debug manufactured-row test also stops on an unallocated
`poly_kind` array in its fixture, reached by Fortran's non-short-circuit
evaluation of the shared eligibility expression. Normal input initializes
this metadata. A copied-fixture adapter initializes it for all two species;
it changes no production source, integrand, expected value or tolerance.
Original and initialized-fixture results are retained separately.

Across the optimized matched relativistic matrix the largest non-extrapolated
raw-determinant change is approximately `1.008e-4`, with largest relative
change `1.513e-8`. Some pre-existing Chebyshev extrapolation diagnostics have
chi of order `1e138` and raw determinant `1e263` in **both** versions. Their
large absolute determinant changes (up to `1.38e253`) accompany small relative
changes: an unchanged huge fitted Landau term amplifies the small corrected
direct-integral contribution in the determinant. These finite extrapolated
values were already inaccurate and are not newly validated physical responses.
All failure-status fields match. Full rows and non-extrapolated summaries are
retained rather than dropping this important limitation.

## 7. Nonrelativistic physical results and independent attribution

The complete optimized and debug four-way ordinary matrices passed the matched
corner/cofactor, unchanged-other-components, allocation and signed-frequency
symmetry checks. All eight optimized 64×128 growing/damped root files reproduce
the **preserved C2e four-way files byte for byte**. These are independent
pre-existing references, not recalculated expected values from the new code.

The four-way 64×128 positive-frequency roots are:

| Configuration | Growing real frequency | Growth | Isotropic real frequency | Damping gamma |
|---|---:|---:|---:|---:|
| I | .395411724579199 | .086220209553379 | .405337298402431 | -.001709481504305 |
| II | .395411725171727 | .086220209204739 | .405337296518827 | -.001709481527499 |
| III | .395411724579199 | .086220209553379 | .405337298402473 | -.001709481504305 |
| IV | .395411725171728 | .086220209204739 | .405337296518808 | -.001709481527499 |

Raw determinant residuals are below the existing `1e-13` threshold. Negative
real-frequency partners were independently solved and satisfy drift-free
reflection within `1e-10`, usually much better. The corner-only root shift is
small because these nearly parallel modes couple weakly to the longitudinal
entry; the full tensor and determinant cofactor checks still resolve its effect.

The corrected default C2b backend retains convergence toward the independent
analytic roots `.39666286088749+.08598637345950i` and
`.404601381588-.001854305370i`:

| Grid | Growing real frequency | Growth | Isotropic real frequency | Damping gamma |
|---|---:|---:|---:|---:|
| 32×64 | .390581513785533 | .086311784441269 | .407435769576102 | -.001312181382559 |
| 64×128 | .395411725171728 | .086220209204739 | .405337296518808 | -.001709481527499 |
| 96×192 | .396168901299062 | .086138999470945 | .404931365196624 | -.001772382115442 |

The residual grid/resonance/continuation errors are much larger than the Bessel
correction. In particular the 64-grid damping error remains appreciable; this
task does not repair the existing finite-window Tlim approximation or adopt
the C2d formulation. Integer winding agreement is not used as its acceptance
criterion, and no new numerical-VDF Nyquist analysis is introduced.

At the oblique fixed wavevector `(.3,.3)`, the combined positive root is
`.428085353535213+.028336683961839i`, compared with
`.428085353348337+.028336683827522i` before either correction. Its independently
solved negative partner and ordinary residuals remain consistent. The harmonic
limit is 43 and all four variants retain identical assignments and halo bounds.

At `kperp=0` the Bessel-only change leaves the ordinary root unchanged. At
`kperp=1e-8`, the original coefficient bias produces a complex root offset
of approximately `6.88e-10` from the exact parallel limit; corrected values
agree with that limit within roughly `2e-16`. The finite-`kperp=1e-6` change
then has the expected small physical scale, not the old nonzero limiting floor.

Numerical electrons use the normal mass-scaled generated grid, analytic
protons and unchanged continuation. Their corner prediction is
`delta_E=+1.460948876232422e-6`, hence the same positive `delta_chi_zz` for
unit density and negative charge. Matched growing/real/damped frequencies
confirm the cofactor relationship with maximum absolute discrepancy below
`6.3e-17`. The corresponding proton checks are below `3e-17`.
Corner-only Landau streams remain byte-identical on every rank, including the
stable, electron and oblique controls. Bessel-only changes can legitimately
change the real Bessel factors in the nonrelativistic Landau contribution;
its formula, fitted model and causal multiplier are unchanged.

The original electron case uses the existing default solver and default
`D_threshold=1e-5`, so its observed `|D|~2.6e-8` is valid under those controls
but insufficient for a high-precision correction comparison. A separate,
explicit stricter run uses the same existing controls as the primary proton
fixture (`secant_method=0`, `D_threshold=1e-13`) for **all** four variants.
Default-control data are retained; no production tolerance is changed.

The stricter electron positive root changes from
`.396817117215722+.085923615643942i` (I) to
`.396817116718169+.085923615909987i` (IV); all four variants and both builds
meet `|D|<1e-13`. This separates the small correction from the looser default
solver's stopping error.

## 8. Completed relativistic and legacy compatibility gates

The full newer relativistic susceptibility, controlled-root, Chebyshev MPI,
parallel susceptibility and fitted-continuation suites completed in both
optimized and debug configurations. The optimized analytic-parallel suite
also passes, retaining byte-identical values for its 864 supported Landau
records. The original debug analytic-parallel fixture fails in both frozen
and corrected sources because it omits allocation of `poly_kind`. The copied
fixture's metadata-only initialization makes both pass; the repository fixture
and its expected scientific values are unchanged.

All 72 controlled relativistic roots in each build retain their existing
scaled-determinant qualification. Across these 144 before/after comparisons,
the largest complex-frequency shift is `1.176803099412033e-10`; the largest
corrected scaled residual is `1.154e-16`. The full ordinary main-program path
also recovers the isotropic transverse root
`1.098159581309793-2.299832948659579e-7i` after correction, versus
`1.098159581424544-2.299832949812484e-7i` before correction.

For the anisotropic relativistic fixture, a generic initial guess near 1
selects a longitudinal near-real solution in **both** versions. It fails the
existing transverse polarization qualification and is retained as a negative
control. Seeding from the independently qualified pre-change controlled root
recovers the proper ordinary transverse branch: the positive root changes
from `1.117122761309876-2.539412191344676e-7i` to
`1.117122761192209-2.539412191838508e-7i`. Both signed partners and both builds
meet the existing residual/polarization criteria. The small signed-frequency
asymmetry of these relativistic fits was already present before correction
(roughly `1.08e-9` in real frequency and `1.66e-10` in damping for the
anisotropic pair); this task does not claim to repair it.

Three inherited harness limitations remain, with original statuses preserved:

* The debug manufactured fixture described above exits 2 before its copied
  metadata initialization; the corrected fixture passes for both sources.
* `run_rel_production.py` exits 1 at its archived mixed-nonrelativistic map
  byte comparison for both sources. The preceding 19 activation/error cases
  agree exactly, as do the pure method-2/mixed-analytic archive comparisons.
  Before/after activation CSVs and all 28 exported model binaries are
  byte-identical. The four-way ordinary map/cofactor checks supplement this
  aborted archive gate without rewriting its fixtures.
* The historical numerical-baseline runner first fails on its stale source
  instrumentation anchor. With only that locator adapted in copied sources,
  both versions reach and fail the same inherited relativistic normalization
  check (`rel4`, measured normalization approximately `.97476`). No physics
  tolerance or normalization was changed to make it pass.

All twelve optimized/debug analytic Nyquist suite runs (A1/A2, growth,
logarithmic, mirror, line and plane) pass, including their established legacy
scientific-output byte comparisons. No assertion is made that every original
repository regression command exits successfully: the original failures above
remain open and are distinct from the completed numerical compatibility gates.

## 9. Measured computational cost

Three sequential four-rank repeats use identical frequencies with diagnostic
I/O disabled. Times are rank-0 collective determinant time, excluding startup.

| Backend and workload | Original mean seconds/call | Corrected mean seconds/call | Change |
|---|---:|---:|---:|
| Numerical proton, 257 matched frequencies | .0605636 | .0609804 | +0.69% |
| Relativistic, ordinary 81-point upper-half-plane map | .0633429 | .0752525 | +18.8% |

The nonrelativistic repeat ranges overlap. The relativistic cost increase is
substantial and is retained as a tradeoff: its existing tensor quadrature
already calls real Bessel functions inside determinant evaluation. Replacing
those calls improves accuracy but is not solely an initialization cost.
No caching, call relocation or additional optimization is introduced here.

Mean harmonic-limit determination times are approximately `120.1/115.2`
microseconds (original/corrected) for the numerical-proton case and
`33.1/45.5` microseconds for the relativistic case. Mean local per-species
array-construction times across participating workers are `48.4/39.4` and
`12.1/11.5` microseconds respectively. These very short timings are noisy;
the compact tables retain every rank, repeat, allocation bound and actual
call count. Stored arrays are still constructed during initialization or
wavevector updates, never newly reconstructed per frequency.

## 10. Production disposition and reproduction

The accepted implementation consists of the three public real-Bessel
wrappers around the Fortran intrinsics `BESSEL_JN`, `BESSEL_J0`, `BESSEL_J1`,
and the independent single-factor `int_ee` correction. Both exact patches are
retained in `validation/bessel_intee/`. The installed source was checked
against the qualified isolated variants; the actual production modules were
built in optimized and debug modes. The production public API passes the
616-row independent reference test, and uninstrumented ordinary production
nonrelativistic/relativistic smoke outputs match the qualified instrumented
variants byte for byte. Protected source checks confirm no changes to other
production files or the C2c/C2d/C2e experimental implementations.

Reproduction entry point:

```sh
bash validation/run_bessel_intee.sh
```

It creates isolated scratch builds, generates distributions through the normal
ALPS generator, executes MPI jobs sequentially, retains failed original gates,
and exports compact evidence into its `compact/` directory. GNU Fortran 11.4.0
and Open MPI 4.1.2 were used with four ranks and one BLAS/OpenMP thread.
Optimized flags were `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`;
debug flags were `-O0 -g -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow`. Relativistic harnesses additionally use
`-ffree-line-length-none`. The wrapper block also compiles under strict
Fortran 2008 intrinsic checking.

The study was run component by component in `/tmp/alps-c2f1/`, including all
physical, mathematical, failure, comparison and final smoke stages. The
assembled convenience entry point was syntax-checked; a second fresh complete
end-to-end replay was not performed. Compact measured evidence, status files,
inputs and the exact two pre-change source files are retained under
`tests/data/bessel_intee/`; those are not regenerated expected results from
the corrected implementation. See the validation README for artifact paths
and the distinction between original failed gates and supplemental checks.

The corrections are ready for review independently. The remaining production
Tlim approximation, finite-grid/continuation errors, inherited harness defects
and relativistic performance penalty are explicit limitations. No broader
numerical-VDF Nyquist validation or production fixed-cell adoption was attempted.
