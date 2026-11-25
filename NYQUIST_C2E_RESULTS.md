# C2e — Bessel accuracy and the residual low-frequency zero

## Decision and scope

The C2d residual growing zero near `1.38e-5` is attributable to the original
Bessel values, rather than the C2d fixed-cell integration error or harmonic
truncation. Replacing **only the nonrelativistic numerical-species Bessel
array values** in copied validation builds removes that zero from the tested,
resolved upper-half-plane windows. The independently predicted static tensor
error agrees with the measured tensor change within `6.8e-19` absolute.
The physical cyclotron roots and their momentum-resolution convergence remain.

This is a qualified success for the targeted artifact, not a certification of
the entire upper half-plane: double-precision wave-matrix conditioning limits
the evidence below approximately `1e-8`. The apparent ordinary-solver roots
near `2e-9` fail the scale-aware checks and are not accepted as replacement roots.

The `int_ee` first-corner momentum factor is independently incorrect. Its
isolated correction has the predicted frequency-independent tensor effect.
The coherent C2d formulation already uses the intended moment and must **not**
receive an additional corner correction.

**No production source, default, namelist, shared BESSJ function, relativistic
routine, continuation fit, Landau algorithm, MPI decomposition, root solver,
or Nyquist algorithm was changed by C2e.** The pre-existing working-tree
changes were preserved. The source-hash manifest also protects the C2c Cauchy
kernel and C2d model, reusable-weight implementation and cached helpers.

## Reproduction and retained evidence

Run `bash validation/run_nyquist_c2e.sh`. It prepares the normal generated VDFs
and derivative inputs, builds isolated source variants, runs sequential
four-rank jobs, checks results, and leaves its scratch directory printed.
`C2E_SCRATCH` selects a fresh scratch directory;
`C2E_REFERENCE_STUDY` optionally reuses an identified C2c/C2d preparation.
No pre-existing `/tmp` directory is necessary for a fresh run.

The resumed measurements used `/tmp/alps-nyquist-c2e`, normal equilibrium data
from `/tmp/alps-nyquist-c2c/study_optimized`, GNU Fortran 11.4.0 and Open MPI
4.1.2. Optimized flags were `-O2 -g -fcheck=bounds -fbacktrace
-ffpe-trap=overflow`; debug flags were `-O0 -g -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow`. OpenBLAS/OMP threads were one. MPI jobs ran
sequentially with four ranks and the existing resonance assignments.

Compact evidence and runnable inputs are in `tests/data/nyquist_c2e/`.
The complete raw tensor, per-harmonic, constructor/halo and timing streams
remain in the isolated study directory. Fixed-frequency PROBE builds deliberately
use the existing validation shutdown; their status is not an adaptive A1
convergence result. Only independent ALPS/SCALED contour runs establish the
reported accepted windings. One exploratory exact-origin point was rejected
and its output retained separately; final four-way comparisons contain 34
finite, nonzero frequencies, including nonzero real frequencies at zero gamma.

New implementation files are confined to `validation/c2e/`: scalar Fortran
audits (`bessel_audit`, `tensor_limits`, `harmonic_limits`, `test_corner`),
tensor/result checks (`static_check`, `check_four`, `conditioning`, `compare`,
`contractions`), copied-source construction (`variants.awk`, `build*.sh`,
`constructor_audit.inc`, `stats.f90`, `scaled.sh`), and sequential execution /
collection scripts. `validation/run_nyquist_c2e.sh` is the public validation
entry point. No routine under `src/` was edited during this stage.

| Variant | Numerical response | Stored Bessel values | Old int_ee corner |
|---|---|---|---|
| coh_old | Preserved coherent C2d | Original | Coherent moment already correct |
| coh_acc | Identical coherent C2d | Intrinsic | Coherent moment already correct |
| coh_extra | Identical coherent C2d, +8 orders diagnostic | Intrinsic | Coherent moment already correct |
| I | Production C2b | Original | Original |
| II | Production C2b | Intrinsic | Original |
| III | Production C2b | Original | One-factor correction |
| IV | Production C2b | Intrinsic | One-factor correction |

## 1. Source audit and independent Bessel reference

The shared real Bessel implementation is `BESSJ` in
`src/ALPS_fns_rel.f90:1845`, with `BESSJ0:1916` and `BESSJ1:1962`.
The numerical-array constructor is `src/ALPS_fns.f90:4332`; harmonic stopping
uses BESSJ in `determine_nmax:4092`. `int_T:1707` uses stored values and
neighbor-order derivative recurrences. Real relativistic tensor construction
also calls BESSJ; complex continuation uses the separate CBESSJ path. This
shared scope is the reason no global replacement was made.

Three distinguishable defects were found:

* The small-argument J0 rational approximation has intercept
  `57568490574/57568490411`, exceeding one by `2.831409995351919e-9`.
  The special exact-zero return is one. At tiny nonzero arguments the squared
  normalization therefore has a defect close to `5.66282e-9`.
* The J1 small-argument slope differs from one half by
  `7.600609031044314e-11`. Its large-argument phase `2.35619491` is a
  default-real literal: after promotion its error relative to `3*pi/4` is
  `4.827995985223765e-7`. The J0 phase and amplitude are also truncated.
* The recurrence selector tests `X > FLOAT(N)`, not `ABS(X) > N`.
  Negative arguments always use the finite Miller recurrence, even when its
  chosen starting order is inadequate. This can be much worse than a few
  parts per billion: at `x=-32,n=7`, the original value is approximately
  `137.245`, whereas the independent value is `0.07621035`.

The audit does not extrapolate intrinsic accuracy without testing it.
`bessel_audit.f90` compares original and intrinsic BESSEL_JN against independent
real128 power series and periodic Fourier integration. The power series is
used through `|x|=16` and at sufficiently high order; otherwise a 1024-node
periodic Fourier integral is used. The reference derivative is differentiated
independently, rather than generated using the recurrence being tested.
See the [Bessel series](https://dlmf.nist.gov/10.2.E2),
[integer-order integral](https://dlmf.nist.gov/10.9.E2), and
[recurrences](https://dlmf.nist.gov/10.6.E1).

Coverage includes the actual loaded proton momentum nodes at all three
resolutions, `kperp=.001` and `.3`, negative arguments, stored halo orders,
arguments near zero and eight, and order/argument transitions through `|x|=63`
and orders through 128. Intrinsic maximum absolute value error was
`3.12858e-16` in both configurations. This is evidence for the tested domain,
not a universal compiler-library guarantee.

`bessel_identities_*.dat` distinguishes the finite sum, reference truncation
tail, function approximation and summation roundoff. Tails are summed directly
in real128, not inferred by subtracting a nearly-unit sum from one.
The symmetric finite `sum(n*Jn**2)` cancels algebraically when parity is imposed;
its small residual alone cannot establish Bessel accuracy. Likewise the finite
derivative sum is `-J_N*J_(N+1)`, rather than identically zero. Neighbor-value
derivatives telescope even for inaccurate values, so the independent derivative
reference is necessary. No values or sums are renormalized.

Representative identity results (the full argument/limit/method tables are retained):

| x | N | Original sum J² − 1 | Intrinsic sum J² − 1 | Independent omitted tail |
|---:|---:|---:|---:|---:|
| .0003 | 15 | +5.66282e-9 | +2.22e-16 | 1.97111e-149 |
| .01 | 15 | +5.66013e-9 | -1.11e-16 | 1.06373e-100 |
| 1 | 15 | -4.38137e-9 | 0 | 1.03378e-36 |
| 8 | 15 | -4.00853e-7 | -1.29677e-8 | 1.29677e-8 |
| 8 | 63 | -3.87886e-7 | -1.11e-16 | 8.81714e-102 |
| 16 | 63 | +3.24518e-7 | 0 | 6.83618e-64 |
| -32 | 63 | +3.77344e4 | -7.77e-16 | 5.25526e-28 |

The difference between the two x=8 rows illustrates why an insufficient
harmonic limit and a Bessel approximation defect must be measured separately.

The actual harmonic cutoff stays **15** at `.001` and **43** at `.3` for all
three grids with either stopping-rule evaluation. The experiment changes only
stored values, keeping original cutoff selection, indices, halos and workload.
The separate enlarged-cutoff diagnostic adds eight orders. At the primary
cutoff the independent tail is negligible compared with `5.66e-9`; the defect
is predominantly function approximation, not omitted harmonics. Actual
`int_T` special `kperp=0` limits and the approach from small nonzero arguments
were tested without changing the formulas for `Jn'` or `n*Jn/x`.

## 2. Static tensor attribution

Use the C2d notation `u=p_perp`, `p=p_parallel`, `A=df_perp`, `B=df_parallel`,
and the exact PWL moment
`E(u)=q/m * [u*M1(B)-M2(A)]`. C2d already combines the n=0 longitudinal terms
before integration. For the finite harmonic tensor, the remaining static
longitudinal contraction defect is predicted by

```
Rz = 2*pi*delta_u*kpar*ns*q * sum_i weight_i*E(u_i)
                                      * (1-sum_n Jn(x_i)**2).
```

This uses the actual prepared derivatives and unchanged perpendicular weights.
It does not fit a determinant offset. The constructor evaluates the finite
Bessel sum in real128 using the actual double values; `static_check.f90`
compares the prediction with `chi_zx*kperp+chi_zz*kpar` from the assembled
collective tensor at `omega=i*1e-8`.

| Grid | Predicted original-minus-accurate contraction | Measured difference | Absolute discrepancy |
|---|---:|---:|---:|
| 32×64 | 8.47160886909e-10 | 8.47160886301e-10 | 6.09e-19 |
| 64×128 | 8.48843348051e-10 | 8.48843348023e-10 | 2.85e-20 |
| 96×192 | 8.49151928799e-10 | 8.49151928127e-10 | 6.73e-19 |

For 64×128, the largest original finite sum defect over perpendicular nodes
is `5.66281931539e-9`; the intrinsic defect is `1.07704527194e-16`.
Its independently predicted residual contraction is `-7.73e-19`.
At `omega=i*1e-8` the raw determinants are:

| Grid | C2d original Bessel | Accurate Bessel |
|---|---:|---:|
| 32×64 | +2.91857100245e-10 | -1.51996502199e-16 |
| 64×128 | +2.94759360660e-10 | -1.54258665812e-16 |
| 96×192 | +2.95044412772e-10 | -1.54205130683e-16 |

The analytic reference is approximately `-1.54618e-16`. Enlarging the accurate
harmonic limit does not restore the original static defect or the extra zero.
All six tensor modes and individual harmonic streams are retained; the
agreement above establishes attribution more strongly than a smaller determinant.

## 3. Independent corner correction and four-way control

At `src/ALPS_fns.f90:1592`, the first usable perpendicular/parallel corner
uses `u*df_perp` instead of `p*df_perp` in
`p*(u*df_parallel-p*df_perp)`. The surrounding corners use the latter expression.
The factor two at this corner and the final trapezoid normalization give

```
delta_E = E_original-E_corrected
        = pi*q/m * delta_u*delta_p * p*(p-u)*df_perp.
```

Only that momentum factor changes in copied configurations III/IV. Original
and corrected `int_ee` are directly compared on prescribed constant,
polynomial and Gaussian derivative arrays for growing, real and damped
frequencies. The measured physical extra original contributions are:

| Grid | delta_E |
|---|---:|
| 32×64 | -5.05089351866e-9 |
| 64×128 | -4.00403595412e-10 |
| 96×192 | -1.01397106636e-10 |

This tail-corner error decreases rapidly under refinement; changing sampled
tail positions and the first perpendicular derivative prevents treating the
three numbers as a pure universal power law.

Configurations I/II/III/IV use C2b integration with original/accurate Bessel
and original/corrected corner respectively. At identical frequencies their
`delta chi_zz=-delta_E` and
`delta D=-delta_E*(Wxx*Wyy+Wxy**2)` agree with direct tensor/determinant
differences; the maximum determinant cofactor discrepancy is below `9e-18`.
Other tensor entries are unchanged by the corner correction. Landau streams
are byte-identical for I versus III and II versus IV on every rank, including
negative gamma. Bessel changes can change the factors multiplying a residue;
the residue algorithm and fitted representation remain unchanged.

The additional full 257-frequency four-way contour comparisons have maximum
absolute cofactor discrepancies below `4.7e-15` (larger-magnitude determinants
make subtraction roundoff larger than in the 34-frequency local probes).
All four C2b configurations give accepted winding two on both standard boxes,
but have substantial residual value errors:

| Configuration | Median relative, original box | Maximum relative, original box | Median relative, lower box | Maximum relative, lower box |
|---|---:|---:|---:|---:|
| I | .00430580415 | .137388653 | .00456227379 | 117.111107 |
| II | .00430580593 | .137386713 | .00456227580 | 117.111297 |
| III | .00430580419 | .137388928 | .00456227383 | 117.111080 |
| IV | .00430580597 | .137386988 | .00456227585 | 117.111270 |

The lower-box extreme is near the small analytic determinant (`1.53982e-6`),
with absolute discrepancy about `1.80330e-4`. Neither the small corner fix nor
the accurate Bessel constants cures the original two-branch integration's
low-frequency approximation error. This is separate from the C2d artifact;
matching an integer count would conceal it. Full absolute/phase statistics and
matched trajectories are retained in `controls/`.

The coherent C2d method uses the intended E moment in the upper half-plane,
so a four-way corner correction there would be mathematically misleading.
`coh_acc_corner` confirms that editing the unreachable old `int_ee` body has
no effect on its 257-frequency trajectory. No second correction is added.
C2b four-way contour counts are boundary-trajectory diagnostics, not a claim
that its remaining Tlim switches form a holomorphic dispersion function.

At 64×128 the C2b positive growing root changes from
`0.395411724579199+0.086220209553379i` (I) to
`0.395411725171728+0.086220209204739i` (IV). The corner-only change is at
root-solver rounding; the Bessel-only change accounts for the visible shift.
The isotropic damped roots are approximately
`0.405337298402431-0.001709481504305i` (I) and
`0.405337296518808-0.001709481527499i` (IV). These are controlled
finite-grid comparisons, not a new validation of the deferred C2d damped
continuation. No experimental coherent damped backend is introduced.

## 4. Residual zero: reproduction and independent exclusion

The original C2d candidate reproduces the three positive imaginary roots:

| Grid | Original extra gamma |
|---|---:|
| 32×64 | 1.38702602862192e-5 |
| 64×128 | 1.38359488691285e-5 |
| 96×192 | 1.38296116456261e-5 |

With accurate values, the ordinary absolute-threshold solver can return near
`2.2e-9`, with raw `|D|<1e-17`. This is false evidence of a resolved zero:
`D/omega**2` remains order one and the appropriately scaled wave matrix is
nonsingular. Additional guesses around and beyond the original zero do not
establish an independently resolved replacement.

For the independent checks the production determinant is unchanged. A copied
validation callback additionally permits `D/omega**2`. Since zero lies outside
each closed rectangle, this has exactly the same closed winding. Both raw
and scaled A1 calculations independently apply the original local refinement,
strict global passes and independently seeded convergence tests.

| Window: real bounds; imaginary bounds | Original Bessel W | Accurate W |
|---|---:|---:|
| ±5e-5; [5e-6,2e-5] | 1 | 0 |
| ±5e-5; [1e-6,1e-5] | 0 | 0 |
| ±5e-5; [2e-5,5e-5] | 0 | 0 |
| ±2e-4; [1e-7,1e-3] | 1 | 0 |
| ±2e-4; [1e-8,1e-3] | 1 | 0 |

All five pairs were tested raw and scaled at 64×128. The narrow and expanded
windows were repeated at all three grids; the narrow test also passed in
debug for both methods. The expanded accurate window remains zero after
doubling the initial sampling density (896 calls, three passes, no local
refinements), and after adding eight harmonics. Analytic reference windows
also return zero. The evidence rules out a moved root in these resolved
windows; it does not cover all other complex frequencies.

`conditioning.f90` independently validates its real128 determinant/adjugate
on diagonal and rank-deficient fixtures. It rotates the actual saved wave
matrix into transverse/longitudinal coordinates and divides its longitudinal
row and column by omega. It reports a Frobenius reciprocal condition measure,
`D/omega**2`, and a conservative perturbation bound for one-ulp-sized **input
matrix** changes. This is not a bound on all susceptibility errors.
At 64×128 and gamma `1e-8`, the normalized determinant is `1.5425867`, the
reciprocal measure is about `.0154`, and the input perturbation ratio is `.047`.
At `1e-9` that ratio exceeds one. At 96×192 and `1e-10` raw D can change sign
at rounding scale. Such points cannot certify a tiny new zero or prove its
absence. Higher precision or a more intrinsically conditioned wave-matrix
evaluation is required for that question; neither is silently introduced here.

The C2b four-way diagnostics I–IV all give zero on the narrow and expanded
low-frequency windows at 64×128. The C2d artifact is specifically exposed by its
restored static structure with inaccurate Bessel coefficients; absence in a
different integration method does not make the Bessel defect harmless.

## 5. Physical roots and standard contours

Positive and negative frequency partners are solved independently with
`kpar=.3>0`. Their reflection errors are below `1e-10` (typically `1e-14`).
The accurate-Bessel positive roots are:

| Grid | omega_real | gamma | Complex error from analytic root |
|---|---:|---:|---:|
| 32×64 | .389163073298398 | .085611504181345 | 7.50915e-3 |
| 64×128 | .394779307977957 | .085883952806634 | 1.88634e-3 |
| 96×192 | .395825247533816 | .085939983884691 | 8.38897e-4 |

The analytic target is `.39666286088749+.08598637345950i`. Residuals are of
order `1e-14`. Accurate Bessel values move the 64×128 C2d root by about
`+5.93e-10-3.46e-10i`; the much larger ordinary momentum-discretization error
continues to decrease approximately quadratically.

Both prescribed rectangles, real bounds `[-1,1]`, imaginary bounds
`[.01,.2]` and `[.001,.2]`, return accepted winding **2**, status **0**, with
both Bessel constructors and all three grids. All isotropic numerical-proton
controls return accepted zero on those finite rectangles. Matched comparisons
use the same 257 frequencies, not separately adapted trajectories.

| Grid | Median relative D error | Max relative, original box | Max relative, lower box | Max absolute, original box |
|---|---:|---:|---:|---:|
| 32×64 | .0127795 | .274638 | .528795 | .0940136 |
| 64×128 | .00322331 | .0675131 | .129543 | .0235663 |
| 96×192 | .00143553 | .0298145 | .0586804 | .0104674 |

Accurate-Bessel maximum accumulated phase errors on the original box decrease
from `.215799` to `.0580640` to `.0259886` radians. Total loop phases agree
with `4*pi` to rounding. Absolute errors at the maximum-relative-error
locations are `.00291079`, `.000715546`, `.000315993`. These values distinguish
ordinary grid error near small determinants from a Bessel-induced instability.
`errors.dat` retains lower-box absolute errors and minimum magnitudes too.

Accurate-Bessel adaptive costs, using unchanged A1 controls:

| Grid | Box gamma bounds | Calls | Passes | Refinements | Maximum depth |
|---|---|---:|---:|---:|---:|
| 32×64 | [.01,.2] | 616 | 3 | 84 | 6 |
| 32×64 | [.001,.2] | 702 | 3 | 127 | 9 |
| 64×128 | [.01,.2] | 630 | 3 | 91 | 6 |
| 64×128 | [.001,.2] | 690 | 3 | 121 | 9 |
| 96×192 | [.01,.2] | 630 | 3 | 91 | 6 |
| 96×192 | [.001,.2] | 686 | 3 | 119 | 9 |

## 6. Oblique and parallel-limit controls

The same equilibrium was tested at the single oblique wavevector `(.3,.3)`;
no numerical-VDF scans or growth searches are used. With cutoff 43, the accurate
Bessel maximum finite J-square defect is `7.07e-16`, versus `2.18e-7` original.
The predicted static z contraction changes from `2.93502e-12` to `3.11e-18`.
An independent +8-order run checks truncation without altering the primary
assignments. Tensor assembly retains the source's symmetric/antisymmetric
component structure.

At 64×128 the oblique matched-contour median relative error is `.00365893`;
maximum relative error is `3.02361` where the analytic determinant is only
`4.0289e-7` (absolute discrepancy `1.21819e-6`). Maximum absolute discrepancy
over the contour is `.0203066`. Thus accurate Bessel values are necessary for
static consistency, but do not make a finite-grid VDF exact. No claim of a
newly identified physical oblique mode is made.

Independent fixed-wavevector runs at `kperp=0,1e-8,1e-6,1e-5,1e-4,.001`
check the approach to the existing exact parallel special case. No change to
harmonic assignments or limit formulas is needed. Coefficients remain fixed
with respect to frequency, so the C2d fixed-domain PWL/Cauchy analyticity
argument is unchanged by replacing their Bessel constants. Roundoff near the
origin remains distinct from mathematical holomorphy of the represented model.

## 7. Performance and regressions

Timing is split into the one-time array constructor, unchanged harmonic-limit
determination, reusable Cauchy-weight construction and repeated determinant
calls. The physical matrix is evaluated through the unchanged collective disp;
all four rank call counts are checked at every run. The constructor executes
once on each applicable worker, never inside the frequency loop.

The scalar constructor-like benchmark (100 repeats, 129 perpendicular nodes,
49 orders) measured about `.133 s` original and `.167 s` intrinsic on this
machine. The physical constructor is much smaller; repeated determinant timing
must be compared using the I/O-disabled three-repeat benchmark rather than
diagnostic-heavy contour wall times. Detailed timing tables are retained with
the compact results. No new cache, wavevector parallelism or communicator is used.

The C2d reusable weights remain unchanged. In the representative 257-frequency
64×128 run, the n=0 owner reports 3855 species/harmonic weight sets, 971460
kernel-cell operations, 357788718 derivative-node reads and 64132 peak cache
bytes. Bessel constructor selection does not change these work counts.

The three I/O-disabled 257-frequency repeats measured:

| Constructor | Collective seconds per 257 calls (range) | Mean seconds/call | n=0-worker constructor | n=0-worker cutoff determination |
|---|---:|---:|---:|---:|
| Original | 11.3183–11.4863 | .0443231 | 34.3–35.1 microseconds | 99.3–111.5 microseconds |
| Intrinsic | 11.6301–11.7027 | .0453850 | 39.7–41.9 microseconds | 99.3–103.3 microseconds |

Mean total time increased by about 2.4% in these sequential repeats. No Bessel
call is added in the frequency loop, so this does not establish an intrinsic
per-frequency Bessel cost; the same integration/work counts are retained.
Weight construction took roughly `.353–.358 s` on that worker per trajectory.
Startup array audits and per-integral I/O are suppressed in this benchmark;
the timing constructor reads the existing validation control before startup.
Earlier diagnostic-heavy runs are not used to estimate this overhead.

The following checks passed:

* Optimized and debug independent Bessel value/derivative/identity audits,
  harmonic stopping replay, actual tensor parallel limits, and prescribed
  derivative corner tests.
* Repeated C2d actual-source parallel/static identity and full-tensor identity /
  complex-derivative tests in both configurations. The latter's zz CR residual
  decreased to approximately `4.0e-11` at perturbation `1e-6`; no frequency
  selector is introduced by changing Bessel constants.
* Optimized and debug physical growing partners, both standard contours at all
  three grids, near-origin narrow counts, isotropic controls, representative
  matched 257-frequency trajectories and the four-way tensor checks. Debug and
  optimized oblique trajectories agree exactly at all 257 frequencies.
* Independent raw/scaled low-frequency contours, enlarged cutoff, stricter
  initial sampling, and analytic control windows; near-origin matrix diagnostics
  are interpreted only within their conditioning limit.
* All six existing suites (`run_nyquist`, `run_nyquist_growth`, `run_nyquist_log`,
  `run_nyquist_mirror`, `run_nyquist_scan`, `run_nyquist_plane`) in optimized and
  debug configurations: **12 suite passes**. Their ordinary scientific-output
  comparisons used `ALPS_REFERENCE_EXE` pointing to the identified unchanged
  baseline executable and retain byte identity, including root/map/scan outputs.
  Existing budget, invalid-input, numerical fault and collective shutdown tests
  passed. The interrupted final debug plane run was restarted and completed.
* `collect.sh` acceptance gates, debug/optimized root comparisons, exact corner
  Landau-stream comparisons and the protected-source SHA256 manifest.

The historical entire C2b/C2c/C2d long-running study matrices were not all
regenerated from scratch. Their relevant mathematical tests, retained baseline
root/determinant comparisons and current physical controls were repeated;
the new full C2e entry point assembles the tested components into one fresh run.
This distinguishes tests actually run from a claim of a second full clean replay.

## 8. Adoption recommendation

1. Independently accept the narrowly scoped **int_ee first-corner factor fix**
   for a subsequent production patch. It has a direct continuum derivation,
   isolated prescribed-array tests and matched growing/real/damped checks.
2. The **nonrelativistic numerical-species constructor** can use independently
   validated intrinsic values while retaining the harmonic selection and halo
   layout. This addresses the demonstrated low-frequency defect without
   touching relativistic or analytic calculations. A separate adoption change
   should retain the conditioning limits and physical grid-convergence tests.
3. A replacement of **shared BESSJ** needs its own relativistic and signed-charge
   compatibility review. The severe negative-argument recurrence defect and
   asymptotic phase literals deserve that review; they are not globally patched
   in this study.
4. C2d's coherent UHP susceptibility still needs an explicit production adoption
   decision. Its damped causal continuation remains separate from the direct
   fixed-cell integral. Accurate Bessel constants do not resolve that design
   requirement or the C2b finite-window Tlim mismatch.

The targeted extra finite-growth root is removed in independently converged,
unchanged low-frequency windows. There is no justified claim of unconditional
physical stability, complete global root exclusion, or reliable double-precision
root finding arbitrarily close to omega zero.
