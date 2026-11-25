# C2b — numerical susceptibility continuity correction

The three interior continuity corrections are implemented and tested. The
tables below distinguish their validation from the remaining finite-window
`Tlim`, support-edge and finite-grid accuracy limitations. Final regression
builds are identified in the validation inventory below.

## Implementation and scope

Only `src/ALPS_fns.f90` changes production mathematics. `funct_g` now uses a
continuous cell interpolant; the original reconstruction is retained as public
`funct_g_legacy`. Private helpers `funct_g_node` and `integrate_cell_part` assemble
the original nodal numerator and integrate fractional background cells.
`integrate_res` keeps its moving symmetric window, with a constant width and
consistent fractional tails instead of the previous sawtooth width/remainder.

No changes are made to `integrate`, `full_integrate`, `landau_integrate`,
`determine_resonances`, `resU`, `int_T`, numerical differentiation, continuation
fits, analytic bi-Maxwellian response, `disp`, MPI communication/decomposition,
A1, A3, root solvers or wavevector scans. The `Tlim` threshold/default, selector
`abs(denomI).GT.Tlim`, regular symmetric-core formula, small-imaginary formula
and limiting i*pi*G terms are unchanged. Function extraction verifies byte
identity of `integrate`, `full_integrate` and `landau_integrate` against the
pre-correction source. Renaming the legacy interpolator back to `funct_g` and
removing its added documentation gives byte identity with the original routine.

The correction was made and tested in three checkpoints. Validation-only source
variants reproduce each checkpoint through the same collective `disp`; there is
no production legacy switch or second susceptibility implementation.

## Reproduction and retained baseline

The original `ALPS_fns.f90` SHA-256 was
`567e05bad770489d578c0024b4230376f1e000c2f04b7e84fc2f7d065d446be5`.
Before any production edits, the C2a constant-numerator tests and full collective
half-cell probes were repeated. Four ranks, P=3, N=100, Tlim=.01, method-2
log-Chebyshev continuation/order8, normal numerical differentiation, and the C1
beta_parallel=3, alpha_p=2 equilibrium were retained. Proton support is
0<=p_perp<=4 sqrt(6), |p_parallel|<=4 sqrt(3); electron response is analytic,
isotropic beta=3 with mass .000544617021487, k=(.001,.3), vA=1e-4.

At omega_c=-.0162379763209582, across ±1e-12:

| gamma | original relative determinant difference | interpolation + remainder only | combined correction |
|---:|---:|---:|---:|
| .001 | 2.38779629123 | .01787189508 | 3.3860e-10 |
| .01 | .009428216625 | .00568599050 | 1.9834e-10 |
| .03 | .002575655213 | .00190010423 | 1.1130e-10 |

These final differences are finite-separation frequency variation, not estimates
of a nonzero limiting jump. Decreasing delta by 100 decreases the differences
by approximately 100 until roundoff becomes appreciable. Integer windings were
not used to select or accept the correction.

`tests/data/nyquist_c2b/legacy_halfcell.txt` retains the original complex values;
`legacy_roots.txt` retains full-precision growing/damped ordinary-solver roots
and residuals on all three grids. C2a fixtures and observations are unchanged.
The original standalone normalized constant-numerator transferred-cell jump
reproduces -0.002643860880366 -0.000160817581835i. The original omitted-remainder
example reproduces -0.00994529849310 -0.000227332356028i.

## Checkpoint 1: positive remainder

The first correction used the actual rounded endpoints:

    start = real(p_res) + Delta
    R = p_upper - start
    ntiny = 0
    if (R > 0) ntiny = max(1,int(R/smdelta))

No epsilon cutoff is introduced. For positive remainders the original panel
convention is retained whenever it already yields a panel. The final endpoint
is assigned directly from p_upper, avoiding accumulated endpoint drift. The
valid interior window has Delta>0 and N>0; no division is performed for R<=0.
The subsequent partition correction supersedes this loop with a cell primitive
that includes every positive fractional interval without a panel-count switch.
The isolated intermediate implementation remains reproducible in validation.

`remainder.f90` exercises the actual intermediate loop with constant, linear and
quadratic numerators; R/step from 1e-10 through 2, including either side of 1;
positive/negative imaginary displacement; and both Tlim branches. It agrees
with an independent dense integral to the expected trapezoidal error. The old
ntiny=1->0 omission jump becomes -2.017e-12 -0.210e-12i at delta=1e-12.
The final cell primitive separately tests a zero-length interval and a positive
one-ULP interval. Perpendicular endpoint weights are unchanged.

At 64×128, the positive growing root is unchanged by this isolated correction;
the damped root changes only at rounding level. All recorded Landau values at
13 identical growing/real/damped frequencies are byte-identical.

## Checkpoint 2: continuous interpolation

At neighboring usable parallel nodes j and j+1, interpolate the **complete
complex numerator**, including the existing frequency, species and tensor
factors:

    t = (p-p_j)/(p_(j+1)-p_j)
    G(p) = (1-t) G_j + t G_(j+1).

The usable nodes are 1..npar-1; the first/last usable cell supplies endpoint
extrapolation. A uniform-grid prediction and neighboring-node correction give
constant-time indexing. Bounds are checked before FLOOR; there are no per-call
allocations or linear searches. Each numerator is assembled from the existing
numerical derivatives, not a separately interpolated VDF.

Tests cover constant, linear, quadratic, cubic, Gaussian and sine numerators,
including complex factors; nodes, quarter-cells, midpoints, both boundaries,
and one-sided limits. Constants and linear functions are exact within rounding;
all nodal values are recovered. Gaussian maximum interpolation errors decrease
by factors 3.96 and 3.99 under successive spacing halvings. The legacy cubic
jump remains exactly -1.5*h^3 times its complex factor (within finite-probe
rounding); the corrected differences tend to zero.

Interpolation alone does not remove the independently diagnosed quadrature
mismatch, as the checkpoint table demonstrates. It was therefore validated
before changing the partition.

## Checkpoint 3: consistent moving partition

The proposal in `validation/c2b/PARTITION_PROPOSAL.md` was written and the copied
implementation passed its standalone tests before production adoption.

Keep a symmetric core [r-Delta,r+Delta], with
**Delta=max(P,2)*h independent of r** in the interior (P*h for the
reference/default P=3), bounded by the available usable momentum support.
This removes the previous sawtooth reset at grid nodes; merely repairing the
half-cell transfer would leave that reset. The existing core formula and N
subintervals are retained, for both Tlim branches and both signs of Im(p_res).

The two tails retain the original whole-cell trapezoid. In a fractional cell,
use the primitive of that same linearly reconstructed **full integrand**
H_j=G_j/(p_j-p_res), rather than a different subdivided G(p)/(p-p_res) rule:

    int_a^b H dp = h [(t_b-t_a) H_j
                     + (t_b-t_a)(t_b+t_a)(H_(j+1)-H_j)/2].

Two fractional pieces reproduce exactly the original whole-cell trapezoid,
up to rounding. A cell changing ownership therefore does not change its
quadrature representation. Left tail, core and right tail have disjoint
interiors and cover the original usable momentum domain. No cell is omitted or
duplicated, no smoothing/blending is applied, and no phase contribution is
spliced. The old remainder-panel count is absent; every positive tail is
included, with zero returning zero. The original perpendicular weights and
2*pi normalization are retained.

The prototype tests the actual routines for constants, polynomials and a
Gaussian, both imaginary signs and both Tlim branches. It checks cell splitting,
zero/one-ULP intervals, domain length, and one-sided limits at the half-cell,
grid node, quarter-cell, old ntiny threshold and a moving internal-node event.
All pass in optimized and debug configurations, with bounds/FP checks.

The legacy special returns near the finite support edges remain unchanged.
These were not the C2a central discontinuities and are not certified continuous
by this correction. This limitation is distinct from the corrected interior
partition and requires a separate support-edge audit.

## Full susceptibility continuity and causal assembly

For each 32×64, 64×128 and 96×192 grid, five central event types are probed at
delta=1e-6,1e-8,1e-10,1e-12 and gamma=±.001,±.003,±.01,±.03,0. This is 360
frequencies per equilibrium, repeated for anisotropic and isotropic protons.
All six numerical runs pass the shrinking-separation test. The largest relative
differences at delta=1e-12 are:

| grid | anisotropic | isotropic |
|---|---:|---:|
| 32×64 | 2.134e-9 | 1.329e-9 |
| 64×128 | 1.868e-9 | 1.251e-8 |
| 96×192 | 7.534e-9 | 4.545e-9 |

These maxima include near-origin/near-small-determinant cases, not just the
original half-cell. Each sequence decreases with separation; no finite plateau
is observed. Analytic raw determinants near z=0 become extremely small and
poorly conditioned; relative comparisons there are not a reliable continuity
metric. Absolute values are retained rather than interpreting their relative
roundoff as a numerical-VDF jump.

Per-rank files retain all harmonic/mode `full_integrate`, `integrate_res` and
`landau_integrate` returns, and assembled species tensors. The C2a analyzer
checks unique harmonic ownership and tensor closure. The original dominant
proton n=0, zz contribution is retained explicitly in those files.

At identical 13 frequencies on each grid, every recorded Landau contribution
is byte-identical before/after correction on all four ranks. The independent
assembly checker verifies actual returned values against I_res, I_res+I_L and
I_res+2*I_L for growing, real and damped frequencies: respectively 126, 21 and
126 checked harmonic/component records per grid. This tests causal assembly,
not merely successful ordinary root convergence.

## Growing and damped roots

The unchanged ordinary solver (two independent guesses, secant_method=0,
D_threshold=1e-13) gives analytic references:

    growing: +0.39666286088749 + 0.08598637345950 i
    damped:  +0.40460138158784 - 0.00185430537032 i.

Negative-frequency partners were solved independently. Corrected results:

| grid | growing positive-frequency root | isotropic damped positive-frequency root |
|---|---|---|
| 32×64 | .390581513199 + .086311784789i | .407435771470 - .001312181370i |
| 64×128 | .395411724579 + .086220209553i | .405337298402 - .001709481504i |
| 96×192 | .396168900706 + .086138999817i | .404931367075 - .001772382098i |

Maximum complex root errors decrease from .00609005 to .00127280 to .000517002
for growing roots, and .00288577 to .000750032 to .000340003 for damped roots.
All residuals are below 2e-12; the measured values are normally below 1e-13.
The negative-frequency solutions recover reflection symmetry to roughly
1e-12 or better, unlike the small legacy interpolation/partition asymmetry.

At 64×128, isolated checkpoints give:

| approximation | growing root | damped root |
|---|---|---|
| legacy | .395485626 + .086310360i | .405307785 - .001742334i |
| remainder only | same within rounding | same within rounding |
| + interpolation | .395339447 + .086164369i | .405324262 - .001712554i |
| + partition | .395411725 + .086220210i | .405337298 - .001709482i |

Continuity does **not** imply a smaller error at every fixed resolution. The
new linear interpolation can have larger ordinary truncation bias than the
legacy centered tangent; in particular the finite-grid damped rate is slightly
farther from the analytic value. Both corrected root sequences move toward
the physical references over the three tested grids. Remaining Tlim/finite-window and
continuation-fit errors have not been repaired or hidden. This is an observed
resolution trend, not proof of asymptotic convergence: with fixed Tlim and P,
Delta=max(P,2)*h shrinks under refinement while the retained small-imaginary formula
still adds its infinite-window residue. C2c must address that limit before
claiming arbitrarily fine-grid convergence of the damped response.

## Reproducibility

```sh
bash validation/run_nyquist_c2b.sh
FFLAGS='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow' \
  bash validation/run_nyquist_c2b.sh
```

See `validation/c2b/README.md`. Generated distributions, inputs, sources,
compiler settings, counters, full-precision roots, matched ASCII values and
per-harmonic/tensor traces remain in the printed external scratch directory.
No Python or new numerical dependency is used. GNU Fortran 11.4.0 and Open MPI
4.1.2 are the identified toolchain; all physical evaluations use four ranks.

## Existing regression gate (completed)

All six existing entry points passed sequentially in **both** optimized and debug
builds, without modifying their analytic fixtures:

| script | optimized | debug |
|---|---|---|
| `validation/run_nyquist.sh` | PASS | PASS |
| `validation/run_nyquist_growth.sh` | PASS | PASS |
| `validation/run_nyquist_log.sh` | PASS | PASS |
| `validation/run_nyquist_mirror.sh` | PASS | PASS |
| `validation/run_nyquist_scan.sh` | PASS | PASS |
| `validation/run_nyquist_plane.sh` | PASS | PASS |

The A2 script includes A1 standalone tests and compares nine ordinary root/map/
line/plane scientific files byte-for-byte against the identified pre-C2b
executable. B1/B2 checks also preserve option-3/4 scientific outputs, validate
scan-state reconstruction and aligned compact/detailed files, and exercise
collective nonfinite/budget failures. A3 checks cover cache/strip equivalence,
both growth spacings, and failed searches. The growth suite reports 10,056
mathematical/cache checks. No MPI worker groups, communicators or resonance
assignment routines change.

Logs and exit statuses:
`/tmp/alps-nyquist-c2b/regressions/status.txt` and adjacent per-suite logs.
Optimized flags: `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`.
Debug flags: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`.
Standalone correction tests additionally passed with invalid/zero/overflow traps.
These analytic regression passes do not substitute for the correction-specific
numerical and C1 checks identified separately below.

## Measured remaining Tlim jump and checkpoint cost

At omega=.4 and gamma=.003±1e-12 (default Tlim), the one-sided relative
**positive-gamma** differences are .0205536797543 (legacy), .0205536797543
(remainder only), .0185282542929 (interpolation), and .0209079976959 (combined).
The corresponding negative-gamma switch differences are .0240773738607,
.0240773738607, .0214158612710 and .0219239782057. These are nonzero limiting
jumps. C2b does not remove them: the finite-window small-imaginary approximation
and its branch switch remain C2c work. Changing the core partition can change
their magnitude even though the formulas/threshold remain unchanged.

An I/O-disabled four-rank matched-contour benchmark evaluates the **same 257
frequencies** at 64×128 for every checkpoint. All ranks report 257 physical
calls. Growth/caching stays disabled for numerical species as in C1; no new
cache policy is introduced. Times exclude equilibrium preparation and are single
instrumented measurements, not general performance predictions:

| variant | G interface calls (all ranks) | collective disp wall seconds | resonance-integrator seconds |
|---|---:|---:|---:|
| legacy | 62,207,649 | 14.0998 | 9.2531 |
| remainder only | 62,220,501 | 14.1057 | 9.2266 |
| + continuous interpolation | 62,220,501 | 15.0629 | 9.1252 |
| + consistent partition | 56,999,754 | 15.7258 | 9.2570 |

Collective time is the maximum summed `disp` time across ranks, not the sum of
four simultaneous wall times. The integrator column sums per-rank time inside
the actual `integrate_res` calls; it excludes observer I/O and the Landau routine.
No Landau contribution is needed on this upper-half-plane benchmark. Ordinary
damped-root runs also time the unchanged continuation evaluation. Counting adds
small common instrumentation overhead. The new interpolator assembles two
nodal numerators instead of three; the new partition removes subdivided tail
requests. These are consequences of the tested correction, not a reduction of
N, contour checks or requested physical evaluations. No allocation/search is
added to the inner loops.

These are the repeated optimized run's measurements. Despite fewer numerator
calls, its combined collective wall time is about 12% higher than legacy;
resonance-integrator time is essentially unchanged. Earlier isolated repeats
gave 12.48 versus 10.45 seconds collectively, so a reliable speedup is not
established. Scheduling and MPI waiting contribute to these short timings.
The correction is accepted for continuity, not a performance claim.

The standalone constant-numerator value check also compares the complete
integral with log(p_max-a)-log(p_min-a). At h=.125, r=-h/2 and P=3, the regular
b=.02 result differs by .00215125 per perpendicular factor; the small-imaginary
b=.005 result differs by .0272066. Negative b gives the conjugate result.
Continuity is therefore not being confused with exact quadrature or a repair
of the deliberately retained finite-window approximation.

The primary 64×128 anisotropic adaptive contour independently accepts W=2,
phase=4*pi, in three passes with 90 local refinements, deepest level 6 and
maximum final increment .4840813. It makes 628 collective determinant calls
(25.3533 seconds for this instrumented run, samples enabled). The pre-C2a
original executable retained W=2 with 626 calls on the same input. Changes in
adaptive requests follow from the changed numerical trajectory, not A1 changes.

A coarse-grid warning emerged in the matched **fixed** contour: the corrected
32×64, 257-point trajectory accumulates 6*pi rather than the analytic 4*pi.
This fixed pass alone is not an accepted adaptive winding; the subsequent
independent adaptive C1 calculation also accepts W=3 (details below). At its bottom-edge
center the relative determinant discrepancy exceeds one, so equal physical
root counts cannot be inferred. The independently converged primary 64×128
count above is unaffected. This discrepancy is retained rather than clipped or
corrected by changing the contour.

## Damped-response accuracy and a remaining real-axis matching issue

The 13 matched stable-control frequencies span negative, zero and positive
imaginary parts at Re(z)=.4046. Near the independently determined analytic damped
root, relative determinant errors are large because the reference determinant
is almost zero, not because an evaluation is nonfinite:

| grid | worst relative error | absolute error at that point |
|---|---:|---:|
| 32×64 | 2131.76 | .00122983 |
| 64×128 | 539.505 | .000311246 |
| 96×192 | 243.128 | .000140263 |

All maxima occur at .4046-.0018543i, where |D_analytic|=5.76910e-7. Both absolute
errors and independently solved root errors decrease with refinement. This
finite-grid response is not claimed accurate to the solver's residual tolerance.
At 64×128, the existing derivative_type=2 diagnostic (same method-2 fit/order)
changes the maximum relative error to 534.681; Tlim=.005 changes it to 565.220.
Neither input change is adopted as a production default or used to force a root.
The corresponding before/after Landau terms remain byte-identical at identical
frequencies for each diagnostic input.

An additional near-real-axis probe, z=.4046±i*eta, eta=1e-4..1e-12, exposes a
remaining **value mismatch between the two causal limits**. At eta=1e-12 the
corrected determinant differences are:

| grid | absolute difference | relative to the negative-side value |
|---|---:|---:|
| 32×64 | .000433869 | .327658 |
| 64×128 | .0000971167 | .124145 |
| 96×192 | .0000404484 | .0530789 |

This mismatch already exists in the legacy response (e.g. absolute .000316177
at 32×64); it is not an altered Landau contribution. It decreases with grid
refinement, but must not be hidden behind successful damped-root convergence.
Raw per-harmonic data and tensor closure confirm the actual causal assembly.
The largest raw component difference is n=0, zz; the largest determinant-weighted
individual contribution is n=1, xy, with substantial cancellations.

Mathematical deduction: as Im(z) tends to zero, the small-imaginary core supplies
±i*pi times the **numerically reconstructed** numerator. The unchanged Landau
routine evaluates the **fitted continuation** and its derivatives, with its
own perpendicular endpoint rule. Exact matching requires these two residues to
agree at finite resolution. They use different derivative/interpolation and
endpoint approximations, so that identity is not enforced. The observed
mismatch is consistent with this difference; separating the derivative and
endpoint effects further is follow-up work. No Landau weights, fitted
representation or causal addition is changed here. This is another reason not
to equate the corrected discrete trajectory with a holomorphic dispersion
function without further convergence evidence.

For a fixed Tlim branch and an interior window, the corrected real-frequency
continuity has a different mathematical basis: each fixed core node samples a
continuous G(r+u), while the two tails share one continuous cell primitive.
Changing cell ownership therefore cannot create a finite value jump in exact
uniform-grid arithmetic. The tests check that implementation down to floating
precision. This local C0 property does not resolve Tlim or causal-limit matching,
nor establish the Cauchy–Riemann conditions.

## Matched-boundary convergence and the physical numerator

At the same 257 ordered contour coordinates in every case:

| grid | legacy median/max relative D error | corrected median/max | corrected max cumulative phase difference |
|---|---|---|---:|
| 32×64 | .00924375 / 1.16099 | .0101841 / 1.25825 | 6.49816 rad |
| 64×128 | .00408394 / .212333 | .00430580 / .137389 | .0582262 rad |
| 96×192 | .00234069 / .0650747 | .00267309 / .0329546 | .0269974 rad |

Median and worst corrected errors decrease with resolution. The worst sampled
points are near the bottom-edge center where the analytic determinant is small.
At 64×128 and 96×192 the corrected worst errors improve on legacy values; the
median remains slightly worse, consistent with the interpolation's truncation
bias. The 32×64 fixed phase is 6*pi; the other corrected fixed phases and all
analytic references are 4*pi. These fixed paths do not themselves supply an
accepted adaptive winding or a proof of a uniform boundary error bound.

The primary isotropic adaptive control independently returns W=0 with 632
collective calls, three passes, 92 refinements, deepest level 6 and 27.0881
seconds (samples enabled). A valid zero remains distinct from an evaluator or
sampling failure.

At the physical 64×128 half-cell event, gamma=.001, proton n=0/zz,
p_perp index 8, the new numerator tends to
-5.185798423e-6 + 1.065451204e-7 i from both sides. The legacy one-sided values
remain +6.710923443e-6 + 2.482576516e-9 i and -3.455721280e-6 + 0i.
For delta=1e-6,1e-8,1e-10,1e-12, the corrected numerator differences are
8.51916e-10, 8.51916e-12, 8.51916e-14 and 8.51914e-16. Legacy differences tend
to 1.01666450e-5. The physical test therefore verifies the actual reconstructed
numerator, not only prescribed polynomial examples.

The corresponding raw n=0/zz `full_integrate` values (before the separate,
frequency-independent ee term and species assembly) are:

| approximation | left value | right value |
|---|---|---|
| legacy | .498642948124 + .000232641570i | .498537727289 - .000396425424i |
| corrected | .501114145295 - .000511477387i | same to about 4.84e-13 |

All remaining harmonics/components, assembled tensors and determinant values
are retained in the per-rank diagnostic files. The original contribution's
finite jump is removed; ordinary finite-grid value changes are not concealed.

## Final positive-width roundoff safeguard

A late standalone corner test found a weakness in the first implementation of
the cell primitive: subtracting independently normalized t values can give zero
for distinct physical endpoints. For h=sqrt(3)/16,
lo=.00010825317547305482 and hi=nearest(lo,+1), the positive span is
1.3552527156068805e-20 but both lo/h and hi/h round to .001. The actual first
primitive fails the added positive-interval test; that implementation was not
retained as the completed default.

The final code evaluates the **same** linear-cell integral as

    (hi-lo) * [(1-tmid) H_j + tmid H_(j+1)],
    tmid = .5*tlo + .5*thi.

It never derives the measure by subtracting normalized coordinates. The test
covers 999 adjacent-endpoint pairs on this nonbinary grid for both imaginary
signs; 93 have coalesced normalized coordinates on the identified compiler.
Every positive interval contributes for the constant-numerator fixture, while
zero intervals return zero. Both optimized and invalid/zero/overflow-checked
debug standalone runs pass. This is a roundoff-safe evaluation of the same
primitive, not a new quadrature or an accuracy/performance tuning change.
The final regression inventory below identifies the repeated correction suite.

The final-source rerun caught an instrumentation issue: a retained worker counter
reported 147 Landau evaluations while its stream had 111 rows. The deliberately
nonzero fixed-probe exit can let the MPI launcher terminate a worker before
implicit Fortran unit cleanup. No Landau values changed; an incomplete file
failed the strict comparison. The copied observer now closes every diagnostic
unit **before** the existing MPI_Finalize call. It introduces no collective or
production change. Per-stream row counts are additionally checked against
counters; the byte-equality gate is retained. Final validation is rerun with
complete streams rather than relaxing that gate.

The debug rerun independently detected the same observer-buffer issue in a full
harmonic stream: all ranks counted 360 determinant calls, but rank 2 had 34,542
full-integral rows instead of 360*96=34,560. The final frequency was missing 18
component records. The analyzer rejected the inconsistent inventory. Thus this
was incomplete output, not different collective determinant-call sequences.
The new full-call counter and pre-finalization closes cover this stream too.
In these reference runs rank 1 retains n=0..7 numerical proton work (87 full
components per frequency), rank 2 retains n=8..15 (96 components), and the
analytic electron work remains on its existing assignment. No assignment,
communicator or frequency-worker group was changed.

## Related C1 regression expectations

Both final-source builds independently accept **W=3 at 32×64**, with 628 calls,
three agreeing passes and deepest level 7. This confirms that the fixed-path
coarse warning is not solely a 257-point alias. The old C1 checker failed its
historical W=2 expectation. That input and the original report are retained;
only this directly affected numerical expectation is updated to a strict W=3
regression. It does not identify three physical unstable modes. Finer grids
remain required to approach the analytic response and W=2.

The unchanged low-gamma input (gamma_min=.001) now sampling-converges to W=2
in both builds: 690 calls, three passes, 121 refinements and deepest level 8.
Its former failure was caused by the corrected real-frequency jumps. The
remaining Tlim discontinuity means this result is not an analyticity
certificate. The related C1 expectation is updated to that measured value.
A separate numerical-VDF max_depth=0 case preserves a controlled refinement-
failure regression: status 5, **no** accepted winding, three collective calls,
normal shutdown. Budget and nonfinite-failure checks are retained unchanged.
No original physical input, contour, derivative, support or tolerance is shifted
to make a numerical count agree with the analytic one.

Ordinary damped-root sensitivity at 64×128 is also measured, using the same
method-2/order8 continuation representation and unchanged physical input:

| control | legacy positive-frequency root | corrected root |
|---|---|---|
| default derivative_type=1, Tlim=.01 | .405307785-.001742334i | .405337298-.001709482i |
| derivative_type=2, Tlim=.01 | .405313952-.001934968i | .405343387-.001901732i |
| derivative_type=1, Tlim=.005 | .405369162-.001750402i | .405367638-.001721746i |

Both signed-frequency roots converge with small determinant residuals. The
corrected partners recover reflection symmetry within rounding; the legacy
partners retain their finite-grid asymmetry. Damping is sensitive to how the
existing fitted derivative is evaluated, even though the fit itself is the
same. These controls are diagnostics, not new defaults or proposed remedies
chosen to obtain a preferred root/count. Their input generators and root
comparisons are included in the C2b suite.

## One-cell-window real-frequency safeguard

Testing a valid nondefault P=1 at a real resonant node exposed another corner of
the first partition prototype: a tiny fractional background interval can borrow
the pole as its other interpolation node. The background fragment itself does
not contain the pole, but H_j at that borrowed node is undefined. This is an
introduced quadrature-stencil problem, not a physical susceptibility pole.
The original integrator returns a finite (biased) PV integral for the constant
numerator; the first prototype triggers a divide-by-zero check.

A copied prototype with **Delta=max(P,2)*h** passes before production adoption.
The two-cell minimum places both background stencil nodes outside the resonance.
The original P-based support-edge returns remain unchanged. The width is bounded
by the usable support and the fractional-cell indices remain clamped; NaN
sentinels at the unused derivative endpoints test those limits. The default P=3 and every physical reference result above are
unchanged; the nondefault P=1 now uses the minimum safe two-cell core. This
behavior is documented under positions_principal in input.md.

`real_node.f90` checks P=1,2,3, parallel grids 64/128/192, real frequencies ±1
with an exact resonant node, gamma=0 and both signs of small/regular imaginary
displacement. It verifies finite values, the zero PV of a constant numerator on
a symmetric interval, and complex-conjugation symmetry. Both optimized and
invalid/zero/overflow-checked debug prototypes pass. This guard leaves the core
formulas, N, Tlim, Landau addition and physical momentum-domain bounds unchanged.

## Completed numerical gate and build inventory

The final sequential queue completed with all four exit statuses zero:
optimized C1, optimized C2b, debug C1, debug C2b. The numerical correction suite
includes all three grids, both signed ordinary roots, matched causal/Landau
records, six 360-frequency interior-event runs, and independent adaptive full
and stable rectangles. The C1 gate also checks unsupported backends, invalid
inputs, budget/depth/nonfinite failures and collective termination.

Retained results: `/tmp/alps-nyquist-c2b/regressions/closed-status.txt`,
`closed_{optimized,debug}_{c1,c2b}.log`, and
`/tmp/alps-nyquist-c2b/closed_suite_{optimized,debug}/`.
The optimized complete numerical suite used the positive-physical-width source
`d2c5b9f8eb15fcbc9950c83daef52e8901476f6309d84ed30817162bff85f1c4`.
The debug complete suite used final source
`905416da5e1a2dfdb03974bb49e167595b1f11ab70e7c1b77ca1a3c1f0058435`.
The later minimum-two-cell and support/stencil bounds safeguards leave default
P=3 interior physical results unchanged. Final-source partition and real-node
standalone tests additionally passed in both optimized and debug builds with
invalid/zero/overflow traps; logs are
`/tmp/alps-nyquist-c2b/final_math_{optimized,debug}_{partition,real_node}.log`.
The earlier complete analytic regression matrix and repeated A2 byte-comparison
passes are identified above; those analytic paths do not invoke these numerical
helpers. No claim is made that every earlier executable had the final hash.

C2b's tested interior continuity corrections are complete. This is not approval
to trust arbitrary coarse numerical VDFs: the accepted coarse W=3, retained
Tlim jump, finite-support fallback and fitted/grid residue mismatch remain
explicit limitations. C2c should compare the exact finite-window linear core
and an independently tested fixed-cell Cauchy representation, preserving this
C2b default until a method is selected and separately authorized.
