# C2c — experimental nonrelativistic Cauchy integration

This is a method-selection study, **not a production-default change**. The
completed C2b implementation remains the production reference. All experimental
susceptibility changes are applied to copied modules under isolated scratch
builds. No namelist method selector is introduced.

**Decision C — continue mathematical development before production adoption.**
A removes the known infinite-window residue error but retains a finite branch
mismatch. B supplies an exact, holomorphic upper-half-plane PWL Cauchy kernel,
but its naive coupling to unchanged tensor quadrature fails a static
cancellation and introduces a grid-sensitive low-growth zero. Develop a
structure-preserving fixed-cell response next; retain C2b defaults meanwhile.

## Baseline gate and implementation boundaries

The outstanding C2b gate completed before interpreting C2c physical results:
C1 and C2b passed in both optimized and debug builds. The six existing analytic
Nyquist/root/map/scan suites had also passed both configurations, including
nine byte-identical legacy scientific outputs. See `NYQUIST_C2B_RESULTS.md` for
identified source versions, the coarse-grid W=3 warning and retained limitations.
The C2c production baseline SHA-256 of `src/ALPS_fns.f90` is
`905416da5e1a2dfdb03974bb49e167595b1f11ab70e7c1b77ca1a3c1f0058435`.

Experimental code lives in `validation/c2c/`; the entry point is
`bash validation/run_nyquist_c2c.sh`. `variants.awk` changes only copied
nonrelativistic numerical integration paths. `build_variant.sh` links them to
the ordinary collective `disp`, MPI setup, harmonic assignments, tensor
assembly, normal derivatives and method-2 continuation. The C1 callback supplies
identical fixed frequencies; adaptive calculations use unchanged A1. The
single-wavevector numerical restrictions and disabled numerical growth scans
remain in force. `funct_g`, `funct_g_legacy`, `landau_integrate`, fitting,
relativistic routines, normalization and production inputs are untouched.

## Signs and actual numerator

Source identities, rather than an assumed sign convention, define

    a = (m_s*omega - n*q_s)/kpar,
    G_j = -q_s/kpar * [omega*df_perp,j
          + (kpar/m_s)*(p_perp*df_parallel,j-p_parallel,j*df_perp,j)] * T_j.

Then `resU*int_T = G_j/(p_j-a)`. `int_T` contains the existing Bessel/tensor
factors and has no frequency dependence at fixed wavevector. Thus G_j is affine
in complex omega for fixed numerical derivatives. Neither species sign nor
harmonic changes this conclusion. In particular, b=Im(a)=m_s*Im(omega)/kpar:
positive growth maps to b<0 when kpar<0. The direct real-momentum transform and
the physical damped continuation must be distinguished.
The existing fitted residue uses `abs(kpar)` in its charge factor. With a
consistent representation it would be i*pi*sign(kpar) times the perpendicularly
integrated grid numerator. Thus physical upper-/lower-omega boundary differences
carry sign(kpar), while upper-/lower-a differences carry the usual positive
2*pi*i factor. Neither candidate substitutes sign(Im(omega)) for sign(b).

## A: finite-window linearized core

For t=p-r and G(r+t)=G0+G1*t,

    J0 = integral[-Delta,Delta] dt/(t-i*b) = 2*i*theta,
    theta = sign(b)*atan2(Delta,abs(b)),
    J1 = integral[-Delta,Delta] t*dt/(t-i*b) = 2*(Delta-b*theta).

Therefore I=J0*G0+J1*G1, with either sign of b. At b=0 the principal value is
2*Delta*G1; the one-sided limits are that value plus/minus i*pi*G0. The existing
one-/two-Landau additions at real/damped omega remain separate and unchanged.
No ratio Delta/b is formed. For Delta much smaller than abs(b), the slope
coefficient uses its convergent alternating series to avoid cancellation.

`finite_core` implements this identity. `integrate_core_small_legacy` preserves
the actual C2b infinite-window constant residue and N-panel trapezoidal slope
calculation, including its b=0 endpoint convention. `test_alps.f90` compares
that scalar reference with the actual C2b `integrate_res` core after subtracting
its unchanged tails, for constants/linear numerators, both kpar signs, both
imaginary signs and PV. The original production small-core body is also retained
unchanged in the baseline and forced-small copied builds.

The A candidate retains the C2b G1 estimate
[G(r+h)-G(r-h)]/(2*h), the moving core, tails, Tlim selector and all other rules.
Only the selected small core is replaced; forced-small, forced-finite and
forced-regular copies compare approximations at identical frequencies. Exact
integration of one linearization does not make that linearization exact over
several PWL cells. Curvature and regular-quadrature errors remain separate.
In particular this moving-window construction has no general holomorphy proof.

Standalone tests span abs(b)/Delta=1e-6,1e-4,1e-2,.1,1,10, both signs and complex
coefficients. Independent polynomial antiderivatives verify constants/linears;
quadratic/cubic/Gaussian rows retain approximation errors. The regular N=100
rule becomes badly under-resolved for narrow resonances; forcing it is a
diagnostic, not a proposed default. Both optimized and debug/FP-checked tests
pass. `core_errors.dat` records legacy, exact-linear and regular errors.

## B: fixed-cell exact PWL transform

For each fixed usable parallel cell, h=p_(j+1)-p_j and s=(G_(j+1)-G_j)/h:

    I_cell = s*h + [G_j+s*(a-p_j)] * L,
    L = log(p_(j+1)-a)-log(p_j-a).

The prototype covers the complete C2b usable domain p_1..p_(npar-1), including
all nonresonant harmonics and resonances outside support. There is no found_res
switch for its direct integral. The same perpendicular rule is retained:
2*pi*dp_perp times weights 1 at the first usable node and interior nodes,
and 1/2 at the last usable node. Distribution endpoints/derivatives are not
invented or extended to change that historical physical domain.

For b!=0 both log arguments lie in the same open half-plane, so their principal
log difference is the continuous branch along the real cell. This statement
is verified by affine antiderivatives, signed tests and boundary-jump tests,
not just assumed from the compiler. Near q=h/(p_j-a)=0 the implementation uses

    L = sum[k>=1] (-1)^(k+1)*q^k/k,
    E = 1-L/q = sum[k>=1] (-1)^(k+1)*q^k/(k+1),
    I_cell = G_j*L + (G_(j+1)-G_j)*E.

Forty terms for abs(q)<1/8 avoid both log subtraction and far-field slope
cancellation. Elsewhere the algebraically identical log difference is used.
There is no alternative quadrature across this numerical arithmetic switch.
Far-field tests retain a 1e-150-sized constant transform accurately.

At real internal nodes the domain PV subtracts G(r) before integrating and adds
G(r)*log[(p_max-r)/(r-p_min)]. This avoids individually divergent cell endpoint
terms and 0*log(0). A support endpoint with nonzero numerator is genuinely
logarithmically singular; the prototype returns NaN there instead of inventing
a finite value. Approaching endpoints with nonzero imaginary displacement is
finite and tested down to abs(b)=1e-100. A production version would need a
clear status-returning interface for such singular real endpoints.

## Analyticity and causal boundary values

**Mathematical result:** fixed support, fixed interpolation weights and affine
G_j(omega) imply a holomorphic integral on either open half-plane of a, by
integrating holomorphic functions over a finite real interval with denominators
uniformly separated from zero on compact subsets. Finite sums, the existing
frequency-independent perpendicular/harmonic factors and determinant assembly
preserve this property in the upper omega half-plane (kpar!=0). The fixed
truncation changes accuracy, not this analyticity statement. For kpar<0 the
same argument uses the lower a half-plane. The numerical branch arithmetic
approximates this one mathematical integral; finite-precision error is tested
separately.

Standalone independent real/imaginary perturbations check Cauchy–Riemann
residuals at cell nodes, mid/quarter cells, support edges and outside support,
for both k signs. Perturbations stay in their chosen open half-plane. Early
exploratory steps crossing the real axis were rejected as invalid analyticity
tests, not used to relax tolerances. Shrinking-separation probes check all
internal nodes and finite-b support crossings. The direct upper/lower limiting
jump agrees with 2*pi*i*G_PWL(r); these scalar tests pass debug and optimized.
Collective determinant CR results and physical accuracy are recorded below.

The direct lower-a transform is **not** the damped causal response. The isolated
B build adds the unchanged fitted ALPS Landau term using the original real/
damped-omega convention for the same interior-resonance fixtures. Its direct
integral is independent of found_res; that flag remains only in the preserved
experimental damped residue addition. No lower-half-plane holomorphy is claimed
for that combination or for arbitrary support crossings.

A diagnostic matching residue can make the two real-axis limits agree if it
uses the same PWL boundary value and perpendicular rule. It is not a global
analytic continuation: selecting a complex linear polynomial by Re(a) can
switch to another polynomial at a cell boundary. PWL slopes generally differ,
and their complex extensions need not agree away from the real axis. A valid
complete damped backend requires an explicitly consistent analytic
representation/continuation, or a separate justified continuation strategy.
No such cell selector or fitted-residue correction is adopted in production.

## Standalone accuracy and numerical conditioning

For G0=1+.3i, G1=0 and Delta=.75, representative absolute errors are:

| abs(b)/Delta | C2b small core | exact finite core | regular N=100 |
|---:|---:|---:|---:|
| 1e-6 | 2.08806e-6 | within rounding | 1.04370e4 |
| 1e-4 | 2.08806e-4 | within rounding | 1.01157e2 |
| .01 | .0208799 | within rounding | .0122727 |
| .1 | .208114 | within rounding | 3.41146e-6 |
| 1 | 1.63996 | within rounding | 8.70026e-6 |
| 10 | 3.07180 | within rounding | 3.41153e-8 |

Both b signs pass. The linear cases also agree with their independent
antiderivatives to the stated 2e-12 acceptance tolerance. The table distinguishes
finite-window bias from regular under-resolution; decreasing Tlim to force the
regular rule is not an accuracy remedy. Curvature remains important: in the
Gaussian example at abs(b)/Delta=1e-6, exact integration of the single estimated
linearization has error .16413, versus .15871 for the old approximation. The
independent Gaussian reference subtracts its analytic linear part and resolves
the smooth remainder; 100,000/200,000-panel results agree within 5e-10. Thus
errors in an under-resolved reference are not mistaken for candidate errors.

For the fixed-cell transform of (1+.3i)*exp(-p^2) on [-2,2] at a=.23+.2i,
independent resolved Simpson comparison gives:

| parallel cells | absolute transform error | reduction on halving h |
|---:|---:|---:|
| 16 | .0214160 | — |
| 32 | .00540488 | 3.96235 |
| 64 | .00135389 | 3.99211 |
| 128 | .000338640 | 3.99801 |
| 256 | .0000846706 | 3.99950 |

This verifies second-order convergence of the PWL representation, separately
from exact integration of that representation. For a fixed positive separation
from the real axis, a uniform numerator error is bounded in the transform by
that error times integral dp/abs(p-a); the bound deteriorates toward the axis.
Relative determinant errors can be much larger near a small analytic
determinant, even when absolute susceptibility errors decrease.

The standalone CR tests sweep six step sizes while staying in the open
half-plane. The largest best-step normalized residual over the selected points
is 1.28e-8. Truncation dominates coarse perturbations near endpoints; rounding
dominates very small perturbations. This numerical evidence accompanies the
analytic derivation, rather than replacing it. Near-node and support-endpoint
affine values remain accurate within 2e-12 relative/absolute scaling at
abs(b)=1e-100. Far-field series, PV node subtraction and both logarithm branches
pass with bounds and invalid/zero/overflow checks in both builds.

For physical grid refinement the generator's four-width support is fixed,
but the usable domain p_1..p_(npar-1) approaches that support as h decreases.
The ordinary finite-difference derivatives and perpendicular quadrature also
change with resolution. Physical convergence therefore includes all these
finite-grid effects; the standalone fixed-support Gaussian test isolates the
parallel PWL interpolation/integration effect.

## Physical measurements

Measurements are generated in isolated builds; the final tables identify actual
accepted windings and failures without inferring a physical inventory merely
from an integer or solver residual. The reference uses numerical protons,
beta_parallel=3, anisotropy=2, analytic isotropic electrons beta=3, no drift,
physical electron mass, k=(.001,.3), vA=1e-4. Supports, normal derivatives,
method-2 fits, P=3/N=100/Tlim=.01 and all A1 controls remain at the C2b values.

### Matched original contour

The same 257 ordered frequencies are used for every method and grid. A and C2b
are byte-identical here: gamma>=.01 places all numerical proton evaluations on
the regular branch, so the proposed small-core change is inactive.

| grid | C2b/A median/max relative D error | B median/max relative D error | C2b/A max cumulative phase error | B max phase error |
|---|---|---|---:|---:|
| 32×64 | .0101841 / 1.25825 | .0106446 / 1.66419 | 6.49816 | 6.52012 |
| 64×128 | .00430580 / .137389 | .00269395 / .397650 | .0582262 | .0646049 |
| 96×192 | .00267309 / .0329546 | .00120016 / .174624 | .0269974 | .0290739 |

Phase errors are in radians, accumulated from the common starting point using
unchanged A1 phase accumulation. The fixed paths give phase/(2*pi)=3 at 32×64
and 2 on the finer grids for both methods. Independent adaptive B calculations
accept 3,2,2 in three passes, with respectively 620,628,628 collective calls.
These are computed finite-grid results, not an imposed expectation of two.
A's adaptive original-contour results coincide with C2b.

B has smaller median error on the finer grids but larger worst relative error
and slightly larger phase error. It is not uniformly more accurate at fixed
resolution. The worst sampled relative errors occur near a small analytic
determinant; absolute values are retained in `contours.dat`. A sampled maximum
below one is evidence, not a proof of a uniform boundary error bound.

For comparison, preserved pre-C2b integration gives median/max errors
.00924375/1.16099, .00408394/.212333 and .00234069/.0650747 on the three grids.
Its apparently favorable coarse winding cannot erase the independently
established discontinuities. Accuracy, continuity and analyticity are assessed
separately.

### Lower contour and grid sensitivity

All methods are also evaluated at the same lower-contour coordinates, with
bottom gamma=.001. Only the bottom samples/endpoints are lowered; vertical
interior samples retain their original positions. This is a nonuniform sampling
of the specified rectangle. Independent adaptive calculations reconstruct the
actual rectangle from its bounds and are separate from this fixed pass.

| grid | C2b median/max relative D error | A median/max | B median/max |
|---|---|---|---|
| 32×64 | .00962039 / 826.124 | .00946926 / 825.654 | .0106446 / 192.900 |
| 64×128 | .00456227 / 117.111 | .00398955 / 117.001 | .00269395 / 46.1587 |
| 96×192 | .00285410 / 35.5954 | .00219951 / 35.5457 | .00120016 / 20.0862 |

B reduces the large near-origin value error but does not eliminate it at these
grid spacings. Its fixed-contour phase differs from the analytic reference by
an additional 2*pi; the coarse original W=3 baseline is retained. A lower
integer closer to the analytic reference is not enough to prefer C2b/A when
their computed boundary curves still have switching/analyticity defects and
even larger value errors. The extra B winding requires separate physical/grid
interpretation rather than automatic identification as an instability.

### Ordinary roots and continuation limits

The independent all-analytic ordinary solver reproduces
+.39666286088749+.08598637345950i and, for the isotropic control,
+.40460138158784-.00185430537032i; negative-frequency partners are solved
separately. Numerical results for the positive-frequency partner are:

| grid | C2b growing root (A identical) | B growing root |
|---|---|---|
| 32×64 | .390581513+.086311785i | .390048496+.086275453i |
| 64×128 | .395411725+.086220210i | .394997861+.086052448i |
| 96×192 | .396168901+.086139000i | .395922144+.086015080i |

A is inactive at these growing roots. B's growth-rate error is smaller but its
real-frequency error is larger; both its real and imaginary errors decrease
under refinement. A successful solver residual is not used to claim continuum
accuracy. All signed partner residuals pass the 2e-12 check and recover the
expected drift-free reflection symmetry without imposing it in the algorithm.

| grid | C2b isotropic damped root | A damped root | B plus unchanged fitted residue (diagnostic) |
|---|---|---|---|
| 32×64 | .407435771-.001312181i | .407436773-.001322050i | .408093876-.001300822i |
| 64×128 | .405337298-.001709482i | .405337551-.001731405i | .405459045-.001720366i |
| 96×192 | .404931367-.001772382i | .404931763-.001805366i | .404986947-.001797267i |

A improves the measured damping-rate error, with little real-frequency change,
over these grids. B's fitted-combination roots also approach the analytic
reference over the tested range, but are **not accepted as validation of a
physical damped continuation**. That requires independent matching and
analyticity of the continuation, not just small determinant residuals.
Continuation fits, derivative methods and Landau formulas are identical in
these comparisons. The preserved pre-C2b root baselines remain in the C2b
report and `tests/data/nyquist_c2b/legacy_roots.txt`.

## Additional consistency gate: fixed-cell term versus unchanged int_ee

The exact cell kernel passes its integration identities, but a drop-in change
of only the resonant integral need not preserve algebraic cancellations with
other, unchanged tensor terms. `test_static.f90` demonstrates this using the
**actual** C2b nodal numerator constructor and `int_ee`, with kperp=0, B0=1,
n=0, df_perp=0 and p_perp*df_parallel=p. These are smooth prescribed derivative
arrays for an identity test, not a new physical VDF fixture.

At omega=0 the parallel Cauchy numerator is -p^3. The accompanying ee integrand
is +p^2, so their continuum sum vanishes. Applying the common original nodal
rule to both terms also cancels within rounding. Integrating PWL(-p^3)/(p-a)
exactly while leaving ee on its original nodal rule does not cancel:

| parallel intervals | h | common nodal-rule residual | fixed-PWL plus unchanged ee | residual/h^2 |
|---:|---:|---:|---:|---:|
| 64 | .03125 | 0 | -.00064181436 | -.657218 |
| 128 | .015625 | 0 | -.00016160909 | -.661951 |
| 192 | .0104167 | 4.11e-15 | -.00007199703 | -.663525 |

Both optimized and debug/FP-checked runs reproduce this **compatibility
failure**. A zero exit for this diagnostic means the predicted defect was
reproduced; it is not acceptance of B as a complete susceptibility replacement.
The result is not an error in the logarithm formula: it analytically integrates
its prescribed PWL numerator. It is an inconsistency between representations
of terms that cancel before discretization. The full C2b moving-core method
has its own local approximation error; the exact cancellation column refers
to the common nodal rule, not a claim that the complete C2b method is exact.

A useful mathematical deduction at kperp=0 is obtained by writing
S=p_perp*df_parallel-p*df_perp. Before discretization,

    q*p^2*(omega*df_perp + kpar*S/m)/(m*omega-kpar*p) + q*p*S/m
      = q*omega*p*p_perp*df_parallel/(m*omega-kpar*p).

Combining the terms first preserves an explicit omega factor. Interpolating
and integrating them independently need not preserve that identity. This
suggests a structure-preserving fixed-cell formulation as follow-up work;
it is **not implemented here**. At finite kperp the other harmonics and tensor
identities must be treated coherently. Adjusting an unrelated constant merely
to recover an expected winding would not be a justified correction.

The physical near-origin data and optional root-location probes are retained
to assess the connection with B's extra low-growth zero. The polynomial test
establishes the compatibility defect; attributing every physical discrepancy
solely to it requires further tensor-level analysis. All production ee,
susceptibility, differentiation and continuation formulas remain unchanged.

## Signed-k and MPI controls

Both signs of kpar are evaluated independently. For the 64×128 original fixed
contour the maximum relative difference between kpar=+.3 and -.3 is
6.13927e-5 for C2b/A, versus 1.13004e-14 for B. The largest C2b absolute
difference is 2.16833e-4 at omega=.90625+.01i. This is an inherited support-edge
limitation, not a failed debug/optimized comparison or a new finite-core defect.

For n=-1 there, the resonance momenta are ±6.3541666667, in parallel cells
122 and 5. With P=3, the existing lower guard 5-P<=2 takes its special fallback,
while the reflected upper guard 122+P>=npar-2 is false. Thus the original
edge-handling logic is not reflection-symmetric. B's complete fixed-domain
integral has no such branch. C2b explicitly left support-edge returns unchanged;
this study records the consequence and does not alter the production guards.

For all three grids and both matched contours, exact metadata comparisons show
identical per-rank harmonic/component/frequency call order for legacy, C2b, A
and B. Every actual determinant call remains collective; independent per-rank
counters agree. The startup assignment depends normally on the species/support
(e.g. a stable control can have a different harmonic maximum); no assignment
routine, communicator or work-scheduling policy changes. Landau comparisons
also include negative-k causal probes rather than assuming sign(b)=sign(gamma).

## Branch matching and frequency differentiability

At the default switch gamma=.003, pairs separated by 2e-12 give these relative
D differences (ranges over a grid node, quarter-cell and half-cell near Re=.4):

| grid | C2b | A finite window | B fixed cells |
|---|---|---|---|
| 32×64 | .02735–.03169 | .02702–.03477 | 2.65–2.84e-11 |
| 64×128 | .01970–.02353 | .00356–.00527 | 2.80–2.92e-11 |
| 96×192 | .03041–.03427 | .00177–.00246 | 2.88–2.92e-11 |

The A jump plateaus as separation decreases; B's differences scale with
frequency separation. A removes most of the fine-grid switch error, but is not
an exact match between branches, and can worsen the coarse-grid jump.

At the 64×128 node probe omega=.389711431703+.003i, the n=0, zz absolute
finite-core/regular discrepancy is .00298134 for N=100, .00298187 for N=400,
and .00298190 for N=1600. The old-small/regular differences are .0179849,
.0179867 and .0179868. The full determinant's relative finite/regular differences
are .00527756, .00527511 and .00527503. Increasing regular quadrature density
does not remove this residual: the single linearization differs from the PWL
numerator over the full core.

At gamma=.0003 the same finite/regular component discrepancy changes from
.286748 (N=100) to .00282752 (400) and .00277384 (1600). Here under-resolution
of the narrow peak dominates N=100. All three momentum grids have N=100/400/1600
comparisons and independent analytic determinant data in the retained supplement.
Curvature, finite-window residue excess and regular quadrature resolution are
therefore distinguished rather than conflated with Tlim selection.

For independent real/imaginary perturbations, the largest best-step normalized
CR residuals in the collective B determinant are 3.35e-9, 2.97e-9 and 3.45e-9
on the three grids. The negative-k 64×128 result is 2.97e-9. At the same 64×128
points C2b and A retain residuals as large as .5267 and .1151, including points
inside the original real-frequency interval. Those observational controls are
not marked as passing analyticity. Perturbations remain within a fixed Tlim
branch, so its switch alone cannot explain the nonzero CR plateau.

## Grid residue versus fitted Landau residue

At Re(omega)=.4046 and imaginary displacements ±1e-12 in the isotropic control,
the actual direct B integral satisfies the grid Plemelj jump. Adding a matching
grid residue would close that boundary jump to rounding/displacement accuracy.
The unchanged fitted ALPS residue does not exactly equal that grid residue:

| grid | n=0, zz direct-jump error | remaining fitted causal jump | relative fitted/grid residue difference |
|---|---:|---:|---:|
| 32×64 | 3.77e-12 | .0276019 | .0153941 |
| 64×128 | 3.85e-12 | .00800260 | .00440626 |
| 96×192 | 3.98e-12 | .00332440 | .00182506 |

All 171 recorded harmonic/component direct boundary values are checked on each
grid. Actual Landau streams are byte-identical to C2b at the same frequencies;
positive and negative-k controls are also retained. The fitted/grid differences
decrease but remain finite at these resolutions. They reflect the different
representations, derivatives and perpendicular endpoint treatments already
identified in C2b; this task does not change them.

A grid-based real-axis matching diagnostic is not a damped continuation rule.
The Gaussian PWL counterexample at a=-.01i gives neighboring complex cell
extensions differing by about .00130, although their values match at the real
shared node. Selecting between them using Re(a) would introduce an off-axis
value jump. No such selector is used in the candidates. A consistent analytic
representation and its continuation/branch structure require separate review.

## Near-origin physical consistency and additional computed zeros

At omega=i*1e-8, the analytic determinant is -1.54618e-16. The numerical
near-origin limits are already effectively constant over the smallest probe
heights; the unchanged ee term is verified identical across all methods:

| grid | C2b Draw | A Draw | B Draw | B proton scaled chi_zz |
|---|---:|---:|---:|---:|
| 32×64 | -.00127059 | -.00126993 | .000303125 | .00290990 |
| 64×128 | -.000179660 | -.000179415 | .0000726613 | .000693422 |
| 96×192 | -.0000544587 | -.0000543680 | .0000316270 | .000300749 |

At 64×128, the B-minus-C2b chi_zz change is .0024207555, while the n=0 integral
change at the same unchanged ee is .0024207545. At 96×192 those changes are
.00082567066 and .00082567071. The other harmonics account for a very small
part of that difference. This supports the static-compatibility diagnosis;
a full finite-kperp tensor-identity audit is still required before claiming
that one isolated term explains every determinant error.

Independent ordinary-solver probes locate the additional B zeros at:

The optional probes use guesses 1e-6+i*(.02,.005,.002) and D_threshold=1e-17
with the existing ordinary solver. This tighter residual target is confined to
the near-origin diagnostic; the established physical-root controls stay at
1e-13, and no Nyquist controls are changed.

| grid | real frequency | growth | raw determinant residual |
|---|---:|---:|---:|
| 32×64 | 5.18e-18 | .0126394925346 | 1.49e-18 |
| 64×128 | 2.41e-17 | .00647826428106 | 3.17e-18 |
| 96×192 | 1.13e-17 | .00435308676425 | 4.49e-18 |

Their growth divided by parallel spacing is .05838, .05984 and .06032:
the extra zero approaches the real axis approximately linearly with h.
This explains why B's original count changes from 3 to 2 under refinement,
while its lower count remains 3 on these grids. Both analytic rectangles
independently accept W=2 (630/690 collective calls); their two known growing
roots remain enclosed. The extra zeros are therefore grid-model discrepancies,
not newly established physical instabilities or missed-loop failures in A1.
Neither numerical method has demonstrated physical convergence of the very
low-frequency determinant on the tested grids. Counts alone would hide that.

## Performance and implementation cost

Repeated four-rank timings use identical frequencies, full-stream diagnostic I/O
disabled, one BLAS/OpenMP thread, and the final counters including background
resU numerator constructions. They exclude equilibrium preparation. Means and
ranges are from three runs; matched point output remains enabled outside the
disp timer.

| path | method | calls | mean collective seconds (range) | numerator constructions | fixed-cell integrals |
|---|---|---:|---|---:|---:|
| original contour | C2b | 257 | 10.2234 (10.1362–10.3893) | 489,256,551 | 0 |
| original contour | A | 257 | 10.2406 (10.1749–10.3522) | 489,256,551 | 0 |
| original contour | B | 257 | 72.4350 (72.2049–72.6781) | 376,295,031 | 373,332,078 |
| bottom gamma=.001 | C2b | 65 | 1.62443 (1.57485–1.68535) | 96,470,730 | 0 |
| bottom gamma=.001 | A | 65 | 1.62700 (1.60156–1.64056) | 95,451,426 | 0 |
| bottom gamma=.001 | B | 65 | 18.2979 (18.2329–18.3641) | 95,171,895 | 94,422,510 |

The bottom rows are raw fixed-frequency timing diagnostics, not horizontal
Nyquist counts. B is about 7.09 times slower on the original contour and 11.26
times slower on this bottom-edge workload. Fewer numerator constructions do
not produce a speedup here. A's small reduction on the bottom edge gives no
resolved collective timing improvement; the existing far-harmonic work and MPI
waiting remain. No resonance reassignment or alternative parallelism is used.

B constructs 1,464,183 nodal numerators and evaluates 1,452,654 cells per
64×128 determinant in this equilibrium. The implementation recomputes each
kernel for each perpendicular node and tensor component. Each cell uses two
complex logs or a conservative 40-term series; there is no kernel reuse. These
are **instrumented prototype costs**, not an intrinsic lower bound for the
fixed-cell method. No performance optimization was implemented to improve the
comparison artificially.

Joint grid-size measurements with full diagnostics enabled give B about
.0739, .2892 and .6500 seconds per determinant at 32×64, 64×128 and 96×192.
Those follow the expected O(nperp*npar*harmonics*components) work; npar and
nperp were not varied independently for timing. Adaptive original-contour B
times were 45.83,181.62,408.23 seconds for 620,628,628 calls. The corresponding
C2b times were 8.73,27.76,55.13 seconds. Scalar/core controls, full tables and
individual repeats remain available rather than a single favorable timing.

The fixed-cell prototype uses automatic p and G work arrays of size npar-1,
about 24*(npar-1) bytes per active integral, plus scalar temporaries. Its current
code has no per-cell allocation. A uses scalar work only. Reusing exact cell
weights per species/harmonic/frequency could remove redundant log/series work
without changing MPI decomposition, but should follow a coherent tensor
representation and independent validation. No general numerical-VDF performance
claim follows from this near-parallel Gaussian benchmark.

## Validation inventory and reproduction

The completed C2b prerequisite includes the six Nyquist entry points
`run_nyquist`, `run_nyquist_growth`, `run_nyquist_log`, `run_nyquist_mirror`,
`run_nyquist_scan`, `run_nyquist_plane`, plus C1 and C2b, in optimized and debug
builds. The existing analytic scientific fixtures, including nine ordinary
root/map/line/plane outputs, remain byte-identical. Logs/statuses are in
`/tmp/alps-nyquist-c2b/regressions/`. No production source changes during C2c;
the source hashes in the C2c common-build provenance are checked unchanged.

C2c execution coverage:

- Constant/linear/polynomial/Gaussian, signed-b/k, PV/endpoint, far-field,
  logarithm-branch, CR and PWL-convergence tests pass in optimized and debug.
  The actual C2b core agrees with its preserved scalar legacy helper.
- The static identity diagnostic reproduces a whole-discretization failure in
  both builds. It is a reason to reject the naive B replacement, not a passing
  physical-accuracy gate.
- The complete three-grid physical matrix runs optimized: 18 independent
  numerical adaptive contours, both analytic reference contours, matched
  original/lower trajectories, signed ordinary roots and all core controls.
- The 64×128 debug physical gate repeats matched contours, forced branches,
  causal/real-axis/signed-k probes and roots. All 3,064 values in 26 matched
  datasets and six two-root files agree byte-for-byte with optimized output.
  The full three-grid C2c adaptive physical matrix was not repeated in debug.
- Boundary and near-origin analyzers also run with debug bounds and
  invalid/zero/overflow traps on the retained physical streams. Their tables
  agree with optimized analyzers. All actual disp-call counts agree across ranks.

GNU Fortran 11.4.0/Open MPI 4.1.2 are the identified tools. Physical optimized
flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags
are `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`. Arithmetic standalones
add invalid/zero/overflow traps in both configurations. Harmless underflow of
negligible far-field series terms is not trapped. No new numerical dependency
or Python is introduced.

The primary runner was edited while active, causing a parse interruption after
its complete three-grid matrix. The remaining controls were run from a frozen
`resume_footer.sh`; completed scientific data were retained. Optional near-origin
root inputs initially hit the existing line-oriented namelist reader and were
rerun in its normal multiline format. Neither recovery changed physics,
thresholds, bounds, guesses or numerical tolerances. The corrected reproducible
runners pass shell syntax checks; the retained logs include both incidents.

Identified artifacts live under `/tmp/alps-nyquist-c2c/`:
`study_optimized`, `study_debug`, `supplement`, `math_optimized`, `math_debug`,
`static_optimized`, `static_debug`, `build_agreement`, `tables_complete` and their
logs. Each physical build retains copied source and executable hashes, inputs,
matched values, per-rank streams and metrics. Compact measured data are kept
under `tests/data/nyquist_c2c/`; the large tensor streams remain external.

From the checkout, reproduce the complete primary/supplementary study with:

```sh
export C2C_SCRATCH=$(mktemp -d /tmp/alps-c2c.XXXXXXXX)
bash validation/run_nyquist_c2c.sh
bash validation/c2c/extras.sh "$C2C_SCRATCH" "$C2C_SCRATCH/supplement"
bash validation/c2c/report_tables.sh "$C2C_SCRATCH" "$C2C_SCRATCH/supplement" "$C2C_SCRATCH/tables"
```

Use `C2C_DEBUG_SMOKE=T` and the debug FFLAGS to reproduce the representative
physical debug gate; omit that variable for the longer full matrix. The scripts
reuse normal VDF generation, differentiation and fitting. `compare_debug.sh`
compares completed optimized/debug trees. Production defaults, inputs, winding,
growth, scans, susceptibility formulas and Landau continuation remain untouched.

## Recommended next step

**Choose C: continue mathematical development.** The finite-window A identity
is correct and removes its targeted error, but A retains finite switch jumps
and significant CR residuals. It is insufficient as the general Nyquist
susceptibility remedy.

B is the preferred integration kernel for further work: exact cell integration,
fixed-domain upper-half-plane analyticity, stable signed branches and restored
support-edge/signed-k behavior are established. Its present drop-in application
is not ready for production: the static cancellation fails, the low-growth
winding has not converged to the analytic reference, the fitted causal boundary
mismatch remains, and the unoptimized cost is substantially higher.

The next controlled prototype should discretize the coupled Cauchy and
nonresonant tensor terms coherently, preserving their algebraic moment/static
identities before assessing weak-growth counts. Test the parallel-limit
identity first, then the finite-kperp tensor and the existing physical boxes.
Treat damped continuation with a consistent analytic representation as a
separate requirement. Reuse of exact kernel weights can be benchmarked after
those identities pass. Do not alter Tlim or frequency bounds to conceal the
additional zeros, and do not introduce grid-cell selection into complex
momentum continuation. Numerical-electron studies, growth searches and
wavevector scans remain deferred. This study authorizes no default change.
