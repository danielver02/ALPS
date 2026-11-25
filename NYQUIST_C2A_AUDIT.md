# C2a — numerical susceptibility discontinuity audit

This audit identifies **two independent half-cell value mismatches**, and a third
remainder-length defect. The floor/cell lookup is not responsible for the
identified half-cell jump: the selected cell does not change there.

1. `funct_g` switches between nearest-node, centered-linear Taylor reconstructions
   whose values generally do not match at their common midpoint.
2. At that same boundary, one whole parallel-momentum cell moves from background
   trapezoidal integration to a subdivided local integral. Those two numerical
   integrals are different, even for a constant numerator.
3. A positive remainder shorter than one nominal substep is omitted when
   `ntiny=0`; this creates additional value jumps away from the midpoint.

The `Tlim` discontinuity is separate. The small-imaginary formula uses an
infinite-limit residue for a **finite** symmetric window and linearizes the
numerator. It therefore does not equal the regular calculation at finite
imaginary resonance displacement. Making `Tlim` extremely small can instead
expose under-resolution of the regular quadrature's narrow resonance peak.

**No production source, susceptibility formula, default, interpolation,
quadrature, continuation, MPI assignment, determinant, root search, scan, A1 or
A3 implementation is changed.** All instrumentation and forced-branch tests are
in copied validation sources. A lower `Tlim` is an explicit diagnostic input,
not a production correction. The numerical-electron result remains winding 3;
no physical root inventory is inferred from that disagreement.

## 1. Reproduction, sources and execution contract

Read with `NYQUIST_C1_RESULTS.md` and `NYQUIST_A0_AUDIT.md`. The equilibrium is
unchanged: numerical protons with beta_parallel=3, alpha=2; isotropic analytic
electrons with beta=3, mass=0.000544617021487, densities 1, charges ±1, no drift,
k=(0.001,0.3), vA/c=1e-4. Proton support is four Gaussian e-folding widths,
0<=p_perp<=4 sqrt(6), |p_parallel|<=4 sqrt(3). The 32×64, 64×128 and 96×192
figures refer to intervals, with endpoints included. The normal ALPS generator,
`read_f0`, numerical derivatives, method-2 log-Chebyshev continuation (order 8),
model transfer, harmonic assignment and Bessel initialization are retained.
P=`positions_principal`=3; N=`n_resonance_interval`=100 unless a test explicitly
specifies 400/1600. Bessel_zero=1e-45; all A1 controls remain at C1 defaults.

The first action, before adding susceptibility instrumentation, reproduced the
C1 limiting jumps with the existing C1 fixed-frequency adapter and original
numerical routines. Across 2e-12 frequency separations:

| Probe | relative determinant jump |
|---|---:|
| Half-cell, gamma=.001 | 2.3877962912 |
| Half-cell, gamma=.01 | .009428216625 |
| Tlim=.01, omega=.4, gamma=.003 | .02055367975 |

The analytic determinants vary smoothly (relative changes ~1e-10 or smaller).
Subsequent instrumented default-mode determinants are **byte-identical** to
these six original complex values. The reproduction is also the first gate in
`bash validation/run_nyquist_c2a.sh`, before compiling the instrumented routines.

Rank 0 retains the existing callback/request protocol. All four ranks execute
every physical `disp` call through the unchanged reductions, determinant
broadcast and barrier. Instrumentation writes per-rank files; it introduces no
MPI calls, communicator, reassignment or sample parallelism. Per-rank determinant
counts are checked by the existing driver and validation harness. Full harmonic
records check unique ownership and reproduce changes in the assembled species
tensor. Ordinary full-contour tests use the **unmodified executable**, not the
forced-branch copy.

Source preservation is checked against `/tmp/alps-nyquist-c2a/snapshot/`.
Relevant unchanged SHA-256 values:

| Source | SHA-256 |
|---|---|
| `src/ALPS_fns.f90` | `567e05bad770489d578c0024b4230376f1e000c2f04b7e84fc2f7d065d446be5` |
| `src/ALPS_com.f90` | `989cfbb8fd3a8252e195d868a93a09552cfbce3e4f1172f93c6a127e947e5cc8` |
| `src/ALPS_nyquist.f90` | `d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31` |
| `src/ALPS_nyquist_growth.f90` | `5eda3a5b6faf58ca54417b5bf26f4bd89abe6a2e47d1f8ac9978a73434d1a212` |

## 2. Actual integration algorithm and switch locations

All line numbers below refer to the unchanged production `src/ALPS_fns.f90`.
Let a=r+ib=(m_s z-n q_s)/k_parallel, h=delta_p_parallel, and j be the grid node
immediately to the left of r. For fixed p_perp and tensor mode define

    G(p;z) = -(q_s/k_parallel) * [z df/dp_perp
                  + (k_parallel/m_s)(p_perp df/dp_parallel - p df/dp_perp)] * T_n(p).

Then `resU*int_T = G(p;z)/(p-a)` at grid nodes. For b>0 there is no real-momentum
pole; with a smooth numerator and finite support this integral changes
continuously as a crosses any real grid node or midpoint. A valid change of
partition cannot change its value. Derivative kinks in an interpolant are a
different issue from finite one-sided value differences.

| Routine / lines | Role and discrete choices |
|---|---|
| `disp`, 369–493, 497–520 | Numerical harmonic loop; modes 1/2/3/4/5/6 are xx/yy/zz/xy/xz/yz. n=0 uses yy/zz/yz only. `int_ee` adds a frequency-independent zz term. Species multiplier is n_s q_s; original MPI reductions are unchanged. |
| `determine_resonances`, 626–730 | Tests **real** resonance positions against the parallel support and an extended edge band. The identified central half-cell pair stays `found_res=T`. |
| `full_integrate`, 735–777 | Direct fixed-grid integration if no resonance is found; `integrate_res` for the strictly upper-half-plane nonrelativistic response. No Landau continuation residue is added in these tests. |
| `integrate`, 784–849 | Tensor-product grid trapezoids of `resU*int_T`, including existing doubled first perpendicular endpoint weight. Its bounds are integer grid indices. |
| `integrate_res`, 957–974 | `floor((real(p_res)-p0)/h)` with neighboring-node corrections identifies j. **At the identified half-cell j=63 on both sides.** |
| `integrate_res`, 995–1004 | Special edge returns integrate only a truncated remote portion. These are not entered by the n=0 central pairs and do not explain their jumps. Edge behavior remains a separate accuracy concern. |
| `integrate_res`, 1010–1020 | Left background ends at j-P. Right background starts at j+P+1 for r-p_j<h/2 and j+P+2 otherwise. This is a real half-cell switch. |
| `integrate_res`, 1031–1039 | Delta=r-p_(j-P), local step=Delta/N; branch is `abs(denomI).GT.Tlim`. Equality goes to the small-imaginary branch. |
| `integrate_res`, 1041–1100 | Regular paired integral over [r-Delta,r+Delta], using `funct_g` on moving subnodes. |
| `integrate_res`, 1105–1180 | Linearized small-imaginary formula, centered finite-difference G_prime, and sign(b)*i*pi*G(r) limiting residue. |
| `integrate_res`, 1188–1247 | Positive-side remainder to the selected upper grid node, `ntiny=int(R/step)`, integrated **only if ntiny>0**. Its substep is adjusted to land on the endpoint. |
| `integrate_res`, 1250–1251 | Perpendicular/parallel quadrature normalization, then addition of the background integral. |
| `funct_g`, 1320–1335 | Predicts nearest node using `nint`, reproduces the old half-cell distance test and clamps boundary stencils. |
| `funct_g`, 1343–1362 | Returns G_j+(G_(j+1)-G_(j-1))/(2h)*(p-p_j), a nearest-node Taylor line rather than a common cell interpolant. |
| `resU`, 1653–1689; `int_T`, 1694–1793 | Original resonant denominator and Bessel tensor factors. No selector changes in these functions explain the central half-cell event. |

### Intended moving partition

Write r=p_j+alpha h, 0<=alpha<1. The symmetric window has half-width
Delta=(P+alpha)h. Its left end p_(j-P) is fixed within the cell; its right end is
p_j+(P+2alpha)h. The positive remainder length is

    R = (1-2 alpha)h    for alpha<1/2,
        (2-2 alpha)h    for alpha>=1/2.

As alpha crosses 1/2, R jumps from 0 to h, while the right background loses
exactly the cell [p_(j+P+1),p_(j+P+2)]. **The geometric partition is correct in
this limiting sense.** No central interval is off by one or duplicated merely
because `upperlimit` changes. Numerical continuity requires the added local-cell
integral to equal the lost background-cell integral, and the numerator sampled
by the local rule to be continuous. Neither condition is currently satisfied.

### Repository history

`git blame` places the half-cell upper-limit logic in commit `6340a8c4`
(2020-05-23). The current centered-linear reconstruction is already present in
`eef5274e` (2021-06-28). The April 2025 simplifications (`567f659e`, `006b3ae4`)
retain the finite-window/remaining-domain structure; the latter changes only the
order of an addition in an endpoint expression. Commit `8d49f15a` (2026-10-13)
explicitly preserves the old search's midpoint selection and stencil boundaries
while replacing its global scan with predicted indexing. It does not establish
interpolant continuity. These findings are not evidence that the newer index
prediction introduced the value mismatch. The comments describe symmetric
window coverage and a remaining interval, not an equality proof for the two
quadrature rules.

## 3. Half-cell interpolation: a finite value mismatch

At p_mid=(p_j+p_(j+1))/2, the two one-sided reconstructions are

    G_L = G_j     + (G_(j+1)-G_(j-1))/4,
    G_R = G_(j+1) - (G_(j+2)-G_j)/4.

Therefore, exactly on a uniform grid,

    G_R-G_L = -(G_(j+2)-3G_(j+1)+3G_j-G_(j-1))/4.

This is minus one quarter of the third finite difference. It vanishes for a
quadratic numerator, but not for a generic smooth numerator. For a smooth
function the leading jump is -h^3 G'''/4. This is a **value** discontinuity,
not merely differing slopes of a continuous piecewise-linear approximation.

A standalone test calls the actual copied `funct_g`, embedding prescribed
numerators in its input arrays. It changes no physical VDF test derivatives.
With h=.125, G=p^3 and midpoint -h/2, the one-sided limits are
+0.001953125 and -0.0009765625: jump -0.0029296875=-1.5 h^3. Constants and linear
numerators are continuous. Quadratic midpoint values agree but are biased: for
G=p^2 both reconstructions give zero at this midpoint, where G=h^2/4.

In the real proton case, at perpendicular index 8 of the 64 grid
(p_perp=1.224744871...), n=0, mode 3, gamma=.001, the actual selected stencils
switch 63->64. The limiting reconstructed numerators are approximately

    G_L =  6.710923444e-6 + 2.4825765e-9 i,
    G_R = -3.455721280e-6 + 0 i.

The recorded actual value agrees with the appropriate one-sided interpolant.
This mismatch feeds both the limiting residue and the finite-difference slope
used in the small-imaginary calculation. At gamma=.001 its residue jump is the
largest imaginary part of the n=0 zz core jump.

## 4. Half-cell quadrature transfer: independent of interpolation

A second standalone test prescribes G=1, for which `funct_g` is exact and has no
stencil jump. Nevertheless the full original `integrate_res` jumps at the
half-cell. With h=.125, b=.01, P=3 and N=100, after dividing out the known
perpendicular quadrature factor, the observed jump is

    -0.002643860880 - 0.000160817582 i.

This isolates a partition/quadrature mismatch without a numerator interpolation
error. Background integration uses one cell trapezoid; the replacement uses
subdivided quadrature of G/(p-a), including a nonpolynomial resonant denominator.
They need not agree even for constant G. A vanishing remainder on one side does
not make the added full cell on the other side match the removed trapezoid.

### Real tensor decomposition

The following values are frequency-scaled proton tensor contributions, before
any conventional susceptibility division by (z vA)^2. For protons n_s q_s=1, so
the integral changes directly equal changes in the corresponding raw species
tensor entry. Dividing by the smooth nonzero frequency factor does not remove
the discontinuity. `int_ee` is constant and contributes no jump.

At the 64-grid midpoint, j=63, lowerlimit=60, upperlimit=67->68, ntiny=0->28.
Neither the cell index, resonance-found flag, species nor harmonic ownership
changes. Let Q_cell be the background trapezoid of the transferred cell.

| n=0, zz term | gamma=.001 | gamma=.01 |
|---|---|---|
| Delta(background) | -1.222952200e-3 -2.152111054e-5 i | -1.206734032e-3 -2.139495390e-4 i |
| Delta(tiny remainder) | +1.181979204e-3 +2.131206769e-5 i | +1.165861332e-3 +2.118755379e-4 i |
| Replacement minus Q_cell | -4.097299592e-5 -2.090428744e-7 i | -4.087270011e-5 -2.074001340e-6 i |
| Delta(symmetric core) | -6.424783847e-5 -6.288579519e-4 i | +5.583762087e-8 -2.277686550e-5 i |
| Delta(core residue alone) | +1.543287774e-7 -6.295263151e-4 i | 0 (regular branch) |
| Delta(total n=0 zz integral) | -1.052208347e-4 -6.290669947e-4 i | -4.081686284e-5 -2.485086662e-5 i |

Delta(background)+Q_cell is only ~4e-13, consistent with the vanishing physical
probe displacement. Thus the lost background cell is identified directly.
The replacement's mismatch and the core's interpolation mismatch both remain.
There is no missing integer cell in the **geometric** half-cell partition.

The determinant changes are -1.344275420e-5-6.471005232e-5 i and
-4.237590395e-6-2.373368266e-6 i respectively. The n=0 zz change alone, multiplied
by its determinant cofactor, accounts for the full determinant jump to about
1e-6 relative or better. Other components give small corrections. Full harmonic
sums reproduce the assembled species-tensor changes to ~1e-17 in this probe;
background+core+remainder reproduce each traced integral to roundoff.

### A further omission: ntiny=0 for a positive remainder

Lines 1188–1191 omit the remainder whenever 0<R<step. At
alpha=(N-P)/(2N+1), ntiny changes 1->0 although the limiting remainder length is
positive. The constant-numerator test isolates this too: h=.125, b=.01, N=100,
P=3, r=-0.0646766169154 produces the normalized jump

    -0.009945298493 - 0.000227332356 i.

This is separate from the exact half-cell limit, where the disappearing side's
remainder genuinely tends to zero. Raising N reduces the omitted interval but
does not make a discontinuous switch continuous at a fixed N.
Other integer changes of panel count can also change a quadrature approximation;
the 1->0 omission is the one independently isolated here.

## 5. One-sided limits and momentum resolution

Probe separations are 2 delta, delta=1e-6,1e-8,1e-10,1e-12. At the 64-grid
half-cell, relative determinant changes plateau rather than vanish:

| gamma | delta=1e-6 | 1e-8 | 1e-10 | 1e-12 |
|---|---:|---:|---:|---:|
| .001 | 2.382551 | 2.387744 | 2.387796 | 2.387796 |
| .003 | .7372696 | .7375360 | .7375387 | .7375387 |
| .01 | .00949994 | .00942891 | .00942822 | .00942822 |
| .03 | .00258682 | .00257574 | .00257566 | .00257566 |

The independent analytic response changes proportionally to displacement,
reaching ~1e-10 relative variation at the closest real-frequency pair.

At each grid's own midpoint omega_c=-k_parallel h/2, delta=1e-12:

| Proton intervals | omega_c | default jump, gamma=.001 | forced-regular jump, gamma=.001 | regular jump, gamma=.01 | regular jump, gamma=.03 |
|---|---:|---:|---:|---:|---:|
| 32×64 | -.0324759526419 | .610976 | .306611 | .0249847 | .00627543 |
| 64×128 | -.0162379763210 | 2.387796 | .0715300 | .00942822 | .00257566 |
| 96×192 | -.0108253175473 | .304710 | .0308020 | .00477786 | .000987799 |

The relative low-gamma determinant jump is not monotone because its denominator
is small and grid-dependent. The **absolute n=0 zz integral jump** at gamma=.001
falls from .00499568 to .000637806 to .000191039; at gamma=.01 it falls from
.000355775 to .0000477868 to .0000149675. This is close to cubic scaling for
these midpoint probes and consistent with the smooth-numerator third-difference
mechanism, not evidence that j is randomly misidentified.

The feature moves with the grid. Probing the same C1 physical frequency
-.016237976321 on the 32/64/96 grids gives default gamma=.001 relative jumps
.000104695 / 2.387796 / .000089466. It is the true resonance midpoint only on the
64 grid. Other smaller jumps can still occur there because moving quadrature
nodes encounter other stencil boundaries.

Controls at grid nodes and quarter-cells also exhibit value changes. For
example, on the 64 grid at gamma=.01: omega=0 gives .0227370 and omega=k h/4
gives .000453481, versus analytic changes of order 1e-10. The data therefore do
not imply that **only** resonance midpoints are problematic. Cell changes alter
window placement, and local moving quadrature nodes can cross stencil boundaries
before the resonance itself does. Full per-harmonic records are retained; a
unique attribution for every secondary moving-node switch is not claimed here.

## 6. Tlim: finite-window approximation and regular peak resolution

The source gives b=m_s gamma/k_parallel, hence the switch is

    gamma_switch,s = |k_parallel| Tlim / m_s.

The per-species mass is essential; `Tlim` is not a frequency tolerance. With the
proton k_parallel=.3, one-sided delta_gamma=1e-12 gives:

| Tlim | gamma_switch | relative determinant jump | limiting D_small | limiting D_regular |
|---|---:|---:|---|---|
| .02 | .006 | .0407065 | -.02938309744 + .01509852003 i | -.02846356682 + .01411728534 i |
| .01 | .003 | .0205537 | -.03014123767 + .01527658866 i | -.02975930929 + .01469648756 i |
| .005 | .0015 | .0110436 | -.03052144074 + .01537657624 i | -.03043764750 + .01500856979 i |
| 1e-5 | 3e-6 | 30986.8 | -.03088135623 + .01546791857 i | -204.1901527 + 1050.601329 i |

The analytic variation is ~3e-11 in all four closest pairs. The last row is a
warning about **under-resolved regular integration**, not evidence that the
analytic response becomes singular there.

### A mathematical mismatch at finite b

For a constant or linear local numerator G(r+t)=G0+G1 t, the exact integral over
the symmetric finite window, b>0, is

    I_window = 2 i atan(Delta/b) G0
               + 2 [Delta-b atan(Delta/b)] G1.

The current small-imaginary code quadratures the linearized derivative part but
uses **i*pi*G0**, the b->0 residue, for the constant part. It omits the finite-b
correction on that window:

    I_small - I_window = i [pi-2 atan(Delta/b)] G0
                       = 2 i atan(b/Delta) G0,

before numerator-curvature and derivative-approximation errors. The remote
integral is still present, so it does not supply an opposite correction.
This mismatch goes to zero as b/Delta->0, but is finite at an arbitrary Tlim.
This local limiting term is distinct from `full_integrate`'s below-real-axis
Landau-continuation addition, which is not traversed in these positive-gamma tests.

The constant-numerator test compares both **original branches at the identical
frequency** with N=1600. It reproduces this predicted difference to 2.98e-9 after
removing the perpendicular factor. Thus the mismatch is not explained solely by
regular-branch quadrature error or a frequency-coordinate change.

For the actual n=0 zz integral at omega=.4, b=.01, Delta=.3590547541:

- observed small-minus-regular core: .002872378221 + .016639965440 i;
- constant-part finite-window residue excess predicted from the recorded G(r):
  -.000059357699 + .015805636561 i;
- residual after that term: .002931735920 + .000834328879 i.

The finite-window residue explains most of this component's imaginary mismatch.
The remaining difference includes numerator curvature, G_prime approximation,
and finite quadrature errors; it is not attributed uniquely to one Taylor term.
Background and remainder agree at identical frequency. The full forced-branch
files expose this comparison separately for every traced harmonic/component.

### Which components produce the Tlim determinant jump?

At the default switch omega=.4+i*.003, n=0 zz has the largest individual tensor
jump, but n=+1 transverse entries also matter. Below are the measured
regular-minus-small integral changes and their first-order determinant changes:

| Proton harmonic / component | Delta(integral) | determinant cofactor times Delta(integral) |
|---|---|---|
| 0 / zz | -.0028723782 -.0166399654 i | +.0003183044 +.0001091011 i |
| +1 / xx | +.0007998758 +.0019749778 i | +.0002905462 +.0000927632 i |
| +1 / yy | +.0007998710 +.0019749661 i | +.0002905769 +.0000927808 i |
| +1 / xy | -.0019749719 +.0007998734 i | -.0004787157 -.0008727951 i |

Other harmonics/components are retained in the full records; their vector
contributions and these entries partly cancel. The full first-order sum is
.0003868587-.0005839241 i, versus the actual determinant change
.0003819284-.0005801011 i. The relative residual is .00898, consistent with a
finite multicomponent jump's higher-order determinant terms. Magnitudes of
these individual contributions must not be added as percentages. Unlike the
half-cell case, a single zz term does not fully explain the Tlim determinant jump.

### Same-frequency branch and quadrature comparison

`branch_override=1/2` exists **only in the copied audit module**, not in an ALPS
namelist. It selects the original regular or small-imaginary body without moving
z or altering the production condition. The numerical tensor, background,
remainder and final determinant are recorded at identical coordinates.

At omega=.4, 64×128, the relative complex determinant errors against the analytic
response are:

| b | regular N=100 | regular N=400 | regular N=1600 | small N=100 |
|---|---:|---:|---:|---:|
| 1e-5 | 31221.5 | 617.478 | 25.5842 | .0113393 |
| 1e-4 | 69.6674 | 3.93875 | .168468 | .0114854 |
| .001 | .559807 | .0205896 | .0194239 | .0129948 |
| .005 | .0198567 | .0197058 | .0197203 | .0200234 |
| .01 | .0202085 | .0202092 | .0202140 | .0297281 |
| .02 | .0211488 | .0211439 | .0211457 | .0508522 |
| .05 | .0239185 | .0239127 | .0239134 | .124236 |
| .1 | .0294931 | .0294875 | .0294878 | .292105 |

At the default b=.01 the regular branch is already essentially quadrature
converged; its difference from the small-imaginary branch remains. At b=1e-5,
even 1600 substeps do not adequately resolve this peak. The constant-numerator
fixture independently reproduces the mechanism: a quadrature node is exactly at
t=r, and its regular contribution scales like step*G0/b. With a fixed substep
this diverges as b->0 although the continuous integral has a finite causal limit.
Thus the two *discrete* branches need not approach each other smoothly as b->0
unless the regular quadrature is refined appropriately.

The comparison also depends on spatial position. At b=.01 the same-frequency
|D_small-D_regular|/|D_regular| is .02093 at omega=.4, about .898/.905 on the two
sides of the C1 half-cell, and 3.494 at omega=0. These use different denominators
from the one-sided Tlim jump table. Near-small determinants amplify relative
errors; this is not a uniform error estimate.

Resolution matters as well: at omega=.4, b=.01, regular errors for 32/64/96 grids
are .08057/.02021/.009250, whereas small-branch errors are
.01355/.02973/.03347. Refining h shrinks Delta with P fixed, so holding b fixed
can move the approximation **away** from its b/Delta->0 limit. A fixed global
Tlim cannot simultaneously certify numerator linearization and regular-peak
resolution across grids/species.

## 7. Explicit Tlim=1e-5 contour experiments

All contours below use the unmodified production executable and unchanged A1
sampling controls, with no growth search and samples enabled. Only the explicit
Tlim or lower-bound diagnostic input differs.

| Numerical representation / gamma_min / Tlim | A1 status | accepted boundary winding | calls | passes / refinements / depth | wall s |
|---|---:|---:|---:|---|---:|
| Numerical p, analytic e / .01 / .01 | 0 | 2 | 626 | 3 / 89 / 5 | 29.83 |
| Same / .01 / 1e-5 | 0 | 2 | 626 | 3 / 89 / 5 | 30.09 |
| Same / .001 / .01 | 5 | NONE | 73 | 1 / 35 / 24 | 1.92 |
| Same / .001 / 1e-5 | 0 | 2 | 684 | 3 / 118 / 7 | 33.71 |
| Analytic p, numerical e / .01 / .01 | 0 | 3 | 586 | 3 / 69 / 4 | 10.37 |
| Same / .01 / 1e-5 | 0 | 3 | 628 | 3 / 90 / 6 | 13.25 |

These are local four-rank timings from the first completed audit, not general
performance estimates. Logical requests equal physical calls (no growth cache
is enabled). Every rank's physical-call count agrees.

For protons gamma_min=.01 gives b_min=.0333333, already above both thresholds.
The two original-contour determinant tables are **byte-identical at all 257
matched frequencies**, with median/max analytic discrepancies .00408394/.212333
and identical accumulated phase 4*pi. Lowering Tlim cannot improve this contour's
proton response because it does not change its branch.

The lower contour does change branch and sampling status. Its accepted winding
2 at Tlim=1e-5 is evidence about that computed finite trajectory, not proof of an
analytic approximation. The half-cell jump at gamma=.001 decreases from 2.388
to .07153 but does not vanish. At .01 it remains .009428; at .03 .002576.
`tiny_half` values are byte-identical to the forced-regular half-cell values.
The imaginary-frequency switch has moved below the contour, while the
real-frequency partition/interpolation mismatch remains.

### Limited numerical-electron control

The supplied mass-scaled support is correct: p_max=.1616836943893 and
h_e=.0025263077483. At gamma=.01, b_e=.0000181539007162. Default Tlim=.01 selects
the small-imaginary treatment throughout the box; Tlim=1e-5 selects regular
integration throughout, since gamma_switch,e≈.00550846<.01.

The accepted computed winding stays **3**, not the analytic reference's 2.
Matched determinant median/max errors change from .0551063/11.3237 to
.0215556/1.58018. The maximum cumulative phase discrepancy remains approximately
one extra turn (6.4093 versus 6.3738 radians), and total phase remains 6*pi.
The regular branch changes and improves these sampled error statistics, but does
not remove the disagreement. Near the lower boundary b_e can be smaller than a
local regular substep; the proton quadrature tests explain why branch selection
alone is insufficient. No electron convergence study or repair is attempted.

## 8. Validation code, artifacts and regression scope

Entry point:

```bash
bash validation/run_nyquist_c2a.sh
```

It builds in fresh external scratch, generates the standard distributions,
reproduces C1 before compiling new susceptibility instrumentation, runs the
paired/forced-branch cases, checks actual-rank call counts, evaluates the six
unmodified-production contours and postprocesses tensor/integral closure. It
also compiles the actual copied numerical functions with full debug checks for
the standalone mathematical tests. `tests/data/nyquist_c2a/README.md` documents
all file formats and runnable inputs.

New files are confined to this report, `validation/c2a/`, the wrapper,
`tests/data/nyquist_c2a/`, and a validation README entry. In particular no ALPS
source or existing scientific regression fixture changes. `input.md` and
`nyquist.md` retain their existing split; no production parameter is added.

Instrumentation exposes the original private helpers only in a copied module.
It snapshots the symmetric core before adding the remainder and background,
records selected interpolation values/stencils, every numerical harmonic's
returned integral, and the assembled species tensors and wave matrix.
Additional diagnostic evaluations occur after the original sum. It does not
splice phase contributions or calculate a second production susceptibility.
The standalone numerator tests are algebraic fixtures, not physical VDFs; the
physical cases always use the normal numerical derivative machinery.

Postprocessing verifies:

- no duplicate harmonic ownership and the same requested harmonic/component set
  on both sides of every pair;
- background+core+remainder equals the returned integral;
- summed per-harmonic changes equal the assembled proton-tensor change, including
  the zero change of the constant `int_ee` term;
- matched analytic, regular and small-branch coordinates are exactly equal;
- the six default instrumented reproduction values equal the uninstrumented
  values byte-for-byte;
- lowering Tlim on the original proton box changes no matched determinant, and
  tiny-Tlim half-cell values equal the forced regular calculations.

The first completed optimized study is `/tmp/alps-nyquist-c2a/full`, including
`build/provenance.txt`, copied `instrumented/fns.f90`, per-run controls, inputs,
coordinates, tensors and integral traces. The final packaged rerun is identified
as `/tmp/alps-nyquist-c2a/final`: its fresh end-to-end entry point completed with
exit 0, including all closure checks and both mathematical compilations. It
contains 50 fixed-frequency diagnostic runs and six full-contour calculations.
The final observations reproduce the first study's statuses, counts and
determinants; timing differences are normal local-run variation.
Representative raw values and compact tables are retained under
`tests/data/nyquist_c2a/results/`; full temporary traces remain in the scratch
runs. Large VDF arrays are reproducibly generated, not added as permanent dumps.

Compiler is GNU Fortran 11.4.0, Open MPI 4.1.2, four local MPI ranks, one BLAS/OpenMP thread.
Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`.
Standalone debug compilation of the copied integration functions uses
`-O0 -g -fcheck=all -fbacktrace`; it reproduces the cubic, constant-window,
finite-window-residue and discarded-remainder identities. The regular constant
finite-window test's normalized residual is 2.98e-9 (acceptance 1e-6).

`validation/run_nyquist.sh` was rerun after the audit with
`ALPS_REFERENCE_EXE=/tmp/alps-nyquist-c1/baseline/ALPS`: PASS. It includes standalone
A1 debug/optimized checks, integrated contours/normalization, MPI error shutdown,
ordinary growing/damped roots, maps and line/plane scans. All nine ordinary
scientific output files compare byte-for-byte with that identified pre-C1
executable. No source correction needs a broader physics regression here;
complete A1/A2/A3/log/mirror/B1/B2 debug/optimized results from C1 remain applicable
to the unchanged sources and are not claimed as newly rerun C2 suites.
The new A1/A2 regression artifacts are `/tmp/alps-nyquist-a2.kkljvlK6`.

## 9. Ranked corrective experiments — not implemented

1. **Make the existing representation and partition internally consistent.**
   Replace nearest-node Taylor-line switching with a common cell interpolation
   that matches values at shared boundaries, and ensure transferred cells use
   the same numerical integral regardless of whether labeled background or
   local. Account for every R>0, including R<step. Test these separately with
   the cubic and constant fixtures, then together. Fixing only the interpolation
   cannot remove the constant-numerator partition jump; fixing only the window
   transfer cannot remove the cubic interpolation jump. No smoothing/averaging
   is justified as a substitute for consistency.
2. **Match finite-window small-b and regular calculations using error controls.**
   Test the finite-window 2*atan(Delta/b) constant term and the corresponding
   linear-numerator integral before considering higher-order numerator terms.
   Study b/Delta and b/step, not just a global absolute Tlim. Require an overlap
   where both branches have independently controlled accuracy; no choice of a
   switching threshold makes two unequal formulas equal. Preserve causal signs
   and the existing lower-half-plane continuation contract in any later change.
3. **If needed, evaluate a fixed-cell Cauchy-integral formulation in a separate
   design review.** Exact integration of a common piecewise polynomial numerator
   against 1/(p-a), or suitably controlled quadrature, can avoid an unresolvable
   moving central peak. A fixed-grid analytic Cauchy transform offers a clearer
   upper-half-plane analyticity argument than ad hoc blends of real/imaginary
   rules, but would be a more substantial numerical-method change.

Do not choose among these by winding agreement alone. First validate values,
phase, one-sided limits, quadrature convergence and momentum convergence against
analytic integrals/bi-Maxwellians. Continuity is necessary; it does not by itself
prove holomorphy of an approximation depending explicitly on Re(a), Im(a), moving
bounds and stencil selections. Only after that review should numerical-electron
convergence, growth searches or numerical wavevector scans resume.

Confirmed facts are the recorded finite jumps, two independent half-cell
counterexamples, the positive-remainder omission, component attribution and
finite-window/regular-resolution tests. Deductions are the third-difference
identity, partition-transfer condition and finite-window residue correction.
The unique cause of every secondary node/quarter-cell event, and of the remaining
numerical-electron winding-three discrepancy, remains unverified; neither is
silently assigned a physical interpretation or repaired here.
