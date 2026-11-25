# C1 and initial C2 — nonrelativistic numerical VDFs

C1 enables numerical/mixed species at one prescribed wavevector through the
existing collective determinant. Numerical proton Gaussian representations
reproduce the analytic reference's accepted boundary winding of 2. Increasing
momentum resolution improves matched-frequency determinant agreement. This is
initial validation of the computed trajectory, not certification of a physical
root inventory: a measurable discontinuity at the existing resonance-treatment
switch remains. No integration correction is included in this change.

## Source changes and scope

Only two production files change:

- `src/ALPS_nyquist_driver.f90`: `read_nyquist` validates common density, charge,
  mass and relativistic status per species; analytic beta/anisotropy/drift checks
  apply only to analytic species. Option 3 accepts numerical species and mixtures;
  options 4/5 still require all warm analytic species. Numerical growth searches
  are rejected during C1 (`nyq_find_max_growth=F` required). `validate_nyquist_vdf`
  checks finite, nonnegative data, uniform increasing rectangular coordinates,
  positive samples for logarithmic continuation, and sufficient grid dimensions.
  Missing files are detected before ordinary VDF reading. Summary diagnostics
  identify the actual species representation and avoid unused analytic defaults.
- `src/ALPS.f90`: call that validation hook immediately after the ordinary
  `read_f0`, before differentiation. Analytic and ordinary operations retain their
  initialization and dispatch.

`ALPS_io.f90`, `ALPS_analyt.f90`, `ALPS_NHDS.f90`, `ALPS_fns.f90`, `ALPS_com.f90`,
the generator, A1 and A3 are byte-identical to the pre-C1 source snapshot.
In particular no susceptibility, resonance, quadrature, derivative, continuation,
root solver, MPI assignment, collective determinant, normalization, cache or
sampling algorithm changes were made. No new namelist parameter or dependency
was introduced. Documentation respects the separate `input.md` reference and
`nyquist.md` guide.

Validation additions: `validation/run_nyquist_vdf.sh`,
`make_nyquist_vdf_cases.sh`, `nyquist_vdf_moments.f90`,
`observe_nyquist_vdf.f90/.awk`, `compare_nyquist_vdf.f90`, and
`audit_nyquist_vdf.awk`; runnable inputs and retained observations are under
`tests/data/nyquist_vdf/`. `validation/README.md` links the entry point.
`run_nyquist_scan.sh` excludes only the newly authorized species representation
and relativity labels from its old-summary byte comparison, alongside its
existing timing exclusion. Every prior scientific field remains compared;
growth/sample files and ordinary scientific fixtures are not edited.

Unchanged SHA-256 values:

| Source | SHA-256 |
|---|---|
| A1 `ALPS_nyquist.f90` | `d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31` |
| A3 `ALPS_nyquist_growth.f90` | `5eda3a5b6faf58ca54417b5bf26f4bd89abe6a2e47d1f8ac9978a73434d1a212` |
| `ALPS_fns.f90` | `567e05bad770489d578c0024b4230376f1e000c2f04b7e84fc2f7d065d446be5` |
| `ALPS_com.f90` | `989cfbb8fd3a8252e195d868a93a09552cfbce3e4f1172f93c6a127e947e5cc8` |

## Shared initialization and MPI contract

The normal startup reads the numerical arrays, constructs the existing central
finite differences, initializes continuation models/fitting, transfers the
completed equilibrium, determines `nmax`, calls `split_processes` and constructs
worker Bessel arrays. All preparation is outside frequency evaluation. These
fixtures use the existing nonrelativistic method-2 log-Chebyshev representation,
kind 1, order 8 and default derivative setting 1. Native initialization and
`determine_param_fit` remain on their normal paths; no continuation work is
skipped for the upper-half-plane calculation.

Rank 0 calls the same `run_nyquist`; uncached frequencies enter the original
command broadcast and collective `disp`, with the existing reductions,
determinant broadcast and barrier. Workers retain species/harmonic assignments.
No communicator or frequency/wavevector parallelism is added. Every normal,
fixed-frequency and failure test checks equality of determinant counts across
four ranks. Ordinary non-growth contours retain the existing uncached behavior;
the A3 exact raw cache implementation is unchanged.

Copied-adapter observations for the primary mixed fixture retain nmax=(15,1).
Worker ranks 1/2 process proton harmonics 0:7 and 8:15, with Bessel-array bounds
-1:8 and 7:16; rank 3 handles the analytic electron species' assigned 0:1 range,
with bounds -1:2. These are the original assignment routine's results, including
its halo bounds. The audit observes initialized arrays and does not alter them.

Rectangles use raw `D_ALPS`; horizontal paths still use
`D_ALPS/(z*vA)^6` through the existing log-polar adapter. The numerical summary's
interpretation states explicitly that a converged result is an accepted
**numerical boundary winding** pending physical convergence evidence.

## Reproduction and diagnostic utilities

Run `bash validation/run_nyquist_vdf.sh`, optionally with debug `FFLAGS`.
Runnable analytic, mixed, fully numerical and isotropic inputs plus generation
instructions are in `tests/data/nyquist_vdf/`. Arrays are produced by the original
`distribution/generate_distribution.f90`; no custom derivatives or finite-grid
renormalization are used. `validation/make_nyquist_vdf_cases.sh` emits every study
input in a fresh external directory.

- `nyquist_vdf_moments.f90` reads the actual arrays. It reports ALPS's full-node
  rectangular quadrature and an independent trapezoid rule. Temperatures are
  conditional moments normalized by measured density: beta_parallel =
  `2*(<p_parallel^2>-<p_parallel>^2)/m`, beta_perp = `<p_perp^2>/m` for input n=1.
  Density error is reported separately; it is not silently corrected.
- `observe_nyquist_vdf.f90/.awk` replace only adaptive sampling in a **copied
  validation adapter** with fixed requests to its unchanged collective callback.
  Each diagnostic is deliberately marked A1 status 6, no accepted integer, and
  exits 1 after normal collective shutdown. It is not a converged Nyquist count.
- `compare_nyquist_vdf.f90` requires exactly identical frequencies and successful
  evaluator statuses. It computes relative complex-determinant discrepancies and
  calls A1 `nyq_accumulate` for continuously unwrapped phases. ASCII output keeps
  both complex determinants, relative errors and cumulative phases.
- `frequencies.dat`: 64 intervals on each CCW edge, repeated closing corner
  (257 contour samples), followed by separate Re(z)=0.4 gamma probes at the
  resonance transition. Initial optimized observations used nine gamma probes;
  extended observations bracket .005/.02 Tlim thresholds (15 probes), add six
  moving-window samples, and final debug diagnostics repeat jumps at separations
  of 2e-10 and 2e-12. Supplemental coordinates are not a closed contour and have
  NaN cumulative phases in the fixed nine-column comparison table.

Fixed-grid comparisons do not rely on coincident adaptive trajectories. Normal
adaptive calculations independently establish the sampling statuses. Statistics
below concern only the 257 identical contour coordinates, with the repeated
closing corner included. Finite sampled errors below unity do not prove a uniform
boundary bound or continuity.

## Equilibrium and momentum moments

All species have density 1, zero drift, proton mass/charge 1, electron mass
0.000544617021487 and charge -1, vA/c=1e-4. k=(0.001,0.3). Beta_parallel,p =
beta_parallel,e = 3, proton anisotropy 2, electron anisotropy 1. Thus Te equals
Tparallel,p. The isotropic control changes only proton anisotropy to 1.
The established ordinary-solver reference is ±0.39666+0.085986i; legacy root
regressions check these independently of Nyquist.

Let S be the number of Gaussian e-folding widths of support. Explicit generator
bounds are pperp_max = S sqrt(3 m alpha), |ppar_max| = S sqrt(3 m).
The generator's automatic `maxP` convention omits sqrt(beta); the fixtures use
explicit bounds to avoid that ambiguity. Both numerical species use their own
mass-scaled support, not proton-scale electron limits.

Measured proton moments (ALPS nodal quadrature, before susceptibility):

| Intervals / S | Density | beta_parallel | beta_perp | anisotropy |
|---|---:|---:|---:|---:|
| 32×64 / 4 | 0.9973916751 | 2.999999101 | 6.015708584 | 2.005236796 |
| 64×128 / 4 | 0.9993486053 | 2.999998817 | 6.003903388 | 2.001301919 |
| 96×192 / 4 | 0.9997104902 | 2.999998709 | 6.001728015 | 2.000576866 |
| 48×96 / 3 | 0.9992290730 | 2.998956416 | 5.998297761 | 2.000128354 |
| 96×192 / 6 | 0.9993487038 | 3.000000000 | 6.003911854 | 2.001303951 |

Drifts are at roundoff (proton velocity magnitude below 1e-16). The 3/4/6 support
study holds the mesh spacing fixed (S/nperp=1/16); the 32/64/96 study holds S=4
fixed. This separates truncation from resolution. At S=4 the central density
quadrature error dominates the tiny Gaussian tails. The S=3 truncation changes
temperature moments appreciably; extending support without refining spacing does
not eliminate the finite-difference/quadrature error.

ALPS's ordinary charge diagnostic therefore shows small mixed-case quadrature
imbalances Np-1: approximately -0.0026083, -0.0006514 and -0.0002895 on the
32/64/96 grids. These are recorded representation errors; declared densities
remain equal and no normalization is adjusted to hide them. In the fully
numerical case the mass-scaled grids give equal measured density deficits for
both species, so its net charge quadrature error is at roundoff. That case still
has the extra boundary turn; mixed-case charge imbalance is not its sole cause.

## Contour and matched-frequency results

All initial runs use the same [-1,1] × [.01,.2] rectangle, growth disabled,
samples enabled, and **unchanged A1 defaults** (initial intervals 8, depth 24,
evaluation cap 200000, up to 6 global passes). Bessel_zero=1e-45,
positions_principal=3, n_resonance_interval=100, Tlim=.01 are held fixed.

Measured optimized initial-contour results:

| Case | W | disp calls | passes / refinements / depth | wall seconds |
|---|---:|---:|---:|---:|
| Analytic | 2 | 630 | 3 / 91 / 6 | 0.0324 |
| Numerical p, 32×64 / S4 | 2 | 612 | 3 / 82 / 5 | 10.56 |
| Numerical p, 64×128 / S4 | 2 | 626 | 3 / 89 / 5 | 37.14 |
| Numerical p, 96×192 / S4 | 2 | 626 | 3 / 89 / 5 | 60.49 |
| Numerical p, 48×96 / S3 | 2 | 632 | 3 / 92 / 5 | 16.91 |
| Numerical p, 96×192 / S6 | 2 | 626 | 3 / 89 / 5 | 81.75 |
| Analytic isotropic | 0 | 630 | 3 / 91 / 6 | .0253 |
| Numerical p isotropic, 64×128 / S4 | 0 | 630 | 3 / 91 / 6 | 29.35 |
| Both numerical, 64×128 / S4 | 3 | 588 | 3 / 70 / 4 | 27.98 |
| Numerical e only, 64×128 / S4 | 3 | 586 | 3 / 69 / 4 | 10.03 |

The unstable numerical-proton cases have total phase 4π; isotropic controls
give 0 and the numerical-electron cases give 6π. Timings are single local
four-rank measurements, not numerical-VDF performance predictions. Every call
count is also the logical request count; cache hits are zero in these non-growth
calculations. Samples are written after sampling and timed separately.

Shared startup times in the same optimized study were .0012 s analytic,
.0719/.2907/.6080 s for the 32/64/96 mixed grids, .6154 s for S6 and .5437 s
for both numerical species. No startup optimization is attempted; numerical
integration dominates these small fixtures. For example the 64-grid mixed
contour spends about 37.14 s in its 626 physical collective calls. Timings
exclude the separate fixed-frequency probe runs, which have their own counts
and timings in `probe_*` summaries.

Matched-frequency errors, dimensionless:

| Numerical p grid | minimum | median | 90th percentile | maximum | max cumulative phase difference, radians |
|---|---:|---:|---:|---:|---:|
| 32×64 / S4 | .006801 | .009244 | .032303 | 1.160986 | .199002 |
| 64×128 / S4 | .001385 | .004084 | .007717 | .212333 | .054455 |
| 96×192 / S4 | .000997 | .002341 | .003654 | .065075 | .024528 |
| 48×96 / S3 | .002627 | .005190 | .010795 | .209778 | .057498 |
| 96×192 / S6 | .001386 | .004082 | .007718 | .212338 | .054593 |

At fixed spacing, extending S=4 to S=6 barely changes these discrepancies;
resolving the Gaussian more finely has a much larger effect than extending this
already broad support. This is sampled evidence for these fixtures, not a
universal momentum cutoff.

Error-bin counts for thresholds .001/.01/.1/1 (257 contour samples):
32 grid: 0/140/110/6/1; 64 grid: 0/236/20/1/0; 96 grid: 1/245/11/0/0.
The worst error occurs near the small analytic determinant at z=0+.01i:
D_analytic≈-1.4851079692e-4. The coarse value is about -3.20930e-4-9.97e-8i.
The relative error is amplified by cancellation here; matching winding alone
would have hidden that discrepancy.

Fixed-contour edge phases (bottom/right/top/left), radians:

- Analytic: -2.5383774473, 1.2067359563, 12.6912761490, 1.2067359563.
- Numerical p 96×192: -2.5363402164, 1.2049618265, 12.6902347894, 1.2075142149.

Both accumulated totals are 4π. Differences include small signed-frequency
asymmetries from the computed numerical response; no symmetry or count doubling
is imposed.

## Resonance transition and limitations

For protons, Im(p_res)=m gamma/kpar, so Tlim=.01 corresponds to gamma=.003.
At fixed Re(z)=.4, compare gamma=.003 and .00300001. The analytic determinant's
relative change is 1.44354e-7. The numerical changes are **6.92%, 2.06%, and
2.82%** on 32×64, 64×128 and 96×192 grids respectively. Thus the switch is
measurably discontinuous at these finite resolutions, and the size is not
monotone with resolution. Changing interpolation/window positions can also
contribute grid-sensitive error; this study does not attribute all determinant
error to that one switch.

The source identifies the relevant branch in existing `integrate_res`:
`abs(denomI)>Tlim` selects the regular moving-window integral, while the other
branch uses the small-imaginary-part approximation. The tiny gamma probe crosses
that branch without changing real resonance location or the input distribution.
This is evidence of a numerical treatment mismatch, not a physical response
singularity. Representation moments and the above-transition contour improve
with resolution, but this does **not** establish uniform upper-half-plane
analyticity of the computed response.

No `integrate`, `integrate_res`, `determine_resonances`, `funct_g`, Tlim default,
quadrature or analytic-continuation correction was attempted. A focused C2 task
should measure branch matching and moving-stencil continuity at fixed frequency
before changing either formulation. Increasing contour refinement alone cannot
repair a discontinuous determinant, and an accepted integer is not proof of its
physical interpretation.

Tlim sensitivity at the 64×128 grid confirms that the jump moves with the
threshold: Tlim=.005 gives a relative step .01104 at gamma=.0015, while
Tlim=.02 gives .04071 at gamma=.006. At the original gamma=.003 these modified
cases are smooth to the probe scale (~1e-7). All three above-transition original
rectangles retain the same count 2 and 626 requests. Their fixed-contour
determinants agree because gamma>=.01 exceeds all three proton thresholds.

## Regression and failure validation

The optimized C1 study is retained in `/tmp/alps-nyquist-c1.ULT1r0wE`;
representative moments, comparison distributions, full matched determinants,
summaries, the unresolved prefix and assignment observations are also retained
under `tests/data/nyquist_vdf/observations/`. Large VDF arrays are reproducibly
generated rather than copied into the repository. Source/build provenance is
retained with the results. The study used separate initial/resumed diagnostic
runs as the secondary discrepancy was characterized; the final entry point also
includes numerical-electron isolation and finer transition probes.

Failure checks distinguish an exhausted budget (A1 status 4, 8 collective calls),
injected nonfinite determinant (status 2, one collective call), and the lower
contour's refinement failure (status 5, 73 calls). All preserve `NONE` rather
than a zero winding and finish collectively. The injected value is assigned in
a copied adapter **after** the actual `disp` call, with original susceptibility
and MPI implementations intact. Input tests reject relativistic species, cold
analytic species, numerical options 4/5, numerical growth searches, missing VDFs
and negative samples. The final suite also tests nonfinite/nonuniform source
data and an identically zero distribution before differentiation.

Standalone comparison checks reject mismatched coordinates, unequal record
counts and a nonfinite value in the supplemental probes. They also check that
comparison tables retain exactly nine numerical fields with NaN supplemental
phases. No independently sampled trajectory is silently paired by row index.

Identified compiler: GNU Fortran 11.4.0, existing MPI/BLAS/LAPACK build. Four local
MPI ranks; `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1`. Optimized flags:
`-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags:
`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`. Full existing suites run
sequentially, followed by C1 debug, so MPI startup does not overlap. The pre-C1
source snapshot is `/tmp/alps-nyquist-c1/original`; its identified executable is
`/tmp/alps-nyquist-c1/baseline/ALPS`.


All six existing suites passed in **both optimized and debug configurations**:

| Entry point | Optimized artifacts | Debug artifacts | Result |
|---|---|---|---|
| `run_nyquist.sh` | `/tmp/alps-nyquist-a2.ST46sxSE` | `/tmp/alps-nyquist-a2.kwmK8OMR` | PASS / PASS |
| `run_nyquist_growth.sh` | `/tmp/alps-nyquist-a3.nfgYOLSJ` | `/tmp/alps-nyquist-a3.QupINb7s` | PASS / PASS |
| `run_nyquist_log.sh` | `/tmp/alps-nyquist-log.0inuqDXc` | `/tmp/alps-nyquist-log.4O3EyVSN` | PASS / PASS |
| `run_nyquist_mirror.sh` | `/tmp/alps-nyquist-mirror.zbhPRpzc` | `/tmp/alps-nyquist-mirror.6lmtjSaF` | PASS / PASS |
| `run_nyquist_scan.sh` | `/tmp/alps-nyquist-b1.CadbxmqN` | `/tmp/alps-nyquist-b1.bC8nQDam` | PASS / PASS |
| `run_nyquist_plane.sh` | `/tmp/alps-nyquist-b2.fr1M6I0y` | `/tmp/alps-nyquist-b2.NYl5HRCB` | PASS / PASS |

A1 additionally runs 306 checks each in standalone debug and optimized builds;
A3 reports 10056 growth/cache checks and the log suite 72 mathematical checks.
The original linear reference still requires exactly 11 steps and 5005 physical
calls. The logarithmic, weak-growth, mirror-root/bracket, strip, budget,
nonfinite-value, MPI shutdown, B1/B2 coordinate, state-rebuild and aligned-output
assertions remain intact.

Both configurations compare the nine ordinary root/map/line/plane scientific
files byte-for-byte with the identified pre-C1 executable. B1 also compares
option-3 growth and samples byte-for-byte and every pre-existing non-timing
summary field. B2 compares the unchanged option-4 compact data byte-for-byte,
detailed scientific columns and aggregate summaries excluding only existing
timings. Source control fixtures and legacy output files were not changed.
The deliberate single-point summary additions are the species response and
relativity labels; numerical species receive applicable VDF fields instead of
unused analytic parameters.

## Additional C2 findings

The fully numerical 64×128 case accepts computed boundary winding **3**, not the
analytic reference's 2. This is a physical disagreement, not a failed MPI job or
an analytic root-count expectation imposed on numerical data. Its measured
numerical-electron density is 0.9993486053, beta_parallel=2.999998817,
beta_perp=3.001951694, drift velocity 2.74e-16; support is
|ppar|max=pperp_max=0.1616836943893. Electron moments are therefore close to target,
but alone do not establish accurate complex susceptibility.

An additional analytic-proton/numerical-electron control also accepts winding 3.
It localizes the extra boundary turn to replacing the electron response at this
resolution. Fully numerical matched median/max determinant errors are
5.24%/1109.54%; the electron-only control gives 5.51%/1132.37%. Both differ most
strongly near the small analytic determinant around z=0+.01i, and both fixed
trajectories accumulate 6π instead of 4π. Numerical-electron refinement has not
been established in this initial study. The electron resonance imaginary
momentum m_e*gamma/kpar is well below Tlim throughout this contour, so the
small-imaginary approximation is a plausible contributor alongside finite-grid
cancellation and derivative error. No specific electron tensor defect is claimed
without a separate susceptibility-level comparison. This secondary case does
not validate fully numerical physical root counting.

The mixed 64×128 contour lowered to gamma_min=.001 returns **unresolved A1
status 5**, no accepted winding, after 73 calls and depth 24. Its accepted prefix
ends near Re(z)=-0.01623799 on the bottom edge; the reported partial phase is not
a complete-contour winding. This is below the Tlim transition and near a moving
window's half-cell boundary. For this grid, delta_ppar=8 sqrt(3)/128, and
Re(z)=-kpar*delta_ppar/2=-0.0162379763209582. Fixed-frequency probes differing by
1e-8 in real frequency show analytic relative changes ~1.22e-6, but numerical
changes of order unity (including 5.52). Even at gamma=.01 the corresponding
numerical steps are ~.0201 and .0107 versus analytic ~1.03e-6. These are confirmed
finite-grid discontinuities; the `integrate_res` half-cell window-boundary change
and `funct_g`'s nearest-node centered-linear stencil changes are the likely
mechanisms. This is distinct from the imaginary-frequency Tlim
switch. Isolating every harmonic/tensor contribution is deferred to focused C2
work, with the input and matched data preserved.

These findings establish why contour sampling convergence and physical
susceptibility convergence must be reported separately. The A1 failure safeguard
correctly prevents assigning a count to the problematic lower contour. The
original upper contour can pass A1 despite smaller discontinuities, so its
integer agreement is encouraging evidence rather than certification.


## Final C1 debug validation

`bash validation/run_nyquist_vdf.sh` completed end-to-end with full debug checks
in `/tmp/alps-nyquist-c1.JoB4Z3i1` (exit 0). The optimized study and resumed
diagnostics described above also completed; the final debug entry point covers
all their cases plus the final validation guards. All source-data, backend,
collective-failure and comparison-format checks passed. Secondary winding 3 is
recorded as a physical discrepancy, not asserted to equal the analytic count.

Final production hashes recorded by that build:

- `ALPS.f90`: `acb9e48ab671964ad537c73651d413b1e6e34f2af43ab3174545de6a5109097b`.
- `ALPS_nyquist_driver.f90`: `ce5f67f48552963bcee7de38774a316fa28d36f21cde919eb89f02bd428d4d8d`.

Debug initial winding, sampling status, logical/physical call counts and cache
hits match optimized for every study case, including all three expected
unresolved outcomes (lower contour, budget, injected determinant). An additional
Fortran comparison of the 257 identical boundary frequencies in **all twelve**
normal/probe case pairs found exactly equal printed complex determinant values
between optimized and debug builds: maximum relative discrepancy zero. Thus
the analytic/numerical disagreement is reproduced across compiler settings.

Representative debug contour wall times, seconds: 32/64/96 mixed grids
25.89/73.90/142.85; S3/S6 support cases 41.31/194.91; isotropic mixed 73.09;
both numerical 69.14; electron-only 23.77. Every contour's physical count remains
the one in the optimized table. Full diagnostics and current build provenance
are retained in `tests/data/nyquist_vdf/observations/debug/`.

At 64×128 resolution the finer two-sided probes confirm nonvanishing jumps.
Here delta is the displacement from the switch coordinate; values compare
`D(+delta)` with `D(-delta)` relative to the latter magnitude:

| Switch / fixed coordinate | delta | analytic relative change | numerical relative change |
|---|---:|---:|---:|
| Tlim: Re(z)=.4, gamma=.003 | 1e-10 | 2.8871e-9 | .02055368 |
| Same | 1e-12 | 2.8894e-11 | .02055368 |
| Half-cell: Re(z)=-.0162379763209582, gamma=.001 | 1e-10 | 2.4480e-8 | 2.387796 |
| Same | 1e-12 | 2.4480e-10 | 2.387796 |
| Same real coordinate, gamma=.01 | 1e-10 | 2.0513e-8 | .00942822 |
| Same | 1e-12 | 2.0515e-10 | .00942822 |

These are resolved coordinate differences in working precision. The numerical
changes plateau while analytic changes scale with separation. This strengthens
the continuity finding without assigning a unique faulty tensor contribution.

## Recommended C2 scope

C1 execution and the primary mixed-equilibrium comparison are established for
these single-wavevector fixtures. The initial C2 representation study moves the
above-transition determinant toward the analytic response and retains the
expected primary windings. **Continuity and fully numerical physical root
counting are not yet validated.**

Keep numerical maximum-growth estimation and numerical options 4/5 disabled.
A focused follow-up should first isolate individual susceptibility/harmonic terms
across the small-imaginary/regular integration switch and moving-window/stencil
boundaries, with fixed frequencies and decreasing probe separations. Review
matching and interpolation continuity before changing either algorithm. Then
repeat electron resolution/support convergence and re-establish closed-contour
agreement, including lower gamma boundaries, before extending the physical scope.
The current data do not justify a universal VDF resolution or stability cutoff.
No corrective integration change is part of C1.
