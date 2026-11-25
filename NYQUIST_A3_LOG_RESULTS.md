# ALPS Nyquist A3 extension — logarithmic strips and oblique mirror validation

Date: 2026-11-21 (America/Phoenix).

**Final tolerance update:** Section 8 supersedes the historical logarithmic
absolute-tolerance results in Sections 1–7. Linear searches remain absolute;
logarithmic searches now require relative bracket-width convergence.

**Implemented and tested.** `nyq_growth_spacing='linear'` remains the default and
reproduces the A3 reference bracket and 11 iterations. Optional `'log'` uses a
geometric trial and the nearest original gamma boundary, with the correct
lower-strip complement count. The original absolute width tolerance is unchanged.

The new physical fixture identifies a proton mirror branch at the single oblique
wavevector `(kperp,kpar)=(0.3,0.05)`, with growth approximately
`0.0381697286246`. Existing ordinary root solving and polarization diagnostics,
anisotropy controls and harmonic convergence establish the reference before
Nyquist assertions. The measured winding is one; both linear and logarithmic
brackets contain the independent ordinary-solver growth.

Nearest-edge logarithmic strips are **not an established cost improvement** over
cached full-height logarithmic rectangles. The existing reference with lower
boundary 0.01 requires two more physical calls; with lower boundary 1e-4 their
sampling work is identical. No sampling safeguard was relaxed to improve this
comparison. Options 4 and 5 remain unimplemented.

## 1. Implementation and preserved components

| File | Extension |
|---|---|
| `src/ALPS_nyquist_growth.f90` | Optional trailing `spacing` argument, safe/testable `nyq_log_trial`, original-edge selection, inferred `N_above`, consistency checks, history metadata and appended columns. Existing call sites without the argument remain linear. |
| `src/ALPS_nyquist_driver.f90` | New namelist string defaulting to `'linear'`, validation when growth estimation is enabled, argument forwarding, explicit spacing/original-bound summary fields. |
| `validation/check_nyquist_growth.f90` | Shared observation readers now also read spacing, original real bounds and appended history fields; legacy 13-column histories remain readable. Existing linear assertions are retained. |
| `validation/test_nyquist_log.f90` | Synthetic analytic counts, independent full rectangles, weak/distinct growth rates, both anchors, coordinate and failure checks. |
| `validation/check_nyquist_log.f90`, `run_nyquist_log.sh` | Linear reference preservation, log searches, every accepted trial checked uncached, cache accounting, and a copied-source cached full-height baseline. |
| `validation/observe_mirror.f90` | Validation-only calls to existing `disp`/`calc_eigen` after ordinary roots, high-precision frequency, residual, field and density diagnostics. |
| `validation/make_nyquist_mirror_cases.f90`, `check_nyquist_mirror.f90`, `run_nyquist_mirror.sh` | Fixed-wavevector mirror fixtures, ordinary references/controls, contour counts and growth comparisons. |
| `input.md`, `validation/README.md`, `tests/data/nyquist_log/`, `tests/data/nyquist_mirror/` | Input/output documentation, examples and regression entry points. |

A1 is unchanged, SHA-256
`d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.
The cache implementation, exact-match rules and determinant normalization are
unchanged. No susceptibility, `disp`, ordinary root solver, analytic continuation,
harmonic assignment, MPI request protocol or production main was modified.
Warm nonrelativistic analytic bi-Maxwellian restrictions remain in force. No
numerical dependency, physical backend, wavevector scan or user-facing optimization
switch was added beyond `nyq_growth_spacing`.

The parameter applies only to an enabled maximum-growth calculation. Accepted
values are lowercase `'linear'` and `'log'`; unsupported enabled values are rejected
with a specific input diagnostic. Existing tolerance, iteration limit, bounds,
normalization and `nyq_sampling` defaults are retained.

## 2. Counting and numerical conventions

Let the original gamma bounds be `[gamma_0,gamma_top]` and its accepted count N0.
For log spacing, the current positive finite bracket `[L,H]` gives

`M=exp(0.5*log(L)+0.5*log(H))`.

This avoids multiplying endpoints. If M is not strictly inside the bracket in
working precision, the search returns the existing coordinate-resolution status
6 without changing the last valid bracket. The stop criterion remains the
**absolute** `H-L < nyq_growth_tolerance`. The output midpoint estimate remains
the arithmetic midpoint of the final bracket; it is not a dispersion root.

The distances `abs(log(M)-log(gamma_0))` and
`abs(log(gamma_top)-log(M))` choose an original boundary. Equal distances in
working precision select the original upper boundary. Near-equalities follow the
computed distances deterministically. In particular, the first 0.01–0.2 trial
selects the lower boundary after rounding; this choice is recorded and tested,
not adjusted to obtain a favorable benchmark.

All rectangles are independently evaluated counterclockwise by unchanged A1:

| Anchor | Actual strip | Accepted strip winding | Count above M in original domain |
|---|---|---|---|
| Original upper | `[M,gamma_top]` | W_upper | `N_above=W_upper` |
| Original lower | `[gamma_0,M]` | W_lower | `N_above=N0-W_lower` |

These equalities follow from root-count additivity for a pole-free response with
no boundary zeros: the two strips partition the original rectangle and their
counts sum to N0. The lower strip may contain roots already below current L;
subtracting its count from the **original** N0 accounts for them correctly.

Only accepted, converged integer windings are used. Strip counts must lie in
`[0,N0]`, and inferred counts must also be consistent with the retained positive
lower count. Negative or inconsistent values are rejected, never clipped. A
positive inferred count sets L=M; zero sets H=M. Failed strips leave the bracket
unchanged and do not become stable results.

Every strip retains A1's local adaptation, stricter global passes and independent
seed. Only exact raw determinant values are reused; no edge phases are spliced.
Normalization follows retrieval. Cache hits make no MPI call, while each miss
uses the existing collective `disp` on all ranks. The per-contour evaluation
budget continues to count **logical** requests, including hits.

The `.nyquist.summary` adds `growth_spacing` and
`growth_original_gamma_bounds`. Existing initial-contour and search status scopes
remain distinct. The first 13 `.nyquist.growth` data columns retain A3 order;
`anchor`, `N_above` and `spacing` are appended. A comment records original gamma
bounds and mode. Linear rows correctly say `current_upper`; log rows say
`original_lower` or `original_upper`. Readers should skip all `#` comment lines.

## 3. Linear/logarithmic validation

Use the established A2/A3 equilibrium at `(0.001,0.3)`, real bounds `[-1,1]`,
top gamma 0.2 and absolute tolerance 1e-4. Independent ordinary roots remain near
`±0.39666+0.085986i`. Complete legacy root checks reproduce them; their printed
frequency precision is respected.

| Lower gamma | Method | Initial W | Steps | Final bracket | Width |
|---:|---|---:|---:|---|---:|
| 0.01 | Omitted/default linear | 2 | 11 | `[0.0859814453125,0.08607421875]` | 9.27734375e-5 |
| 0.01 | Explicit linear | 2 | 11 | Same as default | Same |
| 0.01 | Nearest-edge log | 2 | 12 | `[0.0859334813901669,0.0859963544046832]` | 6.28730145e-5 |
| 0.0001 | Linear | 2 | 11 | `[0.085896923828125,0.08599453125]` | 9.76074219e-5 |
| 0.0001 | Nearest-edge log | 2 | 13 | `[0.0859683182979603,0.0860481205488949]` | 7.98022509e-5 |

All successful brackets contain the independent growth near 0.085986. Omitting
spacing preserves the original 5005 physical / 8690 logical requests as well as
its bracket and iteration count. The log searches need not reproduce the linear
bracket or number of steps.

The first log trial for lower=0.01 illustrates the complement explicitly:
`M=0.0447213595499958`, lower strip `[0.01,M]`, W_strip=0 but **N_above=2**.
The next trial uses original upper `[0.0945741609003176,0.2]`, W_strip=0 and
N_above=0. Both anchors are exercised without confusing the two counts.
With lower=1e-4, this reference's accepted trials all select the original upper
boundary, making its sampled geometry identical to the full-height baseline.

For **every accepted log trial**, the suite independently evaluates an ordinary,
uncached rectangle from M to the original upper boundary. Integer counts agree
exactly and converged phases agree with `2*pi*N_above` to 1e-8 radians. This
includes successful, iteration-limited and pre-failure trials: 30 such rectangles
per build for the reference/control suite, plus 11 for the mirror fixture. These
ordinary calculations have zero cache hits; all actual rank call counts agree.

Additional outcomes:

| Case | Outcome |
|---|---|
| Isotropic reference control | Initial W=0; zero iterations, no bracket/estimate. |
| Raised lower=0.1 | Initial W=0; known growing pair excluded, no estimate. |
| Iteration limit 1 | Status 4; retains `[0.04472135955,0.2]`, no successful estimate. |
| Depth 8, tolerance 1e-14 | Attempt 5 fails after four accepted updates; growth status 3/A1 status 5, retains approximately `[0.07842564913,0.09457416090]`. |
| Unsupported spacing | Specific input error before determinant sampling. |
| Gamma boundary zero | Specific positive-coordinate input error, all-rank abort. |
| Adjacent representable positive endpoints | Direct coordinate-helper check returns nonprogress status 6. |
| Nonfinite/nonpositive helper coordinate | Rejected before taking a logarithm. |

The standalone new math driver passes **58 checks** in debug/full-trap and -O3
builds. It includes three polynomial zeros at growth rates 0.00011, 0.002 and
0.003 in an original `[0.0001,0.2]` domain: a weak root near the lower edge,
distinct rates, both anchors, and independently computed full-height counts.
It also tests degenerate growth, iteration/refinement limits, invalid spacing,
adjacent endpoints and the safe geometric mean over `[1e-300,1e300]`. A1/A3
mathematical suites retain their 306 and 10056 checks, respectively.

## 4. Caching and performance

One serial measurement series, four MPI ranks, same physical equilibrium and
sampling controls; samples disabled. Cache misses equal physical collective
frequency evaluations. The full-height comparison is a validation-only copied
search module that always selects the original upper anchor while retaining
caching and the **same log trial points**. It is not a user-facing switch.

| Original lower | Strategy | Steps | Logical requests | Physical/misses | Hits | Refinements | Max depth | Global passes | Determinant seconds | Calculation seconds |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | Linear shrinking strips | 11 | 8690 | 5005 | 3685 | 1657 | 16 | 36 | 0.243736 | 0.254973 |
| 0.01 | Log nearest-original edge | 12 | 8374 | 5114 | 3260 | 1275 | 15 | 39 | 0.158566 | 0.166875 |
| 0.01 | Log cached full height | 12 | 8264 | 5112 | 3152 | 1220 | 15 | 39 | 0.149599 | 0.157144 |
| 0.0001 | Linear shrinking strips | 11 | 9108 | 4853 | 4255 | 1866 | 15 | 36 | 0.154799 | 0.165311 |
| 0.0001 | Log nearest-original edge | 13 | 8976 | 5584 | 3392 | 1352 | 14 | 42 | 0.162684 | 0.172633 |
| 0.0001 | Log cached full height | 13 | 8976 | 5584 | 3392 | 1352 | 14 | 42 | 0.154545 | 0.164672 |

`MPI_Wtime` calculation duration is the maximum rank time through worker
termination, including the initial rectangle, cache/sampling/bisection work. It
excludes equilibrium preparation, process launch and subsequent history writing.
The existing driver-elapsed field separately includes output overhead. Determinant
time sums actual miss communication and `disp` calls on rank 0. These short
one-series timings vary with host scheduling and frequency-dependent evaluation
cost; call/refinement counts provide the more reproducible algorithm comparison.

**Nearest-edge is slightly worse than cached full height for lower=0.01**:
110 extra logical requests, two extra physical calls and 55 extra refinements.
For lower=1e-4 the work is identical; timing differences do not imply an
algorithmic improvement. No alternate tie policy or sampling optimization was
introduced to favor the benchmark.

Comparing log with linear also includes different numerical convergence paths:
12 versus 11 steps for lower=0.01, and 13 versus 11 at the default lower boundary.
For the 0.01 case, excluding the common 630 logical initial requests, average
logical requests per strip are about 645.3 for log and 732.7 for linear. Lower
per-contour cost therefore does not guarantee lower total cost or fewer physical
calls when more trials are required. No general speedup is claimed, particularly
for numerical-VDF backends that remain outside this validation.

## 5. Independently identified oblique mirror branch

Physical exploration began only after the logarithmic mathematical/reference and
established regression checkpoint passed. The selected single wavevector and
equilibrium were retained after actual root/polarization results, not treated as
predetermined valid fixtures:

| Parameter | Value |
|---|---|
| `(kperp*d_p,kpar*d_p)` | `(0.3,0.05)`, obliquity about 80.54 degrees |
| Proton beta_parallel | 10 |
| Proton Tperp/Tparallel | 2 |
| Electron beta_parallel, anisotropy | 10, 1 |
| Densities, drifts | Both density 1; both drift 0 |
| Charges/masses | Proton 1/1; electron −1/5.44617021487e-4 |
| vA/c | 1e-4 |
| Backend | Both `use_bM=T`, `relat=F`; analytic warm response |
| Harmonic controls | Cap 500, Bessel cutoff 1e-50 |

The unmodified ordinary secant solver, not the Nyquist search, found the growing
mode. Initial guesses were `0.001+0.04i`, `-0.001+0.04i`, and `-0.005+0.08i`;
all converge to the **same** root, not three distinct modes. Its three observed
frequencies are approximately:

- `+5.86956e-17 + 0.038169728624597589i`;
- `-6.88102e-17 + 0.038169728624597575i`;
- `+1.18374e-17 + 0.038169728624597631i`.

They converge in 14, 12 and 14 iterations using `D_threshold=1e-18`, iteration
cap 200 and the unchanged secant method 0. Their raw determinant residuals are
below 1e-18. A copied validation main calls an observer after ordinary root solving;
production root/susceptibility routines are not patched. The observer collectively
refreshes `disp(root)` then calls existing `calc_eigen` with eigen diagnostics on
and heating off. It records high-precision roots and fields in a distinct
`.mirror_reference` file.

Define magnetic compressibility
`C_B=|delta B_z|^2/(|delta B_x|^2+|delta B_y|^2+|delta B_z|^2)` and relative tensor
residual `||wave*E||_2/(||wave||_F*||E||_2)`. The latter guards against relying
only on a small frequency-scaled determinant. The proton density/parallel-field
phase correlation is `Re(delta n_p*conjg(delta B_z))/(|delta n_p|*|delta B_z|)`.

| Proton anisotropy | Ordinary Im(z) | Max abs(Re(z)) over seeds | C_B | Relative tensor residual | Density/B_parallel correlation |
|---:|---:|---:|---:|---:|---:|
| 2 | +0.0381697286246 | 6.89e-17 | 0.829507497 | <=6.85e-17 | −1 |
| 1.5 | +0.0184096806751 | 2.93e-16 | 0.909057310 | <=4.01e-16 | −1 |
| 1 | −0.0196497825507 | 5.43e-17 | 0.759338477 | <=4.85e-16 | −1 |

The magnetic field also satisfies `k dot delta B=0` to the measured relative
roundoff level (below 1e-10, typically about 1e-17 or smaller). Increasing the
harmonic cap to 1000 and reducing the Bessel cutoff to 1e-60 reproduces the growing
root/polarization to the printed precision; tests allow 1e-9 in growth and 1e-10
in relative tensor residual.

The combination of a high-beta, strongly oblique, nearly aperiodic mode,
dominantly parallel magnetic fluctuations, density anticorrelation, proton
perpendicular-anisotropy drive and stabilization under its removal strongly
supports identification as the proton mirror branch. This is an inference from
multiple diagnostics, not from positive growth or small real frequency alone.
Mirror theory treats temperature-anisotropy-driven compressional responses;
see [Hellinger (2007)](https://doi.org/10.1063/1.2768318) and
[Pokhotelov et al. (2002)](https://doi.org/10.1029/2001JA009125).
We do not use their asymptotic threshold formulas as a quantitative growth
reference at this finite perpendicular wavevector. Full perturbed-pressure balance
was not calculated; density anticorrelation is supporting evidence, not a claim
that the full pressure-balance relation was measured. No separate-code numerical
cross-validation is claimed; the independent method here is ordinary root solving
with existing polarization diagnostics.

One exploratory guess approached a very small near-origin frequency rather than
the well-resolved branch. It was not adopted as a physical reference: a small raw
scaled determinant alone is insufficient. The selected roots have much stronger
tensor residual checks, and the Nyquist lower boundaries exclude the origin and
that near-origin iterate.

## 6. Mirror Nyquist checks and growth estimates

All counts below are measured, signed CCW windings of the raw determinant.

| Case | Real interval | Gamma interval | W | Calls |
|---|---|---|---:|---:|
| Selected mirror box | `[-0.05,0.05]` | `[0.005,0.08]` | 1 | 468 |
| Same, stricter sampling | Same | Same | 1 | 922 |
| Expanded box | `[-0.2,0.2]` | `[0.0001,0.15]` | 1 | 566 |
| Isotropic control, same selected box | `[-0.05,0.05]` | `[0.005,0.08]` | 0 | 474 |
| Raised lower boundary | `[-0.05,0.05]` | `[0.05,0.08]` | 0 | 452 |

The stricter case uses 16 initial intervals, phase limit 0.35, curvature 0.075 and
phase/integer tolerances 1e-10. Every count converges under A1's complete checks.
Thus no extra unstable root is indicated inside either prescribed unstable box,
subject to the numerical convergence assumptions; this is not an inventory outside
these frequency/growth bounds. No count was doubled from the nearly zero signed
real frequencies of repeated root guesses.

Both searches use the selected box and absolute growth tolerance 1e-4:

| Method | Steps | Accepted bracket | Width | Logical | Physical/misses | Hits | Refinements | Depth | Passes | Determinant seconds | Calculation seconds |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Linear | 10 | `[0.03810546875,0.0381787109375]` | 7.32421875e-5 | 5656 | 3759 | 1897 | 364 | 10 | 33 | 0.523911 | 0.528348 |
| Log | 11 | `[0.0381488779615716,0.0382005590054078]` | 5.16810438e-5 | 5836 | 3900 | 1936 | 230 | 9 | 36 | 0.527480 | 0.532226 |

Both brackets contain the independently measured growing mirror root. Log trials
have less local refinement here but need another complete contour; total logical
and physical calls are higher. These are single-wavevector pilot timings, not
a performance claim for arbitrary distributions. All eleven accepted mirror log
trials additionally agree with independent uncached full-height count calculations.
The isotropic log search returns a valid initial zero count and no estimate.

## 7. Regression execution, provenance and limitations

Complete suites were run on the final implementation in both configurations:

```sh
bash validation/run_nyquist.sh
bash validation/run_nyquist_growth.sh
bash validation/run_nyquist_log.sh
bash validation/run_nyquist_mirror.sh
```

Repeat with `FFLAGS='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow'` for debug.
Optimized defaults are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`.
The A2 runner includes the complete unchanged A1 suite (306 checks in each of its
debug/-O3 builds); A3 retains 10056 math/cache checks and all its physical/failure
cases. New log math passes 58 checks, also separately compiled with `-std=f2008`,
full invalid/zero/overflow traps and with -O3/warnings. Mirror root, polarization,
control, strict/expanded count and bracket tests pass in debug and optimized builds.

Existing legacy root, stable-root, mapping, consecutive-scan and plane-scan output
comparisons also passed: nine scientific files remain byte-identical to the
pre-extension executable in both suite configurations. The A3 returned-NaN,
budget, cache-hit, invalid-input and coordinated-shutdown tests remain intact.
The initial log checkpoint passed before mirror physical exploration began.

GNU Fortran 11.4.0, Open MPI 4.1.2, four local ranks with `--oversubscribe`, existing
BLAS/LAPACK and one library/OpenMP thread. A clean parallel Automake build with the
existing production -O4 flags also passed the logarithmic reference test. No
in-tree executable or protected scientific output was overwritten; all generated
results/builds used fresh external scratch directories.

Final complete log suite artifacts:

- Optimized: `/tmp/alps-nyquist-log.gJp9E2fL`, executable SHA-256
  `7e6266726a819cd851a38842ad6b1600e528ce09343c483ed3624c3d2db47e74`.
- Debug: `/tmp/alps-nyquist-log.J7l7U0DB`, executable SHA-256
  `1b56b946551a99853cce9ca6fa8a7110a7a504ddbb85cbf8306f941c12618904`.
- Cached full-height optimized validation variant:
  `ae8450d75da1eb6e83df720ffc1948275af6a2f86c35b55e08990d98efc513bd`.

Mirror artifacts: `/tmp/alps-nyquist-mirror.xZzIzWl6` (optimized) and
`/tmp/alps-nyquist-mirror.0VTGdvgV` (debug), with executable hashes
`7defc5cd0b93e4829a1bbf1c0077fd4d9c7bd6ccce1e1aed7523db7e182ce388` and
`a41c83b2f2d9221d20a130aee4c582fb9e5c8dd0dd8df33a570a34ef2d9d7cc2`.
Each build has a source/executable manifest; copied observer/full-height variants
have separate provenance. Final matrix logs are under
`/tmp/alps-nyquist-log/verified/`. Temporary artifacts are not durable archives;
the new Fortran tests, input fixtures and measured tables here retain the evidence.

Existing limitations remain: finite sampling cannot rigorously exclude arbitrary
unsampled loops, roots too close to boundaries can be unresolved, and determinant
accuracy is separate from phase-sampling convergence. The A2 high-frequency
conditioning issue and unpropagated WOFZ overflow flag were not modified. No cold,
relativistic or numerical-VDF Nyquist validation is claimed. All growth estimates
are bounded by the original user domain, not certified whole-plane maxima.

Use linear spacing as the unchanged default and log spacing when its different
trial placement is useful. Do not assume the nearest-edge strategy lowers cost;
measure refinements and misses at the actual equilibrium. Any later edge-phase
reuse or sampling optimization requires separate convergence validation. Future
wavevector-line/plane work must still update context collectively, reset the
single-wavevector cache and provide scan-specific outputs. No option 4/5 scan
capability was implemented in this extension.

## 8. Final A3 update: relative tolerance for logarithmic searches

The same finite positive `nyq_growth_tolerance` now means strict absolute width
`H-L < epsilon` for linear spacing and strict relative width
`(H-L)/(0.5*L+0.5*H) < epsilon` for log spacing. The default remains `1.d-4`;
in log mode this requires relative width below 0.01%. The arithmetic midpoint
avoids forming `L+H`, and remains an estimate rather than a dispersion root.

Only `src/ALPS_nyquist_growth.f90` (shared width calculation and the two stopping
checks) and `src/ALPS_nyquist_driver.f90` (additive summary fields) changed in
production source for this update. `growth_tolerance_type` records `absolute`
or `relative`; `growth_relative_width` supplements the unchanged absolute
`growth_width` for valid brackets. History columns are unchanged. Geometric
trials, original-edge selection, complement counting, exact raw caching, MPI,
failure statuses, susceptibility, root solver and A1 sampling are unchanged.

Measured final brackets, with tolerance `1.d-4` in the selected mode:

| Case / spacing | Final bracket | Relative width | Steps | Logical requests | Physical calls | Hits |
|---|---|---:|---:|---:|---:|---:|
| Reference, linear, lower 0.01 | [0.0859814453125, 0.08607421875] | 1.078411897e-3 | 11 | 8690 | 5005 | 3685 |
| Reference, log, lower 0.01 | [0.08598063183978559, 0.08598849276288539] | 9.142249364e-5 | 15 | 10606 | 6490 | 4116 |
| Reference, log, default lower 1e-4 | [0.08598327558047909, 0.08598826191959801] | 5.799028364e-5 | 17 | 11956 | 7414 | 4542 |
| Mirror, linear | [0.03810546875, 0.0381787109375] | 1.920245791e-3 | 10 | 5656 | 3759 | 1897 |
| Mirror, log | [0.03816825015433743, 0.03817147980940077] | 8.461269289e-5 | 15 | 8098 | 5356 | 2742 |

The reference brackets enclose the independently measured growth near 0.085986.
The mirror bracket encloses the ordinary-solver growth 0.0381697286245976 at
`(kperp,kpar)=(0.3,0.05)`. Initial windings remain two and one respectively.
The linear reference reproduces its previous bracket, 11 steps and evaluation
counts exactly. Linear and log tolerances now have different meanings: these
rows are not equal-accuracy performance comparisons.

Compared with the historical absolute-tolerance log searches above, the reference
with lower 0.01 adds 3 steps and 1376 physical calls; the default-lower reference
adds 4 steps and 1830 calls; the mirror adds 4 steps and 1456 calls. Final maximum
refinement depths are 17, 17 and 13 respectively. No A1 control was relaxed.
The cached full-height log baseline still agrees at every trial and final bracket:
6488 calls versus 6490 for nearest-edge strips with lower 0.01, and identical
7414 calls with the default lower bound. No new optimization was introduced.

The synthetic weak root at gamma=0.00011 needs 15 steps and gives
`[0.00010999725209639309, 0.00011000498180012830]`, relative width
`7.026932044e-5` and absolute width approximately `7.73e-9`. It therefore does
not stop at absolute width 1e-4. The distinct-root analytic test uses relative
tolerance 1e-4 (a tighter absolute result than its former absolute 1e-6 test).
Requesting relative 1e-6 on that same domain exhausts the unchanged default A1
refinement depth: this is now an explicit unresolved-strip test retaining the
last valid bracket, not a successful growth estimate. Tests also check strict
inequality at exact tolerance equality and safe midpoint arithmetic near 1e308.

The complete A1/A2/A3, log and mirror suites passed in optimized and debug builds
using the configurations in Section 7. A1 retains 306 checks per debug/-O3 build,
A3 retains 10056 checks, and log math now has 72 checks. Four-rank accounting,
uncached full-height counts above every accepted log trial, stable/threshold
controls, iteration/refinement/evaluation-budget failures, invalid coordinates
and tolerances, nonfinite determinant injection and coordinated shutdown passed.
The nine legacy scientific files remain byte-identical to the pre-update build
in both configurations. The A1 source SHA-256 is unchanged:
`d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.

All runs used isolated temporary directories. Logs and the pre-update source/build
snapshot are under `/tmp/alps-nyquist-relative/`. Final regression artifacts:

| Suite | Optimized | Debug |
|---|---|---|
| A1/A2, including legacy comparisons | `/tmp/alps-nyquist-a2.Gf51uEWd` | `/tmp/alps-nyquist-a2.JVJ1xaCR` |
| A3 | `/tmp/alps-nyquist-a3.ExYqFdMz` | `/tmp/alps-nyquist-a3.uebbhccv` |
| Log | `/tmp/alps-nyquist-log.Q8emVsal` | `/tmp/alps-nyquist-log.MX2VF8gn` |
| Mirror | `/tmp/alps-nyquist-mirror.OgSXeNC7` | `/tmp/alps-nyquist-mirror.gyCbtj0L` |

Build manifests identify sources, flags and executable hashes. Tight relative
tolerances can require more boundary refinement or fail near roots; failures
retain the bracket without an estimate. All estimates remain restricted to the
original frequency domain and supported warm analytic backend. A3 is ready for
B1; no wavevector scans or new input controls were implemented in this update.
