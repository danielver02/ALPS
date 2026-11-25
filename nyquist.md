# ALPS Nyquist Stability Analysis

This guide describes the Nyquist contour methods, growth-rate estimates,
wavevector scans, output formats, numerical interpretation, and validation
for the current ALPS implementation. For the concise namelist parameter/default
reference and general ALPS input options, see [input.md](input.md).

## Operation selection

**`scan_option`**  
Select operation (default 1):

- 1: Consecutive scans along input paths in wavevector space,  
- 2: Double scan over wavevector plane.
- 3: Nyquist calculation at the input `kperp,kpar`; bypass root solving and scans.
- 4: Nyquist calculations along the existing `scan_input_N` wavevector paths.
- 5: Nyquist wavevector-plane scan using exactly two different `scan_input_N` directions.

**`n_scan`**  
Number of wavevector scans.  
0 turns off wavevector scans.  
Must be 1 or larger for `scan_option`=1.  
Must be set to 2 for `scan_option`=2.  
For `scan_option=3`, `n_scan`, `nroots`, and `use_map` are ignored; no guesses,
maps, or scan instructions are required. Options 1 and 2 retain their existing
root-search/map behavior, including when `n_scan=0`.
Option 4 requires `n_scan>=1`; it ignores `nroots,use_map`, reads species and
`scan_input_1` through `scan_input_n_scan`, and bypasses guesses, maps and root tracking.
Option 5 requires `n_scan=2`; it likewise bypasses guesses/maps/root tracking.
The two scan types must differ and both must be in 1–4 (type 0 is unsupported).

## Nyquist method and controls

### Nyquist configuration

Required for `scan_option=3`, `4` or `5`. All species require `relat=F`, positive
density/mass, nonzero charge, finite applicable parameters, `vA>0`, and `kpar /= 0`.
Analytic species (`use_bM=T`) require positive beta and anisotropy and finite
drift. Options 3/4/5 also accept nonrelativistic numerical species (`use_bM=F`)
and mixed equilibria through the existing production numerical susceptibility.
Maximum-growth searches continue to require warm analytic species.
Cold analytic and relativistic Nyquist responses remain unsupported; existing
root/scan physics is unchanged.
Create the usual `solution/` directory before running.

```fortran
&nyquist
 nyq_geometry='rectangle'
 nyq_omega_min=-1, nyq_omega_max=1
 nyq_gamma_min=0.01, nyq_gamma_max=0.2
 nyq_sampling%initial_intervals=8
 nyq_sampling%max_evaluations=200000
 nyq_save_samples=F
/
```

The `&nyquist` parameters and the defaults for each `nyq_sampling` control
are listed in [input.md](input.md).

For example, `nyq_sampling%max_depth=32` permits additional local refinement.
The kernel also requires global sampling agreement, including an independently
seeded pass. A small integer residual alone is insufficient. Finite sampling is
numerical evidence, not proof that arbitrary phase loops cannot be missed.

Rectangles use the raw determinant by default and report a signed accepted
winding only after convergence. The warm analytic response has no upper-half-plane
poles, so this counts roots in the stated finite rectangle. Numerical VDFs return
the accepted winding of the computed boundary trajectory; its physical root-count
interpretation additionally requires representation and continuity convergence.
It does not count all
unstable roots outside those bounds.

Horizontal paths always use `D_ALPS/(z*vA)**6`, implemented in the adapter with
log-polar arithmetic. Supply the same frequency bounds and `nyq_gamma_min`, omit
`nyq_gamma_max`, and set `nyq_geometry='horizontal'`. The result is an unrounded
finite-interval phase. **Sampling convergence does not establish asymptotic
closure or an integer instability count.** Repeated interval expansion and
independent closure assessment are required; no automatic frequency cutoff is used.

Summary output goes to the console and `solution/<runname>.nyquist.summary`.
Optional samples go to `solution/<runname>.nyquist.samples`, with columns
`Re(z) Im(z) Re(F) Im(F) cumulative_phase`. Rerunning the same name replaces
these Nyquist files; ordinary root/map/scan files are untouched. Rectangles with
retained samples also report the bottom-edge and other-three-edge phases.
`nyq_normalize_rectangle=T` retains these diagnostics without requiring a sample file.

Summaries contain geometry, normalization, equilibrium parameters, status,
phase, accepted winding or `NONE`, per-rank call-count checks, MPI timing, and
refinement diagnostics. `wall_seconds` measures the collective calculation after
initialization; sample writing and total driver time before summary writing are
reported separately. Failed sampling reports an unaccepted latest path or
partial prefix; its endpoints may not reach the requested bounds.

Status codes preserve A1: 0 converged, 1 invalid controls/path, 2 zero/nonfinite
value, 3 evaluator failure (adapter code 101: unrepresentable normalized amplitude),
4 evaluation cap, 5 refinement/coordinate-resolution cap, 6 pass/integer disagreement.
Sampling/output failures return a nonzero process exit after collective shutdown;
invalid input uses communicator-wide abort, consistent with existing fatal input
handling. Unresolved never means stable.

Complete example inputs are in `tests/data/nyquist/`. Build and regression
instructions and measured limitations are in `NYQUIST_A2_RESULTS.md`.

### Numerical-VDF calculations and scans

Use `scan_option=3`, `4` or `5`, the normal numerical species inputs and VDF files, and
`nyq_find_max_growth=F`. No numerical-VDF flag or separate susceptibility is used.
The equilibrium is read, differentiated, fitted and transferred once through the
ordinary startup. Harmonic assignments and worker Bessel arrays are initialized
before the first collective determinant call. Do not omit continuation inputs
simply because a contour lies strictly above the real axis.

The numerical grid must be finite, uniform and rectangular, with increasing
momentum axes and nonnegative perpendicular momentum and VDF values. Logarithmic
continuation with methods 1/2 requires strictly positive samples. Method 0 retains
its supplied analytic-continuation semantics. The VDF cannot be identically zero.
Integration/resonance controls
retain their ordinary meaning. No new numerical integration, differentiation,
analytic continuation or determinant normalization is introduced.

For these Nyquist modes, coordinate validation occurs before derivatives and fits.
Each written decimal coordinate defines a rounding interval of half a unit in its
last recorded decimal place, with a small IEEE arithmetic allowance. A common
positive-spacing affine axis must fit every coordinate interval, repeated row,
endpoint and local increment. Genuinely nonuniform or nonmonotonic grids are
rejected. Accepted rounded axes are replaced **in memory only** by that uniform
model; VDF values and indices are untouched. Already uniform high-precision axes
retain their stored coordinate bits. Canonical coordinates are used for all
derivatives, continuation fits and harmonic/Bessel initialization and broadcast
to every rank. `NYQUIST_GRID` startup diagnostics report the inferred spacing,
residual and serialization bound. This is not VDF interpolation or normalization.

Numerical scans retain the existing harmonic-based MPI decomposition, initialize
the equilibrium once, rebuild harmonic assignments/Bessel arrays when `kperp`
changes, and use fresh per-wavevector determinant caches. Pure `kpar` changes
retain the Bessel arrays. Exact `kperp=0` is supported by the existing parallel
tensor limits. Numerical-VDF winding remains subject to representation and
integration-continuity validation; enabling scans does not certify every measured
distribution as a physical analytic dispersion function.

The summary identifies each species as `analytic_bimaxwellian` or
`numerical_VDF_validation`. Analytic physical fields are retained. Numerical
species instead report density/charge/mass, source filename, grid intervals and
support, continuation settings, harmonic limit, `Tlim`, Bessel cutoff, resonance
window controls and the ordinary density/current integrals. Unused `bM_*`
defaults are not reported as properties of numerical species.

Run `bash validation/run_nyquist_vdf.sh` in a system with Fortran, MPI and the
existing BLAS/LAPACK dependencies. It creates isolated builds and generates VDFs
with the unmodified `generate_distribution` utility. Runnable input examples and
generation instructions are in `tests/data/nyquist_vdf/`; measured comparisons
and limitations are in `NYQUIST_C1_RESULTS.md`. Matched-frequency probes use a
copied validation adapter around the same collective evaluator, not a second
production determinant implementation.

Sampling convergence alone cannot detect a small discontinuity in the computed
susceptibility. Check momentum resolution and support independently and inspect
the resonance-treatment transitions before interpreting a numerical winding as
a physical root inventory. C1 deliberately leaves those integration routines
unchanged.

The subsequent [C2b correction](NYQUIST_C2B_RESULTS.md) uses continuous nodal
interpolation and consistent fractional-cell quadrature in the nonrelativistic
numerical response. It applies to ordinary numerical-VDF roots as well as
Nyquist evaluations; finite-grid numerical results can change. The original
interpolation remains callable as `funct_g_legacy` for validation. Run
`bash validation/run_nyquist_c2b.sh` for the correction-specific tests. The
`Tlim` finite-window mismatch and the legacy support-edge fallback remain
limitations; interior continuity does not establish global holomorphy.

The [C2c method study](NYQUIST_C2C_DESIGN.md) compares an exact finite-window
linear core with a fixed-cell Cauchy prototype in isolated validation builds.
Production remains C2b; no new selector is available. The study identifies
static tensor-cancellation and damped-residue consistency requirements before
adopting a replacement. Reproduction instructions and measured data are in
`validation/c2c/` and `tests/data/nyquist_c2c/`.
For the reference cyclotron case, the corrected coarse 32×64 grid has computed
winding 3, while the 64×128 and 96×192 matched trajectories approach the
analytic winding 2. The coarse value is a numerical convergence warning, not a
physical third-mode identification; see the C2b report and adaptive checks.

#### Maximum-growth estimation with closed strips

This optional calculation uses the same `scan_option=3`, input wavevector and warm
analytic restrictions. Explicit `nyq_omega_min`, `nyq_omega_max` and (for a
rectangle) `nyq_gamma_max` remain mandatory; none is inferred from roots or scans.
For example:

```fortran
&nyquist
 nyq_geometry='rectangle'
 nyq_omega_min=-1, nyq_omega_max=1, nyq_gamma_max=0.2
 ! nyq_gamma_min defaults to 1.d-4
 nyq_find_max_growth=T
 nyq_growth_tolerance=1.d-4
 nyq_growth_max_iterations=32
/
```

The full original rectangle is evaluated first. An accepted initial zero count
means no roots above the selected lower boundary **within this frequency
domain**; no growth estimate or artificial bracket is produced. A positive count
starts a bracket `[L,H]=[nyq_gamma_min,nyq_gamma_max]`. Here the initial upper
boundary bounds the search domain, not all modes outside it.

In default `'linear'` mode each trial `M=(L+H)/2` evaluates the closed strip from `M` to the **current H**,
over the original real-frequency bounds. A positive accepted winding sets `L=M`;
zero sets `H=M`. A change directly from 2 to 0 is supported. Every strip uses
the unchanged A1 adaptive algorithm, stricter global passes and independent seed.
Negative/inconsistent counts or a failed strip stop the search and preserve the
last valid bracket; an unresolved count is never treated as zero.

For `nyq_growth_spacing='log'`, use the geometric trial
`M=exp(0.5*log(L)+0.5*log(H))`, with finite strictly positive boundaries.
Compare the logarithmic distances from M to the **original** gamma_min and
gamma_max, rather than to the changing bracket. A closer original upper boundary
selects `[M,gamma_top]`, whose accepted winding directly gives `N_above`.
A closer original lower boundary selects `[gamma_0,M]`, whose winding counts
roots below M; then `N_above=N0-W_strip`, using the accepted original count N0.
Exact equal distances in working precision choose the original upper anchor;
near ties follow the computed distances deterministically. Set `L=M` for positive
N_above and `H=M` for zero. Negative, excessive or inconsistent counts are rejected,
never clipped. Every chosen strip still receives the complete A1 convergence test.

Linear spacing retains the absolute stopping rule `H-L < nyq_growth_tolerance`.
Log spacing uses `(H-L)/(0.5*L+0.5*H) < nyq_growth_tolerance` instead.
Thus a logarithmic tolerance of `1.d-4` requires relative bracket width below
0.01%. The midpoint in this criterion is arithmetic; trial points remain geometric.
For example, `nyq_growth_spacing='linear', nyq_growth_tolerance=1.d-4`
requires absolute width below `1.d-4` in normalized frequency, while
`nyq_growth_spacing='log', nyq_growth_tolerance=1.d-4` requires width below
`1.d-4` times the arithmetic midpoint. The two modes need not take the same
number of steps and their tolerance values have different meanings.
A geometric trial that rounds to an endpoint
returns coordinate-resolution status 6 without changing the last valid bracket.
Neither original-edge anchoring nor logarithmic splitting guarantees fewer
evaluations. Near-root edges may need more refinement; cache/normalization rules
below are unchanged.

The reported midpoint is a bracket estimate of the fastest enclosed growth,
not a dispersion root. It remains the arithmetic midpoint of the final bracket;
the spacing option changes the trial points, not this output convention.
Success requires the selected absolute or relative width criterion before
the iteration cap. No successful estimate is emitted on a failed or limited
search. The result excludes no roots outside the original real-frequency or
upper-growth bounds, or below the lower threshold. Horizontal maximum-growth
requests and invalid tolerances/iteration limits are rejected.

A rank-0 cache persists across the initial contour, strips and all A1 passes.
It matches exact complex frequencies, stores raw determinants, and applies any
requested normalization after lookup. Hits require no MPI communication; misses
call the existing collective `disp` on every rank. The cache is cleared for each
operation and is not reused at another wavevector. Completed edge phases are
not spliced: A1 still unwraps and validates each complete closed strip. Values
on unchanged edges can be reused, while different floating-point coordinates
remain separate cache entries. No interpolation or approximate matching occurs.

The A1 `max_evaluations` limit remains a **logical per-contour budget**, including
cache hits. Caching does not relax the sampler's limits or convergence tests.
The total search can exceed that limit across multiple contours; the iteration
cap bounds the number of additional strip calculations. Memory grows with unique
cached frequencies.

The existing summary/sample fields describe the **initial contour**. Additional
fields report `growth_status`, `growth_initial_winding`, `growth_converged`,
`growth_no_roots`, `growth_bracket_valid`, bracket/width, midpoint or `NONE`,
iterations/attempts, total passes/refinements/depth, `logical_requests`, actual
`disp_calls`, cache hits/misses/fraction and `determinant_seconds`. `wall_seconds`
includes initial sampling and the entire bisection; it excludes equilibrium
initialization and process launch. `driver_seconds_before_summary` also includes
history/sample writing. The separate A1 `status` still describes the initial
contour; **use `growth_status` for the overall search**.

`growth_spacing` and `growth_original_gamma_bounds` identify the search mode and
fixed domain. `growth_tolerance_type` is `absolute` for linear and `relative` for
log spacing. `growth_width` retains the absolute width; `growth_relative_width`
additionally reports width divided by the arithmetic midpoint for valid brackets.
The first 13 growth-history columns retain the A3 layout; appended
columns are `anchor`, `N_above`, and `spacing`. Linear rows use `current_upper`;
logarithmic rows use `original_lower` or `original_upper`. The history also
includes a comment with original gamma bounds and spacing.

| `growth_status` | Interpretation |
|---:|---|
| 0 | Converged bracket within the specified frequency domain. |
| 1 | Valid initial zero winding; no search or estimate needed. |
| 2 | Initial contour unresolved; no accepted count or bracket. |
| 3 | Strip failed; last valid bracket retained, no successful estimate. |
| 4 | Iteration limit; last valid bracket retained. |
| 5 | Negative or inconsistent count; unresolved search. |
| 6 | Growth-coordinate resolution exhausted; last valid bracket retained. |
| 7 | Search storage allocation failed. |

For strip/initial numerical failures `growth_failure_A1_status` gives the A1
reason. The compact `solution/<runname>.nyquist.growth` file records trial,
strip bounds, winding or `NONE`, A1 status, updated bracket and sampling work for
each attempted iteration. It has only a header when no bisection is attempted.
`nyq_save_samples=T` writes the initial contour only, not a large file for every
strip. Operation exit is nonzero for failed/limited searches, while a valid
initial zero count completes normally.

Example inputs and failure cases: `tests/data/nyquist_growth/`. New tests:
`bash validation/run_nyquist_growth.sh`. Results and caching benchmarks:
`NYQUIST_A3_RESULTS.md`. Option 4 extends these calculations to wavevector paths;
option 5 evaluates wavevector planes using the same per-point calculation.
Logarithmic and oblique mirror-mode fixtures are in `tests/data/nyquist_log/` and
`tests/data/nyquist_mirror/`; validation and measured performance are documented
in `NYQUIST_A3_LOG_RESULTS.md`. Test entry points are
`bash validation/run_nyquist_log.sh` followed by
`bash validation/run_nyquist_mirror.sh`.

### Nyquist wavevector paths (`scan_option=4`)

Use the same species input and `&nyquist` controls as option 3, together with
`n_scan>=1` and the existing indexed `&scan_input_N` groups. Rectangle geometry
is required because the scientific output reports accepted full-box windings.
Horizontal diagnostics remain available in option 3. Frequency bounds and
sampling/growth controls are read once and remain fixed for the entire path.

For example, with the warm mirror species described in `NYQUIST_A3_LOG_RESULTS.md`:

```fortran
&system
 ! Supply the usual species/allocation controls as well.
 scan_option=4, n_scan=1, kperp=0.2, kpar=0.05
/
&nyquist
 nyq_omega_min=-0.05, nyq_omega_max=0.05, nyq_gamma_max=0.08
 nyq_find_max_growth=T, nyq_growth_spacing='log', nyq_growth_tolerance=1.d-4
/
&scan_input_1
 scan_type=3, swf=0.4, swlog=F, ns=2, nres=1
/
```

The initial wavevector is output once, followed by `ns` points per segment,
including its endpoint: `1+sum(ns)` rows. Coordinates use the ordinary scan's
arithmetic at tracking indices `nres,2*nres,...,ns*nres`. Intermediate tracking
steps are not evaluated. Consecutive segments start at the preceding endpoint;
the shared boundary is not output twice and keeps the index of the segment that
ended there. Identical explicitly requested output
coordinates (for example a zero-length segment) still receive their own rows.
`swlog` controls wavevector spacing independently of growth spacing. `nres`
must be positive but creates no extra Nyquist evaluations; `heating` and `eigen`
are root diagnostics and are not used here. Omitted `nres` defaults to 1 and
omitted `swlog,heating,eigen` to F for option 4 only.

All existing types are supported: 0 changes both components to `(swi,swf)`,
1 changes angle to `swf` degrees at fixed magnitude, 2 changes magnitude to
`swf` at fixed angle, 3 changes `kperp`, and 4 changes `kpar`. The existing
angle convention is `atan(kperp/kpar)`. `swi` is used for type 0, not as an
independent starting point for the other types. Logarithms require positive
arguments. Every output wavevector must be finite, with `kperp>=0,kpar/=0`.
The option-4 reader corrects the historical type-2 endpoint bookkeeping when
preparing subsequent segments; legacy option-1/2 bookkeeping is unchanged.

Equilibrium initialization is performed once. At a changed perpendicular
wavevector every rank calls the existing harmonic-limit/assignment routines,
and workers rebuild their Bessel arrays. Pure parallel changes skip this work.
Every point starts with an empty determinant cache; growth iterations and A1
passes at that point share exact raw values. All actual `disp` calls remain
collective. The A1 and growth algorithms and their tolerances are unchanged.

Option 4 automatically creates:

- `solution/<runname>.nyquist.scan`: detailed ASCII diagnostics. Columns are
  `point segment kperp kpar initial_winding A1_status growth_status gamma_lower
  gamma_upper gamma_midpoint gamma_absolute_width gamma_relative_width
  logical_requests disp_calls cache_hits wall_seconds update_seconds
  determinant_seconds passes refinements max_depth growth_iterations
  rank_calls_min rank_calls_max point_status evaluator_status`. Growth status
  `-1` means disabled; other codes are listed above. Point status 0 means the
  requested calculation succeeded; a nonzero value marks a numerical/output
  failure. A valid initial winding can coexist with a failed growth search.
  Real unavailable values are `NaN`; failed searches retain valid bracket
  endpoints but have no midpoint. Header comments give equilibrium, fixed
  frequency bounds, normalization, sampling controls and units.
- `solution/<runname>.nyquist.dat`: exactly four numerical columns,
  **`kperp kpar Wn gamma_max`**, and short `#` comments. `Wn` is the signed,
  accepted winding of the original full frequency rectangle, never a strip
  count. `gamma_max` is the arithmetic midpoint of a converged maximum-growth
  bracket, in reference gyrofrequency units, not a solved root or the input
  upper frequency boundary. Valid zero counts have `Wn=0,gamma_max=NaN`;
  failed initial contours have both `NaN`; valid positive counts with disabled,
  failed or incomplete growth searches retain `Wn` and use `gamma_max=NaN`.
- `solution/<runname>.nyquist.summary`: aggregate requested/evaluated points,
  positive/zero/negative/unresolved counts, successful/failed growth
  calculations, actual determinant calls, rebuild count, shared initialization,
  scan-preparation and total scan times. `scan_status=0` means all requested
  calculations completed; nonzero marks unresolved points or output failure.

Both tables have exactly one aligned row per requested point, including
recoverable failures. Such numerical failures do not stop later points; the
completed scan returns nonzero if any requested calculation remains unresolved.
Invalid input and global initialization failures abort; output errors stop the
scan collectively. A valid zero winding excludes roots only in the prescribed
frequency rectangle. No result certifies roots outside that fixed domain.

Per-point summaries and growth-history files are suppressed. Optional
`nyq_save_samples=T` writes only each point's initial contour to
`<runname>.nyquist.point_NNNNNN.samples`. Option-3 filenames remain unchanged.
All output directories must exist. New validation:
`bash validation/run_nyquist_scan.sh`; results: `NYQUIST_B1_RESULTS.md`.

### Nyquist wavevector planes (`scan_option=5`)

Use `n_scan=2`, the existing `&scan_input_1` and `&scan_input_2` groups, and the
same `&nyquist` controls as options 3/4. Nonrelativistic warm analytic and/or
numerical species and rectangle frequency geometry are required. No new plasma backend or MPI
decomposition is introduced. For example, with the complete mirror species
parameters in `tests/data/nyquist_plane/mirror_plane.in`:

```fortran
&system
 ! Supply the usual plasma/allocation controls too.
 scan_option=5, n_scan=2, kperp=0.2, kpar=0.04
/
&nyquist
 nyq_omega_min=-0.05, nyq_omega_max=0.05, nyq_gamma_max=0.08
 nyq_find_max_growth=T, nyq_growth_spacing='log', nyq_growth_tolerance=1.d-4
/
&scan_input_1
 scan_type=3, swf=0.4, swlog=F, ns=2, nres=1
/
&scan_input_2
 scan_type=4, swf=0.06, swlog=F, ns=2, nres=1
/
```

Scan 1 is the outer loop and scan 2 varies fastest. The output grid has
`(ns_1+1)*(ns_2+1)` points, including initial and endpoint coordinates in each
direction. Here the nine rows are `(0.2,0.04),(0.2,0.05),(0.2,0.06),...,(0.4,0.06)`.
Each outer point starts its own inner scan; starts are retained rather than
deduplicated as consecutive line-segment boundaries. Only output indices
`0,nres,2*nres,...,ns*nres` are evaluated; intermediate root-tracking steps are
skipped. `heating,eigen` remain unused. The two types must differ, type 0 is
disallowed, and all coordinates must satisfy B1's finite/nonzero-parallel rules.

The existing reader calculates each scan's increments **once**, preparing scan 2
after the outer endpoint. `om_double_scan` then applies those same increments
from each outer wavevector. B2 retains this convention: for mixed non-Cartesian
types, an inner scan need not end at the same absolute target on every row.
For example, a kperp outer scan followed by an angle scan translates the fixed
angle increment to each row's initial angle. Do not assume arbitrary scan types
commute when their order is reversed. The existing angle convention remains
`atan(kperp/kpar)`. B1's magnitude-endpoint bookkeeping correction applies to
option 5 as well; legacy options 1 and 2 retain their historical reader behavior.
Wavevector `swlog` is independent of growth spacing. Linear growth uses absolute
tolerance and log growth uses relative width, exactly as at one wavevector.

Options 4/5 share wavevector updates, the authoritative `run_nyquist` calculation,
per-point cache reset, result handling, output and finalization. Rank 0 broadcasts
each wavevector; all ranks refresh the existing harmonic assignments/Bessel arrays
when kperp changes. Pure kpar changes skip rebuilding. Every physical determinant
evaluation remains collective over the original resonance assignments; no rows,
columns or wavevectors are assigned to separate worker groups. Equilibrium
initialization and MPI finalization occur once per plane. Frequency bounds remain
fixed and the exact-frequency cache never spans different grid points.

The B1 filenames and four-column scientific conventions are retained:

- `.nyquist.scan` keeps the first 26 B1 diagnostic columns. For a plane, column 2
  (`segment`) denotes the outer row. Columns **27 and 28** append `outer_index`
  and `inner_index`, both **one-based**, in outer-major order.
- `.nyquist.dat` still has **exactly four numerical columns**,
  `kperp kpar Wn gamma_max`, with no blank row separators. A short comment gives
  grid dimensions and ordering. Wn is the original full-box winding; gamma_max
  is only a converged bracket midpoint, not the input top boundary. The B1 NaN
  conventions apply unchanged to zero, unresolved, disabled and failed searches.
- `.nyquist.summary` adds `grid_dimensions`, `scan_types` and `grid_order` to
  the existing requested/evaluated counts, positive/zero/unresolved outcomes,
  growth-search outcomes, physical calls, rebuilds and timings.

Recoverable numerical failures preserve their rows and valid count/bracket,
continue to later grid points, and yield a nonzero aggregate status. Input,
communication and output failures terminate appropriately. A zero count excludes
roots only in the fixed original frequency rectangle, not outside it. Optional
initial-contour samples retain B1's flat point-index filenames; per-point full
summaries and growth histories are not generated automatically.

Runnable mirror, reverse-order and isotropic inputs are in
`tests/data/nyquist_plane/`. Run `bash validation/run_nyquist_plane.sh` for the
independent pre-B2 coordinate gate, per-point/row/column comparisons and failure
checks. Results, performance and limitations are in `NYQUIST_B2_RESULTS.md`.
