# ALPS Input

This is a reference for the key input parameters used by ALPS.

## Namelists in execution input files.

The following namelists and associated input parameters are read in by ALPS from the input file.

### *&system*  
General system parameters.

**`profile_startup`**
Optional logical in `&system`, default `.false.`. 
When enabled, rank 0 reports elapsed wall time in seconds using `MPI_Wtime`.
This includes 

- VDF input, 
- original-grid processing, 
- relativistic preparation (if applicable), 
- spline matrix/solve/evaluation, 
- postprocessing,
- diagnostic output, 
- analytic-continuation fitting, and total 
- VDF startup before the first dispersion calculation. 

Species timings are summed. 
Detailed reporting is independent of `writeOut`.

**`kperp`**  
Initial perpendicular wavevector $k_{\perp} d_{ref}$.

**`kpar`**  
Initial parallel wavevector $k_{\parallel} d_{ref}$.

**`nspec`**   
Number of plasma species.

**`nroots`**   
Number of dispersion solutions to find and follow.

**`use_map`**   
Choice of:  

- True: Searching for roots over a map in complex frequency space (see &maps_1 namelist).  
- False: Input `nroots` guesses for solutions (see &guess_1 namelist).

**`writeOut`**  
Write or suppress output to screen.

**`nperp`**  
Perpendicular momentum space resolution, $N_{\perp}$.
The input file must have $N_{\perp}+1$ values spanning parallel momentum space.

**`npar`**  
Parallel momentum space resolution, $N_{\parallel}$.
The input file must have $N_{\parallel}+1$ values spanning parallel momentum space.

**`vA`**  
Alfven velocity using the mass density of the reference species, 
normalized to the speed of light, $v_{A,ref}/c$.

**`arrayName`**  
Name of input array, located in 'distribution' folder.

**`Bessel_zero`**  
Maximum amplitude of Bessel function to determine `nmax`.

**`secant_method`**  
Selection for root finding method.

- 0: *secant method from NHDS.* Uses full secant steps and declares convergence when the **magnitude of the dispersion determinant** falls below `D_threshold`.
- 1: *rtsec method from PLUME.* Applies half of each calculated secant step and stops when the calculated **step magnitude** falls below `D_tol`, or the determinant is exactly zero.
- 2: *Secant method with oscillation control (default):* Adds adaptive damping, caps steps at 10% of the current frequency magnitude, and switches to a regularized finite-difference Newton step after repeated detections of revisited frequencies; convergence uses `D_threshold`.

All three use `D_prec` to construct a second initial frequency and `numiter` to limit iterations.

**`numiter`**  
Maximum number of iterations in secant method.

**`D_threshold`**  
Minimum threshold for secant methods 0 and 2.

**`D_prec`**  
Size of bounding region for secant method.

**`D_gap`**  
Size of allowable difference between roots.

**`D_tol`**  
Tolerance for secant method = 1, rtsec.

**`positions_principal`**  
Number of parallel momentum steps distant from the resonant momentum
included in the numerical calculation of Verscharen et al 2018, Eqn 3.5, $M_{I}$.

**`n_resonance_interval`**  
How many steps should be used to integrate around the resonance,
$M_{P}$, used for integrating near poles (see section 3.1 of Verscharen et al 2018).

**`Tlim`**  
Threshold for analytical principal-value integration, $t_{\mathrm{lim}}$.

**`determine_minima`**  
If true, after map search, determine minima and refine solutions.

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

### *&nyquist*

Required for `scan_option=3`, `4` or `5`. Nyquist supports warm analytic species:
all `use_bM=T`, `relat=F`, positive beta/anisotropy/density/mass, nonzero charge,
finite parameters, `vA>0`, and `kpar /= 0`. Numerical-VDF and cold Nyquist
interpretations await separate validation; their existing root/scan modes are unaffected.
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

| Parameter | Meaning / default |
|---|---|
| `nyq_geometry` | `'rectangle'` (default) or `'horizontal'`, lowercase. |
| `nyq_omega_min`, `nyq_omega_max` | Required finite ordered real-frequency bounds, in reference gyrofrequency units. |
| `nyq_gamma_min` | Finite strictly positive imaginary frequency; default `1.d-4`. |
| `nyq_gamma_max` | Required only for rectangles; finite and greater than `nyq_gamma_min`. |
| `nyq_sampling` | The A1 `nyq_options` object; set individual components with `%` syntax below. |
| `nyq_save_samples` | Default `F`; opt into ordered sample storage and ASCII output. |
| `nyq_normalize_rectangle` | Default `F`; diagnostic normalized rectangle for comparison with horizontal phases. Applies only to rectangles. |
| `nyq_find_max_growth` | Default `F`; enable bounded maximum-growth estimation with cached closed strips. Rectangles only. |
| `nyq_growth_tolerance` | Default `1.d-4`; finite positive absolute bracket-width tolerance for `'linear'`, dimensionless relative bracket-width tolerance for `'log'`. Both use strict inequality. |
| `nyq_growth_max_iterations` | Default 32; positive integer limit on bisection steps. |
| `nyq_growth_spacing` | Default `'linear'`; `'log'` enables geometric trials and nearest-original-boundary strips. Used/validated only when maximum-growth estimation is enabled. |

The `nyq_sampling` components use A1 defaults directly; there is no second set
of sampling controls or implementation:

| Component | Default |
|---|---:|
| `initial_intervals` (per segment) | 8 |
| `max_depth` (per initial interval) | 24 |
| `max_evaluations` (logical requests over all passes of each contour) | 200000 |
| `max_passes` | 6 (must be at least 3) |
| `phase_limit` (radians) | 0.7 |
| `curvature_limit` | 0.15 |
| `phase_tolerance` (absolute radians across passes) | `1e-8` |
| `integer_tolerance` (closed winding units) | `1e-8` |
| `log_magnitude_limit` | 1.5 |
| `local_ratio` | 0.05 |

For example, `nyq_sampling%max_depth=32` permits additional local refinement.
The kernel also requires global sampling agreement, including an independently
seeded pass. A small integer residual alone is insufficient. Finite sampling is
numerical evidence, not proof that arbitrary phase loops cannot be missed.

Rectangles use the raw determinant by default and report a signed accepted
winding only after convergence. The warm response has no upper-half-plane poles,
so this counts roots in the stated finite rectangle. It does not count all
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
same `&nyquist` controls as options 3/4. Warm nonrelativistic analytic species and
rectangle frequency geometry remain required. No new plasma backend or MPI
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


### *&fitting_params*

Optional global Levenberg–Marquardt (LM) settings shared by all numerical species using `AC_method=1`, 
for both nonrelativistic and relativistic cases. 
Omission of the full list or individual keys retains the defaults below.

```fortran
&fitting_params
    maxsteps_fit       = 500
    lambda_initial_fit = 1.d0
    lambdafac_fit      = 1.d1
    epsilon_fit        = 1.d-8
/
```

**`maxsteps_fit`** (integer): 
maximum LM fitting iterations.

**`lambda_initial_fit`** (double precision): 
initial LM damping parameter.

**`lambdafac_fit`** (double precision): 
damping multiplier when a trial increases the sum of squared residuals, and divisor otherwise.

**`epsilon_fit`** (double precision): 
convergence threshold for the absolute change in the sum of squared residuals, not the residual itself.

LM parameter files are written during fitting; fit-result diagnostics are
controlled by `writeOut`.

Methods 0 and 2 do not require this name list. 
Method 2 uses `&poly_spec_j`; mixed method-1
and method-2 species use their respective configurations independently.

### *&relativistic_params*

Settings shared by all relativistic species. 
ALPS reads this group when any `&spec_*` block sets `relat=T`, 
after reading all species and before sizing fit arrays. 

```fortran
&relativistic_params
  ngamma = 500
  npparbar = 500
  rel_cache = 'auto'
  rel_cache_dir = 'distribution/'
/
```

The cache is a reusable snapshot of the completed relativistic interpolation. 
ALPS still reads the original ASCII VDF, but with rel_cache='auto' it first checks whether it has already calculated the transformed grids, normalized distribution, and derivatives for those inputs. 
These arrays are saved in a binary file, which preserves every floating-point bit without rounding values into ASCII text. 
The long filename is an automatically generated fingerprint of the VDF contents, interpolation settings, and computational build. 
Matching inputs can reuse the same file across runs or species; changing relevant inputs produces a different fingerprint. 
Before loading, ALPS checks compatibility and a checksum to detect damage. 
If the cache is missing or invalid, ALPS simply calculates the interpolation again and saves a fresh copy. 
The original VDF files remain the authoritative inputs, and fitting and dispersion calculations continue normally.

An absent or empty group uses `ngamma=100`, `npparbar=200`, `rel_cache='off'`, and
`rel_cache_dir='distribution/'`. 
Partial groups retain defaults for omitted fields.
When no species is relativistic, the entire group is ignored, including validation.

**`ngamma`**
Relativistic momentum space resolution, $N_{\Gamma}$.
Requires `ngamma >= 2` when a relativistic species is present.

**`npparbar`**
Relativistic parallel momentum space resolution, $N_{\bar{p}_{\parallel}}$.
Requires `npparbar >= 4` to keep the existing three-point resonance stencil in the interior.

**`rel_cache`**
Optional `&relativistic_params` character parameter, default `'off'`. Modes are case-insensitive:

- `'off'`: construct the relativistic representation normally, without cache identity or file access.
- `'auto'`: load a fully validated entry, or construct and save one on a miss.
- `'refresh'`: construct every relativistic species and atomically replace its entry.

**`rel_cache_dir`**
Optional `&relativistic_params` character parameter, default `'distribution/'`.
The directory must already exist.
Only rank 0 accesses entries.
Missing/unwritable storage, incompatible entries, and integrity failures fall back to normal calculation. Invalid modes are input errors.
Existing relativistic diagnostic files are still written, including on cache hits and when `writeOut=.false.`.
Analytic-continuation fitting is always performed normally.

Entries are identified by numerical VDF contents, coordinates, construction settings, and strict executable/library compatibility, independently of filenames and species
indices.

### *&guess_m*  
Initial guess of complex frequency for $m$th solution.  
Only used when `use_map`=.false.  
Need to have number of name lists equal to `nroots`.

**`g_om`**  
Guess for real solution $\omega_{r}/\Omega_{ref} $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_{ref} $.


### *&maps_1*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{r}/\Omega_{ref}$ map search.
Spacing automatically calculated between `omi` and `omf`.  

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{ref}$ map search.
Spacing automatically calculated between `gami` and `gamf`  

**`omi`**  
Smallest $\omega_{r}/\Omega_{ref}$ value for complex map search.

**`omf`**  
Largest $\omega_{r}/\Omega_{ref}$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_{ref}$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_{ref}$ value for complex map search.

**`ni`**  
Number of $\gamma/\Omega_{p}$ points in frequency grid.

**`nr`**  
Number of $\omega_{r}/\Omega_{ref}$ points in frequency grid.


### *&spec_j*  
Species parameters list for distribution $f_{j}$.
The first species is set as the reference.

**`nn`**  
Relative density $n_{j}/n_{ref}$.

**`qq`**  
Relative charge $q_{j}/q_{ref}$.

**`mm`**  
Relative mass $m_{j}/m_{ref}$.

**`ff`**  
Number of fitted functions for analytical continuation calculation, 
when using `AC_method`=1.

**`relat`**  
Treat $f_{j}$ as non-relativistic or relativistic. 
Analytic bimaxwellian species (`use_bM=T`) require `relat=F`.

**`log_fit`**  
Use linear or $\log_{10}$ fitting routine. 
Relativistic numerical Chebyshev
continuation (`relat=T`, `AC_method=2`) supports only `log_fit=T`.

**`use_bM`**  
Use actual numerical integration (F) or 
bi-Maxwellian/cold-plasma proxy via NHDS routines (T),
with parameters read in from &bM_spec_j namelist.

**`AC_method`**  
Choose the method for the evaluation of the analytic continuation:

- 0: Use the function that is defined analytically in distribution/distribution_analyt.f90
- 1: Use the fit routine as defined in the &ffit_j_k namelist.
- 2: Use a polynomial basis representation as defined in the &poly_spec_j namelist. 
*The limits of this analytic continuation for 
$|\gamma|\gtsrim |\omega_{r}|$ are being investigated*.


Relativistic numerical method 2 uses native SVD/Clenshaw models at fixed Lorentz factor, with `kind=1`, positive requested `order`, and `log_fit=T`. 
Original source samples must be finite and strictly positive on a finite, increasing
rectangular momentum grid.
Zeros are rejected before the log-TPS transformation,
including when caching is enabled. 
No positive floor is substituted.

These models represent the **completed transformed distribution** on the existing Gamma/parallel-momentum grid. Its global natural-log thin-plate spline extends
beyond the supplied source momentum rectangle at some cone-admissible nodes.
Those values are TPS extensions, not source observations. ALPS reports source and
completed bounds, interior/exterior counts, affected Gamma rows and fresh/cache
origin once per selected species. Counts are recomputed from the actual geometry;
cache hits verify identity against the current source arrays and coordinates.
Fit status, degree, rank and residual do not certify source support or accurate
complex continuation. Numerical range or unavailable-row failures remain explicit.
Original-grid fit files and fitted-moment totals omit these fixed-Gamma models;
original VDF moments are unchanged.


### *&ffit_j_k*
Initial Fit Values for species $j$, function $k$.

**`fit_type_in`**  
Kind of fit function:

- 1: Maxwellian,  

$$F_M(\hat{p}\_{\parallel})=u_{1}\mathrm{exp}[-y{\hat{p}}^2\_{\perp}-u_{2}(\hat{p}\_{\parallel}-u_{3})^2]$$

- 2: Kappa,  

$$F_{\kappa}(\hat{p}\_{\parallel})=u_{1}[1+u_2({\hat{p}}\_{\parallel}-u_{3})^2+yu_{5} {\hat{p}}^2\_{\perp}]^{u_{4}}.$$

- 3: Juettner with $p_{\perp},p_{\parallel}$,  

$$F_{J}(\hat{p}\_{\perp},\hat{p}\_{\parallel})=
u_{1}\mathrm{exp}\left[-u_{2}\sqrt{1+\frac{\hat{p}^2\_{\perp}+(\hat{p}^2\_{\parallel}-u_3)^2 v_A^2}{m_{j}^2 c^2}}\right].$$

- 4: Juettner with variable $\Gamma$, constant $\bar{p}_{\parallel}$,  

$$F_{J}(\Gamma)= u_{1} \mathrm{exp}[-y \Gamma].$$

- 5: Juettner with $p_{\perp},p_{\parallel}$; variable $\bar{p}_{\parallel}$,  

$$F_{\kappa}(\hat{p}\_{\perp},\hat{p}\_{\parallel})=u_{1}\mathrm{exp}[-y \hat{p}\_{\perp}]\mathrm{exp}[-u_{2}*(\hat{p}\_{\parallel}-u_{3})^2].$$

- 6: Bi-Moyal distribution

$$F_{bMo}(\hat{p}\_{\perp},\hat{p}\_{\parallel})= u_{1} \mathrm{exp}[0.5 (y u_4 \hat{p}^2\_{\perp} + u_{2} (\hat{p}\_{\parallel} -u_{3})^2 -\mathrm{exp}(y u_{4} \hat{p}^2\_{\perp} + u_{2} (\hat{p}\_{\parallel}-u_{3})^2) )].$$

**`fit_1`-`fit_5`**  
Fit parameters, $u_{1}$ - $u_{5}$, defined in the above equations for each of the types of fit functions.
Not all parameters will be used for all functions.  
Suggested values for parameters generated by generate_distribution.

**`perpcorr`**  
This parameter, $y$ in Eqn. B1, compensates for the strong
$p_{\perp}$ dependence of $u_1$, making the fit more reliable.


### *&bM_spec_j*
Bi-Maxwellian/cold-plasma parameters; for species j.
Only used if `use_bM=T`.

**`bM_nmaxs`**  
Maximum number of resonances to consider.

**`bM_Bessel`**  
Precision threshold for $I_n$.

**`bM_betas`**  
$\beta_{\parallel,j}$ of bi-Maxwellian distribution $f_{j}$. If this variable is set to 0.d0, then the code will treat the given species with the susceptibility from cold-plasma theory.

**`bM_alphas`**  
$T_{\perp,j}/T_{\parallel,j}$ of bi-Maxwellian distribution $f_{j}$.

**`bM_pdrifts`**  
Relative drift of bi-Maxwellian distribution $f_{j}$ or the cold plasma species in units of $m_{ref} v_{A,ref}$.


### *&poly_spec_j*
Input for the polynomial representation of the input distribution for the analytical continuation.
Only used if `AC_method=2`.

**`kind`**  
Type of the basis polynomial:

- 1: Chebychev

**`order`**  
Maximum order of the basis polynomial.

Native LOG10 continuation reconstructs the physical distribution using the shared
Chebyshev kernel. Overflow, underflow and invalid evaluations have explicit
numerical statuses; there is no adjustable logarithmic suppression threshold.
The former `log_max` namelist setting is unsupported and must be removed from inputs.

**`derivative_type`**  
Controls the **nonrelativistic** parallel derivative in the Landau integral
when using polynomial continuation. 
Relativistic AC method 2 always uses the shared
analytic derivative with respect to `u=P_parallel/(m_s*c)` at fixed Gamma; its
Gamma finite-difference stencil remains unchanged.

- 0: Finite Difference Derivative
- 1: Analytic Derivative of Polynomial Basis Function (default)

Omitting `derivative_type` in `&poly_spec_j` selects `1` for that species.
Set `derivative_type=0` explicitly to use the finite-difference derivative.

### *&scan_input_l*
Inputs for scanning parameter space for $l$th scan.  

**`scan_type`**  
Type of parameter scan:

- 0: Current value of $\textbf{k}$ to $k\_{\perp}$=`swi` and $k\_{\parallel}$ =`swf`.   
- 1: $\theta_0 \rightarrow \theta_1$ at fixed $|k|$ from current value of $\theta=\mathrm{atan}(k\_{\perp}/k\_{\parallel})$ to `swf`.  
- 2: Wavevector scan at fixed angle $\theta_{k,B}$ to $|k|$ =`swf`.  
- 3: $k\_{\perp}$ scan with constant $k\_{\parallel}$ to $k\_{\perp}$=`swf`.  
- 4: $k\_{\parallel}$ scan with constant $k\_{\perp}$ to $k\_{\parallel}$=`swf`.  

**`swi`**  
Target perpendicular component for `scan_type=0`; unused for other scan types.

**`swf`**  
Scan variable to define end of scan through wavevector space.

**`swlog`**  
Use $\log_{10}$ (T) or linear (F) spacing.

**`ns`**  
Number of output steps after the initial wavevector. An ordinary segment file
has `ns+1` rows; option 4 writes shared segment boundaries only once.

**`nres`**  
Root-tracking steps per output interval. Option 4 evaluates only output points,
using the same coordinate indices without the intermediate tracking work.

**`heating`**  
Calculates heating rates if true.

**`eigen`**  
Calculates eigenfunctions if true.     

## Namelists in distribution input files.

The following namelists and associated input parameters are read in by the routine `generate_distribution` from input files. This generates `nspec` arrays of distributions that are read in by ALPS

### *&system*
General system parameters.

**`nspec`**

Number of plasma species.

**`beta`**

Reference plasma beta.

**`vA`**

Reference Alfven velocity, normalized to the speed of light, $v_{A,ref}/c$.

**`nperp`**  

Perpendicular momentum space resolution, $N_{\perp}$.
The input file will have $N_{\perp}+1$ values spanning parallel momentum space.
Set the same values of `nperp` in both the distribution input file and ALPS execution input file.

**`npar`**  

Parallel momentum space resolution, $N_{\parallel}$.
The input file will have $N_{\parallel}+1$ values spanning parallel momentum space.
Set the same values of `npar` in both the distribution input file and ALPS execution input file.

**`maxP`**

Maximum value of reference Alfven momentum, $m_{ref} v_{A,ref}$.

**`writeName`**

Name of output arrays.

### *&spec_j*

**`ms_read`**

Mass ratio $m_j/m_{ref}$.

**`taus`**

Temperature ratio $T_{\parallel,j}/T_{\parallel,ref}$.

**`alphs`**

Temperature anisotropy $T_{\perp,j}/T_{\parallel,j}$.

**`ps`**

Normalized drift momentum $m_j v_{drift}/m_{ref} v_{A,ref}$.

**`kappas`**

$\kappa$ index.
Only used if `distribution`=2.

**`distributions`**

Type of distribution:

- 0: Use distribution from distribution_analyt.
- 1: Bi-Maxwellian distribution.
- 2: Bi-Kappa distribution.
- 3: Anisotropic Juettner distribution.
- 4: Bi-Moyal distribution.

**`autoscaleS`**

If `autoscaleS`=.true., define maximum momentum for each array automatically using global value `maxP`.

**`maxPperpS`**

If `autoscaleS`=.false., force maximum perpendicular momentum for species array.

**`maxPparS`**

If `autoscaleS`=.false., force maximum parallel momentum for species array.




