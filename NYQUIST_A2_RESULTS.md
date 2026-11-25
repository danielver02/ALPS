# ALPS Nyquist Stage A2 — integration and bi-Maxwellian validation

Date: 2026-11-20 (America/Phoenix).

**A2a integration and A2b bounded-domain validation passed.** Option 3 now uses
the existing collective determinant at one prescribed wavevector. Rectangles
resolve the growing pair, their individual inclusion/exclusion tests, the stable
control, the growth threshold, and selected negative-wavevector tests. The
normalized horizontal calculation is sampling-converged over the tested intervals
through R=1e5, but its closure is substantial. Larger intervals remain unresolved
near the plasma-frequency scale. **No complete horizontal instability count has
been established.**

The A1 mathematical kernel is unchanged. No susceptibility, `disp`, root-search,
analytic-continuation, numerical-VDF, or harmonic-assignment routine was modified.
Pre-existing unrelated working-tree changes and scientific outputs were preserved.
All builds, MPI calculations, and generated scientific outputs used fresh external
scratch directories. No Python or new numerical dependency was introduced.

## 1. Implementation and build

| File / routine | Change |
|---|---|
| `src/ALPS.f90`, program dispatch/finalization | Broadcast operation selection after input; reject reserved options collectively; read Nyquist input on rank 0; enter the driver after shared initialization and before root search; return to common finalization. |
| `src/ALPS.f90`, startup Bessel guard | Use `.not.proc0` instead of the undefined rank-0 `sproc` sentinel. Worker harmonic assignment is unchanged. |
| `src/ALPS_io.f90`, `init_param` | For option 3 set `nroots=0,use_map=F,n_scan=0` and skip guess/map instructions. For reserved 4/5 return after system input so unused species/scan inputs are not required. Legacy input paths are otherwise unchanged. |
| `src/ALPS_nyquist_driver.f90`, new module | `read_nyquist`, `run_nyquist`, synchronous evaluator, log-polar normalization, compact summary, optional ASCII samples, collective shutdown/error handling. |
| `src/Makefile.am`, `src/Makefile.in` | Register A1 and adapter sources and explicit Fortran module dependencies for clean parallel builds. Only the src Automake template was regenerated. |
| `input.md` | Operation selection, namelist, controls/defaults, output, statuses and interpretation. |
| `validation/build_nyquist.sh` | Isolated identified ALPS build with existing MPI/BLAS/LAPACK dependencies. |
| `validation/make_nyquist_cases.f90` | Fortran generator for integration, legacy, physical, expansion and failure cases. |
| `validation/check_nyquist_summary.f90` | Fortran checks of ordinary roots, summaries, rank counts, finite-phase regressions and independent closure sums. |
| `validation/run_nyquist.sh` | A2a gate, physical tests, optional pre-change comparison, copied-adapter fault injection, output-error checks. |
| `tests/data/nyquist/` | Unstable/stable rectangles, normalized horizontal and ordinary-root inputs, measured finite-phase regression values, fixture notes. |
| `validation/run_numerical_baselines.sh`, `run_rel_production.py`, `run_b01_diagnostics.py` | Add the two new module dependencies to existing scratch build lists. Existing test logic and physical tolerances are unchanged. The existing Python harnesses are not dependencies of the new Fortran/shell suite. |
| `validation/README.md` | New standalone build/test entry points. |

The adapter keeps MPI and ALPS state outside `ALPS_nyquist`. No A1 changes were
needed: its source SHA-256 before/after A2 is
`d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.

Compiler: GNU Fortran 11.4.0. MPI: Open MPI 4.1.2, four ranks, local
`--oversubscribe`. Existing `-llapack -lblas` libraries are used; thread counts
were set to one in the regression suites. The final optimized suite used
`-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; the debug suite used
`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`.

Final validated executable SHA-256 values:

- Optimized: `6c927164ed8355069c0bf2e10127a3b42e6c8a174a53ab5a50cdfb1902302deb`.
- Debug: `eb204d7093aec7d8ae5b191e2446454e08ca9e4232af556ce95413ae9500928c`.
- Pre-change compatibility build: `f8b2957947e79a1a864831a10797f4b10c5d7a9d7b6db7dc76fae1529c9e6841`.

The pre-change build used a source snapshot taken before A2 edits, identical
compiler flags to the optimized build, and the current baseline
`5ee0d335a037744a5ce298c355b84fc316af5f9c` plus pre-existing workspace content.
Source and executable manifests are retained in each scratch build's
`provenance.txt`. Final complete-suite artifacts for this session are
`/tmp/alps-nyquist-a2.aJ53gNW6` (optimized) and
`/tmp/alps-nyquist-a2.DZhnR0rv` (debug). Temporary directories are not durable
archives; the fixtures and measured results here are the retained deliverables.

A separate fresh Automake source/build tree passed configure and `make -j4 -C src`
with the existing production `-O4 ...` flags; that executable also passed the
unstable-rectangle smoke test. No in-tree executable or object was overwritten.

Reproduce the complete isolated test:

```sh
bash validation/run_nyquist.sh
FFLAGS='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow' \
  bash validation/run_nyquist.sh
```

For an optional byte comparison, set `ALPS_REFERENCE_EXE` to an identified
pre-change executable. For a smaller isolated build:

```sh
bash validation/build_nyquist.sh /tmp/new-alps-build
mkdir -p /tmp/new-alps-run/solution /tmp/new-alps-run/distribution
cd /tmp/new-alps-run
mpirun --oversubscribe -np 4 /tmp/new-alps-build/ALPS \
  /absolute/path/to/ALPS_git/tests/data/nyquist/unstable_rectangle.in
```

The normal source build also includes the modules through `src/Makefile.am` and
its regenerated template. The new validation runner requires no plotting tools.

## 2. Input, normalization and output contract

`scan_option=3` always uses the input `kperp,kpar`; no root guesses, maps, or scan
instructions are needed. The old selectors `nroots,use_map,n_scan` are ignored
for this operation. Integration tests deliberately set `use_map=T,n_scan=2` while
omitting every map, guess and scan group, demonstrating the explicit bypass.
Archived examples use simpler `use_map=F,n_scan=0` settings. Existing options 1/2
retain their root/map operation followed by the requested scans. Options 4/5
print a reserved-operation error, finalize all ranks, and return nonzero without
initializing an incomplete operation.

Example:

```fortran
&system
 ! Existing equilibrium, species, wavevector and allocation controls also supplied.
 scan_option=3
/
&nyquist
 nyq_geometry='rectangle'
 nyq_omega_min=-1, nyq_omega_max=1
 nyq_gamma_min=0.01, nyq_gamma_max=0.2
 nyq_save_samples=F
/
```

For a horizontal path set `nyq_geometry='horizontal'` and omit `nyq_gamma_max`.
Frequency bounds must be finite/ordered, the lower imaginary boundary strictly
positive, and a rectangle's upper boundary larger. No domain extension or
origin evaluation is implicit. A2 restricts this operation to finite physical
parameters, positive `vA`, nonzero `kpar`, and warm nonrelativistic analytic
bi-Maxwellian species. Cold/numerical-VDF Nyquist use awaits separate validation;
legacy calculations with those backends retain their previous input behavior.

`nyq_sampling` is the **actual A1 `nyq_options` instance**, read directly as a
namelist derived type. Set, for example,
`nyq_sampling%initial_intervals=16` or `nyq_sampling%max_evaluations=10000`.
Defaults are inherited without a duplicate control structure: 8 intervals per
segment, depth 24, 200000 total calls, 6 passes, phase limit 0.7, curvature 0.15,
log-magnitude limit 1.5, local ratio 0.05, phase/integer tolerances 1e-8.
The validation inputs explicitly use depth 32. Strict tests use 16 intervals,
phase limit 0.35, curvature 0.075, and phase/integer tolerances 1e-10. All other
criteria remain A1's, including the independently seeded convergence pass.

Rectangles use raw `D_ALPS` by default. `nyq_normalize_rectangle=T` enables a
validation-only normalized rectangle and retains its edge diagnostics without
requiring an ASCII file. Horizontal paths always use `D_ALPS/(z*vA)**6`.
The adapter computes that division in log-polar form using component scales,
retaining normalized magnitude as well as continuously correct complex phase.
It never forms an extreme sixth-power divisor. Nonfinite/zero raw values are
passed to A1 for rejection; normalization never hides them. An unrepresentable
normalized amplitude returns evaluator code 101 instead of magnitude clipping.
This cannot repair loss of accuracy or overflow inside `disp` itself.

Default summary: console and `solution/<runname>.nyquist.summary`. Optional
`nyq_save_samples=T`: `solution/<runname>.nyquist.samples`, with ordered columns
`Re(z), Im(z), Re(F), Im(F), cumulative_phase`. Existing root/map/scan filenames
are not used. Rerunning a name replaces only that name's Nyquist files.

Summaries record wavevector, equilibrium parameters, geometry/bounds,
normalization, signed phase/2pi, accepted integer or `NONE`, A1 status and a text
explanation, evaluator code, per-rank call extrema, MPI timing, refinements,
passes/depth, maximum increment, endpoints, and the effective namelist.
Open output always says `asymptotic_closure=UNVERIFIED`. Failed sampling explicitly
labels phase/endpoints as an unaccepted latest path or partial prefix, not a
completed integral. Root counts are reported only for converged finite rectangles.

## 3. Collective execution and A2a checkpoint

Shared ALPS initialization runs once: input, equilibrium/derivatives, continuation
setup, distribution/model transfer, harmonic limits/assignment, arrays and Bessel
setup. The A2 branch follows the startup report and precedes `map_search` and
`refine_guess`. No initialization is repeated for each frequency.

Rank 0 calls A1 and supplies the synchronous callback. Each callback broadcasts
an evaluate command and frequency, then all ranks call the unmodified `disp`.
Its two reductions, determinant broadcast and barrier are unchanged. Rank 0
checks/normalizes the returned value and resumes A1. Workers wait for the next
command. A1 returns normally even on a recoverable failure; rank 0 always sends
the stop command. All ranks then compare their actual local call counts, reduce
wall time, receive the final operation status, and use common finalization.
There is no new frequency/wavevector parallelization, and ranks with zero local
susceptibility contribution still make every collective call.

A2a passed before physical assertions:

| Check | Result |
|---|---|
| Independent A1 compile/tests | 306 checks, zero failures; debug and optimized printed results agree. |
| New modules compile/link | Manual isolated builds and clean parallel Automake build passed. |
| Four-rank debug smoke | Requested wavevector/rectangle used; all four ranks made 630 calls, three passes; summary and ordered samples written. Its winding was provisional at this checkpoint. |
| No root solver in option 3 | No `.roots`/map/scan result created, no guesses/maps/scans required, common finalization reached. |
| Existing root/map/scans | Nine scientific output files byte-identical to the pre-change optimized build. |
| Existing input suite | All ten `validation/run_input_validation.sh` cases passed. |
| Reserved options 4/5 | Clear error, no species or scan inputs required, all-rank MPI finalization, nonzero launcher status. |

A debug error-path test exposed a finalization race: workers returning nonzero
immediately after `MPI_Finalize` could cause the launcher to kill rank 0 before
its last diagnostics were flushed. The final implementation lets rank 0 report
the nonzero process exit after writing/closing diagnostics; workers finish
normally. Operation status is still broadcast to all ranks before finalization.
Both complete suites passed after this correction.

## 4. Reference equilibrium and ordinary roots

Both species are analytic (`use_bM=T,relat=F`), density 1, zero drift. Proton
charge/mass are 1; electron charge −1 and mass `5.44617021487e-4`.
Both parallel betas are 3, so isotropic electrons have `Te=Tparallel_p`.
Proton anisotropy is 2 (unstable) or 1 (stable control). Wavevector is
`(kperp,kpar)=(0.001,0.3)`, `vA/c=1e-4`; selected repeats use `kpar=-0.3`.
Harmonic cap 500 and Bessel threshold 1e-50 are the A0 controls. The 4×4 momentum
arrays are unused initialization placeholders for this analytic backend.

The identified ordinary solver build reproduced both signs before regression
fixtures were established. Initial guesses were `±0.5+0.08i`; `secant_method=0`,
100 iteration cap, residual threshold 1e-13. Existing root files print limited
frequency precision:

| Equilibrium | Root z | Printed Re(D), Im(D) |
|---|---|---|
| Anisotropic | `+0.39666+0.085986i` | `+1.7577e-14, +4.6759e-15` |
| Anisotropic | `-0.39666+0.085986i` | `+1.9785e-14, -5.4347e-15` |
| Isotropic | `+0.40460-0.0018543i` | `-3.5680e-14, -3.2335e-15` |
| Isotropic | `-0.40460-0.0018543i` | `-6.9250e-15, +5.3652e-14` |

All residual magnitudes are below the requested threshold. Root-regression
frequency tolerances are 1e-5 real and 1e-6 imaginary, respecting the printed
precision. The solver reproduces candidates, not an exhaustive inventory.
The independently converged rectangles below establish counts in their regions.

## 5. Closed-contour physical results

All contours are counterclockwise, entirely above the real axis. The warm analytic
response has no upper-half-plane poles, so accepted windings count enclosed
roots with multiplicity. No absolute value or root-count doubling is used.

| Case | Real bounds | Imaginary bounds | W | Calls | Passes | Refinements | Max depth | Final max increment (rad) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Positive root | `[0.3,0.5]` | `[0.04,0.13]` | 1 | 448 | 3 | 0 | 0 | 0.08050 |
| Negative root | `[-0.5,-0.3]` | `[0.04,0.13]` | 1 | 448 | 3 | 0 | 0 | 0.08067 |
| Positive displaced | `[0.52,0.7]` | `[0.04,0.13]` | 0 | 448 | 3 | 0 | 0 | 0.05866 |
| Negative displaced | `[-0.7,-0.52]` | `[0.04,0.13]` | 0 | 448 | 3 | 0 | 0 | 0.05583 |
| Trial, both roots | `[-1,1]` | `[0.01,0.2]` | 2 | 630 | 3 | 91 | 6 | 0.45822 |
| Same, strict sampling | `[-1,1]` | `[0.01,0.2]` | 2 | 1176 | 3 | 140 | 6 | 0.23770 |
| Expanded, strict | `[-2,2]` | `[0.005,0.4]` | 2 | 1288 | 3 | 196 | 8 | 0.24498 |
| Isotropic control | `[-1,1]` | `[0.01,0.2]` | 0 | 630 | 3 | 91 | 6 | 0.48790 |
| Raised threshold | `[-1,1]` | `[0.1,0.3]` | 0 | 570 | 3 | 61 | 4 | 0.44421 |
| Negative kpar, both | `[-1,1]` | `[0.01,0.2]` | 2 | 630 | 3 | 91 | 6 | 0.45822 |
| Negative kpar, positive root region | `[0.3,0.5]` | `[0.04,0.13]` | 1 | 448 | 3 | 0 | 0 | 0.08050 |
| Trial, normalized determinant | `[-1,1]` | `[0.01,0.2]` | 2 | 788 | 3 | 170 | 7 | 0.46214 |
| Expanded, strict, normalized | `[-2,2]` | `[0.005,0.4]` | 2 | 1506 | 3 | 305 | 9 | 0.24451 |

The individual-root contour's nearest edge is about 0.044 away from its measured
root. The displaced contour starts about 0.123 beyond that root's real frequency.
The trial bottom is about 0.076 below its growth rate; the threshold bottom is
about 0.014 above it. These margins are much larger than the ordinary root file's
printing uncertainty. Counts were measured before being pinned by the tests.
Raw and normalized rectangular windings agree, as expected because z=0 is outside.
Both complete debug/optimized suites reproduce these integers.

These are finite-region counts, including a growth threshold. They exclude no
roots outside the tested rectangles or below their lower boundaries.

## 6. Horizontal expansion and finite closure

Use normalized determinant, `gamma=0.01`, and symmetric real bounds `[-R,R]`.
No result is rounded into an accepted integer. All successful entries passed
three A1 global passes including the independent seed. The measured finite
phases are archived as regression observations, not analytic asymptotic counts.
Debug and optimized values satisfy the 1e-8-radian regression tolerance.

| R | Finite phase (rad) | Phase / 2pi | Calls | MPI seconds | Refinements | Max depth | Sampling |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 16.191182474015 | 2.576906725242 | 399 | 0.022867 | 142 | 7 | Converged |
| 10 | 20.424154357934 | 3.250605124537 | 511 | 0.029997 | 198 | 10 | Converged |
| 100 | 25.131904580301 | 3.999866843269 | 605 | 0.039528 | 245 | 14 | Converged |
| 1000 | 25.132678082080 | 3.999989949900 | 713 | 0.037491 | 299 | 17 | Converged |
| 10000 | 31.415914389745 | 4.999998066880 | 1027 | 0.033792 | 456 | 20 | Converged |
| 100000 | 31.415925266686 | 4.999999797999 | 1177 | 0.047364 | 531 | 24 | Converged |
| 1000000 | No completed phase | No count | 115 | 0.006282 | 65 | 32 | Unresolved, status 5, pass 1 |
| 10000000 | No completed phase | No count | 122 | 0.010547 | 65 | 32 | Unresolved, status 5, pass 1 |

Endpoint values for the six completed paths (rounded here; summaries retain full
precision) are:

| R | Fhat(-R+0.01i) | Fhat(+R+0.01i) |
|---:|---|---|
| 1 | `-1.6634732111e24 - 6.7504684870e24 i` | `-1.6634732111e24 + 6.7504684870e24 i` |
| 10 | `-4.8333522108e22 + 4.8517641700e22 i` | `-4.8333522108e22 - 4.8517641700e22 i` |
| 100 | `2.0032573598e19 + 8.3801109916e15 i` | `2.0032573598e19 - 8.3801109916e15 i` |
| 1000 | `2.6178634920e15 + 8.2654639678e10 i` | `2.6178634920e15 - 8.2654639678e10 i` |
| 10000 | `-6.4072125438e9 - 3.8911491215e4 i` | `-6.4072125438e9 + 3.8911491215e4 i` |
| 100000 | `-5244.5108786 - 0.0033281971 i` | `-5244.5108786 + 0.0033281971 i` |

Endpoint phases can look almost identical while endpoint magnitudes are far from
one and the accumulated phase is near four or five windings. Neither observation
justifies omitting the upper closure.

Independently sampled normalized rectangles directly measure that contribution:

| Real bounds | Top gamma | Bottom phase | Other three edges | Sum / 2pi | Calls |
|---|---:|---:|---:|---:|---:|
| `[-1,1]` | 0.2 | 16.191182474015 | −3.624811859656 | 2 | 788 |
| `[-10,10]` | 10 | 20.424154357934 | −7.857783743575 | 2 | 1526 |
| `[-1000,1000]` | 1000 | 25.132678082080 | −12.566307467721 | 2 | 1838 |
| `[-100000,100000]` | 100000 | 31.415925266686 | −18.849554652327 | 2 | 2610 |

Bottom gamma is 0.01 throughout. The latter three rectangles use strict controls;
all converge in three passes. Their refinements/depths are 315/10, 471/17 and
857/23 respectively. Each bottom phase agrees with the separate horizontal run
to 1e-8 rad, and bottom plus other-three equals the accepted total to 1e-10 rad.
The closure is becoming more negative over this sequence, not vanishing.
Expanding to `[-1e7,1e7] × [0.01,1e7]` fails at depth 32 after 170 calls; no
complete rectangular count is returned for that calculation.

A separate far-tail interval `[1e6,1e8]+0.01i` converges in 127 calls, three
passes, 6 refinements, depth 3. Its endpoints are approximately
`0.5438947554 + 7.34496108e-9 i` and
`0.9999448449 + 1.10246385e-14 i`; its finite phase is
`-1.35043682795e-8` rad. This supports the expected normalization trend at truly
large frequencies. It does not resolve the intervening path, certify uniform
upper-closure behavior, or establish a complete root count.

## 7. High-frequency unresolved behavior

Depth-48 diagnostic repeats at R=1e6 and R=1e7 also remain unresolved. The former
uses 136 calls and the latter 144 calls. The failure occurs near
`Re(z)=-4.2862028e5`, a scale comparable to the electron/combined plasma frequency,
well beyond proton-mode frequencies. Retained depth-48 samples show non-monotone
small changes in the normalized real determinant over extremely small frequency
steps, for example:

| Re(z), gamma=0.01 | Re(Fhat) | Im(Fhat) |
|---:|---:|---:|
| −428620.27883421659 | −5.2112327531e-12 | 8.5519704868e-13 |
| −428620.27883415972 | −7.4236255702e-12 | 8.5519717746e-13 |
| −428620.27883415617 | −5.2112368089e-12 | 8.5519698274e-13 |
| −428620.27883415262 | −2.9988380964e-12 | 8.5519639242e-13 |

This is consistent with evaluation precision/cancellation limiting local phase
resolution near a very small determinant. The precise source has not been
isolated or repaired; it is not attributed conclusively to WOFZ. No invalid or
overflowing return was observed in these traces. The existing WOFZ overflow flag
still is not propagated, as identified in A0. Finiteness does not guarantee an
accurate determinant.

Increasing sampling depth alone did not solve this case. No tolerance was relaxed,
no smoothing or determinant formula change was made, and the A1 kernel correctly
withheld acceptance. A separate focused determinant-conditioning investigation
is needed before using this low-gamma horizontal route through that frequency
range for production counts.

## 8. Timing and computational comparison

Timing uses `MPI_Wtime`: maximum rank time from the post-initialization driver
barrier through stop-command completion. This includes callbacks and sampling,
but excludes equilibrium preparation and MPI process launch. The summary also
records ASCII writing time and total driver elapsed time before summary output.
All direct comparisons use four ranks and the same final optimized executable.
The short serial benchmark measurements are recorded below.

| Calculation | Calls | MPI calculation (s) | Sample write (s) | Driver before summary (s) | Passes | Refinements | Depth | Status | Sample file |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Raw trial rectangle, output off | 630 | 0.029092 | 0 | 0.029127 | 3 | 91 | 6 | 0 | Off |
| Same rectangle, output on | 630 | 0.031218 | 0.001236 | 0.032486 | 3 | 91 | 6 | 0 | On |
| Horizontal R=10 | 511 | 0.028152 | 0 | 0.028202 | 3 | 198 | 10 | 0 | Off |
| Normalized rectangle R=10, top=10 | 1526 | 0.052572 | 0 | 0.052605 | 3 | 315 | 10 | 0 | Off |
| Horizontal R=1e5 | 1177 | 0.036402 | 0 | 0.036434 | 3 | 531 | 24 | 0 | Off |
| Normalized rectangle R=1e5, top=1e5 | 2610 | 0.091969 | 0 | 0.092002 | 3 | 857 | 23 | 0 | Off |

All rows are converged for their finite prescribed paths. The rectangles with
sample files off still retain edge diagnostics internally when normalized mode
is selected. The two horizontal rows have no accepted integer or verified
asymptotic closure.

Evaluation counts are the most reproducible comparison. On the R=10 finite region,
the horizontal bottom costs 511 calls versus 1526 for its independently sampled
normalized rectangle. At R=1e5 the comparison is 1177 versus 2610. Those savings
cannot replace the omitted edges: they supply essential closure phase. The raw
small rectangle costs only 630 calls and already provides a bounded count of two.

The diagnostic-on/off rectangle uses the same 630 calls, 3 passes, 91 refinements,
depth 6 and accepted phase. ASCII writing is a separate small cost, but these
short wall times fluctuate with host/MPI scheduling, so a precise speedup or
output-overhead percentage is not established. Initial bi-Maxwellian work is not
distributed like the numerical-VDF harmonic response; these are pilot timings,
not general performance estimates. No A1 optimization was attempted.

## 9. Error handling and regression coverage

Both final complete A2 suites passed. The combined validation also includes:

- Ordinary unstable/stable root frequency and residual checks; mapping and both
  scan modes. The final optimized build's nine scientific files compare byte-for-byte
  with the pre-change build. This also checks the narrow rank-0 Bessel guard fix.
- The existing ten-case input-validation suite (invalid first/second method-2
  settings, valid mixed/rel4/rel5 inputs, obsolete logmax, supported relativistic,
  invalid polynomial kind, invalid bi-Maxwellian relativity, invalid grid).
- Unmodified A1: 306 checks per debug/optimized build, including normalization,
  signed poles, repeated roots, coarse-loop aliasing, invalid values and budgets.
- A real returned-value fault: a copied adapter replaces the result **after the
  collective `disp` call** with quiet NaN. All four ranks make exactly one call;
  A1 returns status 2, no integer is accepted, and common finalization completes.
  Production sources are not patched for this test.
- Evaluation cap 10: status 4, all ranks make exactly 10 calls, no integer.
- Invalid sampling controls (budget 2): status 1, zero collective evaluations,
  workers receive the stop command and finalize normally.
- Invalid gamma=0 input: clear Nyquist-specific input error and communicator-wide
  abort before frequency evaluation, matching existing fatal-input practice.
- Missing output directory: a converged calculation reports summary I/O failure,
  releases all workers, finalizes, and returns a nonzero launcher status.
- Reserved operations, optional sample output, archived example inputs, finite
  horizontal reference phases and independently sampled closure agreement.

A1 status meanings are preserved: 0 converged, 1 bad input, 2 bad value, 3 evaluator
failure, 4 budget, 5 depth/coordinate resolution, 6 no pass/integer agreement.
Sampling failures return a nonzero operation exit after coordinated termination;
open sampling success always retains unverified closure and no integer. No
unresolved result is reported as a stable equilibrium.

The broad numerical-VDF/relativistic scientific suites were not rerun as physical
Nyquist validation; A2 does not claim those backends are certified for winding.
The existing documented broad-suite relativistic limitation is not described as
passing. Existing scratch harnesses that compile main have only had their module
build lists updated.

## 10. Limitations and next steps

The preferred production strategy remains normalized horizontal evaluation **if**
a reliable closure can be demonstrated for the physical determinant. A2 has not
established that prerequisite. For these fixtures, use the converged finite
rectangle as the reference capability. Two roots in a rectangle above gamma=0.01
is not a proof of the total number of unstable modes at every growth rate or
frequency.

Finite deterministic sampling provides numerical evidence, not a theorem excluding
arbitrary missed loops. Analyticity/conditioning of the evaluator and physical
boundary coverage remain separate from A1 sampling convergence. The tests use
A0's fixed harmonic controls; no general harmonic-truncation error bound, complete
root inventory, numerical-VDF holomorphy assessment, cold-pole treatment or
special-function error certification was added.

Recommended next work is a focused high-frequency determinant precision study
around the observed plasma-scale failure, followed by an explicitly justified
closure protocol. Preserve the existing backend unless that study supports a
specific reviewed correction. Separately benchmark evaluation reuse across A1
passes, less redundant sampling of smooth edges and closure-aware bounds, without
dropping independent convergence evidence. These are opportunities, not A2 changes.

### Forward design for options 4 and 5

The existing `scan(:)` stores `type_s`, ranges, linear/log steps, output/substep
counts and the derived increments; `scan_read` computes them and
`pass_instructions` already broadcasts them. Future option 4 can read the same
`scan_input_N` groups and follow their existing wavevector-path semantics, calling
the single-wavevector driver at each selected point. Future option 5 can use the
same two scan structures and plane ordering as `om_double_scan`, with the Nyquist
calculation as the per-wavevector operation. Neither should invoke root tracking
merely to obtain the next wavevector.

Those future drivers must update wavevector state on every rank and retain the
current `determine_nmax`, `split_processes` and Bessel rebuild sequence when kperp
changes. They can reuse the same frequency-path namelist and collective request
protocol. The present early reserved-operation rejection must be replaced only
when the relevant reader/loop is implemented and tested. No new scan reader,
wavevector path definition, plane loop, or frequency/wavevector parallelization
has been implemented in A2.
