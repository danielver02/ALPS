# ALPS Nyquist B1 — wavevector-line scans

Date: 2026-11-21 (America/Phoenix).

**Implemented and validated:** `scan_option=4` evaluates prescribed wavevector
paths using the existing Nyquist contour and optional growth-search algorithms.
The detailed `.nyquist.scan` and four-column `.nyquist.dat` files preserve valid
counts, converged estimates, retained unresolved brackets and explicit missing
values. Recoverable numerical failures do not prevent subsequent wavevectors.
Option 5 remains reserved. No physical backend was added.

## 1. Implementation and preserved components

| Source / routine | B1 change |
|---|---|
| `src/ALPS.f90` | Recognize option 4, read Nyquist controls once, enter the scan after shared initialization and before root/map execution, return to common MPI finalization; measure shared initialization. |
| `src/ALPS_io.f90`, `init_param`, `scan_read`, `alps_error` | Option 4 ignores guesses/maps, requires scans and uses the existing indexed scan reader. New-operation defaults and validation, finite/log-coordinate checks, initialized unused diff2; scoped magnitude-endpoint correction; diagnostic code 20 identifies invalid wavevector paths. |
| `src/ALPS_fns.f90`, `scan_wavevector`, `update_scan_wavevector`, `om_scan` | Extract the ordinary coordinate arithmetic unchanged into a shared routine. Reuse the original determine_nmax/split_processes/Bessel sequence through a helper. The touched Bessel guard uses defined `.not.proc0`, avoiding rank-0 sproc. Root algorithms and tracking decisions are unchanged. |
| `src/ALPS_nyquist_driver.f90` | Optional internal per-point result from `run_nyquist`, suppression of per-point summaries/history in scan mode, collective path driver, two aligned tables, aggregate summary and continuation after numerical failures. Existing option-3 call signature remains valid and filenames are unchanged. |
| `validation/make_nyquist_scan_cases.f90`, `check_nyquist_scan.f90` | Fortran path fixtures and independent table/aggregate/option-3/coordinate checks. |
| `validation/run_nyquist_scan.sh`, `run_nyquist_scan_state.sh` | Coordinate checkpoint, physical equivalence, scientific-output compatibility, allocation stress and failure gates; all copied variants/builds/outputs are external. |
| `validation/run_nyquist.sh` | Replace obsolete reserved-option-4 expectation with a missing-required-scan rejection; retain reserved-option-5 and all existing physics/legacy assertions. |
| `input.md`, `validation/README.md`, `tests/data/nyquist_scan/` | Input/output documentation, complete physical inputs, retained independent coordinate observations and reproduction instructions. |

No changes to `ALPS_com.f90`, `ALPS_var.f90`, susceptibility, `disp`, root-search
algorithms or analytic continuation were needed. Existing `pass_instructions`
already broadcasts all scanner fields. No new build dependency or library was
introduced; the existing module/build ordering suffices.

The A1 and A3 mathematical source hashes are identical to the pre-B1 snapshot:

- `ALPS_nyquist.f90`: `d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.
- `ALPS_nyquist_growth.f90`: `5eda3a5b6faf58ca54417b5bf26f4bd89abe6a2e47d1f8ac9978a73434d1a212`.

Linear growth searches retain strict absolute width tolerance; logarithmic growth
searches retain strict width/arithmetic-midpoint relative tolerance. Their trial
geometry, complement counting, independent A1 passes and exact caching are unchanged.

## 2. Wavevector-path conventions and coordinate checkpoint

One combined scan has `1+sum(ns)` points: the input wavevector and `ns` output
points per segment. Ordinary `om_scan` evaluates `ns*nres` tracking steps and
outputs steps divisible by `nres`; B1 computes exactly those output coordinates
without evaluating intermediate tracking steps. Segment endpoints are included.
The start of the next segment is already present and is not repeated. Explicitly
requested duplicate coordinates, including zero-length segments, retain their rows.

Types 0–4 retain their existing coordinate meanings: both components, angle at
fixed magnitude, magnitude at fixed angle, perpendicular component, and parallel
component. Angle uses the legacy `atan(kperp/kpar)` convention and degrees for the
input target. `swi` is the target perpendicular component of type 0. Linear and
log10 wavevector spacing use `swlog`, independently of growth spacing. Heating,
eigenfunction output and intermediate root tracking do not apply to Nyquist scans.

Before physical assertions, copied **pre-B1** ordinary scans were observed after
coordinate generation but before rebuilding/root tracking. All five types, both
wavevector spacings, ns=4/nres=3 and a seven-point consecutive path agree with B1
to 2e-15. Full-precision pre-B1 observations are retained as `*.coordinates` and
are checked independently of the shared helper in both compiler configurations.

Inspection found an existing type-2 reader inconsistency: after calculating a
magnitude endpoint it retains the old magnitude for preparing subsequent segments.
B1 corrects this only for option 4; ordinary option-1/2 bookkeeping is preserved.
A separate regression starts at (0.2,0.1), changes magnitude to 0.4, then changes
kperp to 0.5. The final kperp is 0.5 rather than the displaced target produced by
that old bookkeeping. This is an endpoint correction, not a new path definition.

## 3. Collective state and cache lifetime

Equilibrium input/preparation, distribution transfer, arrays and initial harmonic
setup run once. At every output point rank 0 broadcasts the wavevector to all
ranks. A perpendicular change invokes the original `determine_nmax`,
`split_processes`, then worker `determine_bessel_array`, which reallocates from
the current nlim. Pure parallel changes skip that sequence. A barrier precedes
contour evaluation. The unchanged `disp` reconstructs k-dependent wave terms and
passes current kperp/kpar into the analytic susceptibility on every call.

Each `run_nyquist` invocation constructs an empty raw-determinant cache. Its
lifetime spans that point's initial contour, growth strips and independent passes;
it never spans wavevectors. Exact hits avoid MPI commands. Misses invoke the
existing collective determinant on every rank. Actual per-rank calls are compared
at every point; logical requests equal physical calls plus cache hits. The
point-end command returns workers to the outer collective scan loop, not to
MPI_Finalize. MPI finalizes once after the complete path.

For the supported warm backend the outer ALPS harmonic limit is always one, so
physical k changes do not naturally change that assignment. A copied-source
allocation stress changes the outer limits to 3 and 5 at successive points while
retaining the same warm susceptibility. All three workers change harmonic ranges
and Bessel shapes; bounds are asserted and both scientific rows and evaluation
counts remain identical to the unmodified warm run. This validates the update
plumbing, **not** numerical-VDF physics or a new MPI decomposition.

## 4. Physical and single-wavevector equivalence

Complete fixtures are in `tests/data/nyquist_scan/`. Each scan point in the
reference, mirror-perpendicular, mirror-parallel and isotropic paths is independently
rerun with option 3 and identical controls. Accepted initial counts, A1/growth
statuses, bracket endpoints, midpoint, logical requests, physical calls and hits
agree exactly. No guesses, map inputs or root tracking are consumed; fixtures
intentionally set use_map=T/nroots=2 while omitting all guesses/maps.

The reference plasma remains beta_parallel=3 for both species, proton anisotropy
2, isotropic electrons, physical mass ratio, no drifts and vA/c=1e-4. At initial
(0.001,0.3), real bounds [-1,1], gamma bounds [0.01,0.2], linear tolerance 1e-4:
W=2; the bracket is `[0.0859814453125,0.08607421875]`, 11 steps,
8690 logical requests, 5005 determinant calls and 3685 hits. The following
kpar=0.31 point also has W=2 and agrees with its independent option-3 calculation.

The mirror plasma has beta_parallel=10 for both species, proton anisotropy 2,
isotropic electrons, physical mass ratio, zero drift and vA/c=1e-4. All rows below
use fixed real bounds [-0.05,0.05], gamma bounds [1e-4,0.08], depth 32 as in A2/A3,
logarithmic growth spacing and relative tolerance 1e-4. Each has W=1 and 17 steps:

| kperp | kpar | Final lower | Final upper | Midpoint estimate | Logical | Calls | Hits |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.05 | 0.04926554764026217 | 0.04926806022460729 | 0.04926680393243473 | 9142 | 6034 | 3108 |
| 0.3 | 0.05 | 0.03816898622648722 | 0.03817093287687182 | 0.03816995955167953 | 9208 | 6070 | 3138 |
| 0.4 | 0.05 | 0.02403493693052050 | 0.02403616273250910 | 0.02403554983151480 | 9354 | 6146 | 3208 |
| 0.3 | 0.04 | 0.03068412055757314 | 0.03068568547351640 | 0.03068490301554477 | 9264 | 6212 | 3052 |
| 0.3 | 0.06 | 0.04552341086866524 | 0.04552573260088915 | 0.04552457173477720 | 9170 | 6040 | 3130 |

Relative widths are approximately 5.099954015e-5. The kperp=0.3 point is an exact
segment endpoint and occurs once. Its bracket contains the independent ordinary
root growth `0.0381697286245976`. The complete existing mirror root/polarization,
reduced-anisotropy, isotropic and harmonic tests also pass, retaining the A3 mirror
identification rather than identifying it solely from growth. These brackets differ
from the earlier lower-bound-0.005 fixture because B1 uses the requested 1e-4 bound.

Expanding the mirror real bounds to [-0.2,0.2] leaves all three perpendicular-scan
counts equal to one; the expanded results also satisfy output-consistency gates.
This is evidence for adequacy of the selected real bounds, not an inventory
outside the contours. The isotropic control has accepted zero counts at all three
perpendicular points, no growth search and no midpoint; it uses 1356 calls total.

## 5. Output contract and required consistency gates

Option 4 automatically writes `.nyquist.scan`, `.nyquist.dat` and an aggregate
`.nyquist.summary` under `solution/<runname>`. The detailed table has 26 columns,
documented in `input.md`, including full-box winding/status, retained bracket,
converged midpoint, both widths, per-point logical/physical/cache accounting,
timing, sampling diagnostics, rank call extrema and evaluator status.

The scientific table has **exactly four numerical columns**:

```text
kperp  kpar  Wn  gamma_max
```

Wn is the signed accepted original-box winding, never a strip winding. gamma_max
is the arithmetic midpoint of a successful bracket, not a dispersion root or the
input top boundary. A valid zero count has `0 NaN`; an unresolved initial contour
has `NaN NaN`; a positive count with disabled/failed/incomplete growth retains its
integer Wn and has NaN for gamma_max. Available unresolved bracket endpoints stay
in the detailed table. All numerical values use sufficient double precision.

The Fortran checker parses the two files independently at every tested physical,
coordinate and failure case. It checks exactly four compact columns, equal row
counts/order/coordinates, original winding, midpoint/NaN semantics, rank and cache
accounting, and aggregate counts/calls. The selected physical, failure and segment
compact files are also byte-identical between debug and optimized builds.

Nyquist line scans write no per-point summary/history. Optional samples have distinct
`.nyquist.point_NNNNNN.samples` names and describe only each initial contour.
The three-point sample test checks all files and the absence of growth histories.
Option 3 retains its filenames: growth history and ordered samples remain
byte-identical to the pre-B1 executable, and summaries are identical after excluding
only timing fields.

## 6. Failure recovery and aggregate status

The copied-adapter nonfinite injection at the interior mirror point produces:

```text
0.2  0.05   1    0.04926680393243473
0.3  0.05  NaN   NaN
0.4  0.05   1    0.02403554983151480
```

Both tables retain all three rows. The first and last calculations converge;
all four ranks make identical actual determinant calls at every point. The
aggregate reports three evaluated points, one unresolved initial contour and a
nonzero overall status, followed by coordinated finalization. It never reports
an entirely successful scan or substitutes stability for failure.

A one-iteration limit preserves positive original counts and finite valid brackets
at every mirror point, with NaN midpoints and failed growth statuses. A ten-request
A1 budget leaves unresolved counts, not zeros, and proceeds through all points.
Tests reject ns=0, nres=0, unsupported types, absent scans/groups, nonfinite/zero
parallel coordinates and incompatible horizontal geometry. Missing output
directories fail collectively before determinant work. Option 5 still reports
reserved/unsupported and finalizes cleanly without requesting species or scans.

The aggregate distinguishes positive/zero/negative accepted windings, unresolved
initial contours, successful growth searches, failed growth calculations and failed
searches after a valid positive count. Numerical point failures continue; input,
initialization and output failures stop appropriately. A completed path with any
unresolved requested calculation returns nonzero. Diagnostic-output I/O failure
cannot guarantee complete aligned tables and is explicitly fatal to the scan.

## 7. Performance observations

Sequential four-rank runs used the identified optimized build, one BLAS/OpenMP
thread and sample storage disabled. Initial input/MPI launch are excluded from the
reported initialization interval, which measures shared preparation after input;
scan wall time includes path preparation, point calculations and table writing.

| Path | Points | Logical requests | Actual calls | Hits | Scan seconds | Sum determinant seconds | Shared init seconds | Rebuilds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Mirror kperp 0.2–0.4 | 3 | 27704 | 18250 | 9454 | 2.2232 | 2.2025 | 0.000616 | 2 |
| Same range, denser | 5 | 46156 | 30396 | 15760 | 3.9092 | 3.8775 | 0.000860 | 4 |
| Same range, denser | 9 | 83032 | 54942 | 28090 | 6.7769 | 6.7244 | 0.000853 | 8 |
| Mirror kpar 0.04–0.06 | 3 | 27642 | 18322 | 9320 | 2.3838 | 2.3636 | 0.000778 | 0 |
| Isotropic kperp control | 3 | 1610 | 1356 | 254 | 0.1546 | 0.1526 | 0.000990 | 2 |

For the three perpendicular points, determinant calls are 6034/6070/6146 and
roughly one third of logical requests are cache hits. Aggregate point times are
2.2224 seconds; additional scan preparation/output costs are small here. Recorded
update/barrier time sums are about 12–44 microseconds across these short scans.
These tiny warm-placeholder rebuilding costs are timer/load dependent; their
reliable structural distinction is two/four/eight rebuilds versus zero for pure
parallel changes. Cost grows approximately with point count in this modest test.

Per-point timings and sampling diagnostics from that same sequential run:

| kperp at kpar=0.05 | Point seconds | Determinant seconds | Passes | Refinements | Maximum depth |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.5956 | 0.5887 | 54 | 539 | 13 |
| 0.3 | 0.7311 | 0.7245 | 54 | 572 | 14 |
| 0.4 | 0.8958 | 0.8892 | 54 | 645 | 20 |

These are observations, not numerical-VDF performance forecasts: the warm analytic
response does not distribute work over harmonics like numerical distributions.
No A1, strip or cache optimization was added. Completed edge phases are not spliced;
independent A1 sampling remains mandatory at every strip and every point.

Benchmark artifacts are `/tmp/alps-nyquist-b1/benchmark/`, executable SHA-256
`15327da7de02f2cccb31641abcf6c70d06d65adbb23d556b7fd656f074b0abde`.
Only subsequent invalid-input diagnostic wording changed after that timed build;
physical code, paths, sampling and accounting are identical.

## 8. Regression execution and provenance

Complete suites were run with four local ranks and `--oversubscribe`:

```sh
bash validation/run_nyquist.sh
bash validation/run_nyquist_growth.sh
bash validation/run_nyquist_log.sh
bash validation/run_nyquist_mirror.sh
bash validation/run_nyquist_scan.sh
```

Both optimized `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow` and debug
`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow` configurations passed.
A1 retains 306 checks in each independent debug/-O3 build, A3 retains 10056
math/cache checks, and logarithmic math retains 72. Existing failures, physical
root/count/normalization/closure tests and mirror identification remain intact.
Nine ordinary scientific root/stable-root/map/consecutive-scan/plane-scan files
remain byte-identical to the identified pre-B1 build in both configurations.

Set `ALPS_REFERENCE_EXE` for the legacy comparison and
`ALPS_B1_REFERENCE_SOURCE` for an actual pre-B1 coordinate/source comparison.
This session used `/tmp/alps-nyquist-b1/original` and its separately built executable.
Retained coordinate observations additionally permit independent comparisons
without that temporary snapshot. GNU Fortran 11.4.0, Open MPI 4.1.2 and existing
BLAS/LAPACK were used. Builds include source/flags/executable manifests and copied
coordinate, fault and allocation-stress variants have separate provenance.

Final B1 artifacts:

- Optimized: `/tmp/alps-nyquist-b1.iBD0rzMS`, executable SHA-256
  `308233d31c1a4203eec87ff81ce4b736f97bd31d8f2732c588e8b976793848a8`.
- Debug: `/tmp/alps-nyquist-b1.gMfJqD2x`, executable SHA-256
  `9bca0316d953be9be27c766df19aaba922add681e935c144d54a855e28f08225`.
- Final in-suite allocation stresses: `/tmp/alps-nyquist-b1-state.kZtMdnmY`
  and `/tmp/alps-nyquist-b1-state.Uynw2QFi`.
- A1/A2: `/tmp/alps-nyquist-a2.avywyj3Z` and `/tmp/alps-nyquist-a2.CTCML7va`.
- A3: `/tmp/alps-nyquist-a3.Zt8pw0Le` and `/tmp/alps-nyquist-a3.3x6DcqSM`.
- Log: `/tmp/alps-nyquist-log.m0lxkeUU` and `/tmp/alps-nyquist-log.VjPuUkSV`.
- Mirror: `/tmp/alps-nyquist-mirror.leCWJgh5` and `/tmp/alps-nyquist-mirror.UNcBEHsa`.

Matrix logs are under `/tmp/alps-nyquist-b1/`; each suite prints its own artifact
path. Temporary files are not durable archives; fixtures, independent coordinate
observations, test sources and measured results are retained in the repository.
All generated builds and scientific outputs were isolated; unrelated working-tree
modifications and existing scientific results were preserved.

## 9. Remaining limits and option-5 implications

Counts and growth brackets apply only within the fixed input real-frequency and
upper-growth boundaries, above the selected lower threshold. A valid zero is not
an unconditional stability certificate. Tight strips may still exhaust A1 limits;
finite sampling is numerical convergence evidence, not exclusion of arbitrary
unresolved loops. The existing high-frequency determinant-conditioning and WOFZ
flag limitations are unchanged. Validation remains warm, nonrelativistic analytic
bi-Maxwellian only; real numerical-VDF harmonic changes need physical validation
before that support is extended.

The per-point adapter, collective wavevector update and explicit cache lifetime
provide a basis for option 5, but plane semantics must separately reuse/validate
the existing two-scan input/grid conventions, row ordering and aggregate output.
No plane execution, new input switch or wavevector-level parallelization was
implemented. B1 preserves the existing single-wavevector mathematics and is ready
for subsequent plane-scan design within those limits.
