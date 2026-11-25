# ALPS Nyquist B2 — wavevector-plane scans

Date: 2026-11-21 (America/Phoenix).

**Implemented and validated.** Option 5 traverses the existing two-scan geometry
through the shared B1 execution machinery. The nine-point mirror grid, its reversed
ordering, every independent option-3 point and all six option-4 rows/columns agree
in counts, brackets and evaluation accounting. Isotropic controls give accepted
zero counts within the prescribed box. Recoverable failures retain all grid rows
and do not corrupt later calculations.

**The MPI resonance decomposition is unchanged.** `ALPS_fns.f90` and
`ALPS_com.f90` are byte-identical to the pre-B2 snapshot, including `disp`,
`determine_nmax`, `split_processes`, `determine_bessel_array` and the B1 update
helper. No communicator, wavevector-worker group, scheduling scheme or additional
parallelism was introduced. Validation remains warm analytic bi-Maxwellian only.

## 1. Source changes and shared implementation

| File / routine | B2 change |
|---|---|
| `src/ALPS.f90` | Read Nyquist controls for option 5 and dispatch both options 4/5 to `run_nyquist_scan` after the existing shared preparation; bypass root/map tracking and use common finalization. |
| `src/ALPS_io.f90`, `init_param`, `scan_read` | Require n_scan=2, distinct types 1–4, bypass guesses/maps, reuse existing input/broadcast structures and extend B1 defaults/validation/magnitude bookkeeping to option 5. Legacy option-1/2 branches are unchanged. |
| `src/ALPS_nyquist_driver.f90` | Extract B1 coordinate preparation into `nyq_scan_coordinates(plane,coordinates,segments,indices)`, add plane coordinates/indexing and conditional metadata/appended indices. One shared update/evaluate/record/finalize loop remains authoritative for options 4/5. |
| `validation/make_nyquist_plane_cases.f90` | Complete physical and coordinate fixture generator using existing scan inputs. |
| `validation/observe_nyquist_plane.f90`, `observe_nyquist_plane.awk`, `check_nyquist_plane_coordinates.f90`, `run_nyquist_plane_coordinates.sh` | Independent coordinate-only gate using copied ordinary plane traversal and retained pre-B2 observations. |
| `validation/check_nyquist_plane.f90`, `run_nyquist_plane.sh` | Plane indices, tables, independent points/lines, order, coverage, assignment auditing and recoverable/fatal failures. |
| `validation/check_nyquist_scan.f90` | Allow the established option-3 input emitter to consume option-5 templates as well as option-4 templates; existing B1 assertions are retained. |
| `validation/run_nyquist.sh` | Replace the obsolete reserved-option-5 test with rejection of its missing required two scans; all physics and legacy assertions remain intact. |
| `input.md`, `validation/README.md`, `tests/data/nyquist_plane/` | Current input/output documentation, runnable mirror/reverse/isotropic inputs and captured coordinate observations. |

No new production module, numerical dependency, namelist parameter or build-order
change was needed. `run_nyquist` is unchanged and remains the single-wavevector
calculation. A1 phase sampling, A3 bisection, exact raw caching, normalization,
susceptibilities and analytic continuation were not modified.

Unchanged source SHA-256 values:

- `ALPS_fns.f90`: `567e05bad770489d578c0024b4230376f1e000c2f04b7e84fc2f7d065d446be5`.
- `ALPS_com.f90`: `989cfbb8fd3a8252e195d868a93a09552cfbce3e4f1172f93c6a127e947e5cc8`.
- `ALPS_nyquist.f90`: `d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.
- `ALPS_nyquist_growth.f90`: `5eda3a5b6faf58ca54417b5bf26f4bd89abe6a2e47d1f8ac9978a73434d1a212`.

## 2. Existing geometry and pre-integration checkpoint

Option 5 requires exactly two different scan types in 1–4. Scan 1 is outer;
scan 2 varies fastest. There are `(ns_1+1)*(ns_2+1)` output points, including
both directions' initial and final indices. Output tracking indices are
`0,nres,2*nres,...,ns*nres`; intermediate root-tracking steps are not evaluated.
Starts are retained for every outer row rather than deduplicated as line-segment
boundaries. Coordinates use the existing `scan_wavevector` arithmetic, including
the legacy `atan(kperp/kpar)` angle convention and base-10 wavevector spacing.

Before enabling option-5 determinant evaluation, **22 coordinate comparisons**
passed against actual copied pre-B2 `om_double_scan` traversal. These cover:

- Cartesian and reversed Cartesian order, including all endpoints/row starts.
- Linear, logarithmic, reversed-logarithmic and mixed spacing.
- Angle–magnitude, logarithmic polar, magnitude–angle and perpendicular–angle combinations.
- nres=2/3 and a magnitude–component transition.

The copied observer bypasses root/Bessel work but retains both original coordinate
loops and row-reset expressions. It does not call the new coordinate helper.
Observations were captured before option 5 was integrated, at full double precision,
and retained under `tests/data/nyquist_plane/*.coordinates`. The standalone gate
compares both live observations and the new list with those independent records
to 2e-15 before the physical suite evaluates any determinant.

Two existing geometry details are explicit:

1. The reader computes each direction's increments once, preparing direction 2
   after direction 1's endpoint. The same inner increment is then translated from
   every outer wavevector. Non-Cartesian rows can consequently have different
   absolute inner endpoints. B2 preserves this; it does not silently recalculate
   every inner increment to a common absolute target.
2. The B1 magnitude-endpoint correction applies to option 5, as requested. Both
   pre-B2 legacy-option-2 and pre-B2 B1-corrected-reader plane observations were
   captured independently. For the magnitude–perpendicular example, legacy
   bookkeeping ends at kperp approximately 0.6422322703, while the corrected
   endpoint is the requested 0.45. Both observations are retained, and the
   corrected geometry is the option-5 reference. Options 1/2 are unchanged.

An additional 2-by-3 rectangle verifies non-square dimensions and ordering.
Arbitrary non-Cartesian scan directions are not assumed to commute.

## 3. MPI state, assignments and cache lifetime

Equilibrium/distribution preparation and initial resonance setup occur once.
At each selected point, the same B1 loop broadcasts the wavevector on
MPI_COMM_WORLD, invokes `update_scan_wavevector` if kperp changed, synchronizes,
and calls `run_nyquist`. The existing helper calls the unchanged limit,
assignment and worker Bessel routines. Pure parallel changes skip rebuilding.
No additional equilibrium initialization is performed per frequency or point.

All actual `disp` calls remain collective with the original reductions,
determinant broadcast and barrier. Rank 0 controls each whole point and its growth
strips; workers return from each point's request loop to the existing outer scan
loop. MPI_Finalize is reached once after the full grid. Cache hits send no MPI
request. Every point constructs a fresh exact-frequency raw-value cache; only its
own initial box, growth strips and independent A1 passes share cached values.
No phase splicing, approximation or cross-wavevector reuse was added.

The Cartesian kperp-outer grid requires two post-initialization rebuilds; its six
pure-parallel transitions require none. The reversed grid requires eight because
kperp changes within each row and resets between rows. Independently evaluated
kpar line scans require zero rebuilds. Pointwise min/max actual calls agree on
all four ranks in every physical and failure test.

A copied-driver audit writes each rank's assignment at every point and checks
worker Bessel allocation bounds against nlim. All three workers have nine audit
records, retaining their existing warm-backend species/harmonic assignments.
Every rank observes the same two rebuild transitions in the forward grid.
The production assignment routines are untouched even in this audit; the injected
invalid value is applied on rank 0 only after the existing collective `disp`.
B1's existing copied-source allocation-stress regression also passes unchanged.

## 4. Nine-point mirror-equilibrium results

Both species are warm analytic (`use_bM=T,relat=F`), beta_parallel=10, density 1,
physical electron/proton mass ratio, no relative drift, vA/c=1e-4. Proton anisotropy
is 2 and electrons are isotropic. The fixed rectangle is
`[-0.05,0.05] x [1e-4,0.08]`; log growth spacing uses relative tolerance 1e-4.
Fixtures retain explicit depth 32 from A2/A3 and all remaining A1 defaults.

All nine measured original-box windings are one. Every search converges in
17 steps, 54 whole-contour passes, with relative width about 5.099954015e-5.
The other eight counts/growth rates were calculated, not prescribed as targets.
Full-precision values and diagnostics are in the scan files; rounded here:

| kperp | kpar | Wn | Lower growth | Upper growth | Midpoint estimate | Actual calls |
|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.04 | 1 | 0.03982137883949 | 0.03982340976329 | 0.03982239430139 | 6094 |
| 0.2 | 0.05 | 1 | 0.04926554764026 | 0.04926806022461 | 0.04926680393243 | 6034 |
| 0.2 | 0.06 | 1 | 0.05835490863881 | 0.05835788478821 | 0.05835639671351 | 6058 |
| 0.3 | 0.04 | 1 | 0.03068412055757 | 0.03068568547352 | 0.03068490301554 | 6212 |
| 0.3 | 0.05 | 1 | 0.03816898622649 | 0.03817093287687 | 0.03816995955168 | 6070 |
| 0.3 | 0.06 | 1 | 0.04552341086867 | 0.04552573260089 | 0.04552457173478 | 6040 |
| 0.4 | 0.04 | 1 | 0.01938292292639 | 0.01938391147175 | 0.01938341719907 | 6192 |
| 0.4 | 0.05 | 1 | 0.02403493693052 | 0.02403616273251 | 0.02403554983151 | 6146 |
| 0.4 | 0.06 | 1 | 0.02855365428634 | 0.02855511054671 | 0.02855438241652 | 6248 |

The central count is explicitly checked as one. Its bracket contains the
independently established ordinary-root growth 0.0381697286245976. The existing
mirror root, polarization, harmonic and anisotropy-identification suite was rerun
in both configurations; B2 does not infer branch identity solely from positive
growth. The table reports bounded estimates, not newly located dispersion roots
or independent polarization identifications at the other eight points.

## 5. Scalar, line, reverse-order, stable and coverage comparisons

Every point in the mirror, reversed and isotropic grids is independently run with
option 3: 27 scalar comparisons. Full-box winding, initial A1 status, growth status,
bracket endpoints, converged midpoint, logical requests, actual calls and hits
agree exactly. The fixtures intentionally set ordinary root/map selectors while
providing no guesses/maps, verifying their bypass.

All three Cartesian rows and all three columns are also evaluated through option 4.
All 18 line-point results match the plane by physical wavevector, including exact
brackets and evaluation accounting. Reversing the Cartesian scan definitions
produces the same nine physical results, matched by coordinates rather than row
index. This statement is not generalized to arbitrary mixed scan types.

The isotropic proton control has nine accepted zero counts, zero growth searches,
NaN estimates, zero unresolved points and scan_status=0. Every point agrees with
its independent scalar result. This excludes roots only in the selected box.

Expanded rectangles `[-0.2,0.2] x [1e-4,0.16]` were evaluated at points 3, 5 and 7:
the largest and smallest measured growth rates and the central reference.
All three expanded counts remain one and their growth brackets overlap the original
brackets. No count transition or apparent stability boundary occurs in the unstable
nine-point grid. The selected expansion supports coverage there; it does not
establish a complete upper-half-plane inventory or justify automatic bounds.

## 6. Output format and independent parsing

The B1 filenames are retained under `solution/<runname>`:

- `.nyquist.scan`: first 26 B1 columns unchanged; append columns 27/28,
  `outer_index,inner_index`, one-based. For a plane, column 2 denotes the outer
  row. Headers describe dimensions, types, order, fixed frequency domain,
  equilibrium, normalization and growth/sampling controls.
- `.nyquist.dat`: exactly four numerical columns `kperp kpar Wn gamma_max`,
  outer-major/inner-fastest, with short comment headers and no blank separators.
- `.nyquist.summary`: existing aggregate counters/timings plus grid dimensions,
  scan types and one-based index/order interpretation.

Wn is the accepted **original** full-box winding, never a strip count. gamma_max
is the arithmetic midpoint of a successfully converged growth bracket, not the
input top boundary. Zero counts use `0 NaN`; unresolved initial contours use
`NaN NaN`; valid positive counts with disabled, failed or incomplete growth use
`Wn NaN`. Valid unresolved brackets remain in detailed output without a midpoint.
Optional samples retain flat `.point_NNNNNN.samples` names. Per-point full summaries
and growth-history files are not automatically created.

Fortran independently parses both tables, verifies exactly 4/28 numerical columns,
row order/count, one-based grid indices, dimensions, original counts, midpoint/NaN
conventions, cache/rank accounting and aggregate counters. The mirror, reversed,
isotropic and fault compact tables are byte-identical between debug/optimized
builds. Option-4 row/column scientific `.dat` files remain byte-identical to the
pre-B2 executable; detailed/summary files agree after excluding only timing values.
The full B1 suite also verifies option-3 growth histories and samples byte-for-byte.

## 7. Failure handling and recovery

The copied-driver NaN injection at central point (outer,inner)=(2,2) occurs only
after its actual collective determinant call. The entire nine-row grid is retained:
8 accepted positive counts/growth estimates, 1 unresolved initial contour, no
spurious zero count, overall status 1 and normal finalization. The failed point
has NaN winding and estimate. All unaffected cells, including every subsequent
cell, match the normal grid exactly in statuses, count, bracket and call accounting.

The failed point makes one actual call per rank. Aggregate calls are 49025,
exactly 55094 minus the normal central 6070 plus that one failed call. This also
checks that a preceding point's cached frequency could not hide the injection.
The assignment audit retains all nine records per rank, with Bessel bounds valid.

A one-iteration limit retains positive initial counts and valid brackets at all
nine points, reports failed searches and emits no midpoint. A ten-request A1
budget retains all rows with unresolved counts, never manufactured zeros.
A rectangular six-point disabled-growth grid retains valid original windings and
NaN estimates. Tests reject wrong scan counts, equal types, type 0, ns/nres=0,
nonfinite/zero-parallel coordinates and horizontal geometry. Missing output
storage terminates collectively before determinant evaluation.

Recoverable numerical failures continue; invalid global input, MPI failure and
output errors terminate appropriately. Existing B1/A3 status definitions and
shared finalization are unchanged. Unwritable/broken output cannot guarantee
complete tables and is fatal to the scan, rather than reported as physical stability.

## 8. Isolated performance measurements

All timings below use the identified optimized executable, four local ranks,
`--oversubscribe`, one BLAS/OpenMP thread and samples disabled. Jobs were run
sequentially after the regression jobs ended. The physical range and frequency/
growth controls are identical; only grid density/order or proton anisotropy changes.

| Grid / case | Points | Logical requests | Actual calls | Cache hits | Hit fraction | Scan seconds | Determinant seconds | Rebuilds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 x 2 mirror | 4 | 36966 | 24592 | 12374 | 0.33474 | 3.8778 | 3.8499 | 1 |
| 3 x 3 mirror | 9 | 83104 | 55094 | 28010 | 0.33705 | 8.9192 | 8.8605 | 2 |
| 5 x 5 mirror | 25 | 230584 | 152640 | 77944 | 0.33803 | 23.8670 | 23.7057 | 4 |
| 3 x 3 reversed | 9 | 83104 | 55094 | 28010 | 0.33705 | 8.4076 | 8.3413 | 8 |
| 3 x 3 isotropic | 9 | 4838 | 4074 | 764 | 0.15792 | 0.5034 | 0.4988 | 2 |

All these grids converge. Cost is approximately proportional to point count in
this modest warm test. Reversing order changes rebuild frequency but not physical
calls or results. Its slightly lower measured time despite more rebuilds is not
an optimization claim: warm rebuilding is tiny and timing/load effects dominate
that difference. Sum update/barrier times are 10.6, 23.4, 59.4, 66.2 and 32.6
microseconds respectively. These include synchronization, not isolated Bessel
routine timings. The shared initialization intervals are about 0.00062, 0.00105,
0.00074, 0.00111 and 0.00070 seconds, excluding MPI launch/initial input parsing.

The nine-point pointwise timings, sampling counts and cache statistics are retained
in the benchmark `.nyquist.scan`; total point time is 8.9178 seconds, close to the
8.9192-second full scan. All 17-step mirror searches use the unchanged independent
A1 checks. No sampling, edge-phase or cross-point-cache optimization was added.

Per-point details in the science table's outer-major order:

| Point | Nyquist seconds | Determinant seconds | Logical requests | Hits | Refinements | Max depth |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.7678 | 0.7615 | 9220 | 3126 | 578 | 16 |
| 2 | 0.8008 | 0.7945 | 9142 | 3108 | 539 | 13 |
| 3 | 0.8008 | 0.7939 | 9100 | 3042 | 518 | 13 |
| 4 | 0.9311 | 0.9251 | 9264 | 3052 | 600 | 15 |
| 5 | 0.9448 | 0.9389 | 9208 | 3138 | 572 | 14 |
| 6 | 0.9545 | 0.9485 | 9170 | 3130 | 553 | 14 |
| 7 | 1.3752 | 1.3677 | 9380 | 3188 | 658 | 15 |
| 8 | 1.1399 | 1.1338 | 9354 | 3208 | 645 | 20 |
| 9 | 1.2029 | 1.1966 | 9266 | 3018 | 601 | 14 |

The 5-by-5 grid also has 25 accepted count-one, converged-growth points; maximum
depth is 20, below the unchanged fixture limit 32.

Warm timings are not forecasts for numerical distributions: their susceptibility
and resonance workload differ. Benchmark artifacts: `/tmp/alps-nyquist-b2/benchmark/`.

## 9. Complete regressions and provenance

The following complete suites passed in optimized and debug configurations:

```sh
bash validation/run_nyquist.sh
bash validation/run_nyquist_growth.sh
bash validation/run_nyquist_log.sh
bash validation/run_nyquist_mirror.sh
bash validation/run_nyquist_scan.sh
bash validation/run_nyquist_plane.sh
```

Optimized flags: `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`.
Debug flags: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`.
A1 additionally runs its independent debug/-O3 comparison: 306 checks each.
A3 retains 10056 math/cache checks and log math 72. Nine ordinary root, stable-root,
map, line-scan and plane-scan scientific files remain byte-identical to the
identified pre-B2 build in both configurations. Existing option-3 physics,
normalization, horizontal closure diagnostics and linear 11-step reference
remain unchanged; option-4 formats and scientific results are preserved.

Two initial overlapping debug-suite attempts stopped during MPI startup before
producing Nyquist results. Their same identified executables passed isolated case
retries, and fresh complete mirror/B1 debug suites subsequently passed sequentially.
No production source, numerical tolerance or MPI/decomposition setting was changed
to obtain those reruns. These incomplete multi-job attempts are not counted as
successful suites; the isolated completed artifacts are listed below.

| Suite | Optimized artifacts | Debug artifacts |
|---|---|---|
| A1/A2, including legacy byte comparisons | `/tmp/alps-nyquist-a2.7kr2i4x1` | `/tmp/alps-nyquist-a2.jigal91G` |
| A3 | `/tmp/alps-nyquist-a3.Hx1pA5Ry` | `/tmp/alps-nyquist-a3.w7ixPqyI` |
| Log | `/tmp/alps-nyquist-log.zlsmYHwP` | `/tmp/alps-nyquist-log.tnaBrzv4` |
| Mirror | `/tmp/alps-nyquist-mirror.zn66NDNv` | `/tmp/alps-nyquist-mirror.fISavryJ` |
| B1 | `/tmp/alps-nyquist-b1.LDjIk2mD` | `/tmp/alps-nyquist-b1.Jt3kIeD0` |
| B2 | `/tmp/alps-nyquist-b2.KccMsYAh` | `/tmp/alps-nyquist-b2.fIo7VI93` |

Source/compiler/executable manifests are in each build. Final B2 executable hashes:

- Optimized/benchmark: `bad038ec4beb470963872776b73db57d0eed8f8b07fc98f37192a099bcf42c98`.
- Debug: `bf5fc4433b957fcda7e3e01955fe75235b218c05d3e55f3657bf77737df31570`.
- Pre-B2 reference: `403d44e50899e5fe0a6ac3187dec07031b7f21a9017393b14760b9dc2e73f3fb`.

GNU Fortran 11.4.0, Open MPI 4.1.2 and existing BLAS/LAPACK were used. The actual
pre-integration coordinate checkpoint is `/tmp/alps-nyquist-b2/checkpoint/` and
its log is `/tmp/alps-nyquist-b2/checkpoint.log`. The pre-B2 source snapshot/build
is `/tmp/alps-nyquist-b2/original` / `baseline`. Set `ALPS_B2_REFERENCE_SOURCE`
for copied-source coordinate observations and `ALPS_REFERENCE_EXE` for legacy/B1
byte comparisons. The retained coordinate files also provide an independent
reference when temporary snapshots are unavailable.

Logs, isolated retry records and final recovery assertions are under
`/tmp/alps-nyquist-b2/`. Temporary directories are not durable archives; fixtures,
captured coordinates, validation code and measured tables are repository deliverables.
All builds/calculated outputs were external; existing scientific outputs and
unrelated working-tree modifications were preserved.

## 10. Remaining limitations

Results describe the original finite frequency domain, above its lower threshold,
not all unstable roots or the global fastest mode. No automatic frequency-bound
selection or instability-boundary interpretation was added. Non-Cartesian geometry
retains legacy translated increments and angle conventions rather than introducing
new grid semantics. Finite phase sampling remains convergence evidence rather than
a completeness proof for arbitrary functions; near-boundary roots can remain
unresolved. Determinant conditioning and the unpropagated WOFZ overflow flag remain
existing limitations, separate from A1 convergence.

Only warm nonrelativistic analytic bi-Maxwellians were validated. Natural warm
outer harmonic assignments remain fixed; physical numerical-VDF assignment changes
require separate validation before extending backend support. No change to the
resonance decomposition was needed, and none was made. B2 provides the requested
plane capability without adding wavevector/frequency parallelism or further
winding/growth-search optimizations.
