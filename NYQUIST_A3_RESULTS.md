# ALPS Nyquist Stage A3 — maximum growth with cached closed strips

Date: 2026-11-20 (America/Phoenix).

**Implemented and validated at one prescribed wavevector.** A converged initial
count of two for the A2 unstable equilibrium leads to the maximum-growth bracket

`0.0859814453125 < gamma_mode,max < 0.08607421875`

within the original real-frequency bounds `[-1,1]` and imaginary bounds
`[0.01,0.2]`. Eleven bisections give width `9.27734375e-5`, smaller than the requested
`1e-4`, and contain the independently reproduced ordinary-solver growth rate
near `0.085986`. The midpoint `0.08602783203125` is an estimate, not a calculated
dispersion root. The two roots at signed real frequencies have equal growth and
are handled by counts jumping directly between two and zero.

The cache reduces the search from 8690 logical A1 requests to **5005 collective
frequency evaluations**, all executed by every MPI rank. Against an uncached
full-original-rectangle bisection baseline, it saves 2653 of 7658 evaluations
(34.6%). Wall-time benefits are less conclusive for this cheap analytic backend;
the measured comparison is reported below without assuming uniform cost per call.

## Implementation and unchanged behavior

| File | A3 work |
|---|---|
| `src/ALPS_nyquist_growth.f90` | New Fortran module: closed-strip bisection, search/result/history types, exact-frequency raw-value cache, growth status descriptions and compact history export. No ALPS state or MPI. |
| `src/ALPS_nyquist_driver.f90` | Read three new controls; default gamma_min=1e-4; enforce explicit real/upper bounds; keep workers available throughout the search; cache-aware synchronous callback; separate physical/logical counts and search diagnostics. |
| `src/Makefile.am`, `src/Makefile.in` | Add growth/cache module and explicit compile dependencies. Only the src template was regenerated. |
| `validation/build_nyquist.sh` | Compile the new module before the adapter; also supports the saved pre-A3 sources. |
| Existing numerical-baseline, relativistic-production and B01 scratch builders | Add the new Fortran module to existing object/compile lists; their test logic is unchanged. |
| `validation/test_nyquist_growth.f90` | Standalone analytic growth, failure and exact-cache tests. |
| `validation/make_nyquist_growth_cases.f90`, `check_nyquist_growth.f90` | Fortran physical fixture generation, bracket/history/accounting verification and independent strip generation/comparison. |
| `validation/run_nyquist_growth.sh` | Four-rank A3 suite, independent raw/normalized strips, uncached baselines, copied-adapter fault injection and invalid controls. |
| `input.md`, `validation/README.md`, `tests/data/nyquist_growth/` | Controls, interpretation, examples and test instructions. |

`src/ALPS_nyquist.f90` is **unchanged**, SHA-256:
`d6f54b3f98c2c9f3eefe9b3d19c0b410de19319fe72c7527f4b2e4f031e0fe31`.
The A1 sampling options and defaults remain unchanged. Direct comparisons with
the pre-A3 source snapshot also confirm no changes to `ALPS.f90`, `ALPS_io.f90`,
`ALPS_var.f90`, `ALPS_com.f90`, `ALPS_fns.f90`, `ALPS_fns_rel.f90`, `ALPS_NHDS.f90`
or `ALPS_analyt.f90`. Thus startup, susceptibility, determinant, root search,
continuation and harmonic/MPI assignment remain the A2 implementations.

Option 3 and the existing input wavevector are retained. Ordinary single-contour
operation remains the default; its caching is disabled, its existing summary and
optional sample fields retain their meanings, and its A2 evaluation counts remain
unchanged. Horizontal diagnostics remain available, but maximum-growth requests
with horizontal geometry are rejected. Options 4/5 remain reserved. All physics
restrictions remain warm, nonrelativistic analytic bi-Maxwellians; no other backend
is claimed as validated.

The new namelist controls are:

```fortran
nyq_find_max_growth = F
nyq_growth_tolerance = 1.d-4
nyq_growth_max_iterations = 32
```

`nyq_gamma_min` now defaults to `1.d-4`. The original real bounds and rectangular
upper-growth bound are mandatory, finite and ordered. Missing real bounds are
explicitly detected rather than silently defaulting to zero; a supplied real
bound of zero is valid. Growth tolerance must be finite/positive and the iteration
limit positive. The normal horizontal mode still needs no gamma_max.

## Algorithm and interpretation

For the supported pole-free upper-half-plane response, an accepted rectangular
winding is its enclosed zero count, including multiplicity. The full user-defined
rectangle is evaluated first by unmodified A1. If that calculation fails, there
is no accepted initial count or growth bracket. If its accepted winding is zero,
there are no roots above the selected lower boundary in this domain; the search
is skipped and no estimate is printed. Negative counts are rejected, never
converted by an absolute value.

For a positive count, initialize `[L,H]=[gamma_min,gamma_max]`. The upper endpoint
has zero roots above it **within the original search domain**; it is not an
independently established bound on all physical modes. At each step:

1. Set `M=0.5*L+0.5*H`.
2. Construct and completely validate the closed rectangle
   `[omega_min,omega_max] × [M,H]`, using the **current** H.
3. Positive accepted count sets L=M; zero sets H=M.
4. Stop successfully only when `H-L < tolerance`.

Because every accepted removed upper strip has zero roots, the remaining strip's
count is also the count above M inside the original domain. No part of a
previous positive strip is incorrectly discarded. The implementation checks that
accepted counts remain nonnegative and cannot exceed the retained lower count.
Both degenerate rates (2→0) and distinct rates (2→1→0) are tested.

A failed strip does not advance either endpoint. The history records its failure
and unchanged bracket. An iteration or coordinate-resolution limit similarly
retains the last bracket but emits no successful scalar estimate. The initial
upper boundary itself is validated by the initial contour; a zero-area rectangle
at that boundary is not evaluated.

All estimates refer only to roots between the user's real-frequency bounds and
below the original upper-growth boundary, above the selected lower threshold.
They do not exclude roots elsewhere. Finite sampling remains numerical convergence
evidence rather than a proof against arbitrary missed loops or unsampled boundary
zeros; A3 does not weaken or replace A1's safeguards.

## Cache and MPI contract

The rank-0 cache lives for exactly one `run_nyquist` operation and is explicitly
cleared at entry. It stores complex frequency and **raw** determinant, using an
open-addressed hash table with exact Fortran complex equality. Adjacent
representable frequencies remain distinct; signed zero is canonicalized for
hashing to agree with numerical equality. There is no interpolation, approximate
matching, eviction, or phase contribution stored in the cache. Storage grows with
unique successful physical evaluations, at load factor no greater than one half.

On a hit the callback retrieves raw D and then applies the selected normalization,
if any, with the existing A2 log-polar routine. On a miss rank 0 broadcasts the
command/frequency and every rank calls the existing collective `disp`. Invalid
raw values are rejected rather than stored or hidden by normalization. Raw storage
prevents mixing normalized and unnormalized conventions. The normalized growth
search and every independently reevaluated normalized strip give the same bracket
and counts as the raw calculation.

Workers remain in the same request loop across the initial contour and every
strip. Hits perform no MPI call; misses keep the existing determinant collectives
unchanged. On success or recoverable numerical failure, rank 0 sends one stop
command, all ranks reduce their actual call counts/timing, receive operation
status, and follow the existing coordinated finalization. Every rank executed
5005 actual calls in the reference search; logical requests were made only by the
rank-0 sampler. No frequency or wavevector parallelism was added.

A3 preserves each complete `nyq_evaluate` invocation, including stricter passes,
the independent seed, curvature/magnitude criteria, depth limits and integer
consistency checks. `nyq_sampling%max_evaluations` is still the logical budget
across all passes **of each contour**, including cache hits. It is not a new
search-wide physical-call budget. The iteration cap limits the number of extra
contours. Caching cannot turn an exhausted logical budget into acceptance.

### Edge reuse and remaining redundancy

Unchanged top-edge points, vertical segments, and a previous bottom becoming a
later top can all supply exact cache hits. Their orientation is irrelevant to
raw-value reuse; A1 still accumulates phase in the current CCW order. Completed
edge phases are **not** spliced. Avoiding substantial A1 API changes preserves
its established whole-contour convergence procedure.

The initial rectangle makes 630 logical requests but only 490 physical calls;
140 hits save repeated pass work. A copied-adapter failure after the first strip
begins provides a direct cross-contour communication check: the first request
`(-1,0.105)` reuses the initial left-edge midpoint, then the next new frequency
fails after an actual collective call. Totals become 632 logical requests,
491 physical calls and 141 hits. Both debug and optimized artifacts verify that
one additional logical edge request generated **no** collective call.

Across the complete search there are 3685 hits. However, A1 still reconstructs
and validates every contour and revisits cached points logically. The shrinking
top can lie close to the roots and require more refinement than the original
fixed top. Hence uncached strips need 8690 requests versus 7658 for fixed-top full
rectangles: 1032 additional logical requests and 516 extra local refinements.
Persistent value caching yields a physical-call saving despite this redundancy.
Floating-point coordinates that differ even slightly are intentionally not merged.

## Physical validation

Use the A2 equilibrium: `kperp=0.001`, `kpar=0.3`, `vA/c=1e-4`, both parallel betas
3, equal densities, zero drift, physical electron/proton mass ratio, isotropic
electrons, proton anisotropy 2 (unstable) or 1 (control). Both species use the
existing analytic backend. The ordinary solver reference remains
`z≈±0.39666+0.085986i`; ordinary unstable/stable root and residual checks were
rerun through the complete A2 suite. Frequency printing limits still apply;
no extra precision is inferred from that output.

All successful search rows below have initial real bounds `[-1,1]`, top 0.2 and
tolerance 1e-4 unless indicated. Counts were checked before establishing these
new fixtures.

| Case | Initial winding | Steps | Final bracket / outcome | Logical requests | Physical calls |
|---|---:|---:|---|---:|---:|
| Unstable, lower=0.01, raw | 2 | 11 | `[0.0859814453125,0.08607421875]` | 8690 | 5005 |
| Omitted lower (default 1e-4) | 2 | 11 | `[0.085896923828125,0.08599453125]` | 9108 | 4853 |
| Normalized rectangles | 2 | 11 | Same bracket as raw | 10282 | 5737 |
| Initial sample output enabled | 2 | 11 | Same bracket/count/work as raw | 8690 | 5005 |
| Stable isotropic, lower=0.01 | 0 | 0 | No search/bracket/estimate | 630 | 494 |
| Raised lower=0.1 | 0 | 0 | No search/bracket/estimate | 582 | 458 |

The omitted-lower width is `9.7607421875e-5`, with midpoint
`0.0859457275390625`; both successful brackets contain the known growth near
0.085986. All accepted physical strip windings are either two or zero. There is
no root-count doubling imposed by symmetry.

Reference raw-search history:

| Step | Trial lower | Strip upper | W | Updated L | Updated H | Logical requests |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.105 | 0.2 | 0 | 0.01 | 0.105 | 572 |
| 2 | 0.0575 | 0.105 | 2 | 0.0575 | 0.105 | 618 |
| 3 | 0.08125 | 0.105 | 2 | 0.08125 | 0.105 | 648 |
| 4 | 0.093125 | 0.105 | 0 | 0.08125 | 0.093125 | 644 |
| 5 | 0.0871875 | 0.093125 | 0 | 0.08125 | 0.0871875 | 702 |
| 6 | 0.08421875 | 0.0871875 | 2 | 0.08421875 | 0.0871875 | 728 |
| 7 | 0.085703125 | 0.0871875 | 2 | 0.085703125 | 0.0871875 | 778 |
| 8 | 0.0864453125 | 0.0871875 | 0 | 0.085703125 | 0.0864453125 | 772 |
| 9 | 0.08607421875 | 0.0864453125 | 0 | 0.085703125 | 0.08607421875 | 828 |
| 10 | 0.085888671875 | 0.08607421875 | 2 | 0.085888671875 | 0.08607421875 | 852 |
| 11 | 0.0859814453125 | 0.08607421875 | 2 | 0.0859814453125 | 0.08607421875 | 918 |

Each of the 11 raw strips and 11 normalized strips was separately evaluated as
an **ordinary uncached closed rectangle** with identical bounds and normalization.
All accepted integers and logical A1 request counts agree, and phases agree within
1e-10 rad. The uncached strip search and the fixed-original-top baseline also
produce exactly the same bracket, trials and count decisions.

## Failure and mathematical tests

| Test | Result |
|---|---|
| Initial budget 10 | Growth status 2, A1 status 4; no initial count/bracket/estimate. |
| Per-contour budget 650 | Initial W=2 succeeds; attempt 5 exhausts A1 budget after four updates. Growth status 3/A1 status 4; retains `[0.08125,0.093125]`, no estimate. |
| Near-root strip, depth 8, tolerance 1e-14 | Initial count succeeds; attempt 7 hits A1 depth limit after six updates. Retains `[0.08421875,0.0871875]`, growth status 3/A1 status 5, no estimate. |
| Iteration cap 1 | Growth status 4; retains `[0.01,0.105]`, no successful estimate. |
| Quiet NaN injected after initial contour | Initial W=2 retained; first strip fails with A1 status 2; no bracket update or stable/success label; all ranks terminate together. |
| Zero/negative/NaN tolerance, zero iteration limit | Clear input error, communicator-wide abort before determinant sampling. |
| Maximum growth with horizontal geometry | Clear incompatible-geometry error. |
| Missing upper or real-frequency bound | Rejected; no silent physical-domain defaults. |

The returned-NaN test modifies only a copied adapter, after the real collective
`disp` call. It and all recoverable failures reach common finalization and return
nonzero operation status without MPI deadlock. A valid initial zero count returns
normal completion, explicitly distinct from every numerical failure.

The standalone Fortran module tests passed **10056 checks** in a full-check/trap
`-O0` build and an `-O3` build. Coverage includes analytic equal/distinct growth
rates, negative winding rejection, an exact sampled zero on the first strip,
invalid tolerance, initial/strip failures, retained brackets, cache-on/off phase
and bracket equivalence, 5000-key rehashing and retrieval, adjacent representable
frequency misses and signed-zero matching. Explicit cross-strip-hit checks also
passed against both completed physical-suite artifact sets.

## Benchmarks

Four ranks, same final optimized build objects, raw determinant, samples disabled;
serial runs on the same host. The comparison variants are copied validation
sources: one disables caching, the other also keeps each trial top at the original
0.2. Production always uses smaller strips and has no added public cache/baseline
selector.

| Method | Steps | Logical | Physical | Hits | Hit fraction | Collective determinant time (s) | MPI calculation time (s) | Passes | Refinements | Max depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cached closed strips | 11 | 8690 | 5005 | 3685 | 42.405% | 0.198736 | 0.209968 | 36 | 1657 | 16 |
| Uncached same strips | 11 | 8690 | 8690 | 0 | 0 | 0.273574 | 0.282484 | 36 | 1657 | 16 |
| Uncached full-original-top rectangles | 11 | 7658 | 7658 | 0 | 0 | 0.200831 | 0.206498 | 36 | 1141 | 16 |

Times use `MPI_Wtime`. Determinant time sums actual miss broadcasts plus collective
`disp` calls on the controller. Calculation time includes the initial contour,
all bisection/sampling/cache work, and worker termination, taking the maximum rank
duration. It excludes equilibrium preparation and process launch. Driver elapsed
before summary writing was 0.210241, 0.282743 and 0.206683 seconds respectively,
including compact history export. Approximate remaining calculation time outside
physical evaluations was 0.011232, 0.008911 and 0.005667 seconds.

The established gain is **42.4% fewer physical calls than uncached strips** and
**34.6% fewer than uncached full rectangles**. The additional strips use 4515
physical calls beyond the initial cached 490, averaging about 410 per step.
Cached strips were faster than uncached strips in this measurement, but were
similar to the full-rectangle baseline in wall time. Cost per frequency and host
scheduling vary; these short timings do not establish a general full-baseline
speedup. No A1 convergence test was sacrificed to improve the numbers. Numerical-VDF
performance must not be extrapolated from these warm analytic timings.

## Builds, regressions and provenance

All builds/outputs were isolated under `/tmp`; existing scientific results and
unrelated working-tree modifications were preserved. GNU Fortran 11.4.0, Open MPI
4.1.2, four ranks, `--oversubscribe`, existing BLAS/LAPACK, one library/OpenMP thread.

Actually run:

- Complete `validation/run_nyquist.sh` in optimized and debug configurations.
  This also runs the complete unchanged A1 suite: 306 checks per debug/optimized
  A1 build and identical printed mathematical results.
- Complete new `validation/run_nyquist_growth.sh` in optimized and debug builds,
  plus the final standalone analytic additions and cross-strip-hit checks.
- Legacy root, stable-root, map, consecutive-scan and plane-scan scientific outputs:
  nine files byte-identical to pre-A3. Comparison passed for the optimized build,
  and separately against a pre-A3 executable with matching debug flags.
- Clean configure and `make -j4 -C src` with the existing production `-O4 ...`
  flags; that executable passed the four-rank growth bracket/history/accounting test.
- Syntax, free-form line length, whitespace and unchanged-source/hash checks.

Optimized flags: `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`.
Debug: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow`.
Standalone math additionally checked `-std=f2008`, warnings, and
`-ffpe-trap=invalid,zero,overflow`; optimized math used `-O3`.

Identified executable SHA-256:

- A3 optimized: `16931d8f5666bdc791cb83de365350974b56d04ad273cdc3982bd84eddce0133`.
- A3 debug: `16f4c5973a9f18ae0469c5d09820b90abc930713d11e366d8a510c0e07075deb`.
- Pre-A3 optimized: `6bd4de8c9b4bc8089a9d5b273f7d0aba28f7282b1e13082518311b37376dcb61`.
- Pre-A3 debug: `1eafeaf34545ca4b434d33eaad2ba28cc781b44954a317cc5e240c04965ddbdc`.

Per-source manifests are in each build's `provenance.txt`; copied baseline
variants have their own source/executable hashes. Complete A3 artifacts:
`/tmp/alps-nyquist-a3.k9SL4Ml3` (optimized),
`/tmp/alps-nyquist-a3.Q5NOPg15` (debug). A2 reruns:
`/tmp/alps-nyquist-a2.TPXJ9Snu` (optimized),
`/tmp/alps-nyquist-a2.NWhCbG2M` (debug). Temporary storage is not a durable archive;
the new fixture inputs, this report and test drivers retain the essential evidence.
The pre-A3 snapshot was taken before edits from the existing workspace, not by
resetting unrelated work. No broad numerical-VDF physical validation is claimed.

Reproduce:

```sh
bash validation/run_nyquist.sh
bash validation/run_nyquist_growth.sh
FFLAGS='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow' \
  bash validation/run_nyquist.sh
FFLAGS='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=overflow' \
  bash validation/run_nyquist_growth.sh
```

For an external legacy comparison, set `ALPS_REFERENCE_EXE` to the identified
pre-change executable. New examples are in `tests/data/nyquist_growth/`.

## Output, limitations and subsequent scans

The existing summary and optional ASCII samples describe the **initial** contour.
New fields distinguish logical requests, actual collective calls, hits/misses,
initial winding, search status, bracket validity, completed updates versus attempted
strips, width, midpoint or `NONE`, aggregate sampling diagnostics and timing.
`growth_status` is the overall search result; the existing A1 `status` remains the
initial-contour status. Valid zero count has growth status 1; numerical/iteration
failures have separate statuses and never produce a successful estimate.
`solution/<runname>.nyquist.growth` records compact iteration history. No per-strip
large sample files are generated; sample output opt-in retains the initial path.

The primary remaining costs are repeated A1 pass work and complete phase-validation
of each strip. No direct edge-phase reuse was attempted because it would require
more invasive convergence/API changes. Reuse of certified edge contributions,
less repeated workspace/sampling work and canonical geometry arithmetic could be
investigated separately, preserving independent evidence against missed loops and
exact-value semantics. They are not required for the validated cached search.

Too-close boundary roots, limited logical budget/depth, and finite-precision
determinants can still make the search unresolved. Last valid brackets survive
these recoverable failures. A3 neither repairs the A2 high-frequency conditioning
problem nor propagates the pre-existing WOFZ overflow flag. A raw determinant
failure inside the backend can still take its existing fatal path. The reported
bracket is numerical evidence under the analytic backend assumptions, not a
certificate excluding roots beyond the original domain or all missed phase loops.

Future option 4/5 loops may reuse the current wavevector input structures and call
this driver after each properly synchronized wavevector update. Cache invalidation
at every operation is essential: frequencies alone are valid keys only while all
plasma/wavevector context stays fixed. The present implementation explicitly resets
that cache, preserves harmonic assignment and requires no frequency parallelism.
Future scans still need wavevector-loop dispatch, Bessel/harmonic rebuild handling
and indexed/aggregate output filenames. None of that scan functionality was added
in A3; options 4 and 5 remain unsupported.
