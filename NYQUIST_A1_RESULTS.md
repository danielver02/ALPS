# ALPS Nyquist Stage A1 — implementation and validation

Date: 2026-11-20 (America/Phoenix). Stage A0 source baseline:
`5ee0d335a037744a5ce298c355b84fc316af5f9c`.

Stage A1 implements the standalone mathematical kernel in
[`src/ALPS_nyquist.f90`](src/ALPS_nyquist.f90). The analytic test driver is
[`validation/nyquist/test_nyquist.f90`](validation/nyquist/test_nyquist.f90), with
[build instructions](validation/nyquist/README.md) and a separate test runner.
No existing production source, build configuration, input handling, determinant,
susceptibility, root solver, or MPI decomposition was changed. Pre-existing
working-tree modifications were preserved. No physical test cases were run.

## Validation outcome

GNU Fortran 11.4.0 completed **306 checks with zero failures in each build**:

- Debug: `-std=f2008 -O0 -g -Wall -Wextra -Wno-compare-reals -fcheck=all
  -fbacktrace -ffpe-trap=invalid,zero,overflow`.
- Optimized: `-std=f2008 -O3 -Wall -Wextra -Wno-compare-reals`.

Both runs produced identical printed results, including phase values to twelve
decimal places and evaluation/refinement counts. Each build separately checks
analytic finite-path phases with absolute tolerance `1e-10` radians and elementary
sample accumulation with `1e-13` radians. This is a standalone regression test;
production ALPS was not rebuilt or executed. Reproduce with:

```sh
bash validation/nyquist/run_tests.sh
```

## Public interfaces and separation of responsibilities

All interfaces use `nyq_dp=kind(1.d0)`, matching ALPS double-precision real/complex
arithmetic without importing `alps_var` or MPI.

| Interface/type | Responsibility |
|---|---|
| `nyq_accumulate(values,phase,status[,cumulative,max_increment])` | Signed phase sum for ordered nonzero finite complex values; no sampling or geometry. |
| `nyq_path` | Ordered polyline vertices plus `closed`; each segment has affine parameter `t` in `[0,1]`. Arbitrary prescribed polylines are permitted. |
| `nyq_rectangle(xmin,xmax,ymin,ymax,path,status)` | Finite ordered bounds; bottom left → bottom right → top right → top left → closure, counterclockwise. No restriction to positive imaginary coordinates. |
| `nyq_horizontal(xmin,xmax,y,path,status)` | Open left-to-right finite horizontal segment. |
| `nyq_evaluate(evaluate,path,options,result[,save_samples])` | Serial adaptive sampling and independent convergence passes on the supplied geometry. |
| `nyq_evaluator` | Abstract subroutine interface: double complex `z` in, double complex `value` out, integer `status` out (zero means success). |
| `nyq_options` | Initial segment density, local criteria, depth/evaluation/pass caps, phase/integer tolerances. |
| `nyq_result` | Phase, convergence/acceptance flags, winding, statuses, diagnostics, endpoints, optional ordered samples. |
| `nyq_write_ascii(filename,result,iostat)` | Explicit optional ASCII export; never called by the sampler. |

Reverse a contour by reversing its vertex order. `nyq_accumulate` closes nothing
implicitly: include the repeated first value to accumulate a closing segment.
The path sampler supplies that repeated value automatically for a closed path,
reusing its first evaluation. Each corner is evaluated once within a pass.
Separate global passes deliberately reevaluate the function; evaluations are
counted across all passes, including probes subsequently discarded by refinement.

The evaluator is neither pure nor concurrent. It can retain context or call
collectives synchronously. A2 can run this sampler on rank 0, broadcasting every
requested frequency to workers through an adapter and collectively calling the
existing `disp(om)`. Worker coordination, termination, initialization, and ALPS
normalization remain A2 responsibilities. The numerical algorithms need no MPI
changes for that integration.

`result%phase` is signed radians; `raw_winding=phase/(2*pi)` is diagnostic.
`winding` is meaningful **only when `integer_accepted` is true**. Open results
never set that flag, even when their phase happens to equal an integer multiple
of `2*pi`. `converged` for an open path means only sampling convergence on the
prescribed finite interval. `endpoints(1:2)` expose evaluated complex function
values at the ends of the accepted path.

Diagnostics include cumulative evaluation count, attempted pass count, total
refinements, deepest recursive level, maximum remaining adjacent phase increment,
and `pass_phases(1:passes)` (completed pass totals; unfinished entries remain zero).
With `save_samples`, `z`, `values`, and `cumulative_phase` contain the final ordered
samples, including midpoint probes and the closing value for closed paths.
On failure these arrays, phase, and endpoints describe only the successfully
completed prefix of the last attempted pass; they must not be interpreted as an
accepted full-path calculation. Empty prefixes have zero phase/endpoints.

## Phase arithmetic and convergence contract

Each value is divided first by its largest absolute real/imaginary component,
then by its safe scaled modulus. The increment is
`atan2(Im(u_next*conjg(u_prev)), Re(u_next*conjg(u_prev)))`.
This avoids products of the original magnitudes and even handles finite complex
components whose Euclidean modulus overflows. Compensated summation preserves
the signed accumulated phase and its unwrapped cumulative history.

There is no absolute magnitude cutoff. Exact zero, nonfinite components, and
nonzero evaluator status reject an evaluation. Subnormal nonzero values can be
used if their phase is still reliable; the kernel cannot assess evaluator
conditioning or recover information already lost before it receives a value.

Default controls:

| Control | Default |
|---|---:|
| Initial intervals per geometry segment | 8 |
| Maximum recursive depth per initial interval | 24 |
| Total evaluator call budget across all passes | 200000 |
| Maximum passes | 6 |
| Local half-interval phase limit | 0.7 rad (0.49 after pass 1) |
| Relative midpoint curvature limit | 0.15 |
| Adjacent midpoint/endpoint log-magnitude change limit | 1.5 |
| Minimum local magnitude ratio | 0.05 |
| Absolute pass phase tolerance | `1e-8` rad |
| Closed integer residual tolerance | `1e-8` winding units |

Every initial interval gets a midpoint probe. Refinement bisects intervals when
half-step phase changes are large, midpoint-resolved and direct phase changes
disagree, scaled complex midpoint curvature is large, or local magnitude changes
are excessive. Curvature uses a common component scale; magnitude comparisons
use logarithms. Accepted intervals retain their midpoint, so all phase probes
contribute to the reported phase. Smooth intervals get no additional local
subdivision beyond that probe.

Successive global passes double the initial density and tighten the phase
criterion after pass 1. Odd passes starting with pass 3 independently perturb
the grid interior using
`t_j = (j + 0.2718281828459*sin(pi*j/n))/n`, with exact segment endpoints.
These nonuniform points do not follow the previous dyadic subdivision sequence.
Acceptance requires agreement of the last **three** completed phase totals;
each such group includes an independently seeded pass. Closed paths must also
pass the integer residual check. Agreement is absolute in accumulated radians,
not agreement of rounded integers. Bounds, paths, and numerical options are
validated before evaluator calls.

This is numerical convergence evidence, **not rigorous exclusion of arbitrary
unsampled loops**. The independent seed is deterministic, not a probabilistic
certificate. A specially constructed function can alias any finite set of probes.
Users can increase density, tighten criteria, or vary geometry, but certification
requires additional analytic information about the supplied function.

## Closed-contour results

Let `i` denote the imaginary unit. Unless indicated otherwise, the contour is
`[-1,1] × [0.5,2.5]`. Expectations follow from analytic zeros minus poles with
multiplicity, and reverse sign for clockwise traversal.

| Function / contour | Expected = accepted W | Evaluations |
|---|---:|---:|
| `z-i`, root included | +1 | 448 |
| `z-i`, imaginary bounds `[1.5,2.5]`, excluded | 0 | 448 |
| `(z-i)(z-2i)`, both included | +2 | 448 |
| Same, imaginary bounds `[0.5,1.5]`, one included | +1 | 448 |
| `(z-i)**2`, repeated root | +2 | 448 |
| `1/(z-i)` | −1 | 448 |
| `(z-i)/(z-2i)`, both included | 0 | 448 |
| `z-i`, reversed vertex order | −1 | 448 |
| `(z-i)(z-i-1e-4)(z-i+1e-4)` | +3 | 452 |
| `z-i`, imaginary bounds `[-1,1.5]` | +1 | 448 |
| `1e280*(z-i)`, imaginary bounds `[0.5,1.5]` | +1 | 448 |
| `1e-280*(z-i)`, same contour | +1 | 448 |
| `exp(650*z)`, same contour, no zeros/poles | 0 | 24576 |
| `z-i`, lower edge `1-1e-4`, upper edge 2 | +1 | 560 |
| `z-i`, lower edge `1+1e-4`, upper edge 2 | 0 | 560 |
| `z-i`, lower edge 1, zero at sampled edge midpoint | Unresolved: bad value | — |
| `z-i`, lower edge `1-1e-12`, upper edge 2 | Unresolved: depth cap | — |

The smooth one-root case requires three passes and no local subdivision; the
final maximum phase increment is `0.06256` rad. The near-boundary resolved cases
require 56 local subdivisions across three passes and have final maximum
increment `0.41708` rad. The exponential spans about 564 decades of magnitude,
requires 12064 subdivisions, and leaves maximum increment `0.32584` rad. Its
large cost reflects both magnitude gradients and rapid mapped phase rotation.

A constant-function rectangle with 1 initial interval per edge verifies exact
corner/closure reuse: passes with 1, 2, and 4 intervals use `8*(1+2+4)=56`
evaluations, including every midpoint and no duplicate closing evaluation.

## Deliberate full-loop undersampling

On the open segment `[0,1]`, use `exp(16*pi*i*z)`. Its analytically known phase
change is `16*pi`: eight complete counterclockwise loops, with identical endpoint
values. Begin with one interval per segment:

| Pass | Grid | Computed phase / `2*pi` |
|---|---|---:|
| 1 | 1 interval, midpoint probe | approximately 0, incorrect |
| 2 | 2 intervals, midpoint probes | approximately 0, incorrect |
| 3 | 4 intervals, independently perturbed grid and adaptation | 8 |

A three-pass budget returns `nyq_no_agreement`, despite the first two apparently
consistent integers. Allowing six passes leads to acceptance after three resolved
passes agree. This demonstrates detection of this constructed failure, not a
universal loop-detection guarantee. Open-path integer acceptance remains false.

## Open horizontal paths: finite sampling versus closure

For `z=x+i*g`, define the exact continuous phase of the linear factor along
`[-R,R]` by

- `P(a) = -2*atan(R/(g-a))` for `g>a`;
- `P(a) = +2*atan(R/(a-g))` for `g<a`.

Thus the exact finite phase of `(z-i*a)/(z-i*b)` is `P(a)-P(b)`. Tests compare
against these formulas, never against a winding inferred from sampled values.
For `g=0.25`, `F1=(z-i)/(z+i)` has limit `2*pi`, while
`F2=(z-i)(z-2i)/((z+i)(z+2i))` has limit `4*pi`. Both tend to unity at infinity,
have their zeros above the horizontal path and poles below it, and have vanishing
upper closure contribution as `R` tends to infinity.

Every entry below is sampling-converged; the closure error is the analytic
asymptotic phase minus the finite phase, not a sampling error.

| R | F1 finite phase (rad) | F1 closure error | Calls | F2 finite phase (rad) | F2 closure error | Calls |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 4.448445335951 | 1.8347 | 115 | 7.605662671661 | 4.9607 | 119 |
| 10 | 5.884755622665 | 0.39843 | 135 | 11.378820712244 | 1.1875 | 163 |
| 100 | 6.243186890381 | 0.039998 | 175 | 12.446383361265 | 0.11999 | 219 |
| 1000 | 6.279185308763 | 0.0040000 | 215 | 12.554370627109 | 0.012000 | 251 |

The tests explicitly reject the inference that sampling convergence at `R=2`
means an adequate asymptotic interval. Even `R=1000` has finite closure errors
far above the sampling comparison tolerance. All open results retain their
unrounded phase and function endpoints.

Additional finite phases agree with the analytic formulas to `1e-10` rad:

| Function / line | Asymptotic phase | Additional validation |
|---|---:|---|
| `1`, any finite horizontal interval | 0 | Exact zero finite phase. |
| `(z+i)/(z+2i)`, `g=0.25` | 0 | No zeros/poles above; repeated R=2,10,100,1000 approaches zero. |
| `(z-i)/(z-2i)`, `g=0.25` | 0 | Equal upper zero/pole counts; same R sequence approaches zero. |
| F2, `g=3` | 0 | Line above all zeros/poles; same R sequence approaches zero. |
| `F1**20`, `g=0.999`, R=1000 | `40*pi` | Finite phase `125.583706250098`, 2209 calls, over 60 decades of magnitude variation. Finite phase is not the asymptotic limit. |
| `exp(650*z)`, `g=0`, R=1 | No unity asymptote | Supplemental generic finite-path check: exact zero phase despite about 564 decades of magnitude variation; no root-count interpretation. |

Similar endpoint phases cannot establish resolved interior variation: the
constructed eight-loop example fails initial sampling despite equal endpoints.
The raw normalization example below likewise has nearly matching principal
endpoint phases but a continuous phase near `-4*pi`.

## Explicit sixth-power normalization and continuous correction

Normalization lives only in the test evaluator. The kernel never mentions `z**6`,
`vA`, plasma tensors, or high-frequency limits. Define `Dtest=z**6*F1`.

On `[-1,1] × [0.25,1.5]`, zero is outside, the zero of F1 at `i` is inside, and its
pole at `-i` is outside. Both F1 and Dtest therefore have accepted winding **+1**:
448 and 514 calls respectively. The extra six zeros at the origin are not enclosed.

On `[-R,R]+0.25*i`, the continuously accumulated phase of `z` is
`-2*atan(R/0.25)`. Consequently

`phase(Dtest)-phase(F1) = -12*atan(R/0.25) -> -6*pi`.

The test evaluates F1, Dtest, and z independently through the same sampler,
checks this exact finite difference to `1e-10` rad, and checks
`phase(Dtest)-6*phase(z) = phase(F1)` to the same tolerance.

| R | Direct F1 phase | Raw Dtest phase | Continuously corrected raw phase |
|---:|---:|---:|---:|
| 2 | 4.448445335951 | −12.908850651026 | 4.448445335951 |
| 10 | 5.884755622665 | −12.664862775447 | 5.884755622665 |
| 100 | 6.243186890381 | −12.576369093657 | 6.243186890381 |
| 1000 | 6.279185308763 | −12.567370612838 | 6.279185308763 |

The phase difference tends to **minus three windings**, as prescribed. Each path
avoids z=0. Principal endpoint argument subtraction is not used for the correction.
Neither direct normalization nor phase correction can repair inaccurate,
overflowing, or underflowing raw determinant evaluations.

## Failure statuses and robustness coverage

| Status | Meaning / demonstrated outcome |
|---|---|
| `nyq_ok=0` | Sampling convergence established; closed integer additionally accepted. |
| `nyq_bad_input=1` | Invalid bounds/path/options, rejected before evaluation. |
| `nyq_bad_value=2` | Exact sampled zero, NaN, or infinity. |
| `nyq_evaluator_failed=3` | Evaluator failure; test verifies original code 37 retained in `evaluator_status`. |
| `nyq_budget=4` | Total budget exhausted; test limits to exactly 10 calls. |
| `nyq_depth=5` | Required local subdivision exceeds depth cap, or midpoint cannot be represented distinctly. Depth-zero and extremely near-boundary tests exercise the cap. |
| `nyq_no_agreement=6` | Allowed passes lack three agreeing totals, or closed integer consistency cannot be accepted. Deliberately aliased loop case exercises this. |

All failure tests assert that integer acceptance is false. Ordered sample tests
verify positive real-coordinate progression on horizontal paths and consistency
between retained values, cumulative phase, and total phase. ASCII export is
checked for success and correct row count, then removed. Elementary accumulation
also tests both branch-crossing signs, `1e-300` samples, and complex values with
finite components of magnitude `0.8*huge(double)` whose direct modulus would
overflow. No arbitrary magnitude cutoff is introduced.

## Numerical limitations and A2 recommendation

Finite sampling cannot prove that an arbitrary function has no unresolved loops
or unsampled boundary zeros. Very close boundary zeros can remain unresolved,
and exact poles must be reported reliably by the evaluator rather than performing
an unchecked division. No derivative bounds, evaluator error estimates,
analyticity certification, or automatic pole inventory are supplied. Memory is
linear in the configured evaluation budget, and global verification has a real
cost even when local intervals are smooth. The tests establish scalar mathematical
behavior, not the accuracy or analytic properties of the ALPS determinant.

Proceed to A2 with this synchronous evaluator contract. Use raw `D_ALPS` for
strictly upper-half-plane rectangles as the finite-region reference. Prefer the
normalized determinant for horizontal production investigations, or apply the
independently resolved continuous six-power correction outside the kernel.
Investigate interval extension and physical determinant accuracy separately from
sampling convergence. Endpoint agreement alone must never justify a closure or
an instability count. The lower-frequency threshold and outer contour bounds
remain explicit physical choices for A2; A1 supplies no universal cutoff.
