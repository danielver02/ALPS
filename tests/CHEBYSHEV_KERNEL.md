# Standalone Chebyshev kernel

The numerical kernel is in `alps_analyt`; production fitting and evaluation
continue to use their existing routines. Run its independent Fortran tests with:

```bash
bash tests/run_chebyshev_kernel.sh --matrix
```

This builds from current source in temporary directories and tests `-O0` and
`-O2`, each with invalid/zero/overflow traps disabled and enabled.

## Interfaces

```fortran
call cheb_fit(u, samples, maximum_degree, model, status, rcond, representation)
call cheb_evaluate(model, u, value, status, derivative, value_status, derivative_status)
call cheb_evaluate_polynomial(model, u, polynomial, status, polynomial_derivative)
```

`rcond` and `representation` are optional; use keywords to omit `rcond`.
`CHEB_DIRECT` is the default, preserving existing callers. `CHEB_LOG10` fits
`log10(samples)` with the same SVD, adaptive-degree, coefficient, and diagnostic
machinery. `CHEB_SUCCESS` and `CHEB_REDUCED` indicate valid fits for either
representation; `model%representation` distinguishes them. Coefficients have
full-weight `c(0)`. Stored derivative coefficients differentiate the raw series
with respect to normalized coordinate x; returned derivatives are with respect
to physical coordinate u.

Both evaluators accept real or complex u and optional derivatives. The default
evaluator returns physical values (`10**P` for LOG10), whereas
`cheb_evaluate_polynomial` always returns P. Derivatives are computed only when
requested. Neither evaluator imposes physical support or a fitting-interval
cutoff. Invalid representation metadata returns `CHEB_INVALID_REPRESENTATION`.

## Sample selection and diagnostics

LOG10 validates all coordinates and values before transforming positive samples,
including subnormals. Exact-zero samples strictly below/above the positive
coordinate interval are excluded. Unsorted and repeated coordinates are allowed.
Negative values return `CHEB_NEGATIVE_SAMPLES`; any zero within the interval,
including one coincident with a positive sample, returns `CHEB_INTERIOR_ZERO`.
An empty or all-zero logarithmic row returns `CHEB_NO_POSITIVE_SAMPLES`.
There are no artificial positive floors.

The model records `original_sample_count`, `retained_sample_count`,
`lower_tail_zeros`, `upper_tail_zeros`, `retained_min`, and `retained_max`.
Existing `sample_count` is the number supplied to the solver after selection.
Degree/rank selection uses that retained count. Residuals for LOG10 measure
errors in log10 space. The existing independent per-diagnostic availability
statuses remain authoritative, and unavailable diagnostics do not invalidate fits.

## Reconstruction and failure semantics

Complex reconstruction uses `exp(real(P)*ln(10))` and the sine/cosine of
`imag(P)*ln(10)`. Exponentiation uses a small binary remainder followed by
checked power-of-two scaling. Derivative reconstruction combines logarithmic
scale factors before exponentiation; it does not first multiply potentially
overflowing `ln(10)*f*P_u` factors. Value and derivative ranges are assessed
independently.

* `CHEB_SUCCESS`: every requested result is finite, with no detected underflow.
* `CHEB_LOG_UNDERFLOW`: a nonzero magnitude/component lies below the least
  subnormal. Affected outputs become zero; other components/results survive.
  Representable nonzero subnormals are ordinary successful results.
* `CHEB_LOG_OVERFLOW`: the distribution magnitude exceeds the supported range.
* `CHEB_LOG_VALUE_FAILURE`: reconstruction produced an invalid result.
* `CHEB_LOG_DERIVATIVE_FAILURE`: a requested reconstructed derivative failed;
  its optional per-output status distinguishes overflow from other failure.
* `CHEB_PHASE_FAILURE`: the phase cannot be reconstructed to the documented
  accuracy bound.

Optional `value_status` and `derivative_status` distinguish independent range
conditions. For example, a zero/UNDERFLOW value can accompany a finite/SUCCESS
derivative. An unrequested derivative has `CHEB_NOT_REQUESTED` status. Underflow
is a nonfatal condition, rather than ordinary success or a failed fit. Hard
failures return quiet NaNs for all requested outputs, even when one intermediate
result was valid. They never invalidate the fitted model.

## Numerical limits

No routine changes the caller's floating-point environment. Finite inputs are
supported; quiet NaNs/infinities are rejected by status. Signaling NaNs may raise
IEEE invalid before a status is returned, especially with invalid traps enabled.
Underflow trapping itself is outside the tested configuration.

The phase must satisfy `abs(imag(P)*ln(10)) <= 1/sqrt(epsilon)` radians, checked
before multiplication. This bounds the absolute phase uncertainty introduced by
working-precision multiplication to O(sqrt(epsilon)); library trigonometric
functions perform argument reduction. Larger phases return an explicit failure,
regardless of physical support. Accuracy near that limit is weaker than for
small phases.

Rounded logarithms can conservatively reject reconstruction within a few ulps
of the overflow boundary; values are never clamped to HUGE. Existing Clenshaw
intermediate and raw physical-derivative range limitations still apply, even
when a later mathematical cancellation/scaling could produce a finite result.

## Relativistic row collection (Task 3B)

```fortran
call build_rel_cheb_rows(gamma, parallel, distribution, species, degrees, status, &
                        selected, diagnostic_unit)
call get_rel_cheb_row(is_rel, igamma, model, info, status)
call rel_cheb_collection_info(nspecies, last_gamma, fit_calls, coefficient_bytes)
call release_rel_cheb_rows
```

The builder is an explicit, test-accessible initialization entry in `alps_analyt`.
Task 3C now calls the conditional `initialize_rel_cheb` hook after
`derivative_f0` returns in `ALPS.f90`, before continuation setup: both fresh
transformation and cache loading have then completed. The relativistic
`AC_method=2` input guard is unchanged; normal production input cannot select
this representation yet. A dedicated MPI harness explicitly sets internal state. Do not call the
builder for every relativistic species unconditionally. Construction must be
explicitly selected; this stage's tests exercise it independently of input-method
selection. The builder has no MPI logic or polynomial cache. Task 3C distributes
completed coefficients as described below.

Inputs are completed, normalized tables with bounds
`(1:nspec_rel,0:ngamma,0:npparbar)`. `species` maps packed slots to distinct,
increasing original species indices. Supply the existing setting as
`degrees=poly_order(species)`; no new input parameter or degree heuristic exists.
The optional `selected` mask retains packed indexing but leaves unselected slots
unfitted. Calling the explicit builder without this mask selects all supplied
slots. Original source-grid arrays, target coordinates, values, derivatives,
normalization, physical support, and integration are not changed.

Only exact `-1` distribution sentinels are removed. Other negative/nonfinite
values and exact zeros pass to kernel validation/filtering. Gamma coordinates
must be finite, at least one, and constant across each row. No new coordinates,
interpolation, curved-boundary extension, or row renormalization is performed.
Models fit the normalized distribution itself, not Gamma times the distribution.
All fitted models use `CHEB_LOG10`; analytic derivatives are with respect to the
supplied dimensionless `u=pparbar_rel`, without another mass or vA factor.

Storage is private, with just packed-species and zero-based Gamma-row dimensions.
The accessor deep-copies the model and `rel_cheb_row_info`; changing the copy
cannot mutate stored coefficients. It never fits. `coefficient_bytes` is an
`integer(int64)` counting value/derivative coefficient payload, excluding model
metadata/descriptors and temporary fitting workspace. `fit_calls` counts calls
to `cheb_fit` in the current construction, including rows rejected by the kernel;
each model separately records its SVD `fit_attempts`. Rebuilding/releasing frees
the entire old collection, including nested arrays, and resets the count. Invalid
collection arguments also release old storage to prevent stale-model access.

`info` records original/packed species, Gamma index/value, input/unmasked/masked/
positive sample counts, row state, and kernel status:

| State | Meaning |
|---|---|
| `REL_CHEB_FITTED` | Valid `CHEB_SUCCESS` model |
| `REL_CHEB_REDUCED` | Valid `CHEB_REDUCED` model; degree/rank reasons remain in model |
| `REL_CHEB_NO_GEOMETRY` | All samples masked; `CHEB_NO_SAMPLES` |
| `REL_CHEB_ZERO` | All unmasked samples zero; `CHEB_NO_POSITIVE_SAMPLES` |
| `REL_CHEB_INVALID_SAMPLES` | Coordinate/value/nonpositive-sample validation failure |
| `REL_CHEB_FIT_FAILURE` | Solver, degree, allocation, or coefficient failure |
| `REL_CHEB_NOT_SELECTED` | Explicitly unselected; no model |
| `REL_CHEB_NOT_BUILT` | Collection absent or accessor indices invalid |

Unavailable rows never acquire invented coefficients. Collection `CHEB_SUCCESS`
means processing completed without unexpected failures; empty/zero rows still
require callers to inspect their state. Otherwise the first unexpected kernel
error is returned, every row is still processed, and each unexpected failure is
reported with species IDs, Gamma index, kernel code, and sample counts. Optional
`diagnostic_unit` redirects those messages from `error_unit`. Structural shape/
mapping errors return `CHEB_INVALID_VALUES` before any row can be constructed.

Numerical zero or empty sampled geometry does not establish zero analytic
continuation. The accessor does not apply physical support restrictions, and a
neighboring Gamma-row evaluation may lie outside its sampled interval. No new
fallback or production interpretation of unavailable rows is introduced here.
An accepted model means the numerical kernel succeeded; residual diagnostics
and independent reconstruction checks still determine approximation quality.
A singleton's constant model has a zero analytic derivative and cannot recover
the nonzero derivative of a general underlying distribution from one sample.

The runner adds independent analytic row tests (isotropic, anisotropic, boosted
thermal), row classifications/lifecycle, degree-30 setup benchmarks at 100x200
and 500x500 intervals, and a separate current-transformation/cache-hit fixture.
All generated distributions, diagnostics, and cache entries are in scratch.
The cache fixture uses the existing Linux/OpenBLAS identity machinery and fails
explicitly if that existing cache backend is unavailable. It does not use or
recapture historical numerical baselines. Its coarse spline anisotropy has a
separate 1% row-reconstruction and 3% peak-scaled derivative budget; exact
analytic-row checks are much tighter and interpolation error is reported apart.

## MPI and continuation adapter (Task 3C, still guarded)

`initialize_rel_cheb(status)` is called on rank 0 in `ALPS.f90` immediately
after `derivative_f0`, covering transformation and cache loading alike. It
selects only `relativistic .and. ACmethod==2 .and. .not.usebM`, preserves all
packed relativistic slots, and calls the row builder once with
`poly_order(species)`. With no selection it allocates/fits nothing. The normal
input guard still rejects relativistic method 2 before this point. There is no
new input option or exposed production activation path.

`pass_distribution` calls `pass_rel_cheb_models(status)` once at initialization.
The latter is in `alps_com`, keeping MPI out of the numerical kernel. Its
collection-specific export/import helpers in `alps_analyt` handle three arrays:

* Integer metadata: original species, row state, kernel status, requested degree,
  effective degree, numerical rank, solver sample count, and representation.
* Real mapping metadata: coordinate center, scale, and row Gamma.
* Real coefficients: padded value and derivative blocks, sized to the largest
  accepted degree, with one column per row in Gamma/packed-species order.

No derived-type descriptors, workspaces, or residual/SVD diagnostics are sent.
On receivers, coefficients get their actual zero-based effective-degree bounds.
Fit diagnostics are unavailable (`fit_attempts=0`, diagnostic availability flags
UNAVAILABLE, residuals/SVD values NaN); row sample-selection counts are -1.
Reduction-reason flags and other fitting-only counts are not meaningful there.
Root retains the complete fit diagnostics. Required evaluation metadata and all
coefficient values are preserved exactly.

A successful nonempty transfer uses six collectives: header broadcast,
allocation-status agreement, three buffer broadcasts, and import-status
agreement. Root errors are communicated in the header before any payload. A
failure on any rank releases the collection everywhere; the production caller
aborts consistently rather than continue with partially initialized data.
An inactive collection uses only the header and no payload transfer. No
frequency-, wavevector-, harmonic-, or row-dependent communication occurs.

```fortran
call eval_rel_cheb(is_rel, igamma, complex_u, value, status, derivative)
value = eval_fit(original_species, igamma, complex_u, status=status)
```

The direct accessor evaluates the private stored model in place: no model copy,
new allocation, fitting, MPI, or transformed-table access. Derivatives are
optional. The copy-returning diagnostic accessor is unchanged.
`eval_fit` has an optional trailing `status` argument and a new relativistic
method-2 branch using the existing packed-species convention. Its inputs are
already Gamma row and normalized complex u. It applies no additional coordinate
transformation or physical-support restriction. Nonrelativistic evaluation
arithmetic remains unchanged. Relativistic method-2 rows are skipped by the
legacy original-grid fitting branch because initialization has already built
them once.

With an explicit status output, continuation callers receive kernel statuses,
including unavailable rows and underflow. Without one, any non-SUCCESS result
uses `alps_continuation_error` to report species/row/status/u and abort the MPI
communicator. In particular, underflow is not silently accepted by legacy
continuation callers. The existing `landau_integrate_rel` source, neighboring-row
evaluations, Gamma/parallel finite differences, resonant numerator, and physical
support tests are unchanged. Analytic derivatives are available for independent
checks but are not substituted into production stencils.

For lifecycle checks, `rel_cheb_collection_info` optionally returns
`build_calls` and `evaluations` (integer/int64); imported workers have zero
builder/fit calls. `rel_cheb_distribution_info` returns transfer count, collective
count, and approximate buffer bytes. These counters are local and reset per
initialization, and reading them performs no communication. They are small
validation diagnostics, not a profiling subsystem.

```bash
bash tests/run_rel_chebyshev_mpi.sh
bash tests/run_chebyshev_kernel.sh --matrix
bash tests/run_input_validation.sh
```

The isolated MPI harness runs with 1/2/4 ranks, synthetic and actual transformed
tables, reduced/unavailable/unselected rows, serial-reference construction,
bitwise cross-rank coefficient/metadata comparisons, and local continuation
equivalence. It removes the transformed tables before repeated evaluation,
checks unchanged finite-difference expressions, and measures construction,
transfer, and 10,000 evaluations separately. Test-only MIN/MAX reductions compare
rank contents; they are outside the production transfer/evaluation counters and
timings. A separate two-rank negative case verifies the status-less failure path.
MPI tests default to invalid/overflow traps; the serial matrix separately tests
invalid/zero/overflow, since existing MPI initialization can be incompatible
with zero-division traps. All files remain in temporary directories; no
historical references are captured or modified.
