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
