# C2f.4a.1 — Focused weak-mode resolution diagnostic

**Completed on 2026-11-24: successful focused resolution diagnostic.** The
two finer-grid roots, tensor errors and determinant errors confirm systematic
near-quadratic convergence. The controlled derivative substitution removes
about 82.7% of the 96×192 growth bias while leaving a material fixed-grid
residual. Both requested lower contours accept winding three. The corrected
frozen entry completed start to finish; no production change or
integration-method adoption is performed.

## 1. Source and fixture provenance

The baseline is Git `99a087271928846cb8dbc49507d5a3a0f4d68d84` with the retained
prior working-tree changes. Current production sources match the qualified
C2f.4a sources. In particular:

| Source | SHA-256 |
|---|---|
| src/ALPS_fns.f90 | 0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488 |
| src/ALPS_fns_rel.f90 | c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6 |
| src/ALPS.f90 | acb9e48ab671964ad537c73651d413b1e6e34f2af43ab3174545de6a5109097b |

The entry reuses the C2f.4a source builder, compact tensor observer, C2f.3a
restricted extended Cauchy arithmetic, C2f.4a sign-bracket/root qualification,
and existing scale-aware wave-matrix reader. The coherent formulation,
fixed-cell/PWL kernel, weights, Bessel values, corrected int_ee, harmonic
stopping, MPI assignment, root stopping and A1 mathematics are unchanged.
Production C2b remains the default.

Both species are numerical, with beta_parallel,p=beta_parallel,e=3,
proton anisotropy 2, isotropic electrons, equal unit densities, zero drift,
electron mass `5.44617021487e-4`, vA/c=`1e-4` and kperp=kpar=.3. Normal ALPS
generation and differentiation use the original four-width support:
proton pperp_max=4*sqrt(6), |ppar_max|=4*sqrt(3), and both electron extents
4*sqrt(3*m_e). N denotes perpendicular intervals; parallel intervals are 2N.
No finite-grid normalization is applied.

The 64/96 roots and tensors are extracted from the immutable C2f.4a archives,
not rerun as matrices. The exact fixed analytic frequency is the archived
ordinary-solver return

```
z_A = 1.3706337381161255e-12 + 0.009424312472916148 i.
```

Its tiny real part is retained to make frequency matching exact; this is the
reference normally quoted as i*0.0094243124, not a claim of an exactly known
16-digit physical growth rate. The old root's stopping uncertainty is much
smaller than the grid/extrapolation differences reported here.

A preliminary harness attempt copied this analytic input before sourcing an
old emitter, which overwrote it with a nearly parallel input. The mismatch
was detected from the actual analytic tensor, not hidden by changing data.
That attempt was stopped after its current MPI job finalized. Preparation now
copies and byte-checks the canonical oblique input after emission; an actual
analytic tensor/wave comparison also rejects the retained bad control at the
existing component-scale bound. The corrected frozen entry reruns the focused
cases. Rejected setup evidence is under `preflight/`, separate from acceptance.
The preliminary parent EXIT trap's signal-status issue was also corrected;
completion requires both exit=0 and COMPLETE.

Reproduce with `bash validation/run_nyquist_c2f4a1.sh`, optionally setting a
fresh external `C2F4A1_SCRATCH`. The corrected frozen run uses
`/tmp/alps-c2f4a1-acceptance-final-20261124`; reproduction has no dependency on
that directory. GNU Fortran 11.4.0 / Open MPI 4.1.2 and existing BLAS/LAPACK are
used, with four ranks and one BLAS/OpenMP thread. Physical flags are the prior
`-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`. Read-only analysis is
also compiled with -O0, all checks and invalid/zero/overflow traps. No full
physical debug matrix, isotropic control, original contour, harmonic audit,
cache benchmark, trajectory matrix or analytic compatibility suite is repeated.
The runner retains its assertions; only its external wall limit grows from
900 to 3600 seconds for the finest contour, not a scientific tolerance.

## 2. Two new bracketed weak roots

Every new ordinary solve follows a measured real-D sign change on the positive
imaginary axis, using thirteen frequencies from gamma=.004 to .016. The
unchanged ordinary residual/positive-trajectory qualifications apply, and the
returned root must remain in the observed bracket with negligible real part.
All roots are far from the unresolved near-origin regime.

| Grid | Root real frequency | Growth rate | Sign bracket | Raw residual magnitude | Scaled residual magnitude |
|---|---:|---:|---|---:|---:|
| 128×256 | 8.83e-14 | .008665637513400 | [.008,.009] | 4.03e-16 | 5.37e-12 |
| 192×384 | 2.81e-13 | .009087488939534 | [.009,.010] | 6.54e-15 | 7.92e-11 |

The growth errors are about 8.05% and 3.57% below the analytic reference,
respectively, versus 14.3% at 96×192. Conditioning of the existing rotated,
longitudinally scaled wave matrix gives reciprocal Frobenius indicators
`2.09e-13` and `3.09e-12`, consistent with finite-frequency singularities.
Raw versus real128 determinants of the *same rounded wave matrix* differ by
`1.46e-19` and `1.25e-20`. This is not a real128 susceptibility calculation.

`resolution.dat` scales the raw returned determinant directly by z².
The existing `conditioning_*.dat` scaled column uses the independently
recomputed real128 determinant of the recorded double matrix; small relative
differences between these near-zero residuals are expected and retained.
Neither diagnostic changes production normalization or the A1 determinant.
Full six-component species tensors and complete 3×3 wave matrices at z_A are
in `fixed_components.dat`, `fixed_wave.dat` and the original tensor records.

## 3. Four-resolution convergence and Richardson analysis

| N | Measured growth | Error from analytic |
|---:|---:|---:|
| 64 | .006369754970337 | -.003054557502579 |
| 96 | .008075085805121 | -.001349226667795 |
| 128 | .008665637513400 | -.000758674959516 |
| 192 | .009087488939534 | -.000336823533382 |

Observed orders against the independent analytic growth are 2.01522, 2.00121
and 2.00267 for successive refinement pairs. These measurements support
approximately quadratic convergence without presuming it in the fit.

The free fit minimizes complex least-squares residuals for
`A+C64*(64/N)^p`, equivalent to `gamma_inf+C*N^(-p)` with
C=C64*64^p. Centred linear least squares determines A/C64 at each p; a coarse
search and golden-section minimization determine p independently. Checked
synthetic complex power laws with p=1.3, 2 and 3.2 test this reader in both
optimized/debug builds. Both readers produce identical physical-data tables.

| Growth model | Fitted p | Extrapolated gamma_inf | Difference from analytic | RMS residual |
|---|---:|---:|---:|---:|
| Free order, all four grids | 2.019770 | .009418495292 | -5.8172e-6 | 6.8668e-7 |
| Fixed order two, all four grids | 2 | .009430192755 | +5.8803e-6 | 2.7276e-6 |

A pure quadratic model describes the leading error well, but it is not an
exact four-point model. Its maximum residual is 4.35e-6 (about .14% of the
coarsest growth error), consistent with higher-order and remaining fixed-domain
terms. The free fit reduces that residual but adds a parameter; four grids
provide little information about model uncertainty. Its fitted intercept is
not automatically more accurate than a recent-pair extrapolation.

| Pair | p=2 Richardson limit | Difference from analytic |
|---|---:|---:|
| 64,96 | .009439350473 | +1.50380e-5 |
| 96,128 | .009424918281 | +6.05808e-7 |
| 128,192 | .009424970080 | +6.57608e-7 |

The two independently added grids therefore support the coarse extrapolation's
near-analytic limit and confirm its approximate order, while showing that its
particular intercept is not final. Recent-pair limits are much closer to the
analytic growth, but the actual 192×384 root remains 3.57% low. No extrapolated
value is substituted for an independently solved root or treated as a rigorous
infinite-resolution error bound.

## 4. Proton/electron yy and determinant errors at the fixed frequency

The table gives signed real parts of numerical minus analytic values at the
identical z_A. Full complex values, including signs and tiny imaginary parts,
are retained in `resolution.dat` and the component tables. The determinant
comparison includes absolute error because analytic D(z_A) is nearly zero.

| N | delta proton yy | delta electron yy | delta D |
|---:|---:|---:|---:|
| 64 | -.00362930524 | -.00176622496 | -1.87852069e-6 |
| 96 | -.00160914359 | -.000783773003 | -8.36252955e-7 |
| 128 | -.000905914443 | -.000440461479 | -4.71169611e-7 |
| 192 | -.000402719710 | -.000195557238 | -2.09526239e-7 |

Successive magnitude orders are 2.00594/1.99704/1.99944 for proton yy,
2.00382/2.00324/2.00256 for electron yy and
1.99600/1.99426/1.99862 for delta D. At N=192 the imaginary errors are
`-2.417e-15`, `-2.836e-17`, and `-1.722e-16`, respectively. They are retained
but do not establish a separate imaginary-error convergence order. The common
complex-fit order is essentially determined by the dominant real error.

| Error quantity | Free fitted p | Free-fit asymptotic real error | Fixed-p=2 asymptotic real error |
|---|---:|---:|---:|
| Proton yy | 2.008368 | -3.9690e-6 | +1.9447e-6 |
| Electron yy | 2.004248 | -2.5690e-7 | +1.2078e-6 |
| D | 1.995225 | +2.4909e-10 | -1.5085e-9 |

Latest-pair p=2 limits of these errors are -1.6392e-7, +3.6615e-7 and
-2.1154e-10. Their reduction supports convergence of the tensor/determinant,
not only of a root-finder return. Intercept differences between fit models and
pairs remain visible; fixed four-width support has not been taken to infinity.

## 5. Controlled 96×192 derivative substitution

For the actual ALPS reference momentum, the unit-density Gaussian is

```
F = pi^(-3/2)/[(3*m)^(3/2)*A] * exp[-p_parallel²/(3*m)-p_perp²/(3*m*A)]
dF/dp_perp = -2*p_perp/(3*m*A) * F
dF/dp_parallel = -2*p_parallel/(3*m) * F.
```

Here A=2 for protons and 1 for electrons. The isolated worker hook checks mass,
signed charge, density, support, every coordinate and the original file values
before using these derivatives. The analytic formulas are evaluated in double
precision; "exact" distinguishes them from finite differences, not floating-point
roundoff. Maximum file/Gaussian difference normalized by
its peak is exactly zero for both species on this machine. The incoming nodes,
f0 values and normal derivatives exactly reproduce the archived 96×192 data.
Normal derivative relative L2 errors are .00100182 (perpendicular) and .00111997
(parallel) for each species. These small global derivative errors can have a
much larger effect on a weak, cancellation-sensitive tensor response.

The replacement occurs after normal fitting and distribution to workers.
Only active interior df0 entries for the assigned species change; original
files, normalization, moment diagnostics, endpoint conventions, continuation
fits, nodes, PWL representation, quadrature and harmonic ownership do not.
No derivative is recomputed as frequency changes. The production derivative
routine is untouched.

| 96×192 variant | Weak growth | Error from analytic | Growth bias removed |
|---|---:|---:|---:|
| Normal numerical derivatives | .008075085805121 | -.001349226667795 | — |
| Exact nodal derivatives | .009190838793755 | -.000233473679162 | 82.7% |

The exact-derivative ordinary root has a measured [.009,.010] sign bracket,
real part 8.16e-14, |D|=2.66e-15 and |raw D/z²|=3.14e-11. Its scaled reciprocal
wave condition indicator is 1.23e-12. All original root/trajectory qualifications
pass; it is not an apparent near-origin root from absolute determinant stopping.

At z_A the signed real errors are:

| Quantity | Normal error | Exact-derivative remaining error | Absolute error removed |
|---|---:|---:|---:|
| Proton yy | -.00160914359 | -.000337240425 | 79.0% |
| Electron yy | -.000783773003 | -.0000817015101 | 89.6% |
| D | -8.36252955e-7 | -1.45089751e-7 | 82.7% |

This is a controlled effect of the supplied derivative arrays, not a claim of
an additive universal physical-error budget. The remaining weak-growth bias is
2.48% of the analytic rate. It combines the remaining PWL representation,
perpendicular quadrature and fixed support; these are not subdivided by a new
scheme. Exact derivatives do not make the finite-grid susceptibility exact.

Matched observations at z_A, the retained normal root and the independently
solved exact-derivative root compare both variants and the analytic response
at all three identical frequencies. `derivative_components.dat` retains all
six components for each species, and `derivative_determinants.dat` retains raw
complex D and signed changes. At the normal root, the exact-derivative D is
about +5.1148e-7; at the exact-derivative root, the normal D is about -6.5825e-7.
Thus the variants genuinely have different discrete-model zeros, not only a
reporting or absolute-threshold difference. Normal/analytic fixed-point
reproduction against the preserved C2f.4a tensors and determinants is exact.

## 6. Two lower-rectangle Nyquist checks

Only [-1,1] x [.001,.2] at the two new normal-derivative grids is evaluated
with A1. No original rectangle or derivative-control A1 calculation is added.
The independent analytic reference count remains three.

| Grid | Accepted winding | A1 status | Passes | Calls on each rank | Refinements | Maximum depth | Wall seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| 128×256 | 3 | 0 | 3 | 622 | 87 | 9 | 563.31 |
| 192×384 | 3 | 0 | 3 | 624 | 88 | 9 | 930.10 |

Both total phases are `18.849555921538759` (6*pi). Accepted winding, global
sampling convergence and physical call equality are checked independently of
the ordinary root returns. Raw samples and complete summaries are retained.
These counts concern the prescribed finite rectangle, not all frequencies.

The final study has **11 MPI jobs, 44 rank records and 1,324 collective disp
calls, with zero rank-call mismatches**. Of these, 1,246 calls are the two A1
calculations; 78 are focused probes, ordinary solves and matched comparisons.
Harmonic limits remain 43/19, with the original proton 0..21 and 22..43 worker
groups and electron 0..19 group. There are no new communicator groups or
frequency/wavevector parallelism.

The corrected frozen entry has exit=0 and the final COMPLETE line. All nine
recorded gates pass. Optimized/debug read-only tables and free-order fits are
byte-identical. New executed validation sources and entry match their frozen
copies. Preservation checks pass for 121 original source/report files, 805
recent historical files and 409 frozen source files. Production source hashes
are unchanged; the pre-existing working-tree changes are preserved. The
earlier reports and C2f.4a reference artifacts are not rewritten.

Compact evidence, inputs, original/exact nodal arrays, generated copied-source
variants, provenance, scalar/component/matrix tables, root qualification,
conditioning, A1 samples and per-rank counters are under
`tests/data/nyquist_c2f4a1/`. Its README explains the formats and separate
rejected preflight evidence; SHA256SUMS checks retained artifact integrity.

## 7. Remaining accuracy limitations and disposition

The measurements confirm approximately second-order weak-mode, yy and
determinant convergence and show that numerical differentiation accounts for
most of the measured 96×192 bias. The remaining fixed-grid error is material:
3.57% in the normal 192×384 growth and 2.48% in the 96×192 exact-derivative
control. Extrapolation improves agreement but is model dependent, with only
four resolved grids and fixed finite support.

**Disposition: success for the requested characterization**, with finite-grid
and extrapolation limitations retained. Numerical differentiation is the
largest measured contributor to this fixture's weak-mode error, and the
unchanged remaining representation/quadrature still matters. The data do not
establish an exact continuum limit, prescribe a universal grid size, or justify
adopting the analytic-derivative diagnostic as a numerical-VDF method.

No production derivative improvement, quadrature change, new susceptibility
formulation, benchmark, non-Maxwellian distribution, scan, growth maximization
or damped causal continuation is implemented. The exact-derivative control is
strictly a diagnostic, not a production method. Stop after this study.
