# C2d — Structure-preserving fixed-cell susceptibility

## Decision and scope

The coupled parallel formulation eliminates the C2c static discretization
residual. Its full-tensor extension converges toward the analytic growing-mode
response and gives independently accepted winding two on both requested boxes
at all three momentum resolutions. **This is a validated experimental
upper-half-plane formulation, not a production-default replacement.**

A smaller, distinct near-origin defect remains with the existing Bessel
approximation. Its failure of the Bessel sum identity leaves a static tensor
residual and a growing numerical zero near 1.38e-5. An isolated attribution
build using intrinsic Bessel values removes the measured static determinant
floor. Production Bessel routines are untouched. Damped causal continuation
also remains a separate unresolved requirement. Therefore the adoption decision
is **partial success; retain the prototype and the C2b production default**.

All new code is under `validation/c2d/`, with entry point
`validation/run_nyquist_c2d.sh` and compact results under
`tests/data/nyquist_c2d/`. No production source or input default was changed by
C2d. In particular the pre-existing C2b edits in `src/ALPS_fns.f90` were
preserved, not reset. SHA-256 comparisons verify preservation of that source,
`ALPS_fns_rel.f90`, `ALPS_com.f90`, A1 and A3 throughout these experiments.
Comparison with the identified C2c build manifest also verifies all 13
production Fortran sources in that manifest; the C2c kernel hash is unchanged.
The existing MPI decomposition, determinant assembly, differentiation, fits,
Landau integration, root algorithms and scan machinery remain authoritative.

## Coupled parallel identity and source verification

Write u=p_perp, p=p_parallel, A=df/dp_perp, B=df/dp_parallel,
S=uB-pA, d=m*w-k*p-n*q, and D_n=n*q*A+k*u*B. The actual `resU` is

```
q*(w*A+k*S/m)/d = q/m*(A+D_n/d).
```

`int_T` has form C_n(u)*p^l: l=0 for xx, yy, xy; l=1 for xz, yz;
l=2 for zz. Its Bessel and imaginary factors remain those of the actual
`int_T`. `disp` multiplies the resulting species tensor by n_s*q_s after
summing the existing worker harmonic ranges. All these normalizations remain
unchanged in the copied builds.

For n=0, the zz resonant term and the nonresonant `int_ee` integrand obey

```
q*p^2*(w*A+k*(u*B-p*A)/m)/(m*w-k*p) + q*p*(u*B-p*A)/m
    = q*w*p*u*B/(m*w-k*p).
```

This equality is pointwise and does not require integration by parts, a
Maxwellian, or a new derivative approximation. The combined integral is

```
-q*w*u/k * integral[p*B(p)/(p-m*w/k)] dp.
```

The explicit w factor is retained. There is no fitted or measured static
subtraction. At w=0 the combined result is exactly zero in the scalar test.
The actual `resU`/`int_T` equality has maximum absolute residual 2.84e-14 in
the prescribed parallel fixtures and 5.12e-16 in the full tensor fixtures.
Positive and negative real frequencies were tested at positive k=0.3;
imaginary frequencies span 1e-2 through 1e-10.

Before the correction, the actual `funct_g_node`/`int_ee` synthetic static
experiment was rerun. Its unchanged-nodal cancellation is within rounding,
but complete-numerator PWL integration plus unchanged ee gives:

| Parallel intervals | Reproduced naive static residual | Coupled residual |
|---:|---:|---:|
| 64 | -6.418143596134129e-4 | 0 |
| 128 | -1.616090865725894e-4 | 0 |
| 192 | -7.199702955351661e-5 | 0 |

`math.sh` runs both the actual C2c baseline and the coupled checkpoint.
`test_convergence.f90` independently compares constant, complex affine,
quadratic, cubic and Gaussian data to a fine Simpson reference; affine cases
also have an independent complex-log primitive check in `test_parallel`.
The smooth nonlinear fixtures exhibit second-order PWL convergence. Neither
the numerical VDF nor its differentiation is replaced with analytic data in
the physical calculations.

## Full-tensor discretization

A coherent interpolation basis is essential. The experimental model uses
fixed PWL interpolants of the **existing nodal derivative arrays** A and B.
Multiplication by p and p^2 is then performed algebraically and integrated
exactly. Thus it does not use a PWL interpolant of each already-multiplied
complete tensor numerator, as the naive C2c B prototype did. No derivative
is recomputed, fitted or replaced; the change is in how those nodal data
represent the coupled integrands. This distinction matters for static moment
identities and is part of the proposed numerical formulation.

Let M_l(F)=integral p^l*PWL(F) dp and
C_l(F;a)=integral p^l*PWL(F)/(p-a) dp. Then

```
I_n,l = q*C_n/m * [M_l(A) - C_l(D_n;a)/k]
a = (m*w-n*q)/k
C_1 = M_0 + a*C_0
C_2 = M_1 + a*M_0 + a^2*C_0.
```

M_l is evaluated exactly cell by cell in local coordinates. C_0 uses the
unchanged C2c `cell_cauchy`/`domain_cauchy` kernel. This is exact integration
of polynomial factors multiplying the chosen PWL basis; it is not an assertion
that the underlying continuous VDF is represented exactly.

For zz, split E=q/m*[u*M_1(B)-M_2(A)] into J0^2*E and
(1-J0^2)*E. Combine the former with the n=0 resonant term before evaluation;
retain the latter as the complementary ee term. This rearrangement is valid
for the finite harmonic sum, without imposing an infinite Bessel sum on a
truncated calculation. The unchanged perpendicular quadrature uses weight
1 at the first usable and interior nodes, 1/2 at the last usable node,
and the common factor 2*pi*dp_perp. Parallel support remains the usable
nodes 1 through npar-1. No cell or support tail is added silently.

`FORMULATION.md` records this proposal before implementation. The copied
`full_integrate` calls the prototype only for nonrelativistic species with
Im(w)>0. Real and damped production evaluation remains C2b. No Tlim,
resonance-window, real-resonance support or stencil selection enters the
experimental upper-half-plane integral.

### Tensor identities and their limits

The tensor can be written v v^dagger with
v=(n*Jn/z, -i*u*Jn', p*Jn), z=k_perp/q. Its contraction with
(k_perp,0,k) is (m*w-d)*Jn. For the complete harmonic sum, the Bessel
identities sum Jn^2=1, sum n*Jn^2=0 and sum Jn*Jn'=0 imply the static
longitudinal cancellation with E. A common linear integration functional
preserves these algebraic cancellations; independently approximated products
generally do not. The same argument does not say every tensor component
vanishes separately, or justify replacing finite harmonic tails by zero.

`test_tensor` checks every nodal rearrangement against actual ALPS `resU`
and `int_T`, and integrates all six modes. With accurate prescribed Bessel
values its static norm |chi*k|_max is 1.00614e-16. At w=i*gamma the residual
scales from 2.09592e-2 at gamma=1e-2 to 2.21226e-10 at gamma=1e-10.
This checks the finite-k_perp extension independently of a successful winding.

### Additional source inconsistency, retained in production

At `int_ee`'s first usable perpendicular/first usable parallel corner,
the original code uses `u*df_perp` where its other corners and its continuum
formula use `p*df_perp`. Its extra contribution relative to the continuum
corner is

```
pi*q/m*dp_perp*dp_parallel*p*(p-u)*df_perp.
```

The experimental moment formulation follows the continuum expression at
that corner. `diagnose_structure` measures the difference from the actual
prepared derivative arrays. No perpendicular endpoint weight is changed.
This is a separately identified source issue, not a static offset correction.
Its measured extra ee contributions are -5.05089e-9, -4.00404e-10 and
-1.01397e-10 on the three reference grids.

## Physical calculations

The equilibrium, finite support, normal differentiation and method-2
continuation are the existing C2b/C2c fixtures: beta_parallel,p=3, A_p=2,
isotropic analytic electrons with beta_parallel,e=3, mass ratio
m_e/m_p=5.44617021487e-4, equal densities, no drift, vA=1e-4,
(k_perp,k_parallel)=(0.001,0.3). Grids are 32x64, 64x128 and 96x192,
with four Gaussian widths of support. Four MPI ranks retain their original
harmonic assignments. No negative-k development or new physical backend is
introduced.

### Growing roots

The ordinary solver independently locates both signed-frequency partners;
the root algorithm and residual observer are unchanged. The analytic positive
reference is 0.39666286088749 + 0.08598637345950 i.

| Grid | Corrected Re(w)>0 | Corrected Im(w) |
|---|---:|---:|
| 32x64 | 0.389163072707080 | 0.085611504527144 |
| 64x128 | 0.394779307385325 | 0.085883953153049 |
| 96x192 | 0.395825246940874 | 0.085939984231271 |

The real-frequency error decreases approximately quadratically with grid
spacing. The negative partners are solved separately at positive k_parallel;
reflection errors are a few e-15 or smaller. Returned determinant residuals
are recorded, not inferred from the printed root precision.

For comparison, production C2b positive roots are approximately
(0.390581513199,0.086311784789), (0.395411724579,0.086220209553),
(0.396168900706,0.086138999817). Naive C2c B roots are
(0.390048496173,0.086275452821), (0.394997860684,0.086052447992),
(0.395922143763,0.086015080415). The coherent formulation is not uniformly
closer for every root component at every resolution. Its advantages must be
judged together with static consistency and determinant-trajectory convergence,
not by selecting the closest isolated root.

### Matched contours and independent A1 calculations

The 257 matched coordinates are the retained C2c coordinates; no trajectories
with different adaptive coordinates are compared pointwise. The lower fixture
moves the bottom edge to gamma=0.001 while retaining the original vertical
interior coordinates, a valid nonuniform sampling of the lower rectangle.
Independent adaptive evaluations use the normal rectangle constructor.

| Grid | Original median/max relative D error | Lower median/max relative D error |
|---|---|---|
| 32x64 | 0.0127795 / 0.274639 | 0.0127795 / 0.528795 |
| 64x128 | 0.00322331 / 0.0675132 | 0.00322331 / 0.129543 |
| 96x192 | 0.00143553 / 0.0298146 | 0.00143553 / 0.0586808 |

Errors decrease systematically with refinement. Large relative errors are
not hidden: absolute errors, minimum |D| and cumulative phase comparisons
are retained in the compact diagnostics. For context, the naive C2c B maximum
errors were 1.66419, 0.397650, 0.174624 on the original contour and
192.900, 46.1587, 20.0862 on the lower contour.

All six adaptive calculations accept winding **2**, A1 status **0**, with
three independent/global convergence passes. Original-box calls are 616,
630 and 630; lower-box calls are 702, 690 and 686, respectively. Maximum
depths are 6 on the original box and 9 on the lower box.
All four ranks execute exactly the same number of physical determinant calls.
These are finite-domain results of the represented numerical susceptibility,
not a complete global physical instability inventory.

### Low-frequency response and additional zeros

Before evaluating the corrected model, the original C2c extra roots were
reproduced through the ordinary solver:

| Grid | Reproduced C2c additional growth | Corrected smaller zero |
|---|---:|---:|
| 32x64 | 0.0126394925345873 | 0.0000138702602862192 |
| 64x128 | 0.00647826428106220 | 0.0000138359488691285 |
| 96x192 | 0.00435308676425105 | 0.0000138296116456261 |

The O(h) growing zero is not simply hidden by a new frequency cutoff: the
low-frequency probes extend through gamma=1e-8, and the smaller zero is
independently located with residuals around 1e-19. It has a different,
nearly resolution-independent scaling. Its persistence prohibits claiming
that the complete inherited numerical tensor is now artifact-free.

At w=i*1e-8 the 32x64 naive C2c B determinant was about +3.03125e-4;
the coherent result is +2.91857e-10. With accurate Bessel values only in the
isolated constructor attribution experiment, it is -1.51997e-16, approaching
the analytic -1.54618e-16. The Bessel-control solver can stop near the
origin because D itself vanishes there; such near-zero solver returns are
not interpreted as a resolved additional unstable mode.

The existing `BESSJ` approximation gives sum Jn^2-1 about 5.6628e-9 at
small nonzero argument, despite returning the exact special value at x=0.
The corresponding static longitudinal defect is proportional to
k*q/m*integral[(1-sum Jn^2)*(u*M_1(B)-M_2(A))] with the existing
perpendicular weights. `structure_defect.dat` compares this prediction with
the assembled tensor. The attribution experiment changes only a copied
constructor, not `ALPS_fns_rel.f90`, harmonic limits, assignments or the
candidate default. It is diagnostic evidence for a separate Bessel-accuracy
review, not permission to substitute a new production backend.

The predicted static |chi*k|_z residuals are 8.4716089e-10, 8.4884335e-10
and 8.4915193e-10. The assembled values at w=i*1e-8 are 8.4715773e-10,
8.4884029e-10 and 8.4914889e-10: agreement within about 3.2e-15, including
the small finite-frequency response. This identifies the remaining tensor
defect quantitatively, rather than only correlating it with the Bessel control.

## Analyticity and continuation

For fixed derivative arrays and momentum nodes, A, B, D_n and their moments
are frequency-independent. Each fixed-support Cauchy transform is analytic
for Im(a)>0, and a is affine in w with positive k_parallel. Polynomial
frequency factors and finite tensor/harmonic sums preserve that analyticity.
Thus the represented upper-half-plane model is holomorphic, including near
internal grid nodes and outside the real momentum support. This mathematical
property does not prove convergence to the continuous VDF susceptibility.

`test_tensor` independently perturbs real and imaginary frequency. Its zz
CR residual decreases from 8.56e-5 at step 1e-3 to 8.56e-9 at step 1e-5
and 4.04e-11 at step 1e-6. Physical determinant CR tests use the retained
C2c groups, including support-edge and near-axis points, at all three grids;
finite-difference truncation decreases before roundoff dominates. Additional
near-origin and growing-root tensor tests are retained separately. No claim
of exact differentiability is made for floating-point arithmetic or for the
arithmetic-series branch inside the stable logarithm evaluator.

At w=i*1e-5 the full proton-tensor CR residual is 2.79e-10 at step 1e-6;
the determinant residual decreases from 4.23e-7 to 4.30e-9 to 3.51e-10 as
the step decreases from 1e-6 to 1e-8. At the analytic growing-root coordinate,
the tensor residual decreases by roughly 100 per tenfold step reduction
(6.52e-4, 6.53e-6, 6.53e-8). Near-origin tensor differences eventually become
roundoff dominated; decreasing the step indefinitely is not an acceptance test.

For causal damped continuation the direct Cauchy transform is insufficient.
The represented boundary residue is

```
G_n(r) = -q*C_n*r^l/(m*k) * [n*q*A_PWL(r)+k*u*B_PWL(r)].
```

The n=0 combined zz residue has the equivalent explicit-w form. Its boundary
jump is 2*pi*i*G_n(r), with the unchanged perpendicular weights. This differs
from both the complete-numerator PWL residue of C2c and, in general, the
existing fitted analytic-continuation residue. Boundary-only copied builds
compare direct integrals on both sides, this grid residue, and the unchanged
`landau_integrate` at the same frequencies near the isotropic control
Re(w)=0.4046. They are diagnostics, not a new physical damped backend.

For n=0 zz, the direct/grid jump discrepancy at separation 2e-12 is
3.83e-12, 3.89e-12 and 3.99e-12 across the three resolutions. All 171 recorded
harmonic/component jumps pass the independent residue check. The unchanged
fitted residue nevertheless differs from the corrected grid residue by
1.4113%, 0.40629% and 0.17092% in this component. Corresponding finite
causal-boundary mismatches are 0.0253369, 0.00738149 and 0.00311376 in the
integral normalization. Their convergence does not establish global analyticity
of the finite-resolution continued model.

The ordinary copied candidate deliberately retains C2b at real/damped
frequencies. Exact comparison of its fitted Landau logs with the unchanged
C2b logs verifies preservation. A future causal extension needs an analytic
momentum representation consistent with the upper-half-plane discretization;
a cell selector based on Re(a) would not establish that property. No damped
roots of the direct fixed-cell model are claimed as physical roots here.

## Weight reuse and performance

After the mathematical checkpoints, the candidate was evaluated both with
and without reusable Cauchy weights. For fixed species, harmonic and frequency,
the same a applies to every perpendicular node and tensor mode. Integrating
the two endpoint hat functions once per parallel cell gives a linear Cauchy
weight vector. Applying it with `sum(weights*values)` (without complex
conjugation) reproduces the original integral. Exact moments are unchanged.

The worker-local cache stores only kernel weights. Frequency, species,
parallel coordinates and k_parallel changes invalidate it; it introduces
no determinant interpolation, new MPI requests, or changed harmonic assignment.
It is separate from the untouched A3 exact-frequency determinant cache.
No per-cell or per-node allocation is introduced. The full npar=128 cache
occupies 64,132 bytes per numerical worker, plus automatic numerator arrays.

Matched 32x64 weights-versus-unoptimized determinants differ by at most
2.09e-14 relative on the original contour. Comparisons also cover the lower
contour and the 64x128 grid, with an absolute tolerance scaled by determinant
magnitude. Debug and optimized numerical comparisons are recorded separately.
Counters distinguish derivative-node reads from Cauchy kernel calls; a
node read is not a complete numerator construction or a collective call.
The derivative counter covers A/B vector construction in the six resonant
tensor modes; it excludes the complementary ee moment's additional reads.

The repeated benchmark uses 257 identical frequencies, four ranks, one
OpenBLAS/OpenMP thread, disabled diagnostic I/O, and sequential launches.
It compares C2b, naive C2c B, coherent without reuse and coherent with reuse.
Timing and operation-count tables accompany this report. These analytic-electron,
numerical-proton timings are not general numerical-VDF scan performance claims.

| Method | Mean seconds / 257 determinants | Range | Seconds / determinant |
|---|---:|---:|---:|
| Production C2b | 10.2372 | 10.0919–10.4301 | 0.0398335 |
| Naive C2c B | 71.9555 | 71.7375–72.1261 | 0.279983 |
| Coherent, no reuse | 76.0647 | 75.9418–76.1310 | 0.295972 |
| Coherent, weight reuse | 11.7084 | 11.1954–12.0445 | 0.0455580 |

Weight reuse improves this formulation by 6.50 times and is 6.15 times faster
than naive C2c B; it remains 14.4% slower than C2b. All methods make exactly
257 collective determinant calls. Across numerical workers the cell-kernel
calls fall from 373,332,078 to 2,007,684; resonant derivative-node reads remain
752,590,062. Reuse changes the work inside a determinant, not the logical
frequency requests or MPI evaluation sequence. The cached implementation
constructs 7,967 harmonic/frequency weight sets in this benchmark.

The polynomial-division evaluation of C_1/C_2 can lose relative accuracy when
|a| is extremely large because of moment cancellation. The C2c scalar kernel's
extreme far-field tests do not certify the entire rearranged tensor at such
frequencies. This study covers the stated physical boxes, harmonics and
support-edge tests; wide-frequency horizontal asymptotics need a separate
conditioning check before using this experimental formulation there.

## Validation and reproducibility

Compiler: GNU Fortran 11.4.0, Open MPI 4.1.2, existing LAPACK/BLAS linkage.
The measured host is an Intel Xeon W-2133 at 3.60 GHz (12 logical CPUs).
Optimized flags
are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags are
`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
Scalar mathematical tests also trap invalid/zero/overflow in optimized mode.
All MPI physical tests use four ranks and one collective determinant path.

The 64x128 debug and optimized matched trajectories, growing-root records and
adaptive sample files are byte-identical. Accepted counts, calls, passes,
refinements and maximum depths also agree exactly. The independent optimized
uncached-versus-cached trajectory comparisons have maximum relative differences
of 2.09e-14 (32x64 original), 1.47e-14 (64x128 original) and 5.08e-14
(64x128 lower); these small regrouping effects do not change accepted counts.

Run the self-contained entry point:

```bash
bash validation/run_nyquist_c2d.sh
```

It prepares normal generated VDF fixtures and analytic reference coordinates
in an isolated scratch directory. `C2D_SCRATCH` selects that directory.
For an already identified C2c study, set `C2D_REFERENCE_STUDY` and
`C2D_REFERENCE_VARIANTS`; this reuses VDF inputs/common reference objects,
not determinant values in a new candidate. No permanent VDF dumps or second
production susceptibility were added.

The mathematical gates cover actual static cancellation, actual tensor source
identities, constant/complex affine/polynomial/Gaussian integrals, source Bessel
identities, complex derivatives, and weight equivalence. The physical matrix
includes both contours and signed-real-frequency roots at all three grids,
low-frequency components/determinants, causal-boundary diagnostics, and a
representative debug 64x128 physical gate. Legacy analytic A1/A2/A3/log/mirror/
line/plane suites and ordinary scientific-file comparisons are recorded below.

| Regression gate | Optimized | Debug |
|---|---|---|
| A1 standalone (306 checks, no failures) | PASS | PASS |
| C2b remainder, interpolation, partition, real-node tests | PASS | PASS |
| C2c legacy-core identity and scalar Cauchy tests | PASS | PASS |
| C2d parallel/tensor identities, convergence, weight invalidation | PASS | PASS |
| A2 integration/root/map regressions | PASS | PASS |
| A3 growth/cache/strip/failure regressions | PASS | PASS |
| Logarithmic growth and weak-growth regressions | PASS | PASS |
| Mirror-mode regressions | PASS | PASS |
| B1 line-scan and output regressions | PASS | PASS |
| B2 plane-scan, MPI fault recovery and output regressions | PASS | PASS |

The six warm suites were launched sequentially in each configuration. Their
existing checks were not weakened. Nine ordinary scientific files per build
(roots, map, line and plane outputs) are byte-identical to the identified
pre-C2d executable: eighteen successful file comparisons. New sources introduce
no communicators, frequency workers or wavevector workers. Each Nyquist job
checks identical physical dispersion-call counts across all four ranks.

An additional numerical isotropic 64x128 control accepts winding zero in both
boxes, using 630 and 714 calls. Replaying the unmodified production numerical
64x128 unstable contour gives winding two in 628 calls. Production numerical
growing and damped root records are byte-identical to their C2c baseline records;
the positive damped root remains approximately
0.405337298402431 - 0.001709481504305 i. It is not a new fixed-cell damped result.
The analytic damped reference remains 0.40460138158784 - 0.00185430537032 i.

The mathematical/scientific matrix was executed through its component scripts
with the identified C2c optimized/debug common builds. A separate fresh run of
`prepare.sh` regenerated all three proton VDF files and both analytic matched
reference curves byte-identically. Its regenerated 32x64 candidate trajectory
also matches the recorded trajectory byte-for-byte. The standalone entry point
packages this sequence for reproduction; it was not redundantly rerun as a
second full monolithic study. Dedicated relativistic physics and numerical-
electron convergence suites were not rerun or claimed as new validation.

An initial exploratory root replay inherited a one-line `&guess_1` rejected
by the custom reader. Those empty-root attempts were not counted as evidence;
the checked driver writes multiline guesses, requires `roots_full.dat`, and
reproduced all three baselines before accepting results. One long unoptimized
matrix was stopped after retaining its completed comparison points; the full
physical matrix was then run with the checked weight-reuse version. Repeated
unoptimized timings still use the identical coherent formulation.

## Recommendation

Retain the coherent upper-half-plane prototype for further development. It
corrects the principal C2c quadrature incompatibility, has controlled momentum
convergence and an analytic fixed-domain mathematical representation, and
supports substantial kernel reuse without changing MPI work assignments.
Before production adoption, review the inherited Bessel sum-identity error
and the `int_ee` corner inconsistency, verify the resulting full static tensor
with the chosen Bessel implementation, and separately design a consistent
causal continuation. Do not normalize away the small numerical zero or select
bounds to exclude it silently. Numerical electrons, growth searches and
numerical-VDF wavevector scans remain outside this experiment.
