# C2f.4b — one tabulated two-component proton distribution

**Result: success for this one-mixture validation.** The final frozen entry completed from start to finish with exit zero: 50 scientific gates passed, 62 sequential four-rank jobs completed, and 21,188 collective determinant calls had no per-rank count mismatches. Optimized and debug results reproduce the same accepted inventories and matched tensors. Production sources and historical fixtures remain unchanged.

## Construction and independent normalization

The physical numerical calculation has `nspec=2`: one unit-density numerical proton species read from a single ordinary `.1.array` table, and analytic isotropic electrons. The proton table is formed externally as `0.8*f1 + 0.2*f2`, with unit-density, zero-drift component shapes `(beta_parallel, A)=(3,2)` and `(1.5,1)`. There is no post-construction renormalization. Electrons retain density one, beta three, anisotropy one, mass `5.44617021487e-4`, and negative unit charge. The wavevector is `(0.001,0.3)` and `vA=1e-4`.

The ordinary unmodified ALPS distribution generator produces the two Gaussian component files on common nodes; `validation/c2f4b/mix.f90` writes their weighted sum in the ordinary input format. Its two generator indices are not separate physical proton species. ALPS reads the resulting single table and constructs its normal numerical derivatives and continuation fit.

Common support is `0 <= pperp <= 4*sqrt(6)` and `|ppar| <= 4*sqrt(3)`, fixed across 32×64, 64×128 and 96×192 grids. Four widths mean four Gaussian exponent scales of the widest component in each direction. The cooler component has eight perpendicular widths and `4*sqrt(2)` parallel widths. Infinite-domain mixture moments are density one, zero drift, parallel beta 2.7, perpendicular beta 5.1, anisotropy 17/9, and parallel kurtosis 3.14814814814815. The latter independently distinguishes the mixture from a single bi-Maxwellian. Analytic finite-support density is 0.999999897638055; moment diagnostics retain both rectangle and trapezoidal estimates separately from this truncation loss.

The independent analytic reference alone uses `nspec=3`, all analytic: `(density,beta_parallel,A)=(0.8,2.4,2),(0.2,0.3,1),(1,3,1)`. In `ALPS_NHDS`, the parallel thermal speed satisfies `vtherm**2=beta/(density*mass)`, so partial betas are density times the requested unit-density shape parameters. `spec_read` retains the specified densities. This reference is not computed by the coherent integrator.

At identical frequencies, native partial proton tensors are compared independently with 0.8/0.2 times unit-density analytic responses. The vacuum matrix is reconstructed once from a unit-component run; the second electron response is not added. No determinants are added. A reporting-only observer retains all three native analytic tensors in `populations.dat` and joins the proton tensors only in the existing two-channel diagnostic format. Native collective `disp`, its wave matrix, determinant, ordinary roots and A1 evaluations remain unchanged.

## Source construction and preservation

New code is confined to `validation/c2f4b/` and its shell entry point. The existing C2f.4a source builder supplies the coherent C2d formulation, reusable weights and C2f.3a restricted extended-precision Cauchy arithmetic. The new builder changes only read-only startup/tensor observers: the old single-Gaussian diagnostic is inappropriate for a mixture, and the analytic three-species reference needs an additional reporting channel. The numerical observer explicitly rejects any numerical configuration other than one unit-density proton at species one with `nspec=2`.

The C2f.1-qualified production files are independently reconstructed from preserved before-sources and the two reviewed patches, then compared byte for byte with the current source. The entry freezes the source/builders/fixtures before execution and records revision, working tree, compiler flags, linked libraries, source hashes and immutable historical-report checks. No production source, mathematical kernel, derivative, Bessel routine, harmonic stopping rule, MPI assignment, root qualification, A1 acceptance criterion or scientific tolerance is modified.

## Reproduction

Run `bash validation/run_nyquist_c2f4b.sh`; optionally set `C2F4B_SCRATCH` to a fresh directory outside the checkout. The entry generates its own VDFs and references and runs four-rank MPI jobs sequentially. It does not depend on an earlier scratch study. The PROBE program's deliberate diagnostic exit one/status six is distinguished from accepted A1 results and final driver completion.

## Measured VDF and representation errors

The table below gives full-support composite-trapezoid moments, normalized by the measured density only when forming temperature moments. The input table itself is never rescaled.

| Grid | Integrated density | Parallel beta | Perpendicular beta | Relative L2 error df/dpperp | Relative L2 error df/dppar |
|---|---:|---:|---:|---:|---:|
| 32×64 | 0.995816736452 | 2.70189647633 | 5.12145429509 | 0.0226213 | 0.0142853 |
| 64×128 | 0.998957208419 | 2.70046888558 | 5.10531681309 | 0.00572415 | 0.00359901 |
| 96×192 | 0.999536732533 | 2.70020735068 | 5.10235490828 | 0.00254992 | 0.00160185 |

Derivative errors use cylindrical-weighted L2 norms over the actual interior derivative nodes against analytic derivatives of the prescribed mixture. These derivatives are diagnostic references only, never supplied to ALPS. The mixture formula matches the stored table within `1.4e-17`; positivity and reflection symmetry hold to rounding. The density deficit at these resolutions is dominated by radial-grid quadrature, not by the `1.02e-7` finite-support loss. Boundary exclusion in the actual susceptibility is retained; density diagnostics also record the production integral separately.

On the ordered 257-point original contour, median relative determinant errors are 0.0121741, 0.00308115 and 0.00137267. Maximum relative errors are 0.335778, 0.0802385 and 0.0353335. Maximum absolute raw determinant errors are 0.0991996, 0.0250983 and 0.0111747. Maximum cumulative phase discrepancies fall from 0.184910 to 0.0518156 to 0.0233312 radians, while the closed total phase agrees. These finite-grid errors are retained even though the independently adapted counts agree.

## Six-component tensor convergence

On the original ordered contour, median meaningful relative proton component errors are:

| Grid | xx | yy | zz | xy | xz | yz |
|---|---:|---:|---:|---:|---:|---:|
| 32×64 | 0.0130278 | 0.0130278 | 0.0104241 | 0.0106443 | 0.00832998 | 0.00794970 |
| 64×128 | 0.00324195 | 0.00324195 | 0.00258177 | 0.00261876 | 0.00208458 | 0.00199791 |
| 96×192 | 0.00144070 | 0.00144070 | 0.00114890 | 0.00116027 | 0.000925949 | 0.000888306 |

The six-component Euclidean relative error has median 0.0119384, 0.00295398 and 0.00131334, with observed orders 2.0149 and 1.9991. Maximum six-component norm errors are 0.0296629, 0.00724350 and 0.00322140. Lower-contour median orders are 2.0034 and 1.9991; maximum-error orders are 2.0754 and 1.8562. Thus quadratic convergence describes most measured errors, without imposing it on every local maximum. The unchanged analytic electron tensors match exactly across these reference/numerical comparisons. Near-zero reference components are excluded only from relative-error statistics using the established reader rule; their full complex values and absolute errors remain in the evidence.

The observers retain raw scaled ALPS tensors and wave matrices. Conventional susceptibility is obtained by dividing a species tensor by `(omega*vA)**2`. The common factor cancels from relative component and six-component norm comparisons at identical frequencies; absolute errors in the retained files refer to the raw scaled convention, not conventional susceptibility. Six complex components, the complete independent-entry wave matrix and complex determinant are retained at every shared coordinate, including analytic-root probes.

## Growing roots and contour inventory (optimized)

The native analytic three-species response independently gives
`+0.36324249734063563 + 0.07461168480752349 i` and
`-0.36324249734063574 + 0.074611684807523476 i`.
Ordinary root solving follows separate accepted winding-one boxes on the positive and negative branches. These local boxes do not change either prescribed full rectangle.

| Grid | Positive root real part | Growth | Complex root error | Raw determinant magnitude |
|---|---:|---:|---:|---:|
| 32×64 | 0.354728547179528 | 0.0748647267255221 | 0.00851770964 | 5.04e-15 |
| 64×128 | 0.361099403534481 | 0.0746691365598714 | 0.00214386375 | 2.47e-15 |
| 96×192 | 0.362288955267324 | 0.0746366409201302 | 0.000953868593 | 3.21e-15 |

Observed complex-root convergence orders are 1.9903 and 1.9973. The negative partners are solved independently and agree with drift-free reflection to roughly `4e-15` in frequency. Root acceptance uses the existing ordinary solver, matrix/residual qualification and independent local A1 enclosures; it is not inferred from small raw determinants alone.

| Response | Original winding / calls | Lower winding / calls | Refinements original / lower | Maximum depth original / lower |
|---|---:|---:|---:|---:|
| Independent analytic mixture | 2 / 646 | 2 / 712 | 99 / 132 | 6 / 9 |
| Numerical mixture 32×64 | 2 / 644 | 2 / 722 | 98 / 137 | 6 / 9 |
| Numerical mixture 64×128 | 2 / 646 | 2 / 714 | 99 / 133 | 6 / 9 |
| Numerical mixture 96×192 | 2 / 646 | 2 / 712 | 99 / 132 | 6 / 9 |

All entries have A1 status zero, evaluator status zero, three global convergence passes, and total phase `4*pi`. Both local branch boxes have winding one. No additional or missing zero is indicated inside either tested rectangle. These are finite-domain accepted numerical inventories, not a global stability classification or certification outside the contours.

The lower ordered contour has median determinant errors 0.0121741, 0.00308115 and 0.00137267, with maximum errors measured separately in the compact evidence. At 64×128 and 96×192 those maxima are 0.105356 and 0.0458628, and maximum cumulative phase discrepancies are 0.0277219 and 0.0140774 radians. A matching integer does not erase these finite-grid errors.

## Focused controls (optimized)

At nine exactly matched component frequencies, native analytic partial-density weighting agrees with the independently weighted unit-density construction to about `1e-15` relative determinant difference. The numerical mixture agrees with the sum of independently initialized endpoint responses to `2.07e-15` relative determinant difference. The maximum component/matrix discrepancy is `1.13e-15` in the raw scaled convention. Actual central derivative arrays are linear to `2.8e-17` absolute. Endpoints alpha one and zero exactly recover all 8,385 parsed component rows, including the original node coordinates. The alpha-one response also reproduces the preserved C2f.3 single-component reference to `1.28e-15` relative determinant difference at these frequencies. The alpha-zero endpoint is independently compared with its cooler analytic unit-density response, with its expected finite-grid error retained.

Cached/direct comparison uses 39 shared points spanning the original contour, both analytic-root frequencies and four additional probes. Maximum relative determinant difference is `1.315e-15`, passing the unchanged `2e-12*max(1e-6,abs(Dcached))` gate. There is no rounding, tolerance adjustment or different susceptibility formula to enforce agreement.

The original actual-source parallel and full-tensor fixtures pass. Exact static longitudinal residual is `9.89e-17`; maximum nodal algebraic identity residual is `5.12e-16`. At the positive analytic mixture root, the determinant Cauchy–Riemann residual decreases from `5.10e-5` to `5.10e-7` to `4.81e-9` as the central-difference step decreases by decades; the proton tensor norm residual decreases from `1.00e-5` to `1.00e-7` to `1.00e-9`. Signed-root, low-growth and imaginary-axis checks are retained too, including roundoff-limited small individual analytic-electron components. This supports numerical evaluation of the fixed-representation holomorphic response; finite-difference tests are not a proof that the tabulated VDF model equals the physical continuum response.

All three numerical grids retain proton harmonic limit 15, rank-one core 0–7 with halo -1–8, rank-two core 8–15 with halo 7–16, and analytic electron assignment on rank three. Stored Bessel values agree exactly with corrected production evaluation. Rank zero coordinates the unchanged collective determinant. No wavevector/frequency work distribution is introduced.

The reused convergence reader labels its six-independent-component norm as “Frobenius”; this report calls it the six-component Euclidean norm. It does not double-weight off-diagonal entries as a full 3×3 Frobenius norm would.

## Validation files

- `inputs.sh` emits generator, physical nspec=2, and independent analytic nspec=3 inputs with separate names; it does not source the older generic-input emitters.
- `mix.f90` performs only external table construction and independent Gaussian moment checks.
- `build.sh`, `startup.inc`, and `tensor.inc` reuse the existing experimental builder and add read-only mixture-aware observers. Actual derivative arrays, population tensors, harmonic ownership and Bessel halos are retained.
- `analyze.f90` checks analytic weighting, numerical superposition, endpoints and derivative diagnostics, and writes exact shared coordinate decks. It does not evaluate a second susceptibility.
- `prepare.sh`, `tools.sh`, and `math.sh` reuse ordinary generation, prior comparison readers and the actual-source mathematical fixtures.
- `driver.sh` and `collect.sh`, invoked by `validation/run_nyquist_c2f4b.sh`, perform the frozen sequential MPI run and aggregate scientific gates separately from expected diagnostic exits.

No other validation infrastructure is altered. The explicit imaginary-frequency floor and strictly growing root checks prevent interpreting this as a new damped coherent response.

## Scope and limitations

This is a genuinely non-Maxwellian **single numerical species** test, not a sum of two physical numerical-species runs. Separate component runs are controls only. Positivity, symmetry, moments, independently additive NHDS tensors, and numerical superposition provide complementary checks of input construction and response linearity.

The finest root is still about `9.54e-4` away from the analytic complex frequency, predominantly in its real part; the growth difference is `2.50e-5`. The fixed 96×192 representation is therefore not exact merely because its contour inventory agrees. Full-grid errors include normal finite differences, PWL representation and the unchanged perpendicular quadrature together. No analytic derivatives are injected and no new scheme is introduced to partition that residual further. Finite-support truncation remains fixed and nonzero, though small compared with the measured discretization errors.

The study covers this one smooth, zero-drift mixture, at one positive-kparallel nearly parallel wavevector and strictly positive-imaginary physical frequencies. It establishes neither general non-Maxwellian accuracy nor a consistent coherent damped continuation. There are no new wavevector scans, growth maximization, non-Maxwellian families, production selector or default change. Additional compatibility suites were not repeated: the construction changes reporting hooks only, while native analytic weighting, current-production source qualification, static/tensor fixtures and the historical numerical-proton endpoint provide the focused controls requested here.

## Final optimized/debug acceptance and reproducibility

The final entry was frozen before acceptance and executed cleanly from start to finish in `/tmp/alps-c2f4b-acceptance-final-20261124`. No running driver was edited. There were no acceptance retries, unexpected diagnostic exits, input failures or evaluator failures. The preceding pilot checked source construction, density conventions, endpoints, actual derivative observation and numerical linearity; it did not supply expected numerical mixture results to the acceptance gates.

Both builds pass all 25 recorded scientific gates apiece, plus the separately fatal mathematical, coordinate, input and source-protection checks. There are 32 independently accepted A1 contours across the two builds (the two standard rectangles and two local root enclosures for the analytic reference and each of three grids), eight ordinary paired-root jobs and 22 fixed-frequency diagnostic jobs. All 248 rank metric records agree on the 21,188 physical collective calls. Deliberate PROBE exit one/status six is retained as diagnostic behavior; A1 results have status zero and the complete driver has exit zero.

All three complete 326-coordinate matched tensor/wave-matrix/determinant streams and every ordinary paired-root file are byte-identical between optimized and debug builds. The parallel/full-tensor fixture outputs are also byte-identical. Both cache/direct comparisons have maximum relative difference `1.3145132896914457e-15`; both numerical-linearity gates pass. All standard contour status, winding, phase, pass count, refinement depth and physical-call results agree between builds. Final-pass sample minima for the original/lower rectangles are:

| Response | Original minimum abs(D) | Lower minimum abs(D) |
|---|---:|---:|
| Analytic mixture | 1.19215809701e-4 | 1.25293397732e-6 |
| 32×64 numerical mixture | 1.16055397609e-4 | 1.21976450020e-6 |
| 64×128 numerical mixture | 1.18435354844e-4 | 1.24516462612e-6 |
| 96×192 numerical mixture | 1.18864347847e-4 | 1.24909182089e-6 |

These adaptive sample minima differ appropriately from the 257-coordinate trajectory minima: the coordinate sets are distinct and are not compared pointwise. The original and lower ordered decks each contain exactly 257 rows, including the closing endpoint; their exact union plus six common probes contains 326 unique coordinates. Every physical numerical input remains `nspec=2`.

The build uses GNU Fortran 11.4.0, Open MPI 4.1.2, LAPACK/BLAS and four ranks, with single-threaded BLAS/OpenMP. Optimized flags are `-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=overflow`; debug flags are `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`. Linked production real Bessel symbols resolve to intrinsic `j0`, `j1`, `jn`. Revision is `99a087271928846cb8dbc49507d5a3a0f4d68d84`; the dirty working tree is preserved and recorded rather than treated as a clean commit baseline.

Protected production hashes:

```text
0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488  src/ALPS_fns.f90
c38f7ddb276f98a209f180f9d14d81e9bd9e85b6d4ac465b05cb93ff7904c2d6  src/ALPS_fns_rel.f90
acb9e48ab671964ad537c73651d413b1e6e34f2af43ab3174545de6a5109097b  src/ALPS.f90
```

All frozen-source, live-production and historical-development preservation manifests pass. The current validation files are separately compared with the executed frozen copies. Documentation and evidence packaging were completed after the acceptance execution; they do not alter the executed driver.

At 64×128 the optimized cached 326-frequency job spends 14.1401 seconds in collective determinants (about 43.4 ms/call); the 39-point uncached control spends 18.9959 seconds (487 ms/call). Debug equivalents are retained separately. These are measured costs on different frequency sets, not a matched performance-speedup claim. No optimization was attempted. The cached rank-one/rank-two peak weight-storage counters are each 127,124 bytes, and the recorded weight construction/work counts remain active. No extra per-frequency Bessel evaluation or derivative reconstruction is added.

## Retained evidence and disposition

`tests/data/nyquist_c2f4b/results/` retains the plain gate, root, contour, determinant-error, six-component convergence, density-weighting, numerical-linearity, derivative, CR and MPI summaries. `inputs.tar.gz` includes the exact generated component and pre-summed VDF tables, generator/physical inputs and shared frequency decks. `cases.tar.gz` retains native population tensors, full matched complex tensors/wave matrices/determinants, A1 samples and summaries, root qualification, assignment observations, work/timing metrics and standard output. To avoid redundant large dumps, actual interior derivative arrays are archived once per optimized mixture grid and endpoint in `derivatives.tar.gz`. `experimental_sources.tar.gz`, `math.tar.gz`, build logs and provenance preserve the executed variants and mathematical evidence. SHA-256 manifests and compressed-archive integrity checks accompany the evidence. Original C1–C2f.4a.1 reports and reference artifacts are unchanged.

**Recommendation:** accept this controlled validation of the existing coherent upper-half-plane method for one smooth, genuinely non-Maxwellian tabulated proton mixture. Independent analytic additivity, numerical superposition, six-component and root convergence, accepted inventories, arithmetic checks and MPI behavior agree. This is not production adoption, a new integration method, or authorization to broaden the study. Stop here as requested; the production C2b default and all physical tolerances remain intact.
