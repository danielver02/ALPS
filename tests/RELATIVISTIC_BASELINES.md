# Relativistic Chebyshev preparation: Tasks 0 and 1

## Scope and review boundaries

Audit baseline and current HEAD: `465a8afd7a0cd1fa6900ae6fd409cb417cf1cac4`.
The audit in the preceding conversation was reviewed against this source. Its
coordinate mismatch, original-grid polynomial allocation, and relativistic
fixed-u gamma stencil remain present.

Task 0 consists of the guard in `src/ALPS_io.f90` and the independent
`test_input_validation.f90` / `run_input_validation.sh` regression. Task 1 is the
remaining test-only Fortran/shell infrastructure, this report, the reference
data, and the optional entry in `run_test_suite.sh`. No commit was made.

Pre-existing modifications, including `distribution/distribution_analyt.f90`,
were preserved. The SHA-256 of the pre-existing tracked diff was
`ba431f983f0cd7e2ede547547012b356cb490a2b941cc6e3bf11680ea29a1cd5`;
excluding the Task 0 file and the Task 1 suite entry reproduces this digest.
The baseline cases do not invoke the user's analytic-distribution function.

## Task 0: species-specific input rejection

`init_param` checks `relativistic(is).and.ACmethod(is)==2` immediately after
`spec_read(is)`. This is before relativistic setup, fit allocation, distribution
input, and coefficient construction. `ACmethod` is an array indexed by species.
The guard also rejects that explicit selection when a proxy is selected.

The diagnostic identifies species and `AC_method=2`; error 18 explains the
original `(pperp,ppar)` versus required `(Gamma,pparbar)` rows and says that
relativistic Chebyshev support is under development. It uses the existing
`alps_error` / `MPI_ABORT` path and writes both stdout and the run log.

The independent regression compiles only the current `ALPS_var`, `ALPS_io`, and
a Fortran fixture/driver. Five four-rank cases test rejection at species 1 and
species 2, supported mixed selections, and relativistic types 4 and 5. Invalid
fixtures deliberately have no distribution tables: successful rejection must
precede any attempt to read them. Tests require a nonzero, non-timeout exit,
the species-specific log, and no bounds/runtime-error diagnostic. Supported
cases must complete actual `init_param` and retain species selections.

## Build and run

From the repository root:

```sh
bash tests/run_input_validation.sh
bash tests/run_test_suite.sh --relativistic-baselines --check
bash tests/run_numerical_baselines.sh --convergence
ALPS_TEST_NP=8 bash tests/run_numerical_baselines.sh --check
```

Requirements: Bash, a GNU-compatible MPI Fortran compiler, MPI launcher,
LAPACK/BLAS (default `-lopenblas`), `rg`, `timeout`, and standard Linux shell
utilities. `pkg-config` version reporting is optional. No Python is used.

Each invocation creates and prints a fresh `mktemp` directory outside the
source tree. `TMPDIR` inside the repository is rejected. Artifacts are retained
for inspection; re-running needs no cleanup and never reuses another run's
outputs. Sources and the two input templates are copied into the scratch
snapshot. Existing `.o`, `.mod`, executable, distribution, and cache files are
excluded. Compilation and all execution occur inside scratch.

Both the real `ALPS` main program and a test driver are freshly linked. The
test driver follows the production initialization order, calls production
`derivative_f0`, `determine_param_fit`, MPI distribution, harmonic setup, and
`disp`, and observes `chi0` and the returned determinant. It then attempts a
secondary root with the existing `rtsec`. Production code has no test hooks.

The existing distribution generator is compiled from source and run against
`test_relativistic_dist.in` and `test_bimax_dist.in` in scratch. The Fortran
fixture derives runs from `tests/test_relativistic.in` / `test_chebyshev.in`.
Generated distributions, logs, fit parameters, relativistic diagnostics, and
solution directories are isolated. Cache mode is off for these comparisons;
no existing cache is read or overwritten. Later cache-enabled runs can use
the same isolated working-directory convention.

Overrides: `MPIFC`, `MPIEXEC`, `MPIEXEC_FLAGS`, `FFLAGS`, `ALPS_TEST_LIBS`,
`ALPS_TEST_NP`, `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and `TMPDIR`.
Compiler/launcher variables are executable paths; flag variables are split
into arguments without shell evaluation. Default build flags are
`-O2 -g -fcheck=bounds -fbacktrace -ffree-line-length-none -ffpe-trap=overflow`.
Default MPI flags are `--oversubscribe`; default BLAS/OMP thread count is 1
unless already specified. Local MPI process communication must be available.

`provenance.txt` records commit, working-tree state, compiler and MPI versions,
wrapper expansion, flags, thread settings, library paths, and binary hashes.
`source.sha256` records the consumed snapshot, and `source.patch` records the
tracked diff. These distinguish a dirty source build from a historical binary.

`--capture` writes candidate observations only in scratch. It never updates the
reviewed reference files automatically. `--check` compares all records against
`tests/data/relativistic_baselines`. `--convergence` captures the baseline cases
plus joint gamma/u refinements and resonance-subinterval refinements, producing
`convergence.txt`. Convergence output is diagnostic, not an accuracy pass.

## Reference cases and data format

| Case | Source table / continuation | Parameters |
|---|---|---|
| `rel4` | Existing pair-plasma Juttner example; functional type 4 | 30x60 source, bounds 5, vA=1, equal masses, Theta=0.5 |
| `rel5` | Test-only `C exp(-2 Gamma-0.15 u^2)`; functional type 5 | Same source coordinates; source-grid normalization |
| `nr_fd` | Existing electron/proton bi-Maxwellian example; Chebyshev | 150x300 source, order 30, derivative_type=0 |
| `nr_analytic` | Same | derivative_type=1 |

Relativistic base target grids are 100x200; the convergence grids are 200x400
and 500x500. The original example uses 500x500. Base resonance subinterval
count is 100, compared separately with 200 and 400. Other existing integration
controls are retained. Scans and screen fit-output diagnostics are disabled.

Each record is `label Re(value) Im(value) absolute_tolerance relative_tolerance`
at round-trip double precision. Labels include species and sample/frequency
index. All six independent susceptibility entries for each species are stored
at each of three fixed complex frequencies, along with reconstruction,
derivatives, moments, and determinant: 66 records per relativistic case and 56
per nonrelativistic case, 244 total. Missing, reordered, nonfinite, and changed
records fail comparison; adversarial test files exercise that failure path.

Frequency samples are `(0.062713,-0.0001)`, its conjugate, and `(1,-0.0001)` for
relativistic runs; `(0.03,-0.00001)`, its conjugate, and `(0.045,-0.019)` for
nonrelativistic runs. Both use kperp=0.001; kpar=0.1 and 0.03 respectively.

`chi` records are exactly production `chi0`, including ALPS normalization by
omega and vA. `determinant` is the raw return from `disp` using its scaled wave
tensor, not a newly rescaled determinant. Momentum units are m_ref*vA_ref;
u=ppar/(m_s*c), frequencies are in the reference gyrofrequency, and wavevectors
are in inverse reference inertial length. Analytic amplitudes in relativistic
shape comparisons are normalized with the same finite transformed grid and
mask, not an assumed infinite-domain density.

Capture configuration: GNU Fortran 11.4.0, Open MPI 4.1.2, OpenBLAS 0.3.20,
one BLAS/OMP thread, four ranks, flags above. Consumed production/input hashes
are in `data/relativistic_baselines/source.sha256`.

Regression tolerance is `1e-12 + 1e-10*abs(reference)` for ordinary records and
`1e-24 + 1e-10*abs(reference)` for raw determinants. These bound floating-point
regression differences, not discretization/physics error. Four-rank fresh-build
repetition was exact for all 244 records. Eight-rank verification also reproduced
all 244 records exactly. Tolerances are deliberately much smaller
than the measured grid errors below; those errors are not hidden by widening
regression tolerances.

## Measured numerical results and unresolved accuracy

Relative distribution errors use `max(abs(reference),1e-12*peak)` in the
denominator. Gamma-grid derivative error is absolute error divided by twice
the completed distribution peak, compared against the exact centered stencil
of the analytic exponential, excluding intentional outer-edge placeholders.

| Metric, base grid | Measured maximum |
|---|---:|
| Isotropic interpolation vs analytic shape | 5.93e-14 |
| Isotropic functional fit vs completed grid | 3.19e-14 |
| Isotropic grid gamma derivative vs exact stencil, scaled | 4.10e-14 |
| Anisotropic interpolation, within original momentum bounds | 2.26e-3 |
| Anisotropic interpolation, including transformed extrapolation | 2.95e-1 |
| Anisotropic functional fit vs completed grid | 7.22e-2 |
| Anisotropic grid gamma derivative vs analytic stencil, scaled | 2.67e-4 |
| Nonrelativistic reconstruction, species 1 / 2 | 4.00e-14 / 3.12e-14 |

The anisotropic maxima include tails; the in-domain/all-domain split exposes
the existing spline extrapolation beyond the original perpendicular bound.
No interpolation or fit algorithm was adjusted. Isotropic affine-log behavior
and the nonrelativistic analytic Maxwellian value/derivative are independently
asserted, in addition to golden-value comparison.

The original generated source density is 0.9955210832998 in the Juttner case
and approximately 0.9997332906504 in the bi-Maxwellian cases. Relativistic
construction renormalizes its completed table to 1 (tested within 2e-12).
Original-grid integration does not renormalize its table. Finite grid/domain
quadrature and the differing normalization paths are preserved.

For six-entry, two-species susceptibility comparisons, `compare_convergence`
reports the Frobenius difference, counting off-diagonal entries twice, divided
by max(1, refined tensor norm). At frequency sample 1:

| Change | rel4 | rel5 |
|---|---:|---:|
| Grid 100x200 -> 200x400 | 0.13745 | 0.14247 |
| Grid 200x400 -> 500x500 | 0.08855 | 0.09085 |
| Resonance intervals 100 -> 200, fixed base grid | 2.95e-5 | 2.31e-5 |
| Resonance intervals 100 -> 400, fixed base grid | 7.32e-5 | 6.31e-5 |

All three frequency results, determinant differences, and other refinement
comparisons are preserved in `data/relativistic_baselines/convergence.txt`.
The observed grid changes are nonmonotone and substantial. Even the isotropic
fit is essentially exact, so fit error alone cannot explain that sensitivity.
These are reproducible implementation baselines, not an independent converged
dielectric reference. No empirical physical-accuracy tolerance is justified
yet, and an independent analytic-derivative quadrature/reference is still
needed. Source-grid, separate gamma/u, Tlim, harmonic, and Bessel convergence
are not established by this joint-grid/subinterval study.

Root attempts are retained separately in `root_status.dat`; they are not used
to mask fixed-frequency changes or automatically become golden roots. The
base-grid relativistic attempts stop with scaled residuals 1.91e-5 (type 4)
and 3.61e-4 (type 5), despite small-step termination, and are marked unresolved.
At 500x500 the returned roots are approximately
`0.03962301059 - 3.12562069e-6 i` (type 4, scaled residual 7.82e-11) and
`0.04573855005 - 3.04974111e-6 i` (type 5, 7.51e-18). They are algebraically
converged; their branch identity and integration accuracy remain unqualified.
The nonrelativistic attempts approach a very low-frequency root rather than
the intended initial branch; they likewise do not supply a trusted physical
root reference. No root-finding or dielectric code was changed to resolve this.

## Coordinate and derivative contracts

`relativistic_references.f90` provides test-only real-Gamma / complex-u values
and both partial derivatives for isotropic Juttner, anisotropic exponential,
and boosted thermal distributions. C is explicit. Tests cover Theta=0.01,
0.1,1,10, real and complex u, and drift beta=0.1,0.3,0.5. Independent amplitude/
phase evaluation checks the complex exponential.

With a_s=vA/ms, tests check r=a_s*Pperp, u=a_s*Ppar,
Gamma=sqrt(1+r^2+u^2), the inverse map, and
`Pperp*abs(d(Pperp,Ppar)/d(Gamma,u)) = Gamma/a_s^3`.
The factor 2*pi is identical on both sides. Non-unit vA=0.6 and mass=2 avoid
hiding scale errors. The same measure is used in `construct_f0_rel`.

Measured maximum analytic complex-value relative error is below 3e-14;
centered-stencil factor error is below 4e-13. Momentum round-trip error is zero
for the chosen point; finite-difference Jacobian relative error is 6.36e-11.
Thresholds are 3e-13 for analytic values, 2e-12/3e-12 for centered/forward
stencil factors, 2e-14 for the round trip, and 2e-10 for the numerically
differentiated Jacobian. The latter includes O(h^2) error with h=1e-5.

Gamma differences evaluate the same complex u on both rows. The exact decay
factors are sinh(lambda*h)/(lambda*h) for a centered difference and
(1-exp(-lambda*h))/(lambda*h) for a forward difference. Halving h verifies
second- and first-order convergence. Parallel differences verify convergence
to the supplied derivative, including the identically zero isotropic result.

A dedicated fixture uses Gamma=1.05, h=0.01, Re(u)=0.3: the point lies inside
the current row's cone and outside the lower neighbor's cone. The analytic
function remains defined there. Another check distinguishes the required
partial derivative from the total derivative along u(Gamma), which adds
u'(Gamma)*F_u. Production derivative calculations are unchanged.

## Sparse gamma rows

Each relativistic run produces all rows in `sparse_rows.dat`: index, Gamma,
sample count, sampled interval, active residue flag, and neighbor extrapolation
flag. `sparse_summary.dat` summarizes counts. Only actual production
`determine_resonances` selections and the per-row real-projection cone test
contribute to active counts. Species 1 is reported; both fixture species have
identical geometry. Frequency is 0.062713-0.0001i, harmonics -1,0,1; only n=0
passes the global resonance dispatcher for this example.

At 100x200 there are no empty rows, one singleton, and five rows with fewer
than 31 samples. Ninety-five rows have active residues; rows 5 and 6 require
neighbor evaluations outside sampled intervals.

At the original 500x500 resolution, rows 0--4 have sample counts
1,15,23,27,31. Their intervals are [0,0], [-0.14,0.14], [-0.22,0.22],
[-0.26,0.26], [-0.30,0.30]. There are no empty rows, one singleton, and four
rows with fewer than 31 samples. The first active residue row is 24
(Gamma=1.29478856457); 476 rows are active and row 24 needs extrapolation into
its lower neighbor. All integrated gamma rows are visited by the Landau loop;
the cone/dispatcher determines actual residues. Row 0 is not used by its
gamma derivative; the top row is needed as a neighbor of ngamma-1.

Degree 30 cannot be fitted unconstrained to 15 samples. Requested degree must
eventually be a maximum, and neighboring-row evaluations cannot be zeroed at
real-data boundaries. No adaptive fitting or continuation status mechanism
was implemented here. The accepted future design decisions remain adaptive
degree, rank-aware SVD (preferably DGELSD), and explicit evaluation status.

## Next step

The review inventory is:

| Files | Change |
|---|---|
| `src/ALPS_io.f90` | Task 0 guard in `init_param` and diagnostic case 18 in `alps_error` |
| `tests/test_input_validation.f90`, `run_input_validation.sh` | Independent Task 0 input-reader/MPI regression |
| `tests/relativistic_references.f90`, `test_coordinate_contract.f90` | Test-only analytic functions, coordinate and derivative contracts |
| `tests/make_baseline_case.f90`, `test_numerical_baselines.f90` | Fortran fixtures and observations of actual production routines |
| `tests/compare_baselines.f90`, `test_comparison_fixtures.f90` | Numerical regression reader and adversarial input checks |
| `tests/compare_convergence.f90` | Grid/integration sensitivity report |
| `tests/run_numerical_baselines.sh`, `run_test_suite.sh` | Isolated source build/run and optional established-suite entry |
| `tests/data/relativistic_baselines/` | Four numerical snapshots, convergence results, consumed source/input hashes |
| `tests/RELATIVISTIC_BASELINES.md` | Reproduction, measured results, and limitations |

The tests and isolated builds are ready to support a separately approved
generalized polynomial kernel task. Physical end-to-end acceptance still needs
an independent relativistic response reference and resolution of the observed
integration/root limitations. Decisions about sparse-row reliability thresholds,
nonpositive data, and the continuation acceptance region remain for maintainers.
Tasks 0 and 1 make no change to relativistic fitting, interpolation/cache format,
contour geometry, Bessel evaluation, or nonrelativistic polynomial arithmetic.
