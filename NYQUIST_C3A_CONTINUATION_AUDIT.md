# C3a — Causal continuation of the coherent numerical susceptibility

**Mathematical/source audit, 2026-11-25. No implementation or physical calculation was performed.** Only this report is added. Production and experimental sources, prior reports, fixtures, tolerances and defaults are preserved.

**Recommendation:** retain C2d/C2f.3a as the validated upper-half-plane reference. For a causal response over a useful damped-frequency domain, develop strategy **B: one analytic momentum representation used consistently in both the coherent integral/moments and its residue**. Exact continuation of the current PWL model (A) is a valuable local reference, but is generally multivalued around its grid nodes and need not converge to the physical damped response at fixed damping. A PWL integral plus the existing fitted residue (C) is an approximation with a measurable boundary defect, not an exact continuation of that upper-half-plane model. The first proposed implementation is a small standalone analytic-polynomial/branch reference, not an ALPS backend or damped solver; its scope is specified below. **C3b has not begun.**

## 1. Evidence, provenance and terminology

Reviewed: `validation/c2d/FORMULATION.md`, `NYQUIST_C2D_RESULTS.md`, `NYQUIST_C2F3A_RESULTS.md`, `NYQUIST_C2F4B_RESULTS.md`, their coherent builders/helpers, the C2c kernel, and the current nonrelativistic production call chain. Historical C2d comments about the original Bessel evaluation and `int_ee` corner describe their then-current baseline; both were subsequently corrected by C2f.1. They are not outstanding corrections proposed here.

The audited checkout is revision `99a087271928846cb8dbc49507d5a3a0f4d68d84` with its existing working-tree changes. Relevant current SHA-256 identities are:

```text
acb9e48ab671964ad537c73651d413b1e6e34f2af43ab3174545de6a5109097b  src/ALPS.f90
0c3dae67819cde17657a7376715d2ad607e30cf4a06f22a4a27378c6e1191488  src/ALPS_fns.f90
2651e19bc7574d692cc9303786a387630f3429264ebe3ff871b572c043e78da2  src/ALPS_analyt.f90
989cfbb8fd3a8252e195d868a93a09552cfbce3e4f1172f93c6a127e947e5cc8  src/ALPS_com.f90
```

Source locations below refer to this checkout. “Direct” means integration over the fixed real momentum interval without a Landau contour deformation. “Causal” means analytic continuation from **Im(omega)>0**, not a fresh direct integral below the axis. “Residue” denotes the numerator coefficient of `1/(p-a)`; it includes tensor/species factors only when explicitly indicated. Upper/lower signs on Cauchy boundary values refer to **Im(a)**, not automatically Im(omega).

A before/after SHA-256 manifest check passed for all 2,916 existing files selected from `src/`, `validation/`, `tests/data/`, and the prior Nyquist/Bessel reports. No compiler, MPI job, numerical validation executable, root search or regression suite was run. The equations below are deductions; numerical values explicitly attributed to earlier reports are historical observations, not new measurements.

General complex-log conventions are those in [NIST DLMF §4.2(i)](https://dlmf.nist.gov/4.2#i): a tracked logarithm can differ from its principal value by integer multiples of `2*pi*i`, and opposite boundary values on the negative real axis differ by that amount. All ALPS-specific expressions and the node analysis below are derived here from the cited source, rather than inferred from a root count.

## 2. Source-level map

| Step | Exact source location | Actual representation/decision |
|---|---|---|
| Read numerical VDF | [read_f0](src/ALPS_io.f90#L1058), lines 1058–1123 | Reads `(pperp,ppar,f0)` at indices `0:nperp,0:npar` from the ordinary species table. Analytic species bypass it. |
| Construct real derivatives | [derivative_f0](src/ALPS_fns.f90#L103), lines 103–122 | Interior centered differences of the original tabulated f0, stored as `df0(:,:,:,1:2)`. Usable derivative nodes exclude both outer endpoints. |
| Initialize continuation once | [ALPS startup](src/ALPS.f90#L98), lines 98–138 | Reads f0, constructs derivatives, calls native NR Chebyshev initialization and the applicable other fits while f0 exists, releases f0, broadcasts data/models, determines harmonic limits/assignments and Bessel arrays. |
| Transfer existing representations | [pass_distribution](src/ALPS_com.f90#L367), especially 430–439; [complete_nr_cheb_initialization](src/ALPS_com.f90#L32), 32–105 | Sends real derivative arrays, nodes and legacy fit parameters; separately exports/imports native NR Chebyshev coefficients/metadata and enables evaluation. No per-frequency fitting. |
| Harmonic susceptibility | [disp](src/ALPS_fns.f90#L337), lines 337–505 | Analytic species use NHDS. Numerical workers integrate their assigned n, including negative-n partners, six tensor modes and one nonresonant zz term. |
| Resonance presence | [determine_resonances](src/ALPS_fns.f90#L691), lines 691–729 | Uses `real(a)` relative to usable parallel support and a nearby exterior margin. This is a real-frequency-dependent selection, not a complex-analytic condition. |
| Production integration/continuation | [full_integrate](src/ALPS_fns.f90#L737), lines 737–779 | If no resonance is found, ordinary `integrate`, with **no Landau addition**. Otherwise NR Im(omega)>0 uses `integrate_res`; Im(omega)=0 adds one `landau_integrate`; Im(omega)<0 adds two. |
| Production direct quadratures | [integrate](src/ALPS_fns.f90#L786), 786–851; [integrate_res](src/ALPS_fns.f90#L857), 857–1190 | Nodal trapezoidal background; moving core/tails, edge returns, and `abs(Im(a))>Tlim` regular versus small-imaginary treatment. C2b is the default. |
| Production interpolated numerator | [funct_g](src/ALPS_fns.f90#L1227), 1227–1269 | Continuous PWL interpolation of the **complete frequency-dependent nodal numerator**, nearest usable cell extrapolation outside the nodes. `funct_g_legacy`, 1272–1379, is preserved separately. |
| Continuation value/derivatives | [eval_fit](src/ALPS_analyt.f90#L1705), 1705–1926; [landau_integrate](src/ALPS_fns.f90#L1385), 1385–1558 | Evaluates a separately represented f0 at complex momentum, then obtains perpendicular differences and the selected parallel derivative. This is not evaluation of `df0`'s PWL interpolants. |
| Tensor numerator/prefactors | [resU](src/ALPS_fns.f90#L1666), 1666–1702; [int_T](src/ALPS_fns.f90#L1707), 1707–1806; [int_T_res](src/ALPS_fns.f90#L1810), 1810–1908 | NR denominator `m*omega-k*p-n*q`; real Bessel coefficients; p powers are evaluated at complex a for the residue, without complex conjugation. |
| Nonresonant term | [int_ee](src/ALPS_fns.f90#L1563), 1563–1661 | Frequency-independent NR longitudinal quadrature, with the C2f.1 corner correction. It has no Cauchy pole/residue of its own. |
| Collective assembly | [disp](src/ALPS_fns.f90#L499), 499–619 | Multiplies each numerical species contribution by `n_s*q_s`, MPI-reduces tensors, adds vacuum once, forms wave matrix and raw determinant, broadcasts D and barriers. Conventional species chi is raw chi divided by `(omega*vA)^2`, lines 535–536. |

### 2.1 Production C2b versus experimental coherent C2d

C2b does not compute an exact single fixed-grid Cauchy transform. Its background interpolates/integrates the nodal quotient, while its core uses the interpolated complete numerator. The core width and tail ownership use `real(a)` and grid-cell choices ([lines 952–1024](src/ALPS_fns.f90#L952)); near-edge returns precede the core ([991–1000](src/ALPS_fns.f90#L991)). The regular/small-imaginary switch is at [1038](src/ALPS_fns.f90#L1038). The small branch estimates a numerator slope and uses signed limiting `i*pi*G(real(a))` for nonzero Im(a), with no such term at exactly zero ([1100–1186](src/ALPS_fns.f90#L1100)). Thus its Landau prescription is a physical/numerical approximation; selection by damping sign alone does not prove analyticity of the entire approximate backend. This audit neither changes nor requalifies those quadratures.

The C2d substitution is explicit in [variant.awk](validation/c2d/variant.awk#L1), lines 1–16: **only nonrelativistic Im(omega)>0** calls `integrate_coherent` and uses `coherent_ee`. It bypasses `found_res` there; real and damped queries retain C2b. This composite experimental executable is therefore **not claimed to be one holomorphic susceptibility across the real axis**.

The coherent formula is in [helpers_cached.inc](validation/c2d/helpers_cached.inc#L1), lines 1–64, and [model.f90](validation/c2d/model.f90#L6), lines 6–39. It uses PWL interpolants of **A=df0/dpperp and B=df0/dppar**, then integrates their p and p² factors exactly. This differs from C2b's PWL interpolation of the complete numerator. Its parallel support is stored nodes `1:npar-1`; its perpendicular rule is `2*pi*du` times weights one at `i=1,...,nperp-2` and one-half at `i=nperp-1`.

The current path from [C2f.4b build](validation/c2f4b/build.sh#L5) through [C2f.4a](validation/c2f4a/build.sh#L13), [C2f.3](validation/c2f3/build.sh#L18) and [C2f.2](validation/c2f2/build.sh#L21) reuses that formulation. [C2f.3a build_wide.sh](validation/c2f3a/build_wide.sh#L5) and [wide.awk](validation/c2f3a/wide.awk#L1) promote only the Cauchy kernel/intermediates/weight storage to `selected_real_kind(18)`, retaining double APIs and cache keys. This is restricted extended precision, not a claim that all tensor arithmetic is real128. Its independent real128 accuracy references are described in C2f.3a; no arithmetic study is repeated here.

Cached weights implement the same linear functional; [work.f90](validation/c2d/work.f90#L13), lines 13–36, invalidate on frequency, species, grid, k or storage dimensions. They currently encode the **direct** sheet only. A future multivalued continuation would require an explicit sheet/domain contract; merely reusing frequency-only direct weights cannot silently encode path information.

### 2.2 What the existing Chebyshev machinery represents

For eligible NR numerical species, `AC_method=2`, `kind=1`, [nr_cheb_species_eligible](src/ALPS_analyt.f90#L679), 679–705, selects native row models. [initialize_nr_cheb_svd](src/ALPS_analyt.f90#L740), especially 769–827, fits **f0 itself**, separately at each perpendicular row on the original parallel nodes, not A, B, D_n or the tensor numerator. The direct representation is a finite Chebyshev polynomial

\[
\widehat f_i(p)=P_i(x),\quad x=(p-c_i)/s_i,
\]

or, for log fits,

\[
\widehat f_i(p)=10^{P_i(x)}=\exp[\ln(10)P_i(x)].
\]

The former is polynomial and the latter entire in p for fixed finite coefficients and nonzero real scale. Fitting is rectangular SVD on actual nodes, with rank/degree diagnostics ([964–1024](src/ALPS_analyt.f90#L964)); LOG10 operates on positive samples and rejects interior zeros, excluding only exact exterior zero samples, not inventing a physical support cutoff ([cheb_fit, 517–592](src/ALPS_analyt.f90#L517)). These are finite-fit models, not proof of analytic accuracy away from the real data.

[cheb_polynomial_complex](src/ALPS_analyt.f90#L1487), 1487–1540, uses complex affine mapping/Clenshaw with no unit-disk or momentum-support clipping. [cheb_evaluate_complex](src/ALPS_analyt.f90#L1269), 1269–1324, and [cheb_reconstruct](src/ALPS_analyt.f90#L1327), 1327–1466, reconstruct the physical f and parallel derivative, including the `ln(10)*P'` chain rule. Range/phase/underflow statuses limit reliable floating-point evaluation; these are not physical domain boundaries. Status-less native evaluation failures call the existing continuation error handler ([eval_fit 1780–1787](src/ALPS_analyt.f90#L1780)). Missing method-2 models do not fall back to old coefficient storage ([1910–1924](src/ALPS_analyt.f90#L1910)).

In `landau_integrate`, interior

\[
\widehat A_i(a)=\frac{\widehat f_{i+1}(a)-\widehat f_{i-1}(a)}{2\Delta u},
\]

and either

\[
\widehat B_i(a)=\widehat f'_i(a)\quad(\texttt{derivative_type=1}),
\qquad
\widehat B_i(a)=\frac{\widehat f_i(a+h)-\widehat f_i(a-h)}{2h}\quad(\texttt{type=0}).
\]

See [1461–1505](src/ALPS_fns.f90#L1461). Native type one skips unused parallel FD queries; perpendicular FD queries remain. The outer perpendicular rows use one-sided perpendicular differences ([1514–1546](src/ALPS_fns.f90#L1514)). There is no general equality between these analytic functions and the PWL A/B of C2d, even if a fit reproduced f at every tabulated node. Type zero might reproduce interior nodal differences in that special case; type one need not even do that. Matching nodal values alone does not match complex values or a PWL function between nodes.

The other existing methods remain distinct: method zero calls `distribution_analyt` at complex parallel p ([eval_fit 1796–1818](src/ALPS_analyt.f90#L1796)); method one evaluates configured fitted analytic functions ([1821–1908](src/ALPS_analyt.f90#L1821)), whose powers/square roots may have their own branch domains. Methods zero/one preserve a zero-neighbor early return in `landau_integrate` ([1487–1495](src/ALPS_fns.f90#L1487)); native method two instead reports failures and allows successfully evaluated zeros. A Chebyshev-specific analyticity statement must not be generalized to every method-one formula or numerical failure policy.

## 3. Coherent tensor and Cauchy conventions

Suppress the species label temporarily. Let k=kparallel be real and nonzero, m>0, q signed, u=pperp, p=ppar,

\[
a_n=\frac{m\omega-nq}{k},\qquad d=m\omega-kp-nq=-k(p-a_n),
\quad A=\partial_u f,\ B=\partial_p f,\ D_n=nqA+kuB.
\]

From actual `resU`, pointwise

\[
\frac{q[\omega A+(k/m)(uB-pA)]}{d}
=\frac qm\left(A+\frac{D_n}{d}\right).
\tag{1}
\]

Write `int_T=T_nu(p)=K_nu(u)*p^ell`, with ell=0 for xx,yy,xy; ell=1 for xz,yz; ell=2 for zz. For `z0=kperp/q` and `J=J_n(z0*u)`, the finite-kperp prefactors are

| nu | ell | K_nu |
|---|---:|---|
| xx | 0 | n² J²/z0² |
| yy | 0 | u² (J')² |
| zz | 2 | J² |
| xy | 0 | i u n J J'/z0 |
| xz | 1 | n J²/z0 |
| yz | 1 | -i u J J' |

Exact parallel limits are the special cases in `int_T`, not divisions by zero in this table. The Bessel factors are frequency independent. The i factors and antisymmetric off-diagonal signs must be retained; there is no Hermitian conjugation of complex p in the continued tensor.

For fixed usable support [P0,PN], define

\[
M_\ell(F)=\int_{P_0}^{P_N}p^\ell F(p)\,dp,
\quad C_\ell(F;a)=\int_{P_0}^{P_N}\frac{p^\ell F(p)}{p-a}\,dp.
\]

The per-node coherent integral is

\[
I_{n\nu}=\frac{qK_{n\nu}}m\left[M_\ell(A)-\frac1k C_\ell(D_n;a_n)\right].
\tag{2}
\]

The perpendicular functional is \(\mathcal Q[X]=2\pi\Delta u\sum_i w_i X_i\) with the coherent weights above; raw species tensors acquire the further **n_s*q_s** in `disp`. No density or charge factor may be added a second time.

The nonresonant moment is \(E_i=(q/m)[uM_1(B)-M_2(A)]\). For n=0,zz, coherent integration combines J0²E with (2) **before evaluation**:

\[
I_{0,zz}^{\rm combined}=-\frac{q\omega uJ_0^2}{k}C_1(B;m\omega/k),
\quad E_i^{\rm complementary}=(1-J_0^2)E_i.
\tag{3}
\]

This preserves the explicit omega factor and does not require an infinite harmonic identity. Any future representation must use the same A/B and moments in both (2) and (3). A residue-only change cannot repair an independently inconsistent moment discretization. No new residue belongs to E itself.

## 4. Exact cell primitive and the direct sheet

On one cell [L,U], h=U-L>0 and F(p)=alpha+beta*p, polynomial division gives

\[
\frac{\alpha+\beta p}{p-a}=\beta+\frac{\alpha+\beta a}{p-a},
\]

hence

\[
C_{0,j}=\beta h+(\alpha+\beta a)\Lambda_j(a),
\quad \Lambda_j(a)=\int_L^U\frac{dp}{p-a}
=\operatorname{Log}_{\Gamma}(U-a)-\operatorname{Log}_{\Gamma}(L-a).
\tag{4}
\]

Writing `log((U-a)/(L-a))` is shorthand for this **specified branch**, not permission to discard a multiple of `2*pi*i` on a continued sheet. For a nonreal a and a straight increasing real cell, U-a and L-a lie in the same open half-plane. Their principal-argument difference lies in (-pi,pi); thus the two principal logs, and the principal log of their ratio, **do give the correct direct integral under these conditions**. This justifies the current kernel; the warning concerns continuation to another sheet, real-axis limits/signed zeros, or differently tracked paths, not a defect in the direct formula.

For Im(a)>0, the difference in arguments is positive; for Im(a)<0 it is negative. For a real interior r,

\[
\Lambda_j(r+i0)=\ln\frac{U-r}{r-L}+i\pi,
\qquad
\Lambda_j(r-i0)=\ln\frac{U-r}{r-L}-i\pi.
\tag{5}
\]

Outside the cell the direct upper/lower limits agree. Equation (5) follows directly from the endpoint logarithms: L-r approaches the negative real axis from below/above, respectively. Equivalently, for t=p-r the imaginary part of `1/(t-i*epsilon)` integrates to +pi as epsilon decreases to zero. This fixes the Plemelj sign for ALPS's denominator convention without relying on a remembered sign rule.

The actual [cell_cauchy](validation/c2c/cauchy.f90#L41), lines 41–78, uses this endpoint-log difference, with a small-ratio series to avoid far-field cancellation; both are approximations to the same analytic direct function. [domain_cauchy](validation/c2c/cauchy.f90#L79), 79–99, subtracts the common resonant numerator to obtain a finite **principal value** at real interior nodes. That exact-real branch is not already a causal boundary value. Individual cell endpoint evaluations with a nonzero endpoint numerator are singular; their cancellation in an interior-node domain PV must not be replaced by evaluating divergent cells separately. Extended arithmetic preserves these mathematical definitions; it does not choose a causal sheet.

## 5. Boundary jumps, higher moments and signed k

For continuous PWL F, including its internal nodes where the PV exists,

\[
C_\ell^\pm(F;r)=\operatorname{PV}\int\frac{p^\ell F(p)}{p-r}\,dp
\ \pm i\pi r^\ell F(r),
\qquad C_\ell^+-C_\ell^-=2\pi i\,r^\ell F(r).
\tag{6}
\]

Here a superscript denotes approach in the a-plane. For a **specified cell j** crossed away from its endpoints, F_j(a)=alpha_j+beta_j*a is an entire polynomial. Local continuation from Im(a)>0 through that cell is

\[
C_{\ell,j}^{\rm continued}(a)=C_\ell^{\rm direct}(F;a)+2\pi i\,a^\ell F_j(a)
\quad(\operatorname{Im}a<0).
\tag{7}
\]

This refers to the full-domain direct integral plus the one cell's crossed residue, not a separate change in all un-crossed cells. For continuation from Im(a)<0 to Im(a)>0, the sign is reversed.

Polynomial division yields, on every consistently chosen sheet,

\[
C_1(F;a)=M_0(F)+aC_0(F;a),
\quad C_2(F;a)=M_1(F)+aM_0(F)+a^2C_0(F;a).
\tag{8}
\]

Thus their corrections in (7) are respectively `2*pi*i*a*F_j(a)` and `2*pi*i*a**2*F_j(a)`; the moments have no jump. More generally \(C_\ell=M_{\ell-1}+aC_{\ell-1}\).

Let sigma=sign(k). Since `Im(a)=m*Im(omega)/k`, the physical starting half-plane is Im(a) of sign sigma. The causal real-frequency boundary and the local damped continuation are

\[
C_\ell^{\rm causal}(r)=C_\ell^{\rm PV}(r)+i\pi\sigma r^\ell F(r),
\]
\[
C_\ell^{\rm causal}(a)=C_\ell^{\rm direct}(a)+2\pi i\sigma a^\ell F_j(a)
\quad(\operatorname{Im}\omega<0),
\tag{9}
\]

provided that the chosen continuation path crosses cell j and has no additional winding around branch points. If it crosses the real a-axis **outside the support**, the corresponding direct-sheet crossing has no jump instead. A later move in the damped half-plane must stay on that chosen sheet; it cannot re-decide the crossing from the current real resonance position.

For a general harmonic/component of (2), the full local damped correction, before perpendicular summation and the final n_s*q_s, is therefore

\[
\Delta I_{n\nu}=-\frac{2\pi i qK_{n\nu}}{m|k|}\,a_n^\ell D_{n,j}(a_n).
\tag{10}
\]

For the n=0 combined zz term it is equivalently

\[
\Delta I_{0,zz}=-\frac{2\pi i q\omega uJ_0^2}{|k|}\,a_0 B_j(a_0).
\tag{11}
\]

Using `k*a0=m*omega` makes (10) and (11) identical; no extra longitudinal correction should be added. These signs are compatible with the source's `abs(kpar)` in the residue while its direct core uses signed Im(a). They must be checked at both k signs, not by conjugating a computed response. Negative q remains signed in a_n and D_n, while the final species normalization remains n_s*q_s.

## 6. Interpolation nodes: continuity is not global analytic continuation

Write the full direct C0 as

\[
\sum_j\{\beta_j h_j+F_j(a)[\log(p_{j+1}-a)-\log(p_j-a)]\}.
\tag{12}
\]

At an interior node x=p_j the coefficient of `log(x-a)` from the adjacent cells is

\[
F_{j-1}(a)-F_j(a)
=(\beta_{j-1}-\beta_j)(a-x),
\tag{13}
\]

because F is continuous at x. Constant logarithmic singularities cancel. Unless the slopes also match, the remaining term is proportional to **(a-x)*log(x-a)**: its value tends to zero at the node, but its derivative contains a logarithmic singularity and its analytic continuation has nonzero monodromy. One counterclockwise circuit of a about x changes this term by

\[
2\pi i(\beta_{j-1}-\beta_j)(a-x).
\tag{14}
\]

For C_ell multiply the log coefficient by a^ell; polynomial moments cannot remove its monodromy. Special zeros or cancellations can weaken/remove a particular component's singularity. There is no general cancellation of every slope jump across all tensor components/harmonics. Continuous A/B therefore imply continuous Plemelj values at internal nodes, **not** a holomorphic extension through every node.

Two paths crossing different cells can reach the same lower-half-plane a on different sheets. Each local germ has unique analytic continuation along a specified nonsingular path; the multi-sheet continuation is not ambiguous once that path is specified. There is generally **no single-valued holomorphic domain containing the upper half-plane and continuing across all interpolation nodes and cell interiors at once**. A useful local domain crossing one cell (and avoiding the nodes) does exist. A fixed sheet can also be defined over a useful lower domain with appropriate cuts; it need not match the physical boundary germ at every other real cell.

The direct real-line Cauchy transform is itself a well-defined single-valued function off its integration segment. Reaching its lower value by going around a support endpoint is a valid direct-sheet path, but is generally **not** the path through the resonant interval that supplies Landau continuation. Neither the existence of that direct function nor node branch points make all continuation impossible.

A particularly unsafe shortcut is to evaluate (7) with `j=cell(real(a))` at every damped query. At `a=x-i*eta`, neighboring affine continuations differ by `(beta_left-beta_right)*(-i*eta)`, which is nonzero for eta>0. Consequently the residue correction has a finite value jump along that vertical switching line. Real-axis continuity of PWL interpolation does not prevent this off-axis jump. Tracking a sheet by the *previous solver query* is also insufficient: a root solver may visit frequencies out of order. The evaluated susceptibility must be a deterministic function on a declared sheet/domain, not depend on query history. Current caches have no such sheet key.

### 6.1 Why grid refinement alone does not settle damped convergence

For a smooth target F, PWL direct integrals converge on compact sets a fixed distance off the real support; slope jumps shrink under refinement. This is consistent with C2d upper-half-plane validation. It does **not** imply convergence of the local analytic continuation at a fixed negative imaginary displacement.

A concrete counterexample is F(p)=p². Its cell interpolant on [L,U] is

\[
F_j(a)=(L+U)a-LU.
\]

As cells shrink around a crossing point r, this tends to `2*r*a-r**2`, whereas the analytic target residue is a². At `a=r-i*eta`, their limiting difference is **eta²**, independent of the grid width. The direct PWL integral converges, but the continued residue leaves a `2*pi*i*eta**2` discrepancy for k>0. Thus a fixed-degree local PWL germ need not converge to the physical damped response at fixed eta. It can be controlled in a shrinking near-axis domain with curvature/depth bounds, or for globally affine data, but that is a separate approximation statement. Small node-switch jumps or a successful damped root solve would not remove this obstruction.

There is a stronger scalar restriction: on any fixed sheet, changing logarithm branches in (12) adds `2*pi*i` times a fixed integer linear combination of the affine F_j(a). The sheet difference for C0 is therefore affine in a, however fine the grid; for C_ell it is a^ell times that affine function. If the direct integrals converge locally uniformly, a locally uniform limit of those affine C0 sheet differences is still affine. It cannot reproduce the nonlinear analytic residue of a generic smooth F on a fixed open lower-half-plane domain. Choosing different sheets separately at each frequency evades that premise only by abandoning one analytic response on that domain. This does not preclude controlled near-axis or exceptional-component approximations, but it rules out a general grid-refinement argument for strategy A.

## 7. Finite-support endpoints and derivative boundary data

At the left endpoint the log coefficient is `-F_0(a)` and at the right endpoint it is `+F_last(a)`. A nonzero endpoint value gives a logarithmic divergence. If F vanishes at the endpoint but its first nonzero one-sided derivative does not, the value can remain finite while a weaker `(a-P)*log(P-a)` branch point remains. For PWL data, zero endpoint value alone does not remove that branch; both value and slope of that cell must vanish to eliminate that cell's endpoint contribution. Higher p powers can soften an origin endpoint singularity but need not remove it.

The relevant support here is **usable derivative support**, not the outermost f0 samples: C2d integrates from node 1 to node npar-1, whose A/B values are ordinary centered differences using outer f0 data. There is no imposed analytic taper there. The native fit uses f0 rows on 0:npar, potentially excluding exact exterior zeros for LOG10 fitting; its fit support and data quantity thus differ from C2d's integral support. Inaccurate endpoint derivatives change the residue/log coefficients and all direct moments. A jump introduced by extending a nonzero endpoint value by zero is an endpoint feature, not an internal interpolation discontinuity.

A nonzero globally analytic function cannot be identically zero outside a finite support interval and analytic through its endpoints. Even a globally analytic fitted numerator integrated only over a finite interval retains endpoint branch points. Strategy B can remove **interior grid-node** singularities while retaining a causal branch on an explicitly bounded domain away from the endpoints; it does not automatically produce a full-real-line physical susceptibility. No tail extension or clipping is authorized by this audit.

In frequency coordinates, branch points occur at

\[
\omega_{n,j}=\frac{k p_j+nq}{m}.
\tag{15}
\]

The full tensor has the union of the contributing species/harmonic branch loci, unless a coefficient cancels. Negative k reverses their ordering. A useful damped domain must account for every relevant harmonic, including those whose real resonance is outside the finite support. Reapplying `found_res` as a lower-plane switch is not a branch prescription.

## 8. Existing Landau term and representation mismatch

Set \(\widehat D_n(a)=nq\widehat A(a)+ku\widehat B(a)\). On resonance,

\[
[\omega-ka/m]\widehat A+(ku/m)\widehat B=\widehat D_n/m.
\]

From [landau_integrate 1507–1552](src/ALPS_fns.f90#L1507), the source returns the **half-jump**

\[
L_{n\nu}=i\pi\,\mathcal Q_L\left[-\frac{qK_{n\nu}}{m|k|}a_n^\ell\widehat D_n(a_n)\right].
\tag{16}
\]

The source multiplier is `i*dpperp*pi*2*pi`; equivalently it is `i*pi` times the usual `2*pi*dpperp` perpendicular measure. `full_integrate` adds L at a real frequency and 2L below the frequency axis. This matches the signs/factors in (9)–(10) **when integral and residue describe the same numerator, weights, support and continuation chart**. It does not establish those equalities for the actual approximations.

There is a concrete perpendicular coverage difference: the Landau loop at [1468–1474](src/ALPS_fns.f90#L1468) covers 1:nperp-1, with weight one at the first usable row and one-half at the last. It additionally includes row 0 with weight one-half and row nperp with weight one-half ([1514–1550](src/ALPS_fns.f90#L1514)). Coherent Q includes neither outer row. The loop's test `iperp==0` is unreachable inside its 1:nperp-1 loop; the surrounding comment is not an alternative weighting rule. Smooth-axis/tail factors can make endpoint contributions small or exactly zero in special components, but they must be measured/derived, not silently removed or assumed equal. No weighting correction is made here.

On a fixed local chart j, using identical Q for clarity, a hybrid lower value differs from the exact PWL continuation by

\[
\delta I_{n\nu}(a)=-\frac{2\pi i qK_{n\nu}}{m|k|}a_n^\ell
[\widehat D_n(a)-D_{n,j}(a)].
\tag{17}
\]

At a real interior r, the mismatch between the lower hybrid limit and the validated growing boundary is the same expression with a=r, plus any Q_L-versus-Q contribution. A directly useful bound is

\[
|\delta I|\le\frac{2\pi|qK|}{m|k|}|a|^\ell
\bigl(|nq|\,|\delta A|+|ku|\,|\delta B|\bigr),
\tag{18}
\]

before perpendicular/harmonic sums and the final density/charge normalization. This is a componentwise absolute bound, not a relative determinant-error estimate. Determinants can amplify tensor errors near a singular wave matrix.

**Previously observed, not recomputed here:** [C2d continuation evidence](NYQUIST_C2D_RESULTS.md#analyticity-and-continuation), lines 291–313, reports n=0,zz direct/grid jump errors about `3.83e-12,3.89e-12,3.99e-12` at probe separation `2e-12` on 32×64,64×128,96×192 grids, respectively. All 171 recorded harmonic/component jump checks passed. The fitted residue differed from the coherent grid residue by `1.4113%,0.40629%,0.17092%`; measured causal-boundary mismatches were `0.0253369,0.00738149,0.00311376` in that integral normalization, near the isotropic control Re(omega)=0.4046. These are archived measurements of that study's setup, not a new validation of every current Chebyshev option or a universal error constant. They illustrate (17), not an analyticity certification from convergence of three numbers.

With an entire fitted residue, `C_PWL_direct + 2L_fit` is analytic **within an open lower half-plane region** where the fit and direct integral are analytic. Its failure is generally in gluing to the chosen upper response across the real axis. Smooth lower-plane values alone do not make it the continuation of that response. A finite mismatch is not removed by choosing the real-axis value to be the average, or by switching formulas exactly at Im(omega)=0.

## 9. Strategies A, B and C

| Strategy | Consistency and convergence | Conditioning/complexity | Disposition |
|---|---|---|---|
| **A: exact PWL sheets** | Exactly continues the represented C2d model locally on a specified sheet. Keeps all upper results and coherent identities. Interior slope jumps create nontrivial monodromy; fixed-depth physical damped convergence is not guaranteed, as the quadratic example proves. | Reuses cell primitives/weights, but needs explicit sheet/domain metadata, endpoint/node handling, deterministic solver queries and cache semantics. Small damping can still be ill-conditioned. | Retain as a mathematical/local diagnostic, or a narrowly bounded near-axis approximation with explicit error estimates. Not the default general damped solution. |
| **B: one globally analytic momentum model** | Use analytic A/B consistently in Cauchy integrals, polynomial moments, coupled zz and residue. Removes interior PWL node branches. Can define a single causal branch across an interior real-frequency domain; finite support and model singularities still constrain it. Physical convergence requires controlling fit, derivatives and complex extrapolation, not merely quadrature. | A direct Chebyshev polynomial permits exact polynomial Cauchy recurrences. LOG10 gives exp(polynomial): no general finite polynomial primitive; controlled quadrature or another justified integration representation is required. Complex extrapolation/degree can magnify noise; no blanket full-real-line integrability. | **Recommended formulation for a useful damped backend**, after a small scalar proof/accuracy stage and renewed upper-half-plane comparison. |
| **C: coherent PWL integral plus analytic fitted residue** | Can approximate the physical result if both representations converge with uniform error in the required complex domain and boundary mismatch is bounded. At fixed grid it generally is not analytic continuation of the validated upper model. A good fitted residue could improve physical damping while leaving the gluing defect. | Smallest code change, but two error models, endpoint/weight discrepancies and an axis seam remain. Existing production behavior supplies compatibility evidence, not a proof. | Useful as an explicitly labeled comparison; not recommended as the definitive coherent causal response. |

For B, suitable consistent derivative operators can preserve existing nodal-difference conventions: apply the same perpendicular FD and, if chosen, parallel FD to an analytic row model throughout both the real integral and complex residue. Finite differences of analytic functions are analytic. Alternatively use analytic parallel derivatives throughout, but this changes the represented data and requires separate validation; it is not a harmless residue-only option. Independent fits to signed A/B are another possible model, but positivity/log fitting no longer applies and compatibility/convergence to derivatives of a single f must be assessed. No particular choice is implemented or silently selected here.

For a fixed analytic pair A*,B*, recommended tensor equations are (2)–(3) with every M and C evaluated from that pair, and

\[
C_\ell^{*,\rm causal}=C_\ell^{*,\rm direct}+2\pi i\sigma a^\ell F^*(a)
\]

on a declared branch domain reached through the interior support interval, with the appropriate PV half-jump on the real boundary. Off-support harmonic paths need their declared sheet rather than a moving resonance switch. This retains the continuum tensor identity, signed-k convention, species factors and finite harmonic MPI organization.

B cannot be exactly identical to the present nontrivial PWL response on an upper-half-plane open set **and** remove its intrinsic node monodromy: uniqueness of analytic continuation forbids that. Its upper response must instead be independently compared with C2d and the analytic physical reference, with finite-grid differences reported. Preserve C2d as an immutable comparator; do not recast such differences as an unnoticed production correction.

A small real-axis fit residual also does not bound the complex residue. High-order Chebyshev terms grow off the fitted interval; finite noisy samples do not uniquely determine a physically accurate analytic distribution. Model degree, regularity assumptions, usable complex strip and support/tail errors need explicit acceptance evidence. The native evaluator's finite output/status is necessary for numerical reliability, not sufficient for physical accuracy.

## 10. Proposed first implementation task — specification only

**Suggested subsequent task: “C3b: standalone causal scalar Cauchy reference with a global analytic polynomial.”** This is a prerequisite mathematical/implementation checkpoint for B, not a damped ALPS backend.

1. Add only an isolated Fortran validation utility. Reuse the C2c direct kernel conventions and C2d moment identities; leave all production and existing experimental sources unchanged. No MPI, root solver, fit replacement or physical simulation in this first checkpoint.
2. Accept fixed support and **one global analytic polynomial**, initially constants through quartics and low-degree complex Chebyshev coefficients. Do not select a polynomial from `real(a)`. Use the same polynomial in the direct integral, moments and residue. A supplied native `CHEB_DIRECT` model can be a later input adapter without refitting it; defer LOG10 integration until a quadrature design is reviewed.
3. Define the direct sheet, a particular interior crossing interval, and an explicit continuation chart/path. Implement C0,C1,C2 and the signed-k causal/PV corrections. For a monomial basis, exact recurrence is `C_j(1)=M_(j-1)(1)+a*C_(j-1)(1)` with `C_0(1)=Lambda`; this supplies an independent analytic reference for the global polynomial, without a second plasma susceptibility.
4. Include the current PWL interpolation of the same quadratic as a **reference diagnostic**, showing its correct local jump, node monodromy and the fixed-depth eta² residue defect. Do not try to smooth or remove those facts in the PWL kernel.
5. Verify both k signs, both a half-planes, real PV/causal limits, interval interiors, node limits, support endpoints and paths of the same/different winding. Query-order independence within a fixed chart is mandatory. Same-chart paths must agree; different monodromy classes must differ by the analytically predicted amount.
6. Use debug bounds/FPE checks and independent higher-precision scalar primitives where needed. Report absolute errors and conditioning, retaining existing mathematical tolerances where the same fixtures apply. For new scalar identities justify thresholds from roundoff/conditioning; do not change physical tolerances.
7. Stop with a compact scalar report and a decision on the analytic row/integration representation. Do not connect to `full_integrate`, add a selector, compute damped roots, or claim C2d upper-half-plane equivalence.

This checkpoint tests the signs/branch/data-model contract needed for B. It is intentionally not presented as completing LOG10 integration, choosing a physically accurate analytic fit, or solving the full tensor continuation problem.

## 11. Minimum subsequent physical validation — designed, not executed

After that scalar gate and approval of a consistent analytic representation/integration method, a separate isolated full-tensor prototype would require the following. These are future tests, not results of C3a.

| Gate | Required evidence |
|---|---|
| Analytic bi-Maxwellian reference | Start with the established isotropic proton/electron case and its independently qualified damped branch near `0.40460138158784 - 0.00185430537032i` (historical analytic reference in C2d). Recalculate the analytic reference under the eventual fixed setup rather than copy a numerical production root as truth. Retain a growing cyclotron control. |
| Upper-half-plane preservation | At identical frequencies compare all six species components, coherent moments/static identities, wave matrices and complex D against both immutable C2d/C2f.3a and analytic NHDS. For B, quantify representation changes; require convergence rather than exact finite-grid identity that is mathematically unavailable. |
| Axis matching and signs | Probe decreasing positive/negative imaginary omega at fixed real frequencies. Compare direct jumps to `2*pi*i*r^ell*F(r)` in a, signed causal limits to (9), and the full prefactors to (10). Test k>0 and k<0 explicitly at fixed q and both charge signs. Distinguish direct jump, causal boundary equality and floating-point PV evaluation. |
| Coupled tensor structure | Check every mode, n=0 combined zz's explicit omega factor, complementary ee exactly once, common A/B moments and perpendicular weights, harmonic truncation, species density/charge normalization and vacuum once. Verify any required endpoint-weight change independently before using it. |
| Nodes/support/branches | Use fixed real frequencies mapping onto cell interiors, nodes and neighborhoods of support endpoints for several harmonics. Same-chart loop/path tests must be consistent, and known endpoint monodromy must not be mistaken for a code discontinuity. Exclude singular points explicitly from an admissible solver domain, not by clipping a result. |
| Analyticity/conditioning | Independent real/imaginary perturbations of scalar integrals, all tensor components and D must show CR convergence before roundoff dominates. Include physical roots, near-axis points and support-sensitive frequencies; record matrix singular values or existing scale-aware diagnostics near small D. No root acceptance based solely on a tiny raw determinant. |
| Damped roots | Use the unchanged ordinary ALPS solver and physical qualifications. Establish a complex analytic domain containing its iterates; do not allow query-history-dependent sheets. Compare growing/damped frequency errors and residuals under independent grid, analytic-model degree and support refinement. A converged root of an inaccurate representation is not validation. |
| Causal model versus residue | Compare the new representation's own real boundary residue with the old fitted residue at identical frequencies, retaining A/B and endpoint contributions separately. Demonstrate that the proposed model's integral and residue agree with each other, not only with a fitted expected root. |
| MPI/compatibility | Keep the original four-rank species/harmonic decomposition and collective disp sequence. Compare per-rank harmonic contributions/call counts, shared-model initialization, failure propagation and deterministic cached/uncached results. A sheet must be part of the declared evaluation domain or an explicit cache key if more than one is supported. No wavevector/frequency parallelism. |
| Subsequent mixture control | Only after the bi-Maxwellian gates, use C2f.4b's **one tabulated proton mixture plus analytic electrons**, with its independently additive analytic reference. This tests model accuracy beyond fitting a Gaussian logarithm exactly; it does not authorize drifting/new VDF families. |

No lower-half-plane Nyquist counting is proposed. The purpose of any later root experiments would be causal accuracy, not developing a new contour capability.

## 12. Questions requiring resolution before a physical damped-root claim

- Which analytic object will define the candidate: f with declared derivative operators, or consistent analytic approximants of A/B? Which physical regularity assumptions justify it for tabulated data?
- How will a native LOG10 exp(polynomial) be integrated accurately with the same model in moments and residue, without reintroducing a frequency-dependent quadrature/stencil switch? Is the initial restricted domain small enough for controlled analytic extrapolation?
- What fixed support, tail-error bound and per-harmonic continuation chart are admissible around the target damped roots? An analytic row model is not automatically an integrable whole-real-line VDF.
- How will the actual perpendicular endpoint mismatch between Q and Q_L be handled and independently qualified? Small tails cannot be assumed zero, and existing method-zero/one zero-return behavior is not an analytic theorem.
- Can the selected model/degree represent the derivative and complex residue accurately enough for electrons as well as protons? Real-grid fit diagnostics alone cannot establish this.
- How will a solver request outside its approved chart or a failed/underflowed fit evaluation terminate collectively, without inventing a zero susceptibility or a history-dependent sheet?
- What conditioning limit permits a meaningful static/near-axis comparison after raw/conventional normalization is distinguished? The prior nearly parallel `~1e-8` conditioning observation is not a universal damped-domain limit.

**Disposition:** the fixed-cell continuation is locally well defined and reviewable; it is not mathematically impossible. Its node structure, global sheet ambiguity and fixed-depth convergence limitation rule out treating a cell-selected residue as a generally satisfactory causal backend. A consistent analytic representation is the strongest route to a controlled damped susceptibility, but it must earn new upper/lower accuracy and tensor-identity validation. The hybrid remains an explicitly approximate comparator. No source, solver, integration, fit, MPI or Nyquist change is made by this audit, and no C3b work is started.
