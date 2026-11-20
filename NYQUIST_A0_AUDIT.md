# ALPS Nyquist stability analysis — Stage A0 implementation audit

Audit date: 2026-11-20 (America/Phoenix). Source revision: `5ee0d335a037744a5ce298c355b84fc316af5f9c`, with the pre-existing working-tree changes left intact. Line references below describe the source inspected during this audit.

**Scope:** implementation audit and a small ordinary bi-Maxwellian root calculation only. No production source, numerical routine, build configuration, or existing input was modified. No Nyquist routine or winding calculation was implemented. The two supplied papers were read as scientific references, not as instructions.

## 1. Executive assessment

**Yes: Nyquist analysis can reuse the public collective function `alps_fns:disp(om)` after the existing initialization, without substantial restructuring or invoking a root solver.** Keep the present decomposition over species/cyclotron harmonics. Add a separate Fortran Nyquist module, a narrow input/dispatch branch, and a rank-0-controlled collective evaluation loop.

The mathematical qualification is essential:

- The returned scalar is the determinant of a **frequency-scaled wave tensor**, not the determinant of the usual dimensionless dielectric wave tensor. In proton reference units, the relation is \(D_{\rm ALPS}=(\omega v_A)^6D_{\rm conventional}\), where both \(\omega\) and the input \(\mathrm{vA}=v_{Ap}/c\) are dimensionless.
- This factor has zero winding on a rectangle whose **entire interior and boundary lie at strictly positive imaginary frequency**. For finite-temperature bi-Maxwellians, fixed real \(k_\parallel\ne0\), and valid physical species parameters, the underlying scaled tensor is analytic there and has no physical upper-half-plane susceptibility poles. Thus the existing determinant is the appropriate initial observable; **no additional frequency normalization is needed**.
- Do not extend that conclusion to contours enclosing/touching \(\omega=0\), the real axis, cold-species resonances, or unsupported \(k_\parallel=0\) evaluations.
- The numerical-VDF path is callable through the same interface, but its moving resonance windows and threshold approximations are not an exactly holomorphic frequency representation. Physical analyticity does not certify the implemented approximation. Numerical-VDF counts require a later convergence/continuity assessment, not merely the same phase sampler.
- A finite rectangle counts zeros **inside that rectangle, with multiplicity**. It does not establish the total number of unstable modes at that wavevector without evidence excluding roots outside its real/imaginary limits and below \(\gamma_{\min}\). A zero count with \(\gamma_{\min}>0\) is a thresholded, bounded result.

A small existing-solver run located useful growing reference roots for the requested equilibrium at \((k_\perp d_p,k_\parallel d_p)=(0.001,0.3)\):
\[
\omega/\Omega_p\simeq \pm0.39666+0.085986i.
\]
The isotropic control has corresponding damped reference roots near \(\pm0.40460-0.0018543i\). These are reference candidates, not an exhaustive root inventory or an assumed winding number (§6).

## 2. Code map and repeated-evaluation contract

### 2.1 Initialization and persistent data

| Responsibility | Exact location and routines | Data/dependency |
|---|---|---|
| MPI startup and process-count restriction | [src/ALPS.f90:58](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS.f90:58), lines 58–87 | `mpi_init`, size/rank, `proc0`; currently requires even `nproc>2`; rank 0 calls `init_param`, all call `pass_instructions`. |
| Input, species and analytic species parameters | [src/ALPS_io.f90:82](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_io.f90:82), `init_param`, especially 120–263; `spec_read` 647–688; `bM_read` 734–768 | `kpar,kperp,vA`; `ns,qs,ms`; `usebM,relativistic`; `bMbetas,bMalphas,bMpdrifts,bMnmaxs,bMBessel_zeros`. |
| Broadcast inputs and allocate determinant/VDF arrays | [src/ALPS_com.f90:193](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_com.f90:193), `pass_instructions` 193–363 | Broadcasts scalars/species arrays/guesses, sets `pi` (254), allocates `nmax`, root-only `wave,chi0,chi0_low` (344–349), and all-rank `pp,df0` (359–361). |
| Read equilibrium and differentiate | `ALPS.f90` 89–112; [src/ALPS_io.f90:989](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_io.f90:989), `read_f0` 989–1054; [src/ALPS_fns.f90:34](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_fns.f90:34), `derivative_f0` 34–254 | Rank 0 reads numerical distributions and computes derivatives/diagnostics. Bi-Maxwellian species use zero placeholder `f0,pp,df0` and analytic density/current diagnostics; no VDF files are read for them. Original `f0` is deallocated before determinant evaluation. |
| Continuation models, once per equilibrium | [src/ALPS_analyt.f90:142](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_analyt.f90:142), `initialize_rel_cheb` 142–177, `initialize_nr_cheb_production` 691–706, `determine_param_fit` 2013–2357 | Called on rank 0 before releasing `f0`. Bi-Maxwellian species skip fitting (2093–2099). Numerical species retain fit/model data. |
| Distribute equilibrium and continuation models | `ALPS.f90` 116–123; `ALPS_com.f90:pass_distribution` 367–448; `complete_nr_cheb_initialization` 32–51; `pass_nr_cheb_models` 53–107; `pass_rel_cheb_models` 117–180 | All-rank initialization; broadcasts `pp,df0,param_fit,fit_type,perp_correction`, relativistic grids/derivatives when applicable, and selected Chebyshev models. Allocation/status agreement includes `MPI_ALLREDUCE`. |
| Harmonic limits and work assignment | [src/ALPS_fns.f90:4054](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_fns.f90:4054), `determine_nmax` 4054–4158; `split_processes` 4162–4290 | `nmax(nspec)`, per-rank `sproc,nlim(2)`. Numerical limits depend on \(k_\perp\), momentum support, and `Bessel_zero`; allocation can add numerical harmonics when ranks exceed nominal work. |
| Cached ordinary Bessel values | `ALPS.f90` 125–132; `ALPS_fns.f90:determine_bessel_array` 4297–4338 | Worker-local `bessel_array(nlim(1)-1:nlim(2)+1,0:nperp)`; deallocated/reallocated when rebuilt; independent of frequency. |
| Global data definitions | [src/ALPS_var.f90:100](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_var.f90:100), 100–134, 218–310, 361–384 | MPI state; normalized wavevector and `vA`; equilibrium grids/derivatives; harmonic assignment; root-only tensor/physical susceptibility snapshots; analytic-species settings. |

The original `f0` is not needed for ordinary subsequent `disp` calls. This does not make `disp` self-initializing: the arrays, species settings, model transfer and harmonic assignment above must already exist.

### 2.2 Susceptibility, tensor, determinant

| Responsibility | Exact routine/location | Implemented behavior |
|---|---|---|
| Public determinant entry | [src/ALPS_fns.f90:258](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_fns.f90:258), `disp(om)` 258–621; public declaration 26 | Complex input and complex return, both double precision. No root-search call inside. |
| Analytic bi-Maxwellian/cold response | `disp` 335–364; [src/ALPS_NHDS.f90:59](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_NHDS.f90:59), `calc_chi` 59–233 | Only the worker assigned a species interval starting at \(n=0\) calls `calc_chi`; it computes that species' whole analytic harmonic sum. Other workers assigned that analytic species contribute zero. |
| Analytic harmonic tensor | `ALPS_NHDS.f90:calc_ypsilon` 241–368 | \(\Omega_j=q_j/m_j\), thermal/drift speeds, resonance argument, \(A_n,B_n\), scaled tensor entries. Handles exactly \(k_\perp=0\) explicitly at 340–366. |
| Plasma dispersion function | `ALPS_NHDS.f90:dispfunct` 479–520; `WOFZ` 528–731 | Faddeeva evaluation, with a fixed-\(k_\parallel\) sign convention. `besselI` 457–473 and `BESSI/BESSI0/BESSI1` 738–846 supply exponentially scaled real Bessel functions. |
| Cold response | `ALPS_NHDS.f90:calc_chi_cold` 372–451 | Selected by exactly `bMbetas(j)==0` at 135–140. Explicit real-axis resonance denominators. |
| Numerical-VDF harmonic response | `disp` 369–493 | Assigned nonnegative \(n\), explicitly adds both signs for \(n>0\); includes longitudinal `int_ee`/`int_ee_rel` once on the \(n=0\) worker. |
| Resonance detection/integration selection | `ALPS_fns.f90:determine_resonances` 626–730; `full_integrate` 735–777 | Detects resonances using real parts of complex resonance momenta; chooses direct or resonance integration. Additional Landau terms occur only at \(\operatorname{Im}\omega\le0\). |
| Nonrelativistic numerical integration | `ALPS_fns.f90:integrate` 784–849; `integrate_res` 855–1255; `funct_g` 1260–1366; `resU` 1653–1689; `int_T` 1694–1793 | Fixed-grid trapezoids or moving resonance windows, local interpolation and a small-imaginary-part approximation. |
| Numerical continuation and longitudinal term | `ALPS_fns.f90:landau_integrate` 1372–1545; `int_ee` 1550–1648; `int_T_res` 1797–1895; `ALPS_analyt.f90:eval_fit` 1705–1926 | Fitted/analytic/Chebyshev continuation is used for residue terms, not added on the strict upper-half-plane path. `int_ee` accepts `om` but does not use it. |
| Relativistic numerical integration | [src/ALPS_fns_rel.f90:688](/home/kgklein/Codes/ALPS_git/ALPS/src/ALPS_fns_rel.f90:688), `integrate_res_rel` 688–739; `integrate_resU_rel` 744–947; `principal_integral_rel` 952–1141; `funct_g_rel` 1146–1257 | Transformed \((\Gamma,\bar p_\parallel)\) integration with resonance windows and threshold treatment. `landau_integrate_rel` 1263–1359, `int_ee_rel` 1364–1482, `resU_rel` 1487–1516, `int_T_rel` 1523–1635, `int_T_res_rel` 1640–1760 complete this branch. |
| Species normalization and MPI assembly | `disp` 345–360, 497–520 | Analytic result divided by `ns*qs`, then multiplied back by the same factor; numerical result multiplied by `ns*qs`. Two complex sum reductions to rank 0. |
| Physical susceptibility snapshots | `disp` 533–542 | Stores `chi0=chi/(om*om*vA*vA)` and `chi0_low`, with tensor symmetries, on rank 0 for eigenfunction/heating use. These divided snapshots do not feed the determinant assembly. |
| Full tensor construction | `disp` 554–600 | Sums species into local `eps`, fills non-Hermitian tensor symmetries, adds \((om*vA)^2 I\), then vacuum wavevector terms into global `wave`. |
| Determinant and return | `disp` 603–617 | Explicit symmetric \(3\times3\) determinant expansion; broadcasts completed scalar; world barrier. |

### 2.3 Frequency and wavevector state

For a fixed equilibrium and wavevector, every `disp(om)` invocation recomputes the frequency-dependent response. Its local work arrays are reset at 318–322. Analytic species recompute \(A_n,B_n,Z\) and the harmonic sum; the Bessel stopping rule depends on \(k_\perp\) and species data, **not on \(\omega\)**. Numerical species redetermine resonance positions/windows for each frequency. There is no last-frequency cache, stored root iterate, or dependency on the order of previous frequencies in this entry point.

The call has side effects: it overwrites rank-0 `wave,chi0,chi0_low` and shared `ierror`, and performs collectives. It is not pure, thread-safe, or suitable for simultaneous independent frequency calls. Copy any needed diagnostic tensor before evaluating another frequency. No ordinary per-frequency file output occurs inside `disp`; called error paths can print and abort.

One startup portability issue is visible: `sproc` is declared without initialization (`ALPS_var.f90` 188), and `split_processes` assigns it only to workers (4273–4282), while the main program tests `sproc==0` to skip rank-0 Bessel setup (130). `pass_instructions` initializes `nlim` (257), but not `sproc`. The present executable successfully treats rank 0 as the zero sentinel, but this is not an explicit Fortran initialization contract. A future entry path should use the already-defined `proc0` for its worker guard or explicitly establish that sentinel; no A0 correction was made.

Persistent fits and Bessel arrays are equilibrium/wavevector state, not determinant caches. The relativistic disk cache is used during `derivative_f0_rel` (`ALPS_fns_rel.f90` 34–97); its `cache_identity` (`ALPS_rel_cache.f90` 336–375) uses source grids/VDF, mass, `vA`, grid resolution and build identity, not frequency. Chebyshev model arrays/activation flags/counters persist in `ALPS_analyt.f90` 70–106, but the strict upper-half-plane path does not call the Landau continuation evaluator.

Wavevectors are read once and broadcast in `pass_instructions` (224–225). Existing scans update them on every rank; when \(k_\perp\) changes they repeat `determine_nmax`, `split_processes`, and worker `determine_bessel_array` (`om_scan`, 2548–2555). A change only in \(k_\parallel\) does not require rebuilding ordinary Bessel arrays; its frequency/resonance dependence is read afresh. Retain this dependency if scans are added later.

**Repeated-evaluation contract:** initialize MPI and all shared state once; hold the equilibrium and wavevector fixed; have every rank call `disp` with the same prescribed complex frequency in the same collective order. No call to `secant`, `secant_osc`, `rtsec`, `refine_guess`, or `map_search` is required.

## 3. Mathematical assessment

### 3.1 The implemented tensor and its scaling

Let \(x=\omega_{\rm physical}/\Omega_p\), \(a=v_{Ap}/c=\texttt{vA}>0\), and \(\boldsymbol{\kappa}=(k_\perp d_p,0,k_\parallel d_p)\). Use proton reference units initially, so \(\Omega_p>0\) and positive \(\operatorname{Im}x\) means growth. In this section \(x\) is the code argument `om`.

Define \(\widetilde\chi_j\) as the completed species tensor named local `chi` after the two reductions, and \(S=\sum_j\widetilde\chi_j\), with the symmetry completion used in the source. The actual stored `wave` is
\[
 A(x)=a^2x^2I+S(x)+\boldsymbol{\kappa}\boldsymbol{\kappa}^{T}
                 -\kappa^2I
 =
 \begin{pmatrix}
 a^2x^2+S_{xx}-\kappa_\parallel^2&S_{xy}&S_{xz}+\kappa_\perp\kappa_\parallel\\
 -S_{xy}&a^2x^2+S_{yy}-\kappa_\parallel^2-\kappa_\perp^2&S_{yz}\\
 S_{xz}+\kappa_\perp\kappa_\parallel&-S_{yz}&a^2x^2+S_{zz}-\kappa_\perp^2
 \end{pmatrix}.
\]
This follows directly from `disp` 557–600, including the signs; no conjugation is taken. For nonzero \(x\), the code's physical susceptibility snapshot is
\[
 \chi^{\rm phys}_j=\frac{\widetilde\chi_j}{a^2x^2}.
\]
The conventional wave tensor with this susceptibility is
\[
 M=I+\sum_j\chi^{\rm phys}_j+
    \frac{\boldsymbol{\kappa}\boldsymbol{\kappa}^{T}-\kappa^2I}{a^2x^2},
 \qquad A=a^2x^2 M.
\]
Consequently,
\[
 D_{\rm ALPS}(x)=\det A=a^6x^6\det M,\qquad
 W[D_{\rm ALPS};C]-W[\det M;C]=6\,\operatorname{Ind}_C(0).
\]

For the proposed rectangle with \(\gamma_{\min}>0\), the origin is outside, \(\operatorname{Ind}_C(0)=0\), and the windings agree. Multiplication by the nonzero constant \(a^6\) has no winding effect. This is a contour-specific justification, not a claim that algebraically equivalent root equations always have equal winding.

For a meromorphic determinant with no zero or pole on a counterclockwise contour, the argument principle is \(W=(2\pi)^{-1}\Delta_C\arg D=N_{\rm zeros}-N_{\rm poles}\). The warm upper-half-plane assessment below supplies the no-poles condition needed to interpret this winding as a zero count.

The formula does **not** establish six extra zeros at the origin: \(\det M\) may itself have a pole or zero there, and orders can cancel. A contour including zero requires the actual local order and a separate analysis. Do not divide out guessed “trivial” roots.

There is a misleading normalization comment: `ALPS_NHDS.f90` 57 says elements were multiplied by \(k_\perp^2d_p^2\), and `disp` 326 mentions an additional \(k_\perp^2\). **The active expressions do not apply a common \(k_\perp^2\) factor to the wave tensor/determinant.** In `calc_ypsilon`, the nonzero-\(k_\perp\) terms still contain \(1/z\), where \(z\propto k_\perp^2\) (318–355); the exact parallel limit has finite nonzero transverse entries (341–348). `zp` is assigned but not used. The vacuum entries in `disp` are simply \(a^2x^2\) and quadratic wavevector terms. The audit therefore uses the executable expressions above, not those comments. Even a genuine nonzero fixed-\(k\) scalar would have zero frequency winding; a zero factor at \(k_\perp=0\) would be a different issue.

Other conditioning/scaling operations are frequency independent: \(\ell_j^{-2}=n_jq_j^2/m_j\), cancellation of the temporary `ns*qs` division/multiplication for analytic species, and \(I_n(z)e^{-z}\) evaluated directly. The Bessel exponential is the physical harmonic coefficient, not a new frequency zero/pole. Require positive finite density/mass, nonzero charge, nonzero finite `vA`, positive temperature/anisotropy for the warm branch; the analytic wrapper still explicitly divides by `ns*qs` and cannot safely represent a zero-density species by setting `ns=0`.

### 3.2 Warm bi-Maxwellians: exact formula versus numerical evaluation

For each species, `calc_chi` and `calc_ypsilon` use
\[
 \Omega_j=\frac{q_j}{m_j},\quad
 \ell_j^2=\frac{m_j}{n_jq_j^2},\quad
 v_{\parallel j}=\sqrt{\frac{\beta_{\parallel j}}{n_jm_j}},\quad
 U_j=\frac{p_{{\rm drift},j}}{m_j},\quad
 z_j=\frac{k_\perp^2v_{\parallel j}^2\alpha_j}{2\Omega_j^2},
\]
and
\[
 \zeta_{nj}=\frac{x-k_\parallel U_j-n\Omega_j}{k_\parallel v_{\parallel j}}.
\]
Here masses/densities/charges and velocities are in the code's reference units. For fixed real \(k_\parallel\ne0\), all denominators outside the plasma dispersion function are constants in \(x\); \(A_n,B_n\) are polynomials in \(x\) times that function plus polynomials (`calc_ypsilon` 316–326).

The exact function underlying `WOFZ` is entire:
\(w(z)=e^{-z^2}\operatorname{erfc}(-iz)\). This definition and its entirety are given by [NIST DLMF §7.2(i)](https://dlmf.nist.gov/7.2#i). The implemented sign convention is
\[
 Z_{\rm code}(\zeta)=
 \begin{cases}
 i\sqrt{\pi}\,w(\zeta),&k_\parallel>0,\\
 -i\sqrt{\pi}\,w(-\zeta),&k_\parallel<0 .
 \end{cases}
\]
This is `dispfunct` 513–519. The choice is fixed by the wavevector, not by the frequency. For a strict upper-half-plane frequency, the argument passed to `WOFZ` is in its upper half-plane for **either** sign of \(k_\parallel\). No frequency branch cut or conjugation is introduced by this convention.

Thus the finite harmonic sum, evaluated with exact \(w\), gives an entire scaled warm-species response at fixed supported wavevector. Its continuation across a cyclotron resonant real frequency is finite for \(k_\parallel\ne0\), unlike the cold model. The physical velocity-integral explanation is consistent: denominators \(x-k_\parallel v_\parallel-n\Omega_j\) cannot vanish for real particle velocities and \(\operatorname{Im}x>0\). With converged thermal/harmonic sums, there are no physical susceptibility poles in this open half-plane. Instability appears as a zero of the wave determinant, not a pole of a species response.

Numerical qualifications:

- `WOFZ` uses different power-series, Taylor/continued-fraction regions, variable truncation orders, and quadrant transformations (591–724). Its finite-precision map is not literally an entire mathematical function. These are approximations to one entire function, rather than physical frequency branch cuts, but joins/roundoff must be small compared with the determinant's distance from zero along the contour.
- The embedded accuracy statement is about 14 significant digits in the upper argument half-plane (543–547). It is **not** an error bound on a cancellation-dominated susceptibility or determinant.
- `WOFZ` can signal overflow with `FLAG=.TRUE.` (599, 712, 728); `dispfunct` does not inspect this returned flag. An early exit need not provide a valid new \(U,V\). This is a confirmed missing error propagation, not evidence that the initial contour will encounter it.
- Real Bessel approximations use modest-length coefficient tables (800–844), so full determinant accuracy need not reach the precision claimed for `WOFZ`. These coefficients/truncation limits are constant with frequency and do not themselves introduce frequency branch switches.
- Harmonic truncation in `calc_chi` (162–186) is frequency independent, but uniform accuracy over an enlarged contour has not been established. A coefficient cutoff alone does not bound every response element after multiplication by frequency-dependent factors.
- Large \(|\zeta|\) causes cancellation in combinations such as constant terms plus \(\zeta Z\). An analytic exact formula can still be evaluated inaccurately there.

**Assessment:** the existing returned determinant is mathematically suitable for the argument principle on a strict upper-half-plane rectangle for these warm bi-Maxwellian cases. Its numerical winding is trustworthy only after sampling and response accuracy checks. A new frequency-dependent normalization is not a remedy for susceptibility error.

### 3.3 Zero frequency, zero wavenumbers, and cold species

At \(x=0\), `disp` unconditionally forms `chi0` and `chi0_low` by division by \(a^2x^2\) (533–534). Even when the scaled tensor has a finite analytic limit, the current function is not a clean zero-frequency evaluation interface: those diagnostic assignments can yield NaN/Inf or traps. The determinant expansion itself does not divide by \(x\). Initial contours should exclude the origin and should not attempt to evaluate it.

Crossing **\(\operatorname{Re}x=0\) at \(\operatorname{Im}x>0\)** is valid: \(x\ne0\), so it is not the zero-frequency singularity. Negative real frequencies are equally legitimate. No special half-contour replacement or root-count doubling is justified by the determinant interface.

Warm `calc_ypsilon` has divisions by \(k_\parallel\) and \(k_\parallel^2\), including the drift contribution to \(\chi_{zz}\) in `calc_chi` 196/207. The numerical resonance paths also divide by \(k_\parallel\). **Exactly perpendicular propagation is not supported by this warm entry path as written.** A nearly parallel initial wavevector has nonzero \(k_\parallel\), so it avoids this issue. Exactly parallel \(k_\perp=0\) has explicit limits in the warm analytic path and in numerical `int_T`; small nonzero \(k_\perp\) is nevertheless a useful initial oblique cross-check.

The cold branch has denominators \(x-k_\parallel U_j\) and \(x-k_\parallel U_j\pm\Omega_j\), including a squared denominator (`calc_chi_cold` 430–438). For real inputs its potential poles are on the real axis; actual determinant pole orders can be reduced by cancellations. The strict positive-\(\gamma_{\min}\) rectangle excludes them. At zero frequency and zero drift, the expression \(x^2/(x-k_\parallel U_j)^2\) is literally \(0/0\) even if a limiting expression exists. Do not infer an everywhere-entire determinant for cold cases or use a real-axis contour without a separate pole/limit treatment. Initial validation uses positive-beta warm species only.

### 3.4 Numerical VDFs: physical analyticity is not implementation analyticity

The exact causal susceptibility for a real gyrotropic distribution with sufficient decay/regularity has no resonant velocity denominators vanishing in the strict upper frequency half-plane. At a fixed grid, the direct `integrate` path is a finite weighted sum of `resU` rational functions; its poles are real:
\[
 x=\frac{k_\parallel p_\parallel+nq_j}{\Gamma m_j}
\]
(`resU` 1680–1685, real positive \(\Gamma,m_j\)). Hence a fixed-grid direct sum is holomorphic in the open upper half-plane, although it approximates the continuum response.

**The actual dispatcher does more than this fixed sum.** Confirmed sources of nonholomorphic numerical dependence are:

1. `determine_resonances` chooses paths using \(\operatorname{Re}p_{\rm res}\) and grid-domain membership (693–727).
2. `integrate_res` chooses cells with `floor(real(p_res)...)`, moves integration limits at half-cell boundaries, and has special edge returns (957–1020).
3. Its quadrature nodes and interpolated numerator depend on \(\operatorname{Re}p_{\rm res}\); `funct_g` changes interpolation stencils with real momentum (1319–1362).
4. It switches at \(|\operatorname{Im}p_{\rm res}|=\texttt{Tlim}\) (1031–1101). The small-imaginary-part approximation explicitly separates real and imaginary parts and uses a limiting residue term (1112–1164). For nonrelativistic species the switch occurs at \(\gamma=\texttt{Tlim}|k_\parallel|/m_j\), potentially at a positive height inside a proposed rectangle, separately for each species.
5. The relativistic branch similarly uses real resonance/cone tests (`ALPS_fns_rel.f90` 823–849), moving bounds (875–899), and `Tlim` switching in `principal_integral_rel` (1029–1101). The corresponding imaginary resonance scale is \(\Gamma\gamma\,\texttt{vA}/|k_\parallel|\).

These are finite-resolution approximations to analytic integrals. Their explicit dependence on real/imaginary parts is not merely the rounding inherent in any complex arithmetic; no exact holomorphic representation or exact matching of all switches is established. The size of the resulting mismatch and its effect on a winding are **unmeasured**, not a demonstrated false physical count. Moving only the contour above `Tlim` would not remove all moving-grid concerns.

In the strict upper half-plane, `full_integrate` adds **no** `landau_integrate` or `landau_integrate_rel` residue (756–774). Thus possible complex-power/square-root conventions in fitted kappa/Jüttner continuation (`ALPS_analyt.f90` 1976–1984), complex square roots in `int_T_res_rel` (1682), or continuation overflow/zero-return policies are not traversed by the initial contour. They become relevant if frequencies on/below the real axis are later allowed. Existing continuation setup must still be shared correctly.

A continuous nonanalytic approximate boundary curve has a winding as a curve, but that integer alone is not an argument-principle count of its own interior frequency zeros. To transfer its winding to the physical determinant, establish adequate response convergence and a boundary error margin: for example, if an analytic reference \(D\) and a continuous approximation \(\widehat D\) obey \(|\widehat D-D|<|D|\) everywhere on \(C\), their boundary curves are homotopic without crossing zero and have equal winding. A practical convergence study can support confidence but is not that rigorous error bound.

**Disposition:** retain this backend and interface; defer physical certification of numerical-VDF windings. No multiplicative determinant normalization removes moving-window/interpolation errors. No susceptibility refactor or repair is authorized in A0.

## 4. MPI assessment

### 4.1 Existing decomposition and collective sequence

Rank 0 assembles the wave tensor. Worker ranks receive a species `sproc` and nonnegative harmonic interval `nlim` from `split_processes` (4162–4290). Numerical workers evaluate \(n=0\) once and paired \(\pm n\) for positive \(n\); the longitudinal `ee` term belongs to the worker whose lower bound is zero. Analytic species are a special case: only their worker with a zero lower bound performs the entire internal NHDS harmonic sum. Preserve both behaviors. “All ranks participate” does not mean that every rank performs an equal number of susceptibility operations.

Every rank calls `disp`, including rank 0 and workers with no nonzero contribution. Each invocation executes, in order:

1. `MPI_REDUCE(schi,chi,...,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD)` — 516–517.
2. The analogous reduction of `schi_low` — 519–520.
3. Rank-0 species summation, tensor and determinant assembly — 523–609.
4. `MPI_BCAST(disp,...,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD)` — 613.
5. `MPI_BARRIER(MPI_COMM_WORLD)` — 617.

All ranks receive the scalar determinant. Only rank 0 owns the completed `wave,chi0,chi0_low`; the reduction receive arrays on other ranks must not be used as completed responses. There is no frequency broadcast inside `disp`.

Current callers generate matching frequencies on all ranks: input guesses are broadcast in `pass_instructions` (317–318), `map_search` constructs a common grid on all ranks (3779–3802), and root solvers update local identical iterates from the broadcast determinant. `map_search` broadcasts selected minima/guesses before `refine_guess` (3861–3864). These root-search control paths should be bypassed for Nyquist.

### 4.2 Minimal reliable rank-0 control

Have all ranks enter the future Nyquist driver. On every iteration, rank 0 selects a command such as evaluate/finish/unresolved and broadcasts it. For evaluate, broadcast **one complex frequency** from rank 0; then all ranks call the existing `disp` once. After return, rank 0 records the result, updates the contour sample cache/refinement queue, and chooses the next command. Workers return to the same command broadcast. Finish or unresolved ends this loop collectively and reaches the shared program finalization.

The extra command/frequency broadcasts belong in the new driver or its adapter, leaving the legacy `disp` communication sequence untouched. Do not start a rank-0-only `disp`, launch independent rank-local refinement, or decompose a list of frequencies/wavevectors across ranks.

Constraints for adaptive sampling:

- Broadcast every request/termination decision. Reuse of a cached sample on rank 0 must cause no unmatched worker determinant call; workers simply remain in the command loop until the next request.
- Identify samples by their contour positions and assemble phase in **geometric contour order**, even if the refinement queue evaluates them out of order.
- Keep the closing segment and all corners in the ordered contour. The existing barrier is sufficient to separate completed calls; extra barriers are unnecessary.
- No rank may skip either susceptibility reduction just because its local work is zero. Reductions still include root and idle contributions.
- Recoverable bad values detected after a completed call must lead to a broadcast unresolved/stop decision, never a rank-0-only early return.
- The initial warm path has no structured susceptibility error status. `alps_error` uses `MPI_ABORT` (`ALPS_io.f90` 1461), and overflow traps may terminate before return. A new driver cannot reliably turn every existing hard failure into a recoverable unresolved result without a later, narrowly justified status interface.

The existing executable requires even `nproc>2` (`ALPS.f90` 77–78); use four ranks for the initial two-species case. Numerical harmonic-limit adjustment can depend on rank count (`determine_nmax` 4131–4147), so later MPI reproducibility comparisons must record actual harmonic coverage. No decomposition change is needed.

## 5. Recommended integration design for A1/A2

### 5.1 Dispatch and shared startup

Use a new `src/ALPS_nyquist.f90` module in a future stage. Keep the phase/contour algorithm independent of susceptibility formulas: it should consume prescribed complex-frequency/value samples through an evaluator interface and maintain contour/refinement state. The ALPS-specific adapter should know only the existing collective `disp` contract, MPI control, and optional evaluation diagnostics. A synthetic scalar evaluator can then exercise the mathematical algorithm without plasma initialization.

The smallest production entry point is **after `ALPS.f90` 132–136, before 142**, once input, equilibrium, continuation models, harmonic splitting, and Bessel arrays are ready. Select Nyquist versus the complete existing operation, preserving the present `use_map` versus `refine_guess` branch inside the normal operation. The Nyquist branch must also bypass the subsequent `om_scan/om_double_scan` block at 152–170; it should reach the common finalization at 173–181.

There is an earlier input dependency: `init_param` currently reads map parameters or `nroots` guesses **before reading species** (155–178). Merely adding the late execution branch would still require unnecessary map/guess input. Future input handling should conditionally skip that block for Nyquist. Assign defined neutral values to root/map bookkeeping used by shared initialization; retaining the small existing `wroots` allocation is less disruptive than redesigning it. Likewise, guard the map-only broadcasts in `pass_instructions` (238–252) so Nyquist never broadcasts uninitialized map scalars if `use_map` happens to be true.

Reuse unchanged:

- `init_param`'s species/equilibrium portions, `pass_instructions`' common allocation/broadcast portions, and all equilibrium preparation in `ALPS.f90` 89–130.
- `read_f0`, `derivative_f0`, analytic-continuation initialization/fitting, `pass_distribution`, and `complete_nr_cheb_initialization`.
- `determine_nmax`, `split_processes`, and `determine_bessel_array`, once per initial wavevector.
- The entire `disp` tensor/determinant path, NHDS susceptibility formulas, numerical integrators and harmonic decomposition.
- All root-search/map/scan routines and their default dispatch behavior. The future driver should not call their search or output machinery.

Cleanup is modest: the current program has **no common equilibrium finalizer**. It explicitly deallocates the original `f0` during startup, then relies on process exit for most other arrays; `mpi_finalize` and error-log closing are the reusable common termination path. The Nyquist module should release its own contour/sample storage before that path. Existing `release_nr_cheb` (`ALPS_analyt.f90` 667–671) and `release_rel_cheb_rows` (507–514) can be reused if explicit model release becomes necessary; do not add broad lifecycle refactoring merely for a single calculation.

### 5.2 Minimal input and interface changes

A minimal future input contract is:

| Item | Recommended behavior |
|---|---|
| Operation selector | Add `use_nyquist=.false.` to `&system`, preserving current behavior when absent; broadcast it before conditional root/map input broadcasts. |
| Wavevector/equilibrium | Continue using current `kperp,kpar`, `&spec_j`, `&bM_spec_j`, and numerical-VDF inputs. No second species/equilibrium specification. |
| Contour | A separate `&nyquist` group containing four explicit bounds: `omega_min,omega_max,gamma_min,gamma_max`, in \(\Omega_p\) units for the initial proton-reference use. Require \(\omega_{\min}<0<\omega_{\max}\) for the initial two-sided configuration and \(0<\gamma_{\min}<\gamma_{\max}\). Do not silently choose physical bounds from root-search tolerances. |
| Sampling controls | Document conservative defaults for initial edge sampling, phase-increment target strictly below \(\pi\), refinement/evaluation caps, and a convergence tolerance. Expose only the controls needed to reproduce/refine the calculation; they are numerical controls, not extra plasma parameters. |
| Single wavevector | Require `n_scan=0` in initial Nyquist mode; reject an incompatible scan request explicitly. Root guesses and `&maps_1` must not be required. |
| Result | Report the wavevector, bounds, normalization, raw accumulated phase/integer residual, sample/evaluation count, and a convergence status. Publish an integer count only when the calculation meets its convergence criteria. An unresolved result is not zero. |

`disp` is already public with the needed scalar argument/result, so an initial warm-species prototype needs **no change to its signature**. Mode/configuration declarations and broadcasts, early input guards, late dispatch, and module build dependencies are the small required interfaces.

A structured collective error status or an optional diagnostic interface may later be necessary for robust recovery, especially for numerical VDFs and ignored special-function failures. Finite checks in the adapter can reject invalid returned values, but cannot detect every plausible finite bad response or catch every pre-return trap. Any such interface extension should preserve the old call form and arithmetic when absent, and needs separate review; do not change legacy formulas or claim that the scalar return already provides error certification.

In a future build change, list the module and explicit Fortran dependencies in [src/Makefile.am:15](/home/kgklein/Codes/ALPS_git/ALPS/src/Makefile.am:15) (source list 15–24; dependency rules 41–49). Keep it downstream of the determinant backend so it does not create a circular dependency. No Python runtime or new external numerical library is needed.

## 6. Contour, precision, and initial bi-Maxwellian validation

### 6.1 Rectangular contour and what it counts

The proposed rectangle is appropriate for the initial warm cases. Traverse
\[
(\omega_{\min},\gamma_{\min})
\longrightarrow(\omega_{\max},\gamma_{\min})
\longrightarrow(\omega_{\max},\gamma_{\max})
\longrightarrow(\omega_{\min},\gamma_{\max})
\longrightarrow(\omega_{\min},\gamma_{\min}),
\]
which is counterclockwise. Evaluate all four edges, including the top and closing edge. No PLUMAGE real-axis asymptotic crossing correction belongs in this closed-contour calculation.

The intended result is
\[
N_C=\#\{x:\ D_{\rm ALPS}(x)=0,\ 
\omega_{\min}<\operatorname{Re}x<\omega_{\max},\
\gamma_{\min}<\operatorname{Im}x<\gamma_{\max}\},
\]
with algebraic multiplicities, subject to the no-poles/accuracy qualifications already given. A root on the boundary makes the usual winding count undefined; do not assign half a root or round through an undefined phase.

Frequency/wavevector considerations:

- Use the **same reference units as the existing solver**, \(x=\omega/\Omega_p\), \(\kappa=k d_p\), not the papers' \(k\rho_p\) without conversion. With the supplied proton parameters, the code gives \(v_{\parallel p}/v_A=\sqrt3\) and \(v_{\perp p}/v_A=\sqrt6\), so \(k\rho_{\parallel p}=\sqrt3\,k d_p\), or \(k\rho_{\perp p}=\sqrt6\,k d_p\) if the perpendicular thermal gyroradius is intended. State which definition is used.
- Keep the positive proton gyrofrequency convention. Electron-reference conventions are discussed in `practices.md`, but a signed reference frequency can change what “upper half-plane” means physically. That extension needs its own convention check.
- Include negative and positive real frequencies explicitly. Passing through the imaginary axis above zero is supported. Linear edge parameterization handles the sign change; the existing positive-only logarithmic real map is unsuitable for an edge spanning zero.
- A positive \(\gamma_{\min}\) intentionally excludes marginal modes and all \(0<\gamma\le\gamma_{\min}\). Decrease it systematically only while response accuracy and boundary separation remain adequate. Never call a bounded thresholded zero count an unconditional stability proof.
- Lateral/top limits must be justified from reference root searches and expansion checks. Proton cyclotron frequencies, Doppler shifts, thermal resonance widths, and ordinary dispersion branches provide search scales; none alone bounds every unstable root. Electron and electromagnetic branches are also present in the full tensor.
- At sufficiently large \(x\) the exact warm scaled tensor tends toward its vacuum/high-frequency form, with leading \(D_{\rm ALPS}\sim a^6x^6\). Reaching a regime where this dominance is quantitatively valid requires comparison with both species plasma/gyrofrequency and wavevector scales, and control of cancellation in the numerical evaluation. An ion-scale finite top edge is not automatically an asymptotic edge. No universal finite “large enough” boundary is established here.
- A top/lateral expansion that leaves the count unchanged is useful evidence of coverage, but a finite set of expansions does not prove that no remote instability exists. Record the tested region.

If a zero or a possible pole is near/on an edge, refine and inspect the evaluation quality. If phase remains unreliable, mark the result unresolved. A deliberate boundary shift must be recorded as a **different contour**, with the resulting count conditional on those new bounds. Do not silently shift the user's region or treat displaced results as the same domain.

### 6.2 Precision and phase reliability

The implementation uses `double precision`/`double complex` and `MPI_DOUBLE_COMPLEX`, not arbitrary precision. The configured GNU/MPI build uses `mpifort`; [src/Makefile.am:4](/home/kgklein/Codes/ALPS_git/ALPS/src/Makefile.am:4) specifies `-O4 -fbounds-check -ffpe-trap=overflow` among its flags. IEEE binary64 is the effective precision of the existing executable, but portable Fortran declarations alone are not a cross-platform precision guarantee.

The determinant is an unscaled explicit polynomial in `wave`:
\[
D=A_{11}(A_{22}A_{33}+A_{23}^2)
 +2A_{12}A_{23}A_{13}-A_{13}^2A_{22}+A_{12}^2A_{33}.
\]
For the matrix symmetries in the source this is the ordinary determinant; it does not use \(|D|\), conjugation, a square root, an absolute-value normalization, or a root-solver deflation factor.

Risks and useful diagnostics:

- Cubic products can overflow/underflow, and large terms can nearly cancel. A finite nonzero result alone does not prove a reliable phase.
- Small determinant magnitude can be physically correct near a root, or result from a poorly scaled tensor; an absolute threshold inherited from `D_threshold` cannot distinguish them. Conversely, an enormous \(|D|\) does not imply poor phase if the tensor is accurately evaluated and well separated from singularity.
- A dimensionless diagnostic such as \(|\det A|/\|A\|^3\), computed with safe scaling, or a reciprocal condition estimate can flag difficult evaluations. Its interpretation must account for roundoff and susceptibility error; it is not a universal root-distance measure. Rank 0 can inspect its existing `wave` immediately after a call without altering susceptibility calculations.
- A determinant perturbation is amplified near a singular matrix. In a nonsingular case, \(\delta D/D\) is governed to first order by \(\operatorname{tr}(A^{-1}\delta A)\). Close to a zero, both susceptibility and determinant assembly errors can dominate phase.
- `chi0` diagnostic divisions can overflow at very small \(|x|\) even when the scaled determinant remains representable. The strict positive lower edge avoids exact zero but does not give unlimited precision as \(\gamma_{\min}\to0\).
- Existing `map_search` applies log magnitude, sentinel substitution, and limited-format output (`ALPS_fns.f90` 3806–3829). Those values/files are useful for ordinary root location, **not** as a phase oracle for the new driver. Use the complex scalar returned directly by `disp`.
- If actual contour evaluations exhaust dynamic range, a later optional phase/log-magnitude evaluation of a safely scaled tensor or a pivoted factorization can be considered using existing numerical facilities. It must preserve determinant sign/phase and report failures, and leave legacy evaluation unchanged. It is not required or implemented by this audit. Scaling a scalar after it has already overflowed cannot repair it.

For phase accumulation, work with safely normalized complex samples and an \(\operatorname{atan2}\)-based phase difference; avoid directly forming a potentially overflowing ratio/product of unscaled determinants. Multiplication of each **valid boundary value** by a positive real scale preserves its phase and curve winding. Such a scale need not be analytic, so this is a boundary arithmetic device, not a new analytic dispersion determinant. For a genuine analytic normalization \(g(x)\), the correction is \(W[gD]=W[D]+W[g]\); never assume a zero correction merely because roots away from \(g\)'s singularities agree.

The numerical stopping contract must distinguish:

1. **Converged numerical count:** finite reliable samples, a closed contour separated from zero within the adopted error assessment, resolved phase increments, stable winding under stricter sampling, and an accumulated phase sufficiently close to an integer multiple of \(2\pi\).
2. **Unresolved:** near/on-boundary zero, invalid or unreliable values, persistent phase ambiguity, response discontinuity/error, noninteger residual beyond tolerance, or refinement/evaluation/precision limits reached.

Small endpoint phase increments and a near-integer sum are insufficient by themselves: an entire loop can occur between samples and integer aliasing can hide missed winding. Use initial edge coverage, midpoint/refinement comparisons, tighter repeats and independent seed grids to detect this; assess a maximum unresolved segment scale. Without derivative/error bounds, finite adaptive samples provide numerical convergence evidence, not a rigorous exclusion theorem for arbitrary scalar functions. Report that distinction. For the warm pole-free domain, a negative converged count is a diagnostic failure requiring investigation, not a negative number of modes.

### 6.3 Equilibria expressed with current inputs

Use the current analytic bi-Maxwellian backend for **both** species: `use_bM=T, relat=F`. A Maxwellian represented as a numerical table still takes the numerical integration path and is outside initial physical validation.

| Input | Unstable candidate | Stable control |
|---|---|---|
| `&system` | `nspec=2`, proton reference units, `n_scan=0`; prescribed nonzero \(k_\parallel\), small \(k_\perp\) | Same wavevector and global normalization |
| Proton `&spec_1` | `nn=1,qq=1,mm=1,use_bM=T,relat=F` | Same |
| Proton `&bM_spec_1` | `bM_betas=3,bM_alphas=2,bM_pdrifts=0` | `bM_betas=3,bM_alphas=1,bM_pdrifts=0` |
| Electron `&spec_2` | `nn=1,qq=-1,mm≈1/1836.15267343,use_bM=T,relat=F` | Same |
| Electron `&bM_spec_2` | `bM_betas=3,bM_alphas=1,bM_pdrifts=0` | Same |
| Harmonic controls | Explicit `bM_nmaxs`, `bM_Bessel_zeros` for each species | Same values and convergence checks |

Equal electron/proton densities and \(T_e=T_{\parallel p}\) imply **\(\beta_{\parallel e}=3\)**, not a mass-ratio-scaled beta. Drift inputs are momenta in \(m_pv_A\) units; both zero yields charge neutrality and zero relative/net drift. Thermal speed follows the square root of inverse mass automatically.

Use the current [tests/test_bimax.in](/home/kgklein/Codes/ALPS_git/ALPS/tests/test_bimax.in) as an input-format reference, with explicit changed beta/anisotropy and scans disabled. Some older local inputs, for example `inputs/AC_test_bimax.in`, contain controls no longer accepted in `&system`; do not copy them verbatim into new validation cases. Use the reader's actual `bM_Bessel_zeros` spelling (`ALPS_io.f90` 761–767), despite the abbreviated label in `input.md`.

The finite-temperature, isotropic, drift-free pair is the intended stable physical control. The small root calculation below checks particular damped modes; it is not an exhaustive numerical demonstration of stability of the implemented determinant throughout an unbounded half-plane.

### 6.4 Small ordinary-solver check performed in A0

Used the existing `src/ALPS` executable without rebuilding or editing it. SHA-256:
`49b318f2ab0b2824a3d168b1eb4d66366a8ecc0435a052482e73a711d2afa298`.

All inputs, logs and solver outputs were isolated in `/tmp/alps_nyquist_a0/`, with local `solution/` and `distribution/` directories. The initial sandboxed MPI launch could not open local PMIx sockets; the same local calculation completed with sandbox escalation. No production output was overwritten. The existing executable's build provenance was not independently reconstructed; results are a small runtime check, not a newly built regression baseline.

Shared settings: \(k_\perp d_p=0.001\), \(|k_\parallel|d_p=0.3\), \(\texttt{vA}=10^{-4}\), both species warm and analytic, mass ratio \(m_e/m_p=5.44617021487\times10^{-4}\), harmonic cap 500, Bessel threshold \(10^{-50}\), `secant_method=0`, `numiter=100`, and initial guesses \(\pm0.5+0.08i\). Ordinary grid allocations used `nperp=npar=4`; these are unused zero placeholders for the analytic response, **not** a four-cell VDF approximation. Four MPI ranks were used.

| Case | \(k_\parallel d_p\) | Reference roots \(x\), as printed by the existing solver | Final root status / \(|D|\) |
|---|---:|---|---|
| \(\beta_{\parallel p}=3,\alpha_p=2;\ \beta_e=3,\alpha_e=1\) | \(+0.3\) | \(+0.39666+0.085986i\), \(-0.39666+0.085986i\) | Both converged in 10 iterations; \(|D|\approx1.02\times10^{-16},1.12\times10^{-16}\); tolerance \(10^{-15}\). |
| Same, isotropic protons \(\alpha_p=1\) | \(+0.3\) | \(+0.40460-0.0018543i\), \(-0.40460-0.0018543i\) | Converged in 12/11 iterations; \(|D|\approx3.58\times10^{-14},5.41\times10^{-14}\); tolerance \(10^{-13}\). |
| Same unstable equilibrium | \(-0.3\) | Same growing pair to printed precision | Both converged in 10 iterations; same printed determinants as \(+0.3\). |

Exploratory overly tight absolute determinant tolerances reached 100 iterations without satisfying the requested residual despite similar printed frequencies. The final checks above used tolerances above the observed residual fluctuations. This is evidence against equating a very small arbitrary `D_threshold` with extra physical accuracy; the residual is normalization dependent and the supplied root output has limited digits.

The angle is approximately \(0.191^\circ\), so this is nearly parallel. Finding the growing pair at this prescribed \(k\) supplies a concrete starting wavevector within an unstable region, consistent with the desired proton-cyclotron test. Polarization/mode tracking and uniqueness of each root remain to be checked; this audit does not identify every growing mode solely from those two solves.

The following is the complete temporary unstable input, reproduced here for review/reproducibility, not added as a production test:

```fortran
&system
 kperp=1.d-3, kpar=0.3d0, nspec=2, nroots=2,
 use_map=F, writeOut=T, nperp=4, npar=4, vA=1.d-4,
 arrayName='audit_bimax', Bessel_zero=1.d-45,
 secant_method=0, numiter=100, D_threshold=1.d-15,
 D_prec=1.d-5, D_gap=1.d-5, D_tol=1.d-10,
 positions_principal=3, n_resonance_interval=100,
 Tlim=0.01d0, determine_minima=F, scan_option=1, n_scan=0
/
&guess_1
 g_om=0.5d0, g_gam=0.08d0
/
&guess_2
 g_om=-0.5d0, g_gam=0.08d0
/
&spec_1
 nn=1.d0, qq=1.d0, mm=1.d0, AC_method=1, ff=1,
 relat=F, log_fit=T, use_bM=T
/
&bM_spec_1
 bM_nmaxs=500, bM_Bessel_zeros=1.d-50, bM_betas=3.d0,
 bM_alphas=2.d0, bM_pdrifts=0.d0
/
&spec_2
 nn=1.d0, qq=-1.d0, mm=5.44617021487d-4, AC_method=1, ff=1,
 relat=F, log_fit=T, use_bM=T
/
&bM_spec_2
 bM_nmaxs=500, bM_Bessel_zeros=1.d-50, bM_betas=3.d0,
 bM_alphas=1.d0, bM_pdrifts=0.d0
/
```

For the stable run change only proton `bM_alphas` to 1 and `D_threshold` to \(10^{-13}\); for the reverse run change only `kpar` to \(-0.3\) in the unstable input. With directories created and the working directory set to the temporary directory, the successful local command was:

```sh
mpirun --oversubscribe -np 4 /home/kgklein/Codes/ALPS_git/ALPS/src/ALPS unstable.in
```

`solution/unstable.roots`, `solution/stable.roots`, `solution/reverse.roots`, their inputs, and `*.stdout` remain in that temporary directory for this session; temporary storage is not a durable test archive. The report retains the input and printed results.

### 6.5 Minimal next physical validation sequence

1. Re-establish the reference roots with an identified build, independent starting guesses, and direct residual checks; inspect eigenfunctions/polarization with existing facilities if needed to label the ion-cyclotron branch. Vary analytic harmonic controls to establish response/root convergence. A root solver finding two roots cannot establish completeness.
2. Use ordinary `use_map=T`, `loggridw=F`, `loggridg=F` searches spanning both real-frequency signs and positive growth to look for additional roots at the chosen wavevector. Refine minima with existing routines. Use adequate grid resolution and overlapping/expanded searches: map minima and `nroots/numroots` limits can miss roots, so map output is reference evidence, not a certified count. Keep all ranks in every evaluation.
3. In A1's scalar tests verify orientation, zero-minus-pole sign, multiplicity, branch wrapping and unresolved outcomes. In A2 first use local contours around a well isolated reference root (and displaced contours excluding it) after root isolation; then use a two-sided rectangle covering the measured unstable region. Counts are acceptance criteria only after boundary separation and reference coverage are established.
4. Set \(\gamma_{\min}\) below measured growing rates with a margin exceeding their uncertainty and evaluation error. Set upper/lateral bounds beyond the reference-root region, with margins verified by boundary maps/evaluations. Expand those bounds and lower \(\gamma_{\min}\) in controlled repeats. No numerical contour limits are prescribed by A0.
5. Repeat the identical contour/wavevector protocol on the isotropic control. A resolved zero count is expected for the stable warm equilibrium, but it must be measured and qualified by the tested bounds.
6. Exercise \(k_\parallel<0\) and both signs of \(\operatorname{Re}\omega\) explicitly. For real physical fields, conjugate solutions are related between \((\omega,\mathbf k)\) and \((-\omega^*,-\mathbf k)\); the no-drift, even parallel distributions in this test supply additional symmetry. The runtime comparison supports that symmetry here. With drifting/asymmetric distributions, do not impose same-\(k\) symmetry or simply multiply a positive-frequency count by two. At fixed \(\mathbf k\), distinct signed-frequency zeros are counted separately.
7. Probe a small neighborhood in \(k_\parallel\) around the selected unstable point using the existing ordinary workflow if needed to establish branch continuity and avoid a marginal endpoint. This is selection of reference cases, not authorization for a Nyquist wavevector scan. Subsequent Nyquist scans must retain the same collective harmonic decomposition.

## 7. Risks and unresolved questions

### Confirmed implementation facts/limitations

| Finding | Consequence / disposition |
|---|---|
| Returned determinant differs from \(\det M\) by \(a^6x^6\) | Safe winding equivalence only on contours with zero index about the origin; initial positive lower boundary satisfies this. Retain raw `disp`. |
| Unconditional \(1/x^2\) susceptibility snapshots | Exact \(x=0\) is unsafe even if a scaled tensor limit exists. Exclude it. |
| Warm \(k_\parallel=0\) divisions | Do not use an exactly perpendicular wavevector for initial validation. |
| Cold rational resonance denominators | Real-axis poles/limits need a separate analysis; not initial warm validation. |
| Numerical-VDF moving windows and `Tlim` approximation | The computed backend is not certified analytic merely because the physical response is. Later boundary continuity/response convergence is required. |
| `WOFZ` overflow flag ignored by caller | No complete special-function failure status currently reaches `disp`. Extreme-argument robustness is incomplete. |
| `disp` is collective, mutates root-only tensors, and does not broadcast its input | New driver must broadcast requests/frequencies and coordinate termination; no independent rank sampling. |
| Existing map output transforms/clips values | Use direct complex determinant returns for winding. |
| No common equilibrium teardown; hard error paths abort MPI | Reuse program finalization; new module cleans its own storage. Do not claim all failures are recoverable. |
| Rank-0 `sproc` sentinel is not explicitly initialized before startup's Bessel guard | Existing run succeeds, but future startup/adapter guards should use defined rank state; record this portability issue without changing the harmonic decomposition. |
| Source comments mention \(k_\perp^2\) scaling absent from active expressions | Use the derived active tensor; resolve documentation discrepancy in a later reviewed change if desired. No A0 code edit. |

These do not demonstrate an incorrect initial warm-species winding calculation; no such calculation has been implemented or run.

### Potential concerns requiring evidence

- How close can the lower edge approach marginal roots while preserving reliable susceptibility and determinant phase? There is no measured universal \(\gamma_{\min}\).
- What lateral/top bounds establish the desired physical coverage at a selected \(k\)? Two located roots and a few finite expansions are not a global exclusion theorem.
- What response/determinant error margin is sufficient near a boundary root, including `WOFZ` region joins, large-\(\zeta\) cancellation, Bessel approximations, and harmonic truncation? The special-function accuracy comment is not a determinant error certificate.
- Does the raw determinant remain representable and sufficiently conditioned over the selected contour? No phase/log-determinant extension should be added before a demonstrated need.
- Do numerical-VDF cell crossings, edge returns and `Tlim` switches create measurable boundary jumps or spurious winding at practical resolutions? Their existence is confirmed; their quantitative harm has not been measured. Include both nonrelativistic and relativistic paths when that stage is reached.
- How should recoverable backend failure be reported collectively without changing legacy behavior? A finite-value adapter is a minimum check, not full propagation of special-function/status information.
- Are multiple roots or nearly coincident branches present in the reference rectangle? Winding counts multiplicity, whereas a root list often deduplicates positions.
- Can sampling limits expose missed fast loops rather than return an aliased integer? Numerical stability under refinement must be part of the convergence contract.
- Existing executable results were not regenerated from a clean, independently identified build. They establish a practical candidate case, not preserved regression baselines.
- Electron-reference normalization, real-axis contours, exactly perpendicular warm propagation, and cold pole subtraction are outside the initial supported mathematical contract.

### Relationship to the supplied papers

[Klein et al. (2017), JGR, §2, pp. 9816–9818](https://doi.org/10.1002/2017JA024486) describes the historical inverse-determinant Nyquist diagram with handed real-axis crossings and a large-frequency contribution. Its notation \(|D|\) denotes a determinant, not the positive real modulus used in ordinary \(\log|D|\) diagnostics. That presentation is historical context, not a prescription to integrate the positive scalar \(|D|^{-1}\), which contains no complex phase.

[Klein et al. (2018), PRL, p. 205102-2](https://doi.org/10.1103/PhysRevLett.120.205102) discusses counting modes above a minimum growth threshold. Its observational electron temperature is based on the scalar proton temperature; the present requested test explicitly uses \(T_e=T_{\parallel p}\) instead. Its \(k\rho_p\) conventions and illustrated counts must not be transplanted to the ALPS test.

For a closed counterclockwise curve, \(W[1/D]=-W[D]\); orientation and inverse choice must therefore be tracked explicitly when comparing conventions. A contour integral of \(1/D\) alone is not generally a zero count: the argument principle uses \(D'/D\), or equivalently the continuous phase change of \(D\). No historical crossing algorithm, asymptotic half-crossing correction, or numerical PLUMAGE reproduction is needed for the proposed implementation.

## 8. Ordered implementation checklist

All items below are recommendations for subsequent authorized stages; none has been implemented in A0.

1. **Freeze the initial mathematical contract.** Proton reference units; warm positive-beta analytic species for physical validation; real fixed \(k_\parallel\ne0\); rectangle entirely at \(\operatorname{Im}\omega>0\); use \(D_{\rm ALPS}=\det A\) directly; count only the stated finite region with multiplicity and no enclosed susceptibility poles.
2. **Define the scalar winding interface and statuses (A1).** A separate Fortran module with no susceptibility formulas; ordered closed-contour samples, independent evaluation adapter, explicit converged/unresolved outcomes, no integer result on failure.
3. **Validate the scalar algorithm independently (A1).** Known analytic/rational examples with zeros and poles, reversed orientation, repeated zeros, loops missed by coarse sampling, and a zero on/near an edge. Verify sign conventions and that sample/refinement caps lead to unresolved. These mathematical checks require no physical VDF cases.
4. **Add narrow opt-in input guards.** Default-off selector, explicit four-bound `&nyquist` group, documented sampling controls, validation of finite/ordered bounds and supported wavevector/species settings, `n_scan=0`, no required guesses/maps. Define neutral shared bookkeeping without changing normal input behavior.
5. **Share the complete current startup.** Equilibrium/model construction and transfer, harmonic limits/splitting, and Bessel setup as mapped in §2. Ensure the worker Bessel guard uses defined rank state rather than relying on an uninitialized rank-0 `sproc`. No per-frequency reinitialization.
6. **Add the all-rank driver and rank-0 request protocol (A2).** Command broadcast, complex frequency broadcast, one collective existing `disp` per request, rank-0 sample storage/refinement/phase, coordinated termination. Preserve the two reductions, determinant broadcast and barrier; add no frequency/wavevector decomposition.
7. **Dispatch after initialization and bypass root/scan execution.** Keep normal `map_search/refine_guess/om_scan/om_double_scan` code unchanged. Reuse MPI finalization/log closing. Register only necessary module build/dependency changes.
8. **Implement evaluation-quality checks before accepting counts.** Reject invalid/zero/ambiguous values; inspect available tensor scaling/conditioning when needed; distinguish backend failure from unresolved phase and from zero modes. Decide whether explicit special-function failure propagation is needed before claiming robust convergence. Do not alter default determinant arithmetic or legacy error behavior incidentally.
9. **Establish reference-root coverage for the supplied physical cases.** Reproduce §6's candidate roots with identified build, response convergence, additional guesses/maps and signed-frequency coverage; select documented contour boundaries from evidence.
10. **Validate contour convergence and physical controls.** Local inclusion/exclusion, two-sided unstable rectangle, isotropic control, negative \(k_\parallel\), stricter/independent sampling and boundary expansion/lower-edge studies. Record bounds and status for every result; do not assume the integer from the two current root solves.
11. **Check legacy preservation.** Run existing relevant bi-Maxwellian/root/map checks with Nyquist omitted/default-off; compare outputs using established numerical tolerances. Check that collective ordering and rank-0 diagnostics remain correct. No broad refactor or extensive unrelated test campaign.
12. **Gate numerical-VDF physical interpretation separately.** Reuse the same adapter, but study grid/velocity-support/harmonic/threshold convergence and continuity across resonance switches against analytic bi-Maxwellian references before certifying numerical-VDF counts. Keep future wavevector scans and any deeper precision/pole treatment outside the initial single-wavevector implementation.

A0 deliverable verification: production `src/` had no working-tree diff before or after the audit; existing unrelated working-tree changes were left in place. The sole repository deliverable created by this task is this Markdown report.
