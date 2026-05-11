# FastEarth.jl — Design Document

**Status:** initial design, pre-implementation
**Date:** 2026-05-11
**Scope:** global 2D Glacial Isostatic Adjustment (GIA) model, generalising FastIsostasy
(Swierczek-Jereczek et al. 2024) from a regional Cartesian domain to the whole sphere,
with gravitationally self-consistent sea level and an interface suitable for interactive
coupling to ice-sheet models at 10–50 km resolution.

---

## 1. Goals and design constraints

1. **Global domain.** Whole sphere; no regional boundary conditions.
2. **2D effective viscosity (layer-stacking).** Retain the core FastIsostasy simplification:
   the 3D radial viscosity profile $\eta(r;\theta,\phi)$ is collapsed into a laterally
   variable 2D field $\eta_{\text{eff}}(\theta,\phi,\ell)$ by analytic layer-stacking
   (Cathles 1975; Swierczek-Jereczek et al. 2024, §2.3). The dimensional reduction is
   radial, not lateral, so the philosophy carries over to the sphere directly.
3. **Gravitationally self-consistent sea level.** Ocean loading must respond to the
   evolving geoid; rotational feedback included at degree-2.
4. **Ice-sheet coupling target.** Interactive coupling at $\Delta t_{\text{couple}} \in
   [10, 100]$ yr on grids of 10–50 km equivalent resolution, over paleo and
   future-projection timescales ($10^2$–$10^5$ yr). Speed matters; ~1% accuracy vs.
   3D viscoelastic codes (SELEN, ISSM-SESAW, VILMA, CALSEA) is the target.
5. **Implementation language:** Julia.

Non-goals (out of scope, at least initially): power-law / Burgers rheology,
3D viscosity solvers, full normal-mode / propagator-matrix machinery, compressibility
corrections, ellipsoidal-Earth or precession effects.

---

## 2. Governing equations

### 2.1 Notation

| Symbol | Meaning |
|---|---|
| $(\theta,\phi)$ | colatitude, longitude on sphere of radius $R$ |
| $t$ | time |
| $u(\theta,\phi,t)$ | radial bedrock displacement [m] |
| $u_{\text{visc}}, u_{\text{el}}$ | viscous and elastic parts; $u = u_{\text{visc}} + u_{\text{el}}$ |
| $N(\theta,\phi,t)$ | geoid perturbation [m] |
| $S(\theta,\phi,t)$ | relative sea level [m] |
| $C(\theta,\phi,t)$ | ocean function $\in\{0,1\}$ |
| $H_i, H_w$ | grounded ice thickness, ocean column thickness [m] |
| $L = \rho_i H_i + \rho_w H_w$ | surface mass load [kg m⁻²] |
| $\eta_{\text{eff}}(\theta,\phi,\ell)$ | layer-stacked, wavelength-dependent effective viscosity [Pa s] |
| $D(\theta,\phi) = E T_e^3 / [12(1-\nu^2)]$ | flexural rigidity [N m] |
| $\rho_m, g$ | mantle density, surface gravity |
| $\Delta_S$ | surface Laplace–Beltrami operator |
| $\lvert\nabla_S\rvert$ | spherical half-derivative pseudo-differential operator |
| $h_\ell^L, k_\ell^L$ | load Love numbers (elastic) |

### 2.2 Core PDE (spherical LV-ELVA)

FastIsostasy's regional ELVA equation

$$2\eta_{\text{eff}}(x,y)\,\lvert\nabla\rvert\,\partial_t u + \rho_m g\, u + D\,\nabla^4 u = \sigma^{zz}$$

uses the **half-derivative operator** $\lvert\nabla\rvert = \sqrt{-\Delta}$ (in Fourier
space: $\lvert k\rvert$) — the viscous half-space normal-stress kernel of Cathles (1975).
This is the FastIsostasy heritage operator and we keep it.

Its spherical analogue uses the surface half-derivative $\lvert\nabla_S\rvert$, defined
through its action on spherical harmonics:

$$\lvert\nabla_S\rvert\, Y_{\ell m} = \frac{\sqrt{\ell(\ell+1)}}{R}\, Y_{\ell m}, \qquad
\Delta_S\, Y_{\ell m} = -\frac{\ell(\ell+1)}{R^2}\, Y_{\ell m}.$$

The global LV-ELVA equation is then

$$\boxed{\;
2\,\eta_{\text{eff}}(\theta,\phi,\ell)\,\lvert\nabla_S\rvert\,\partial_t u_{\text{visc}}
\;+\; \rho_m g\, u
\;+\; \Delta_S\!\bigl(D(\theta,\phi)\,\Delta_S u\bigr)
\;=\; -g\bigl(L + \rho_m N\bigr).
\;}\tag{FE-1}$$

Notes on each term:

- **Viscous term.** $\lvert\nabla_S\rvert$ replaces FastIsostasy's $\lvert\nabla\rvert$.
  Modal relaxation time at constant viscosity is
  $\tau_\ell = 2\eta\sqrt{\ell(\ell+1)}/R\,/\,[\rho_m g + D\,\ell^2(\ell+1)^2/R^4]$,
  recovering the Haskell-style spectrum — short wavelengths relax faster, long
  wavelengths slower, no ad-hoc relaxation times.
- **Effective viscosity.** $\eta_{\text{eff}}$ is **wavelength-dependent**: short loads
  probe the asthenosphere, long loads probe the lower mantle. Layer-stacking provides
  $\eta_{\text{eff}}(\theta,\phi,\ell)$ from a 1D column integral at each grid point
  using a depth set by $\ell$. Applied per-degree in spectral space.
- **Buoyancy.** $\rho_m g u$ provides the Winkler restoring force, dominant at long
  wavelengths.
- **Flexure.** Self-adjoint form $\Delta_S(D\Delta_S u)$ to handle laterally variable
  $D(\theta,\phi)$. Membrane / thin-shell curvature terms ($\sim E T_e/[R^2(1-\nu^2)]$)
  are dropped — they are $\lesssim 1\%$ at the wavelengths that matter for ice-sheet
  loading, and the Winkler term already provides finite stiffness at $\ell = 0, 1$.
- **Forcing.** $-g(L + \rho_m N)$ is the effective load: surface mass plus the
  $\rho_m$-equivalent column representing geoid-uplift-induced buoyancy
  ("self-gravity" coupling).

### 2.3 Elastic displacement (diagnostic)

The viscous PDE (FE-1) is first-order in time and cannot capture instantaneous elastic
response. We add an elastic displacement diagnostically, applied per-mode in spectral
space:

$$u_{\text{el},\ell m}(t) = \frac{3\, h_\ell^L}{\rho_E (2\ell + 1)}\, L_{\ell m}(t),
\qquad u = u_{\text{visc}} + u_{\text{el}}.\tag{FE-2}$$

This is a one-line spectral multiplication with tabulated elastic load Love numbers
$h_\ell^L$ (e.g. PREM-based, Farrell 1972). Cost is negligible.

### 2.4 Self-gravity and the sea-level equation

We use a **kinematic** self-gravity formulation (no viscoelastic Love-number
convolution): the geoid is the Newtonian potential of the surface load plus the
$\rho_m$-equivalent column from the deformed bedrock, with elastic load Love-number
boosting:

$$\boxed{\;N_{\ell m}(t) = \frac{3}{\rho_E(2\ell+1)}\Bigl[(1 + k_\ell^L)\,L_{\ell m}(t)
\;+\; \rho_m\, u_{\text{visc},\ell m}(t)\Bigr].\;}\tag{FE-3}$$

Two SLE solvers are supported (selectable at runtime — see §4.4):

**SLE-A (production, "emulate SELEN within ~1%"):** Kendall et al. (2005) ocean-function
iteration. Per coupling step:

$$S(\theta,\phi,t) = \bigl[N(\theta,\phi,t) - u(\theta,\phi,t) + \Delta\Phi(t)/g\bigr]\,C(\theta,\phi,t),
\tag{FE-4a}$$

$$\Delta\Phi(t)/g \;=\; -\frac{(M_{\text{ice}}(t)-M_{\text{ice}}(t_0))/\rho_w + \int_\Omega (N-u)\,C\,dA}
{\int_\Omega C\,dA},\tag{FE-4b}$$

iterated to convergence on $(S, C)$, with topography $T = T_0 + u - N$ used to update
$C$ each inner sweep.

**SLE-B (fast, "qualitatively self-consistent"):** tabulated single-pass approximation

$$N_\ell \approx \alpha_\ell u_\ell + \beta_\ell L_\ell, \qquad S = N - u + c(t),\tag{FE-5}$$

where $c(t)$ is the spatially uniform mass-conservation constant and $(\alpha_\ell,
\beta_\ell)$ are precomputed coefficients consistent with (FE-3). Skips the
ocean-function iteration; suitable for fast ensemble runs where the spatial pattern
of $S$ matters more than $\sim 1\%$ pointwise accuracy.

### 2.5 Rotational feedback (degree-2)

Polar motion $\mathbf{m}(t)$ is computed from the degree-2 surface-load inertia tensor
via the linearised Liouville equation (Milne & Mitrovica 1998), and adds a degree-2
contribution to the centrifugal potential which feeds back into $N_{2m}$. Two extra
ODEs per timestep; negligible cost. Always on.

### 2.6 Closure on load and topography

$$L = \rho_i H_i + \rho_w H_w, \qquad H_w = \max(0, -T)\cdot C, \qquad T = T_0 + u - N,$$

where $T_0$ is the prescribed reference bedrock topography (e.g. ETOPO / BedMachine
restricted to ice-free conditions). $C$ is updated from $T$ and the grounded-ice mask
inside the SLE iteration (SLE-A only).

### 2.7 Assumptions stated explicitly

- Newtonian linear rheology, incompressible mantle.
- Laterally variable $\eta_{\text{eff}}(\theta,\phi,\ell)$, $D(\theta,\phi)$.
- No 3D viscosity solver; layer-stacking is the dimensional reduction.
- Self-gravity is kinematic (no viscoelastic Love-number time-convolution); the
  viscous response is carried entirely by $u_{\text{visc}}$.
- Elastic response added diagnostically via $h_\ell^L$.
- Rotational feedback at degree-2 only; no ellipsoidal-Earth or precession.
- No spatial boundary conditions; only regularity at the poles (automatic in the
  SH basis).

---

## 3. Numerical strategy

### 3.1 Spatial discretisation

**Pseudo-spectral spherical harmonics on a Gauss–Legendre grid**, truncation
$\ell_{\max}$ tunable (default target 512 for $\sim$40 km equatorial resolution).

Why SH wins for this problem:

- $\Delta_S$, $\Delta_S^2$, $\lvert\nabla_S\rvert$ all **diagonal** in $(\ell,m)$.
- The SLE Love-number convolutions in (FE-3) are **diagonal**. This is the killer
  feature — every alternative grid pays extra cost here.
- No coordinate singularity at the poles.

Laterally variable coefficients ($\eta_{\text{eff}}, D, C$) are handled
pseudospectrally: synthesise the operand to the Gauss–Legendre grid, multiply
pointwise, transform back. 2/3 dealiasing rule applied to the variable-coefficient
products. Smooth ocean function via a brief spectral filter at $\ell \approx 0.9\,
\ell_{\max}$ to control Gibbs ringing.

Rejected alternatives and why:

| Grid | Why rejected |
|---|---|
| Lat-lon FD | Pole singularity in $\Delta_S^2$; CFL collapse; SLE convolution becomes dense $O(N^2)$. |
| Cubed sphere (ClimaCore) | Excellent for hyperbolic systems, but SLE wants SH — extra round-trip negates the win. |
| Icosahedral / geodesic | Same critique as cubed sphere here. |
| HEALPix | Equal-area is attractive for ocean integrals, but thinner Julia ecosystem and we'd rebuild SH machinery anyway. Possible future swap if `Healpix.jl` matures further. |

### 3.2 Time integration

Equation (FE-1) is **stiff**: modal relaxation times $\tau_\ell$ span $\sim 10^2$ yr
(short wavelengths, asthenosphere) to $\sim 10^5$ yr (long wavelengths, deep mantle).
Explicit treatment is infeasible.

**IMEX with implicit homogeneous-Earth diagonal:**

- *Implicit:* the spectral-diagonal linear operator at a reference viscosity $\bar\eta$
  and rigidity $\bar D$. In SH coefficients this is a per-mode scalar solve.
- *Explicit:* the lateral-coefficient residuals
  $(\eta_{\text{eff}} - \bar\eta)\lvert\nabla_S\rvert\partial_t u$,
  $\Delta_S((D-\bar D)\Delta_S u)$, evaluated pseudospectrally; and the SLE load
  update.

Preconditioner / per-mode solve (from the integrated-proposal synthesis):

$$P^{-1}_\ell \;=\; \Bigl[\rho_m g \;+\; 2\bar\eta\, k_\ell/\Delta t \;+\; \bar D\, k_\ell^4\Bigr]^{-1},
\quad k_\ell = \sqrt{\ell(\ell+1)}/R.\tag{FE-6}$$

Default integrator: an exponential time differencing scheme (ETD1) applied to the
diagonal part — exact for the linear stiff piece, unconditionally stable, no Krylov
inner iteration when lateral contrasts are modest. Fallback: preconditioned GMRES with
$P^{-1}_\ell$ as preconditioner, for the case of strong lateral viscosity contrasts.

Typical timestep: $\Delta t = 10$–$100$ yr (coupling-driven, not stability-driven).
SLE iteration nested inside each step.

### 3.3 SLE inner solver

Per outer step, SLE-A runs the Kendall et al. (2005) fixed-point loop:

1. From current $u, H_i$: compute load $L^{(k)}$ on grid.
2. Spectral transform: apply (FE-3) for $N^{(k)}_{\ell m}$.
3. Inverse transform; update $S^{(k+1)}$ from (FE-4a) with $\Delta\Phi$ from (FE-4b).
4. Update $H_w$ and ocean function $C$ from topography $T = T_0 + u - N$.
5. Repeat until $\|S^{(k+1)} - S^{(k)}\|_\infty < \varepsilon$ (typically 3–6
   iterations; Anderson acceleration via `NLsolve.jl` if needed).
6. Outer rotation update (degree-2 Liouville), re-sweep 1–2 times.

SLE-B applies (FE-5) once per step with $c(t)$ from global mass conservation.

### 3.4 Ice-sheet coupling interface

Contract is deliberately minimal so FastEarth stays ice-model-agnostic:

**Input (per coupling step):**
- $H_i(\theta,\phi)$ on the FastEarth Gauss–Legendre grid (regridded conservatively
  from the ice model's native grid by the coupler).

**Output:**
- $u(\theta,\phi)$ — bedrock displacement, hence bedrock topography $b = T_0 + u$.
- $N(\theta,\phi)$ — geoid (needed if the ice model uses geoid-referenced sea level).
- $S(\theta,\phi)$ — relative sea level (the field marine-ice-sheet grounding-line
  tests actually want).

Internal $\Delta t \le \Delta t_{\text{couple}}$ with checkpointable state for
restartability.

### 3.5 Julia package stack

| Layer | Package | Role |
|---|---|---|
| SH transforms (primary) | **FastTransforms.jl** | Synthesise/analyse on Gauss–Legendre grid |
| SH transforms (fallback) | **SHTns_jll** / `SHTnsSpheres.jl` | High-$\ell_{\max}$ performance / GPU |
| Reference SH+IMEX template | **SpeedyWeather.jl** | Borrow `SpectralTransform`, `LowerTriangularMatrix` patterns |
| Time integration | **OrdinaryDiffEq.jl** (`SplitODEProblem`, ETD, IMEX-RK) | Linear/nonlinear splitting |
| Implicit linear solves | **LinearSolve.jl** + **SciMLOperators.jl** | Per-mode diagonal operator wrapper |
| Fixed-point SLE | **NLsolve.jl** (Anderson) | Inner SLE iteration |
| Regridding ice ↔ FastEarth | **Interpolations.jl** + custom conservative remap | Coupling boundary |
| GPU portability | **KernelAbstractions.jl** | Grid-space pointwise ops |
| I/O | **NCDatasets.jl**, **TOML.jl** | Config and output |
| Diagnostics / plotting | **Makie.jl**, **GeoMakie.jl** | Maps on the sphere |

Deliberate non-uses: `ClimaCore.jl` (wrong grid primitive for this problem),
`Oceananigans.jl` (Cartesian-oriented).

---

## 4. Architecture and code organisation

### 4.1 Module structure

```
FastEarth.jl
├── src/
│   ├── FastEarth.jl              # top-level module, public API
│   ├── grid.jl                   # GaussLegendreGrid, SH transform plumbing
│   ├── operators.jl              # |∇_S|, Δ_S, Δ_S², variable-coefficient applications
│   ├── viscosity.jl              # layer-stacking, η_eff(θ,φ,ℓ) construction
│   ├── lithosphere.jl            # D(θ,φ), flexure operator
│   ├── elastic.jl                # u_el diagnostic from h_ℓ^L
│   ├── sle.jl                    # SLE-A (Kendall) and SLE-B (tabulated)
│   ├── rotation.jl               # degree-2 Liouville polar motion
│   ├── timestepping.jl           # ETD / IMEX integrators, preconditioner
│   ├── coupling.jl               # ice-sheet input/output interface
│   └── io.jl                     # NetCDF, TOML config
├── test/                         # unit tests, benchmarks vs. analytic / SELEN
├── docs/                         # this design doc, derivations, ADRs
├── benchmarks/                   # Farrell, disc-load, SELEN comparisons
└── Project.toml
```

### 4.2 Public API sketch

```julia
using FastEarth

grid = GaussLegendreGrid(ℓmax = 512)
earth = EarthModel(grid;
    η_profile = load_viscosity_field("data/viscosity_3d.nc"),  # → η_eff(θ,φ,ℓ)
    Te        = load_lithosphere_thickness("data/Te.nc"),
    T0        = load_topography("data/etopo.nc"),
    sle       = :kendall,         # or :tabulated
    rotation  = true,
)

state = initialise(earth; H_i_initial = zeros(grid))

# Coupling step
for t in 0:Δt:tmax
    H_i = get_ice_thickness_from_coupler(t)
    advance!(state, earth, H_i, Δt)
    push_to_coupler(state.u, state.N, state.S)
end
```

### 4.3 Type design notes

- `GaussLegendreGrid` carries Legendre nodes/weights, $\ell_{\max}$, and FastTransforms
  plans. One canonical grid object passed everywhere.
- Spectral state stored as `LowerTriangularMatrix`-style coefficient arrays (borrowed
  from SpeedyWeather idiom; one real-valued triangle per real scalar field).
- Operators (`|∇_S|`, `Δ_S`, `Δ_S²`, the preconditioner) are diagonal arrays of length
  $\ell_{\max}+1$, broadcast over the $m$-stride.
- `EarthModel` is immutable after construction; `state` holds the time-varying fields.

### 4.4 Configuration

TOML-driven with explicit sections:

```toml
[grid]
ℓmax = 512

[rheology]
viscosity_file = "data/viscosity_3d.nc"
Te_file        = "data/Te.nc"
wavelength_dependent_η = true   # use η_eff(θ,φ,ℓ); false → single η(θ,φ)

[sle]
solver         = "kendall"      # "kendall" | "tabulated"
tolerance      = 1.0e-4
max_iterations = 20

[rotation]
enabled = true

[time]
Δt_internal = 25.0              # years
integrator  = "ETD1"            # "ETD1" | "IMEX_GMRES"
```

---

## 5. Validation hierarchy

Progressive benchmarking before exposing the coupling interface. Each level gates the
next.

0. **Farrell (1972) elastic point-load Green's function.** Pure spectral test — no
   time integration, no SLE, only the SH transform plumbing and elastic Love-number
   tables. Verifies $u_{\text{el}}$ and the geoid response (FE-2)–(FE-3).

1. **Homogeneous spherical ELVA, analytic.** Constant $\eta$, constant $D$, single
   spherical-cap load. Each $(\ell,m)$ mode decouples; relaxation should match the
   analytic $\tau_\ell$ spectrum to numerical precision.

2. **Homogeneous viscoelastic — TABOO / SELEN.** Match a SELEN run for an idealised
   disc-load deglaciation on a layered radial-1D Earth. Target: 1% in peak uplift,
   spatial pattern, and far-field subsidence.

3. **Lateral viscosity anomalies.** Synthetic checkerboard or single-anomaly $\eta_{\text{eff}}$
   fields. No external reference exists, but cross-check between SLE-A and SLE-B and
   between ETD and IMEX-GMRES integrators.

4. **Regional comparison vs. FastIsostasy v1.0.** Restrict FastEarth to a hemisphere
   and compare against the original Cartesian FastIsostasy on a shared Antarctic
   deglaciation problem. Should agree to within the curvature corrections.

5. **ESMIP-style benchmark.** Match Spada et al. and/or the Martinec et al. (2018)
   SLE intercomparison cases, in both SLE-A and SLE-B modes.

6. **Coupled AIS / GIS deglaciation.** End-to-end run coupled to an ice-sheet model
   over a glacial cycle.

---

## 6. Known risks and open questions

1. **Spherical layer-stacking derivation.** FastIsostasy's $\eta_{\text{eff}}$
   derivation is Cartesian (Cathles 1975). The spherical analogue with wavelength
   $\to$ degree mapping is plausible at leading order in $h/R$ but should be written
   down explicitly. **Action:** dedicate a `docs/derivation_layer_stacking.md` note
   before benchmark level 2.

2. **Kinematic vs. viscoelastic Love-number self-gravity.** The framework in (FE-3)
   omits the viscoelastic time-convolution kernel that SELEN-class codes use. The
   claim is that $u_{\text{visc}}$ already carries the viscous gravitational signal.
   This should validate against level 2 / level 5. If it fails, the cheapest upgrade
   is to add a small number of viscous Love-number modes (Maxwell relaxation
   spectrum), not a full normal-mode solve.

3. **Variable-coefficient biharmonic stability.** Sharp $D$ contrasts at
   ocean–continent boundaries may produce Gibbs ringing in the pseudospectral
   evaluation of $\Delta_S(D \Delta_S u)$. Mitigation: smooth $D$ to the resolved
   scale and apply 2/3 dealiasing. Fallback: Schwarz-style preconditioned solve.

4. **Wavelength-dependent $\eta_{\text{eff}}(\theta,\phi,\ell)$ vs. a single field
   $\eta_{\text{eff}}(\theta,\phi)$.** Building the spectrally adaptive form from
   day one is the recommendation, but it doubles memory for the viscosity field.
   Whether the accuracy improvement justifies it should be settled empirically at
   benchmark level 3.

5. **Rotational-feedback coupling stability.** The degree-2 Liouville feedback is
   known to be only weakly contractive (Mitrovica et al. 2005). Aggressive timesteps
   may stall the SLE/rotation outer loop. Mitigation: under-relax the rotation update
   ($\omega \approx 0.5$), or solve SLE + rotation as a small dense coupled problem
   each iteration.

---

## 7. References

- Cathles, L. M. (1975). *The Viscosity of the Earth's Mantle.* Princeton UP.
- Farrell, W. E. (1972). Deformation of the Earth by surface loads. *Rev. Geophys.*
  10, 761–797.
- Farrell, W. E., & Clark, J. A. (1976). On postglacial sea level. *Geophys. J. R.
  astr. Soc.* 46, 647–667.
- Gomez, N., Mitrovica, J. X., Huybers, P., & Clark, P. U. (2010). Sea level as a
  stabilizing factor for marine-ice-sheet grounding lines. *Nat. Geosci.* 3, 850–853.
- Kendall, R. A., Mitrovica, J. X., & Milne, G. A. (2005). On post-glacial sea level
  II. *Geophys. J. Int.* 161, 679–706.
- Martinec, Z., et al. (2018). A benchmark study of numerical implementations of the
  sea level equation. *Geophys. J. Int.* 215, 389–414.
- Milne, G. A., & Mitrovica, J. X. (1998). Postglacial sea-level change on a rotating
  Earth. *Geophys. J. Int.* 133, 1–19.
- Mitrovica, J. X., & Peltier, W. R. (1991). On post-glacial geoid subsidence over
  the equatorial oceans. *J. Geophys. Res.* 96, 20053–20071.
- Mitrovica, J. X., Wahr, J., Matsuyama, I., & Paulson, A. (2005). The rotational
  stability of an ice-age Earth. *Geophys. J. Int.* 161, 491–506.
- Spada, G., & Stocchi, P. (2007). SELEN: a Fortran 90 program for solving the
  sea-level equation. *Comp. Geosci.* 33, 538–562.
- Swierczek-Jereczek, J., Robinson, A., Blasco, J., Alvarez-Solas, J., & Montoya, M.
  (2024). FastIsostasy v1.0 — a regional, accelerated 2D glacial isostatic adjustment
  (GIA) model accounting for the lateral variability of the solid Earth.
  *Geosci. Model Dev.* 17, 5263–5290.
- Turcotte, D. L., & Schubert, G. (2014). *Geodynamics*, 3rd ed. Cambridge UP.

---

## 8. Provenance

This design synthesises three independent agent proposals plus an external
ChatGPT-authored alternative. Key decisions and where they came from:

| Decision | Source |
|---|---|
| Half-derivative operator $\lvert\nabla_S\rvert$ (not $\Delta_S$) for the viscous term | External proposal — corrected a unanimous agent error |
| Kinematic self-gravity, no viscoelastic Love-number convolution | Agent C + external proposal |
| Elastic correction $u_{\text{el}}$ added diagnostically via $h_\ell^L$ | Agent C |
| Kendall et al. (2005) SLE-A as production solver; tabulated SLE-B as fast option | Synthesis (user-requested both options) |
| Wavelength-dependent $\eta_{\text{eff}}(\theta,\phi,\ell)$ from day one | Agent B + external proposal |
| ETD1 integrator default, GMRES fallback | Agent A + external proposal preconditioner formula |
| Reference template: SpeedyWeather.jl | Agent B |
| Drop spherical-shell membrane terms | Agents A and C |
| Degree-2 rotational feedback always on | All three agents |
