# Input Parameters

This page is the reference for every `[ARES-*]` parameter. It mirrors the output of the `DocGen` tool (`src/app/docgen.f90`), which builds the registry from the same `Register_*` routines the solver uses at run time. To regenerate it, build and run `DocGen` (it writes `input-parameters.md`).

!!! note
    Only the `[ARES-*]` solver sections are registered here. The pre-processing blocks (`[GRIB-*]`, `[GPB-*]`, `[ICB-*]`, `[BCB-*]`) are interpreted by ATLAS and are documented in the [Input File](input.md) page and the [ATLAS documentation](https://github.com/open-hydra/ATLAS).

## ARES-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false | no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 | no | Residual convergence threshold |
| time-threshold | 1e30 | > 0 | no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 | no | Maximum number of iterations |

## ARES-Numerics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| time-scheme | euler | euler, RK2, RK3 | **yes** | Time integration solver |
| cfl | 0.5 | > 0 | **yes** | CFL number |
| vnn | 0.3 | > 0 | no | Viscous (von Neumann) stability number |
| cfl-rise-threshold | 0 | >= 0 | no | Iterations over which CFL is ramped up |
| time-accurate | .false. | logical | **yes** | Time-accurate switch (false = march to steady state) |
| integration-variables | prim | prim, prec | no | Update variables: primitive or preconditioned |
| irs | .false. | logical | no | Implicit Residual Smoothing |
| irs-beta | 0.0 | >= 0 | no | IRS smoothing coefficient |
| preconditioning-Uref | -1.0 |  | no | Reference velocity for preconditioning (< 0 ⇒ use sound speed) |
| preconditioning-eps-min | -1.0 |  | no | Minimum cut-off $\varepsilon$ for preconditioning (< 0 ⇒ default 0.10) |
| space-reconstruction | none | MUSCL, first-order, none | **yes** | Space reconstruction method |
| flux-limiter | none | vanalbada, minmod, superbee, vanleer, mc, MC, LIMO3, none | no | Flux limiter for reconstruction (`mc` and `MC` are the same limiter) |
| riemann-solver | HLLC | Rusanov, PLLF, HLLE, HLLC, HLLC Prec, HLLC Rotated | no | Riemann (numerical flux) solver |
| riemann-options-Minf | 0.0 | >= 0 | no | Reference Mach number for the Riemann solver |

!!! tip "Preconditioning"
    Setting `integration-variables = prec` activates [low-Mach preconditioning](../theory/preconditioning.md). The `preconditioning-Uref` and `preconditioning-eps-min` keys then control the reference velocity and the lower cut-off of the preconditioning parameter. Pair it with a preconditioned Riemann solver (`PLLF` or `HLLC Prec`).

## ARES-Multigrid

Per-level iteration counts. The number of levels is set by the `levels` key; one `levelN-iter` entry is registered per level.

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| levels | 1 | >= 1 | no | Number of multigrid levels |
| level1-iter | 0 | >= 0 | no | Iterations on multigrid level 1 (finest) |
| level2-iter | 0 | >= 0 | no | Iterations on multigrid level 2 |
| … | 0 | >= 0 | no | (one entry per level) |

## ARES-Physics

The keys of `[MOSE-Physics]` that apply to a single real fluid.

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| equations | *(none)* | euler , navier-stokes | **yes** | Gas dynamics equations |
| turbulence | none | SA, SA-R, SA-RC, SA-QCR2000, SA-rough, SA-rough-QCR2000, SAcomp, SST, SST-RC, Wilcox2006, SSGLRR, SSGLRR-SD, SSGLRR-w2019, SSGLRR-w2019-SD, none | no | RANS turbulence model |
| schmidt | 0.0 | >= 0 | no | Laminar Schmidt number (no effect: ARES solves a single species) |

!!! note
    `navier-stokes` with `turbulence = none` is a laminar run. `euler` takes no turbulence model: naming one stops the run at set-up.

## ARES-Turbulence

These keys are read for every run and used only when `turbulence` names a model.

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| Prt | 0.90 | > 0 | no | Turbulent Prandtl number |
| Sct | 0.90 | > 0 | no | Turbulent Schmidt number |
| point-implicit | .false. | logical | no | Point-implicit (Patankar) treatment of the turbulence destruction source terms (steady runs only) |
| omega-wall-bc | practical | practical, asymptotic | no | Omega wall condition of the k-omega models (SST, Wilcox2006) and of SSGLRR: Menter practical $60\nu/(\beta_1 y^2)$, or asymptotic $6\nu/(\beta_1 y^2)$, the exact $y \to 0$ limit, for wall-resolved meshes |
| sst-production | incompressible | incompressible, compressible | no | SST production (NASA TMR naming): incompressible ($P = \mu_t S^2$, SST-2003m) or compressible ($P = \tau_{ij}\,\partial u_i/\partial x_j$, SST-2003e, as OpenFOAM kOmegaSST: adds $-\tfrac23\mu_t(\nabla\cdot\mathbf u)^2 - \tfrac23\rho k\,\nabla\cdot\mathbf u$, and $-\tfrac23\gamma\rho\omega\,\nabla\cdot\mathbf u$ to the $\omega$ equation) |
| k-energy-coupling | .false. | logical | no | Two-equation models (SST, Wilcox2006): $k$ in the total energy, the full compressible RANS form of NASA TMR: $-\tfrac23\rho k\,\delta_{ij}$ in the stress tensor of the momentum and energy equations, and the $k$-equation source (production − dissipation) taken from the mean-flow energy. `.false.` = TMR m/e variants |
| Prt-correction | .false. | logical | no | Turbulent-Prandtl correction for wall roughness (SA with `-rough`) |
| sst-asymptotic | .false. | logical | no | Impose the asymptotic solution $\omega = 6\nu/(\beta_1 d^2)$ on the near-wall cells (SST only) |
| sst-asymptotic-cells | 1 | >= 1 | no | Number of cells off viscous walls where the asymptotic $\omega$ is imposed |

!!! note "Defaults that differ from MOSE"
    `Prt` defaults to 0.90 in ARES and to 0.85 in MOSE; `point-implicit` to `.false.` in ARES and to `.true.` in MOSE.

## ARES-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ic-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Initial-condition (`INPUT/ic.*`) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Solution (`OUTPUT/field.*`) format |
| sol-diter | 1000000000 | > 0 | no | Solution output frequency (iterations) |
| sol-dtime | 1e30 | > 0 | no | Solution output frequency (physical time) |
| sol-overwrite | true | true, false | no | Overwrite the solution file (false ⇒ numbered snapshots) |
| sol-variables | thermo |  | no | Extra solution variable groups to write (`thermo`, `transport`) |
| wall-variables | mech thermal |  | no | Wall variable groups to write (`mech`, `thermal`) |
| res-diter | 1 | > 0 | no | Residual-history write frequency |
| shell-diter | 1 | > 0 | no | Console (shell) update frequency |
| ini-diter | 10000 | > 0 | no | `input.ini` re-read frequency |
| timer-diter | 0 | >= 0 | no | Wall-clock timing report frequency (0 ⇒ off) |

## ARES-Probes

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| probe1 | probe-placeholder |  | no | Name of the section defining probe 1 (add `probe2`, `probe3`, … for more) |

### Per-probe section

Each probe name listed above refers to a section with these keys:

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| variables | none |  | no | Variables to record at the probe |
| dtime | 1e30 | > 0 | no | Probe output frequency (physical time) |
| diter | 1000000000 | > 0 | no | Probe output frequency (iterations) |
| index-position | 0 0 0 0 | >= 0 | no | Probe location by index (`i j k block`) |
| position | 0.0 0.0 0.0 |  | no | Probe location by coordinates (`x y z`) |
