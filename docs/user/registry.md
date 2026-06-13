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
| simulation-type | euler | euler , laminar , turbulent | **yes** | Flow regime (euler = inviscid; laminar / turbulent = Navier–Stokes) |

## ARES-Numerics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| time-scheme | euler | euler, RK2, RK3 | no | Time integration solver |
| cfl | 0.5 | > 0 | **yes** | CFL number |
| vnn | 0.3 | > 0 | no | Viscous (von Neumann) stability number |
| cfl-rise-threshold | 0 | >= 0 | no | Iterations over which CFL is ramped up |
| time-accurate | .false. | logical | **yes** | Time-accurate switch (false = march to steady state) |
| integration-variables | prim | prim, prec | no | Update variables: primitive or preconditioned |
| irs | .false. | logical | no | Implicit Residual Smoothing |
| irs-beta | 0.0 | >= 0 | no | IRS smoothing coefficient |
| preconditioning-Uref | -1.0 |  | no | Reference velocity for preconditioning (< 0 ⇒ use sound speed) |
| preconditioning-eps-min | -1.0 |  | no | Minimum cut-off $\varepsilon$ for preconditioning (< 0 ⇒ default 0.10) |
| space-reconstruction | none | MUSCL-SD, MUSCL, first-order, none | **yes** | Space reconstruction method |
| flux-limiter | none | vanalbada, minmod, superbee, vanleer, MC, LIMO3, none | no | Flux limiter for reconstruction |
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

## ARES-RANS

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| turbulence-model | *(none)* | SA, SA-R, SA-RC, SA-QCR2000, SA-rough, SA-rough-QCR2000, SAcomp, SST, Wilcox2006, SSGLRR, none | no | RANS turbulence model |
| Prt | 0.90 | > 0 | no | Turbulent Prandtl number |
| Sct | 0.90 | > 0 | no | Turbulent Schmidt number |
| Sc | 0.7 | > 0 | no | Laminar Schmidt number |
| k-coupling | .false. | logical | no | Couple turbulent kinetic energy into the mean-flow energy equation |
| Prt-correction | .false. | logical | no | Turbulent-Prandtl correction for wall roughness |

!!! note
    `turbulence-model` is only consulted when `simulation-type = turbulent`. For `euler` and `laminar` runs the `[ARES-RANS]` section is ignored.

## ARES-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ini-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Initial-condition (`INPUT/ic.*`) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Solution (`OUTPUT/field.*`) format |
| sol-diter | 1000000000 | > 0 | no | Solution output frequency (iterations) |
| sol-dtime | 1e30 | > 0 | no | Solution output frequency (physical time) |
| sol-overwrite | true | true, false | no | Overwrite the solution file (false ⇒ numbered snapshots) |
| sol-variables | thermo |  | no | Extra solution variable groups to write (`thermo`, `transport`) |
| wall-variables | mech thermal |  | no | Wall variable groups to write (`mech`, `thermal`) |
| res-diter | 1 | > 0 | no | Residual-history write frequency |
| shell-diter | 1 | > 0 | no | Console (shell) update frequency |
| ini-diter | 10000 | > 0 | no | `input.ini` re-read frequency |

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
