# Input Parameters

This page is the reference for every `[ARES-*]` parameter. Its tables are the output of the `DocGen` tool (`src/app/docgen.f90`), which builds the registry from the same `Register_*` routines the solver uses at run time, in the order of registration; the notes between the tables are written by hand. To regenerate the tables, build and run `DocGen` (it writes `input-parameters.md`) and copy them over those of this page, keeping the notes.

!!! note
    Only the `[ARES-*]` solver sections are registered here. The pre-processing blocks (`[GRIB-*]`, `[GPB-*]`, `[ICB-*]`, `[BCB-*]`) are interpreted by ATLAS and are documented in the [Input File](input.md) page and the [ATLAS documentation](https://github.com/open-hydra/ATLAS).

## ARES-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false | no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 | no | Residual convergence threshold, relative to the first iteration (steady runs only) |
| time-threshold | 1e30 | > 0 | no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 | no | Maximum number of iterations (re-read at runtime) |

## ARES-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ic-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Initial condition (INPUT/ic.*) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw | no | Solution (OUTPUT/field.*) format |
| decomposition-map | none |  | no | Block-to-rank map in INPUT/ (none: ARES assigns the blocks) |
| sol-diter | 1000000000 | > 0 | no | Solution output frequency (iterations) |
| sol-dtime | 1e30 | > 0 | no | Solution output frequency (physical time) |
| sol-overwrite | true | true, false | no | Overwrite the solution file (false = numbered snapshots) |
| sol-variables | thermo |  | no | Extra solution variable groups to write (thermo, transport) |
| wall-variables | mech thermal |  | no | Wall variable groups to write (mech, thermal) |
| res-diter | 1 | > 0 | no | Residual-history write frequency (iterations) |
| shell-diter | 1 | > 0 | no | Console (shell) update frequency (iterations) |
| ini-diter | 10000 | > 0 | no | input.ini re-read frequency (iterations) |
| ini-error | warn | warn, stop | no | input.ini read again with an error: warn and keep the values held before the reading, or stop (at set-up every error stops) |
| timer-diter | 0 | >= 0 | no | Wall-clock timing report frequency (iterations; 0 = off) |

`decomposition-map`: see [Block-to-rank assignment](using.md#block-to-rank-assignment).

## ARES-Probes

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| probe1 | probe-placeholder |  | no | Name of the section of the probe, and of its file OUTPUT/<name>.txt (16 characters at most; probe1, probe2, ... one per probe) |

## probe-A

Each probe named in `[ARES-Probes]` (`probe1 = <name>`) has a section `[<name>]` with these keys; DocGen registers one probe, named `probe-A`.

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| variables | none |  | no | Probe variables to write |
| dtime | 1e30 | > 0 | no | Probe output frequency (physical time) |
| diter | 1000000000 | > 0 | no | Probe output frequency (iterations) |
| index-position | 0 0 0 0 | >= 0 | no | Probe location by index: block i j k |
| position | 0.0 0.0 0.0 |  | no | Probe location by coordinates: x y z (the cell with the nearest centre) |

## ARES-Numerics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| time-scheme | euler | euler, RK2, RK3 | **yes** | Time integration solver |
| cfl | 0.5 | > 0 | **yes** | CFL number of the convective spectral radii summed over the directions (RK2 <= 1, RK3 <= 1.17) |
| vnn | 0.3 | > 0 | no | Diffusion number of the largest diffusivity summed over the directions (RK2 <= 0.5, RK3 <= 0.62) |
| cfl-rise-threshold | 0 | >= 0 | no | Iterations over which the time step is ramped up from 1/N of its value (0 = off) |
| time-accurate | .false. | logical | **yes** | Time-accurate switch (false = march to steady state) |
| integration-variables | prim | prim, prec | no | Update variables: primitive (prim) or preconditioned (prec, with riemann-solver = HLLC Prec) |
| irs-beta | 0.0 | >= 0 | no | Implicit Residual Smoothing coefficient (> 0 switches IRS on; read at setup only; steady runs only; 0.5 damps the shortest wave the most) |
| irs-variables | all | flow, all | no | Residuals smoothed by IRS: all (also the RANS equations; for Spalart-Allmaras: the k-omega models can diverge with all) or flow (mass, momentum, energy) |
| preconditioning-Uref | -1.0 |  | no | Reference velocity for preconditioning (< 0 = local sound speed; it sets U_r only: alpha takes the speed of sound) |
| preconditioning-eps-min | -1.0 |  | no | Minimum cut-off epsilon of the preconditioning parameter (< 0 = default 0.10) |
| space-reconstruction | none | MUSCL, first-order, none | **yes** | Space reconstruction method |
| flux-limiter | none | vanalbada, minmod, superbee, vanleer, mc, MC, LIMO3, none | no | Flux limiter for MUSCL reconstruction (van Leer if omitted or none; ignored with first-order; mc and MC are the same limiter) |
| diffusive-metric | stencil | stencil, mean | no | Face metric of the viscous fluxes: stencil takes the face Jacobian from the same differences of the cell centres as the face gradient, exact for a linear field on any mesh (at the symmetry faces and the k faces of type 200 for the linear fields that they allow) except at the chimera faces (102), which take the inverse of the mean of the inverse metrics of the two cells with every value, and at the symmetry faces with one cell across them, which take the metric of the cell; mean takes the mean of the two metrics (the faces between blocks and of type 200 the inverse of the mean of the inverse metrics), overstates the gradient by (1+r)^2/(4r) at a cell-size ratio r and does not converge next to a symmetry plane (300) that bounds a sector of a polar mesh |
| riemann-solver | HLLC | Rusanov, PLLF, HLLE, HLLC, HLLC Prec, HLLC Rotated | no | Riemann (numerical flux) solver (HLLC Prec with prec only, and prec with HLLC Prec only; HLLC Rotated on 2-D meshes only) |

!!! tip "Preconditioning"
    Setting `integration-variables = prec` activates [low-Mach preconditioning](../theory/preconditioning.md). The `preconditioning-Uref` and `preconditioning-eps-min` keys then control the reference velocity and the lower cut-off of the preconditioning parameter. It takes the preconditioned Riemann solver `HLLC Prec` (the set-up stops on any other).

`diffusive-metric`: see [Spatial Discretization](../theory/numerics.md#gradient-computation-for-diffusive-fluxes).

## ARES-Multigrid

Per-level iteration counts. The number of levels is set by `levels` (an integer, default 1, at least 1), read at set-up only and outside the registry; one `levelN-iter` key is registered per level (DocGen registers two levels). The finest level runs `level1-iter` iterations, or `iter-threshold` when `level1-iter` is 0; a coarse level whose `levelN-iter` is 0 runs one iteration, with a warning ([Time Integration](../theory/time-integration.md#multigrid-acceleration)).

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| level1-iter | 0 | >= 0 | no | Iterations for multigrid level 1 (re-read at runtime; 0 = iter-threshold on level 1, one iteration on a coarse level) |
| level2-iter | 0 | >= 0 | no | Iterations for multigrid level 2 (re-read at runtime; 0 = iter-threshold on level 1, one iteration on a coarse level) |
| prolongated-diter | 0 | >= 0 | no | Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write; re-read at runtime) |

## ARES-Physics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| equations | *(none)* | euler , navier-stokes | **yes** | Gas dynamics equations |
| turbulence | none | SA, SA-R, SA-RC, SA-QCR2000, SA-rough, SA-rough-QCR2000, SAcomp, SST, SST-RC, Wilcox2006, SSGLRR, SSGLRR-SD, SSGLRR-w2019, SSGLRR-w2019-SD, none | no | RANS turbulence model |
| schmidt | 0.0 | >= 0 | no | Laminar Schmidt number (unused: ARES solves a single species) |

!!! note
    `navier-stokes` with `turbulence = none` is a laminar run. `euler` takes no turbulence model: naming one stops the run at set-up.

## ARES-Turbulence

These keys are read for every run and used only when `turbulence` names a model.

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| Prt | 0.90 | > 0 | no | Turbulent Prandtl number |
| Sct | 0.90 | > 0 | no | Turbulent Schmidt number |
| point-implicit | .false. | logical | no | Point-implicit factor 1/(1 + dt d) on the increment of the turbulence equations, steady runs only (time-accurate runs: explicit sources, stiffness in the time step) |
| wall-distance | exact | exact, face-centre | no | Wall distance of the turbulence models (yn, and the roughness of the nearest wall face): exact (minimum distance from the cell centre to the wall surface, NASA TMR) or face-centre (distance to the nearest centre of a wall face, as ARES did before the exact distance; kept to compare with the runs of that version: larger than the exact distance where the cell centre is not above a face centre, on skewed cells near the wall) |
| omega-wall-bc | practical | practical, asymptotic | no | Omega wall condition for k-omega models and SSGLRR: Menter practical (60*nu/(beta1*y^2), beta1 of SST, for every model) or asymptotic (6*nu/(beta*y^2), the value of the y->0 limit set at the wall face, with the beta of the model, beta1 for SST and SSGLRR and beta0 for Wilcox2006; kept for comparison: practical is the wall condition of the references) |
| sst-production | incompressible | incompressible, compressible | no | SST production (NASA TMR naming): incompressible (P = mu_t*S^2, SST-2003m) or compressible (P = tau_ij*du_i/dx_j, SST-2003e, as OpenFOAM kOmegaSST: adds -2/3*mu_t*divU^2 - 2/3*rho*k*divU, and -2/3*gamma*rho*omega*divU to the omega equation) |
| k-energy-coupling | .false. | logical | no | Two-equation models (SST, Wilcox2006): k in the total energy, the full compressible RANS form of NASA TMR: -2/3*rho*k*delta_ij in the stress tensor of the momentum and energy equations, and the k-equation source (production - dissipation) taken from the mean-flow energy. .false. = TMR m/e variants |
| Prt-correction | .false. | logical | no | Turbulent-Prandtl correction for wall roughness (SA with -rough only) |
| sst-asymptotic | .false. | logical | no | Impose the asymptotic solution omega = 6*nu/(beta_1*d^2) on the near-wall cells (SST only) |
| sst-asymptotic-cells | 1 | >= 1 | no | Number of cells off viscous walls where the asymptotic omega is imposed |

`point-implicit`: see [Numerical treatment of the source terms](../theory/turbulence.md#numerical-treatment-of-the-source-terms).

!!! note "Defaults"
    `Prt` defaults to 0.90 and `point-implicit` to `.false.`.
