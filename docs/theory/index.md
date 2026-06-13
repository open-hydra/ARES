# Theoretical Guide

This section provides the theoretical foundations for the physical models and numerical methods implemented in ARES: the governing equations of compressible real-fluid flow, the finite-volume spatial discretization, Riemann solvers, time integration, low-Mach preconditioning, turbulence closures, and the real-fluid thermodynamic framework.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Compressible Navier–Stokes for a single real fluid: conservative form, the $(p,\mathbf{u},h)$ primitive set, convective and viscous fluxes, source terms

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretization__

    ---

    Cell-centred finite volume, MUSCL / MUSCL-SD reconstruction, six flux limiters, shock detection

    [:octicons-arrow-right-24: Spatial discretization](numerics.md)

-   :material-waves:{ .lg .middle } __Riemann Solvers__

    ---

    Six numerical-flux solvers: Rusanov, PLLF, HLLE, HLLC, and preconditioned / rotated HLLC variants

    [:octicons-arrow-right-24: Riemann solvers](riemann-solvers.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    SSP Runge–Kutta, CFL/VNN stability, implicit residual smoothing, geometric multigrid

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-speedometer-slow:{ .lg .middle } __Low-Mach Preconditioning__

    ---

    Weiss–Smith preconditioning of the pseudo-time derivative for accurate, fast convergence at low Mach number

    [:octicons-arrow-right-24: Preconditioning](preconditioning.md)

-   :material-weather-windy:{ .lg .middle } __Turbulence Modelling__

    ---

    Spalart–Allmaras (R / RC / QCR / rough), Menter SST, Wilcox 2006 $k$–$\omega$, SSG-LRR

    [:octicons-arrow-right-24: Turbulence models](turbulence.md)

-   :material-thermometer:{ .lg .middle } __Real-Fluid Thermodynamics__

    ---

    Tabulated $(p,h)$ equation of state, property interpolation, and the FLINT table back-end

    [:octicons-arrow-right-24: Thermodynamic models](thermo.md)

</div>

---

## Overview

ARES advances the compressible Navier–Stokes equations in conservative form on structured multi-block grids. The numerical pipeline:

| Stage | Method | Page |
|-------|--------|------|
| **Governing system** | Single-fluid Euler / Navier–Stokes with a real-fluid EOS | [Governing Equations](governing-equations.md) |
| **Spatial discretization** | Cell-centred FVM + MUSCL-limited reconstruction | [Spatial Discretization](numerics.md) |
| **Interface fluxes** | Six Riemann solvers (Rusanov, PLLF, HLLE, HLLC and variants) | [Riemann Solvers](riemann-solvers.md) |
| **Time marching** | SSP RK (1–3 stage), IRS, multigrid | [Time Integration](time-integration.md) |
| **Low-Mach** | Weiss–Smith preconditioning | [Preconditioning](preconditioning.md) |
| **Turbulence closure** | SA, SST, Wilcox 2006, SSG-LRR | [Turbulence Modelling](turbulence.md) |
| **Thermodynamics** | Tabulated real-fluid $(p,h)$ properties (FLINT) | [Thermodynamics](thermo.md) |
