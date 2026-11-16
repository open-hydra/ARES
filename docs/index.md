---
title: ARES
hide:
  - navigation
  - toc
---

<div class="ares-hero" markdown>

<div class="ares-hero__text" markdown>

<p class="ares-hero__eyebrow">HYDRA CFD SUITE</p>

<h1 class="ares-hero__title">ARES</h1>

<p class="ares-hero__lead">Open-source solver for compressible real-fluid flows on multi-block structured grids,
from supercritical heat transfer to turbulent boundary layers.</p>

<p class="ares-hero__cmd"><span>$</span> ./install.sh build --compilers=gnu --use-openmp</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Verification](vv/index.md){ .md-button .ares-hero__ghost }

<p class="ares-hero__name"><strong>A</strong> <strong>R</strong>eal-fluid <strong>E</strong>quations
<strong>S</strong>olver</p>

<!-- </div>

<div class="ares-hero__shot" markdown>

<div class="ares-term" markdown>

<div class="ares-term__bar"><span></span><span></span><span></span></div>

```bash
# build
git clone https://github.com/open-hydra/ARES.git
cd ARES && ./install.sh build --compilers=gnu --use-openmp

# run the turbulent flat plate
cd test/Flat_Plate_SST
./ARES.sh solve -p 4
```

</div>

<p class="ares-hero__caption" markdown>From clone to a first run: the [Quick Start](getting-started/quick-start.md)</p> -->

</div>

</div>

<p class="ares-kicker">The solver</p>

## Compressible flow with real-fluid physics { .ares-heading }

<p class="ares-sub" markdown>ARES solves the Euler and Navier–Stokes equations for a real fluid, with properties
read from pressure–enthalpy tables. Read the [Overview](overview.md) for the models in one page.</p>

<div class="ares-tools" markdown>

[:material-weather-windy: <b>Euler and Navier–Stokes</b> <i>Steady or time-accurate flows, inviscid or viscous, on multi-block structured grids.</i>](theory/governing-equations.md){ .ares-tool }
[:material-chart-bell-curve: <b>Shock-capturing numerics</b> <i>MUSCL reconstruction and six Riemann solvers, preconditioned and rotated HLLC included.</i>](theory/riemann-solvers.md){ .ares-tool }
[:material-water-thermometer-outline: <b>Real-fluid thermodynamics</b> <i>Tabulated equation of state and transport properties on a (p, h) grid, through FLINT.</i>](theory/thermo.md){ .ares-tool }
[:material-waves: <b>Turbulence</b> <i>RANS closures: Spalart–Allmaras and its variants, Menter SST, Wilcox k-ω 2006 and SSG/LRR.</i>](theory/turbulence.md){ .ares-tool }
[:material-speedometer-slow: <b>Low-Mach preconditioning</b> <i>Weiss–Smith preconditioning keeps the density-based scheme accurate as the Mach number goes to zero.</i>](theory/preconditioning.md){ .ares-tool }
[:material-chip: <b>Parallel runs and multigrid</b> <i>OpenMP and MPI, alone or together, with geometric multigrid towards steady state.</i>](user/convergence.md){ .ares-tool }

</div>

<p class="ares-kicker"></p>

## Hydra CFD Suite { .ares-heading }

<p class="ares-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. ARES is its real-fluid flow solver.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid .is-here }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>FLINT</b> thermochemistry](https://MarcoGrossi92.github.io/FLINT/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="ares-links" markdown>

[Quick Start](getting-started/quick-start.md)
[Overview](overview.md)
[User Guide](user/index.md)
[Input Parameters](user/registry.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
