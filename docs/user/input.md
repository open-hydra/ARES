# Input File

ARES is configured through a single **INI-format** file, `input.ini`, in the case root directory. It is shared with the ATLAS pre-processor: ATLAS reads the case-definition blocks and generates the data in `INPUT/`, while ARES reads the solver blocks.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Parameters not listed take their **default** value. Example:

```ini
[ARES-Parameters]
simulation-type = turbulent
iter-threshold  = 700000

[ARES-Numerics]
cfl                   = 0.6
time-scheme           = RK3
integration-variables = prim
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[ARES-RANS]
turbulence-model = SST
```

!!! warning "Defaults are silent"
    If a parameter is omitted, the default is used without warning. Check the [Parameter Registry](registry.md) for exact names, defaults, and allowed values — names are case-sensitive.

---

## Solver sections (`[ARES-*]`)

These are read by ARES and validated against the registry.

| Section | Description | Reference |
|---------|-------------|-----------|
| `[ARES-Parameters]` | Simulation type, convergence thresholds, restart | [→](registry.md#ares-parameters) |
| `[ARES-Numerics]` | Time/space discretization, Riemann solver, CFL, IRS, preconditioning | [→](registry.md#ares-numerics) |
| `[ARES-IO]` | Output formats, frequencies, and variables | [→](registry.md#ares-io) |
| `[ARES-RANS]` | Turbulence model and turbulent transport numbers | [→](registry.md#ares-rans) |
| `[ARES-Probes]` | Point probes for time-history recording | [→](registry.md#ares-probes) |
| `[ARES-Multigrid]` | Per-level multigrid iteration counts | [→](registry.md#ares-multigrid) |

!!! note "`simulation-type`, not `equations`"
    The flow regime is selected with a single key, `simulation-type` = `euler` · `laminar` · `turbulent`, in `[ARES-Parameters]`. There is no separate "equations" key: `euler` solves the inviscid system, while `laminar` and `turbulent` solve Navier–Stokes (the latter activating the model named in `[ARES-RANS]`).

---

## Pre-processing sections (ATLAS)

These define the case and are consumed by **ATLAS** to generate the `INPUT/` data. They are documented here because they live in the same file, but their full semantics belong to the [ATLAS documentation](https://github.com/open-hydra/ATLAS).

| Section | Role |
|---------|------|
| `[GRIB-*]` | **GRI**d **B**lock — mesh generation (e.g. `method = gmsh`, a `.geo`/`.msh` input, or analytic surfaces) |
| `[GPB-*]` | **G**as **P**roperty **B**lock — the real-fluid table: `type = real-fluid`, `fluid`, pressure/temperature range (`pmin/pmax/Tmin/Tmax`), and resolution (`NP × NH`) |
| `[ICB-*]` | **I**nitial **C**ondition **B**lock — the initial flow state per block (`p`, `T`, `u`, and turbulence variables `kappa`/`omega` if applicable) |
| `[BCB-*]` | **B**oundary **C**ondition **B**lock — assigns a named boundary type to each of the six block faces (see [Boundary Conditions](boundary-conditions.md)) |

A minimal real-fluid block:

```ini
[GPB-Phase1]
type  = real-fluid
fluid = air
pmin  = 0.80e5
pmax  = 2.0e5
Tmin  = 280.0
Tmax  = 340.0
NP    = 200
NH    = 200
```

This requests a real-fluid air table on a $200\times200$ pressure–enthalpy grid spanning the given $p$ and $T$ ranges. See [Real-Fluid Thermodynamics](../theory/thermo.md) for what the table contains and how it is used.

---

## Parameter Registry

The complete list of all `[ARES-*]` parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code by the `DocGen` tool (`src/app/docgen.f90`) — do not edit it by hand.
