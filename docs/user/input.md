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
iter-threshold  = 700000

[ARES-Physics]
equations  = navier-stokes
turbulence = SST

[ARES-Numerics]
cfl                   = 0.6
time-scheme           = RK3
integration-variables = prim
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC
```

!!! warning "Defaults are silent"
    If a parameter is omitted, the default is used without warning. Check the [Parameter Registry](registry.md) for exact names, defaults, and allowed values — names are case-sensitive.

---

## Solver sections (`[ARES-*]`)

These are read by ARES and validated against the registry.

| Section | Description | Reference |
|---------|-------------|-----------|
| `[ARES-Parameters]` | Convergence thresholds, restart | [→](registry.md#ares-parameters) |
| `[ARES-Physics]` | Equations and turbulence model | [→](registry.md#ares-physics) |
| `[ARES-Turbulence]` | Options of the turbulence model: turbulent Prandtl/Schmidt numbers, omega wall condition, ... | [→](registry.md#ares-turbulence) |
| `[ARES-Numerics]` | Time/space discretization, Riemann solver, CFL, IRS, preconditioning | [→](registry.md#ares-numerics) |
| `[ARES-IO]` | Output formats, frequencies, and variables | [→](registry.md#ares-io) |
| `[ARES-Probes]` | Point probes for time-history recording | [→](registry.md#ares-probes) |
| `[ARES-Multigrid]` | Per-level multigrid iteration counts | [→](registry.md#ares-multigrid) |

!!! note "Same layout as MOSE"
    The flow regime is set in `[ARES-Physics]` with the keys of `[MOSE-Physics]`: `equations` = `euler` · `navier-stokes`, and `turbulence` = the RANS model (default `none`). `navier-stokes` with `turbulence = none` is a laminar run; `euler` takes no turbulence model. The options of the model are in `[ARES-Turbulence]`, as in `[MOSE-Turbulence]`.

!!! warning "Input files of earlier ARES versions"
    `simulation-type` in `[ARES-Parameters]`, the `[ARES-RANS]` section and `ini-format` in `[ARES-IO]` are no longer read. ARES stops at set-up if it finds them, and says where each key went:

    | Earlier key | Now |
    |-------------|-----|
    | `[ARES-Parameters]` `simulation-type = euler` | `[ARES-Physics]` `equations = euler` |
    | `[ARES-Parameters]` `simulation-type = laminar` | `[ARES-Physics]` `equations = navier-stokes` |
    | `[ARES-Parameters]` `simulation-type = turbulent` + `[ARES-RANS]` `turbulence-model = <model>` | `[ARES-Physics]` `equations = navier-stokes`, `turbulence = <model>` |
    | `[ARES-RANS]` `Sc` | `[ARES-Physics]` `schmidt` |
    | `[ARES-RANS]` `Prt`, `Sct`, `point-implicit`, `Prt-correction`, `sst-asymptotic`, `sst-asymptotic-cells` | `[ARES-Turbulence]`, same names |
    | `[ARES-IO]` `ini-format` | `[ARES-IO]` `ic-format` |

    `time-scheme` in `[ARES-Numerics]` is now required: an input without it stops at set-up.

---

## Pre-processing sections (ATLAS)

These define the case and are consumed by **ATLAS** to generate the `INPUT/` data. They are documented here because they live in the same file, but their full semantics belong to the [ATLAS documentation](https://github.com/open-hydra/ATLAS).

| Section | Role |
|---------|------|
| `[GRIB-*]` | **GRI**d **B**lock — mesh generation (e.g. `method = gmsh`, a `.geo`/`.msh` input, or analytic surfaces) |
| `[GPB-*]` | **G**as **P**roperty **B**lock — the real-fluid table: `type = real-fluid`, `fluid`, pressure/temperature range (`pmin/pmax/Tmin/Tmax`), and resolution (`NP × NH`) (see [Gas Properties](gas-properties.md)) |
| `[ICB-*]` | **I**nitial **C**ondition **B**lock — the initial flow state per block (`p`, `T`, `u`, and the density-weighted turbulence variables `rhoK`/`rhoOmega` if applicable) (see [Initial Conditions](initial-conditions.md)) |
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
