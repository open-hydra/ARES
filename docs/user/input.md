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

!!! note "Comments"
    A line that starts with `!`, `;` or `#` is a comment, also when it holds `=` or a bracket. After a value only `;`
    starts a comment (`cfl = 0.9 ; CFL number`): a `#` or a `!` after a value is part of the value. A number followed by
    one does not read and stops the run; a name followed by one stops the run when its parameter takes one of a list of
    names.

!!! warning "Defaults are silent"
    If a parameter is omitted, the default is used without warning. Check the [Parameter Registry](registry.md) for exact names, defaults, and allowed values — names are case-sensitive.

!!! warning "Keys that are not read"
    ARES warns at the start for every line `key = value` of `input.ini` whose key the reader does not give in its section, with its line: such a key is an error in the input, and the run goes on with the default of its parameter in place of the value written. A parameter of the registry written with no value is named here, and stops the run (below). Comment lines are skipped, also when they hold `=`.

!!! warning "An input.ini that cannot be read stops the run at the start"
    If `input.ini` cannot be opened (missing, or without read permission), holds no line that is not blank, or the INI reader returns an error while it loads the file (a file with no section, a line it cannot parse, which the reader names on its error output), the run stops at the start with the reason and exit code 1; with MPI every rank stops. For a reading during the run, see the note «Errors of a reading during the run» below.

!!! warning "Values outside their rules stop the run at the start"
    A value outside the rule of its parameter (column Allowed of the [Parameter Registry](registry.md): a bound such as `> 0`, an interval such as `in (0, 1)`, where `(` and `)` leave the end out and `[` and `]` take it in, each held by every value of a list of numbers; a list of names; `logical`) stops the run with a message at the start. ARES reads `input.ini` again during the run, every `ini-diter` iterations (10000 by default), to take the values changed meanwhile; what an error of such a reading does is in the note «Errors of a reading during the run» below. In an MPI run rank 0 reads `input.ini`, at the start and at every reading, and passes its text to every rank, so that every rank takes the same values. The parameters that `input.ini` leaves out keep their values.

!!! warning "Values that do not read stop the run at the start"
    A value that does not read as the type of its parameter in the [Parameter Registry](registry.md) stops the run at the start with a message that gives the file, the section, the key and the value as written (`levels` of `[ARES-Multigrid]` is read at the start only): a real must be a finite number (not a word, `NaN`, `0,5`, or `1e999`, which gfortran reads as an infinity), an integer has no decimal point and no exponent (not `2.0` or `1e3`), a logical starts with `T` or `F` after an optional period (not `yes` or `1`; a word that starts with one of them reads as true or false), a list holds no more values than its parameter, and a key has a value.

!!! warning "Errors of a reading during the run"
    When ARES reads `input.ini` again during the run, an error of the reading (a file that cannot be read, a value that does not read as the type of its parameter, is out of its rule or is not finite) does what `ini-error` of `[ARES-IO]` says. With `warn`, the default, the reading is not taken: ARES writes a warning with the iteration and the error, which names the key and its value as written, and the run goes on with every value it had before the reading, those of the keys of the same reading that were right too. With `stop` the run stops with the message and exit code 1, as at the start; with MPI every rank stops. The value of `ini-error` in force when the reading starts decides. At the start every error stops the run, whatever `ini-error` says.

!!! warning "What a reading during the run changes"
    When ARES reads `input.ini` again during the run, only the parameters that the solver reads at every step take the new values, such as `cfl`, `vnn`, the intervals of the outputs, `iter-threshold` and `level<n>-iter`. The choices that the set-up makes once from the input are not made again, nor are its checks between keys: a new value of `equations`, `turbulence`, `space-reconstruction`, `flux-limiter`, `riemann-solver` or `time-scheme` is checked against its rule and leaves the run as it was. Residual smoothing cannot be turned on during a run started without it (`irs = true` or `irs-beta > 0`): its work arrays exist only when it was on at the start, so the reading ignores it, with a warning.

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

!!! note "Turbulence sources in time-accurate runs"
    In time-accurate RANS runs (`time-accurate = .true.`, or a time-accurate HYDRA coupling) the turbulence sources
    are explicit and the local time step includes their stiffness, with no parameter: near walls, where the
    destruction of $k$, $\omega$ or $\tilde\nu$ is fast, the step can be smaller than the convective and diffusive
    limits ([theory](../theory/turbulence.md#time-accurate-runs-explicit-sources)). `point-implicit` serves steady
    runs only: a `point-implicit = .true.` written in `input.ini` is not applied there, and a warning at start-up says
    so. With SST or Wilcox 2006, a time-accurate start far from the wall value of $\omega$ (100 times or more below it,
    in wall cells where diffusion sets the step) needs the CFL ramp `cfl-rise-threshold`: without it the first step
    can make $\omega$ negative. At the end of a RANS run in which a stage made $\tilde\nu$, $k$ or $\omega$ negative,
    one line gives how many cell updates had such a stage.

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

## Time step: `cfl` and `vnn`

Each cell takes the largest explicit step that keeps it stable, from one rule
([theory](../theory/time-integration.md#stability-conditions)):

- `cfl` is the Courant number of the **sum over the directions** of the convective
  spectral radii, $|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$, measured on the
  width of the cell across its faces;
- `vnn` is the diffusion number of the **largest diffusivity** of the equations, also summed
  over the directions: in air the energy one, $\gamma\nu/\text{Pr} \approx 1.94\,\nu$;
  with Spalart–Allmaras also $(\nu + \tilde\nu)/\sigma$;
- the convective and the diffusive rates **add up**: a cell where both matter gets a step
  smaller than either limit alone;
- in time-accurate RANS runs the stiffness of the turbulence sources adds up as well, with no
  parameter (note *Turbulence sources in time-accurate runs* above).

Stability limits of the linearised schemes (MUSCL with upwind fluxes, central diffusion);
the additive rule stays stable for any `cfl` and `vnn` up to these values:

| `time-scheme` | `cfl` ≤ | `vnn` ≤ |
|---|---|---|
| `RK2` (Heun) | 1 | 0.5 |
| `RK3` | 1.17 | 0.62 |
| `euler` | 1 with `first-order` reconstruction (unstable with `MUSCL`) | 0.5 |

These are linear limits. With the limiter, positivity (a maximum principle; the scheme is not TVD in more than one
dimension) is guaranteed only up to `cfl = 0.5` and `vnn = 0.5`, with Heun and with RK3
([theory](../theory/time-integration.md#local-time-step)), and the `vnn` limits are sharp: a run
a little beyond them blows up. Suggested values:

- `RK3`: `cfl = 0.8` (above the guaranteed 0.5) and `vnn = 0.5`;
- `RK2` (Heun): `cfl = 0.5` with moving shocks, up to `cfl = 0.7` where the flow is smooth, and
  `vnn = 0.4`; never `cfl` ≥ 0.95 or `vnn` ≥ 0.5.

The defaults are `cfl = 0.5` and `vnn = 0.3`; in a viscous time-accurate run raise `vnn` within the limits of the
table before `cfl`.

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
