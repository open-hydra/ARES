# Initial Conditions

ARES initial conditions are specified in two layers:

1. **The user-facing `[ICB-*]` interface** in `input.ini` — a *name-based* description of the starting state of each block. This is what you write.
2. **The solver-side IC file** — a structured file holding the mesh and the cell-centred flow field for every block, read directly by ARES. ATLAS compiles the `[ICB-*]` blocks (with the mesh and the real-fluid table) into this file.

This page covers both, **for a real-gas (real-fluid) phase** — i.e. a phase declared `type = real-fluid` in its `[GPB-*]` block (internal phase type `RF`). The same builder serves the ideal-gas phase, but the keys below are the real-fluid set; the dispersed-phase (`DP`) and solid-phase (`SP`) initial conditions are out of scope.

!!! info "There is no inline flow field"
    The initial state is **not** written inside the `[ARES-*]` sections of `input.ini`. You describe it in the `[ICB-*]` blocks and let ATLAS materialise the IC file; if ATLAS is not available the file can be prepared by hand as long as it conforms to the [layout below](#the-ic-file-ares-reads).

---

## The `[ICB]` interface

The initial field is described by one `[ICB-Block<b>]` section per block. ATLAS evaluates each block against the mesh and the real-fluid $(p,h)$ table and writes the IC file.

```ini
[ICB-Block1]
phase    = air       ; optional — restricts the block to a named phase (blank = all)
p        = 1.0e5     ; static pressure [Pa]
T        = 555.0     ; static temperature [K]   (give T or h — see below)
u        = 100.0     ; x-velocity [m/s]
rhoK     = 2.25e-3   ; freestream rho*k       (two-equation models)
rhoOmega = 6250.0    ; freestream rho*omega
```

The optional `phase =` key restricts a block to a named phase; `type =` defaults to `homogeneous` (a uniform block). The remaining keys set the **thermodynamic state**, the **velocity**, and the **turbulence** start values.

### Thermodynamic state

Give the **pressure** plus **either** the enthalpy **or** the temperature:

- `p` + `h` → ARES stores $(p,h)$ directly, the natural real-fluid pair; $T$ is recovered from the table.
- `p` + `T` → the enthalpy is obtained by inverting the $(p,h)$ table at the requested $T$.

The bundled HTD case uses `p` + `h`; the flat plates use `p` + `T`.

### Velocity

Give the Cartesian components directly, or a speed magnitude with flow angles:

| Key(s) | Meaning |
|--------|---------|
| `u`, `v`, `w` | Cartesian velocity components [m/s] |
| `vel` + `alpha`, `beta` | Speed magnitude with in-plane / out-of-plane angles [m/s, rad] |

### Turbulence

The transported variables of the active model, under the names of the NASA
[Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/) with a `rho` prefix — they are
**density-weighted**, exactly as ARES advances them and as they are stored in the IC file:

| Key(s) | Model |
|--------|-------|
| `rhoNuTilde` | Spalart–Allmaras (all SA variants) — $\rho\tilde\nu$ (the SA working variable, *not* the eddy viscosity) |
| `rhoK`, `rhoOmega` | SST / Wilcox 2006 — $\rho k$, $\rho\omega$ |
| `rhoRij`, `rhoOmega` | SSG-LRR (RSM) — the Reynolds stresses $\rho R_{ij}$ and $\rho\omega$ |

!!! info "Why density-weighted here and specific at a boundary"
    A real-fluid `[ICB-*]` block takes the density-weighted keys **only**: the $(p,h)$ table carries
    no density, so ATLAS cannot turn a specific value into the weighted one the solver transports.
    The [inlet keys](boundary-conditions.md#inlets) are the opposite — `nuTilde`, `k`, `omega` and
    `Rij` are all *specific*, and ARES weights them with the density of the boundary state. The retired names `mit`
    and `kappa` are refused by ATLAS with a hint; `mut` and `nut` name the **eddy viscosity**, which
    is an output, never an input.

---

## Field from a file

Instead of a constant, any state key can be read from an ASCII table by giving it its `<key>-file` variant together with a `<key>-direction`: it is then interpolated as a **1-D profile** along that coordinate (`x, y, z, r, t`). Without a direction the file is interpolated from a 2-D/3-D field file. As soon as one field comes from a file the block initialization becomes *variable* (cell-by-cell) instead of homogeneous.

```ini
[ICB-Block1]
p           = 48.0e5
h-file      = inlet_profile.dat   ; enthalpy from a table ...
h-direction = y                   ; ... interpolated along y
u           = 100.0
```

This is the IC counterpart of the [spatially-varying boundary data](boundary-conditions.md#spatially-varying-boundary-data) on a face — here it fills a whole block from a measured or precomputed profile.

## Multi-zone initialization

A block can be split into coordinate zones, each pointing at its own section. Adding a `direction` key switches the block to multizone:

```ini
[ICB-Block1]
direction = x
zone1 = plenum
range1 = -0.1 0.0
zone2 = duct
range2 =  0.0 0.5

[plenum]
p = 48.0e5
T = 300.0
u = 0.0

[duct]
p = 46.0e5
T = 310.0
u = 100.0
```

## Initialization from an existing solution

To start from a previous run (a coarser grid, or a different mesh), give `old-solution`; ATLAS interpolates it onto the new mesh:

```ini
[ICB-Block1]
old-solution      = ../coarse/OUTPUT/field.tec
interpolation-law = index          ; outlaw | index | extrude
```

`interpolation-law` selects the mapping — `outlaw` (general distance-weighted), `index` (index-to-index), or `extrude` (2-D → 3-D extrusion, using `theta` and `nz`). `old-block-id` (`0` = auto) picks the source block. This produces the same $(p, h, u, v, w, \dots)$ field as a fresh initialization, ready for a restart-style run.

---

## The IC file ARES reads

ATLAS writes the IC file ARES actually loads. It stores both the **node coordinates** and the **cell-centred flow variables** for every block in a single file — no separate mesh file is used. It is read through the **[ORION](https://github.com/MarcoGrossi92/ORION)** I/O library, which supports several formats:

| Format | Extension | Notes |
|--------|-----------|-------|
| Tecplot ASCII | `.tec` | Human-readable; default in bundled cases |
| Tecplot Binary | `.szplt` | Compact and fast; requires TecIO (`--use-tecio`) |
| VTK structured | `.vts`, `.vtm` | Open standard; good for large datasets |

The format is selected with `ic-format` in `[ARES-IO]` (`tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw`) and must match the file ATLAS actually wrote:

```
INPUT/ic.tec          ← Tecplot ASCII (default)
INPUT/ic.szplt        ← Tecplot binary
INPUT/ic.vtm          ← VTK multi-block, with its .vts block files
```

A vtk initial condition is the `.vtm` file with its block files, which the `.vtm` lists relative to its own folder, as
ORION writes them: a solution of ARES (`OUTPUT/field.vtm` with its `OUTPUT/vtk/` folder) starts a new run when copied
into `INPUT/` as `ic.vtm` with its `vtk/` folder. A `.vtm` whose block files are not there but at
`INPUT/vtk/field<entry>.vts` is read there. The format of the block files is the one they declare (ascii, base64 or raw
appended data), whichever of `vtk ascii` and `vtk raw` `ic-format` names.

The variable list and zone structure are identical across formats.

!!! note "Solution output uses the same formats"
    The solver writes `OUTPUT/field.*` with the twin key `sol-format`, accepting the same four values. See [Running ARES](using.md#io-formats).

### Required variables

ARES uses a **real-fluid primitive set** built around pressure and enthalpy — *not* density. Every IC file stores the node coordinates plus the following cell-centred variables:

| Variable | Symbol | Description | Units |
|----------|:------:|-------------|-------|
| Coordinates | $x, y, z$ | Node coordinates (nodal) | m |
| Pressure | $p$ | Static pressure | Pa |
| Velocity | $u, v, w$ | Cartesian velocity components | m/s |
| Enthalpy | $h$ | Specific static enthalpy | J/kg |

A representative Tecplot ASCII header from a bundled case:

```
 VARIABLES = "x" "y" "z" "p" "u" "v" "w" "h"
 ZONE T = B1-RF, I=129, J=65, K=2, DATAPACKING=BLOCK,
      VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=...
```

Here coordinates `x y z` are nodal and `p u v w h` are cell-centred. From the pressure and enthalpy, ARES recovers the remaining state (temperature, density, speed of sound) by inverting the real-fluid $(p,h)$ table — the *thermo inversion* step reported during loading.

!!! note "Why $(p,h)$ instead of $(\rho, p)$?"
    The real-fluid equation of state is tabulated on a pressure–enthalpy grid. Storing $p$ and $h$ as the primitive thermodynamic pair means the initial state maps directly onto the table without an extra inversion, and it is the natural set for the low-Mach preconditioned update. See [Real-Fluid Thermodynamics](../theory/thermo.md).

### Turbulence variables

When a turbulence model is active, the IC file must also carry the model's transported variables (set in the `[ICB-*]` block by the keys in [Turbulence](#turbulence) above):

| Model | Extra variables | Description |
|-------|-----------------|-------------|
| Spalart–Allmaras (all SA variants) | `rhoNuTilde` = $\rho\tilde\nu$ | Density times the SA working variable |
| SST / Wilcox 2006 ($k$–$\omega$) | `rhoK` = $\rho k$, `rhoOmega` = $\rho\omega$ | Density times the turbulent kinetic energy and the specific dissipation rate |
| SSG-LRR (RSM) | `ru'u' rv'v' rw'w' ru'v' ru'w' rv'w'` = $\rho R_{ij}$, `rhoOmega` | Density times the Reynolds-stress components and $\omega$ |

The band names are those ATLAS ICB writes and ARES writes back in `OUTPUT/field.*`, so an IC file, a
restart and a solution file are interchangeable. A file whose turbulence bands carry any other name
is refused at start-up (the names of earlier files — `mi_t`, `mi_tilde`, `kappa` with its `omega` —
are still accepted), which stops a *specific* field from being read as a density-weighted one.

---

## Restart

When resuming (`newrun = .false.`), ARES reads the initial state from the latest solution file in `OUTPUT/` instead of `INPUT/`: `OUTPUT/field.*`, written at the end of a run, or else the highest-numbered `OUTPUT/field<N>.*` (`sol-overwrite = false`), in the `sol-format` format. The solution files use the same variable layout as the initial condition. The flow field is restored, with the physical time of a time-accurate run and the global iteration count of a steady run: from a vtk solution, `CYCLE`, or `|TIME|` when the file holds no `CYCLE`; from a Tecplot solution, the magnitude of `SOLUTIONTIME`.

A steady run writes its solutions (field, wall, diagnostic) at a negative time, minus the global iteration count: `SOLUTIONTIME` in Tecplot, in zones without `STRANDID` (static zones), `TIME` in vtk. A time-accurate run writes them at its physical time, with the zone of block *b* in strand *b* (`STRANDID` = *b*) in Tecplot. A restart reads a negative time as a steady solution, whose global iteration count is its magnitude. A time-accurate restart from a steady solution starts a new transient at t = 0, with the counters of a new run.

A solution that holds no time, a vtk file without `TIME` or a Tecplot ASCII file without `SOLUTIONTIME`, cannot be restarted from, since a restart goes on from the time of its solution: the run stops at set-up with a message. To start a run from such a solution, use it as the initial condition (`INPUT/ic.*`, `newrun = true`). A solution whose time, or one of the counters it holds (below), is not finite (NaN or infinity) stops the run at set-up as well, with the file and the name of the value.

### Counters stored with every solution

The field solution also stores the counters that decide what the run does and writes next, with the same names in every format: as field data of the `.vts` files of a vtk solution, as data set auxiliary data of a Tecplot solution (`DATASETAUXDATA` records after the zones of a `.tec` file; a script that reads a `.tec` file token by token must skip them):

| Item | Content |
|------|---------|
| time | Physical time; in a steady run minus the global iteration count |
| `CYCLE` | Global iteration count, over all grid levels: first column of `residual-history.dat`, probe rows with `diter` |
| `ARES_ITERATION` | Iteration count of the grid level of the solution: CFL ramp, `sol-diter`, `res-diter` |
| `ARES_OUTPUT_NUMBER` | Number of the last numbered solution (`sol-overwrite = false`) |
| `ARES_GRID_LEVEL` | Grid level being solved (multigrid sequence) |
| `ARES_LAST_OUTPUT_TIME` | Time of the last time-based solution (`sol-dtime`) |
| `ARES_RES_REFERENCE` | Pressure residual of the first iteration of the grid level: the unit of `res-threshold` (steady runs) |

The text formats (`vtk ascii`, `tecplot ascii`) write 17 significant digits, the fewest that give back every 64-bit value, and the binary formats (`vtk raw`, `tecplot binary`) store the bits. A restart from a solution with its counters, written by a run of the same kind (steady or time-accurate), goes on as the uninterrupted run would: the CFL ramp, the output cadence and numbering, the residual history, the probe rows and the multigrid level go on from the counters. The last iteration of a run writes the outputs due there as any other iteration does (the probe rows, and the numbered or time-based solution with its wall and diagnostic files), then the final solution, with the counters as they are after those outputs.

If the run went on after the solution the restart reads (it stopped without writing `field.*`, and the restart reads the last `field<N>`), `residual-history.dat` and the probe files still hold the rows written after that solution: ARES passes every row of these files to the file system when it writes it, so a run stopped before its end leaves in them all the rows it wrote. A restart from a solution with its counters drops those rows first, the rows of later iterations and the time-based probe rows of later times, and then appends its own; a last row cut in the middle (a file cut by hand) is dropped too, and a file that does not end with a newline gets one. With a solution without counters the rows stay and the restart appends its rows after them: remove them before restarting, and leave a newline at the end of the file.

Time-based probe rows (`dtime`) fall on the multiples of `dtime`. A restart that goes on in time writes next the row of the first multiple its time has not reached; a transient started from a steady solution starts again from the first one.

### Limits

- **Solutions without the counters** (written before ARES stored them): the restart restores the field and the time or the global iteration count, and starts the other counters as a new run does. The CFL ramp starts again, `sol-diter` counts from the restart, the numbered solutions start again from `field1` (replacing those of the first part) and a multigrid sequence starts again from the coarsest level. A steady restart from such a solution measures `res-threshold` from the pressure residual of its own first iteration, as a new run, and says so with a warning: to end about where the uninterrupted run would, set `res-threshold` to ε R₁/Rₙ, with ε the threshold of the first part and R₁, Rₙ the pressure residuals of `residual-history.dat` at the first iteration and at the iteration of the solution.
- **Probe variables changed before the restart.** To drop the rows written after the solution, the restart counts the values of each probe row with the `variables` of its own input (a row can take more than one line): if the variables of a probe changed, it may keep or drop the wrong rows.
- **Multigrid.** A run ends on the finest grid level, where `field` and `field<N>` are written: from a solution with its counters, the restart goes on on the finest level, without solving the coarser ones again. The solutions written while a coarser level is solved (`field-level<L>`, on level *L*, and `field-prolongated`, both with the time of their iteration) are not read by a restart.
- **`res-diter` greater than 1, steady runs.** The first iteration of a restart computes the pressure residual, which the uninterrupted run computes only every `res-diter` iterations: until the next multiple the convergence test uses that residual, so the two runs can stop at different iterations.
