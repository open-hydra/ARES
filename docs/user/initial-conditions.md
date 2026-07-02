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
phase = air          ; optional — restricts the block to a named phase (blank = all)
p     = 1.0e5        ; static pressure [Pa]
T     = 555.0        ; static temperature [K]   (give T or h — see below)
u     = 100.0        ; x-velocity [m/s]
kappa = 2.25e-3      ; freestream k       (two-equation models)
omega = 6250.0       ; freestream omega
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

The transported variables of the active model:

| Key(s) | Model |
|--------|-------|
| `mit` | Spalart–Allmaras (all SA variants) — modified eddy viscosity |
| `kappa`, `omega` | SST / Wilcox 2006 — $k$, $\omega$ |
| `rhoRij`, `omega` | SSG-LRR (RSM) — Reynolds stresses + $\omega$ |

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

The format is selected with `ini-format` in `[ARES-IO]` (`tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk raw`) and must match the file ATLAS actually wrote:

```
INPUT/ic.tec          ← Tecplot ASCII (default)
INPUT/ic.szplt        ← Tecplot binary
INPUT/ic.vts          ← VTK structured
```

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
| Spalart–Allmaras (all SA variants) | $\tilde\nu$ | Modified eddy viscosity |
| SST / Wilcox 2006 ($k$–$\omega$) | $k$, $\omega$ | Turbulent kinetic energy, specific dissipation rate |
| SSG-LRR (RSM) | $R_{ij}$, $\omega$ | Reynolds-stress components and $\omega$ |

---

## Restart

When resuming (`newrun = .false.`), ARES reads the initial state from the latest solution file in `OUTPUT/` instead of `INPUT/`. The solution files use the same variable layout and format, and the flow field, physical time, and global iteration counter are all restored.
