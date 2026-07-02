# Boundary Conditions

ARES boundary conditions are specified in two layers:

1. **The user-facing `[BCB-*]` interface** in `input.ini` — a *name-based* description that assigns a boundary type to each block face. This is what you write.
2. **The solver-side `bc.txt`** — a face-by-face table of numeric type codes read directly by ARES. ATLAS compiles the `[BCB-*]` blocks into this file.

This page covers both, **for a real-gas (real-fluid) phase** — i.e. a phase declared `type = real-fluid` in its `[GPB-*]` block (internal phase type `RF`). In ATLAS the real-fluid and ideal-gas phases **share the same boundary builder**, so every boundary type documented here applies unchanged to a real-fluid run. Only the boundary types that the ARES solver actually consumes are documented.

---

## The `[BCB]` interface

Each block has a boundary block that maps its six faces to **named** boundary definitions:

```ini
[BCB-Block1]
phase = air        ; optional — restricts these faces to a named phase (blank = all)
face1 = inflow
face2 = outflow
face3 = symmetry
face4 = wall
face5 = none       ; degenerate in 2-D
face6 = none

[inflow]
type = inlet
g    = 80.6478     ; mass flux [kg/m²s]
T    = 300.0       ; static temperature [K]

[outflow]
type = outlet
p    = 1.0e5       ; back pressure [Pa]

[symmetry]
type = symmetry

[wall]
type = wall
q    = 0.0         ; wall heat flux [W/m²] (0 ⇒ adiabatic)
```

The right-hand side of each `faceN =` line is a **section name you choose** (e.g. `inflow`, `wall`, `symmetry`). That section's `type` key selects the physical boundary condition, and the remaining keys provide its data.

The optional `phase =` key restricts a block's faces to a named phase; for a single real-fluid phase it can be omitted.

### Face numbering

For a block with $n_i \times n_j \times n_k$ cells the six faces are:

| Face | Location | Varying indices |
|:----:|----------|-----------------|
| 1 | $i = 1$ (imin) | $j, k$ |
| 2 | $i = n_i$ (imax) | $j, k$ |
| 3 | $j = 1$ (jmin) | $i, k$ |
| 4 | $j = n_j$ (jmax) | $i, k$ |
| 5 | $k = 1$ (kmin) | $i, j$ |
| 6 | $k = n_k$ (kmax) | $i, j$ |

For 2-D cases ($n_k = 1$) faces 5 and 6 are degenerate (`none` / symmetry); for a 2-D axisymmetric case they carry `axisymmetric`.

---

## Boundary types

The `type` key of a boundary section selects one of the following. They fall into three groups: **geometric / connectivity** boundaries (carry no flow data), **walls**, and **flow** boundaries (inlets / outlets).

| `type` | Group | Meaning | Data keys |
|--------|-------|---------|-----------|
| `symmetry` | geometric | Symmetry plane / inviscid slip wall | — |
| `axisymmetric` | geometric | Axis of a 2-D axisymmetric domain | — |
| `extrapolation` | geometric | Zeroth-order extrapolation (supersonic outflow / far-field) | — |
| `connection` | connectivity | Inter-block abutting interface (usually auto-detected) | — |
| `periodic` | connectivity | Periodic face pair | `blocks`, `faces` |
| `wall` | wall | Viscous no-slip wall (or slip wall if no data) | `q` / `T`, `ks`, `eps` |
| `inlet` | flow | Inflow boundary | see [Inlets](#inlets) |
| `outlet` | flow | Outflow boundary | `p`, `rf` |
| `null` | — | No boundary condition (placeholder) | — |

### Geometric and connectivity boundaries

These carry no data line — just `type =`:

```ini
[symmetry]
type = symmetry      ; symmetry plane / inviscid slip wall

[axis]
type = axisymmetric  ; the axis of a 2-D axisymmetric (2Daxi) domain

[farfield]
type = extrapolation ; zeroth-order extrapolation: supersonic outflow / far-field

[link]
type = connection    ; abutting inter-block interface (normally found automatically)
```

`connection` faces are usually detected automatically by ATLAS from the mesh topology, so you rarely name them by hand. A **periodic** pair does carry data:

```ini
[periodic]
type   = periodic
blocks = 1 1         ; source / destination block
faces  = 3 4         ; source / destination face
```

### Walls

`type = wall` is the viscous no-slip wall. The keys present select the thermal treatment; roughness `ks` (0 = smooth) and emissivity `eps` are optional add-ons:

| Keys present | Wall condition |
|--------------|----------------|
| `q` (`= 0` ⇒ adiabatic) | Prescribed **heat flux** |
| `T` | **Isothermal** wall |
| *(none)* | Inviscid **slip** wall (same as `symmetry`) |

```ini
[adiabatic]
type = wall
q    = 0.0           ; adiabatic (zero heat flux)

[heated]
type = wall
q    = 4.766e6       ; prescribed wall heat flux [W/m²]

[cold]
type = wall
T    = 300.0         ; isothermal wall [K]

[rough]
type = wall
q    = 0.0
ks   = 1.0e-4        ; equivalent sand-grain roughness height [m]
```

### Inlets

`type = inlet` data is given in the **laboratory frame**. ARES supports two inlet flavours, selected by the keys you provide:

| Keys present | Inlet flavour |
|--------------|---------------|
| `g` + `T` | Subsonic, **mass flux + static temperature** |
| `mach` + `T` + `p` | **Supersonic** inflow (Mach + static temperature + static pressure) |

```ini
[inflow]
type  = inlet
g     = 59.29        ; mass flux [kg/m²s]
T     = 555.0        ; static temperature [K]
kappa = 2.25e-3      ; freestream k       (two-equation models)
omega = 6250.0       ; freestream omega
```

Two optional groups apply to every inlet:

- **Flow direction** — `alpha` (in-plane angle), `beta` (out-of-plane angle); with no angles the inflow is face-normal.
- **Turbulence inflow** — `mit` for the SA models, `kappa` + `omega` for the two-equation models, `rhoRij` + `omega` for the RSM. `rf` is a relaxation factor (default 1).

### Outlets

`type = outlet` prescribes a back pressure; `rf` is the relaxation factor:

```ini
[outflow]
type = outlet
p    = 46.11e5       ; back (ambient) pressure [Pa]
```

The inlet and outlet routines fall back to a symmetry wall if the local state would make the boundary ill-posed (e.g. reversed flow). An outlet with `p = 0` reverts to extrapolation.

---

## Splitting a face into patches

A single face can carry different conditions along its length. The face's section is split by coordinate (or index) range:

```ini
[wall]
direction = x          ; coordinate used to split the face (x,y,z,r,t,i,j,k)
patch1    = symmetry   ; section applied on range1
range1    = -2. 0.     ; x ∈ [-2, 0] is a symmetry plane (leading edge)
patch2    = adiabatic  ; section applied on range2
range2    = 0. 2.      ; x ∈ [0, 2] is an adiabatic wall (the plate)

[adiabatic]
type = wall
q    = 0.0
```

This is the standard way to model a flat plate: an inviscid symmetry segment upstream of the leading edge, then a no-slip wall over the plate.

## Spatially-varying boundary data

Each boundary cell gets its own data line in `bc.txt`, so a wall or inlet condition can **vary along the face** instead of being uniform. Give a scalar key its `<key>-file` variant — pointing at an ASCII table — together with a `direction`, and ATLAS interpolates the tabulated values onto the face cells along that coordinate.

The most common use is a **variable wall temperature or heat flux**: pair `T-file` or `q-file` with `direction`, and the wall becomes isothermal with a temperature $T(x)$ — or a heat-flux wall with $q(x)$ — that follows the table along the chosen direction:

```ini
[wall]
type      = wall
direction = x            ; the file is read along x
q-file    = wall_q.dat   ; heat flux q(x)  →  variable heat flux along x
; T-file  = wall_T.dat   ; ... or a variable wall temperature T(x)
```

`wall_q.dat` is a two-column ASCII table `coordinate  value`; ATLAS assigns each wall cell the value interpolated at the cell's coordinate, so `q-file` gives a heat flux that varies along `direction` and `T-file` gives a wall temperature that varies along `direction`. The same `<key>-file` + `direction` mechanism works for the other boundary scalars (`ks-file`, `g-file`, `alpha-file`, …). `direction` accepts `x, y, z, r, t`; with two directions the table is read as a 2-D field and bilinearly interpolated. (For an index-based `direction` of `i, j, k`, add `file-direction` to name the spatial axis the file is tabulated against.)

ARES reads the resulting per-cell `q`/`T` directly (wall codes 301/302), so the variable distribution is honoured by the solver.

---

## The solver-side `bc.txt`

ATLAS writes `INPUT/bc.txt`, the face-by-face table ARES actually reads. It is the ground truth from which the solver builds its boundary data structures.

### File format

Each boundary cell is listed as an identifier line, optionally followed by a data line:

```
<block>  <i>  <j>  <k>  <face>  <type>
<data line, comma-separated, if the type requires it>
```

| Field | Description |
|-------|-------------|
| `block` | 1-based block index |
| `i j k` | Cell indices of the boundary cell |
| `face` | Face index (1–6) |
| `type` | Numeric BC type code |

Types that carry no data (symmetry, axisymmetric, extrapolation) have no second line.

### Numeric type codes

The codes ARES reads (`IO_BC`, `Mod_BC_Fluxes`), with the `[BCB]` type that produces them:

| Code | `[BCB]` `type` | Meaning | Data line |
|:----:|----------------|---------|-----------|
| `0`   | `null`          | Interior / no BC | — |
| `101` | `connection`    | Abutting inter-block interface | — |
| `103` | —               | Multi-solver (HYDRA) coupling: symmetry wall + external flux | — |
| `200` | `axisymmetric`  | Axisymmetric-axis treatment | — |
| `201` | `periodic`      | Periodic face pair (dispatched as a connection) | — |
| `300` | `symmetry` / `wall` (no data) | Symmetry plane / inviscid slip wall | — |
| `301` | `wall` (`q`)    | No-slip wall, prescribed heat flux `qw` (`0` = adiabatic) | `q, ks, eps` |
| `302` | `wall` (`T`)    | No-slip wall, prescribed temperature `Tw` (isothermal) | `T, ks, eps` |
| `400` | `extrapolation` | Zeroth-order extrapolation (supersonic outflow / far-field) | — |
| `404` | `inlet` (`g`,`T`) | Subsonic inlet: mass flux + static temperature | `T, g, α, β, rel` |
| `405` | `inlet` (`mach`,…) | Supersonic inlet: Mach + static `T` + static `p` | `M, T, p, α, β, rel` |
| `406` | `outlet` (`p`)  | Subsonic outlet: back pressure `pAmb` | `p, rel` |

where `α`, `β` are the in-plane / out-of-plane inflow angles (`normal` = face-normal) and `rel` a relaxation factor. The inlet/outlet routines fall back to a symmetry wall if the local state would make the BC ill-posed (e.g. reversed flow).

### Example — mass-flux inlet (type 404)

A mass-flux inlet (`g = 80.65 kg/m²s`, `T = 300 K`) compiles to:

```
       1       1       1       1       1     404
    0.300000E+03,    0.806478E+02,         normal,         normal,    0.100000E+01,
```

The identifier line places the BC on block 1, cell `(1,1,1)`, face 1, with type `404` (inlet: mass flux + static temperature). The data line carries `T`, `g`, the in-plane and out-of-plane flow angles (`normal` = face-normal inflow), and a relaxation factor.

!!! info "Numeric codes and ATLAS"
    In normal use you do not write these codes by hand — you write the `[BCB]` blocks above and let ATLAS translate them. The codes are documented here so you can interpret or hand-edit `bc.txt` when debugging a case.

---

## Multigrid

When geometric multigrid is enabled, ARES needs one boundary table per grid level. ATLAS produces the coarse-level tables alongside the finest one; the per-level iteration counts are set in `[ARES-Multigrid]` (see the [registry](registry.md#ares-multigrid)).
