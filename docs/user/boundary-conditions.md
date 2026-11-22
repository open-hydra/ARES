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
The null type has no flux through the k faces. In a planar case (not an axisymmetric wedge) with symmetry on all the k faces ARES sets the z velocity to zero and does not integrate the z momentum. The two mirror planes one cell apart would otherwise give $w$ the mode $dw/dt = -2\,a\,w/\Delta z$ ($a$ the speed of sound, $\Delta z$ the thickness of the layer), which the explicit schemes cannot follow at the time step of the in-plane cells: round-off in $w$ would grow until the run fails. ARES does so only when every k face is a plane normal to z (the x and y components of its unit normal at most $10^{-10}$). Where the k faces are not normal to z, as in a layer between parallel planes inclined to z, symmetry on them asks the gas to move along them, across the layer, and setting $w$ to zero would change the flow (the mass flow, the pressure): there ARES integrates the z momentum, with a warning, and the mode above limits the time step. Heun follows it for time steps up to $\Delta z/a$ and RK3 up to $1.26\,\Delta z/a$: with the time step of the in-plane cells, about $\text{cfl}\,\Delta x/a$, layers thinner than about $\text{cfl}\,\Delta x$ (Heun) or $0.8\,\text{cfl}\,\Delta x$ (RK3) make the run fail.
In a 1-D case (one cell in j and one in k) the same holds in both directions: with symmetry on all the k faces ARES sets $w$ to zero and does not integrate it, and with symmetry on all the j faces it does the same for $v$, whose mode is $dv/dt = -2\,a\,v/\Delta y$, with the same condition: every j face a plane normal to y for $v$, every k face a plane normal to z for $w$. A channel whose section varies along x, one cell across with its faces along the walls, has faces that are not, and the velocity across is integrated, with a warning; cells thinner across than about $\text{cfl}\,\Delta x$ (Heun) or $0.8\,\text{cfl}\,\Delta x$ (RK3) then make the run fail.
With symmetry on only one of the two faces and extrapolation (`400`) on the other, the velocity across keeps half that mode, $dv/dt = -a\,v/\Delta y$, and is integrated, since the extrapolation face lets the gas through. Heun follows that mode for time steps up to $2\,\Delta y/a$ and RK3 up to $2.51\,\Delta y/a$; as the time step of a 1-D case is about $\text{cfl}\,\Delta x/a$, cells thinner across than about $\text{cfl}\,\Delta x/2$ (Heun) or $0.4\,\text{cfl}\,\Delta x$ (RK3) make the run fail. Give both faces the same type.

---

## Boundary types

The `type` key of a boundary section selects one of the following. They fall into three groups: **geometric / connectivity** boundaries (carry no flow data), **walls**, and **flow** boundaries (inlets / outlets).

| `type` | Group | Meaning | Data keys |
|--------|-------|---------|-----------|
| `symmetry` | geometric | Symmetry plane / inviscid slip wall | — |
| `axisymmetric` | geometric | k faces (5 and 6) of an axisymmetric wedge or sector (the axis itself takes `symmetry`) | — |
| `extrapolation` | geometric | Zeroth-order extrapolation (supersonic outflow / far-field) | — |
| `connection` | connectivity | Inter-block abutting interface (usually auto-detected) | — |
| `periodic` | connectivity | Periodic face pair | `blocks`, `faces` |
| `wall` | wall | Viscous no-slip wall (or slip wall if no data) | `q` / `T`, `ks`, `eps` |
| `inlet` | flow | Inflow boundary | see [Inlets](#inlets) |
| `outlet` | flow | Outflow boundary | `p`, `rf` |
| `manifold` | flow | Faces joined by a lumped plenum (a turn-around, a header) | see [Manifold](#manifold) |
| `null` | — | No boundary condition (placeholder) | — |

### Geometric and connectivity boundaries

These carry no data line — just `type =`:

```ini
[symmetry]
type = symmetry      ; symmetry plane / inviscid slip wall

[wedge]
type = axisymmetric  ; the k faces (5, 6) of a 2-D axisymmetric (2Daxi) wedge; the axis takes symmetry

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

The two records of a connection face (`101`, `201`) in `bc.txt` must point at each other (each one's `bs is js ks fs` is the other's block, cell and face, not its own); ARES stops at set-up otherwise. The block across a connection face must have at least two cells across that face, at every multigrid level: the two ghost cells of the face are its first two cells, and their metrics come from its first two node planes past the face. A record whose block across has one cell stops the run at set-up, with the message that names the record, the block and its cells; blocks one cell thick across a connection are not supported. The flux of the face is computed once, by the record of the cell on the low-index side (between faces of the same parity, by the record with the smaller block, face, k, j, i, whatever the order of the records in the file), and given with the opposite sign to the other cell, also across MPI ranks (see [Spatial Discretization](../theory/numerics.md#faces-between-blocks)).

### Walls

`type = wall` is the viscous no-slip wall. The keys present select the thermal treatment; roughness `ks` (0 = smooth) and emissivity `eps` are optional add-ons:

| Keys present | Wall condition |
|--------------|----------------|
| `q` (`= 0` ⇒ adiabatic) | Prescribed **heat flux**, positive into the fluid on every face |
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

The wall output, `OUTPUT/wall.tec`, writes the same quantities on every face of every wall, oriented into the gas: the
stress `tauX`, `tauY`, `tauZ` is the traction of the gas on the wall, $\boldsymbol{\tau}\cdot\mathbf{n}$ with
$\mathbf{n}$ the unit normal into the gas (along the flow under an attached boundary layer); the heat flux `qw` is
positive into the gas, as `q` (a heat-flux wall writes its `q`, an isothermal wall the heat it conducts into the gas).

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
k     = 2.25e-3      ; freestream k       (two-equation models)
omega = 6250.0       ; freestream omega
```

Two optional groups apply to every inlet:

- **Flow direction** — `alpha` and `beta`; with no angles the inflow is face-normal. The two inlets read them differently: on the mass-flux inlet (`404`) they are the angles of the velocity in the $x$–$y$ and $x$–$z$ planes, $V_y/V_x = \tan\alpha$ and $V_z/V_x = \tan\beta$; on the supersonic inlet (`405`) $\alpha$ is an azimuth and $\beta$ an elevation, direction $(\cos\alpha\cos\beta, \sin\alpha\cos\beta, \sin\beta)$. The two agree when one of the angles is zero.
- **Turbulence inflow** — `nuTilde` for the SA models, `k` + `omega` for the two-equation models, `Rij` + `omega` for the RSM. `rf` is a relaxation factor (default 1).

!!! info "Boundary turbulence values are specific"
    `nuTilde` ($\tilde\nu$, m²/s), `k` (m²/s²), `omega` (1/s) and `Rij` ($R_{11} = R_{22} = R_{33}$,
    m²/s²; the shear stresses are set to $10^{-8}$) are the **specific** quantities of the NASA
    [Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/): ARES multiplies each by the
    density of the boundary state it builds, since the variables it transports are the
    density-weighted ones. The [`[ICB-*]` keys](initial-conditions.md#turbulence) are the opposite —
    `rhoNuTilde`, `rhoK`, `rhoOmega`, `rhoRij`, already weighted — and a boundary section carrying one
    of those is refused by ATLAS with a hint, as are the retired `mit` and `kappa`; `mut`/`nut` name
    the eddy viscosity, which is never an input.

### Outlets

`type = outlet` prescribes a back pressure; `rf` is the relaxation factor:

```ini
[outflow]
type = outlet
p    = 46.11e5       ; back (ambient) pressure [Pa]
```

An inlet (`404`, `405`) takes its face as a slip (symmetry) wall, for that evaluation, when the gas of its boundary cell is not moving into the domain — a cell at rest included: a mass-flux inlet takes nothing from the velocity of the cell, so a case started at rest, with the outlet at the initial pressure, never starts moving. An outlet (`406`) does the same when the gas of its cell moves into the domain. At the end of a run in which it happened one line gives the number of flux evaluations of each kind; start such a case with a velocity into the domain at the inlet. An outlet with `p = 0` reverts to extrapolation.

### Manifold

`type = manifold` joins faces through a lumped plenum: the gas that leaves the domain through the `outlet` faces of the
section comes back through its `inlet` faces, any number of each, on any blocks and ranks, matching or not cell by cell
(the turn-around between the two channels of a counter-flow element, a header feeding several channels). One section
per manifold; the roles are fixed by the section, not by the direction of the flow:

```ini
[turn]
type   = manifold
outlet = 1:2        ; block:face where the gas leaves the domain (comma-separated list)
inlet  = 2:2        ; block:face where it comes back
K      = 1.0        ; loss on the dynamic pressure of the outlet faces (default 0)
```

Once per stage ARES sums, over the cells of the outlet faces whose gas leaves, the mass flow $\dot m = \sum \rho\,
u_n\,A$ and, weighted by it, the total enthalpy $h_0 = h + v^2/2$, the dynamic pressure $q = \rho v^2/2$ and the
turbulence variables; over the cells of the inlet faces, the area $A_\text{in}$, the static pressure and their own dynamic
pressure. Then

- the inlet faces take the mass flux $g = \texttt{mdot-ratio}\;\dot m / A_\text{in}$, uniform, along the face normal,
  the total enthalpy $h_{0,\text{out}} - \texttt{dh0}$ and the mean turbulence; the static state of the face follows
  from the pressure of its cell and $h = h_0 - (g/\rho(p,h))^2/2$ through the fluid tables;
- the outlet faces take the static pressure $p_{0,\text{out}} - q_\text{out}$ of the total pressure
  $p_{0,\text{out}} = p_{0,\text{in}} / \texttt{p0-ratio} + K\,q_\text{out}$, with $p_{0,\text{in}}$ the mean static
  pressure of the inlet cells plus their mass-weighted dynamic pressure, through the outlet `406`.

With the defaults (`K` = 0, `mdot-ratio` = 1, `p0-ratio` = 1, `dh0` = 0) the plenum is loss-free and perfectly mixed:
mass, total enthalpy and total pressure cross it exactly. `K` = 1 loses the dynamic pressure of the gas that leaves, as
a jet into a plenum at rest does (the static pressure is then about the same on the two sides). The three ratios lump
the topology of a flow system on one link (one element simulated out of N, a bleed, the total-pressure ratio and the
specific work `dh0` [J/kg] of a turbine): with `mdot-ratio` ≠ 1 or `dh0` ≠ 0 the link is a mass or energy sink by
construction. `rf` (default 1) blends the value of the manifold with that of the face's own cell: the mass flux and
total enthalpy of the inlet faces, the static pressure of the outlet faces.

An outlet face whose cell has the gas moving inward is a slip wall for that evaluation, as on a `406`. While nothing
leaves through the outlet faces (a case started at rest), the inlet faces are slip walls and the outlet faces take the
mean static pressure of the inlet side. Unlike the inlet `404`, an inlet face of a manifold takes its mass flux whatever
the direction of the gas of its cell. The log has one line per manifold and grid level at the shell interval
(`shell-diter`): `manifold <id>, grid level <l>, iteration <n>: mdot out/in, h0 out/in, p0 out/in, taken out` (the
mass flows, total enthalpies and total pressures of the two sides, and the power $\dot m\,\texttt{dh0}$). The test
`Manifold` checks it ([Testing](../development/testing.md#ci-tests)).

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

All boundary cells are listed — one entry per face cell, face by face, block by block, as ATLAS writes them (another
order gives the same solution and the same wall output).

Types that carry no data (symmetry, axisymmetric, extrapolation) have no second line.

Every value of the data line is required, except those with a default value: `k_rough` (`ks`) and `eps_wall` (`eps`)
(0). An empty value (two commas with nothing between them) keeps that default; an empty required value stops the run at
set-up, with the line of `bc.txt`, the record, its type and the name of the value. So does a line short of a value, or
with a value that is not a number, and a value that is not finite (NaN or infinity), required or with a default, the
flow directions `alpha` and `beta` included. Blank lines between records are skipped. On a coarse level of a multigrid
run the messages name the file of that level.

### Numeric type codes

The codes ARES reads (`IO_BC`, `Mod_BC_Fluxes`), with the `[BCB]` type that produces them:

| Code | `[BCB]` `type` | Meaning | Data line |
|:----:|----------------|---------|-----------|
| `101` | `connection`    | Abutting inter-block interface | — |
| `102` | `chimera`       | Overset interface: the two ghost cells take the volume-fraction mean of the donor-cell primitives | `n1 n2`, then `n1+n2` lines `bD iD jD kD vol_frac` |
| `103` | —               | Multi-solver (HYDRA) coupling: symmetry wall + external flux | `bS iS jS kS fS d11 d12 d21 d22 [ks]` |
| `104` | —               | Multi-solver (HYDRA-AF) coupling on a non-conforming interface: as `103`, with the partner cells found from the donors | `n1 n2 [ks]`, then `n1+n2` lines `bD iD jD kD vol_frac` (FUSS donors) |
| `200` | `axisymmetric`  | Rotational periodicity of a wedge or sector: the two k faces (5 and 6) of one block | — |
| `201` | `periodic`      | Periodic face pair (dispatched as a connection) | — |
| `300` | `symmetry` / `wall` (no data) | Symmetry plane / inviscid slip wall | — |
| `301` | `wall` (`q`)    | No-slip wall, prescribed heat flux `qw`, positive into the fluid on every face (`0` = adiabatic) | `q, ks, eps` |
| `302` | `wall` (`T`)    | No-slip wall, prescribed temperature `Tw` (isothermal) | `T, ks, eps` |
| `400` | `extrapolation` | Zeroth-order extrapolation (supersonic outflow / far-field) | — |
| `404` | `inlet` (`g`,`T`) | Subsonic inlet: mass flux + static temperature | `T, g, α, β, rel[, turb]` |
| `405` | `inlet` (`mach`,…) | Supersonic inlet: Mach + static `T` + static `p` | `M, T, p, α, β, rel[, turb]` |
| `406` | `outlet` (`p`)  | Subsonic outlet: back pressure `pAmb` | `p, rel` |
| `501` | `manifold`      | Face of a lumped plenum: role 1 the gas leaves the domain, 2 it comes back ([Manifold](#manifold)) | `gid, role, K, rel, mdot_ratio, p0_ratio, dh0` |

where `α`, `β` are the inflow angles (`normal` = face-normal; `404`: $V_y/V_x = \tan\alpha$, $V_z/V_x = \tan\beta$; `405`: azimuth and elevation, see [Inlets](#inlets)) and `rel` a relaxation factor: required on the three lines, it is used by the supersonic inlet `405` only, which blends its Mach number and pressure with those of the cell, $\text{rel}\,M_\text{bc} + (1-\text{rel})\,M_\text{cell}$. A type that is not in this table stops the run at set-up with the line of `bc.txt`; so do a block or a face out of range, a cell outside its block or not next to the face it names, two records for the same cell and face, and a temperature or a pressure of a data line (`302` `T`; `404` `T`; `405` `p` and `T`; `406` `p` when not zero) outside the range of the fluid tables. The mass flux `g` of the inlet `404` is positive into the domain on every face: zero, or a negative value on a face 1, 3 or 5, stops the run at set-up (an inlet only lets gas in; give a face where the gas leaves an outlet type, `406` or `400`); on a face 2, 4 or 6 a negative value is taken as its magnitude. The supersonic inlet `405` with `α` and `β` both `normal` lets the gas in along the normal of the face into the domain, on every face; with one angle `normal`, that angle is the one of that normal (`α` its azimuth, `β` its elevation). On an inlet, `turb` is the trailing group of as many values as the model transports (1, 2 or 7), read by position as the **specific** quantities `nuTilde`, `k`, `omega`, `Rij` — see [Inlets](#inlets) above. An inlet or an outlet whose cell has the gas moving the wrong way, or at rest on an inlet, acts as a slip wall for that evaluation, and the end of the run counts the evaluations (see [Outlets](#outlets)).

### Example — mass-flux inlet (type 404)

A mass-flux inlet (`g = 80.65 kg/m²s`, `T = 300 K`) compiles to:

```
       1       1       1       1       1     404
    0.300000E+03,    0.806478E+02,         normal,         normal,    0.100000E+01,
```

The identifier line places the BC on block 1, cell `(1,1,1)`, face 1, with type `404` (inlet: mass flux + static temperature). The data line carries `T`, `g`, the in-plane and out-of-plane flow angles (`normal` = face-normal inflow), and a relaxation factor.

!!! note "Extrapolation (400)"
    A `400` face needs no data: it takes the flux of the state of its boundary cell, $\mathbf{F}(q_1)\cdot\mathbf{n}\,A$,
    whatever the direction of the gas. Where the gas of the boundary cell moves into the domain, it enters with the state of that
    cell: the characteristics that enter take the values of the cell, not data from outside, so a `400` face keeps any uniform
    inflow going by itself, whether the case wants it (the inlet of a truncated domain whose state outside is that of its
    boundary cells) or not (a free boundary through which the flow draws gas in, a plenum fed through its outlet). Where the gas
    outside has a known state, use an inlet (`404`, `405`).

    The run says where the gas enters through `400` faces. The first time it enters through the `400` faces of one face of a
    block (on each grid level), the run prints one line with the block, the face (and the grid level, with more than one),
    the iteration, the faces through which the gas enters at that stage and the largest $|w|/a$ among them, $w$ the velocity
    of the gas of the boundary cell along the normal out of the domain and $a$ its speed of sound. Each time the largest
    $|w|/a$ of that boundary at a stage is more than 10 times that of its last line, one more line follows, with that value
    and the faces above 10 times the last one: a boundary has at most one line for each order of magnitude that its inflow
    grows, not one for each stage. The lines are written at the end of the boundary fluxes of their stage, before the cells
    are updated, so a run that then stops on an error has them in its output. At the end of a run in which the gas entered
    through `400` faces, one more line reports the flux evaluations of `400` faces where it entered, over the stages, grid
    levels and ranks, and the largest $|w|/a$ among them; a run in which no gas entered, and a run that stops on an error, do
    not print it. The lines do not say whether or when the inflow ends; the end line says how many flux evaluations it lasted.
    After a restart the lines and the end line count from the restart, with the iteration of the restarted run. $|w|/a$ of
    about $10^{-9}$ or less: only round-off moves the gas towards the domain, as on the side faces of 1-D cases; from
    $10^{-6}$ to $10^{-3}$: small motions towards the domain, such as the numerical oscillations next to a shock or a junction
    between blocks, or a weak flow across the face; from $10^{-2}$ up: the gas really enters, at the inlet of a truncated
    domain or where the flow draws gas in through a free boundary.

!!! note "Type 200"
    With one cell in $k$, type `200` on faces 5 and 6 makes the mesh a 2-D axisymmetric wedge: the two faces lie at
    $\theta = \mp\delta/2$ and the ghost cell across each is the cell rotated by $\mp\delta$ about the $x$ axis. The
    inviscid flux gives the pressure and centrifugal terms; with the Navier–Stokes equations the viscous stress of the
    faces is applied as well ($-\tau_{\theta\theta}/r$ in the radial momentum, see
    [Spatial Discretization](../theory/numerics.md#axisymmetric-wedge)). Keep the wedge thin: the azimuthal
    derivatives carry an error of order $\delta^2$. With more than one cell in $k$, type `200` is a rotational
    periodicity: faces 5 and 6 make one interior face of the whole annulus, with its inviscid and viscous fluxes (see
    [Spatial Discretization](../theory/numerics.md#3-d-sector)). A viscous run on a wedge or a sector whose boundary
    faces on the axis have zero area stops at set-up, with the block, the cell and the face: so do axis nodes exactly at
    $r = 0$, and axis nodes at a small radius whose positions on the two $k$ planes coincide. Put the axis nodes at a
    small radius, $10^{-6}$ of the radial size of the first cell, turned with the other nodes of their $k$ plane.

    A face of type `200` takes what lies beyond it from the other $k$ face of its own block, in the same column
    $(i, j)$, turned by the angle of the wedge or of the sector: type `200` joins the two $k$ faces (5 and 6) of one
    block, and a block closed by type `200` holds the whole wedge or sector in $k$. A record of type `200` on a face 1 to
    4, with its cell outside its block, or without a record of type `200` on the other $k$ face of the same column stops
    the run at set-up with a message that names the record: so do a sector cut along $k$ into several blocks, a ring of
    360° in several blocks closed by type `200` between the last block and the first, and a wedge with type `200` on one
    $k$ face only. Where a wedge or a sector reaches the axis, the axis lies on one of the faces 1 to 4 and takes
    symmetry (`300`), as the message of a record of type `200` there says. Cut a sector along $i$ or $j$; join the
    blocks of a ring by block connections (`101`), and keep type `200` for a block that holds the whole sector or the
    whole turn.

!!! info "Numeric codes and ATLAS"
    In normal use you do not write these codes by hand — you write the `[BCB]` blocks above and let ATLAS translate them. The codes are documented here so you can interpret or hand-edit `bc.txt` when debugging a case.

---

## Multigrid

When geometric multigrid is enabled, ARES needs one boundary table per grid level. ATLAS produces the coarse-level tables alongside the finest one; the per-level iteration counts are set in `[ARES-Multigrid]` (see the [registry](registry.md#ares-multigrid)).
