# Gas Properties

ARES gets its thermodynamic and transport properties from a **pre-tabulated real-fluid table**, specified in two layers:

1. **The user-facing `[GPB-*]` interface** in `input.ini` — a *name-based* request for a fluid and a $(p,T)$ range. This is what you write.
2. **The solver-side tables** — `thermo.dat` and `transport.dat`, sampled on a pressure–enthalpy $(p,h)$ grid and read by ARES through the **FLINT** library. ATLAS compiles the `[GPB-*]` block into these tables.

This page covers both, **for a real-gas (real-fluid) phase** (`type = real-fluid`, internal phase type `RF`). For the physics of what the table holds and how it is used, see [Real-Fluid Thermodynamics](../theory/thermo.md).

---

## The `[GPB]` interface

A single `[GPB-Phase<n>]` section requests one real-fluid table:

```ini
[GPB-Phase1]
type  = real-fluid     ; selects the real-fluid phase (RF)
name  = air            ; optional phase label (prefixes the output files)
fluid = air            ; the fluid, passed to the property backend
pmin  = 0.80e5         ; pressure range [Pa]
pmax  = 2.0e5
Tmin  = 280.0          ; temperature range [K]
Tmax  = 340.0
NP    = 200            ; pressure grid points
NH    = 200            ; enthalpy grid points
model = coolprop       ; property backend (optional)
```

| Key | Required | Meaning |
|-----|:--------:|---------|
| `type` | ✓ | `real-fluid` selects the real-gas phase (internal type `RF`) |
| `fluid` | ✓ | Fluid identifier passed to the backend (e.g. `air`, `parahydrogen`) |
| `pmin`, `pmax` | ✓ | Pressure range of the table [Pa] |
| `Tmin`, `Tmax` | ✓ | Temperature range of the table [K] |
| `name` | — | Optional phase label; prefixes the generated files |
| `NP` | — | Pressure grid points (default 200) |
| `NH` | — | Enthalpy grid points (default 200) |
| `model` | — | Property backend (default `coolprop`) |

### Fluid and backend

`fluid` names the substance; `model` chooses how its properties are evaluated:

| `model` | Backend | Use |
|---------|---------|-----|
| `coolprop` (default) | CoolProp reference EOS | High-accuracy properties; the bundled cases use `air` and `parahydrogen` |
| `redlich-kwong` | Redlich–Kwong cubic EOS (Cantera) | Cubic-EOS approximation |
| `peng-robinson` | Peng–Robinson cubic EOS (Cantera) | Cubic-EOS approximation |

### Table range

ATLAS samples the fluid over $[p_\text{min},p_\text{max}]$ and the enthalpy band corresponding to $[T_\text{min},T_\text{max}]$, on an $N_P \times N_H$ grid, and writes the tables. The enthalpy bounds are computed from the corners $(T_\text{min},p_\text{min})$ and $(T_\text{max},p_\text{max})$, so the table is stored directly in $(p,h)$.

!!! warning "Stay inside the table"
    Properties are only valid inside the tabulated $(p,h)$ box, whose range ARES prints at loading (pressure, enthalpy, temperature). An initial condition with a cell outside it, or a temperature or a pressure of `bc.txt` outside it (302, 404, 405, 406), stops the set-up with the cell or the line. During the run a state within one table cell below $p_\text{min}$ or $h_\text{min}$ is extrapolated from the first cell; a state farther out, or at or beyond $p_\text{max}$ or $h_\text{max}$, stops the run with the cell and its state. The wall temperature of a heat-flux wall (301) or of a wall coupled to the solid by hydra (103, 104), and the ghost states of the walls, are brought within the tables instead (for the state of the gas only: the coupled wall keeps its temperature in the heat flux), and counted at the end of the run. Choose the range to comfortably bracket the expected flow states, transients included.

### Resolution

`NP` and `NH` set the number of pressure and enthalpy grid points. Finer grids reduce interpolation error where properties vary steeply — e.g. across the pseudo-critical line of a supercritical fluid — at the cost of a larger table. The default is $200\times200$; the bundled supercritical-hydrogen ([HTD](../vv/htd.md)) case uses $800\times800$.

---

## The tables ARES reads

From the `[GPB-*]` block ATLAS builds two ASCII tables on the **pressure–enthalpy $(p,h)$ grid**:

| File | Contents |
|------|----------|
| `thermo.dat` | Equation of state — the eight variables density $\rho$, temperature $T$, $\partial\rho/\partial T$, $\partial\rho/\partial h$, specific heat $c_p$, entropy, $\partial\rho/\partial p$ and speed of sound |
| `transport.dat` | Transport — dynamic viscosity $\mu$ and thermal conductivity $k$ |

Both are indexed by $(p,h)$, so the solver recovers every state by a direct table lookup with no equation-of-state inversion at run time. This is why the [initial-condition](initial-conditions.md) and solution files store $p$ and $h$ as the primitive thermodynamic pair. ARES reads the tables through the **FLINT** library at start-up.

!!! note "Both tables must share the same $(p,h)$ orientation"
    `thermo.dat` and `transport.dat` are built together on the same grid and must stay consistent — if you change the range or resolution, regenerate both with ATLAS.

See [Real-Fluid Thermodynamics](../theory/thermo.md) for the full table layout and the inverse table $h(p,T)$ ("Thermo inversion") built at loading.
