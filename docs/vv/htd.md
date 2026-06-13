# Supercritical Heat Transfer (HTD)

**Case:** `test/HTD`

This case simulates **heat-transfer deterioration (HTD)** of **supercritical para-hydrogen** flowing through a heated pipe — the flagship real-fluid validation in the suite. Near the pseudo-critical point the properties of hydrogen vary so steeply with enthalpy that the wall heat transfer can *deteriorate*, producing a localised spike in wall temperature. Reproducing that spike requires the full real-fluid equation of state, which is exactly what ARES provides.

---

## Why this case matters

HTD is a stringent test because it stacks the distinguishing ARES features:

| Feature | Activated by |
|---------|--------------|
| Real-fluid EOS over a wide range | `[GPB] fluid = parahydrogen`, $p\in[30,80]$ bar, $T\in[25,800]$ K, $800\times800$ table |
| Low-Mach preconditioning | `integration-variables = prec`, `riemann-solver = HLLC Prec`, `preconditioning-Uref = 600` |
| Wall heating | a heated-wall boundary (`q = 4.766` MW/m²) on the downstream block |
| Multi-block, axisymmetric | `type = 2Daxi`, two blocks joined by a `connection` |

---

## Configuration

```ini
[ARES-Parameters]
simulation-type = turbulent

[ARES-Numerics]
cfl                     = 0.3
vnn                     = 0.3
time-scheme             = RK3
integration-variables   = prec
riemann-solver          = HLLC Prec
preconditioning-Uref    = 600.0
preconditioning-eps-min = 0.1

[ARES-RANS]
turbulence-model = SA
Prt = 0.9

[GPB-Phase1]
type  = real-fluid
fluid = parahydrogen
pmin  = 30.0d5   ;  pmax = 80.0d5     ! [Pa]
Tmin  = 25.0d0   ;  Tmax = 800.0d0    ! [K]
NP    = 800      ;  NH   = 800
```

The geometry is a two-block axisymmetric pipe: an **adiabatic entrance** block develops the flow, and a **heated** block applies the wall heat flux. The cryogenic supercritical inflow and the high back pressure are:

```ini
[inflow]
type = inlet
g    = 2741.0      ; mass flux [kg/m²s]
T    = 31.39       ; inlet temperature [K]  (cryogenic)

[outflow]
type = outlet
p    = 46.11d5     ; back pressure [Pa]  (above the H₂ critical pressure ≈ 12.8 bar)
```

!!! note "Genuinely supercritical"
    At ~46 bar and ~31 K the hydrogen is above its critical pressure and temperature, so it is a single supercritical phase whose density and specific heat change by an order of magnitude across the heated section. This is the regime the real-fluid table is built for.

---

## Reference data and verification

The verification script is `reference/validate_htd.py`. The experimental reference — **NASA test 24-1027** (*“Experimental heat-transfer results for cryogenic hydrogen flowing in tubes at subcritical and supercritical pressures to 800 psia”*) — is embedded directly in the script: measured **bulk temperature** and **wall temperature** at twelve axial stations along the heated tube.

The script reads:

| Input | Content |
|-------|---------|
| `OUTPUT/1d.dat` | Section-averaged 1-D profiles ($x$, $T_w$, $T$, …) extracted from the 2-D solution |
| `INPUT/thermo.dat` | The $(p,h)$ table, used for the optional bulk-temperature cross-check |

`OUTPUT/1d.dat` is produced by the shared **`extract1d`** tool (`test/common/Extract1D.f90`), which extracts wall and bulk 1-D profiles from `field.tec` / `wall.tec` using the real-fluid tables:

```bash
cd test/HTD
./ARES.sh solve -b -p 8                                          # long run; iter-threshold = 6,000,000
../common/extract1d OUTPUT/field.tec OUTPUT/wall.tec OUTPUT/1d.dat INPUT
cd reference && python3 validate_htd.py                          # overlays ARES on the NASA data
```

The script plots the bulk- and wall-temperature distributions along the tube against the experimental points; the HTD wall-temperature peak in the heated section is the feature under test.

---

## What this validates

- The **real-fluid $(p,h)$ table and thermo inversion** under steep property variation near the pseudo-critical line.
- The **low-Mach preconditioned** scheme on a genuine internal flow.
- The wall-heating boundary condition and the model's ability to capture the **wall-temperature peak** characteristic of heat-transfer deterioration.

!!! warning "Long run"
    HTD is a deep-convergence case (`iter-threshold = 6,000,000`, `cfl = 0.3`). Run it in the background and monitor `logfile`.
