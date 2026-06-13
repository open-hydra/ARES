# Turbulent Prandtl Correction

**Case:** `test/Prt-correction`

This case validates the **wall-roughness turbulent-Prandtl correction** — the model in `src/lib/physics/Lib_Prt_Correction.f90`. Physically, a rough wall enhances momentum transfer more than heat transfer, so a constant turbulent Prandtl number $\mathrm{Pr}_t$ over-predicts the wall heat flux on rough surfaces. The correction adds a roughness-dependent increment $\Delta\mathrm{Pr}_t$ near the wall; this case isolates and verifies that effect on a heated, axisymmetric, rough pipe.

---

## Configuration

```ini
[ARES-Parameters]
simulation-type = turbulent

[ARES-Numerics]
cfl                     = 0.2
vnn                     = 0.2
time-scheme             = RK2
integration-variables   = prec
riemann-solver          = HLLC Prec
preconditioning-Uref    = 120.0
preconditioning-eps-min = 0.1

[ARES-RANS]
turbulence-model = SA-rough
Prt              = 0.9
Prt-correction   = .true.        ; <-- the feature under test

[GRIB-meshgen]
type = 2Daxi                     ! axisymmetric pipe

[ICB-Block1]
p   = 18.0e5
h   = 1.1e5
u   = 12.0
mit = 1.0e-10
```

The wall is a **heated rough wall**: a prescribed heat flux plus a sand-grain roughness height.

```ini
[qw]
type = wall
q    = 3.0d6        ; wall heat flux [W/m²]
ks   = 1.60d-4      ; sand-grain roughness height [m]
```

!!! note "Constant-property real-fluid table"
    Unlike the other cases, Prt-correction does **not** request a `[GPB]` real-fluid table from ATLAS. Instead `make-ares-table.py` builds a *constant-property* table (water, $\rho = 1002.4\ \mathrm{kg/m^3}$, constant $c_p$ and transport) in the ARES/ATLAS `(p,h)` table format. Holding the properties constant removes real-fluid variation from the comparison, so the only thing being tested is the **turbulent-Prandtl roughness correction** itself.

---

## Verification

```bash
cd test/Prt-correction
python3 make-ares-table.py          # build the constant-property table
./ARES.sh solve -b -p 4
../common/extract1d OUTPUT/field.tec OUTPUT/wall.tec OUTPUT/1d.dat INPUT
python3 validate_Prt_correction.py
```

`validate_Prt_correction.py` reads the section-averaged profiles in `OUTPUT/1d.dat` (produced with the shared `extract1d` tool) and compares the run with the rough-pipe correlations of the reference paper (Latini, Fiore, Nasuti, *Aerosp. Sci. Technol.* 126, 2022): the non-dimensional Darcy friction factor $f_D/f_{D,theo}$ (Colebrook–White) and Nusselt number $\mathrm{Nu}/\mathrm{Nu}_{theo}$ (Dipprey–Sabersky), in the style of the paper's Figs. 5 and 9. With the correction active, $\mathrm{Nu}/\mathrm{Nu}_{theo} \approx 1$; without it the heat transfer is over-predicted by a factor that grows with the relative roughness.

---

## What this validates

- The `Prt-correction` model (`Lib_Prt_Correction.f90`) and its coupling to the [effective conductivity](../theory/thermo.md#transport-properties).
- The `SA-rough` model's roughness shift via the sand-grain height `ks`.
- The heated-wall (`q`) boundary condition on an axisymmetric geometry.

!!! tip "Pairing with the model"
    `Prt-correction = .true.` only has an effect together with a rough-wall turbulence model (`SA-rough` here) and a non-zero `ks`. On a smooth wall the correction reduces to the standard constant $\mathrm{Pr}_t$.
