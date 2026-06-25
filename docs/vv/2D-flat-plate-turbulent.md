# Turbulent Flat Plate (NASA TMR)

**Cases:** `test/Flat_Plate_SA`, `test/Flat_Plate_SST`, `test/Flat_Plate_Wilcox2006`, `test/Flat_Plate_SGGLRR`

The turbulent flat plate is the standard RANS validation case. ARES is run with each of its turbulence models on the same flow and compared against the **[NASA Turbulence Modeling Resource](https://turbmodels.larc.nasa.gov/) (TMR)** flat-plate reference data.

| Case directory | `turbulence-model` | Verification script | Reference data |
|----------------|--------------------|---------------------|----------------|
| `Flat_Plate_SA` | `SA` | `validate_sa.py` (reads `OUTPUT/1d.dat`) | `cf.dat` |
| `Flat_Plate_SST` | `SST` | `validate_sst.py` | `cf.dat`, `kw.dat`, `mut.dat` |
| `Flat_Plate_Wilcox2006` | `Wilcox2006` | `validate_wilcox2006.py` | `cf.dat`, `mut.dat` |
| `Flat_Plate_SGGLRR` | `SSGLRR` | `validate_sgglrr.py` | `cf.dat`, `mut.dat` |

Running the same geometry through all four models exercises the one-equation, both two-equation, and the Reynolds-stress closures against a common benchmark.

---

## Problem definition

A zero-pressure-gradient turbulent boundary layer over a flat plate, with `simulation-type = turbulent` and a `[ARES-RANS]` block selecting the model. The working fluid is air from a real-fluid $(p,h)$ table; the inflow carries the freestream turbulence values needed by the model — `kappa`/`omega` for the two-equation models, `rhoRij`/`omega` for the RSM, `mit` for SA.

The lower boundary is split into an upstream **symmetry** segment (ahead of the leading edge) and a downstream **adiabatic no-slip wall** (the plate), using the patch mechanism:

```ini
[ARES-Parameters]
simulation-type = turbulent

[ARES-RANS]
turbulence-model = SST     ; or SA / Wilcox2006 / SSGLRR
Prt = 0.85

[wall]
direction = x
patch1 = symmetry   ;  x ∈ [-2, 0]
range1 = -2. 0.
patch2 = adiabatic  ;  x ∈ [0, 2]   (the plate)
range2 = 0. 2.
```

---

## Reference quantities

`validate_sst.py` (and its siblings) reproduce the quantities plotted in the case's `PLOT.lay`, comparing against the NASA TMR data in `reference/`:

1. **Skin-friction coefficient** $c_f(x)$ vs. `reference/cf.dat`. Because the run may not sit exactly at the TMR Reynolds number, the comparison is done in $\mathrm{Re}_x$ space.
2. **Turbulence profiles** at $x \approx 0.97$ vs. `reference/kw.dat`, in wall-scaled form:

$$
k^+ = \frac{k}{a_\infty^2},
\qquad
\omega^+ = \frac{\omega\,\mu_\infty}{\rho_\infty a_\infty^2}.
$$

!!! note "Reynolds consistency"
    The $c_f(\mathrm{Re}_x)$ curve collapses onto the reference regardless of the exact run Reynolds number, but the $x\approx0.97$ profiles only collapse if the run is at the TMR Reynolds number ($\mathrm{Re}=5\times10^{6}$). The verification scripts document this so a profile mismatch is not misread as a model error.

---

## Solution variables

The turbulent `field.tec` carries the real-fluid primitive set plus the model and derived quantities:

```
x y z (nodal) | p u v w h kappa omega T rho sound mil kl mit (cell-centred)
```

where `kappa` = $k$, `omega` = $\omega$, `mil` = $\mu_\ell$, `kl` = $k_\ell$, `mit` = $\mu_t$. Wall output (`wall.tec`) provides the skin friction used for $c_f$.

---

## Running and verifying

```bash
cd test/Flat_Plate_SST
./ARES.sh solve -p 4
python3 validate_sst.py
```

The script reads `OUTPUT/field.tec` and `OUTPUT/wall.tec`, computes $c_f(x)$ and the wall-scaled $k^+,\omega^+$ profiles, overlays them on the NASA TMR reference, and prints the RMS / max relative $c_f$ error with a **PASS/FAIL** verdict (PASS for RMS < 5 %). The shared comparison logic lives in `test/common/_fp_turb_common.py`.

The SA case is the exception: `validate_sa.py` validates $c_f$ from `OUTPUT/1d.dat`, which must first be generated with the shared extractor:

```bash
cd test/Flat_Plate_SA
../common/extract1d OUTPUT/field.tec OUTPUT/wall.tec OUTPUT/1d.dat INPUT
python3 validate_sa.py
```

!!! tip "Constant-property comparison"
    The `Flat_Plate_SST`, `Flat_Plate_Wilcox2006`, and `Flat_Plate_SGGLRR` cases ship `set_constant_cp.py` and `set_constant_transport.py`, which flatten the real-fluid table to constant $c_p$ / transport so the comparison against the constant-property NASA reference is exact rather than approximate.

---

## What this validates

- Each RANS model reproduces the reference skin-friction law of the flat-plate boundary layer.
- The near-wall behaviour of $k$ and $\omega$ (or $\tilde\nu$, or the Reynolds stresses $R_{ij}$) matches the established reference profiles.
- The wall-output machinery (skin friction, $y^+$) is correct.
