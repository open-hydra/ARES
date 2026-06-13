# Verification & Validation

ARES ships with a suite of test cases that exercise the solver on real-fluid, wall-bounded turbulent flows with known reference or experimental data. Each case is a self-contained directory under `test/` with its own `input.ini`, run script (`ARES.sh`), reference data, and (where applicable) a Python verification script.

The current suite is entirely **2-D real-fluid turbulent** and is organised as follows.

<div class="grid cards" markdown>

-   :material-chart-line:{ .lg .middle } __Turbulent Flat Plate__

    ---

    Skin friction and turbulence profiles compared with the **NASA Turbulence Modeling Resource**, across four RANS models.

    [:octicons-arrow-right-24: Turbulent flat plate](2D-flat-plate-turbulent.md)

-   :material-pipe:{ .lg .middle } __Supercritical Heat Transfer (HTD)__

    ---

    Heat-transfer deterioration of supercritical **para-hydrogen** in a heated pipe, validated against experimental wall-temperature data.

    [:octicons-arrow-right-24: HTD case](htd.md)

-   :material-grain:{ .lg .middle } __Turbulent Prandtl Correction__

    ---

    The wall-roughness turbulent-Prandtl correction on a heated axisymmetric rough pipe.

    [:octicons-arrow-right-24: Prt correction](prt-correction.md)

</div>

---

## The current test set

| Case directory | Page | Fluid | Highlights |
|----------------|------|-------|-----------|
| `test/Flat_Plate_SA` | [Turbulent flat plate](2D-flat-plate-turbulent.md) | air | SA vs. NASA TMR |
| `test/Flat_Plate_SST` | [Turbulent flat plate](2D-flat-plate-turbulent.md) | air | SST vs. NASA TMR |
| `test/Flat_Plate_Wilcox2006` | [Turbulent flat plate](2D-flat-plate-turbulent.md) | air | Wilcox 2006 vs. NASA TMR |
| `test/Flat_Plate_SGGLRR` | [Turbulent flat plate](2D-flat-plate-turbulent.md) | air | SSG-LRR (RSM) vs. NASA TMR |
| `test/HTD` | [HTD](htd.md) | para-hydrogen | supercritical, preconditioned, 2-block axisymmetric |
| `test/Prt-correction` | [Prt correction](prt-correction.md) | water (constant-property table) | SA-rough + `Prt-correction` |

!!! note "These are the only cases"
    This list reflects exactly what is documented from the `test/` folder. The flat plates and HTD live directly under `test/`.

---

## What "verification" and "validation" mean here

- **Verification** checks that the equations are solved *correctly* — e.g. each turbulent flat plate reproduces the skin-friction law of the boundary layer on its own grid.
- **Validation** checks that the *models* reproduce reference or experimental data — e.g. the flat plates match the NASA Turbulence Modeling Resource, and HTD matches measured pipe wall temperatures for supercritical hydrogen.

Each case reads its `OUTPUT/`, compares against the `reference/` data (or against reference correlations / experimental points embedded in the script), and reports the error metrics. The turbulent flat plates share helper logic for the Reynolds-number-consistent comparison (`test/common/_fp_turb_common.py`); the cases that work on section-averaged 1-D profiles (`Flat_Plate_SA`, `HTD`, `Prt-correction`) use the shared `test/common/extract1d` tool to produce `OUTPUT/1d.dat` first.

!!! note "Real-fluid tables and the analytic checks"
    The canonical turbulence references (NASA TMR) assume a calorically simple gas. To compare cleanly, the `Flat_Plate_SST`, `Flat_Plate_Wilcox2006`, and `Flat_Plate_SGGLRR` cases ship `set_constant_cp.py` and `set_constant_transport.py`, which flatten the real-fluid $(p,h)$ table to constant $c_p$ and constant transport, isolating the numerics from real-fluid property variation. The [HTD](htd.md) case instead uses the full real-fluid table, since the property variation *is* the physics under test.
