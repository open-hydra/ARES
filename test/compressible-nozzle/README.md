# Non-ideal supersonic nozzle (MDM)

Test M1.5H of Spinelli et al. (2018), the SU2 tutorial *NICFD_nozzle*. Results: `docs/vv/nicfd-nozzle.md`.

## Files

| Path | Content |
|------|---------|
| `input.ini` | ARES, GPB, ICB and BCB setup |
| `custom.py` | MDM for GPB `thermo = custom`: CoolProp HEOS (gas phase imposed) + SU2 constant mu, k |
| `MESH/mesh_SU2.dat` | SU2 tutorial mesh, Tecplot FE quads |
| `MESH/su2tec2ares.py` | rebuilds the 83 x 60 structured grid, writes `MESH/mesh.tec` (planar, one cell in z) |
| `IC/make_ic.py` | quasi-1D isentropic solution: `IC/{p,T,u}.dat` and the inlet `g`, `T` |
| `reference/verify.py` | centreline P/PT and Mach vs experiment and SU2, writes the SVG figures |
| `reference/experimentald_data.txt`, `reference/SU2_data.txt` | experimental and SU2 data |

## Setup notes

- **Planar, not axisymmetric:** the TROVA nozzle is planar. Faces 5 and 6 are `symmetry`.
- **Inlet:** ARES has no total-pressure inlet. `g = 1006.41`, `T = 541.925` are the choked mass flux and static
  temperature from (P_T, T_T) (`IC/make_ic.py`). The inlet total pressure settles 0.3% above 904388 Pa.
- **Tables:** `pmax = 1.3e6` stays below p_c = 1.41 MPa but well above P_T. The wall stagnation at start-up exceeds
  9.6e5. The high-p, low-h corner is metastable vapour, never visited.
- **ICB:** a real-fluid block takes `vel-file`, not `u-file`.

## Run

```bash
source $(conda info --base)/etc/profile.d/conda.sh && conda activate ct-env   # CoolProp, matplotlib
python3 MESH/su2tec2ares.py MESH/mesh_SU2.dat MESH/mesh.tec planar 1e-3
python3 IC/make_ic.py
cp MESH/mesh.tec mesh.tec
ATLAS GPB BCB ICB
for f in bc.txt ic.tec phase.txt thermo.dat transport.dat; do cp fromATLAStoSolver/MDM-$f INPUT/$f; done
./ARES.sh -b -p 4 solve        # ~50k iterations to converge (~7 orders on the pressure residual), ~12 min
python3 reference/verify.py
```
