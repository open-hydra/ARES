# Low-Mach cylinder (CO2): preconditioning

CI test `LowMachCylinder` (`docs/development/testing.md`), the ARES counterpart of MOSE's
`test/2D/euler/low-mach-cylinder`, in a real fluid: CO2 at 8 MPa, 350 K (Z = 0.74).

| Run | `integration-variables` | `riemann-solver` | Mach | `preconditioning-eps-min` |
|---|---|---|---|---|
| `prec-m2` | `prec` | `HLLC Prec` | 0.01 | 1e-3 |
| `prec-m3` | `prec` | `HLLC Prec` | 0.001 | 1e-4 |
| `prim-m2` | `prim` | `HLLC` | 0.01 | — |

All steady, first order, forward Euler at cfl 0.8, 3000 iterations, run side by side in `runs/`. ARES has no
total-pressure inlet: the upstream third of the far field is a mass-flux inlet (404, the normal mass flux
rho U |n_x| per face, flow along x), the rest an outlet (406 at p_inf).

| Path | Content |
|------|---------|
| `run.sh` | the three runs and the checks (listed in its header with the measured values) |
| `make_case.py` | the O-grid, `ic.tec`, `bc.txt`; the checks |
| `tables/` | the CO2 table (committed), its `[GPB]` block and `make_tables.sh` (ATLAS + CoolProp) |
| `plot_vv.py` | the figure of `docs/vv/low-mach-cylinder.md` from `runs/` (matplotlib) |

```bash
ARES_BIN=$PWD/../../bin/ARES ./run.sh
```
