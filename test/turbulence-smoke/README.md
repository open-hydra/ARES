# Turbulence smoke runs

CI tests `TurbulenceSmoke.<model>` (`docs/development/testing.md`): 300 iterations of each turbulence model on a coarse
flat plate in air (40 x 32 cells, 1e5 Pa, 300 K, 70 m/s), started from a turbulent boundary layer. They check that the
model runs (finite field, non-negative turbulence variables and mut, mut/mu between 5 and 1e4, falling residual), not
its accuracy: that is `test/flat-plate-*`. `SAcomp` is left out: ARES refuses it at set-up.

| Path | Content |
|------|---------|
| `run.sh [MODEL ...]` | the runs in `runs/<model>` and the checks; all models without arguments |
| `make_case.py` | grid, `ic.tec` (the boundary layer), `bc.txt`; the checks |
| `tables/` | the air table (committed), its `[GPB]` block and `make_tables.sh` (ATLAS + CoolProp) |

```bash
ARES_BIN=$PWD/../../bin/ARES ./run.sh SST SSGLRR
ITER=1000 ./run.sh SA           # longer runs
```
