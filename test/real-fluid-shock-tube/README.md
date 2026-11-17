# Real-fluid shock tube (CO2)

CI test `RealFluidShockTube` (`docs/development/testing.md`). A Riemann problem in supercritical CO2, checked against
its exact solution from the CoolProp equation of state.

| Path | Content |
|------|---------|
| `run.sh` | runs the case in `runs/muscl/` and the checks; exit 0 when all pass |
| `make_case.py` | writes `ic.tec`, `bc.txt` (1-D, 400 cells); the checks |
| `tables/` | the CO2 table (committed), its `[GPB]` block and `make_tables.sh` (ATLAS + CoolProp) |
| `reference/exact_riemann.py` | the exact solution (CoolProp, ~8 min), written to `reference/exact.txt` |
| `reference/plot_vv.py` | the figure of `docs/vv/real-fluid-shock-tube.md` from `runs/muscl` (matplotlib) |

```bash
ARES_BIN=$PWD/../../bin/ARES ./run.sh
```

Changing the initial states: edit `PL, TL, PR, TR` in `reference/exact_riemann.py` and rerun it (conda env `ct-env`),
then check the table range in `tables/input.ini` and rerun `tables/make_tables.sh`. The isentrope of the left state
must stay out of the two-phase region (10 MPa, 360 K reaches 304 K at 5 MPa, saturation 287 K).
