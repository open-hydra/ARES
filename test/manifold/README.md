# Manifold (BC 501), real fluid

CI test `Manifold`, and `ManifoldMPI` in an MPI build (`docs/development/testing.md`). A duct of supercritical CO2
(7 MPa, 345 K, 2000 kg/m²s) cut in two blocks joined by a manifold: the converged field must be that of the one-block
duct, and the log line of the manifold must show mass, total enthalpy and total pressure crossing it, or the relations
of the record when the ratios of the link are set.

| Path | Content |
|------|---------|
| `run.sh` | the runs in `runs/` and the checks; `run.sh mpi`: the two layouts on 1 and 2 ranks, bit for bit |
| `make_case.py` | writes `ic.tec`, `bc.txt` of the layouts `one`, `manifold`, `reverse` (the turn-around of a counter-flow element); the checks |

The fluid table is that of `../real-fluid-shock-tube/tables` (committed: neither ATLAS nor CoolProp is needed).

```bash
ARES_BIN=$PWD/../../bin/ARES ./run.sh
ARES_BIN=$PWD/../../bin/ARES ARES_MPIEXEC=mpiexec ./run.sh mpi
```
