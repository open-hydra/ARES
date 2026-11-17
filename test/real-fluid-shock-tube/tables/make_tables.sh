#!/usr/bin/env bash
# Rebuilds thermo.dat, transport.dat and phase.txt from input.ini with ATLAS GPB (CoolProp).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
ATLAS GPB
for f in phase.txt thermo.dat transport.dat; do mv "fromATLAStoSolver/CO2-$f" "$f"; done
rm -rf fromATLAStoSolver
