#!/usr/bin/env bash
#===============================================================================
#  Real-fluid thermodynamics: a shock tube in supercritical CO2
#
#  The case of make_case.py: 400 cells on [0, 1], at rest, 10 MPa, 360 K
#  (Z = 0.71) on the left and 5 MPa, 330 K (Z = 0.80) on the right, the table
#  of tables/; Euler, MUSCL (van Leer), HLLC, RK3 at cfl 0.8, up to t = 6e-4 s.
#  Checks against the exact solution from the CoolProp equation of state
#  (reference/exact.txt), measured values with gfortran 13 in brackets:
#    1. p, u, rho, T of the two plateaus (the rarefied and the shocked gas)
#       within 0.5% of the exact ones in every plateau cell (0.02%; an ideal
#       gas with the gamma of the initial states gives u* 16% lower);
#    2. the shock and the contact within one cell of their exact positions
#       (0.11 and 0.00 cells).
#  The binary is $ARES_BIN, or bin/ARES of the source tree.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
BIN="${ARES_BIN:-$HERE/../../bin/ARES}"
[[ -x "$BIN" ]] || { echo "no ARES binary at $BIN (set ARES_BIN)"; exit 1; }

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false

fail=0
check() {
  if [[ $1 -eq 0 ]]; then echo "  [ok  ] $2"; else echo "  [FAIL] $2"; fail=1; fi
}

dir=runs/muscl
rm -rf "$dir"
mkdir -p "$dir/INPUT" "$dir/OUTPUT"
for f in phase.txt thermo.dat transport.dat; do ln -s "../../../tables/$f" "$dir/INPUT/$f"; done
(cd "$dir" && python3 -B ../../make_case.py case 400) || exit 1
cat > "$dir/input.ini" <<EOF
[ARES-Parameters]
newrun = .true.
time-threshold = 6.0e-4

[ARES-Physics]
equations = euler

[ARES-Numerics]
cfl = 0.8
time-scheme = RK3
time-accurate = .true.
integration-variables = prim
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[ARES-IO]
sol-diter = 1000000000
shell-diter = 1000000000
EOF
(cd "$dir" && "$BIN" > log 2>&1)
if [[ ! -f $dir/OUTPUT/field.tec ]] || grep -a -qiE 'NaN|error' "$dir/log"; then
  echo "  [FAIL] the run did not complete (see $dir/log)"
  echo "RealFluidShockTube: FAILED"
  exit 1
fi

python3 -B make_case.py plateaus $dir/OUTPUT/field.tec 5e-3
check $? "the plateaus match the exact real-fluid solution"
python3 -B make_case.py shock $dir/OUTPUT/field.tec 1
check $? "the shock and the contact are where the exact solution puts them"

if [[ $fail -ne 0 ]]; then echo "RealFluidShockTube: FAILED"; exit 1; fi
echo "RealFluidShockTube: all checks passed"
