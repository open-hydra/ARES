#!/usr/bin/env bash
#===============================================================================
#  Preconditioning on a real fluid: inviscid CO2 past a cylinder at low Mach
#
#  The case of make_case.py: an O-grid of 96 x 48 cells, CO2 at 8 MPa, 350 K
#  (Z = 0.74, the table of tables/), first-order steady runs of 3000 iterations
#  (forward Euler, local time steps, cfl 0.8), run side by side:
#    prec-m2  integration-variables = prec, HLLC Prec, Mach 0.01,  eps-min 1e-3
#    prec-m3  the same at Mach 0.001, eps-min 1e-4
#    prim-m2  integration-variables = prim, HLLC, Mach 0.01
#  Checks (measured values with gfortran 13 in brackets):
#    1. the preconditioned runs converge: the pressure residual falls below
#       5e-2 of its first value (2.0e-2; the wake keeps it from falling much
#       further);
#    2. they give the same Cp within 5e-3 for r < 1.5 R (1.1e-3): the
#       incompressible limit, Mach-independent to O(M^2);
#    3. on the front half of the cylinder their Cp is within 0.75 of potential
#       flow (0.55: the dissipation of the first-order scheme separates the
#       flow at the back, and the wake lowers the suction at the shoulders);
#    4. the unpreconditioned run differs from it by more than 2 in Cp (6.6):
#       its dissipation scales with the speed of sound, Cp errors of order 1/M.
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

# setup <tag> <Mach> <integration variables> <riemann solver> <eps-min>: runs/<tag> with its inputs
setup() {
  local tag=$1 dir=runs/$1
  rm -rf "$dir"
  mkdir -p "$dir/INPUT" "$dir/OUTPUT"
  for f in phase.txt thermo.dat transport.dat; do ln -s "../../../tables/$f" "$dir/INPUT/$f"; done
  (cd "$dir" && MACH=$2 python3 -B ../../make_case.py case) || return 1
  cat > "$dir/input.ini" <<EOF
[ARES-Parameters]
iter-threshold = 3000
newrun = .true.
res-threshold = 1e-300

[ARES-Physics]
equations = euler

[ARES-Numerics]
cfl = 0.8
time-scheme = euler
time-accurate = .false.
integration-variables = $3
space-reconstruction = first-order
riemann-solver = $4
preconditioning-eps-min = $5

[ARES-IO]
sol-diter = 1000000000
shell-diter = 1000000000
EOF
}

# solved <tag>: the run of runs/<tag> completed
solved() {
  local dir=runs/$1
  if [[ -f $dir/OUTPUT/field.tec ]] && ! grep -a -qiE 'NaN|error' "$dir/log"; then return 0; fi
  echo "  [FAIL] the run $1 did not complete (see $dir/log)"
  fail=1
  return 1
}

setup prec-m2 0.01  prec "HLLC Prec" 1e-3 || exit 1
setup prec-m3 0.001 prec "HLLC Prec" 1e-4 || exit 1
setup prim-m2 0.01  prim HLLC        0.1  || exit 1
for tag in prec-m2 prec-m3 prim-m2; do
  (cd runs/$tag && "$BIN" > log 2>&1) &
done
wait

M="python3 -B make_case.py"
echo "-- preconditioned runs at Mach 0.01 and 0.001 --"
ok2=0; ok3=0
solved prec-m2 && ok2=1
solved prec-m3 && ok3=1
for tag in prec-m2 prec-m3; do
  [[ -f runs/$tag/OUTPUT/residual-history.dat ]] || continue
  $M converged runs/$tag/OUTPUT/residual-history.dat 5e-2
  check $? "$tag converges"
done
if [[ $ok2 == 1 && $ok3 == 1 ]]; then
  $M same runs/prec-m2/OUTPUT/field.tec 0.01 runs/prec-m3/OUTPUT/field.tec 0.001 5e-3
  check $? "Cp does not depend on the Mach number (incompressible limit)"
fi
for tag in prec-m2 prec-m3; do
  [[ $tag == prec-m2 && $ok2 == 1 || $tag == prec-m3 && $ok3 == 1 ]] || continue
  mach=0.01; [[ $tag == prec-m3 ]] && mach=0.001
  $M front runs/$tag/OUTPUT/field.tec $mach 0.75
  check $? "$tag: the front half is within 0.75 of potential flow"
done

echo "-- unpreconditioned run at Mach 0.01 --"
if solved prim-m2 && [[ $ok2 == 1 ]]; then
  $M differs runs/prim-m2/OUTPUT/field.tec 0.01 runs/prec-m2/OUTPUT/field.tec 0.01 2
  check $? "prim with HLLC loses the incompressible limit (Cp off by more than 2)"
fi

if [[ $fail -ne 0 ]]; then echo "LowMachCylinder: FAILED"; exit 1; fi
echo "LowMachCylinder: all checks passed"
