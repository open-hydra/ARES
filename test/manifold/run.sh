#!/usr/bin/env bash
#===============================================================================
#  The manifold (BC 501): faces joined by a lumped plenum, real fluid
#
#  A duct of supercritical CO2 (make_case.py; the table of
#  ../real-fluid-shock-tube/tables) cut in two blocks joined by a manifold,
#  straight or with the second block turned back (the two face-2 ends joined,
#  the turn-around of a counter-flow element in small). Euler, MUSCL van Leer,
#  HLLC, RK3, steady (local time steps), from the uniform state with the
#  pressure of block 1 raised by 2 %, to a residual of 1e-12. The uniform state
#  of the inlet (T, g) and of the back pressure is the steady solution of every
#  layout, so:
#    1. the manifold runs converge to the one-block run (p, u, h within 1e-9;
#       measured 1e-13), and the log line of the manifold shows the same mass
#       flow, total enthalpy and total pressure on its two sides (1e-10);
#    2. the ratios of the link act: mdot-ratio 0.5 halves the mass flux of
#       block 2, dh0 lowers its total enthalpy by dh0, and the log line shows
#       the mass flow, total enthalpy and total pressure in the relations of
#       the record (mdot-ratio 0.5, K 0.5, p0-ratio 0.9, dh0 20 kJ/kg).
#  With the argument "mpi" (test ManifoldMPI): the two manifold layouts on 1
#  and on 2 ranks (one block each, the two sides of the manifold on two
#  ranks), bit for bit the same; $ARES_MPIEXEC (default mpiexec) runs them.
#  The binary is $ARES_BIN, or bin/ARES of the source tree.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
BIN="${ARES_BIN:-$HERE/../../bin/ARES}"
[[ -x "$BIN" ]] || { echo "no ARES binary at $BIN (set ARES_BIN)"; exit 1; }
MODE="${1:-serial}"
MPIEXEC="${ARES_MPIEXEC:-mpiexec}"
NAME=Manifold; [[ $MODE == mpi ]] && NAME=ManifoldMPI

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false

fail=0
check() {
  if [[ $1 -eq 0 ]]; then echo "  [ok  ] $2"; else echo "  [FAIL] $2"; fail=1; fi
}

# run <tag> <ranks> <layout> [manifold parameters]: the case in runs/<tag>, its log runs/<tag>/log
run() {
  local tag=$1 ranks=$2 layout=$3; shift 3
  local dir=runs/$tag
  rm -rf "$dir"
  mkdir -p "$dir/INPUT" "$dir/OUTPUT"
  for f in phase.txt thermo.dat transport.dat; do ln -s "../../../../real-fluid-shock-tube/tables/$f" "$dir/INPUT/$f"; done
  (cd "$dir" && python3 -B ../../make_case.py case "$layout" "$@") || return 1
  cat > "$dir/input.ini" <<EOF
[ARES-Parameters]
newrun = .true.
iter-threshold = 100000
res-threshold = 1e-12

[ARES-Physics]
equations = euler

[ARES-Numerics]
cfl = 0.8
time-scheme = RK3
time-accurate = .false.
integration-variables = prim
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[ARES-IO]
sol-diter = 1000000000
shell-diter = 500
EOF
  if [[ $ranks -gt 1 ]]; then
    (cd "$dir" && "$MPIEXEC" -n "$ranks" "$BIN" > log 2>&1)
  else
    (cd "$dir" && "$BIN" > log 2>&1)
  fi
  if [[ ! -f $dir/OUTPUT/field.tec ]] || grep -a -qiE 'NaN|error' "$dir/log"; then
    echo "  [FAIL] the run $tag did not complete (see $dir/log)"
    return 1
  fi
  return 0
}

if [[ $MODE == mpi ]]; then
  for layout in manifold reverse; do
    echo "-- $layout on 1 and on 2 ranks --"
    run "mpi-$layout-np1" 1 "$layout" || { fail=1; continue; }
    run "mpi-$layout-np2" 2 "$layout" || { fail=1; continue; }
    python3 -B make_case.py same "runs/mpi-$layout-np1/OUTPUT" "runs/mpi-$layout-np2/OUTPUT" 0
    check $? "$layout: 2 ranks give the field of 1 rank bit for bit"
    python3 -B make_case.py carried "runs/mpi-$layout-np2/log" "runs/mpi-$layout-np2/OUTPUT" 1e-10
    check $? "$layout on 2 ranks: mass, total enthalpy and total pressure cross the manifold"
  done
else
  echo "-- the steady state through the manifold, from block 1 at +2 % pressure --"
  run one 1 one || fail=1
  for layout in manifold reverse; do
    run "$layout" 1 "$layout" || { fail=1; continue; }
    python3 -B make_case.py same runs/one/OUTPUT "runs/$layout/OUTPUT" 1e-9
    check $? "$layout: the field of the one-block duct"
    python3 -B make_case.py carried "runs/$layout/log" "runs/$layout/OUTPUT" 1e-10
    check $? "$layout: mass, total enthalpy and total pressure cross the manifold"
  done

  echo "-- the ratios of the link: mdot-ratio 0.5, K 0.5, p0-ratio 0.9, dh0 20 kJ/kg --"
  if run ratios 1 manifold mr=0.5 K=0.5 pr=0.9 dh0=20000; then
    python3 -B make_case.py split runs/ratios/OUTPUT 0.5
    check $? "the mass flux of block 2 is half that of block 1"
    python3 -B make_case.py drop runs/ratios/OUTPUT 20000
    check $? "the total enthalpy of block 2 is 20 kJ/kg below that of block 1"
    python3 -B make_case.py carried runs/ratios/log runs/ratios/OUTPUT 1e-10 mr=0.5 K=0.5 pr=0.9 dh0=20000
    check $? "the log line holds the relations of the record"
  else
    fail=1
  fi
fi

if [[ $fail -ne 0 ]]; then echo "$NAME: FAILED"; exit 1; fi
echo "$NAME: all checks passed"
