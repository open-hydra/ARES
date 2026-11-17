#!/usr/bin/env bash
#===============================================================================
#  Turbulence models: short smoke runs on a coarse flat plate in air
#
#  ./run.sh [MODEL ...]   the models given, or all of them
#  The case of make_case.py (40 x 32 cells, air at 1e5 Pa, 300 K, 70 m/s),
#  ITER steady iterations (MUSCL, HLLC, RK2, cfl 0.4), one run per model in
#  runs/<model>. Checks: the field is finite, the turbulence variables and the
#  eddy viscosity are not negative, mut/mu exceeds 10 in the boundary layer,
#  the first residual falls. They do not check accuracy (the flat-plate
#  validations of test/flat-plate-* do).
#  The binary is $ARES_BIN, or bin/ARES of the source tree.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
BIN="${ARES_BIN:-$HERE/../../bin/ARES}"
[[ -x "$BIN" ]] || { echo "no ARES binary at $BIN (set ARES_BIN)"; exit 1; }
ITER=${ITER:-300}

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false

# SAcomp is left out: ARES refuses it at set-up (its correction depends on the units)
MODELS=(SA SA-R SA-RC SA-QCR2000 SA-rough SA-rough-QCR2000 SST SST-RC Wilcox2006
        SSGLRR SSGLRR-SD SSGLRR-w2019 SSGLRR-w2019-SD)
[[ $# -gt 0 ]] && MODELS=("$@")

fail=0
for model in "${MODELS[@]}"; do
  echo "-- $model --"
  dir=runs/$model
  rm -rf "$dir"
  mkdir -p "$dir/INPUT" "$dir/OUTPUT"
  for f in phase.txt thermo.dat transport.dat; do ln -s "../../../tables/$f" "$dir/INPUT/$f"; done
  (cd "$dir" && python3 -B ../../make_case.py case "$model") || { fail=1; continue; }
  extra=""
  [[ $model == SA-rough* ]] && extra="Prt-correction = .true."
  cat > "$dir/input.ini" <<EOF
[ARES-Parameters]
iter-threshold = $ITER
newrun = .true.
res-threshold = 1e-300

[ARES-Physics]
equations = navier-stokes
turbulence = $model

[ARES-Numerics]
cfl = 0.4
vnn = 0.4
time-scheme = RK2
time-accurate = .false.
integration-variables = prim
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[ARES-IO]
sol-diter = 1000000000
shell-diter = 1000000000
sol-variables = transport

[ARES-Turbulence]
Prt = 0.9
$extra
EOF
  (cd "$dir" && "$BIN" > log 2>&1)
  if [[ ! -f $dir/OUTPUT/field.tec ]] || grep -a -qiE 'NaN|error' "$dir/log"; then
    echo "  [FAIL] the run did not complete (see $dir/log)"
    fail=1
    continue
  fi
  python3 -B make_case.py check "$dir/OUTPUT/field.tec" "$dir/OUTPUT/residual-history.dat" || fail=1
done

if [[ $fail -ne 0 ]]; then echo "TurbulenceSmoke ${MODELS[*]}: FAILED"; exit 1; fi
echo "TurbulenceSmoke ${MODELS[*]}: all checks passed"
