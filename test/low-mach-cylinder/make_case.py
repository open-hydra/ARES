"""Inputs of the low-Mach cylinder test, and the checks of its output (Python standard library only).

Inviscid flow of CO2 past a circular cylinder at p_inf = 8 MPa, T_inf = 350 K (compressibility factor 0.74), the
real-fluid table of tables/ (ATLAS GPB, CoolProp). A 2-D O-grid of one block:
  i  around the cylinder, NI cells, clockwise from the downstream axis (i, j, k right-handed); faces 1 and 2 meet on
     the downstream axis and are joined by a connection (101) on themselves;
  j  radial, NJ cells from the cylinder (r = R, face 3, slip wall 300) to r = R_OUT (face 4), the first cell as wide
     as it is long and the cells growing geometrically;
  k  one cell, faces 5 and 6 symmetry (300): a planar case.
Face 4 is a mass-flux inlet (404: the normal mass flux rho U |n_x| and T_inf, flow along x) where its face centre has
cos(theta) < -0.5, the upstream third, and an outlet (406 at p_inf) elsewhere. The initial field is the free stream.

At Mach number M the flow is incompressible to O(M^2), whatever the equation of state: the pressure coefficient
Cp = (p - p_inf) / (rho_inf U^2 / 2) of a scheme with an incompressible limit does not depend on M. HLLC Prec with
integration-variables = prec has that limit when the floor of the reference velocity scales with M
(preconditioning-eps-min = M / 10); HLLC with prim has dissipation of order a and Cp errors of order 1/M.

The free-stream state is read from the table (bilinear in p, h, as ARES reads it), so that h_inf, rho_inf and a_inf are
those of the run, and the test needs neither ATLAS nor CoolProp.

Usage (from a run directory holding INPUT/; MACH in the environment, 0.01 by default):
  python3 -B make_case.py case                       writes INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py converged HISTORY FALL     exit 0 if the first residual of the last row of the history is
                                                     below FALL times that of its first row
  python3 -B make_case.py same F1 M1 F2 M2 TOL       exit 0 if the Cp of F1 (Mach M1) and F2 (M2) differ by less than
                                                     TOL in every cell with r < 1.5 R
  python3 -B make_case.py differs F1 M1 F2 M2 MIN    exit 0 if they differ by more than MIN in some cell there
  python3 -B make_case.py front F M TOL              exit 0 if the Cp of the first row of cells of the front half is
                                                     within TOL of potential flow
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from ares_io import Table, read_field as read_tecplot  # noqa: E402
TABLE = os.path.join(HERE, 'tables', 'thermo.dat')

NI, NJ = 96, 48
R, R_OUT = 0.5, 20.0
DZ = 0.01
P_INF, T_INF = 8.0e6, 350.0
MACH = float(os.environ.get('MACH', '0.01'))


def free_stream():
    tab = Table(TABLE)
    h = tab.h_of(P_INF, T_INF)
    rho, a = tab('Density', P_INF, h), tab('SoundSpeed', P_INF, h)
    return h, rho, a, MACH * a


# ---------------------------------------------------------------- grid and inputs
def radii():
    """NJ + 1 node radii from R to R_OUT, the first cell as wide as it is long."""
    d0 = 2.0 * math.pi * R / NI
    lo, hi = 1.0, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        if d0 * (q ** NJ - 1.0) / (q - 1.0) > R_OUT - R:
            hi = q
        else:
            lo = q
    r = [R]
    for j in range(NJ):
        r.append(r[-1] + d0 * q ** j)
    r[-1] = R_OUT
    return r


def write_case():
    h_inf, rho_inf, a_inf, u_inf = free_stream()
    r = radii()
    xs, ys = [], []
    for j in range(NJ + 1):
        for i in range(NI + 1):
            th = -2.0 * math.pi * i / NI
            xs.append(r[j] * math.cos(th))
            ys.append(r[j] * math.sin(th))
    nc = NI * NJ
    with open('INPUT/ic.tec', 'w') as f:
        f.write(' VARIABLES = x y z p u v w h\n')
        f.write(' ZONE  T = B1-RF, I=%d, J=%d, K=2, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n' % (NI + 1, NJ + 1))
        for col in (xs + xs, ys + ys, [0.0] * len(xs) + [DZ] * len(xs)):
            f.write(''.join(' %.16E\n' % v for v in col))
        for v in (P_INF, u_inf, 0.0, 0.0, h_inf):
            f.write(' %.16E\n' % v * nc)
    rec = lambda i, j, face, t: ''.join('%8d' % v for v in (1, i, j, 1, face, t)) + '\n'
    out = []
    for j in range(1, NJ + 1):                 # the seam: face 1 of cell (1, j) on face 2 of cell (NI, j)
        out += [rec(1, j, 1, 101), '%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (1, NI, j, 1, 2, 1, 0, 0, 1)]
        out += [rec(NI, j, 2, 101), '%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (1, 1, j, 1, 1, 1, 0, 0, 1)]
    for i in range(1, NI + 1):
        out.append(rec(i, 1, 3, 300))
    for i in range(1, NI + 1):
        th = -2.0 * math.pi * (i - 0.5) / NI
        if math.cos(th) < -0.5:                # 404: T, the normal mass flux, alpha = beta = 0 (along x), rel
            g = rho_inf * u_inf * abs(math.cos(th))
            out += [rec(i, NJ, 4, 404), '  %.17E, %.17E, 0.0, 0.0, 1.0,\n' % (T_INF, g)]
        else:                                  # 406: p, rel
            out += [rec(i, NJ, 4, 406), '  %.17E, 1.0,\n' % P_INF]
    for face in (5, 6):
        for j in range(1, NJ + 1):
            for i in range(1, NI + 1):
                out.append(rec(i, j, face, 300))
    with open('INPUT/bc.txt', 'w') as f:
        f.write(''.join(out))
    print('  free stream: p = %.6g Pa, T = %.6g K, h = %.8g J/kg, rho = %.6g kg/m3, a = %.6g m/s, U = %.6g m/s'
          % (P_INF, T_INF, h_inf, rho_inf, a_inf, u_inf))


# ---------------------------------------------------------------- checks
def read_field(path):
    """(p, x, y, j) of every cell of the solution, (x, y) the mean of the four nodes of the cell at k = 1."""
    (ni, nj, _), nodes, cells = read_tecplot(path)
    x, y, p = nodes['x'], nodes['y'], cells['p']
    out = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            idx = [i + ni * j, i + 1 + ni * j, i + ni * (j + 1), i + 1 + ni * (j + 1)]
            out.append((p[i + (ni - 1) * j], sum(x[m] for m in idx) / 4.0, sum(y[m] for m in idx) / 4.0, j))
    return out


def cp_field(path, mach):
    global MACH
    MACH = mach
    _, rho_inf, _, u_inf = free_stream()
    q = 0.5 * rho_inf * u_inf ** 2
    return [((p - P_INF) / q, cx, cy, j) for p, cx, cy, j in read_field(path)]


def residual_fall(path):
    rows = [line.split() for line in open(path) if line.strip()]
    return float(rows[-1][2]) / float(rows[0][2])


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    what, args = sys.argv[1], sys.argv[2:]
    if what == 'case':
        write_case()
        return
    if what == 'converged':
        fall = residual_fall(args[0])
        print('  first residual: last row / first row = %.3e (required below %s)' % (fall, args[1]))
        sys.exit(0 if fall < float(args[1]) else 1)
    if what in ('same', 'differs'):
        a = cp_field(args[0], float(args[1]))
        b = cp_field(args[2], float(args[3]))
        d = max(abs(x[0] - y[0]) for x, y in zip(a, b) if math.hypot(x[1], x[2]) < 1.5 * R)
        print('  largest |Cp(M = %s) - Cp(M = %s)| over r < 1.5 R: %.3e' % (args[1], args[3], d))
        if what == 'same':
            sys.exit(0 if d < float(args[4]) else 1)
        sys.exit(0 if d > float(args[4]) else 1)
    if what == 'front':
        worst = 0.0
        for cp, cx, cy, j in cp_field(args[0], float(args[1])):
            if j == 0 and cx < 0.0:
                r, th = math.hypot(cx, cy), math.atan2(cy, cx)
                worst = max(worst, abs(cp - (2.0 * (R / r) ** 2 * math.cos(2.0 * th) - (R / r) ** 4)))
        print('  largest |Cp - Cp_potential| on the first row of the front half: %.3f' % worst)
        sys.exit(0 if worst < float(args[2]) else 1)
    sys.exit(__doc__)


if __name__ == '__main__':
    main()
