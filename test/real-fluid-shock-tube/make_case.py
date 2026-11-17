"""Inputs of the real-fluid shock tube, and the checks of its output (Python standard library only).

A Riemann problem in supercritical CO2: at rest, p = 10 MPa, T = 360 K (Z = 0.71) for x < 0.5 and p = 5 MPa, T = 330 K
(Z = 0.80) for x > 0.5, in a 1-D tube of N cells on [0, 1] (one cell in y and z, symmetry on the faces 3 to 6, slip
walls at the ends, which the waves do not reach). The tables of tables/ (ATLAS GPB, CoolProp); the initial enthalpies
are read from them, so the test needs neither ATLAS nor CoolProp.

The reference is the exact solution of reference/exact.txt (reference/exact_riemann.py, from the CoolProp equation of
state): a left rarefaction, a contact and a right shock. With the gamma = cp/cv of the initial states an ideal gas gives
a velocity u* 16% lower: the plateaus check the real-fluid thermodynamics of the table, not only the scheme.

Usage (from a run directory holding INPUT/):
  python3 -B make_case.py case [N]            writes INPUT/ic.tec and INPUT/bc.txt (N cells, 400 by default)
  python3 -B make_case.py plateaus F TOL      exit 0 if p, u, rho, T of both plateaus of F are within TOL (relative,
                                              the largest deviation of a cell from the exact value) of the exact ones
  python3 -B make_case.py shock F CELLS       exit 0 if the shock of F, where p crosses (p* + p_R)/2, is within CELLS
                                              cells of the exact position; the same for the contact (rho)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from ares_io import Table, read_field, solution_time  # noqa: E402

DY = 0.01


def exact():
    rows = [line.split() for line in open(os.path.join(HERE, 'reference', 'exact.txt')) if not line.startswith('#')]
    pL, TL, pR, TR, x0, t = map(float, rows[0])
    keys = ('p', 'u', 'rhoL', 'TL', 'rhoR', 'TR', 'head', 'tail', 'shock')
    star = dict(zip(keys, map(float, rows[1])))
    return (pL, TL, pR, TR, x0, t), star


def write_case(n):
    (pL, TL, pR, TR, x0, _), _ = exact()
    tab = Table(os.path.join(HERE, 'tables', 'thermo.dat'))
    hL, hR = tab.h_of(pL, TL), tab.h_of(pR, TR)
    xs = [m / n for m in range(n + 1)]
    xn = (xs * 2) * 2
    yn = ([0.0] * (n + 1) + [DY] * (n + 1)) * 2
    zn = [0.0] * 2 * (n + 1) + [DY] * 2 * (n + 1)
    xc = [0.5 * (xs[m] + xs[m + 1]) for m in range(n)]
    with open('INPUT/ic.tec', 'w') as f:
        f.write(' VARIABLES = x y z p u v w h\n')
        f.write(' ZONE  T = B1-RF, I=%d, J=2, K=2, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n' % (n + 1))
        for col in (xn, yn, zn):
            f.write(''.join(' %.16E\n' % v for v in col))
        f.write(''.join(' %.16E\n' % (pL if x < x0 else pR) for x in xc))
        f.write(' %.16E\n' % 0.0 * (3 * n))
        f.write(''.join(' %.16E\n' % (hL if x < x0 else hR) for x in xc))
    rec = lambda i, face: ''.join('%8d' % v for v in (1, i, 1, 1, face, 300)) + '\n'
    out = [rec(1, 1), rec(n, 2)]
    for face in (3, 4, 5, 6):
        out += [rec(i, face) for i in range(1, n + 1)]
    with open('INPUT/bc.txt', 'w') as f:
        f.write(''.join(out))
    print('  left  p = %.6g Pa, T = %.6g K, h = %.8g J/kg, rho = %.6g kg/m3' % (pL, TL, hL, tab('Density', pL, hL)))
    print('  right p = %.6g Pa, T = %.6g K, h = %.8g J/kg, rho = %.6g kg/m3' % (pR, TR, hR, tab('Density', pR, hR)))


def solution(path):
    (ni, _, _), nodes, cells = read_field(path)
    x = nodes['x'][:ni]
    xc = [0.5 * (x[m] + x[m + 1]) for m in range(ni - 1)]
    return xc, x[1] - x[0], cells


def plateaus(path, tol):
    (_, _, pR, _, x0, _), s = exact()
    t = solution_time(path)
    xc, dx, c = solution(path)
    xt, xm, xs = x0 + s['tail'] * t, x0 + s['u'] * t, x0 + s['shock'] * t
    margin = 0.15 * (xs - xt)                  # away from the smeared waves
    worst = 0.0
    for name, lo, hi, ref in (('left', xt, xm, {'p': s['p'], 'u': s['u'], 'rho': s['rhoL'], 'T': s['TL']}),
                              ('right', xm, xs, {'p': s['p'], 'u': s['u'], 'rho': s['rhoR'], 'T': s['TR']})):
        sel = [m for m, x in enumerate(xc) if lo + margin < x < hi - margin]
        for var, v in ref.items():
            dev = max(abs(c[var][m] / v - 1.0) for m in sel)
            worst = max(worst, dev)
            print('  %-5s plateau %-3s exact %.6g, largest deviation %.2f%% (%d cells)' % (name, var, v, 100 * dev, len(sel)))
    print('  largest deviation %.2f%% (required below %.2f%%)' % (100 * worst, 100 * tol))
    return worst < tol


def crossing(xc, f, level, lo, hi):
    """The x in (lo, hi) where f crosses level, linear between cells (the first crossing)."""
    for m in range(len(xc) - 1):
        if lo < xc[m] < hi and (f[m] - level) * (f[m + 1] - level) <= 0.0 and f[m] != f[m + 1]:
            return xc[m] + (level - f[m]) / (f[m + 1] - f[m]) * (xc[m + 1] - xc[m])
    return float('nan')


def shock(path, cells):
    (_, _, pR, _, x0, _), s = exact()
    t = solution_time(path)                    # the last step ends past the time-threshold
    xc, dx, c = solution(path)
    xm, xs = x0 + s['u'] * t, x0 + s['shock'] * t
    ok = True
    for name, var, level, ref in (('shock', 'p', 0.5 * (s['p'] + pR), xs),
                                  ('contact', 'rho', 0.5 * (s['rhoL'] + s['rhoR']), xm)):
        x = crossing(xc, c[var], level, 0.5 * (xm + x0 + s['tail'] * t), 1.0) if name == 'shock' else \
            crossing(xc, c[var], level, x0 + s['tail'] * t, 0.5 * (xm + xs))
        err = abs(x - ref) / dx
        print('  %-7s at x = %.5f, exact %.5f: %.2f cells (required below %s)' % (name, x, ref, err, cells))
        ok = ok and err < cells
    return ok


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    what, args = sys.argv[1], sys.argv[2:]
    if what == 'case':
        write_case(int(args[0]) if args else 400)
    elif what == 'plateaus':
        sys.exit(0 if plateaus(args[0], float(args[1])) else 1)
    elif what == 'shock':
        sys.exit(0 if shock(args[0], float(args[1])) else 1)
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
