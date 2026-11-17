"""Inputs of the turbulence smoke runs, and their checks (Python standard library only).

A coarse turbulent flat plate in air (the tables of tables/): p = 1e5 Pa, T = 300 K, U = 70 m/s, one block of NI x NJ
cells on x in [-X_LE, L], y in [0, H], one cell in z:
  face 1  mass-flux inlet (404, g = rho U, T, the free-stream turbulence of the model);
  face 2  outlet (406 at p);
  face 3  symmetry (300) for x < 0, adiabatic no-slip wall (301) from the leading edge x = 0, rough (ks) for the
          rough Spalart-Allmaras models;
  face 4  symmetry;  faces 5, 6 symmetry (planar).
The cells cluster at the leading edge and at the wall (first cell 2e-6 m, y+ about 0.5 at x = 0.5 m).
Free-stream turbulence after the NASA TMR flat plate: nuTilde = 3 nu; k = 9e-9 a^2 and mut/mu = 9e-3 for the
two-equation models; Rij = 2/3 k on the diagonal (1e-8 off it) and the same omega for SSG/LRR.

The initial field holds a turbulent boundary layer from the leading edge: u = U (y/delta)^(1/7) below
delta = 0.37 x Re_x^-0.2, an eddy viscosity nu_t = NUT_BL nu 4 (y/delta)(1 - y/delta) (at least that of the free stream),
and the turbulence variables of each model from it: nuTilde = nu_t (fv1 is about 1); k = u_tau^2 / 0.3 (Bradshaw),
u_tau^2 = U^2 cf/2 with cf = 0.0592 Re_x^-0.2, omega = k / nu_t; Rij = 2/3 k on the diagonal. A model has no time to
build a boundary layer from the free stream in a short run (mut/mu stays below 1 after 3000 iterations): these runs
are too short to converge, and they check that each model sets up, runs, keeps a finite field with non-negative
turbulence variables and eddy viscosity, and keeps the eddy viscosity of the layer alive and bounded.

Usage (from a run directory holding INPUT/):
  python3 -B make_case.py case MODEL          writes INPUT/ic.tec and INPUT/bc.txt for the model (nrans 1, 2 or 7)
  python3 -B make_case.py check F HISTORY     the checks of the solution F and its residual history
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from ares_io import Table, read_field  # noqa: E402

NI_LE, NI, NJ = 8, 40, 32          # cells ahead of the plate, in all along x, across
X_LE, L, H, DZ = 0.1, 1.0, 0.2, 0.01
DY1 = 2.0e-6
P, T, U = 1.0e5, 300.0, 70.0
KS = 1.0e-5                         # sand-grain roughness of the rough SA models [m]
NUT_BL = 50.0                       # peak nu_t / nu of the initial boundary layer


def nrans(model):
    return 7 if model.startswith('SSGLRR') else 1 if model.startswith('SA') else 2


def stretched(n, length, d1):
    """n + 1 points on [0, length], the first interval d1, growing geometrically."""
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        if d1 * (q ** n - 1.0) / (q - 1.0) > length:
            hi = q
        else:
            lo = q
    pts = [0.0]
    for m in range(n):
        pts.append(pts[-1] + d1 * q ** m)
    pts[-1] = length
    return pts


def grid():
    dx = 0.004
    xp = stretched(NI - NI_LE, L, dx)
    xu = [-v for v in reversed(stretched(NI_LE, X_LE, dx))]
    return xu[:-1] + xp, stretched(NJ, H, DY1)


_FREE = []


def free_stream():
    if _FREE:
        return _FREE[0]
    thermo = Table(os.path.join(HERE, 'tables', 'thermo.dat'))
    transport = Table(os.path.join(HERE, 'tables', 'transport.dat'))
    h = thermo.h_of(P, T)
    rho, a = thermo('Density', P, h), thermo('SoundSpeed', P, h)
    mu = transport('Viscosity', P, h)
    _FREE.append((h, rho, a, mu))
    return _FREE[0]


def turbulence(model, nut=None, k_bl=None, rho=None):
    """Specific (inlet) turbulence values in the order of ARES, from the free stream; with nut, k_bl and rho, the
    density-weighted ones of a cell of the initial field."""
    _, rho_inf, a, mu = free_stream()
    n = nrans(model)
    k_inf = 9.0e-9 * a * a
    nut_inf = 9.0e-3 * mu / rho_inf
    if nut is None:
        nut, k, r = (3.0 * mu / rho_inf if n == 1 else nut_inf), k_inf, 1.0
    else:
        nut, k, r = max(nut, nut_inf), max(k_bl, k_inf), rho
    if n == 1:
        return [r * nut]
    omega = k / nut
    return [r * v for v in ([k, omega] if n == 2 else [2.0 / 3.0 * k] * 3 + [1.0e-8] * 3 + [omega])]


def layer(x, y):
    """u, nu_t / nu and k of the initial boundary layer at (x, y)."""
    h, rho, _, mu = free_stream()
    if x <= 0.0:
        return U, 0.0, 0.0
    rex = rho * U * x / mu
    delta = 0.37 * x * rex ** -0.2
    if y >= delta:
        return U, 0.0, 0.0
    eta = y / delta
    utau2 = 0.5 * U * U * 0.0592 * rex ** -0.2
    return U * eta ** (1.0 / 7.0), NUT_BL * 4.0 * eta * (1.0 - eta), utau2 / 0.3 * (1.0 - eta)


def write_case(model):
    h, rho, _, mu = free_stream()
    spec = turbulence(model)
    xs, ys = grid()
    ni, nj = len(xs) - 1, len(ys) - 1
    xn = [x for _ in range(2) for _ in ys for x in xs]
    yn = [y for _ in range(2) for y in ys for _ in xs]
    zn = [z for z in (0.0, DZ) for _ in ys for _ in xs]
    names = {1: ['rhoNuTilde'], 2: ['rhoK', 'rhoOmega'],
             7: ["ru'u'", "rv'v'", "rw'w'", "ru'v'", "ru'w'", "rv'w'", 'rhoOmega']}[len(spec)]
    nv = 8 + len(spec)
    nc = ni * nj
    cells = []
    for j in range(nj):
        for i in range(ni):
            u, nr, k = layer(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]))
            cells.append((u, turbulence(model, nr * mu / rho, k, rho)))
    with open('INPUT/ic.tec', 'w') as f:
        f.write(' VARIABLES = x y z p u v w h %s\n' % ' '.join(names))
        f.write(' ZONE  T = B1-RF, I=%d, J=%d, K=2, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED), SOLUTIONTIME=0.0\n' % (ni + 1, nj + 1, nv))
        for col in (xn, yn, zn):
            f.write(''.join(' %.16E\n' % v for v in col))
        f.write(' %.16E\n' % P * nc)
        f.write(''.join(' %.16E\n' % c[0] for c in cells))
        f.write(' %.16E\n' % 0.0 * (2 * nc))
        f.write(' %.16E\n' % h * nc)
        for m in range(len(spec)):
            f.write(''.join(' %.16E\n' % c[1][m] for c in cells))
    rec = lambda i, j, face, t: ''.join('%8d' % v for v in (1, i, j, 1, face, t)) + '\n'
    turb = ''.join(' %.10E,' % v for v in spec)
    ks = KS if 'rough' in model else 0.0
    out = []
    for j in range(1, nj + 1):
        out += [rec(1, j, 1, 404), '  %.10E, %.10E, normal, normal, 1.0,%s\n' % (T, rho * U, turb)]
    for j in range(1, nj + 1):
        out += [rec(ni, j, 2, 406), '  %.10E, 1.0,\n' % P]
    for i in range(1, ni + 1):
        if 0.5 * (xs[i - 1] + xs[i]) < 0.0:
            out.append(rec(i, 1, 3, 300))
        else:
            out += [rec(i, 1, 3, 301), '  0.0, %.4E, 0.0,\n' % ks]
    for i in range(1, ni + 1):
        out.append(rec(i, nj, 4, 300))
    for face in (5, 6):
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                out.append(rec(i, j, face, 300))
    with open('INPUT/bc.txt', 'w') as f:
        f.write(''.join(out))
    print('  %s: %d x %d cells, nrans %d, inlet turbulence %s' % (model, ni, nj, len(spec),
                                                               ' '.join('%.4g' % v for v in spec)))


def check(path, history):
    (ni, nj, _), nodes, c = read_field(path)
    ok = True

    def verdict(cond, text):
        nonlocal ok
        print('  [%s] %s' % ('ok  ' if cond else 'FAIL', text))
        ok = ok and cond
    bad = [name for name, v in c.items() if any(not math.isfinite(x) for x in v)]
    verdict(not bad, 'every variable finite' + (' (not: %s)' % ' '.join(bad) if bad else ''))
    turb = [n for n in c if n.startswith('rho') and n != 'rho' or n in ("ru'u'", "rv'v'", "rw'w'")]
    for name in turb + ['mut']:
        lo = min(c[name])
        verdict(lo >= 0.0, '%s non-negative (smallest %.3e)' % (name, lo))
    ratio = max(m / l for m, l in zip(c['mut'], c['mil']))
    verdict(5.0 < ratio < 1.0e4, 'the eddy viscosity of the layer stays alive and bounded: largest mut/mu %.3g '
            '(initially %g, required between 5 and 1e4)' % (ratio, NUT_BL))
    rows = [line.split() for line in open(history) if line.strip()]
    first, last = float(rows[0][2]), float(rows[-1][2])
    verdict(math.isfinite(last) and last < first,
            'the first residual falls: %.3e -> %.3e over %d iterations' % (first, last, len(rows)))
    return ok


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    if sys.argv[1] == 'case':
        write_case(sys.argv[2])
    elif sys.argv[1] == 'check':
        sys.exit(0 if check(sys.argv[2], sys.argv[3]) else 1)
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
