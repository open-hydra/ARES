"""Inputs of the Manifold test (BC 501), and the checks of its output (Python standard library only).

A duct along x, x in [0, 1] m, N = 40 cubic cells of 25 mm, one layer in y and z, supercritical CO2 of the table of
../real-fluid-shock-tube/tables (ATLAS GPB, CoolProp): a mass-flux inlet 404 (T = 345 K, g = 2000 kg/(m2 s)) at x = 0, a
back-pressure outlet 406 (7 MPa) at x = 1 m, symmetry (300) on the j and k faces. Euler: the uniform state of T, g and
the back pressure (rho 141.6 kg/m3, u 14.1 m/s, Mach 0.056) is the steady solution of every layout of the duct:
  one       one block of N cells;
  manifold  two blocks of N/2 cells joined by a manifold (501): face 2 of block 1 is its outlet (the gas leaves the
            domain there), face 1 of block 2 its inlet;
  reverse   the manifold between the two face-2 ends: block 2 is written from x = 1 down to x = 0.5 (j reversed too, so
            that the mesh stays right-handed), its face 1 is the outlet 406 and its face 2 the manifold's inlet; the gas
            still flows toward +x, that is toward decreasing i. The face pairing of a counter-flow element's turn-around.
The parameters of the manifold come after the layout: K=, rf=, mr= (mdot-ratio), pr= (p0-ratio), dh0= (J/kg); the
defaults 0, 1, 1, 1, 0 are a loss-free plenum. The initial state is the uniform one with the pressure of the first block
raised by 2 %: the run has to find the steady state through the manifold.

Checks (exit 0 when they hold):
  same DIR_A DIR_B TOL         p, u, h of the two outputs DIR/field.tec, in x order, equal within TOL (relative; 0 = bit
                               for bit); a reversed second block is put back in x order
  carried LOG DIR TOL [par]    the last manifold line of LOG: mdot in = mr mdot out, h0 in = h0 out - dh0 (to TOL of h0),
                               p0 in = pr (p0 out - K q_out), q_out = rho u^2/2 of the last cell of block 1 of DIR
                               (uniform there), each to TOL (relative)
  split DIR MR                 the mass flux rho u of block 2 is MR times that of block 1, block 1 at g (relative 1e-9)
  drop DIR DH0                 the total enthalpy h + u^2/2 of block 2 is that of block 1 minus DH0 (1e-9 of h0 + 1e-6)

Usage (from a run directory holding INPUT/):
  python3 -B make_case.py case LAYOUT [K=..] [rf=..] [mr=..] [pr=..] [dh0=..]   INPUT/ic.tec and INPUT/bc.txt
"""
import os
import re
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from ares_io import Table  # noqa: E402

N = 40                                   # cells along x
H = 1.0 / N                              # cell size, m
T_IN, G_IN, P_BACK = 345.0, 2000.0, 7.0e6
TABLE = os.path.join(HERE, '..', 'real-fluid-shock-tube', 'tables', 'thermo.dat')


def fresh(path):
    """Remove path first, so that a link left by another case is never written through."""
    if os.path.lexists(path):
        os.remove(path)
    return open(path, 'w')


def exact():
    """p, u, h of the uniform steady state (h at the inlet temperature and the back pressure)."""
    tab = Table(TABLE)
    h = tab.h_of(P_BACK, T_IN)
    return P_BACK, G_IN / tab('Density', P_BACK, h), h


def write_ic(blocks):
    """blocks: (i0, n, reversed) per block: first global cell, cells along x, written from high x to low x."""
    p, u, h = exact()
    with fresh('INPUT/ic.tec') as f:
        f.write(' VARIABLES = x y z p u v w h\n')
        for b, (i0, n, rev) in enumerate(blocks, start=1):
            f.write(' ZONE  T = B%d-RF, I=%d, J=2, K=2, DATAPACKING=BLOCK, '
                    'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n' % (b, n + 1))
            for c in range(3):
                for k in range(2):
                    for j in range(2):
                        for i in range(n + 1):
                            if rev:
                                x, y = (i0 + n - i) * H, (1 - j) * H
                            else:
                                x, y = (i0 + i) * H, j * H
                            f.write(' %.16E\n' % (x, y, k * H)[c])
            pb = p * (1.02 if b == 1 else 1.0)
            for value in (pb, u, 0.0, 0.0, h):
                f.write(' %.16E\n' % value * n)


def rec(b, i, f, t, prop=''):
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, 1, 1, f, t)
    if prop:
        s += prop + '\n'
    return s


def manifold(role, par):
    """The property line of a 501 record of group 1, as the boundary-condition builder of ATLAS writes it."""
    return '       1, %7d, %s' % (role, ''.join(' %.16E,' % par[k] for k in ('K', 'rf', 'mr', 'pr', 'dh0')))


def write_case(layout, par):
    inlet = '%r, %r, normal, normal, 1.0' % (T_IN, G_IN)
    outlet = '%.6E, 1.0' % P_BACK
    sides = lambda b, n: [rec(b, i, f, 300) for f in (3, 4, 5, 6) for i in range(1, n + 1)]
    n = N // 2
    if layout == 'one':
        write_ic([(0, N, False)])
        out = [rec(1, 1, 1, 404, inlet), rec(1, N, 2, 406, outlet)] + sides(1, N)
    elif layout == 'manifold':
        write_ic([(0, n, False), (n, n, False)])
        out = [rec(1, 1, 1, 404, inlet), rec(1, n, 2, 501, manifold(1, par)),
               rec(2, 1, 1, 501, manifold(2, par)), rec(2, n, 2, 406, outlet)] + sides(1, n) + sides(2, n)
    elif layout == 'reverse':
        write_ic([(0, n, False), (n, n, True)])
        out = [rec(1, 1, 1, 404, inlet), rec(1, n, 2, 501, manifold(1, par)),
               rec(2, 1, 1, 406, outlet), rec(2, n, 2, 501, manifold(2, par))] + sides(1, n) + sides(2, n)
    else:
        sys.exit('unknown layout %s' % layout)
    with fresh('INPUT/bc.txt') as f:
        f.write(''.join(out))


def read_zones(path):
    """The zones of an ASCII Tecplot solution of ARES (DATAPACKING=BLOCK, x y z nodal, the rest cell-centred, one cell in
    j and k here): a list of {name: [cell values]} with the cell centre x, in the order of the zones."""
    with open(path) as f:
        names = f.readline().split('=', 1)[1].replace('"', ' ').split()
        text = f.read()
    zones = []
    for part in re.split(r'^\s*ZONE', text, flags=re.M)[1:]:
        head, _, body = part.partition('\n')
        ni, nj, nk = (int(re.search(r'\b%s=\s*(\d+)' % d, head).group(1)) for d in ('I', 'J', 'K'))
        vals = []
        for line in body.split('\n'):
            s = line.strip()
            if not s:
                continue
            if s[0].isalpha():
                break
            vals += [float(t) for t in s.split()]
        nn, nc = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)
        x = vals[:ni]
        z = {'x': [0.5 * (x[m] + x[m + 1]) for m in range(ni - 1)]}
        for m, name in enumerate(names[3:]):
            z[name] = vals[3 * nn + m * nc: 3 * nn + (m + 1) * nc]
        zones.append(z)
    return zones


def in_x_order(directory):
    zones = read_zones(os.path.join(directory, 'field.tec'))
    out = {}
    for z in zones:
        rev = len(z['x']) > 1 and z['x'][0] > z['x'][-1]
        for k, v in z.items():
            out.setdefault(k, []).extend(v[::-1] if rev else list(v))
    return out


def same(a, b, tol):
    fa, fb = in_x_order(a), in_x_order(b)
    worst = 0.0
    for k in ('p', 'u', 'h'):
        if len(fa[k]) != len(fb[k]):
            print('    %s: %d against %d cells' % (k, len(fa[k]), len(fb[k])))
            return False
        for x, y in zip(fa[k], fb[k]):
            if struct.pack('<d', x) != struct.pack('<d', y):
                worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1e-300))
    if worst == 0.0:
        print('    bit for bit the same (p, u, h, %d cells)' % len(fa['p']))
        return True
    print('    largest relative difference of p, u, h %.3e (required at most %.1e)' % (worst, tol))
    return worst <= tol


def last_manifold_line(log):
    line = None
    for ln in open(log, errors='replace'):
        if ln.lstrip().startswith('manifold 1,'):
            line = ln
    if line is None:
        print('    no manifold line in %s' % log)
        sys.exit(1)
    nums = [float(v) for v in re.findall(r'[-+]?\d+\.\d+E[-+]\d+', line)]
    if len(nums) < 7:
        print('    cannot read the manifold line: %s' % line.strip())
        sys.exit(1)
    return dict(zip(('mdot_out', 'mdot_in', 'h0_out', 'h0_in', 'p0_out', 'p0_in', 'power'), nums))


def carried(log, directory, tol, par):
    v = last_manifold_line(log)
    b1 = read_zones(os.path.join(directory, 'field.tec'))[0]
    q = 0.5 * b1['rho'][-1] * b1['u'][-1] ** 2   # dynamic pressure of the outlet cells (uniform across the duct)
    bad = 0
    d = abs(v['mdot_in'] - par['mr'] * v['mdot_out']) / abs(v['mdot_out'])
    print('    mdot out %.10e in %.10e: in = mr out off by %.1e' % (v['mdot_out'], v['mdot_in'], d))
    bad += d > tol
    d = abs(v['h0_in'] - (v['h0_out'] - par['dh0'])) / max(abs(v['h0_out']), 1.0)
    print('    h0 out %.10e in %.10e: in = out - dh0 off by %.1e' % (v['h0_out'], v['h0_in'], d))
    bad += d > tol
    d = abs(v['p0_in'] - par['pr'] * (v['p0_out'] - par['K'] * q)) / abs(v['p0_in'])
    print('    p0 out %.10e in %.10e: in = pr (out - K q_out) off by %.1e' % (v['p0_out'], v['p0_in'], d))
    bad += d > tol
    return bad == 0


def split(directory, mr):
    zones = read_zones(os.path.join(directory, 'field.tec'))
    g = [sum(r * u for r, u in zip(z['rho'], z['u'])) / len(z['u']) for z in zones]
    d1, d2 = abs(g[0] - G_IN) / G_IN, abs(abs(g[1]) - mr * G_IN) / G_IN
    print('    mass flux of block 1 %.10e (g %r), of block 2 %.10e (mr g %r): off by %.1e, %.1e'
          % (g[0], G_IN, g[1], mr * G_IN, d1, d2))
    return max(d1, d2) <= 1e-9


def drop(directory, dh0):
    zones = read_zones(os.path.join(directory, 'field.tec'))
    h0 = [sum(h + 0.5 * u * u for h, u in zip(z['h'], z['u'])) / len(z['u']) for z in zones]
    d = abs((h0[0] - h0[1]) - dh0)
    print('    h0 of block 1 %.8f, of block 2 %.8f J/kg: the drop %.6f against dh0 %.6f, off by %.1e J/kg'
          % (h0[0], h0[1], h0[0] - h0[1], dh0, d))
    return d <= 1e-9 * abs(h0[0]) + 1e-6


def parameters(words):
    par = {'K': 0.0, 'rf': 1.0, 'mr': 1.0, 'pr': 1.0, 'dh0': 0.0}
    for w in words:
        if '=' in w:
            k, v = w.split('=', 1)
            if k not in par:
                sys.exit('unknown parameter %s' % k)
            par[k] = float(v)
    return par


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    what, args = sys.argv[1], sys.argv[2:]
    if what == 'case':
        write_case(args[0], parameters(args[1:]))
        return
    if what == 'same':
        ok = same(args[0], args[1], float(args[2]))
    elif what == 'carried':
        ok = carried(args[0], args[1], float(args[2]), parameters(args[3:]))
    elif what == 'split':
        ok = split(args[0], float(args[1]))
    elif what == 'drop':
        ok = drop(args[0], float(args[1]))
    else:
        sys.exit(__doc__)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
