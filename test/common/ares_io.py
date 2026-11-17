"""Readers shared by the CI tests (Python standard library only): the real-fluid table of ATLAS GPB and the ASCII
Tecplot solution of ARES."""


class Table:
    """thermo.dat of ATLAS GPB: Tecplot DATAPACKING=BLOCK on a regular (p, h) grid, p fastest. Interpolated
    bilinearly in (p, h), as ARES reads it."""

    def __init__(self, path):
        with open(path) as f:
            f.readline()
            names = [s.strip().strip('"') for s in f.readline().split('=', 1)[1].split(',')]
            zone = f.readline()
            ni = int(zone.split('I=')[1].split(',')[0])
            nj = int(zone.split('J=')[1].split(',')[0])
            vals = [float(t) for t in f.read().split()]
        n = ni * nj
        self.var = {name: vals[m * n:(m + 1) * n] for m, name in enumerate(names)}
        self.ni, self.nj = ni, nj
        self.p = self.var['Pressure'][:ni]
        self.h = self.var['Enthalpy'][::ni]

    def __call__(self, name, p, h):
        """The variable name ('Density', 'Temperature', 'SoundSpeed', ...) at (p, h)."""
        def locate(axis, x):
            d = (axis[-1] - axis[0]) / (len(axis) - 1)
            m = min(max(int((x - axis[0]) / d), 0), len(axis) - 2)
            return m, (x - axis[m]) / (axis[m + 1] - axis[m])
        i, s = locate(self.p, p)
        j, t = locate(self.h, h)
        v = self.var[name]
        a, b = v[j * self.ni + i], v[j * self.ni + i + 1]
        c, d = v[(j + 1) * self.ni + i], v[(j + 1) * self.ni + i + 1]
        return (1 - t) * ((1 - s) * a + s * b) + t * ((1 - s) * c + s * d)

    def h_of(self, p, T):
        """The enthalpy at (p, T), by bisection of the interpolated T(p, h)."""
        lo, hi = self.h[0], self.h[-1]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self('Temperature', p, mid) < T:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def read_field(path):
    """The first zone of an ASCII Tecplot solution of ARES (DATAPACKING=BLOCK, x y z nodal, the rest cell-centred):
    (dims, nodes, cells), dims = (I, J, K) node counts, nodes = {'x': [...], 'y': ..., 'z': ...}, cells = {name: [...]},
    i fastest."""
    with open(path) as f:
        names = f.readline().split('=', 1)[1].replace('"', ' ').split()
        zone = f.readline()
        vals = []
        for line in f:
            if not line.strip() or line.lstrip()[0].isalpha():
                break
            vals += [float(t) for t in line.split()]
    ni, nj, nk = (int(zone.split(' %s=' % d)[1].split(',')[0]) for d in ('I', 'J', 'K'))
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    nodes = {name: vals[m * nn:(m + 1) * nn] for m, name in enumerate(names[:3])}
    cells = {name: vals[3 * nn + m * nc:3 * nn + (m + 1) * nc] for m, name in enumerate(names[3:])}
    return (ni, nj, nk), nodes, cells


def solution_time(path):
    """The SOLUTIONTIME of the first zone of an ASCII Tecplot solution of ARES."""
    with open(path) as f:
        f.readline()
        return float(f.readline().split('SOLUTIONTIME=')[1].split(',')[0])
