#!/usr/bin/env python3
"""Centreline pressure ratio and Mach number of the ARES solution against test M1.5H (Spinelli et al. 2018) and
the SU2 solution of the NICFD_nozzle tutorial.

usage: python3 reference/verify.py [OUTPUT/field.tec]
       prints the comparison at the pressure taps and writes reference/nozzle_p.svg, reference/nozzle_mach.svg

The experimental x (origin at the nozzle inlet section) is the x of the mesh. The SU2 tutorial mesh ends at
x = 0.1229 m: the taps and Mach points downstream of it are outside the domain."""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
P0_EXP = 904387.9          # total pressure of test M1.5H [Pa]


def read_field(path):
    """Tecplot ASCII block field of ARES: nodal x, y, z, cell-centred variables, arrays indexed [k, j, i]."""
    with open(path) as fh:
        names = fh.readline().split('=')[1].split()
        z = fh.readline()
        ni, nj, nk = (int(z.split(f'{c}=')[1].split(',')[0]) for c in 'IJK')
        vals = np.array([v for l in fh if 'AUXDATA' not in l for v in l.split()], dtype=float)
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    out, pos = {}, 0
    for n in names:
        size = nn if n in ('x', 'y', 'z') else nc
        shape = (nk, nj, ni) if size == nn else (nk - 1, nj - 1, ni - 1)
        out[n] = vals[pos:pos + size].reshape(shape)
        pos += size
    return out


def read_experiment(path):
    """Pressure taps (x, P, U, ...) and Mach points (x, M, U) of experimentald_data.txt."""
    pres, mach, block = [], [], None
    for l in open(path).read().splitlines():
        if 'Pressure and temperature data' in l: block = pres
        elif 'Mach data' in l: block = mach
        elif block is not None:
            try: row = [float(v) for v in l.split()]
            except ValueError: continue
            if row: block.append(row)
    return np.array(pres), np.array(mach)


def centreline(f):
    """x, p, Mach of the cells on the axis (j = 1)."""
    xc = 0.25 * (f['x'][0, 0, :-1] + f['x'][0, 0, 1:] + f['x'][0, 1, :-1] + f['x'][0, 1, 1:])
    return xc, f['p'][0, 0], np.abs(f['u'][0, 0]) / f['sound'][0, 0]


def plot(xa, ya, exp, xlabel, ylabel, path, su2=None, xmax=0.125):
    """Figure readable on the light and the dark theme of the docs: transparent background, mid-grey text and axes,
    coloured lines, hollow markers. x is cut at xmax (the end of the CFD domain)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    grey = '#8a8a8a'
    plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'legend.fontsize': 11,
                         'svg.fonttype': 'none', 'font.family': 'sans-serif',
                         'text.color': grey, 'axes.labelcolor': grey, 'axes.edgecolor': grey,
                         'xtick.color': grey, 'ytick.color': grey, 'grid.color': grey,
                         'savefig.transparent': True})
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(xa, ya, '-', color='#e8710a', lw=2.0, label='ARES')
    if su2 is not None:
        ax.plot(*su2, '--', color='#3b8bd9', lw=1.8, label='SU2')
    keep = exp[0] <= xmax
    ax.errorbar(exp[0][keep], exp[1][keep], yerr=exp[2][keep], fmt='o', ms=5, mfc='none', mew=1.3,
                color='#d63c3c', capsize=2, label='Experiment (M1.5H)')
    ax.set_xlim(right=xmax)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.grid(alpha=0.25); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print(f'wrote {os.path.relpath(path, CASE)}')


def plot_field(f, path):
    """Mach number field, mirrored about the axis to show the whole nozzle (cell values on the cell corners)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x = f['x'][0]; y = f['y'][0]                       # nodes (j, i)
    M = np.hypot(f['u'][0], f['v'][0]) / f['sound'][0]   # cells (j, i)
    xc = 0.25 * (x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:])
    yc = 0.25 * (y[:-1, :-1] + y[1:, :-1] + y[:-1, 1:] + y[1:, 1:])
    # close the field on the axis (y = 0) and on the wall: those rows take the values of their cells
    yw = np.interp(xc[-1], x[-1], y[-1])
    X = np.vstack([xc[:1], xc, xc[-1:]])
    Y = np.vstack([np.zeros_like(yc[:1]), yc, yw[None, :]])
    MM = np.vstack([M[:1], M, M[-1:]])
    levels = np.linspace(0.0, 1.6, 33)
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for sgn in (1, -1):
        cs = ax.contourf(X, sgn * Y, MM, levels=levels, cmap='viridis', extend='max')
        cs.set_edgecolor('face')                          # no seams between the bands in the SVG
        cs.set_linewidth(0.3)
        ax.contour(X, sgn * Y, MM, levels=[1.0], colors='white', linewidths=1.0, linestyles='--')
        ax.plot(x[-1], sgn * y[-1], color='#8a8a8a', lw=1.2)
    for side in ('top', 'right', 'bottom', 'left'):
        ax.spines[side].set_visible(False)
    ax.set_xlim(x.min(), x.max())
    ax.set_aspect('equal'); ax.set_xlabel(r'$x$ [m]'); ax.set_ylabel(r'$y$ [m]')
    cb = fig.colorbar(cs, ax=ax, shrink=0.8, pad=0.02, ticks=np.arange(0, 1.61, 0.4))
    cb.set_label('Mach number'); cb.outline.set_edgecolor('#8a8a8a')
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print(f'wrote {os.path.relpath(path, CASE)}')


def main():
    field = read_field(sys.argv[1] if len(sys.argv) > 1 else os.path.join(CASE, 'OUTPUT', 'field.tec'))
    xa, pa, Ma = centreline(field)
    pres, mach = read_experiment(os.path.join(HERE, 'experimentald_data.txt'))
    xs, ps = np.loadtxt(os.path.join(HERE, 'SU2_data.txt')).T

    print('   x [m]   P/PT exp   SU2      ARES     ARES-exp   U(exp)')
    for x, P, U in pres[:, :3]:
        if x > xa[-1]:
            print(f'  {x:.4f}   {P / P0_EXP:.4f}     outside the domain'); continue
        a, s = np.interp(x, xa, pa) / P0_EXP, np.interp(x, xs, ps) / P0_EXP
        print(f'  {x:.4f}   {P / P0_EXP:.4f}   {s:.4f}   {a:.4f}   {a - P / P0_EXP:+.4f}    {U / P0_EXP:.4f}')
    d = np.interp(xs, xa, pa) / P0_EXP - ps / P0_EXP
    print(f'ARES - SU2 along the centreline: max |dP/PT| = {np.abs(d).max():.4f}, rms {np.sqrt((d ** 2).mean()):.4f}')
    inside = mach[mach[:, 0] <= xa[-1]]
    dm = np.interp(inside[:, 0], xa, Ma) - inside[:, 1]
    print(f'Mach: {len(inside)} of {len(mach)} points inside, ARES - exp mean {dm.mean():+.4f}, rms {np.sqrt((dm ** 2).mean()):.4f}')
    print(f'exit: Mach {Ma[-1]:.3f}, P/PT {pa[-1] / P0_EXP:.4f}')

    try:
        plot(xa, pa / P0_EXP, (pres[:, 0], pres[:, 1] / P0_EXP, pres[:, 2] / P0_EXP), r'$x$ [m]', r'$P/P_T$',
             os.path.join(HERE, 'nozzle_p.svg'), su2=(xs, ps / P0_EXP))
        plot(xa, Ma, mach.T, r'$x$ [m]', 'Mach number', os.path.join(HERE, 'nozzle_mach.svg'))
        plot_field(field, os.path.join(HERE, 'nozzle_mach_field.svg'))
    except ImportError:
        print('matplotlib not available: no figures')


if __name__ == '__main__':
    main()
