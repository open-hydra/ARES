"""Figure of docs/vv/real-fluid-shock-tube.md: p, u, rho, T of the ARES run against the exact solution (exact.txt).

  python3 plot_vv.py [RUN]     RUN = the run folder (../runs/muscl by default); needs matplotlib
Writes ../../../docs/vv/img/shock_tube.svg.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
from ares_io import read_field  # noqa: E402


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    run = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'runs', 'muscl')
    (ni, _, _), nodes, c = read_field(os.path.join(run, 'OUTPUT', 'field.tec'))
    x = nodes['x'][:ni]
    xc = [0.5 * (x[m] + x[m + 1]) for m in range(ni - 1)]
    ex = [list(map(float, line.split())) for line in open(os.path.join(HERE, 'exact.txt'))
          if not line.startswith('#')][2:]
    grey = '#8a8a8a'
    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'legend.fontsize': 10,
                         'svg.fonttype': 'none', 'font.family': 'sans-serif',
                         'text.color': grey, 'axes.labelcolor': grey, 'axes.edgecolor': grey,
                         'xtick.color': grey, 'ytick.color': grey, 'grid.color': grey,
                         'savefig.transparent': True})
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.2), sharex=True)
    panels = (('p', 1, 1e-6, '$p$ [MPa]'), ('u', 2, 1.0, '$u$ [m/s]'),
              ('rho', 3, 1.0, r'$\rho$ [kg/m$^3$]'), ('T', 4, 1.0, '$T$ [K]'))
    for ax, (var, col, scale, label) in zip(axes.flat, panels):
        ax.plot([r[0] for r in ex], [r[col] * scale for r in ex], '-', color='#3b8bd9', lw=1.6, label='exact (CoolProp)')
        ax.plot(xc, [v * scale for v in c[var]], 'o', ms=2.4, mfc='none', mew=0.8, color='#e8710a', label='ARES')
        ax.set_ylabel(label)
        ax.set_xlim(0.25, 0.8)
        ax.grid(alpha=0.25)
    for ax in axes[1]:
        ax.set_xlabel('$x$ [m]')
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    out = os.path.join(HERE, '..', '..', '..', 'docs', 'vv', 'img', 'shock_tube.svg')
    fig.savefig(out)
    print('wrote', os.path.normpath(out))


if __name__ == '__main__':
    main()
