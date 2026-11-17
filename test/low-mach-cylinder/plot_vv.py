"""Figure of docs/vv/low-mach-cylinder.md: Cp on the cylinder (first row of cells) of the three runs, and potential flow.

  python3 plot_vv.py [RUNS]    RUNS = the folder of the runs (runs/ by default); needs matplotlib
Writes ../../docs/vv/img/cylinder_cp.svg and prints the numbers quoted on the page.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc  # noqa: E402


def wall(path, mach):
    """theta from the front stagnation point [deg] and Cp of the first row of cells, sorted by theta, upper half."""
    pts = []
    for cp, cx, cy, j in mc.cp_field(path, mach):
        if j == 0 and cy >= 0.0:
            pts.append((180.0 - math.degrees(math.atan2(cy, cx)), cp, math.hypot(cx, cy)))
    return sorted(pts)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
    f = lambda tag: os.path.join(runs, tag, 'OUTPUT', 'field.tec')
    m2, m3, pr = wall(f('prec-m2'), 0.01), wall(f('prec-m3'), 0.001), wall(f('prim-m2'), 0.01)
    th = [p[0] for p in m2]
    pot = [2.0 * (mc.R / r) ** 2 * math.cos(2.0 * math.radians(180.0 - t)) - (mc.R / r) ** 4 for t, _, r in m2]
    grey = '#8a8a8a'
    plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'legend.fontsize': 10,
                         'svg.fonttype': 'none', 'font.family': 'sans-serif',
                         'text.color': grey, 'axes.labelcolor': grey, 'axes.edgecolor': grey,
                         'xtick.color': grey, 'ytick.color': grey, 'grid.color': grey,
                         'savefig.transparent': True})
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(th, pot, '-', color=grey, lw=1.4, label='potential flow (cell centres)')
    ax.plot(th, [p[1] for p in m2], '-', color='#e8710a', lw=2.0, label='prec, HLLC Prec, M = 0.01')
    ax.plot([p[0] for p in m3], [p[1] for p in m3], 'o', ms=4, mfc='none', color='#3b8bd9',
            label='prec, HLLC Prec, M = 0.001')
    ax.plot([p[0] for p in pr], [p[1] for p in pr], '--', color='#d63c3c', lw=1.6, label='prim, HLLC, M = 0.01')
    ax.set_xlabel(r'$\theta$ from the front stagnation point [deg]')
    ax.set_ylabel('$C_p$')
    ax.set_xlim(0, 180)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc='lower left')
    fig.tight_layout()
    out = os.path.join(HERE, '..', '..', 'docs', 'vv', 'img', 'cylinder_cp.svg')
    fig.savefig(out)
    print('wrote', os.path.normpath(out))
    for name, d in (('prec M=0.01', m2), ('prec M=0.001', m3), ('prim M=0.01', pr)):
        stag = d[0][1]
        cmin = min(p[1] for p in d)
        print('%-13s Cp stagnation %.3f, minimum %.3f; range of Cp %.2f to %.2f'
              % (name, stag, cmin, min(p[1] for p in d), max(p[1] for p in d)))
    print('potential     Cp first cell %.3f, minimum %.3f' % (pot[0], min(pot)))


if __name__ == '__main__':
    main()
