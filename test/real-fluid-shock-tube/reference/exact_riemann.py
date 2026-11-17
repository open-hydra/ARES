"""Exact solution of the real-fluid shock tube of this test, from the CoolProp reference equation of state of CO2.

Left and right states at rest, (p, T) below; the solution is a left rarefaction, a contact and a right shock:
  rarefaction  u* = u_L + int_{p*}^{p_L} dp / (rho a), along the isentrope s = s_L;
  shock        the Hugoniot h - h_R = (p - p_R)(1/rho + 1/rho_R)/2, mass flux m^2 = (p - p_R)/(1/rho_R - 1/rho),
               u* = u_R + (p - p_R)/m, shock speed S = u_R + m/rho_R;
p* is the root of u*_L(p) = u*_R(p). The ideal-gas solution with the gamma = cp/cv and the density of each initial state
is printed too, to show how far the non-ideal solution lies from it.

Writes exact.txt (read by ../make_case.py at check time): the star state and the wave speeds, and the profile of
p, u, rho, T at the time of the check on 2001 points, for figures. Needs CoolProp (conda env ct-env):
  python3 exact_riemann.py
"""
import math
import CoolProp.CoolProp as CP

FLUID = 'HEOS::CO2'
PL, TL = 10.0e6, 360.0
PR, TR = 5.0e6, 330.0
X0, T_END = 0.5, 6.0e-4          # diaphragm and time of the check (the waves stay inside [0, 1])
N_INT = 400


def prop(out, a, va, b, vb):
    return CP.PropsSI(out, a, va, b, vb, FLUID)


rhoL, rhoR = prop('D', 'P', PL, 'T', TL), prop('D', 'P', PR, 'T', TR)
PC = CP.PropsSI('pcrit', FLUID)
sL, hR = prop('S', 'P', PL, 'T', TL), prop('H', 'P', PR, 'T', TR)


def isentrope(p):
    """rho, a, T on the isentrope s = s_L: T by bisection on s(p, T) (the (p, s) flash of CoolProp fails near the
    pseudo-critical line), between saturation, below p_c, and 600 K."""
    lo = prop('T', 'P', p, 'Q', 1) + 1e-6 if p < PC else 250.0
    hi = 600.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if prop('S', 'P', p, 'T', mid) < sL:
            lo = mid
        else:
            hi = mid
    T = 0.5 * (lo + hi)
    return prop('D', 'P', p, 'T', T), prop('A', 'P', p, 'T', T), T


def u_rarefaction(p):
    """u behind the left rarefaction that ends at p (Simpson on N_INT intervals)."""
    n = N_INT
    dp = (PL - p) / n
    f = lambda q: 1.0 / (isentrope(q)[0] * isentrope(q)[1])
    s = f(p) + f(PL) + sum((4 if k % 2 else 2) * f(p + k * dp) for k in range(1, n))
    return s * dp / 3.0


def hugoniot(p):
    """T, rho behind the right shock to p (bisection on T)."""
    def res(T):
        rho = prop('D', 'P', p, 'T', T)
        return prop('H', 'P', p, 'T', T) - hR - 0.5 * (p - PR) * (1.0 / rho + 1.0 / rhoR)
    lo, hi = TR, TR + 200.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if res(mid) > 0:
            hi = mid
        else:
            lo = mid
    T = 0.5 * (lo + hi)
    return T, prop('D', 'P', p, 'T', T)


def u_shock(p):
    T, rho = hugoniot(p)
    m = math.sqrt((p - PR) / (1.0 / rhoR - 1.0 / rho))
    return (p - PR) / m, m, T, rho


def solve():
    lo, hi = PR, PL
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if u_rarefaction(mid) - u_shock(mid)[0] > 0:
            lo = mid
        else:
            hi = mid
    ps = 0.5 * (lo + hi)
    us, m, TsR, rhosR = u_shock(ps)
    rhosL, asL, TsL = isentrope(ps)
    return dict(p=ps, u=us, rhoL=rhosL, TL=TsL, rhoR=rhosR, TR=TsR,
                head=-prop('A', 'P', PL, 'T', TL), tail=us - asL, shock=m / rhoR)


def ideal(gL, gR):
    """p*, u* of the ideal-gas Riemann problem with gamma gL on the left and gR on the right (same p, rho)."""
    aL, aR = math.sqrt(gL * PL / rhoL), math.sqrt(gR * PR / rhoR)
    def f(p):
        uL = 2 * aL / (gL - 1) * (1 - (p / PL) ** ((gL - 1) / (2 * gL)))
        A, B = 2 / ((gR + 1) * rhoR), (gR - 1) / (gR + 1) * PR
        return uL - (p - PR) * math.sqrt(A / (p + B))
    lo, hi = PR, PL
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    p = 0.5 * (lo + hi)
    return p, 2 * aL / (gL - 1) * (1 - (p / PL) ** ((gL - 1) / (2 * gL)))


def main():
    s = solve()
    gL = prop('Cpmass', 'P', PL, 'T', TL) / prop('Cvmass', 'P', PL, 'T', TL)
    gR = prop('Cpmass', 'P', PR, 'T', TR) / prop('Cvmass', 'P', PR, 'T', TR)
    pi, ui = ideal(gL, gR)
    print('left  p %.6e T %.4f rho %.6f Z %.4f' % (PL, TL, rhoL, prop('Z', 'P', PL, 'T', TL)))
    print('right p %.6e T %.4f rho %.6f Z %.4f' % (PR, TR, rhoR, prop('Z', 'P', PR, 'T', TR)))
    for k in ('p', 'u', 'rhoL', 'TL', 'rhoR', 'TR', 'head', 'tail', 'shock'):
        print('%-6s %.10e' % (k, s[k]))
    print('ideal gas (gamma %.3f / %.3f): p* %.6e (%+.1f%%), u* %.4f (%+.1f%%)'
          % (gL, gR, pi, 100 * (pi / s['p'] - 1), ui, 100 * (ui / s['u'] - 1)))
    with open('exact.txt', 'w') as f:
        f.write('# Exact solution of the CO2 shock tube (exact_riemann.py, CoolProp %s)\n' % CP.get_global_param_string('version'))
        f.write('# pL TL pR TR x0 t\n%.10e %.10e %.10e %.10e %.10e %.10e\n' % (PL, TL, PR, TR, X0, T_END))
        f.write('# p* u* rho*L T*L rho*R T*R head tail shock\n')
        f.write(' '.join('%.10e' % s[k] for k in ('p', 'u', 'rhoL', 'TL', 'rhoR', 'TR', 'head', 'tail', 'shock')) + '\n')
        f.write('# x p u rho T\n')
        # the fan tabulated from p* to p_L: (u - a, p, u, rho, T), u - a increasing with u (decreasing p)
        fan, u, n = [], 0.0, 400
        ps = [PL - (PL - s['p']) * m / n for m in range(n + 1)]
        for m, p in enumerate(ps):
            rho, a, T = isentrope(p)
            if m:
                rho0, a0, _ = isentrope(ps[m - 1])
                u += 0.5 * (ps[m - 1] - p) * (1 / (rho * a) + 1 / (rho0 * a0))
            fan.append((u - a, p, u, rho, T))
        for n in range(2001):
            x = n / 2000.0
            xi = (x - X0) / T_END
            if xi < s['head']:
                row = (PL, 0.0, rhoL, TL)
            elif xi < s['tail']:          # inside the fan: u - a = xi, interpolated in the tabulated isentrope
                k = max(m for m in range(len(fan)) if fan[m][0] <= xi)
                k = min(k, len(fan) - 2)
                w = (xi - fan[k][0]) / (fan[k + 1][0] - fan[k][0])
                row = tuple((1 - w) * fan[k][c] + w * fan[k + 1][c] for c in range(1, 5))
            elif xi < s['u']:
                row = (s['p'], s['u'], s['rhoL'], s['TL'])
            elif xi < s['shock']:
                row = (s['p'], s['u'], s['rhoR'], s['TR'])
            else:
                row = (PR, 0.0, rhoR, TR)
            f.write('%.6f %.10e %.10e %.10e %.10e\n' % ((x,) + row))


if __name__ == '__main__':
    main()
