# -*- coding: utf-8 -*-
"""
Shared helpers for the NASA-TMR turbulent flat-plate validation scripts
(validate_sst.py, validate_wilcox2006.py, validate_sgglrr.py).

All three cases are the NASA Turbulence Modeling Resource 2D zero-pressure
gradient flat plate, defined at Re = 5e6, M = 0.2. The ARES runs use the
real-fluid air table (mu ~ 2.9e-5 at T~555 K), so unless the inflow mass flux
is rescaled the run sits at a LOWER Reynolds than the nominal 5e6:

  * cf(x)  -> compared in a Reynolds-consistent way: the tabulated NASA cf is
             re-expressed at the run's own Re_x (cf_ref(x) = cf_NASA(Re_x(x))),
             exactly as the laminar script evaluates Blasius at the run's Re_x.
  * wall-normal profiles at the NASA station x~0.97 (mut/mu, k+, omega+) are
             Reynolds-sensitive and DO NOT collapse if the run is not at
             Re = 5e6. The scripts print a clear warning in that case.

cf definition (from each PLOT.lay):
    Tau_x = mu * |dU/dy|_wall ,   c_f = Tau_x / (0.5 * rho_inf * U_inf^2)
Here the solver's own wall shear tau (column 25 of 1d.dat) is used.
"""

import os
import sys
import numpy as np

# Column indices in 1d.dat (0-based)
# 0:x 6:U 11:rho 16:mi(=mu) 25:tau(total)
C_X, C_U, C_RHO, C_MU, C_TAU = 0, 6, 11, 16, 25

X_LO      = 0.0125    # cf error-metric lower bound (matches MOSE verify.py)
RE_REF    = 5.0e6     # NASA-TMR reference Reynolds per unit length
X_STATION = 0.97      # NASA profile station


# --------------------------------------------------------------------------- #
def read_wall_tec(path):
    """Read the wall Tecplot BLOCK file (x,y,z nodal; y+,tauX,tauY,tauZ,pw,Tw,qw
       cell-centred; J=1 surface zone). Returns x cell-centres and a dict of
       wall arrays. tauX is the solver's wall shear stress."""
    with open(path) as fh:
        fh.readline()                       # VARIABLES = ...
        zone = fh.readline()
        body = np.array(fh.read().split(), dtype=float)

    def _dim(tag):
        s = zone.upper().split(tag + "=")[1]
        return int("".join(c for c in s.split(",")[0] if c.isdigit()))
    NI, NJ, NK = _dim("I"), _dim("J"), _dim("K")
    nnode = NI * NJ * NK
    ncell = (NI - 1) * (NK - 1)             # J = 1 surface
    off = 0
    xN = body[off:off + nnode]; off += nnode
    off += 2 * nnode                        # skip y, z nodal
    names = ["yplus", "tauX", "tauY", "tauZ", "pw", "Tw", "qw"]
    cv = {nm: body[off + k*ncell: off + (k+1)*ncell] for k, nm in enumerate(names)}
    xN = xN.reshape(NK, NJ, NI)
    xc = 0.5 * (xN[0, 0, :-1] + xN[0, 0, 1:])
    return xc, cv


def field_freestream(xN, yN, cv):
    """Free-stream U, rho, mu from the field edge (row farthest from the wall),
       median over the plate columns x>0. Consistent with wall.tec (same run),
       unlike a possibly stale 1d.dat."""
    u, rho, mu = cv["u"], cv["rho"], cv["mil"]
    fs = -1 if u[0, :].mean() < u[-1, :].mean() else 0   # freestream j-edge
    xc = 0.25 * (xN[0, :-1, :-1] + xN[0, :-1, 1:] + xN[0, 1:, :-1] + xN[0, 1:, 1:])
    cols = xc[0, :] > 0.0
    return (float(np.median(u[fs, cols])),
            float(np.median(rho[fs, cols])),
            float(np.median(mu[fs, cols])))


def read_tecplot_zones(path):
    """Parse a (possibly multi-zone) Tecplot point file. Returns an ordered
       dict {zone_title: ndarray}. Case-insensitive on ZONE/VARIABLES; handles
       commas inside quoted zone titles."""
    zones, cur, rows = {}, None, []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            low = s.lower()
            if low.startswith("zone"):
                if cur is not None:
                    zones[cur] = np.array(rows)
                if "t=" in low:
                    title = s[low.index("t=") + 2:].strip().strip(",").strip().strip('"')
                else:
                    title = "zone%d" % len(zones)
                cur, rows = title, []
            elif s and not s.startswith("#") and "variable" not in low:
                try:
                    rows.append([float(v) for v in s.replace(",", " ").split()])
                except ValueError:
                    pass
        if cur is not None:
            zones[cur] = np.array(rows)
    return zones


def pick_zone(zones, prefer=("CFL3D", "FUN3D")):
    """Pick a reference zone: first whose title contains a preferred tag,
       else the first zone."""
    for tag in prefer:
        for name in zones:
            if tag.lower() in name.lower():
                return name, zones[name]
    name = next(iter(zones))
    return name, zones[name]


def read_field_tec(path, cell_vars):
    """Read a Tecplot BLOCK field.tec with 3 nodal coords (x,y,z) followed by
       `cell_vars` cell-centred variables. Returns NI,NJ,NK, nodal xN,yN
       (shape K,J,I) and a dict of cell arrays (shape K-1,J-1,I-1)."""
    with open(path) as fh:
        fh.readline()                       # VARIABLES = ...
        zone = fh.readline()
        body = np.array(fh.read().split(), dtype=float)

    def _dim(tag):
        s = zone.upper().split(tag + "=")[1]
        return int("".join(c for c in s.split(",")[0] if c.isdigit()))
    NI, NJ, NK = _dim("I"), _dim("J"), _dim("K")
    nnode = NI * NJ * NK
    ncell = (NI - 1) * (NJ - 1) * (NK - 1)

    off = 0
    xN = body[off:off + nnode]; off += nnode
    yN = body[off:off + nnode]; off += nnode
    off += nnode                            # z (unused)
    cv = {}
    for nm in cell_vars:
        cv[nm] = body[off:off + ncell].reshape(NK - 1, NJ - 1, NI - 1)[0]
        off += ncell
    xN = xN.reshape(NK, NJ, NI)
    yN = yN.reshape(NK, NJ, NI)
    return NI, NJ, NK, xN, yN, cv


def wall_profile(xN, yN, cv, varlist, x_station):
    """Wall-normal profile at exactly x_station, linearly interpolated in x
       between the two bracketing cell columns (the x-mesh is coarse downstream,
       so nearest-cell would land ~0.04 m off). Ordered from the wall outward;
       the wall is the j-end with the lowest streamwise velocity 'u'.
       Returns x_sta(=x_station), wall-distance yn, and {var: array}."""
    xc = 0.25 * (xN[0, :-1, :-1] + xN[0, :-1, 1:] + xN[0, 1:, :-1] + xN[0, 1:, 1:])
    yc = 0.25 * (yN[0, :-1, :-1] + yN[0, :-1, 1:] + yN[0, 1:, :-1] + yN[0, 1:, 1:])
    xrow = xc[0, :]
    if x_station <= xrow[0]:
        i, w = 0, 0.0
    elif x_station >= xrow[-1]:
        i, w = len(xrow) - 2, 1.0
    else:
        i = int(np.searchsorted(xrow, x_station)) - 1
        w = (x_station - xrow[i]) / (xrow[i + 1] - xrow[i])

    def lerp(a):                                   # interpolate a column in x
        return (1.0 - w) * a[:, i] + w * a[:, i + 1]

    y_col = lerp(yc)
    y_nod = (1.0 - w) * yN[0, :, i] + w * yN[0, :, i + 1]
    u_col = lerp(cv["u"])
    y_wall = y_nod[0] if u_col[0] < u_col[-1] else y_nod[-1]
    yn = np.abs(y_col - y_wall)
    o = np.argsort(yn)
    out = {v: lerp(cv[v])[o] for v in varlist}
    return x_station, yn[o], out


# --------------------------------------------------------------------------- #
def reference_curves(zones):
    """Reference zones to plot, as MOSE's verify.py does: finest grid (level 1)
       or single-zone references, skipping coarse levels and the TAU solver."""
    out = []
    for name, arr in zones.items():
        low = name.lower()
        if low.startswith("tau"):                       # skip TAU solver
            continue
        if "level" in low and "level 1" not in low:     # only the finest grid
            continue
        out.append((name, arr))
    return out


def validate_cf(f_field, f_wall, f_cf, cell_vars, case_name, ares_color="C0"):
    """Skin friction cf(x), computed exactly like MOSE's verify.py:
         - wall shear tauX read from wall.tec (solver value, AT the wall),
         - free-stream U_inf, rho_inf from the field edge,
         - Cf = tauX / (0.5 * rho_inf * U_inf**2),
         - plotted directly vs x and compared to the NASA cf.dat zones.
       (Valid when the run is at the NASA conditions Re=5e6, M=0.2, i.e. with
       a constant mu=1.18587e-5 in transport.dat.) Returns free-stream state."""
    import matplotlib.pyplot as plt
    print("=" * 74)
    print(" 1) SKIN FRICTION  cf(x)   (wall.tec + field.tec, MOSE verify.py style)")
    print("=" * 74)

    _, _, _, xN, yN, cv = read_field_tec(f_field, cell_vars)
    U_inf, rho_inf, mu_inf = field_freestream(xN, yN, cv)
    Re_per_m = rho_inf * U_inf / mu_inf
    print(f"  Free-stream (field edge): U={U_inf:.3f} m/s  rho={rho_inf:.5f}  mu={mu_inf:.4e}")
    print(f"  Re/m = {Re_per_m:.4e}   (NASA nominal {RE_REF:.1e})")
    if abs(Re_per_m - RE_REF) / RE_REF > 0.10:
        print("  [warn] run NOT at NASA Re=5e6 -> set a constant mu=1.18587e-5")
        print("         (set_constant_transport.py) and g=59.29 to compare vs x.")

    xc, wcv = read_wall_tec(f_wall)
    tau = np.abs(wcv["tauX"])
    m = xc > 0.0
    xp, cf = xc[m], tau[m] / (0.5 * rho_inf * U_inf**2)

    # error vs the primary reference (CFL3D if present, else FUN3D), MOSE-style
    zones = read_tecplot_zones(f_cf)
    zname, ref = pick_zone(zones)
    x_hi = min(xp.max(), ref[:, 0].max())
    me = (xp >= X_LO) & (xp <= x_hi)
    cf_ref = np.interp(xp[me], ref[:, 0], ref[:, 1])
    err = np.abs(cf[me] - cf_ref) / cf_ref * 100.0
    rms, emax = float(np.sqrt(np.mean(err**2))), float(err.max())

    print(f"\n  Reference: '{zname}'   ({X_LO} <= x <= {x_hi:.2f} m)")
    print(f"    RMS relative error : {rms:.2f} %")
    print(f"    Max relative error : {emax:.2f} %")
    verdict = "PASS" if rms < 5.0 else "FAIL"
    print(f"  -> cf consistent with NASA-TMR reference : {verdict}  (tol 5% RMS)")

    # ---- plot Cf vs x (MOSE verify.py style) ----
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(xp, cf, "o", ms=5, mfc="none", color='red', label=f"ARES ({case_name})")
    for name, arr in reference_curves(zones):
        ax.plot(arr[:, 0], arr[:, 1], "-", lw=2, label=name)
    ax.set_xlim(0.0, 1.8)
    ax.set_ylim(0.002, 0.006)
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$C_f$")
    ax.set_title(f"Turbulent flat plate ({case_name}) - skin friction")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    return dict(U_inf=U_inf, rho_inf=rho_inf, mu_inf=mu_inf, Re_per_m=Re_per_m)
