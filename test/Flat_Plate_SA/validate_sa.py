#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the ARES turbulent flat-plate test case (Flat_Plate_SA,
Spalart-Allmaras) against the NASA Turbulence Modeling Resource (TMR)
2D zero-pressure-gradient flat plate.

This is the turbulent counterpart of the laminar 'validate_blasius.py'. The
quantities reproduced here are exactly those of PLOT.lay:

    Tau_x = mu * |dU/dy|_wall ,   c_f = Tau_x / (0.5 * rho_inf * U_inf^2)

1.  Skin-friction coefficient  cf(x)   (source: OUTPUT/1d.dat)
    cf is built from the wall shear tau_w and the free-stream state and is
    compared with the NASA-TMR reference (reference/cf.dat, CFL3D grid level 1).

    IMPORTANT - Reynolds consistency:
    The NASA case is defined at Re = 5e6 (per unit length). This ARES run uses
    the real-fluid air table, whose viscosity at the operating point
    (T~555 K) is mu ~ 2.9e-5 Pa s, i.e. about 2.5x the ideal value the NASA
    case assumes; the run therefore sits at Re/m ~ 2.0e6, NOT 5e6. Comparing
    cf directly versus x would carry a ~18% offset that is pure Reynolds
    effect, not a model error. So - exactly as the laminar script evaluates
    the Blasius law at the run's own Re_x - here the tabulated NASA cf is
    re-expressed at the run's actual Re_x:  cf_ref(x) = cf_NASA( Re_x(x) ),
    with Re_x_NASA = RE_REF * x_ref. In Re_x space the curves collapse.

The files are produced by the run; this script only reads them.

Usage
-----
    python3 validate_sa.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# paths (relative to this script's directory)
# --------------------------------------------------------------------------- #
HERE    = os.path.dirname(os.path.abspath(__file__))
F_1D    = os.path.join(HERE, "OUTPUT", "1d.dat")
F_REF   = os.path.join(HERE, "reference", "cf.dat")

# Column indices in 1d.dat (0-based), see the "variables=" header line:
# 0:x 1:s 2:A 3:P 4:Dh 5:G 6:U 7:p 8:p0 9:T 10:T0 11:rho 12:h 13:s 14:Prw
# 15:cp 16:mi(=mu) 17:k 18:Tm 19:Re 20:Pr 21:tau3 22:tau4 23:tau5 24:tau6
# 25:tau(total) ...
C_X, C_U, C_RHO, C_MU, C_TAU = 0, 6, 11, 16, 25

# The wall (adiabatic flat plate) starts at x = 0; x in [-2, 0] is symmetry
# (input.ini, [wall] patch1=symmetry / patch2=adiabatic). The leading-edge
# region is under-resolved, so cf statistics/plot start at X_MIN.
X_MIN = 0.01

# NASA-TMR reference
RE_REF   = 5.0e6                 # reference Reynolds per unit length of the TMR case
REF_ZONE = "CFL3D, grid level 1"  # finest-grid reference curve


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def read_cf_reference(zone_wanted):
    """Parse the multi-zone Tecplot reference cf.dat and return (x, cf) of the
       requested zone. Zone headers look like  ZONE T="CFL3D, grid level 1"
       or  ZONE, T="FUN3D, grid level 1"."""
    zones, cur, rows = {}, None, []
    with open(F_REF) as fh:
        for ln in fh:
            if "ZONE" in ln:
                if cur is not None:
                    zones[cur] = np.array(rows)
                cur = ln.split("T=")[1].strip().strip(",").strip().strip('"')
                rows = []
            elif ln.strip() and not ln.startswith("#") and "VARIABLES" not in ln:
                try:
                    rows.append([float(v) for v in ln.split()])
                except ValueError:
                    pass
        if cur is not None:
            zones[cur] = np.array(rows)
    if zone_wanted not in zones:
        sys.exit("ERROR: zone '%s' not found in %s. Available: %s"
                 % (zone_wanted, F_REF, list(zones.keys())))
    ref = zones[zone_wanted]
    return ref[:, 0], ref[:, 1]


# --------------------------------------------------------------------------- #
# 1. cf validation from 1d.dat
# --------------------------------------------------------------------------- #
def validate_cf():
    print("=" * 74)
    print(" 1) SKIN-FRICTION COEFFICIENT  cf(x)   (source: OUTPUT/1d.dat)")
    print("=" * 74)

    data = np.genfromtxt(F_1D, skip_header=3)
    x   = data[:, C_X]
    U   = data[:, C_U]
    rho = data[:, C_RHO]
    mu  = data[:, C_MU]
    tau = data[:, C_TAU]

    # Free-stream / edge state from the symmetry strip upstream (x < 0).
    up = x < 0.0
    if up.sum() < 3:
        up = np.zeros_like(x, dtype=bool); up[:5] = True
    U_inf, rho_inf, mu_inf = U[up].mean(), rho[up].mean(), mu[up].mean()
    Re_per_m = rho_inf * U_inf / mu_inf

    print(f"  Free-stream state (x<0 strip):")
    print(f"    U_inf   = {U_inf:10.4f}  m/s")
    print(f"    rho_inf = {rho_inf:10.6f} kg/m^3")
    print(f"    mu_inf  = {mu_inf:10.4e} Pa s   (real-fluid air table)")
    print(f"    Re/m    = {Re_per_m:10.4e} 1/m  (NASA nominal: {RE_REF:.1e})")

    # Developed wall region: x >= X_MIN with a real wall shear.
    pl = (x >= X_MIN) & (tau > 0.0)
    xp, taup = x[pl], tau[pl]
    Re_x   = rho_inf * U_inf * xp / mu_inf
    cf_cfd = taup / (0.5 * rho_inf * U_inf**2)

    # NASA reference re-expressed at the run's own Re_x (Reynolds-consistent,
    # same idea as evaluating Blasius at the run's Re_x in the laminar case).
    xr, cfr = read_cf_reference(REF_ZONE)
    Re_x_ref = RE_REF * xr
    cf_ref   = np.interp(Re_x, Re_x_ref, cfr)
    err      = (cf_cfd - cf_ref) / cf_ref * 100.0
    mae      = np.mean(np.abs(err))

    print("\n  Sample of the comparison (x >= %.2f m):" % X_MIN)
    print("    {:>10s} {:>12s} {:>12s} {:>12s} {:>9s}".format(
        "x[m]", "Re_x", "cf_CFD", "cf_NASA", "err[%]"))
    for i in np.linspace(0, len(xp) - 1, 10).astype(int):
        print("    {:10.4f} {:12.3e} {:12.4e} {:12.4e} {:+9.2f}".format(
            xp[i], Re_x[i], cf_cfd[i], cf_ref[i], err[i]))

    print(f"\n  Mean |error| for x >= {X_MIN:.2f} m : {mae:6.2f} %  "
          f"(vs {REF_ZONE})")
    verdict = "PASS" if mae < 5.0 else "CHECK"
    print(f"  -> cf consistent with NASA-TMR SA reference : {verdict}")

    # ---- plot cf vs x ----
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xp, cf_cfd, "o", ms=4, mfc="none", color="red", label="ARES (SA)")
    xs    = np.linspace(xp.min(), xp.max(), 400)
    Re_xs = rho_inf * U_inf * xs / mu_inf
    ax.plot(xs, np.interp(Re_xs, Re_x_ref, cfr), "-", color="k",
            label="NASA-TMR CFL3D (at run Re$_x$)")
    ax.set_xlabel(r"$x$ [m]   (leading edge at $x=0$)")
    ax.set_ylabel(r"$c_f$")
    ax.set_yscale("log")
    ax.set_title("Turbulent flat plate (SA) - skin friction vs x")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()

    return dict(U_inf=U_inf, rho_inf=rho_inf, mu_inf=mu_inf)


# --------------------------------------------------------------------------- #
def main():
    if not os.path.exists(F_1D):
        sys.exit(f"ERROR: {F_1D} not found - run the case first.")
    validate_cf()
    print()
    print("Done. Showing plot (close the window to exit)...")
    plt.show()


if __name__ == "__main__":
    main()
