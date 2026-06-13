#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the ARES turbulent flat-plate test case Flat_Plate_SST
(SST k-omega) against the NASA Turbulence Modeling Resource (TMR).

Quantities reproduced from PLOT.lay:
  1) cf(x)                          vs reference/cf.dat
  2) k+(y), omega+(y) at x~0.97     vs reference/kw.dat
     k+     = k / a_inf^2
     omega+ = omega * mu_inf / (rho_inf * a_inf^2)

(The mut/mu profile plot was removed: its peak scales with Re, so it is
 misleading unless the run is exactly at the NASA Re=5e6.)

See _fp_turb_common.py for the Reynolds-consistency handling: cf is compared
in Re_x space; the x~0.97 profiles only collapse if the run is at Re = 5e6.

Usage:  python3 validate_sst.py
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))   # ../common
import _fp_turb_common as C

F_FIELD = os.path.join(HERE, "OUTPUT", "field.tec")
F_WALL  = os.path.join(HERE, "OUTPUT", "wall.tec")
F_CF    = os.path.join(HERE, "reference", "cf.dat")
F_MUT   = os.path.join(HERE, "reference", "mut.dat")
F_KW    = os.path.join(HERE, "reference", "kw.dat")

# field.tec cell variables (after x,y,z nodal)
CELL_VARS = ["p", "u", "v", "w", "h", "kappa", "omega",
             "T", "rho", "sound", "mil", "kl", "mit"]


def validate_profiles():
    print()
    print("=" * 74)
    print(" 2) PROFILES at x ~ %.2f m   (source: OUTPUT/field.tec)" % C.X_STATION)
    print("=" * 74)
    if not os.path.exists(F_FIELD):
        print("  field.tec not found - skipping profiles."); return

    _, _, _, xN, yN, cv = C.read_field_tec(F_FIELD, CELL_VARS)
    x_sta, yn, p = C.wall_profile(xN, yN, cv,
                                  ["u", "rho", "mil", "mit", "kappa", "omega", "sound"],
                                  C.X_STATION)
    # free-stream edge (outermost cell of the column)
    rho_inf, mu_inf, a_inf = p["rho"][-1], p["mil"][-1], p["sound"][-1]

    # ARES profiles (kappa, omega in field are conserved: rho*k, rho*omega)
    k_spec = p["kappa"] / p["rho"]
    om_spec = p["omega"] / p["rho"]
    k_nd  = k_spec / a_inf**2
    om_nd = om_spec * mu_inf / (rho_inf * a_inf**2)

    print(f"  station x = {x_sta:.4f} m")

    # ----- k+, omega+ vs kw.dat -----
    have_kw = os.path.exists(F_KW)
    if have_kw:
        zkw = C.pick_zone(C.read_tecplot_zones(F_KW))
        ykw, omkw, kkw = zkw[1][:, 1], zkw[1][:, 2], zkw[1][:, 3]
        print(f"  k+      : ARES peak = {k_nd.max():.3e}   ref peak = {kkw.max():.3e}")
        print(f"  omega+  : ARES max  = {om_nd.max():.3e}   ref max  = {omkw.max():.3e}")

    if abs(C.RE_REF - rho_inf * p['u'][-1] / mu_inf) / C.RE_REF > 0.10:
        print("  [note] run not at NASA Re=5e6 -> profiles are expected to differ")
        print("         (e.g. mut/mu peak scales with Re).")

    # ----- plots -----
    ymax = 0.04
    if have_kw:
        fig2, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.5))
        a1.plot(k_nd, yn, "o", ms=3, mfc="none", color="C2", label="ARES")
        a1.plot(kkw, ykw, "-", color="k", label=f"NASA {zkw[0].split(',')[0]}")
        a1.set_xlabel(r"$k^+ = k/a_\infty^2$"); a1.set_ylabel(r"$y$ [m]")
        a1.set_ylim(0, ymax); a1.grid(True, ls=":", alpha=0.6); a1.legend()
        a1.set_title("SST - turbulent kinetic energy")
        a2.semilogx(om_nd, yn, "o", ms=3, mfc="none", color="C3", label="ARES")
        a2.semilogx(omkw, ykw, "-", color="k", label=f"NASA {zkw[0].split(',')[0]}")
        a2.set_xlabel(r"$\omega^+ = \omega\,\mu_\infty/(\rho_\infty a_\infty^2)$")
        a2.set_ylabel(r"$y$ [m]")
        a2.set_ylim(0, ymax); a2.grid(True, which="both", ls=":", alpha=0.6); a2.legend()
        a2.set_title("SST - specific dissipation rate")
        fig2.tight_layout()


def main():
    if not os.path.exists(F_FIELD):
        sys.exit(f"ERROR: {F_FIELD} not found - run the case first.")
    C.validate_cf(F_FIELD, F_WALL, F_CF, CELL_VARS, "SST")
    validate_profiles()
    print("\nDone. Showing plots (close the windows to exit)...")
    plt.show()


if __name__ == "__main__":
    main()
