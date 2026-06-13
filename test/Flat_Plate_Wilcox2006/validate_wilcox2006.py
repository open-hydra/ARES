#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the ARES turbulent flat-plate test case Flat_Plate_Wilcox2006
(Wilcox k-omega 2006) against the NASA Turbulence Modeling Resource (TMR).

Quantities reproduced from PLOT.lay:
  1) cf(x)                 vs reference/cf.dat
  2) mut/mu(y) at x~0.97   vs reference/mut.dat   (mut/mu = kappa/omega*rho/mu)

See _fp_turb_common.py for the Reynolds-consistency handling: cf is compared
in Re_x space; the x~0.97 profile only collapses if the run is at Re = 5e6.

Usage:  python3 validate_wilcox2006.py
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

CELL_VARS = ["p", "u", "v", "w", "h", "kappa", "omega",
             "T", "rho", "sound", "mil", "kl", "mit"]


def validate_profile():
    print()
    print("=" * 74)
    print(" 2) EDDY-VISCOSITY RATIO  mut/mu(y) at x ~ %.2f m" % C.X_STATION)
    print("=" * 74)
    if not os.path.exists(F_FIELD):
        print("  field.tec not found - skipping profile."); return

    _, _, _, xN, yN, cv = C.read_field_tec(F_FIELD, CELL_VARS)
    x_sta, yn, p = C.wall_profile(xN, yN, cv, ["u", "mil", "mit"], C.X_STATION)
    mut_mu = p["mit"] / p["mil"]

    zmut = C.pick_zone(C.read_tecplot_zones(F_MUT))
    ymut, mutref = zmut[1][:, 1], zmut[1][:, 2]
    print(f"  station x = {x_sta:.4f} m")
    print(f"  mut/mu : ARES peak = {mut_mu.max():7.1f}   "
          f"ref({zmut[0]}) peak = {mutref.max():7.1f}")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.plot(mut_mu, yn, "o", ms=3, mfc="none", color="C1", label="ARES")
    ax.plot(mutref, ymut, "-", color="k", label=f"NASA {zmut[0].split(',')[0]}")
    ax.set_xlabel(r"$\mu_t/\mu$"); ax.set_ylabel(r"$y$ [m]")
    ax.set_ylim(0, 0.04)
    ax.set_title(f"Wilcox2006 - eddy-viscosity ratio at x={x_sta:.2f} m")
    ax.grid(True, ls=":", alpha=0.6); ax.legend()
    fig.tight_layout()


def main():
    if not os.path.exists(F_FIELD):
        sys.exit(f"ERROR: {F_FIELD} not found - run the case first.")
    C.validate_cf(F_FIELD, F_WALL, F_CF, CELL_VARS, "Wilcox2006")
    validate_profile()
    print("\nDone. Showing plots (close the windows to exit)...")
    plt.show()


if __name__ == "__main__":
    main()
