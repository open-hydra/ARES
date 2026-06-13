#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the ARES turbulent flat-plate test case Flat_Plate_SGGLRR
(SSG-LRR Reynolds-Stress Model) against the NASA Turbulence Modeling
Resource (TMR).

Quantity reproduced from PLOT.lay:
  1) cf(x)                 vs reference/cf.dat

The mut/mu(y) profile plot was removed: its peak scales with Re, so it is
misleading unless the run is exactly at the NASA Re=5e6.

See _fp_turb_common.py for the Reynolds-consistency handling.

Usage:  python3 validate_sgglrr.py
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

# RSM field layout: Reynolds stresses replace k; mut (mit) still present.
CELL_VARS = ["p", "u", "v", "w", "h",
             "ruu", "rvv", "rww", "ruv", "ruw", "rvw", "omega",
             "T", "rho", "sound", "mil", "kl", "mit"]


def main():
    if not os.path.exists(F_FIELD):
        sys.exit(f"ERROR: {F_FIELD} not found - run the case first.")
    # ARES drawn in red: the SSG-LRR cf overlaps the (C0/blue) NASA reference
    # curve and was indistinguishable in the default colour.
    C.validate_cf(F_FIELD, F_WALL, F_CF, CELL_VARS, "SSG-LRR", ares_color="red")
    print("\nDone. Showing plots (close the windows to exit)...")
    plt.show()


if __name__ == "__main__":
    main()
