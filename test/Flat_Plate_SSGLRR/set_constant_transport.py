#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overwrite the transport.dat (FLINT real-fluid table) of this case with
CONSTANT molecular transport properties, so the flat plate runs like the
MOSE / NASA-TMR ideal-gas setup:

    Re/m = g / mu = 59.29 / 1.18587e-5 = 5.0e6   (at M = 0.2)

instead of the real-fluid viscosity (~2.9e-5 at 555 K) that would put the run
at Re/m ~ 2.0e6. Only the Viscosity and Conductivity blocks are replaced; the
Pressure/Enthalpy grid is preserved verbatim, so the table stays consistent
with thermo.dat.

transport.dat is a Tecplot BLOCK file:
    VARIABLES = "Pressure","Enthalpy","Viscosity","Conductivity"
    ZONE I=..,J=..,K=1, DATAPACKING=BLOCK
i.e. 3 header lines then 4 blocks of I*J*K values (P, H, mu, k), one per line.

The original is backed up to transport.dat.orig (once). Re-run to re-apply;
restore with:  cp INPUT/transport.dat.orig INPUT/transport.dat

Usage:  python3 set_constant_transport.py
"""

import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
F    = os.path.join(HERE, "INPUT", "transport.dat")

# --------------------------------------------------------------------------- #
# constants to impose  (edit here)
# --------------------------------------------------------------------------- #
VISCOSITY = 1.1858685985e-5        # [Pa s] molecular viscosity (MOSE/NASA value)
PRANDTL   = 0.72                   # [-]    laminar Prandtl number
CP        = 1004.5                 # [J/kg/K] cp = gamma*R/(gamma-1) = 1.4*287/0.4
CONDUCTIVITY = VISCOSITY * CP / PRANDTL   # [W/m/K]  k = mu*cp/Pr  (or set directly)


def main():
    with open(F) as fh:
        lines = fh.read().splitlines()

    header, data = lines[:3], lines[3:]
    zone = header[2]

    def dim(tag):
        m = re.search(tag + r"\s*=\s*(\d+)", zone)
        if not m:
            raise ValueError(f"cannot find {tag} in ZONE line: {zone}")
        return int(m.group(1))

    I, J, K = dim("I"), dim("J"), dim("K")
    nval = I * J * K
    if len(data) != 4 * nval:
        raise ValueError(f"expected {4*nval} data lines (4 x {nval}), got {len(data)}")

    P = data[0:nval]              # Pressure  (kept verbatim)
    H = data[nval:2 * nval]       # Enthalpy  (kept verbatim)
    visc = [repr(VISCOSITY)]    * nval
    cond = [repr(CONDUCTIVITY)] * nval

    bak = F + ".orig"
    if not os.path.exists(bak):
        shutil.copy2(F, bak)
        print(f"  backup written: {os.path.relpath(bak, HERE)}")
    else:
        print(f"  backup already exists: {os.path.relpath(bak, HERE)} (kept)")

    with open(F, "w") as fh:
        fh.write("\n".join(header + P + H + visc + cond) + "\n")

    print(f"  transport.dat updated  (grid {I}x{J}x{K}, {nval} pts):")
    print(f"    Viscosity    = {VISCOSITY:.6e} Pa s   (constant)")
    print(f"    Conductivity = {CONDUCTIVITY:.6e} W/m/K (constant, Pr={PRANDTL}, cp={CP})")
    print(f"    Pressure / Enthalpy grid preserved")
    print(f"  -> with g=59.29: Re/m = g/mu = {59.29/VISCOSITY:.3e}  (target 5.0e6, M=0.2)")


if __name__ == "__main__":
    main()
