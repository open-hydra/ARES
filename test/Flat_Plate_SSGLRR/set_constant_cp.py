#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set a CONSTANT Cp in thermo.dat (FLINT real-fluid table) so that, together with
the constant mu and k imposed by set_constant_transport.py, the laminar Prandtl
number is exactly the MOSE/NASA value:

    Pr_lam = mu * Cp / k = 1.18587e-5 * 1004.5 / 1.65445e-2 = 0.72

thermo.dat is a Tecplot BLOCK file with 10 variables:
    "Pressure","Enthalpy","Density","Temperature","dRho/dT","dRho/dh",
    "Cp","Entropy","dRho/dp","SoundSpeed"
i.e. 3 header lines then 10 blocks of I*J*K values, one per line. Only the Cp
block is replaced; everything else (EOS: density, temperature, sound speed, ...)
is preserved, so T(p,h) etc. stay real-fluid - Cp here only sets the
conductivity/Prandtl coupling.

Original backed up to thermo.dat.orig (once). Restore with:
    cp INPUT/thermo.dat.orig INPUT/thermo.dat

Usage:  python3 set_constant_cp.py
"""

import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
F    = os.path.join(HERE, "INPUT", "thermo.dat")

# --------------------------------------------------------------------------- #
# constant to impose  (edit here)
# --------------------------------------------------------------------------- #
CP = 1004.5      # [J/kg/K]  perfect-gas cp = gamma*R/(gamma-1) = 1.4*287/0.4
                 # gives Pr = mu*Cp/k = 0.72 with the constant transport values


def main():
    with open(F) as fh:
        lines = fh.read().splitlines()

    header, data = lines[:3], lines[3:]
    names = re.findall(r'"([^"]+)"', header[1])      # VARIABLES line
    zone  = header[2]

    try:
        idx = next(i for i, n in enumerate(names) if n.strip().lower() == "cp")
    except StopIteration:
        raise ValueError(f"'Cp' not found in VARIABLES: {names}")

    def dim(tag):
        return int(re.search(tag + r"\s*=\s*(\d+)", zone).group(1))
    I, J, K = dim("I"), dim("J"), dim("K")
    nval = I * J * K
    nvar = len(names)
    if len(data) != nvar * nval:
        raise ValueError(f"expected {nvar*nval} data lines ({nvar} x {nval}), got {len(data)}")

    # report the old Cp range, then overwrite that block
    old = [float(v) for v in data[idx*nval:(idx+1)*nval]]
    data[idx*nval:(idx+1)*nval] = [repr(CP)] * nval

    bak = F + ".orig"
    if not os.path.exists(bak):
        shutil.copy2(F, bak)
        print(f"  backup written: {os.path.relpath(bak, HERE)}")
    else:
        print(f"  backup already exists: {os.path.relpath(bak, HERE)} (kept)")

    with open(F, "w") as fh:
        fh.write("\n".join(header + data) + "\n")

    print(f"  thermo.dat: Cp (variable #{idx+1}) set constant")
    print(f"    old Cp range : {min(old):.2f} .. {max(old):.2f} J/kg/K (real-fluid)")
    print(f"    new Cp       : {CP:.2f} J/kg/K (constant)")
    print(f"    other EOS variables (rho, T, sound, ...) preserved")
    print(f"  -> Pr_lam = mu*Cp/k = 0.72  (with mu=1.18587e-5, k=1.65445e-2)")


if __name__ == "__main__":
    main()
