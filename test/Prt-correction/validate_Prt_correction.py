#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Turbulent Prandtl correction model validation (Prt-correction)
================================================================================

Reference:
  B. Latini, M. Fiore, F. Nasuti,
  "Modeling liquid rocket engine coolant flow and heat transfer in high
   roughness channels", Aerospace Science and Technology 126 (2022) 107672.

Test case (axisymmetric, uniformly heated channel, high roughness) -- Prt-correction2:
  - Real fluid: water (Pr ~ 6.033)
  - Re (nominal, inlet) ~ 2.7356e4   [paper Sec.4.1; local Re grows along x]
  - hs/D = 0.08   (relative roughness, equivalent sand grain; read from 1d.dat)
  - inlet pressure p ~ 50 bar, imposed wall heat flux q_w = 3 MW/m^2
  Corresponds to the calibration point of Fig. 9(a)/Fig.15 of the paper.

NOTE on the Reynolds number:
  The strong wall heating lowers the viscosity downstream, so the LOCAL Re in
  1d.dat is far from constant.  The reference correlations (Colebrook,
  Dipprey-Sabersky) are evaluated at a SINGLE Reynolds number: the value taken
  at the smallest x of the region used for the correlations, i.e. the start of
  the developed region (x >= 0.8 L).  All reference quantities are therefore
  single scalars (horizontal lines); the per-section f_D and Nu coming from the
  1D file are what gets validated against them.  The Darcy friction factor of
  the reference is computed with Colebrook-White (Eq.9) only.

The model ("THRC", thermal high-roughness correction) introduces a local
turbulent Prandtl number Pr_t = Pr_ts + dPr_tr that GROWS approaching the
wall, reducing the turbulent conductivity k_t and hence the heat transfer,
breaking the Reynolds analogy typical of the "equivalent sand grain" approach.
Without correction the Nu is strongly overpredicted (factor that grows with
hs/D: ~2.2x at hs/D=0.04 up to >4x at hs/D=0.21, Fig.5 of the EUCASS paper);
with correction the Nu falls back in line with the experimental correlations
(Dipprey-Sabersky), i.e. Nu/Nu_theo ~ 1.

This script compares the ARES results (file OUTPUT/1d.dat, section-averaged
data along x) with the theoretical correlations of the paper and with the
REFERENCE values digitized from the article figures, producing two
non-dimensional plots in the style of Fig. 5 / Fig. 9 of the paper:

  1) Non-dimensional friction factor  f_D / f_D,theo  vs  x/D   (Fig. 5 style)
        f_D,theo  -> Colebrook-White for rough pipes, Eq.(9)
        ( equivalent to  c_f / c_f,theo , since c_f = f_D/4 )

  2) Non-dimensional Nusselt number  Nu / Nu_theo    vs  x/D   (Fig. 9 style)
        Nu_theo   -> Dipprey-Sabersky, Eq.(11), via Stanton number

The plots are shown on screen (plt.show()), not saved.
================================================================================
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt

try:
    from scipy.optimize import brentq
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

# ------------------------------------------------------------------------------
# Case parameters (read from input.ini where possible)
# ------------------------------------------------------------------------------
HERE     = os.path.dirname(os.path.abspath(__file__))
DATFILE  = os.path.join(HERE, "OUTPUT", "1d.dat")
INI      = os.path.join(HERE, "input.ini")

KARMAN   = 0.41        # von Karman constant (not used directly here)
KF_DS    = 5.19        # k_f coefficient of the Dipprey-Sabersky correlation [8]
PR_TS    = 0.9         # smooth-wall turbulent Prandtl number (info)

# x/L at which the SINGLE reference Reynolds number is evaluated (also the start
# of the developed region used for the validation averages).  Change this knob
# to move the reference station along the channel and check how Re_ref varies.
# (paper design Reynolds, Sec.4.1, is ~2.7356e4; here Re_ref comes from 1d.dat.)
X_REF_FRAC = 0.90

# ------------------------------------------------------------------------------
# REFERENCE values from the paper (digitized from the figures)
# Case: Pr = 6.033, hs/D = 0.08, Re_nominal = 2.7356e4
# ------------------------------------------------------------------------------
# Fully-developed Nu/Nu_theo: the reference is the Dipprey-Sabersky correlation
# itself (Eq.11), i.e. the corrected model should give Nu/Nu_theo ~ 1.
REF_NU_RATIO_EQ11     = 1.0    # theoretical reference Dipprey-Sabersky (Eq.11)
REF_NU_RATIO_CORR     = 1.0    # ARES paper, "SA - Pr_t correction" -> ~1
# Fig. 5(b) / Sec.4 text: developed f_D/f_D,theo (~5% below Colebrook at hs/D=0.08)
REF_FD_RATIO_EXTRAP    = 0.95  # extrapolated value (Richardson), ~5% shift


def read_ks_from_ini(path, default=1.60e-4):
    """Read the equivalent sand-grain height hs (= ks) from input.ini.
    The wall [qw] block contains the line 'ks = ...'."""
    if not os.path.isfile(path):
        return default
    for line in open(path):
        m = re.match(r"\s*ks\s*=\s*([0-9.eEdD+-]+)", line)
        if m:
            val = m.group(1).replace("d", "e").replace("D", "e")
            try:
                return float(val)
            except ValueError:
                pass
    return default


def roughness_from_1d(d, ini_path):
    """Equivalent sand-grain roughness hs [m] along x.

    Preferred source: the 'hs[m]' column now written by ARES into 1d.dat (per
    wall cell).  Falls back to parsing 'ks' from input.ini if the column is
    absent or all-zero (e.g. an old 1d.dat from before the solver change)."""
    if "hs[m]" in d and np.any(d["hs[m]"] > 0.0):
        return d["hs[m]"].copy()
    hs = read_ks_from_ini(ini_path)
    print("   (hs[m] column missing/zero in 1d.dat -> using ks=%.3e from input.ini)" % hs)
    return np.full_like(d["x[m]"], hs)


def load_1d(path):
    """Load OUTPUT/1d.dat (Tecplot POINT) mapping the columns by name."""
    lines = open(path).read().splitlines()
    hdr = next(l for l in lines if l.strip().lower().startswith("variables"))
    names = re.findall(r'"([^"]+)"', hdr)
    rows = []
    for l in lines:
        s = l.split()
        if len(s) == len(names):
            try:
                rows.append([float(x) for x in s])
            except ValueError:
                continue
    a = np.array(rows)
    return {n: a[:, i] for i, n in enumerate(names)}


# ------------------------------------------------------------------------------
# Theoretical correlations
# ------------------------------------------------------------------------------
def colebrook_white(Re, rough):
    """Darcy friction factor f_D from Colebrook-White, Eq.(9):
        1/sqrt(f) = -2 log10( hs/(3.71 D) + 2.51/(Re sqrt(f)) )
    'rough' = hs/D (relative roughness)."""
    def residual(f):
        return 1.0 / np.sqrt(f) + 2.0 * np.log10(rough / 3.71 + 2.51 / (Re * np.sqrt(f)))
    if _HAVE_SCIPY:
        return brentq(residual, 1e-5, 5.0, maxiter=200)
    # fallback: fixed-point iteration
    f = 0.05
    for _ in range(100):
        f = (1.0 / (-2.0 * np.log10(rough / 3.71 + 2.51 / (Re * np.sqrt(f))))) ** 2
    return f


def stanton_dipprey_sabersky(fD, Re, Pr, rough, kf=KF_DS):
    """Stanton number from the Dipprey-Sabersky correlation, Eq.(11):
        St = (f/8) / ( 1 + sqrt(f/8) * { kf [Re (hs/D) sqrt(f/8)]^0.2 Pr^0.44 - 8.48 } )
    Links the heat transfer to the measured friction (breaks the Reynolds analogy).
    """
    fr = fD / 8.0
    roughness_reynolds = Re * rough * np.sqrt(fr)
    bracket = kf * roughness_reynolds ** 0.2 * Pr ** 0.44 - 8.48
    return fr / (1.0 + np.sqrt(fr) * bracket)


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    if not os.path.isfile(DATFILE):
        raise SystemExit("File not found: %s\n(run the ARES solver first)" % DATFILE)

    d  = load_1d(DATFILE)

    # --- section quantities ---
    x    = d["x[m]"]
    D    = d["Dh[m]"]                 # hydraulic diameter (= pipe diameter, axisym.)
    rho  = d["rho[kg/m^3]"]
    U    = d["U[m/s]"]                # bulk velocity
    Tb   = d["T[K]"]                  # bulk (static) temperature
    Tw   = d["Tw[K]"]                 # wall temperature
    qw   = d["qw[W/m^2]"]            # wall heat flux
    kth  = d["k[W/m K]"]             # thermal conductivity
    tau  = d["tau[Pa]"]              # wall shear stress
    Re   = d["Re"]                    # LOCAL section Reynolds number
    Pr   = d["Pr"]
    hs   = roughness_from_1d(d, INI)  # sand-grain roughness hs [m] (from 1d.dat)

    xD       = x / D                  # non-dimensional abscissa
    rough    = hs / D                 # relative roughness hs/D (per section)

    # ------------------------------------------------------------------
    # Single reference Reynolds number.
    # Taken at x/L = X_REF_FRAC (start of the developed region).  ALL reference
    # correlations are evaluated only at this Re -> single scalar reference
    # values.  The 1D-file quantities (f_D, Nu) stay per-section: they are what
    # gets validated.  Move X_REF_FRAC to change the reference station.
    # ------------------------------------------------------------------
    L     = x.max()                   # channel length (mesh starts at x ~ 0)
    xL    = x / L                     # axial fraction x/L
    i_ref = int(np.argmin(np.abs(xL - X_REF_FRAC)))   # station closest to X_REF_FRAC
    dev   = np.arange(len(x)) >= i_ref  # developed region: from the reference station on
    Re_ref    = Re[i_ref]
    Pr_ref    = Pr[i_ref]
    rough_ref = rough[i_ref]

    # ------------------------------------------------------------------
    # 1) Friction
    # ------------------------------------------------------------------
    # Numerical (validated, per section), from wall shear stress:
    #   f_D = 8 tau_w / (rho U^2)
    fD_wall = 8.0 * tau / (rho * U ** 2)
    # Reference (single value): Colebrook-White Eq.(9) at Re_ref
    fD_theo_ref = colebrook_white(Re_ref, rough_ref)

    cf_ratio = fD_wall / fD_theo_ref   # f_D/f_D,theo  (from tau_w)

    # ------------------------------------------------------------------
    # 2) Heat transfer
    # ------------------------------------------------------------------
    # Numerical (validated, per section): h = q_w/(T_w - T_b) Eq.(10), Nu = hD/k
    h     = qw / (Tw - Tb)
    Nu    = h * D / kth
    # Reference (single value): Dipprey-Sabersky Eq.(11) at Re_ref, using the
    # Colebrook-White friction factor (Eq.9) -- NOT the numerical wall friction.
    St_ref      = stanton_dipprey_sabersky(fD_theo_ref, Re_ref, Pr_ref, rough_ref)
    Nu_theo_ref = St_ref * Re_ref * Pr_ref

    nu_ratio = Nu / Nu_theo_ref

    # ------------------------------------------------------------------
    # Plots (shown on screen)
    # ------------------------------------------------------------------
    relabel = r"$Re_{ref}=%.2g$ @ $x/D=%.0f$" % (Re_ref, xD[i_ref])

    # Fig 1 -- non-dimensional friction (Fig. 5 style)
    fig1, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xD, cf_ratio, "b-", lw=1.8, label=r"ARES")
    ax.axhline(REF_FD_RATIO_EXTRAP, color="grey", ls="--", lw=1.4,
               label=r"Paper reference")
    ax.set_xlabel(r"$x/D$")
    ax.set_ylabel(r"$f_D\,/\,f_{D,\,theo}\;\;(\equiv c_f/c_{f,\,theo})$")
    ax.set_title(r"Non-dimensional friction ($h_s/D=%.2f$, %s)"
                 % (rough_ref, relabel))
    ax.set_ylim(0.0, 2.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig1.tight_layout()

    # Fig 2 -- non-dimensional Nusselt (Fig. 9 style)
    fig2, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xD, nu_ratio, "r-", lw=1.8, label=r"ARES ($SA$ - $Pr_t$ correction)")
    ax.axhline(REF_NU_RATIO_EQ11, color="k", ls=":", lw=1.4,
               label="Dipprey-Sabersky ")
    ax.set_xlabel(r"$x/D$")
    ax.set_ylabel(r"$Nu\,/\,Nu_{theo}$")
    ax.set_title(r"Non-dimensional Nusselt ($h_s/D=%.2f$, $Pr=%.2f$, %s)"
                 % (rough_ref, Pr_ref, relabel))
    ax.set_ylim(0.0, 5.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig2.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
