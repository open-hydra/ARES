"""Case-directory custom.py of the NICFD nozzle: siloxane MDM (octamethyltrisiloxane) as a RealFluid.

Thermodynamics: the Helmholtz equation of state of CoolProp (HEOS, the multiparameter model of the paper,
Spinelli et al. 2018), evaluated with the gas phase imposed. The expansion starts 7 K above the dew line
(p0 = 904388 Pa, T0 = 542.13 K, Tsat(p0) = 535.0 K) and stays in the vapour (MDM is a dry fluid), but the
(p, h) rectangle of the table crosses the dew line in its high-p, low-h corner, which the flow never visits:
there the imposed gas phase gives the metastable vapour of the same equation, continuous with the rest of the
table; where even that fails the equilibrium state is taken, with finite-difference derivatives.

Transport: CoolProp has no viscosity or conductivity model for MDM, so the constants of the SU2 tutorial
(NICFD_nozzle.cfg: MU_CONSTANT, THERMAL_CONDUCTIVITY_CONSTANT) are used.
"""
import math

import CoolProp.CoolProp as CP

from custom_classes import RealFluid


class CoolPropVapour(RealFluid):
    """A CoolProp HEOS fluid tabulated as a (metastable) vapour, with constant transport properties."""

    def __init__(self, name, fluid, mu, k):
        super().__init__(name)
        self.mu, self.k = mu, k
        self.gas = CP.AbstractState('HEOS', fluid)
        self.gas.specify_phase(CP.iphase_gas)
        self.eq = CP.AbstractState('HEOS', fluid)

    def reference(self, T, p):
        self.gas.update(CP.PT_INPUTS, p, T)
        return self.gas.hmass(), self.gas.smass()

    def state(self, h, p):
        try:
            return self._state_gas(h, p)
        except ValueError:
            return self._state_equilibrium(h, p)

    def _T_gas(self, h, p):
        """T of the gas at (h, p) by Newton on h(T, p) (PT inputs: CoolProp's (h, p) flash with the phase imposed
        brackets T below the upper bound of the tables)."""
        g = self.gas
        T = 540.0
        for _ in range(50):
            g.update(CP.PT_INPUTS, p, T)
            dT = (h - g.hmass()) / g.cpmass()
            T += max(min(dT, 20.0), -20.0)
            if abs(dT) < 1.0e-9 * T:
                g.update(CP.PT_INPUTS, p, T)
                return T
        raise ValueError('gas temperature not converged')

    def _state_gas(self, h, p):
        g = self.gas
        self._T_gas(h, p)
        rho, T, cp = g.rhomass(), g.T(), g.cpmass()
        drho_dT = g.first_partial_deriv(CP.iDmass, CP.iT, CP.iP)
        drho_dp = g.first_partial_deriv(CP.iDmass, CP.iP, CP.iHmass)
        a = g.speed_sound()
        vals = (rho, T, cp, drho_dT, drho_dp, a)
        if not all(math.isfinite(v) for v in vals) or rho <= 0 or cp <= 0 or a <= 0:
            raise ValueError('gas phase not available')
        return {"T": T, "rho": rho, "drho_dT": drho_dT, "drho_dh": drho_dT / cp, "drho_dp": drho_dp,
                "cp": cp, "s": g.smass(), "mu": self.mu, "k": self.k, "a": a}

    def _state_equilibrium(self, h, p):
        e = self.eq

        def at(hh, pp):
            e.update(CP.HmassP_INPUTS, hh, pp)
            return e.rhomass(), e.T()

        rho, T = at(h, p)
        s = e.smass()
        dh, dp = 1.0e-4 * abs(h) + 10.0, 1.0e-5 * p
        rp, Tp = at(h + dh, p); rm, Tm = at(h - dh, p)
        drho_dh = (rp - rm) / (2 * dh)
        dT_dh = (Tp - Tm) / (2 * dh)
        drho_dp = (at(h, p + dp)[0] - at(h, p - dp)[0]) / (2 * dp)
        # In the dome T is a function of p alone: cp and d(rho)/dT|p are those of the saturated vapour
        if abs(dT_dh) < 1.0e-12:
            e.update(CP.PQ_INPUTS, p, 1.0)
            cp = e.cpmass()
            drho_dT = e.first_partial_deriv(CP.iDmass, CP.iT, CP.iP)
        else:
            cp = 1.0 / dT_dh
            drho_dT = drho_dh * cp
        # 1/a^2 = d(rho)/dp|s = d(rho)/dp|h + d(rho)/dh|p / rho  (dh = dp / rho along an isentrope)
        a2 = 1.0 / max(drho_dp + drho_dh / rho, 1.0e-12)
        return {"T": T, "rho": rho, "drho_dT": drho_dT, "drho_dh": drho_dh, "drho_dp": drho_dp,
                "cp": cp, "s": s, "mu": self.mu, "k": self.k, "a": math.sqrt(a2)}


MDM = CoolPropVapour("MDM", "MDM", mu=1.21409e-5, k=0.030542828)
