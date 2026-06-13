# Real-Fluid Thermodynamics

## Overview

ARES does **not** assume an ideal gas. The equation of state and the transport properties are supplied as **tabulated real-fluid data** on a structured **pressure–enthalpy $(p,h)$ grid**, generated and interpolated by the [FLINT](https://github.com/MarcoGrossi92/FLINT) library. This lets ARES reproduce the genuine real-fluid behaviour of a chosen working fluid (compressibility factor $\neq 1$, variable specific heats, near-critical effects) over the tabulated range, while keeping the per-cell property evaluation as cheap as a bilinear lookup.

The table is requested from a `[GPB-*]` block in `input.ini`:

```ini
[GPB-Phase1]
type  = real-fluid
fluid = air
pmin  = 0.80e5   ;  pmax = 2.0e5      ! pressure range [Pa]
Tmin  = 280.0    ;  Tmax = 340.0      ! temperature range [K]
NP    = 200      ;  NH   = 200        ! resolution: NP × NH points
```

ATLAS / FLINT build the data over $[p_\min,p_\max]\times$ the enthalpy range corresponding to $[T_\min,T_\max]$, sampling reference data (e.g. CoolProp or NASA correlations) at $N_P \times N_H$ points and writing `thermo.dat` (thermodynamics) and the transport table.

!!! warning "Stay inside the table"
    Properties are only valid inside the tabulated $(p,h)$ box. A simulation whose pressure or enthalpy leaves $[p_\min,p_\max]$ or the enthalpy range corresponding to $[T_\min,T_\max]$ will extrapolate, which is inaccurate and can destabilise the run. Choose the `[GPB]` range to comfortably bracket the expected flow states.

---

## The $(p,h)$ Table

The state variable is the pair $(p,h)$. Every other thermodynamic or transport quantity is stored as a 2-D array over the $(p,h)$ grid and retrieved by **bilinear interpolation**, denoted here `ph2vars(p,h,·)`:

| Quantity | Symbol | Table |
|----------|:------:|-------|
| Density | $\rho$ | `rho_tab` |
| Temperature | $T$ | `T_tab` |
| Speed of sound | $a$ | `sound_tab` |
| Specific heat ratio | $\gamma$ | `gamma_tab` |
| $\partial h/\partial T\big|_p$ ($=c_p$) | $c_p$ | `hT_tab` |
| Dynamic viscosity | $\mu_\ell$ | transport table |
| Thermal conductivity | $k_\ell$ | transport table |

For example the speed of sound used by the time-step routine is simply

$$
a = \texttt{ph2vars}(p,\,h,\,\texttt{sound\_tab}).
$$

!!! note "Specific heat lives in the enthalpy derivative"
    At constant pressure $c_p = \partial h/\partial T\big|_p$, so on a $(p,h)$ grid the specific heat is obtained from the temperature/enthalpy derivative table (`hT_tab`), **not** a separate `cp_tab`. This is a common pitfall when reading or post-processing the tables directly.

### Table orientation

`thermo.dat` and the transport table **must share the same $(p,h)$ orientation** (the same ordering of the pressure and enthalpy axes). A mismatch between the two is treated as an error during loading, not silently auto-detected — the two files are assumed consistent because they are produced together by ATLAS/FLINT.

---

## Using the Table in the Solver

The real-fluid EOS enters the solver at three points:

1. **Primitive → conservative.** Given $(p,h)$, density $\rho=\texttt{ph2vars}(p,h,\texttt{rho\_tab})$ closes the conservative vector $\mathbf U = [\rho,\rho\mathbf u,\rho E_0]$ with $e=h-p/\rho$.

2. **Conservative → primitive (thermo inversion).** From $(\rho,e)$ recovered after a time step, ARES inverts the table to find the $(p,h)$ consistent with the stored $\rho(p,h)$ and $e(p,h)=h-p/\rho$. This local two-variable inversion is the real-fluid replacement for the ideal-gas temperature Newton iteration, and is the "Thermo inversion … OK" line in the loading report.

3. **Fluxes and time step.** The speed of sound, $\gamma$, $c_p$, and the transport coefficients needed by the [Riemann solvers](riemann-solvers.md), the [viscous fluxes](governing-equations.md#viscous-fluxes), the [CFL/VNN time step](time-integration.md#stability-conditions), and the [preconditioner](preconditioning.md) are all read from the table at the local $(p,h)$.

---

## Transport Properties

For a single real fluid the molecular viscosity and conductivity are direct table look-ups:

$$
\mu_\ell = \mu_\ell(p,h), \qquad
k_\ell = k_\ell(p,h).
$$

The **effective** transport coefficients used in the viscous fluxes add the turbulent contributions from the active RANS model:

$$
\mu = \mu_\ell + \mu_t,
\qquad
\kappa = k_\ell + \frac{\mu_t\,c_p}{\mathrm{Pr}_t},
$$

with the turbulent Prandtl number $\mathrm{Pr}_t$ taken from `[ARES-RANS]` (default 0.90). An optional [wall-roughness Prandtl correction](turbulence.md) (`Prt-correction`) modifies $\mathrm{Pr}_t$ near rough walls.

---

## Validation of the Tables

The FLINT real-fluid tables are verified independently of ARES: a self-consistency check confirms the tabulated $(\rho,T,a,\dots)$ satisfy the thermodynamic relations to machine precision, and the values are cross-checked against reference databases (CoolProp). Validation scripts live under `lib/FLINT/test/real-fluid/`.

!!! tip "Constant-property checks"
    The bundled flat-plate cases ship helper scripts (`set_constant_cp.py`, `set_constant_transport.py`) that flatten the table to constant $c_p$ / constant transport. These reduce the real-fluid case to a textbook ideal-gas / constant-property problem so the solution can be checked against an analytic result (e.g. Blasius) without table-interpolation effects.

---

## References

1. Bell, I. H., Wronski, J., Quoilin, S., Lemort, V. "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp." *Ind. Eng. Chem. Res.* 53 (2014).
2. Poinsot, T., Veynante, D. *Theoretical and Numerical Combustion*, 3rd ed., 2012.
3. FLINT documentation — [https://github.com/MarcoGrossi92/FLINT](https://github.com/MarcoGrossi92/FLINT).
