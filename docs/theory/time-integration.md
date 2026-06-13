# Time Integration

ARES advances the semi-discrete equations with an explicit strong-stability-preserving (SSP) Runge–Kutta scheme. This page covers the time-stepping algorithm, the CFL and VNN stability conditions, implicit residual smoothing (IRS) for convergence acceleration, and the geometric multigrid strategy. Low-Mach preconditioning of the pseudo-time derivative is described on its [own page](preconditioning.md).

---

## Explicit Runge–Kutta Schemes

Three explicit schemes are available via `time-scheme`, all written in the unified Shu–Osher convex-combination form:

| `time-scheme` | Scheme | Order | Stages |
|---------------|--------|:-----:|:------:|
| `euler` | Forward Euler | 1 | 1 |
| `RK2` | 2-stage SSP Runge–Kutta (Heun) | 2 | 2 |
| `RK3` | 3-stage SSP Runge–Kutta (Shu–Osher) | 3 | 3 |

All three are SSP: each stage is a convex combination of forward-Euler steps, so the positivity and TVD properties of the spatial discretization are inherited under the standard CFL limit.

**Forward Euler**

$$
\mathbf U^{n+1} = \mathbf U^n + \Delta t\,\mathcal L(\mathbf U^n)
$$

**SSP-RK2 (Heun)**

$$
\begin{aligned}
\mathbf U^{(1)} &= \mathbf U^n + \Delta t\,\mathcal L(\mathbf U^n)\\
\mathbf U^{n+1} &= \tfrac12\mathbf U^n + \tfrac12\bigl(\mathbf U^{(1)} + \Delta t\,\mathcal L(\mathbf U^{(1)})\bigr)
\end{aligned}
$$

**SSP-RK3 (Shu–Osher)**

$$
\begin{aligned}
\mathbf U^{(1)} &= \mathbf U^n + \Delta t\,\mathcal L(\mathbf U^n)\\
\mathbf U^{(2)} &= \tfrac34\mathbf U^n + \tfrac14\bigl(\mathbf U^{(1)} + \Delta t\,\mathcal L(\mathbf U^{(1)})\bigr)\\
\mathbf U^{n+1} &= \tfrac13\mathbf U^n + \tfrac23\bigl(\mathbf U^{(2)} + \Delta t\,\mathcal L(\mathbf U^{(2)})\bigr)
\end{aligned}
$$

### State update per stage

At each RK stage:

1. Form the conservative state from primitives via the real-fluid table, $\mathbf U = \text{prim2cons}(\mathbf P)$.
2. Scale the residual, $\mathbf R^\ast = -\mathbf R\,\Delta t / V$.
3. *(Optional)* apply [implicit residual smoothing](#implicit-residual-smoothing-irs) to $\mathbf R^\ast$.
4. Combine according to the RK coefficients.
5. Recover primitives by **thermo inversion**, $\mathbf P^{(k)} = \text{cons2prim}(\mathbf U^{(k)})$, inverting the $(p,h)$ table from the updated $(\rho,e)$.
6. Enforce realizability (positive pressure, in-range enthalpy).

!!! note "Primitive vs. preconditioned update"
    With `integration-variables = prim` the residual is applied in primitive variables directly. With `integration-variables = prec` the update is performed in the preconditioned variables, scaling the pseudo-time term by the preconditioning matrix — see [Low-Mach Preconditioning](preconditioning.md).

---

## Stability Conditions

### CFL (convective) limit

For cell $i$ in direction $d$,

$$
\Delta t^{(d)}_{\text{CFL},i} = \frac{\Delta x_d}{|\mathbf u\!\cdot\!\hat{\mathbf e}_d| + a}\times\text{CFL},
$$

with the real-fluid speed of sound $a = a(p,h)$. When preconditioning is active, $a$ is replaced by the **preconditioned** signal speed (built from $U_r$), which raises the allowable step at low Mach number.

### VNN (viscous) limit

For viscous runs the von Neumann limit is

$$
\Delta t^{(d)}_{\text{VNN},i} = \frac{\rho\,(\Delta x_d)^2}{\mu_\ell + \mu_t}\times\text{VNN}.
$$

### Global time step

$$
\Delta t = \min_{i,d}\Bigl[\min\bigl(\Delta t^{(d)}_{\text{CFL},i},\ \Delta t^{(d)}_{\text{VNN},i}\bigr)\Bigr]
$$

for a steady (`time-accurate = false`) run with a single global step; in time-accurate mode the same minimum sets the physical step.

!!! tip "CFL ramp-up"
    `cfl-rise-threshold` linearly ramps the CFL from a low starting value to the target over the given number of iterations, improving robustness when the initial condition is far from steady state.

---

## Implicit Residual Smoothing (IRS)

IRS raises the effective CFL limit by smoothing the residual with a Laplacian operator before the update. In each direction $d$ the smoothed residual satisfies

$$
\mathbf R^\ast_i = \frac{\mathbf R_i + \varepsilon(\mathbf R^\ast_{i-1} + \mathbf R^\ast_{i+1})}{1 + 2\varepsilon},
$$

approximated with a couple of Jacobi sweeps. The stable CFL limit is amplified by roughly $1+2\varepsilon$; typical $\varepsilon$ (`irs-beta`) is 0.1–0.5. Enable it with `irs = .true.`.

---

## Multigrid Acceleration

ARES supports **geometric multigrid** with $2\!:\!1$ coarsening per direction to accelerate convergence to steady state. The number of levels is set by `levels` in `[ARES-Multigrid]`, and the work per level by the `levelN-iter` counts.

### Grid hierarchy

- Each coarse cell aggregates $2^d$ fine cells ($d$ = active spatial dimensions); for 2-D grids ($n_k=1$) coarsening acts only in $i$ and $j$.
- The fine-grid dimensions must be divisible by $2^{\text{levels}-1}$ in each active direction.
- One boundary-condition table per level is required (ATLAS generates the coarse tables).

### Restriction (fine → coarse)

Conservative volume averaging:

$$
\mathbf U_\text{coarse} = \frac{1}{V_\text{coarse}}\sum_{\text{fine}\in\text{coarse}} \mathbf U_\text{fine}\,V_\text{fine}
$$

(8-cell average in 3-D, 4-cell in 2-D). Coarse-grid primitives are recovered by thermo inversion of the averaged conservative state.

### Prolongation (coarse → fine)

The coarse-grid correction is transferred back with cubic interpolation (3-D) or biquadratic interpolation (2-D); each fine cell receives a weighted contribution from its parent and neighbouring coarse cells according to its position within the coarse cell.

### Cycle

A V-cycle (pre-smooth → restrict → coarse solve → prolong → post-smooth) is the default; the driver can recurse for W-cycles. Multigrid can also accelerate each physical step of a transient run.

---

## References

1. C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.* 77 (1988).
2. S. Gottlieb, C.-W. Shu, E. Tadmor, "Strong stability-preserving high-order time discretization methods," *SIAM Rev.* 43 (2001).
3. A. Jameson, "Solution of the Euler equations for two-dimensional transonic flow by a multigrid method," *Appl. Math. Comput.* 13 (1983).
