# Riemann Solvers

At every cell interface the convective flux is computed by approximately (or, with sufficient dissipation, robustly) resolving the local Riemann problem between the reconstructed left and right states. ARES provides **six** numerical-flux solvers, selected by `riemann-solver` in `[ARES-Numerics]`. They share a common interface: given the left/right primitive states $(p,\mathbf u,h)_{L,R}$ and the face normal $\hat{\mathbf n}$, each returns the numerical flux $(\,F_\rho, F_{\rho u}, F_{\rho v}, F_{\rho w}, F_{\rho E}\,)$.

| `riemann-solver` | Internal name | Family | Preconditioned | Notes |
|------------------|---------------|--------|:--------------:|-------|
| `Rusanov` | Local Lax–Friedrichs (Rusanov) | central + max-eigenvalue dissipation | no | Most robust, most diffusive |
| `PLLF` | Preconditioned Local Lax–Friedrichs | central + preconditioned dissipation | **yes** | Low-Mach companion of Rusanov |
| `HLLE` | HLLE | HLL two-wave | no | Positivity-preserving; smears contacts |
| `HLLC` | HLLC (Batten) | HLL three-wave | no | Default; resolves the contact/shear |
| `HLLC Prec` | Preconditioned HLLC | HLL three-wave | **yes** | HLLC with preconditioned wave speeds |
| `HLLC Rotated` | Rotated HLLC Batten / HLLE | hybrid | no | Carbuncle-resistant rotated-Riemann |

!!! tip "Choosing a solver"
    Start from **HLLC** (the default) for most flows — it resolves contact and shear layers that HLLE smears. For **low-Mach** runs combine the preconditioner (`integration-variables = prec`) with **`HLLC Prec`** or **`PLLF`**, whose dissipation is scaled consistently so accuracy does not collapse as $M\to0$. For strong-shock / blunt-body flows prone to the carbuncle, use **`HLLC Rotated`**. Use **`Rusanov`** when robustness matters more than resolution.

---

## Local Lax–Friedrichs (Rusanov)

The simplest and most dissipative choice: a central flux plus dissipation proportional to the largest signal speed,

$$
\mathbf F = \tfrac12\bigl(\mathbf F_L + \mathbf F_R\bigr)
          - \tfrac12\,|\lambda|_{\max}\,(\mathbf U_R - \mathbf U_L),
\qquad
|\lambda|_{\max} = \max\bigl(|u_n|_L + a_L,\ |u_n|_R + a_R\bigr),
$$

where $u_n = \mathbf u\!\cdot\!\hat{\mathbf n}$ and $a$ is the real-fluid speed of sound from the $(p,h)$ table. It is monotone and very robust, but smears contact discontinuities and shear layers.

---

## HLLE (two-wave HLL)

The HLL solver of Harten–Lax–van Leer with the Einfeldt wave-speed estimates models the interface with two acoustic waves enclosing a single averaged intermediate state:

$$
\mathbf F^{\text{HLL}} =
\begin{cases}
\mathbf F_L & S_L \ge 0\\[4pt]
\dfrac{S_R\,\mathbf F_L - S_L\,\mathbf F_R + S_L S_R(\mathbf U_R - \mathbf U_L)}{S_R - S_L} & S_L < 0 < S_R\\[10pt]
\mathbf F_R & S_R \le 0
\end{cases}
$$

with Einfeldt's $S_L = \min(u_{n,L}-a_L,\ \tilde u-\tilde a)$ and $S_R = \max(u_{n,R}+a_R,\ \tilde u+\tilde a)$ from Roe-averaged $\tilde u,\tilde a$. HLLE is positivity-preserving and entropy-satisfying, but because it omits the contact wave it diffuses contact and shear layers — a drawback for boundary layers.

---

## HLLC (three-wave HLL)

HLLC restores the missing middle wave by inserting a contact discontinuity of speed $S^\ast$ between the two acoustic waves, giving two intermediate states $\mathbf U_L^\ast,\mathbf U_R^\ast$:

$$
S^\ast = \frac{p_R - p_L + \rho_L u_{n,L}(S_L - u_{n,L}) - \rho_R u_{n,R}(S_R - u_{n,R})}
              {\rho_L(S_L - u_{n,L}) - \rho_R(S_R - u_{n,R})}
$$

$$
\mathbf F =
\begin{cases}
\mathbf F_L & 0 \le S_L\\
\mathbf F_L + S_L(\mathbf U_L^\ast - \mathbf U_L) & S_L \le 0 \le S^\ast\\
\mathbf F_R + S_R(\mathbf U_R^\ast - \mathbf U_R) & S^\ast \le 0 \le S_R\\
\mathbf F_R & 0 \ge S_R
\end{cases}
$$

ARES uses the Batten et al. formulation of the intermediate states. HLLC resolves contact and shear waves accurately, making it the default for viscous and boundary-layer flows. (In the loading report it is labelled *HLLC Batten*.)

---

## Preconditioned variants — PLLF and HLLC Prec

For low-Mach flows the acoustic eigenvalues $u_n\pm a$ are enormous compared with the convective speed $u_n$, so the upwind dissipation $\propto|\lambda|$ overwhelms the physical fluxes and accuracy collapses. The **preconditioned** solvers replace the acoustic speed by the *preconditioned* eigenvalues built from the reference velocity $U_r$ (see [Low-Mach Preconditioning](preconditioning.md)):

$$
u_n' = \tfrac12 u_n\,(1+\alpha),\qquad
a' = \sqrt{\tfrac14 u_n^2(1-\alpha)^2 + U_r^2},\qquad
\alpha = \frac{U_r^2}{a^2}\ \ (\text{schematically}),
$$

so the dissipation scales with $U_r$ (the convective scale) instead of $a$. `PLLF` applies this to the Rusanov dissipation; `HLLC Prec` applies it to the HLLC wave speeds $S_L, S_R$. They are the solvers to use whenever `integration-variables = prec`.

!!! warning "Match the solver to the update"
    The preconditioned solvers are only consistent when the preconditioned update is active. Using `HLLC Prec` / `PLLF` without `integration-variables = prec` (or vice-versa) mixes scalings and is not recommended.

---

## Rotated HLLC (carbuncle cure)

Strong grid-aligned shocks can trigger the *carbuncle* instability with single-direction upwind solvers. The **rotated** solver of Nishikawa & Kitamura evaluates the Riemann problem in two locally-defined directions: one aligned with the velocity-difference vector (where it uses the diffusive **HLLE** flux to damp the instability) and one orthogonal (where it uses the sharp **HLLC** flux). Blending the two yields a flux that is carbuncle-resistant at shocks yet accurate elsewhere — useful for blunt-body and hypersonic cases.

---

## References

1. A. Harten, P. D. Lax, B. van Leer, "On upstream differencing and Godunov-type schemes for hyperbolic conservation laws," *SIAM Rev.* 25 (1983).
2. B. Einfeldt, "On Godunov-type methods for gas dynamics," *SIAM J. Numer. Anal.* 25 (1988).
3. E. F. Toro, M. Spruce, W. Speares, "Restoration of the contact surface in the HLL-Riemann solver," *Shock Waves* 4 (1994).
4. P. Batten, M. A. Leschziner, U. C. Goldberg, "Average-state Jacobians and implicit methods for compressible viscous and turbulent flows," *J. Comput. Phys.* 137 (1997).
5. E. Turkel, "Preconditioned methods for solving the incompressible and low speed compressible equations," *J. Comput. Phys.* 72 (1987).
6. H. Nishikawa, K. Kitamura, "Very simple, carbuncle-free, boundary-layer-resolving, rotated-hybrid Riemann solvers," *J. Comput. Phys.* 227 (2008).
