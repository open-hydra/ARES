# Turbulence Modelling

ARES provides Reynolds-Averaged Navier–Stokes (RANS) closures for turbulent flows, selected at run time by `turbulence-model` in `[ARES-RANS]` (active only when `simulation-type = turbulent`). The families range from a one-equation eddy-viscosity model to a full seven-equation Reynolds-stress model.

| `turbulence-model` | Type | Notes |
|--------------------|------|-------|
| `SA` | 1-equation eddy viscosity | Baseline Spalart–Allmaras |
| `SA-R` | 1-equation | + rotation correction |
| `SA-RC` | 1-equation | + Spalart–Shur rotation/curvature correction |
| `SAcomp` | 1-equation | + compressibility correction |
| `SA-QCR2000` | 1-eq + algebraic | + Quadratic Constitutive Relation (stress anisotropy) |
| `SA-rough` | 1-equation | + sand-grain wall roughness |
| `SA-rough-QCR2000` | 1-eq + algebraic | roughness + QCR2000 |
| `SST` | 2-equation $k$–$\omega$ | Menter Shear-Stress Transport |
| `Wilcox2006` | 2-equation $k$–$\omega$ | Wilcox 2006 revision |
| `SSGLRR` | 7-equation RSM | Speziale–Sarkar–Gatski / Launder–Reece–Rodi |
| `none` | — | Laminar (no model) |

---

## Boussinesq Hypothesis

The eddy-viscosity models assume the Reynolds stress is proportional to the mean strain rate:

$$
\tau_{ij}^R = 2\,\mu_t\,S_{ij} - \tfrac23\,\rho\,k\,\delta_{ij},
\qquad
S_{ij} = \tfrac12\Bigl(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\Bigr),
$$

with the isotropic term retained for two-equation models and dropped for SA. The SSG-LRR model abandons this hypothesis and transports the $\tau_{ij}^R$ directly.

---

## Spalart–Allmaras (SA) — One Equation

SA transports a modified eddy viscosity $\tilde\nu$:

$$
\frac{\partial(\rho\tilde\nu)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u\tilde\nu)
= \underbrace{c_{b1}\tilde S\rho\tilde\nu}_{\text{production}}
+ \underbrace{\frac{1}{\sigma}\bigl[\nabla\!\cdot\!((\mu+\rho\tilde\nu)\nabla\tilde\nu) + c_{b2}\rho|\nabla\tilde\nu|^2\bigr]}_{\text{diffusion}}
- \underbrace{c_{w1}f_w\rho\Bigl(\frac{\tilde\nu}{y}\Bigr)^2}_{\text{destruction}}
$$

| Constant | Value |
|:--------:|:-----:|
| $c_{b1}$ | 0.1355 |
| $c_{b2}$ | 0.622 |
| $\sigma$ | 2/3 |
| $\kappa$ | 0.41 |
| $c_{w1}$ | $c_{b1}/\kappa^2 + (1+c_{b2})/\sigma$ |
| $c_{w2}$ | 0.3 |
| $c_{w3}$ | 2.0 |
| $c_{v1}$ | 7.1 |

Auxiliary relations:

$$
\chi = \frac{\tilde\nu}{\nu},\quad
f_{v1} = \frac{\chi^3}{\chi^3+c_{v1}^3},\quad
f_{v2} = 1 - \frac{\chi}{1+\chi f_{v1}},\quad
\tilde S = \Omega + \frac{\tilde\nu}{\kappa^2 y^2}f_{v2},
$$

$$
r = \frac{\tilde\nu}{\kappa^2 y^2\tilde S},\quad
g = r + c_{w2}(r^6-r),\quad
f_w = g\Bigl(\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\Bigr)^{1/6},\quad
\mu_t = \rho\tilde\nu f_{v1},
$$

with wall value $\tilde\nu_\text{wall}=0$, $\Omega$ the vorticity magnitude and $y$ the wall distance.

### SA variants

- **SA-R (rotation).** Adds $c_\text{rot}(\lVert S\rVert-\Omega)$ to the production, sensitising the model to system rotation.
- **SA-RC (rotation/curvature).** The Spalart–Shur correction multiplies the production by a factor $f_{r1}(r^\ast,\tilde r)$ built from the strain/vorticity ratio and the material derivative of $S_{ij}$ — the recommended choice in rotating frames.
- **SAcomp (compressibility).** Paciorri–Sabetta correction scaling production with the turbulent-stress ratio, for high-speed shear layers.
- **SA-rough (roughness).** Modifies $\chi$ and the wall distance with an offset $d_0 = 0.03\,k_s$ ($k_s$ = sand-grain height) so that a non-zero $\tilde\nu$ enters at a rough wall.
- **SA-QCR2000.** Replaces the Boussinesq stress with the Quadratic Constitutive Relation (see [below](#quadratic-constitutive-relation-qcr2000)) to capture stress anisotropy in corner and secondary flows.

---

## Menter SST $k$–$\omega$ — Two Equations

The Shear-Stress Transport model blends $k$–$\omega$ (near walls) with $k$–$\varepsilon$ (freestream):

$$
\frac{\partial(\rho k)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u k)
= P_k - \beta^\ast\rho\omega k + \nabla\!\cdot\!\bigl[(\mu+\mu_t/\sigma_k)\nabla k\bigr]
$$

$$
\frac{\partial(\rho\omega)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u\omega)
= \gamma\frac{\rho P_k}{\mu_t} - \beta\rho\omega^2 + \nabla\!\cdot\!\bigl[(\mu+\mu_t/\sigma_\omega)\nabla\omega\bigr]
+ 2(1-F_1)\frac{\rho}{\sigma_{\omega2}\omega}\nabla k\!\cdot\!\nabla\omega
$$

with the production limiter $P_k = \min(\mu_t S^2, 10\beta^\ast\rho\omega k)$ and eddy viscosity

$$
\mu_t = \frac{\rho k a_1}{\max(a_1\omega, S F_2)},\qquad a_1 = 0.31 .
$$

Blended constants $\phi = F_1\phi_1 + (1-F_1)\phi_2$:

| Constant | Set 1 | Set 2 |
|:--------:|:-----:|:-----:|
| $\sigma_k$ | 0.85 | 1.0 |
| $\sigma_\omega$ | 0.5 | 0.856 |
| $\beta$ | 0.075 | 0.0828 |
| $\gamma$ | 5/9 | 0.44 |

Universal: $\beta^\ast=0.09$, $\kappa=0.41$. Blending functions $F_1 = \tanh(\arg_1^4)$, $F_2=\tanh(\arg_2^2)$ from the wall distance $y$. Wall conditions: $k_\text{wall}=0$, $\omega_\text{wall}=6\nu/(0.075\,y^2)$. ARES implements the **Menter 2003** form.

---

## Wilcox 2006 $k$–$\omega$ — Two Equations

| Constant | Value |
|:--------:|:-----:|
| $\sigma_k$ | 0.6 |
| $\sigma_\omega$ | 0.5 |
| $\beta^\ast$ | 0.09 |
| $\beta_0$ | 0.0708 |
| $\gamma$ | 13/25 |
| $C_\text{lim}$ | 7/8 |
| $\sigma_d$ | 1/8 |

Stress-limited eddy viscosity and vortex-stretching destruction:

$$
\mu_t = \frac{\rho k}{\hat\omega},\quad
\hat\omega = \max\Bigl(\omega, C_\text{lim}\frac{\sqrt{2S_{ij}S_{ij}}}{\beta^\ast}\Bigr),\quad
\beta = \beta_0 f_\beta,\quad
f_\beta = \frac{1+85X_\omega}{1+100X_\omega},
$$

with $X_\omega = |W_{ij}W_{jk}S_{ki}|/(\beta^\ast\omega)^3$. Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega>0$ (coefficient $\sigma_d$), and the production limiter is the more permissive $P_k=\min(\mu_t S^2, 20\beta^\ast\rho\omega k)$.

---

## SSG-LRR Reynolds-Stress Model

The SSG-LRR model transports the six independent Reynolds stresses $R_{ij}$ plus $\omega$ (seven equations total), removing the Boussinesq assumption entirely. The transport of $R_{ij}$ balances production $P_{ij}$, the pressure–strain redistribution $\Pi_{ij}$, dissipation, and diffusion:

$$
\frac{\partial(\rho R_{ij})}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u R_{ij})
= P_{ij} + \Pi_{ij} - \tfrac23\beta^\ast\rho\omega k\,\delta_{ij} + D_{ij},
$$

where the pressure–strain model blends the Speziale–Sarkar–Gatski (SSG) closure away from walls with the Launder–Reece–Rodi (LRR) closure near walls, using the same $\omega$-based length scale and $F_1$-type blending as SST. RSM closures capture stress anisotropy, secondary flows, and strong streamline curvature that eddy-viscosity models miss, at the cost of seven transported fields and stiffer convergence.

---

## Quadratic Constitutive Relation (QCR2000)

QCR2000 is an **algebraic correction**, not a separate transport model: it replaces the linear Boussinesq stress with a nonlinear (quadratic) constitutive relation,

$$
\tau_{ij}^\text{QCR} = \tau_{ij}^\text{Boussinesq} - c_{nl1}\bigl(O_{ik}\tau_{jk} + O_{jk}\tau_{ik}\bigr),
\qquad
O_{ik} = \frac{2W_{ik}}{\sqrt{W_{mn}W_{mn}}},\quad c_{nl1}=0.3,
$$

which reintroduces stress anisotropy. In ARES it is layered on top of SA (`SA-QCR2000`, `SA-rough-QCR2000`) and is the model of choice for the 3-D corner/ablation cases.

---

## General RANS Architecture

All models are accessed through **procedure pointers** in `Mod_RANS`, bound at setup by the `turbulence-model` key:

| Pointer | Purpose |
|---------|---------|
| `Eddy_Viscosity` | Compute $\mu_t$ from the model variables |
| `RANS_Diffusive_Flux` | Turbulent diffusion of the model variables |
| `Stress_Vector` | Viscous + Reynolds stress on a face |
| `RANS_Set_Wall_Values` | Wall boundary conditions for the model variables |

This lets the model be switched at run time without recompilation.

### Turbulent transport numbers and corrections

`[ARES-RANS]` also sets the turbulent Prandtl ($\mathrm{Pr}_t$, default 0.90) and Schmidt ($\mathrm{Sc}_t$) numbers used in the [effective conductivity](thermo.md#transport-properties), an optional turbulent-kinetic-energy coupling into the energy equation (`k-coupling`), and a wall-roughness Prandtl correction (`Prt-correction`) that adjusts $\mathrm{Pr}_t$ near rough walls.

---

## References

1. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for aerodynamic flows," AIAA-92-0439, 1992.
2. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial experience with the SST turbulence model," *Turbulence, Heat and Mass Transfer 4*, 2003.
3. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.
4. P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models to rotation and curvature," *Aerosp. Sci. Technol.* 1 (1997).
5. C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain correlation of turbulence," *J. Fluid Mech.* 227 (1991).
6. S. R. Allmaras, F. T. Johnson, P. R. Spalart, "Modifications and clarifications for the implementation of the Spalart–Allmaras turbulence model," ICCFD7-1902, 2012. *(QCR2000)*
