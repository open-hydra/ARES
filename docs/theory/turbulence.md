# Turbulence Modelling

ARES provides Reynolds-Averaged Navier–Stokes (RANS) closures for turbulent flows, selected at run time by `turbulence-model` in `[ARES-RANS]` (active only when `simulation-type = turbulent`). The families range from a one-equation eddy-viscosity model to a full seven-equation Reynolds-stress model. All models are bound to procedure pointers at setup, so switching model requires no recompilation.

| `turbulence-model` | Type | Notes |
|--------------------|------|-------|
| `SA` | 1-equation eddy viscosity | Baseline Spalart–Allmaras |
| `SA-R` | 1-equation | + rotation correction |
| `SA-RC` | 1-equation | + Spalart–Shur rotation/curvature correction |
| `SAcomp` | 1-equation | + compressibility correction |
| `SA-QCR2000` | 1-eq + algebraic | + Quadratic Constitutive Relation (stress anisotropy) |
| `SA-rough` | 1-equation | + sand-grain wall roughness |
| `SA-rough-QCR2000` | 1-eq + algebraic | roughness + QCR2000 |
| `SST` | 2-equation $k$–$\omega$ | Menter Shear-Stress Transport (2003) |
| `Wilcox2006` | 2-equation $k$–$\omega$ | Wilcox 2006 revision |
| `SSGLRR` | 7-equation RSM | Speziale–Sarkar–Gatski / Launder–Reece–Rodi ($\omega$-2019) |
| `none` | — | Laminar (no model) |

Any of the `-blowcorr` suffix (wall-transpiration correction) and the `Prt-correction` key can be combined with the eligible base models (see [General RANS Architecture](#general-rans-architecture)).

---

## Boussinesq Hypothesis

All eddy-viscosity models assume the Reynolds-stress tensor is proportional to the mean strain rate:

$$
\tau_{ij}^R = 2\,\mu_t\,S_{ij} - \tfrac23\,\rho\,k\,\delta_{ij},
\qquad
S_{ij} = \tfrac12\Bigl(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\Bigr),
$$

where $\mu_t$ is the eddy (turbulent) viscosity and $k$ the turbulence kinetic energy. The isotropic term $-\tfrac23\rho k\,\delta_{ij}$ is retained for the two-equation models and dropped for SA. The SSG-LRR model abandons this hypothesis and transports the Reynolds stresses $\tau_{ij}^R$ directly.

---

## Spalart–Allmaras (SA) — One-Equation Model

### Transport equation

SA transports a single modified eddy viscosity $\tilde\nu$:

$$
\frac{\partial(\rho\tilde\nu)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u\tilde\nu)
= \underbrace{c_{b1}\tilde S\rho\tilde\nu}_{\text{production}}
+ \underbrace{\frac{1}{\sigma}\bigl[\nabla\!\cdot\!((\mu+\rho\tilde\nu)\nabla\tilde\nu) + c_{b2}\rho|\nabla\tilde\nu|^2\bigr]}_{\text{diffusion}}
- \underbrace{c_{w1}f_w\rho\Bigl(\frac{\tilde\nu}{y}\Bigr)^2}_{\text{destruction}}
$$

### Model constants

| Constant | Value | Description |
|:--------:|:-----:|-------------|
| $c_{b1}$ | 0.1355 | Production coefficient |
| $c_{b2}$ | 0.622 | Diffusion coefficient |
| $\sigma$ | 2/3 | Turbulent Schmidt number |
| $\kappa$ | 0.41 | von Kármán constant |
| $c_{w1}$ | $c_{b1}/\kappa^2 + (1+c_{b2})/\sigma$ | Destruction coefficient |
| $c_{w2}$ | 0.3 | Destruction coefficient |
| $c_{w3}$ | 2.0 | Destruction coefficient |
| $c_{v1}$ | 7.1 | Near-wall damping constant |

### Auxiliary functions

$$
\chi = \frac{\tilde\nu}{\nu},\qquad
f_{v1} = \frac{\chi^3}{\chi^3+c_{v1}^3},\qquad
f_{v2} = 1 - \frac{\chi}{1+\chi f_{v1}}
$$

Modified vorticity (with $\Omega$ the vorticity magnitude and $y$ the wall distance):

$$
\tilde S = \Omega + \frac{\tilde\nu}{\kappa^2 y^2}\,f_{v2}
$$

Destruction function:

$$
r = \frac{\tilde\nu}{\kappa^2 y^2\,\tilde S},\qquad
g = r + c_{w2}(r^6-r),\qquad
f_w = g\Bigl(\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\Bigr)^{1/6}
$$

### Eddy viscosity

$$
\mu_t = \rho\,\tilde\nu\,f_{v1}(\chi)
$$

### Wall boundary condition

$$
\tilde\nu_\text{wall} = 0
$$

### SA Variants

- **SA-R (rotation).** Adds $c_\text{rot}\,\min\!\bigl(0,\,\lVert S\rVert-\Omega\bigr)$ to the production ($c_\text{rot}=2$), sensitising the model to system rotation.
- **SA-RC (rotation/curvature).** The Spalart–Shur correction multiplies the production by a factor $f_{r1}(r^\ast,\tilde r)$ built from the strain/vorticity ratio $r^\ast=S/\Omega$ and the material derivative of $S_{ij}$, bounded to $[0,\,1.25]$ — the recommended choice for strongly curved and swirling flows.
- **SAcomp (compressibility).** Paciorri–Sabetta correction scaling the production with the turbulent-stress ratio $S_\tau = \omega\,\tilde\nu f_{v1}/a^2$, for high-speed shear layers.
- **SA-rough (roughness).** Modifies $\chi$ and the wall distance with an offset $d_0 = 0.03\,k_s$ ($k_s$ = sand-grain height) so that a non-zero $\tilde\nu$ enters at a rough wall.
- **SA-QCR2000.** Replaces the Boussinesq stress with the [Quadratic Constitutive Relation](#quadratic-constitutive-relation-qcr2000) to capture stress anisotropy in corner and secondary flows.

---

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear-Stress Transport model blends a $k$–$\omega$ formulation near walls with a $k$–$\varepsilon$-like behaviour in the freestream, through the blending function $F_1$. ARES implements the **Menter 2003** form.

### Transport equations

$$
\frac{\partial(\rho k)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u k)
= P_k - \beta^\ast\rho\omega k + \nabla\!\cdot\!\bigl[(\mu+\sigma_k\mu_t)\nabla k\bigr]
$$

$$
\frac{\partial(\rho\omega)}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u\omega)
= \gamma\frac{\rho P_k}{\mu_t} - \beta\rho\omega^2 + \nabla\!\cdot\!\bigl[(\mu+\sigma_\omega\mu_t)\nabla\omega\bigr]
+ \underbrace{2(1-F_1)\frac{\rho\,\sigma_{\omega2}}{\omega}\,\nabla k\!\cdot\!\nabla\omega}_{\text{cross-diffusion}}
$$

The cross-diffusion term is active only away from walls, where $F_1\to0$.

### Production limiter

$$
P_k = \min\!\bigl(\mu_t S^2,\;10\,\beta^\ast\rho\omega k\bigr)
$$

which prevents unbounded growth of $k$ in stagnation regions ($S=\sqrt{2S_{ij}S_{ij}}$).

### Eddy viscosity

$$
\mu_t = \frac{\rho\,k\,a_1}{\max(a_1\omega,\;S\,F_2)},\qquad a_1 = 0.31
$$

### Model constants

All blended coefficients are computed as $\phi = F_1\phi_1 + (1-F_1)\phi_2$.

| Constant | Set 1 ($\phi_1$) | Set 2 ($\phi_2$) |
|:--------:|:-----:|:-----:|
| $\sigma_k$ | 0.85 | 1.0 |
| $\sigma_\omega$ | 0.5 | 0.856 |
| $\beta$ | 0.075 | 0.0828 |
| $\gamma$ | 5/9 | 0.44 |

Universal: $\beta^\ast = 0.09$, $a_1 = 0.31$, $\kappa = 0.41$.

### Blending functions

$$
F_1 = \tanh\!\bigl(\arg_1^4\bigr),\qquad
\arg_1 = \min\!\left[\max\!\left(\frac{\sqrt{k}}{\beta^\ast\omega y},\;\frac{500\nu}{\omega y^2}\right),\;\frac{4\rho\sigma_{\omega2}k}{CD_{k\omega}\,y^2}\right]
$$

with $CD_{k\omega} = \max\!\bigl(2\rho\sigma_{\omega2}\,\omega^{-1}\nabla k\!\cdot\!\nabla\omega,\;10^{-10}\bigr)$, and

$$
F_2 = \tanh\!\bigl(\arg_2^2\bigr),\qquad
\arg_2 = \max\!\left(\frac{2\sqrt{k}}{\beta^\ast\omega y},\;\frac{500\nu}{\omega y^2}\right)
$$

### Wall boundary conditions

$$
k_\text{wall} = 0,\qquad
\omega_\text{wall} = \frac{60\,\nu}{\beta_1\,y^2} = \frac{800\,\nu}{y^2}\quad(\beta_1=0.075)
$$

the Menter approximate smooth-wall condition for the first cell at distance $y$.

### Energy coupling *(optional)*

When `k-coupling = .true.`, the turbulence kinetic energy is added to the total energy and its production/destruction feeds back into the mean-flow energy equation. It is disabled by default.

---

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

### Model constants

| Constant | Value |
|:--------:|:-----:|
| $\sigma_k$ | 0.6 |
| $\sigma_\omega$ | 0.5 |
| $\beta^\ast$ | 0.09 |
| $\beta_0$ | 0.0708 |
| $\gamma$ | 13/25 |
| $C_\text{lim}$ | 7/8 |
| $\sigma_d$ | 1/8 |

### Eddy viscosity (stress limiter)

$$
\mu_t = \frac{\rho\,k}{\hat\omega},\qquad
\hat\omega = \max\!\left(\omega,\;C_\text{lim}\frac{\sqrt{2\bar S_{ij}\bar S_{ij}}}{\beta^\ast}\right)
$$

The limiter $C_\text{lim}$ caps the eddy viscosity where $\omega$ is small relative to the (deviatoric) strain rate $\bar S_{ij}$.

### Destruction with vortex-stretching correction

The $\omega$-destruction coefficient is modified by the Wilcox stress-limiter function:

$$
\beta = \beta_0\,f_\beta,\qquad
X_\omega = \frac{|W_{ij}\,W_{jk}\,\bar S_{ki}|}{(\beta^\ast\omega)^3},\qquad
f_\beta = \frac{1 + 85\,X_\omega}{1 + 100\,X_\omega}
$$

### Cross-diffusion *(conditional)*

Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega > 0$ (otherwise $\sigma_d=0$):

$$
\text{CD} = \frac{\sigma_d\,\rho}{\omega}\,(\nabla k\!\cdot\!\nabla\omega)
\qquad\text{if } \nabla k\!\cdot\!\nabla\omega > 0
$$

### Production limiter

More permissive than SST:

$$
P_k = \min\!\bigl(\mu_t S^2,\;20\,\beta^\ast\rho\omega k\bigr)
$$

The $\omega$-production uses the $\gamma\,(\omega/k)\,P_k$ form.

---

## SSG-LRR Reynolds-Stress Model

The SSG-LRR-$\omega$ model (NASA RSM-SSGLRR-$\omega$2019) transports the six independent Reynolds stresses $R_{ij}$ plus $\omega$ — seven equations total — removing the Boussinesq assumption entirely. The transport of $R_{ij}$ balances production $P_{ij}$, the pressure–strain redistribution $\Pi_{ij}$, dissipation, and diffusion:

$$
\frac{\partial(\rho R_{ij})}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u R_{ij})
= P_{ij} + \Pi_{ij} - \tfrac23\beta^\ast\rho\omega k\,\delta_{ij} + D_{ij}
$$

The pressure–strain model blends the Speziale–Sarkar–Gatski (SSG) closure away from walls with the Launder–Reece–Rodi (LRR) closure near walls, using the same $\omega$-based length scale and an $F_1$-type blending as SST. The blended LRR/$\omega$ coefficients ($C_1=1.8$, $C_2^{\text{LRR}}=0.52$, …) and SSG/$\varepsilon$ coefficients ($C_1=1.7$, $C_1^\ast=0.9$, $C_2=1.05$, …) follow the NASA TMR reference. RSM closures capture stress anisotropy, secondary flows, and strong streamline curvature that eddy-viscosity models miss, at the cost of seven transported fields and stiffer convergence.

---

## Quadratic Constitutive Relation (QCR2000)

QCR2000 is an **algebraic correction**, not a separate transport model: it replaces the linear Boussinesq stress with a nonlinear (quadratic) constitutive relation,

$$
\tau_{ij}^\text{QCR} = \tau_{ij}^\text{Boussinesq} - c_{nl1}\bigl(O_{ik}\tau_{jk} + O_{jk}\tau_{ik}\bigr),
\qquad
O_{ik} = \frac{2W_{ik}}{\sqrt{W_{mn}W_{mn}}},\quad c_{nl1}=0.3,
$$

which reintroduces stress anisotropy. In ARES it is layered on top of SA (`SA-QCR2000`, `SA-rough-QCR2000`) and is the model of choice for 3-D corner and secondary-flow cases.

---

## General RANS Architecture

### Procedure-pointer architecture

All models are accessed through **procedure pointers** in `Mod_RANS`, bound at setup by the `turbulence-model` key:

| Pointer | Purpose |
|---------|---------|
| `Eddy_Viscosity` | Compute $\mu_t$ from the model variables |
| `RANS_Diffusive_Flux` | Turbulent diffusion of the model variables |
| `Stress_Vector` | Viscous + Reynolds stress on a face |
| `RANS_Set_Wall_Values` | Wall boundary conditions for the model variables |
| `RANS_Set_Blowing_Wall` | Wall-transpiration (blowing/suction) wall values |

This lets the model be switched at run time without recompilation.

### Blowing / wall-transpiration correction

The `-blowcorr` suffix activates a correction (`RANS_Set_Blowing_Wall`) that modifies the wall treatment under wall injection or suction. SA uses its own blowing wall values; the $k$–$\omega$ models (SST, Wilcox 2006) switch between the corrected and uncorrected $\omega$ wall BC accordingly.

### Turbulent transport numbers and corrections

`[ARES-RANS]` also sets the turbulent Prandtl ($\mathrm{Pr}_t$, default 0.90) and Schmidt ($\mathrm{Sc}_t$) numbers used in the [effective conductivity](thermo.md#transport-properties), the optional turbulent-kinetic-energy coupling into the energy equation (`k-coupling`), and a wall-roughness Prandtl correction (`Prt-correction`, requires an `-rough` SA model) that adjusts $\mathrm{Pr}_t$ near rough walls.

---

## References

1. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for aerodynamic flows," AIAA-92-0439, 1992.
2. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial experience with the SST turbulence model," *Turbulence, Heat and Mass Transfer 4*, 2003.
3. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.
4. P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models to rotation and curvature," *Aerosp. Sci. Technol.* 1 (1997).
5. C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain correlation of turbulence," *J. Fluid Mech.* 227 (1991).
6. S. R. Allmaras, F. T. Johnson, P. R. Spalart, "Modifications and clarifications for the implementation of the Spalart–Allmaras turbulence model," ICCFD7-1902, 2012. *(QCR2000)*
