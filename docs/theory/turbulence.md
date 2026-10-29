# Turbulence Modelling

ARES provides Reynolds-Averaged Navier–Stokes (RANS) closures for turbulent flows, selected at run time by `turbulence` in `[ARES-Physics]` (with `equations = navier-stokes`; `none`, the default, is a laminar run), with the options of the model in `[ARES-Turbulence]`: the same keys and sections as MOSE. The families range from a one-equation eddy-viscosity model to a full seven-equation Reynolds-stress model. All models are bound to procedure pointers at setup, so switching model requires no recompilation.

| `turbulence` | Type | Notes |
|--------------|------|-------|
| `SA` | 1-equation eddy viscosity | Baseline Spalart–Allmaras |
| `SA-R` | 1-equation | + rotation correction |
| `SA-RC` | 1-equation | + Spalart–Shur rotation/curvature correction |
| `SAcomp` | 1-equation | + compressibility correction |
| `SA-QCR2000` | 1-eq + algebraic | + Quadratic Constitutive Relation (stress anisotropy) |
| `SA-rough` | 1-equation | + sand-grain wall roughness |
| `SA-rough-QCR2000` | 1-eq + algebraic | roughness + QCR2000 |
| `SST` | 2-equation $k$–$\omega$ | Menter Shear-Stress Transport (2003) |
| `SST-RC` | 2-equation $k$–$\omega$ | + Spalart–Shur rotation/curvature correction (Smirnov–Menter) |
| `Wilcox2006` | 2-equation $k$–$\omega$ | Wilcox 2006 revision |
| `SSGLRR` | 7-equation RSM | Speziale–Sarkar–Gatski / Launder–Reece–Rodi (w2012; `-SD`, `-w2019` variants) |
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
- **SA-RC (rotation/curvature).** The Spalart–Shur correction multiplies the production by a factor $f_{r1}(r^\ast,\tilde r)$, built from the strain/vorticity ratio $r^\ast=S/\Omega$ and from the material derivative of $S_{ij}$. It is the recommended choice for strongly curved and swirling flows. For SA-RC, $f_{r1}$ is not bounded and can be negative (NASA TMR, SA-RC); SST-RC bounds its own to $[0,\,1.25]$.

    ARES computes the factor $2r^\ast/(1 + r^\ast)$ of $f_{r1}$ as $2S/(S + \Omega)$: where the vorticity vanishes it takes its limit, 2. Where strain and vorticity are both below $10^{-30}\ \mathrm{s^{-1}}$, in a fluid at rest or a uniform flow, $r^\ast$ is not defined: the velocity gradient there is the $10^{-40}$ that ARES adds to the velocity differences, times the metric of the cell. ARES then takes the value of simple shear, 1, with $\tilde r = 0$. The correction does not act there, and the production is that of SA.
- **SAcomp (compressibility).** Paciorri–Sabetta correction scaling the production with the turbulent-stress ratio $S_\tau = \omega\,\tilde\nu f_{v1}/a^2$, for high-speed shear layers.
- **SA-rough (roughness).** Equivalent sand-grain roughness via the Aupoix–Spalart (“Boeing”) extension — see [Rough walls](#rough-walls-sa-rough) below.
- **SA-QCR2000.** Replaces the Boussinesq stress with the [Quadratic Constitutive Relation](#quadratic-constitutive-relation-qcr2000) to capture stress anisotropy in corner and secondary flows.

### Rough walls (SA-rough)

The `-rough` suffix activates the Aupoix–Spalart (“Boeing”) equivalent sand-grain extension. The geometric wall distance $y$ is offset by a roughness length, and the viscosity ratio gains a roughness term:

$$
d = y + d_0,\qquad d_0 = 0.03\,k_s,\qquad
\chi = \frac{\tilde\nu}{\nu} + c_{R1}\frac{k_s}{d},\qquad c_{R1} = 0.5
$$

with $d$ replacing $y$ in $\bar S$, in the destruction function $r$, in the destruction term itself and in the point-implicit factor. The wall condition on $\tilde\nu$ becomes $(\partial\tilde\nu/\partial n)_w = \tilde\nu_w/(0.03\,k_s)$, discretised in `Spalart_Set_Wall_Values` as

$$
\tilde\mu_w = \frac{2\,\tilde\mu_c\,m\,(0.03\,k_s)}{1 + 2\,m\,(0.03\,k_s)},\qquad m = \mathbf{M}_{(\text{dir})}\!\cdot\!\mathbf{n}
$$

Note also that $f_{v2}$ is evaluated from the **molecular** ratio rather than the augmented $\chi$:

$$
f_{v2} = 1 - \frac{\tilde\nu}{\nu + \tilde\nu\,f_{v1}(\chi)}
$$

Smooth-wall this is algebraically the same as the $1-\chi/(1+\chi f_{v1})$ form quoted above, since $\chi=\tilde\nu/\nu$; with roughness the two differ, because $\chi$ carries the extra $c_{R1}k_s/d$ term while the numerator must not. ARES uses this form unconditionally, so it is correct in both regimes.

#### Where $\chi$ is evaluated

$\chi$ is a **field function of the local coordinate** $y$, so `Spalart_Eddy_Viscosity` must receive the wall distance *of the point being evaluated* — not a fixed reference length:

| Location | $\tilde\nu$ | $y$ |
|---|---|---|
| Cell centre | cell value | $y_n(i,j,k)$ |
| Interior / connection face | interface value | $\tfrac12\bigl(y_{n,L}+y_{n,R}\bigr)$ |
| **Wall face** | **wall value $\tilde\mu_w$** | **$0$** |

On a wall face the flux itself is formed there, so $y=0$ and the roughness term collapses to a constant independent of $k_s$:

$$
c_{R1}\frac{k_s}{0 + 0.03\,k_s} = \frac{0.5}{0.03} \simeq 16.67
$$

The expression is regular at $y=0$ (the denominator tends to $0.03\,k_s>0$), so this is a genuine field value, not a limit that has to be guarded.

!!! note "Eddy viscosity does not vanish at a rough wall"
    On a **smooth** wall $\tilde\nu_w = 0$, hence $\chi=0$, $f_{v1}(0)=0$ and $\mu_t = 0$ — exact, since all Reynolds stresses vanish under no-slip and $\tau_w$ is purely molecular. On a **rough** wall the Boeing extension deliberately leaves $\tilde\nu_w \neq 0$, so $\mu_t$ at the wall is finite: that non-zero wall eddy viscosity *is* the mechanism producing the $\Delta U^+$ log-layer shift. Evaluating $\chi$ with the first-cell distance instead of $y=0$ therefore under-predicts the wall eddy viscosity, by an amount that grows as the near-wall cell coarsens.

Roughness is implemented for the Spalart–Allmaras family only; with a $k$–$\omega$ model the `-rough` suffix raises a warning and the smooth-wall $\omega$ condition is used.

---

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear-Stress Transport model blends a $k$–$\omega$ formulation near walls with a $k$–$\varepsilon$-like behaviour in the freestream, through the blending function $F_1$.

ARES implements the **SST-2003** version of the model (Menter, Kuntz & Langtry [2]; the original model is Menter 1994 [7]) as specified by the NASA Turbulence Modeling Resource [8]:

- eddy viscosity built from the strain invariant $S$;
- production limiter $10\,\beta^\ast\rho\,\omega k$ in both equations;
- the 2003 constants.

Two options select the remaining details, and together they span TMR's naming of the variants:

| `sst-production` | Production $P$ | `k-energy-coupling = false` (default) | `k-energy-coupling = true` |
|---|---|---|---|
| `incompressible` (default) | $\mu_t S^2$ | SST-2003m | SST-2003s |
| `compressible` | $\tau_{ij}\,\partial u_i/\partial x_j$ (exact) | SST-2003e | SST-2003 |

The suffix records the approximations [8]:

- Without the [energy coupling](#energy-coupling-k-energy-coupling), the $-\tfrac23\rho k\,\delta_{ij}$ part of $\tau_{ij}$ is left out of the momentum and energy equations ("m", "e").
- "m" and "s" approximate the production by $\mu_t S^2$. TMR describes this as exact for incompressible flow and "typically considered a very good approximation, except perhaps for very high Mach number flows".
- With `compressible`, the limiter acts on the deviatoric part of $P$ only, as in OpenFOAM, whereas TMR's SST-2003 limits all of it ([below](#compressible-production-sst-production)).

`turbulence = SST-RC` adds the [rotation/curvature correction](#rotationcurvature-correction-sst-rc).

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

which prevents unbounded growth of $k$ in stagnation regions ($S=\sqrt{2S_{ij}S_{ij}}$). The $\omega$ production, $\gamma\rho P_k/\mu_t$, uses the limited $P_k$, as corrected on the TMR page [8]; the 2003 paper has the unlimited value there.

### Compressible production (`sst-production`)

The production above is the SST-2003m form $\mu_t S^2$, which drops the dilatation terms of $P = \tau_{ij}\,\partial u_i/\partial x_j$. With `sst-production = compressible` (default `incompressible`) ARES keeps them, with $\tau_{ij} = \mu_t\bigl(2S_{ij} - \tfrac23\,\nabla\!\cdot\!\mathbf u\,\delta_{ij}\bigr) - \tfrac23\,\rho k\,\delta_{ij}$:

$$
P_k = \min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf u)^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\rho k\,\nabla\!\cdot\!\mathbf u,
$$

$$
P_\omega = \frac{\gamma\,\rho}{\mu_t}\,\min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf u)^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\gamma\,\rho\,\omega\,\nabla\!\cdot\!\mathbf u.
$$

This is the exact SST-2003 production. TMR calls it SST-2003e, since the $-\tfrac23\rho k$ term still stays out of the momentum and energy equations; with the [energy coupling](#energy-coupling-k-energy-coupling) it is SST-2003. The dilatation terms are arranged as in OpenFOAM's `kOmegaSST` [9]:

- the deviatoric part, $\mu_t[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf u)^2] = 2\mu_t S^{d}_{ij}S^{d}_{ij} \ge 0$, stays under the limiter;
- the dilatation terms are added outside it. TMR's SST-2003 limits the whole of $P$ instead.

In attached, nearly incompressible boundary layers the two forms coincide. Through a compression ($\nabla\!\cdot\!\mathbf u<0$: shock feet, reattachment) the dilatation terms amplify $k$ more than $\omega$ ($\gamma<1$), so $\mu_t$ grows.

### Rotation/curvature correction (SST-RC)

`turbulence = SST-RC` sensitizes the model to streamline curvature and system rotation. It uses the Spalart–Shur correction as adapted to SST by Smirnov and Menter [10]: TMR SST-RC [8], here on SST-2003, so SST-2003RC. The production $P$ of both equations is multiplied by

$$
f_{r1} = \max\!\bigl[\min(f_\text{rot},\,1.25),\,0\bigr], \qquad
f_\text{rot} = (1 + c_{r1})\,\frac{2\,r^\ast}{1 + r^\ast}
  \bigl[1 - c_{r3}\arctan(c_{r2}\,\tilde r)\bigr] - c_{r1},
$$

$$
r^\ast = \frac{S}{\Omega}, \qquad
\tilde r = \frac{2\,\Omega_{ik}\,S_{jk}}{\Omega\,D^3}\,\frac{DS_{ij}}{Dt}, \qquad
D^2 = \max(S^2,\;0.09\,\omega^2),
$$

with $c_{r1} = 1$, $c_{r2} = 2$, $c_{r3} = 1$, $S = \sqrt{2S_{ij}S_{ij}}$ and $\Omega = \sqrt{2\Omega_{ij}\Omega_{ij}}$. $f_{r1}$ multiplies $P$ before the production limiter, so it reaches the $\omega$ equation through the limited $P$. In simple shear $r^\ast = 1$ and $\tilde r = 0$, so $f_{r1} = 1$ and the model reduces to SST.

Implementation notes:

- The factor $2r^\ast/(1 + r^\ast)$ and the numerator of $\tilde r$ are shared with SA-RC (`Lib_SpalartShur`), and so is their value where strain and vorticity vanish ([SA variants](#sa-variants)). In a fluid at rest or a uniform flow, $f_{r1} = 1$.
- Where only the vorticity vanishes, the factor takes its limit 2 and $\tilde r = 0$, so $f_{r1} = 1.25$. The production there is negligible.
- $DS_{ij}/Dt$ is its steady convective part, $\mathbf u\!\cdot\!\nabla S_{ij}$.

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

Both arguments diverge as $y\to0$, so $F_2\to\tanh(\infty)=1$ at a wall face. `compute_F2` returns that limit directly for $y\le0$, which keeps the blending function well defined when the wall boundary condition evaluates transport properties on the face itself. The value is immaterial to the result there — $k_\text{wall}=0$ forces $\mu_t=0$ regardless — but it avoids a division by zero.

### Wall boundary conditions

$$
k_\text{wall} = 0,\qquad
\omega_\text{wall} = \frac{60\,\nu}{\beta_1\,y^2} = \frac{800\,\nu}{y^2}\quad(\beta_1=0.075)
$$

the Menter approximate smooth-wall condition for the first cell at distance $y$. This is the default, `omega-wall-bc = practical` in `[ARES-Turbulence]`. On a wall-resolved mesh ($y^+ < 1$), `omega-wall-bc = asymptotic` uses instead the exact $y\to0$ limit of the $\omega$ equation,

$$
\omega_\text{wall} = \frac{6\,\nu}{\beta_1\,y^2} = \frac{80\,\nu}{y^2},
$$

as in MOSE. The key sets the wall condition of the Wilcox 2006 and SSG-LRR models as well. Separately, `sst-asymptotic = .true.` (SST only, ARES only) imposes the same asymptotic value on the first `sst-asymptotic-cells` cells off each viscous wall, removing them from the time integration of the $\omega$ equation.

### Energy coupling (`k-energy-coupling`)

The compressible RANS equations [11] carry the turbulent kinetic energy in two places:

- the total energy, $E = e + \tfrac12 u_iu_i + k$;
- the Boussinesq Reynolds stress, $\tau_{ij} = 2\mu_t\bigl(S_{ij} - \tfrac13\nabla\!\cdot\!\mathbf u\,\delta_{ij}\bigr) - \tfrac23\rho k\,\delta_{ij}$.

By default ARES drops both, as in the TMR "m"/"e" variants: its energy is that of the mean flow, $E_m = e + \tfrac12 u_iu_i$, and the $-\tfrac23\rho k$ term is not in the stress. With `k-energy-coupling = true` (two-equation models: SST and Wilcox 2006) ARES keeps both:

- The stress of the momentum and energy equations includes $-\tfrac23\rho k\,\delta_{ij}$ (`Stress_Vector_2eq`), so it also enters the work term $\nabla\!\cdot\!(\boldsymbol\tau\cdot\mathbf u)$ of the energy equation.
- The mean-flow energy gives up what the turbulence receives:

$$
\frac{\partial(\rho E_m)}{\partial t} + \dots = \dots - \bigl(P_k - \beta^\ast\rho\,\omega\,k\bigr).
$$

The sink is the same source as in the $k$ equation, including its point-implicit factor when `point-implicit = .true.`. So $E_m + k$ is conserved.

This is the total-energy form of TMR [11] and Wilcox [3], written for $E_m$. Subtracting the $k$ equation from the equation for $E = E_m + k$ leaves exactly this sink, and the diffusion of $k$ cancels between the two equations. The thermodynamic state follows from the mean-flow energy in either case: $k$ is not part of the internal energy.

TMR notes that the $k$ terms are sometimes dropped "for non-supersonic speed flows". They matter where $k$ is not small against the kinetic energy of the mean flow.

---

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

The Wilcox (2006) model [3, 12] uses the exact production $P = \tau_{ij}\,\partial u_i/\partial x_j$ and the TMR $k$-equation limiter (below). In TMR's naming this is Wilcox2006-klim, in its "e" form (Wilcox2006-klim-e) unless [`k-energy-coupling`](#energy-coupling-k-energy-coupling) is on [8].

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

More permissive than SST, and in the $k$ equation only, as in TMR's Wilcox2006-klim:

$$
P_k = \min\!\bigl(\tau_{ij}\,\partial u_i/\partial x_j,\;20\,\beta^\ast\rho\omega k\bigr)
$$

The $\omega$ production uses the unlimited $P$: $\gamma\,(\omega/k)\,P$.

---

## SSG-LRR Reynolds-Stress Model

The SSG-LRR-$\omega$ model (NASA TMR SSGLRR-RSM-w2012) transports the six independent Reynolds stresses $R_{ij}$ plus $\omega$ — seven equations total — removing the Boussinesq assumption entirely. The transport of $R_{ij}$ balances production $P_{ij}$, the pressure–strain redistribution $\Pi_{ij}$, dissipation, and diffusion:

$$
\frac{\partial(\rho R_{ij})}{\partial t} + \nabla\!\cdot\!(\rho\mathbf u R_{ij})
= P_{ij} + \Pi_{ij} - \tfrac23\beta^\ast\rho\omega k\,\delta_{ij} + D_{ij}
$$

The pressure–strain model blends the Speziale–Sarkar–Gatski (SSG) closure away from walls with the Launder–Reece–Rodi (LRR) closure near walls, using the same $\omega$-based length scale and an $F_1$-type blending as SST. The blended LRR/$\omega$ coefficients ($C_1=1.8$, $C_2^{\text{LRR}}=0.52$, …) and SSG/$\varepsilon$ coefficients ($C_1=1.7$, $C_1^\ast=0.9$, $C_2=1.05$, …) follow the NASA TMR reference. RSM closures capture stress anisotropy, secondary flows, and strong streamline curvature that eddy-viscosity models miss, at the cost of seven transported fields and stiffer convergence.

| Model | TMR name | Difference |
|-------|----------|------------|
| `SSGLRR` | SSGLRR-RSM-w2012 | Generalized gradient diffusion, $D\,\rho R_{kl}/(C_\mu\omega)\,\partial R_{ij}/\partial x_l$ |
| `SSGLRR-SD` | SSGLRR-RSM-w2012-SD | Simple diffusion, $(D/C_\mu)\,\mu_t\,\partial R_{ij}/\partial x_k$ |
| `SSGLRR-w2019` | SSGLRR-RSM-w2019 | w2012 + length-scale correction of the $\omega$ equation |
| `SSGLRR-w2019-SD` | SSGLRR-RSM-w2019-SD | w2012-SD + length-scale correction |

The w2019 length-scale correction (Eisfeld & Rumsey, 2020) multiplies the $\omega$ destruction $\beta_\omega\rho\omega^2$ by $1-F^{(LSC)}$:

$$
F^{(LSC)} = \tfrac12\left\{1+\tanh\left[A(\chi-\chi_T)\right]\right\},\qquad
\chi = \max\left[\left(\frac{\mathcal L_t}{\ell_{log}}-1\right)\left(\frac{\mathcal L_t}{\ell_{log}}\right)^2,\,0\right],\qquad
\frac{\mathcal L_t}{\ell_{log}} = \frac{\sqrt{k}}{C_\mu^{1/4}\,\kappa\,\omega\,d}
$$

with $A=31$, $\chi_T=1$, $\kappa=0.41$. It vanishes in the log layer ($\mathcal L_t=\ell_{log}$) and is active only near stagnation and reattachment points, where it removes the non-physical "backbending" of the original model.

At walls $R_{ij}=0$ and $\omega$ follows `omega-wall-bc` (the default `practical` value is the TMR one). On symmetry planes, which must be $x$-, $y$- or $z$-planes, the shear stresses with one index normal to the plane vanish (e.g. $R_{12}=R_{13}=0$ on an $x$-plane) and the other components have zero gradient; the momentum flux through the wall and symmetry faces uses these face values.

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

All models are accessed through **procedure pointers** in `Mod_RANS`, bound at setup by the `turbulence` key:

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

`[ARES-Turbulence]` also sets the turbulent Prandtl ($\mathrm{Pr}_t$, default 0.90) and Schmidt ($\mathrm{Sc}_t$) numbers used in the [effective conductivity](thermo.md#transport-properties), the optional turbulent-kinetic-energy coupling into the energy equation (`k-energy-coupling`), and a wall-roughness Prandtl correction (`Prt-correction`, requires an `-rough` SA model) that adjusts $\mathrm{Pr}_t$ near rough walls.

---

## References

1. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for aerodynamic flows," AIAA-92-0439, 1992.
2. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial experience with the SST turbulence model," *Turbulence, Heat and Mass Transfer 4*, 2003.
3. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.
4. P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models to rotation and curvature," *Aerosp. Sci. Technol.* 1 (1997).
5. C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain correlation of turbulence," *J. Fluid Mech.* 227 (1991).
6. S. R. Allmaras, F. T. Johnson, P. R. Spalart, "Modifications and clarifications for the implementation of the Spalart–Allmaras turbulence model," ICCFD7-1902, 2012. *(QCR2000)*
7. F. R. Menter, "Two-equation eddy-viscosity turbulence models for engineering applications," *AIAA J.* 32(8), pp. 1598–1605, 1994, doi:10.2514/3.12149.
8. NASA Langley Turbulence Modeling Resource, "The Menter Shear Stress Transport turbulence model," <https://tmbwg.github.io/turbmodels/sst.html> (formerly turbmodels.larc.nasa.gov/sst.html): definitions and naming of the SST variants (SST-2003, SST-2003m/s/e, SST-RC) and of Wilcox2006-klim.
9. OpenFOAM Foundation, OpenFOAM-10, `kOmegaSST` model (`src/MomentumTransportModels/momentumTransportModels/Base/kOmegaSST/kOmegaSSTBase.C`).
10. P. E. Smirnov, F. R. Menter, "Sensitization of the SST turbulence model to rotation and curvature by applying the Spalart–Shur correction term," *J. Turbomach.* 131(4), 041010, 2009, doi:10.1115/1.3070573.
11. NASA Langley Turbulence Modeling Resource, "Implementing turbulence models into the compressible RANS equations," <https://tmbwg.github.io/turbmodels/implementrans.html>.
12. D. C. Wilcox, "Formulation of the k-ω turbulence model revisited," *AIAA J.* 46(11), pp. 2823–2838, 2008, doi:10.2514/1.36541.
