# Turbulence Modelling

ARES provides Reynolds-Averaged Navier–Stokes (RANS) closures for turbulent flows, selected at run time by `turbulence` in `[ARES-Physics]` (with `equations = navier-stokes`; `none`, the default, is a laminar run), with the options of the model in `[ARES-Turbulence]`: the same keys and sections as MOSE. The families range from a one-equation eddy-viscosity model to a full seven-equation Reynolds-stress model. All models are bound to procedure pointers at setup, so switching model requires no recompilation.

| `turbulence` | Type | Notes |
|--------------|------|-------|
| `SA` | 1-equation eddy viscosity | Baseline Spalart–Allmaras |
| `SA-R` | 1-equation | + rotation correction |
| `SA-RC` | 1-equation | + Spalart–Shur rotation/curvature correction |
| `SAcomp` | 1-equation | + compressibility correction (stops at set-up: its result depends on the units) |
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
+ \underbrace{\frac{1}{\sigma}\bigl[\nabla\!\cdot\!((\mu+\rho\tilde\nu)\nabla\tilde\nu) + c_{b2}\rho|\nabla\tilde\nu|^2
  - (\nu+\tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu\bigr]}_{\text{diffusion}}
- \underbrace{c_{w1}f_w\rho\Bigl(\frac{\tilde\nu}{y}\Bigr)^2}_{\text{destruction}}
$$

This is the conservation form of the original equation, which has no density and holds for
compressible flows too: it is $\rho$ times the original equation plus $\tilde\nu$ times the
continuity equation (S. R. Allmaras, F. T. Johnson, P. R. Spalart, *Modifications and
Clarifications for the Implementation of the Spalart–Allmaras Turbulence Model*, ICCFD7-1902,
2012, eq. 9; NASA Turbulence Modeling Resource, Spalart–Allmaras page). Since
$\rho\,\nabla\cdot\bigl((\nu+\tilde\nu)\nabla\tilde\nu\bigr) = \nabla\cdot\bigl(\rho(\nu+\tilde\nu)\nabla\tilde\nu\bigr)
- (\nu+\tilde\nu)\nabla\rho\cdot\nabla\tilde\nu$, the last diffusion term is needed wherever density
and $\tilde\nu$ vary together (boundary layers at high Mach number, shocks, heated walls): without it,
as in implementations that put the density inside the derivatives, the predictions of supersonic flows
change (NASA Turbulence Modeling Resource, Note 4). In `SA_Source_Blk` the density gradient takes the
same central differences as the gradient of $\tilde\nu$.

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

The ghost cell behind a viscous wall takes $2\tilde\nu_\text{wall} - \tilde\nu_\text{cell}$
([Wall ghost cells of the turbulence variables](#wall-ghost-cells-of-the-turbulence-variables)). Negative values are
clipped: after every stage ARES keeps $\rho\tilde\nu \ge 10^{-10}$ (the realizability floor) and counts the cell
updates in which a stage made it negative ([Time-accurate runs](#time-accurate-runs-explicit-sources)).

### SA Variants

- **SA-R (rotation).** Adds $c_\text{rot}\,\min\!\bigl(0,\,\lVert S\rVert-\Omega\bigr)$ to the production ($c_\text{rot}=2$), sensitising the model to system rotation.
- **SA-RC (rotation/curvature).** The Spalart–Shur correction multiplies the production by a factor $f_{r1}(r^\ast,\tilde r)$, built from the strain/vorticity ratio $r^\ast=S/\Omega$ and from the material derivative of $S_{ij}$. It is the recommended choice for strongly curved and swirling flows. For SA-RC, $f_{r1}$ is not bounded and can be negative (NASA TMR, SA-RC); SST-RC bounds its own to $[0,\,1.25]$.

    ARES computes the factor $2r^\ast/(1 + r^\ast)$ of $f_{r1}$ as $2S/(S + \Omega)$: where the vorticity vanishes it takes its limit, 2. Where strain and vorticity are both below $10^{-30}\ \mathrm{s^{-1}}$, in a fluid at rest or a uniform flow, $r^\ast$ is not defined: the velocity gradient there is the $10^{-40}$ that ARES adds to the velocity differences, times the metric of the cell. ARES then takes the value of simple shear, 1, with $\tilde r = 0$. The correction does not act there, and the production is that of SA.

    ARES evaluates the Lagrangian derivative $DS_{ij}/Dt = \partial S_{ij}/\partial t + u_k\,\partial S_{ij}/\partial x_k$ with its convective part only. Without the time term the correction holds for steady flows only: «Without the time term, the model is not correct in a time-dependent flowfield» (NASA Turbulence Modeling Resource, Spalart–Allmaras page, SA-RC). Time-accurate runs with a `-RC` correction (SA-RC, SST-RC) print a warning at start-up.
- **SAcomp (compressibility).** Paciorri–Sabetta correction scaling the production with the turbulent-stress ratio $S_\tau = \omega\,\tilde\nu f_{v1}/a^2$, for high-speed shear layers.

    !!! warning "Stopped at set-up"
        A run with `turbulence = SAcomp` stops at set-up. The correction takes
        $c$ from $c^2 f(c) = 25\,\Omega\,\tilde\nu f_{v1}/a$, with
        $f(c) = 0.44/(1 + 14c^5) + 0.56$, and multiplies the production by
        $f_1(c)\,f(c)$, with $f_1(c) = 0.6/(1 + 9c^6) + 0.4$: the left side is a
        number and the right side a speed (the code divides by $a$, not by $a^2$
        as written above), so the factor depends on the units. For
        $\Omega = 10^3$ 1/s, $\nu_t = 10^{-3}$ m²/s and $a = 340$ m/s it is 0.99
        with lengths in metres and 0.22 with lengths in millimetres; and since
        $a$ only divides, at $\Omega = 10^5$ 1/s it cuts the production by 78 %
        whatever the Mach number. Until the correction is checked against its
        source (the article the code cites, doi:10.2514/2.3967) it is not
        available: use `SA` or another model.
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

as in MOSE. The key sets the wall condition of the Wilcox 2006 and SSG-LRR models as well, the asymptotic one with the $\beta$ of the $\omega$ equation of each model at the wall: $6\nu/(\beta_0 y^2) = 84.7\,\nu/y^2$ for Wilcox 2006 (see its section below), $80\,\nu/y^2$ for SSG-LRR, whose $\omega$ equation takes $\beta = 0.075$ there, as SST. Separately, `sst-asymptotic = .true.` (SST only, ARES only) imposes the same asymptotic value on the first `sst-asymptotic-cells` cells off each viscous wall, removing them from the time integration of the $\omega$ equation.

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

The sink is the same source as in the $k$ equation, so $E_m + k$ is conserved. In steady runs with `point-implicit = .true.` the [point-implicit](#steady-runs-point-implicit-factor) factor treats the energy and $k$ equations as one system: the two increments add up to that of the fluxes of $E_m + k$ at every stage, and the steady state does not depend on the time step.

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
X_\omega = \frac{|W_{ij}\,W_{jk}\,\hat S_{ki}|}{(\beta^\ast\omega)^3},\qquad
f_\beta = \frac{1 + 85\,X_\omega}{1 + 100\,X_\omega}
$$

with $\hat S_{ki} = S_{ki} - \tfrac{1}{2}\,(\partial u_m/\partial x_m)\,\delta_{ki}$ (NASA Turbulence Modeling
Resource, Wilcox 2006 page): off the diagonal $\hat S = S$. On that page the formula of $\chi_\omega$ reads
$\Omega_{ij}\Omega_{ij}\hat S_{ki}$, whose indices do not close; ARES takes $W_{ij}W_{jk}\hat S_{ki}$, the form of the Wilcox
1998 model on the same page ($\Omega_{ij}\Omega_{jk}S_{ki}$) and the only one for which the statement of the page that
$\chi_\omega$ is zero in 2-D holds. In 2-D flows the product $W_{ij}W_{jk}\hat S_{ki}$ is zero and $f_\beta = 1$.

### Wall boundary conditions

$k_\text{wall} = 0$. With `omega-wall-bc = asymptotic` the wall takes the near-wall limit of the $\omega$ equation of the
model, $\omega_\text{wall} = 6\nu/(\beta_0 y^2)$ with $\beta_0 = 0.0708$ (NASA Turbulence Modeling Resource, Wilcox 2006
page), that is $84.7\,\nu/y^2$, $y$ the distance of the centre of the wall cell from the wall. With `practical`, the default,
it takes Menter's condition of the [SST section](#wall-boundary-conditions), $800\,\nu/y^2$ with the $\beta_1$ of SST,
as the same page advises for this model.

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
| `RANS_Extrapolate_Wall` | Ghost-cell value of the model variables at walls, $2Z_\text{wall} - Z_\text{cell}$ |
| `RANS_Set_Blowing_Wall` | Wall-transpiration (blowing/suction) wall values |
| `RANS_Stiffness_Bound` | Bound of the stiffness of the destruction terms from the cell alone, for the time step of [time-accurate runs](#time-accurate-runs-explicit-sources) |

This lets the model be switched at run time without recompilation.

### Wall ghost cells of the turbulence variables

Behind a viscous wall (`301`, `302`) the ghost cell of each turbulence variable takes the value that puts the wall
value on the face,

$$
Z_\text{ghost} = 2\,Z_\text{wall} - Z_\text{cell}, \qquad Z = \tilde\nu,\ k,\ \omega,\ R_{ij} ,
$$

with the wall values of the models above (`RANS_Extrapolate_Wall`). It is the rule of CFL3D, «the actual wall boundary
conditions for the turbulence quantities (unlike the primitive variables) are applied at ghost cells» (S. L. Krist,
R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual (Version 5.0)*, NASA/TM-1998-208444, 1998, appendix H, eq. H-186).
The wall flux routines set these values at every stage; the source terms read them in the gradients of the wall cells:
the $c_{b2}$ term of SA and its term $-(\nu+\tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu/\sigma$, $F_1$ and the
cross-diffusion of SST, the cross-diffusion of Wilcox 2006. With $\tilde\nu = a\,y$ next to a wall the $c_{b2}$ term of
the wall cell sees $\partial\tilde\nu/\partial y = a$; a ghost cell that copied the cell would give $a/2$. The other
primitive variables are copied into the wall ghost cells, with the velocity reversed (`Ghost_Wall_Extrapolation`); the
fluxes through the first interior faces use the extrapolation of `Fill_Ghost_Cell`. The walls coupled to another solver
(`103`, `104`) have no wall flux routine of their own: their ghost cells copy the cell, turbulence variables included.

### Wall distance

The distance $y$ of the models (SA, the SST blending functions, the wall
value of $\omega$) is the minimum distance from the cell centre, the mean of
its eight nodes, to the wall surface, as the NASA Turbulence Modeling
Resource requires [13]. It is not the distance to the nearest wall node
or face centre, which differs from it wherever the grid lines are not normal
to the wall; next to a convex corner or at the end of a wall it is the
distance to the corner or the end.

The distance to the centre of the nearest wall face is wrong by a fixed
fraction in the first cell of a skewed mesh. With the grid lines at an angle
$\theta$ to the wall normal, the centre of the first cell, at the height
$d_1$, lies $d_1\tan\theta$ along the wall from the centre of its wall face,
the nearest face centre while $d_1\tan\theta$ is at most half the face, at
the distance $d_1/\cos\theta$. The relative error,
$1/\cos\theta - 1$, does not fall with refinement: 1.5 % at 10°, 6.4 % at 20°,
15.5 % at 30°, 41 % at 45°. It enters squared the destruction of SA,
$c_{w1} f_w (\tilde\nu/d)^2$ (a factor $\cos^2\theta$: 0.97, 0.88, 0.75,
0.50), and the arguments of $F_1$ and $F_2$ and the wall value
$\omega_w = 60\nu/(\beta_1 d_1^2)$ of SST. That distance is exact only on
meshes orthogonal to the wall, with the cell centres above the face centres.

- Each face of a wall boundary (types 103, 104, 301, 302) is the bilinear
  surface through its four nodes, the face of the cells
  ([Cell volume](numerics.md#cell-volume)):
  $\mathbf{x}(s,t) = \mathbf{x}_1 + s\mathbf{a} + t\mathbf{b} + st\mathbf{T}$,
  $\mathbf{a} = \mathbf{x}_2 - \mathbf{x}_1$, $\mathbf{b} = \mathbf{x}_4 - \mathbf{x}_1$,
  with the twist $\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$
  (corners in cyclic order). The distance to it is its minimum over the face,
  exact to the rounding of the size of the face and of the distance: a face plane
  to the resolution of its coordinates (every face of a 2-D mesh) takes the four
  triangles between its centre and its sides, any other face a certified minimum
  ([Code Structure](../development/structure.md#wall-distance) gives the
  algorithm). The four triangles of a twisted face would be wrong by up to
  $(|\mathbf{T}|/16)/y$, 6 % on a unit face twisted by 0.1 at the height 0.1.
- On a 2-D axisymmetric wedge (BC 200) the distance is measured in the
  meridional plane $(x, r)$, $r = \sqrt{y^2+z^2}$, from the cell centre to the
  wall curve. This is the distance to the wall surface of revolution, and it
  does not depend on the angle of the wedge. A wall on the faces of the wedge
  (faces 5–6) stops the run. It is the one exception to the single geometric
  model of the cells: for the wall distance the wall of the wedge is the
  surface of revolution of its meridional curve, for the volumes and the fluxes
  the cells are polyhedra with plane faces; the two differ by
  $r(1 - \cos(\delta/2)) \approx r\delta^2/8$.
- The distance and the face that gives it do not depend on the threads or on
  the ranks.
- Across periodic boundaries the wall goes on: the distance is the minimum over
  the periodic images of the wall too. A periodic connection (type 201) gives
  the translation from its source face to its own face, and on a 3-D mesh
  without sectors so does a block closed in $k$ by type `200` (its last node
  plane is the first one translated). The images are the integer combinations of
  these translations (along $x$ on an axisymmetric wedge). The two sides of a
  periodic boundary match only to the round-off of the coordinates in the mesh
  file: translations that
  differ by at most $10^{-2}$ of the thinnest cell across the periodic
  boundaries are the same, as in the type of the mesh. The thickness of a
  cell across a periodic face is the distance from the centre of the face to
  the centre of the opposite face of the cell, along the normal of the face;
  for a block closed in $k$ by type `200`, the move of its nodes per cell in
  $k$. The round-off of a coordinate grows with its size, so a tolerance
  relative to the size of the coordinates would depend on where the mesh
  lies; relative to the cells it does not. On a 3-D mesh whose
  sectors are closed by the rotational periodicity (type `200` on faces 5–6)
  the images are also the turns by multiples of the sector angle about the
  $x$ axis, up to half a turn either way. A channel periodic in $x$ and turned by
  30 degrees about $z$ gives translations that differ by 1.1e-4 of its cells with
  five significant digits, by 8.8e-4 with four, and by 0.15 with two, which stops
  the run. Translations that do not form a lattice stop the run: one that is not
  an integer combination of the others, more than three independent ones, or one
  not along $x$ on a wedge or with sectors. Only a bounded set of images can hold
  the nearest wall, and only those are searched: a mesh without periodic
  boundaries costs nothing more.
- With SA-rough each cell takes the sand-grain roughness $k_s$ of the face that gives the
  minimum distance; between faces at the same distance, the first one in the
  boundary file.

### Blowing / wall-transpiration correction

The `-blowcorr` suffix activates a correction (`RANS_Set_Blowing_Wall`) that modifies the wall treatment under wall injection or suction. SA uses its own blowing wall values; the $k$–$\omega$ models (SST, Wilcox 2006) switch between the corrected and uncorrected $\omega$ wall BC accordingly.

### Turbulent transport numbers and corrections

`[ARES-Turbulence]` also sets the turbulent Prandtl ($\mathrm{Pr}_t$, default 0.90) and Schmidt ($\mathrm{Sc}_t$) numbers used in the [effective conductivity](thermo.md#transport-properties), the optional turbulent-kinetic-energy coupling into the energy equation (`k-energy-coupling`), and a wall-roughness Prandtl correction (`Prt-correction`, requires an `-rough` SA model) that adjusts $\mathrm{Pr}_t$ near rough walls.

---

## Numerical treatment of the source terms

The turbulence source terms are added to the residual and advanced with the same explicit
Runge–Kutta stages as the mean flow. Their destruction parts are stiff near walls: for the
$k$–$\omega$ models $\omega$ grows like $1/y^2$ towards the wall, and an explicit stage on
$\dot q = -D(q)$ gives $q - \Delta t\,D = q\,(1 - \Delta t\,D/q)$, negative for $\Delta t\,D/q > 1$. Here
$q$ is the conserved turbulence variable; the destruction rate $D/q$ is at most $d$, the derivative
of the destruction ($d/2$ for $\tilde\nu$ and $\omega$, whose destruction is quadratic):

| Model | Equation | $d = \partial D/\partial q$ |
|-------|----------|------------------------------|
| SST / Wilcox | $k$ | $\beta^\ast\,\omega$ |
| SST / Wilcox | $\omega$ | $2\,\beta\,\omega$ |
| SA | $\tilde\nu$ | $2\,c_{w1}\,f_w\,\tilde\nu/y^2$ |

($y$ the wall distance; $f_w$ is held fixed.) Steady and time-accurate runs treat this
stiffness in two different ways.

### Steady runs: point-implicit factor

With `point-implicit = .true.` (the default is `.false.`) the increment of each turbulence equation is divided
by $1 + \Delta t\,d$, with the local time step $\Delta t$ of the cell:

$$
q^{n+1} = q^n + \frac{\Delta t}{1 + \Delta t\,d}\,\bigl( F(q^n) + S(q^n) \bigr),
$$

where $F$ stands for the fluxes and $S = P - D$ for the source (in each Runge–Kutta stage, with the
residual of the stage). For $\omega$ of SST, $d$ also takes the derivative of the cross-diffusion where that term is
positive, $\max(D_\omega^{cd}, 0)/(\rho\omega)$: there it decreases as $\omega$ grows, and it is the only part of the
source that the factor can take besides the destruction. This is the point-implicit form of Bussing and Murman,
$D\,\Delta q = \Delta t\,(F + S)$ with $D = 1 - \Delta t\,\partial S/\partial q$ approximated by its
destruction part $1 + \Delta t\,d$, in its first-order form, «suitable for steady-state
calculations» (H. C. Yee, J. L. Shinn, NASA TM-89415, 1986, eq. 3.1). Like any approximation of the
implicit operator of a delta form, it leaves the steady state where the residual vanishes
(T. H. Pulliam, *Time Accuracy and the Use of Implicit Methods*, AIAA Paper 93-3360, eq. 8): the
converged solution is that of $F + S = 0$, whatever the time step, and the factor only lets
steady runs keep the mean-flow CFL near walls. The residuals of $\tilde\nu$, $k$ and $\omega$ printed during
the run (the change of the variable over the step) are scaled by the same factor, so near walls they
are smaller than without it.

With [`k-energy-coupling = true`](#energy-coupling-k-energy-coupling) the energy of the mean flow has
the sink $-S_k$ of the $k$ equation, whose derivative with respect to $\rho k$ is $+\beta^\ast\omega$, so
the factor is that of the two equations together: the implicit operator of the $k$ and energy rows
is $D = \begin{pmatrix} 1 + \Delta t\,\beta^\ast\omega & 0 \\ -\Delta t\,\beta^\ast\omega & 1 \end{pmatrix}$, and the
energy takes its own increment plus $\Delta t\,\beta^\ast\omega$ times that of $k$,

$$
\Delta(\rho E_m) = \Delta t\,\bigl(F_E - S_k\bigr) + \Delta t\,\beta^\ast\omega\,\Delta(\rho k), \qquad
\Delta(\rho k) = \frac{\Delta t}{1 + \Delta t\,\beta^\ast\omega}\,\bigl(F_k + S_k\bigr).
$$

The two increments add up to $\Delta t\,(F_E + F_k)$, the fluxes of $E_m + k$ alone, at every stage, and
the converged solution, $F_k + S_k = 0$ and $F_E = S_k$, does not depend on the time step.

```ini
[ARES-Turbulence]
point-implicit = .true.   ; point-implicit factor on the turbulence increment, steady runs (default .false.)
```

### Time-accurate runs: explicit sources

In time-accurate runs (`time-accurate = .true.`, or a time-accurate HYDRA coupling) the factor is not applied,
whatever `point-implicit` says; when `input.ini` sets `point-implicit = .true.`, a warning at start-up recalls it (a
run that leaves the key out takes the default and gets no warning). Since
$(F + S)/(1 + \Delta t\,d) = (F + S)(1 - \Delta t\,d) + O(\Delta t^2)$, the factor would make the turbulence equations,
and the mean flow through $\mu_t$, first-order accurate in time. With the $k$–$\omega$ models the error of the factor
would not even vanish when the mesh and the step are refined together, as in an explicit run, whose global step goes as
$h^2$ in the wall cells: there $\omega \approx 6\nu/(\beta_1 y^2)$ makes $\Delta t\,d$ independent of the mesh, and of
order ten at the usual `vnn`. With Spalart–Allmaras $\Delta t\,d$ grows with $f_w\,y_1^+$ ($y_1^+$ the height of the
first cell in wall units): it vanishes as the mesh is refined, but is of order one at the usual heights.
The sources are explicit, like the fluxes, and the local time step of each cell
([Time Integration](time-integration.md#stability-conditions)) includes a bound of their stiffness:

$$
\frac{1}{\Delta t_i} = \frac{1}{\Delta t_i^\text{cfl, vnn}} + \hat d_i ,
\qquad \hat d_i \ge d_i = \max_q \frac{\partial D_q}{\partial q} ,
$$

with $\hat d$ a bound of the rates $d$ of the table that the state of the cell at the start of the step gives
alone, without gradients or ghost cells:

| Model | $\hat d$ | $\hat d / d$ |
|-------|----------|--------------|
| SST | $\max(\beta^\ast, 2\max(\beta_1, \beta_2))\,\omega$ | $\beta = F_1\beta_1 + (1 - F_1)\beta_2$, and $F_1$ takes gradients: $\beta_2/\beta_1 = 1.10$ where $F_1 = 1$ (next to walls), 1 where $F_1 = 0$ |
| Wilcox 2006 | $\max(\beta^\ast, 2\beta_0)\,\omega$ | $\beta = \beta_0 f_\beta$, and $f_\beta \le 1$ takes the velocity gradient: $1/f_\beta$, 1 in 2-D ($\chi_\omega = 0$), at most $1/0.85$ |
| SA | $2\,c_{w1}\,(1 + c_{w3}^6)^{1/6}\,\tilde\nu/y^2$ | $f_w = g\,[(1 + c_{w3}^6)/(g^6 + c_{w3}^6)]^{1/6}$ grows with $g \ge 0$ towards $(1 + c_{w3}^6)^{1/6} = 2.005$: $2.005/f_w$, 1 where $r$ reaches its cap of 10, about 2 where $r \approx 1$ (an equilibrium boundary layer, where the destruction is slow against the diffusion) |

($y$ the wall distance of the source, moved below a rough wall with SA-rough.) The Reynolds-stress model SSG/LRR
takes no bound: none was derived for its sources, and its $\hat d$ is zero. Then
$\Delta t\,d \le \Delta t\,\hat d \le 1$: a forward-Euler stage cannot make
$\tilde\nu$, $k$ or $\omega$ negative through their destruction, and neither can RK2 (Heun) and RK3, whose
stages are convex combinations of forward-Euler stages (W. Hundsdorfer, *Numerical Solution of
Advection–Diffusion–Reaction Equations*, lecture notes, 2000, §6.2, Theorem 6.4; C.-W. Shu,
S. Osher, *J. Comput. Phys.* 77, 1988, Proposition 2.1 of its report version, ICASE Report 87-33).
The destruction of $\tilde\nu$ and of $\omega$ is quadratic, $D = (d/2)\,q$, which leaves them a
margin of 2. Where the destruction is faster than the transport, as in the wall cells of fine
meshes, above all with SST, the step is smaller than the convective and diffusive limits.

The bound covers the destruction only. Three sinks stay outside it: the cross-diffusion of SST,
$2(1-F_1)\,\rho\,\sigma_{\omega 2}\,\nabla k\cdot\nabla\omega/\omega$, negative where $\nabla
k\cdot\nabla\omega < 0$ (Wilcox 2006 switches its cross-diffusion off there); with `sst-production =
compressible`, the dilatation part of the production of SST, $-\tfrac23\rho k\,\nabla\cdot\mathbf{u}$
and $-\tfrac23\gamma\rho\omega\,\nabla\cdot\mathbf{u}$, negative where the flow expands; and the
production of SA-R and SA-RC, which can be negative («the SA-RC production term can be negative»,
NASA Turbulence Modeling Resource). Like the transport, they can make a stage negative, and the count
at the end of the run shows it.

The convex combination keeps the bound of the forward-Euler stage only where that bound holds on the
state each stage starts from: Shu and Osher assume that the forward-Euler step keeps the property,
under one step restriction, for every state it is applied to (ICASE Report 87-33, (2.7) and
Proposition 2.1), while here the restriction depends on the state through $d$, and the rates are
those of the state at the start of the step. A stage from a state where the rate is $d_s$ takes
$q_s - \Delta t\,(d_s/2)\,q_s$ through a quadratic destruction, positive while $\Delta t\,d_s < 2$; with
$\Delta t\,\hat d \le 1$ this holds as long as $d_s < 2\,\hat d$. The margin of 2 of $\tilde\nu$ and $\omega$
therefore covers a rate that grows within the step up to twice the bound at the start, and no more; $k$, whose
destruction is linear, keeps only the margin by which its rate $\beta^\ast\omega$ lies below $\hat d$. An
impulsive start of a $k$–$\omega$ model from a state far from the wall value of $\omega$ goes beyond
that: during the first step the wall condition raises $\omega$ in the wall cells by orders of
magnitude, and with the full step the explicit destruction overshoots ($\omega < 0$, then the floor
and a very large $\mu_t$). Start such runs with the CFL ramp (`cfl-rise-threshold`).

At the end of a RANS run (one- and two-equation models) in which a stage made $\tilde\nu$, $k$ or $\omega$ negative
before the realizability floor, ARES prints how many cell updates had such a stage; a run in which none did prints no
such line:

```text
 Turbulence variables below zero after a stage (set to the floor): <n>
```

In a steady run the count covers the whole run: the start, and the fronts of $k$ that the fluxes
push below zero where $k$ is about zero, which can last to the last iteration; it is a total, which
convergence does not bring back to zero.

### Relation to under-relaxation

Point-implicit damping is a *physically weighted, automatic* form of
under-relaxation: the effective relaxation factor $1/(1+\Delta t\,d)$ is
large (≈1) where the source is mild and small only where the destruction is
stiff, and it scales with the local time step. A constant global
under-relaxation factor on the turbulence update achieves a similar
stabilising effect but is blunter: it slows convergence everywhere and must be
hand-tuned per case and grid.

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
13. NASA Langley Turbulence Modeling Resource, "The Spalart–Allmaras Turbulence Model," <https://tmbwg.github.io/turbmodels/spalart.html> (formerly turbmodels.larc.nasa.gov/spalart.html): minimum distance to the wall.
