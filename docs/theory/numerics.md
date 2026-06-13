# Spatial Discretization

ARES uses a cell-centred finite volume method (FVM) on structured multi-block grids. This page describes the discretization framework, the MUSCL reconstruction that provides second-order accuracy, the slope limiters that ensure monotonicity, and the shock-detection sensor that adaptively reduces the order near discontinuities.

---

## Finite Volume Framework

The integral conservation law over a cell $\Omega_i$ with boundary $\partial\Omega_i$ is

$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{\Omega_i}\!\mathbf{U}\,\mathrm{d}V
+ \oint_{\partial\Omega_i}\!\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
= \oint_{\partial\Omega_i}\!\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
+ \int_{\Omega_i}\!\mathbf{S}\,\mathrm{d}V .
$$

With midpoint quadrature this gives the semi-discrete update

$$
\frac{\mathrm{d}\mathbf{U}_i}{\mathrm{d}t}
= -\frac{1}{V_i}\sum_{f}\bigl(\mathbf{F}^{c}_f - \mathbf{F}^{v}_f\bigr)\,A_f
+ \mathbf{S}_i,
$$

summed over the six faces $f$ of each hexahedral cell. At each face:

- **Convective flux** — left and right primitive states are reconstructed, then a [Riemann solver](riemann-solvers.md) returns the numerical convective flux from those states and the face normal.
- **Diffusive flux** — velocity and temperature gradients are computed at the face with a wide stencil and mapped to Cartesian coordinates by the face-metric tensor.

---

## MUSCL Reconstruction

First- and second-order spatial accuracy are available (`space-reconstruction`). Second order is achieved by piecewise-linear MUSCL reconstruction (Monotone Upstream-centred Schemes for Conservation Laws). At each interface ARES uses a four-point stencil $\{i{-}1, i, i{+}1, i{+}2\}$ of the primitive variables $(p,u,v,w,h,\dots)$:

1. Compute slopes from consecutive cell values and physical spacings $\Delta l_0,\Delta l_1,\Delta l_2$:

$$
s_0 = \frac{P_i - P_{i-1}}{\Delta l_0},\quad
s_1 = \frac{P_{i+1} - P_i}{\Delta l_1},\quad
s_2 = \frac{P_{i+2} - P_{i+1}}{\Delta l_2}.
$$

2. Limit the slopes to suppress oscillations near discontinuities:

$$
\bar{s}_L = \phi(s_1, s_0),\qquad
\bar{s}_R = \phi(s_2, s_1),
$$

where $\phi$ is the chosen limiter.

3. Reconstruct interface values:

$$
P_L = P_i + \beta\,\bar{s}_L\,\delta l_L,\qquad
P_R = P_{i+1} - \beta\,\bar{s}_R\,\delta l_R,
$$

with $\delta l_L,\delta l_R$ the cell-centre-to-interface distances and $\beta$ a blending parameter set by the shock detector.

| `space-reconstruction` | Meaning |
|------------------------|---------|
| `first-order` | Donor-cell ($\beta=0$); most robust, most diffusive |
| `MUSCL` | Second-order MUSCL with the selected limiter |
| `MUSCL-SD` | MUSCL with the **shock-detection** sensor modulating $\beta$ (see below) |
| `none` | No reconstruction (used for some pure-diffusion configurations) |

!!! note "Role of β"
    $\beta = 1$ recovers full second-order MUSCL; $\beta = 0$ reduces to first order (donor-cell) in the immediate vicinity of shocks. The `MUSCL-SD` variant lets the shock sensor smoothly modulate $\beta$ between these limits.

---

## Flux Limiters

The limiter $\phi(a,b)$ is selected by `flux-limiter`. All are expressed via the slope ratio $r = b/a$.

| `flux-limiter` | Character |
|----------------|-----------|
| `minmod` | Most diffusive TVD limiter; strictly monotone |
| `vanleer` | Smooth, differentiable; good dissipation/accuracy balance |
| `vanalbada` | Smooth; less dissipative than van Leer near extrema |
| `MC` | Monotonized-central; steepest gradients inside the TVD region |
| `superbee` | Sharpest TVD limiter; excellent at steep gradients |
| `LIMO3` | Third-order limiter (Čada–Torrilhon) — third-order accurate in smooth regions while remaining non-oscillatory |
| `none` | No limiting (unlimited reconstruction) |

**MinMod**

$$
\phi(a,b)=\begin{cases}\operatorname{sign}(a)\min(|a|,|b|) & ab>0\\ 0 & \text{otherwise}\end{cases}
$$

**Van Leer**

$$
\phi(a,b)=\frac{r+|r|}{1+|r|}\,a
$$

**Van Albada**

$$
\phi(a,b)=\frac{r^2+r}{1+r^2}\,a
$$

**MC (monotonized central)**

$$
\phi(a,b)=\max\!\bigl(0,\min(2r,\tfrac12(1+r),2)\bigr)\,a
$$

**Superbee**

$$
\phi(a,b)=\max\!\bigl(0,\min(2r,1),\min(r,2)\bigr)\,a
$$

| Limiter | Dissipation | Smoothness |
|---------|:-----------:|:----------:|
| minmod | High | $C^0$ |
| van Leer | Medium | $C^\infty$ |
| van Albada | Medium-low | $C^\infty$ |
| MC | Low | $C^0$ |
| superbee | Lowest | $C^0$ |
| LIMO3 | Low (3rd-order smooth) | piecewise |

---

## Shock Detection

Even with TVD limiters, MUSCL reconstruction can produce oscillations near strong shocks (the carbuncle phenomenon). With `MUSCL-SD`, ARES uses a pressure-based sensor to detect shocks and smoothly degrade to first order in their vicinity.

A $3\times3\times3$ stencil of pressures centred on cell $(i,j,k)$ gives second-order undivided differences in each direction:

$$
\kappa_d = \frac{|p_{+}-2p_0+p_{-}|}{p_{+}+2p_0+p_{-}},\qquad d=1,2,3,
$$

and the shock strength is $\sigma=\max(\kappa_1,\kappa_2,\kappa_3)$. A smooth blending with threshold $\Delta$ gives

$$
\beta=\begin{cases}1-\tanh\!\bigl(10(\sigma\Delta)^3\bigr) & \sigma<1/\Delta\\ 0 & \text{otherwise}\end{cases}
$$

so $\beta\to1$ in smooth flow (full second order) and $\beta\to0$ at shocks (first order). This $\beta$ multiplies the limited slope in the reconstruction step.

!!! note "Solver-level use"
    Some Riemann solvers (notably the rotated HLLC) also consult the shock flag to tune their numerical dissipation, adding a second layer of robustness near discontinuities.

---

## Gradient Computation for Diffusive Fluxes

Velocity and temperature gradients for the viscous flux are computed at each interface from a stencil spanning the face-normal and tangential directions. The computational-space gradient $(\xi,\eta,\zeta)$ is mapped to Cartesian $(x,y,z)$ with the **face-metric tensor** $M_{3\times3}$, computed from the grid geometry and stored per face. The same metric provides the physical spacing $\Delta l$ used in the MUSCL reconstruction.

---

## References

1. B. van Leer, "Towards the ultimate conservative difference scheme. V," *J. Comput. Phys.* 32 (1979).
2. P. K. Sweby, "High resolution schemes using flux limiters for hyperbolic conservation laws," *SIAM J. Numer. Anal.* 21 (1984).
3. M. Čada, M. Torrilhon, "Compact third-order limiter functions for finite volume methods," *J. Comput. Phys.* 228 (2009).
4. A. Jameson, W. Schmidt, E. Turkel, "Numerical solution of the Euler equations by finite volume methods using Runge–Kutta time-stepping schemes," AIAA-81-1259, 1981.
