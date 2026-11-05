# Spatial Discretization

ARES uses a cell-centred finite volume method (FVM) on structured multi-block grids. This page describes the discretization framework, the MUSCL reconstruction that provides second-order accuracy, and the slope limiters that ensure monotonicity.

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
P_L = P_i + \bar{s}_L\,\delta l_L,\qquad
P_R = P_{i+1} - \bar{s}_R\,\delta l_R,
$$

with $\delta l_L,\delta l_R$ the cell-centre-to-interface distances.

| `space-reconstruction` | Meaning |
|------------------------|---------|
| `first-order` | Donor-cell; most robust, most diffusive |
| `MUSCL` | Second-order MUSCL with the selected limiter |
| `none` | No reconstruction (used for some pure-diffusion configurations) |

!!! note "Choosing the reconstruction"
    `space-reconstruction = MUSCL` limits the slopes with the function named by `flux-limiter` (see
    [Flux Limiters](#flux-limiters)); if `flux-limiter` is omitted, van Leer is used and a warning is printed.
    `space-reconstruction = first-order` sets every slope to zero (donor-cell); a `flux-limiter` given with it is
    ignored, with a warning.

---

## Flux Limiters

The limiter $\phi(a,b)$ is selected by `flux-limiter`. All are expressed via the slope ratio $r = b/a$.

| `flux-limiter` | Character |
|----------------|-----------|
| `minmod` | Most diffusive TVD limiter; strictly monotone |
| `vanleer` | Smooth, differentiable; good dissipation/accuracy balance |
| `vanalbada` | Smooth; less dissipative than van Leer near extrema |
| `MC` (or `mc`, as in MOSE) | Monotonized-central; steepest gradients inside the TVD region |
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

### TVD Diagram

In Sweby's $(r,\phi)$ diagram every second-order TVD limiter must lie inside the shaded region bounded below by the **minmod** envelope and above by the **superbee** envelope, while passing through the consistency point $\phi(1)=1$:

$$
0 \le \phi(r) \le \min(2r,\,2),\qquad \phi(r)=0\ \text{for } r\le 0 .
$$

A limiter sitting near the lower (minmod) edge is the most dissipative and monotone; one near the upper (superbee) edge is the sharpest but can steepen smooth extrema. Van Leer and van Albada are smooth ($C^\infty$) curves through the middle of the region; MC follows the monotonized-central envelope. LIMO3 is a third-order limiter that leaves the strict second-order TVD region in smooth flow to recover third-order accuracy, while staying non-oscillatory.

| Limiter | Dissipation | Smoothness | Position in TVD region |
|---------|:-----------:|:----------:|------------------------|
| minmod | High | $C^0$ | lower envelope |
| van Leer | Medium | $C^\infty$ | mid (smooth) |
| van Albada | Medium-low | $C^\infty$ | mid (smooth) |
| MC | Low | $C^0$ | monotonized-central edge |
| superbee | Lowest | $C^0$ | upper envelope |
| LIMO3 | Low (3rd-order smooth) | piecewise | third-order (leaves 2nd-order TVD region) |

---

## Gradient Computation for Diffusive Fluxes

Velocity and temperature gradients for the viscous flux are computed at each interface from a **10-point stencil**:

- **2 points** in the face-normal direction — the two cells sharing the face, giving the normal derivative directly;
- **8 points** in the two tangential directions — four cells per direction, each tangential derivative being the average of the two consecutive differences.

The resulting computational-space gradient $(\xi,\eta,\zeta)$ is mapped to Cartesian $(x,y,z)$ with the **face-metric tensor** $M_{3\times3}$, taken as the average $M = \tfrac12(M_1 + M_2)$ of the two adjacent cells' metrics, computed from the grid geometry and stored per face. The same metric provides the physical spacing $\Delta l$ used in the MUSCL reconstruction.

Block-connection (`101`, `201`) and chimera (`102`) faces use the inverse mean instead, $M_f = \left[\tfrac12\left(M_1^{-1} + M_2^{-1}\right)\right]^{-1}$, the inverse of the mean of the two cell edge vectors (`Face_Metric`): the one-index difference across the face is then divided by the centre-to-centre distance. The mean of $M$ divides it by the mean of $1/h$ instead: across a large cell-size jump at a block interface (86:1 makes it 22 times too steep) the normal gradient is far too steep and the explicit diffusion on the large cell becomes unstable, although its own time step is fine.

The viscous flux of a chimera face is the mean of two fluxes: one with the
central differences along the face of the row of the cell, one with those of
the row of the ghost cells, filled from the donors. With a stress linear in the
gradient (laminar, Spalart–Allmaras) the mean is the flux of the mean gradient,
whose differences along the face are the mean of those of the two rows, as on
an interior face (Vinokur 1986, NASA CR-177416, eqs. 90–91, p. 23). Beyond the
edges of the block face the row of the ghost cells is extrapolated linearly, in
the variables whose differences the stencil takes (velocity, temperature and the
turbulence variables per unit mass, at the pressure of the ghost cell): the
difference there is the one-sided one, exact for a linear field (a copy of the
ghost cell would halve it, and the face would take 3/4 of the derivative of a
linear field along the edge). In the k direction of a 2-D mesh the neighbours of
a ghost cell are the ghost cell turned by the wedge angle, as those of a cell
are; across a block one cell wide in another direction they are the ghost cell
itself. The ghost cells of a chimera face take the wall distance and the
roughness of their receiver cell (zeroth order on purpose: a linear extrapolation
can reach $y_n \le 0$ next to walls).

---

## Faces Between Blocks

A face shared by two blocks (block connection `101`, periodic connection `201`)
has one flux, computed once and given to both cells with opposite signs, as an
interior face. Of the two boundary records of the face, the one of the cell on
the low-index side (whose face is 2, 4 or 6 while the other is 1, 3 or 5)
computes the inviscid and the viscous flux. When the two faces have the same
parity, as where blocks meet with their index directions permuted, the record
with the smaller (block, face, $k$, $j$, $i$) computes it, so that the result
does not depend on the order of the records in `bc.txt`. The face therefore
conserves mass, momentum and energy exactly, whatever the Riemann solver and the
preconditioning, on any number of MPI ranks; between ranks the other cell's
share travels in one more message per Runge–Kutta stage. The two records of a
face must point at each other, and the block across a face must have at least
two cells across it at every multigrid level: otherwise the run stops at set-up.

The inviscid flux of a face between blocks is, bit for bit, that of the same
face in the unsplit mesh when the index directions of the two blocks agree
across the face and the two blocks store the nodes of the common face with the
same bits (the ghost cell takes its geometry from the nodes of the other block,
so nodes that differ in the last digits give a slightly different face). This
is an identity of the fluxes, not of the solution: in 2-D and 3-D a cell next to
the cut adds the flux of the face between blocks after those of its interior
faces, so its sum is taken in another order than in the unsplit mesh and the two
solutions differ at round-off.

---

## References

1. B. van Leer, "Towards the ultimate conservative difference scheme. V," *J. Comput. Phys.* 32 (1979).
2. P. K. Sweby, "High resolution schemes using flux limiters for hyperbolic conservation laws," *SIAM J. Numer. Anal.* 21 (1984).
3. M. Čada, M. Torrilhon, "Compact third-order limiter functions for finite volume methods," *J. Comput. Phys.* 228 (2009).
4. A. Jameson, W. Schmidt, E. Turkel, "Numerical solution of the Euler equations by finite volume methods using Runge–Kutta time-stepping schemes," AIAA-81-1259, 1981.
