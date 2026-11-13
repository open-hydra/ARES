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

### Cell volume

A cell is the image of the unit cube by the trilinear map
$\mathbf{x}(\xi,\eta,\zeta)$ of its eight nodes. Each face is the bilinear
(doubly ruled) surface through its four nodes, the same for the two cells that
share it, and its area vector, the one the fluxes take, is half the cross
product of its diagonals [5, pp. 15–16]. The volume of the cell is the integral
of the Jacobian of the map,

$$
V = \int_0^1\!\!\int_0^1\!\!\int_0^1 \det J\,\mathrm{d}\xi\,\mathrm{d}\eta\,\mathrm{d}\zeta
  = \frac{1}{6}\sum_{d=1}^{3}\bar{\mathbf{e}}_d\cdot\left(\mathbf{S}_d^{-}+\mathbf{S}_d^{+}\right),
$$

with $\bar{\mathbf{e}}_d$ the mean of the four edges of the cell along index $d$
and $\mathbf{S}_d^{\mp}$ the area vectors of the two faces across index $d$,
both along $+d$: an identity of polynomials in the 24 coordinates, which can be
checked term by term ($\det J$ has degree at most two in each of $\xi$, $\eta$,
$\zeta$). It is also
$\frac13\sum_f\bar{\mathbf{x}}_f\cdot\mathbf{S}_f$ over the six faces, outward,
with $\bar{\mathbf{x}}_f$ the mean of the nodes of the face: the divergence
theorem for $\nabla\cdot\mathbf{x} = 3$, with
$\int_f\mathbf{x}\cdot\mathrm{d}\mathbf{S} = \bar{\mathbf{x}}_f\cdot\mathbf{S}_f$ on a
bilinear face. So the discrete divergence of a field that grows as $\mathbf{x}$
is exact on every cell and, as two neighbours share their face, the volumes of
the cells sum to the volume enclosed by the boundary faces, with no overlap and
no gap [5, p. 16]. On a cell with plane faces (a 2-D mesh, the wedge, a 3-D
mesh of plane faces) it is the volume of the polyhedron. A sum of tetrahedra,
each face cut along a diagonal, gives that volume only where the faces are
plane: the two cells that share a twisted face cut it along different
diagonals, and each volume misses or counts twice
$\mathbf{T}\cdot\mathbf{S}_f/12$, with
$\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$ the
twist of the face.

The volume has the sign of the block (negative on a left-handed cell) and is
stored in absolute value. A cell whose volume has not the sign of the
determinant at its centre is turned over somewhere and stops the run at setup
(below). The volume is computed from differences of the nodes, so its rounding
is that of the size of the cell, not that of its coordinates; the compiler must
keep the order of these differences, which the code enforces (VOLATILE
variables, also under a fast floating-point model).

### Checks of the cells at setup

A cell whose determinant is zero to $10^{-12}$ of the product of its edge lengths
stops the run at setup, with its block and cell. So does a cell folded at
one of its corners, or turned against its block, with its block, cell and corner
node: a cell is the image of the unit cube by the trilinear map of its eight
nodes, and at a corner the Jacobian of that map is
$J = \mathbf{e}_i\cdot(\mathbf{e}_j\times\mathbf{e}_k)$, with
$\mathbf{e}_i, \mathbf{e}_j, \mathbf{e}_k$ the three edges from that corner along
$i$, $j$ and $k$. At every corner it must have the sign of the block (the sign of
the determinant of the first cell, from which the face normals of the block take
their orientation) and exceed $10^{-12}\,|\mathbf{e}_i||\mathbf{e}_j||\mathbf{e}_k|$;
a corner with an edge of zero length (nodes that coincide on an axis) is not
checked. On a plane 2-D mesh $J = \Delta z\, j(\xi, \eta)$, with $j$ the Jacobian
of the bilinear map of the section, $\mathbf{x} = \mathbf{a} + \mathbf{b}\xi +
\mathbf{c}\eta + \mathbf{d}\xi\eta$: $j = (\mathbf{b} + \mathbf{d}\eta)\times
(\mathbf{c} + \mathbf{d}\xi)$ is affine in $(\xi, \eta)$ because
$\mathbf{d}\times\mathbf{d} = 0$, so its sign in the cell is decided at the
corners. On the 2-D wedge $J = R\,\sin\delta\; j_m(\xi, \eta)$, with $R > 0$ the
bilinear radius when no node is on the axis and $j_m$ the Jacobian of the
meridional section, affine as above: the corners decide again. In 3-D the corners
are a necessary condition only; a cell whose volume (above) has not the
sign of the determinant at its centre is turned over somewhere, and stops the
run at setup as well ($\det J$ has the sign of the block in every point of a
cell that is not turned over, and so has its integral). These checks run on every
multigrid level too, whose nodes are every other node of the level above: the stop of
the corners names the level, the other stops give the block and the cell of that
level. Behind a face of zero
area (an axis) the ghost cells are the reflection of the first cell through the
axis, whose determinant has the opposite sign by construction: their sign is
not checked (a wall on such a face stops, see [Face metric `stencil`](#face-metric-stencil)).
The first layer of ghost cells behind a chimera face is
checked too (see [Metrics of the ghost cells](#metrics-of-the-ghost-cells)).

A degenerate or inverted cell is found inside the OpenMP threads of the set-up:
the first thread that finds one writes it, and every MPI rank stops together,
with exit code 1, after the parallel region.

### Rotational periodicity (type 200)

Type 200 on the $k$ faces closes a mesh that turns about the $x$ axis: an
axisymmetric wedge (one cell in $k$) or a 3-D sector (several). The flow
beyond face 5 is the flow next to face 6 turned back by the angle $\delta$
of the wedge or sector, and the flow beyond face 6 the flow next to face 5
turned forward. With $g = +1$ on face 5 and $g = -1$ on face 6, the ghost
cell takes the state of its source cell (on a wedge the cell itself, on a
sector the cell next to the opposite $k$ face) with the velocity turned
about $x$ by $\varphi = -g\,\delta$:

$$
v_g = v\cos\varphi - w\sin\varphi = v\cos\delta + g\,w\sin\delta, \qquad
w_g = v\sin\varphi + w\cos\varphi = -g\,v\sin\delta + w\cos\delta .
$$

$u$ and the scalars do not change under a turn about $x$, and are copied.
The ghost cells and the fluxes of these faces take $\cos\delta$ and
$\sin\delta$ computed once, with $\delta$, when the mesh is read. The angle
$\delta$ and the type of the mesh (plane, wedge, sector or translated block)
come from the nodes of every block ([Running ARES](../user/using.md), note *Mesh type*).

#### Axisymmetric wedge

A 2-D axisymmetric case is a wedge of the 3-D mesh: one cell in $k$, whose two
$k$ faces lie in the meridional planes $\theta = \mp\delta/2$ and carry the
boundary type `200`. ARES integrates the 3-D equations on the wedge; there are
no axisymmetric source terms. The ghost cell across each wedge face is the cell
itself rotated by $\mp\delta$ about the $x$ axis, so every face flux follows
from the state of the cell:

- inviscid: the pressure and the momentum carried by the cell state rotated to
  the face (no Riemann problem). In the radial momentum the two faces give the
  $+p/r$ and the centrifugal terms of the axisymmetric equations;
- viscous: the flux of a face between the cell and its rotated image, with the
  stencil of an interior face (the $i$ and $j$ neighbours of the image rotated
  in the same way) and, with `diffusive-metric = stencil`, the face metric of
  the same differences of the positions, the row of the image the row of the
  cell turned (exact for the linear fields that the rotation allows, see the 3-D
  sector below); with `mean`, the inverse-mean face metric of the cell and its
  image, as on the faces between blocks (`Face_Metric`). Both are computed once
  at setup (they depend on the mesh only).

On a planar 2-D mesh (faces 5 and 6 parallel, no wedge angle) type `200` makes
the image the cell itself, and the fluxes of the two faces cancel.

For an axisymmetric field the mass, energy and turbulence fluxes of the two
faces are equal and opposite, and only the momentum fluxes are applied. Without
swirl the viscous traction of each face is the normal stress

$$
\tau_{\theta\theta} = 2\mu\,\frac{u_r}{r} + \lambda\,\nabla\!\cdot\mathbf{u},
\qquad \lambda = -\tfrac23\,\mu ,
$$

with $\mu$ including the eddy viscosity, and the two faces together give the
$-\tau_{\theta\theta}/r$ of the radial momentum equation [6, 7]; with swirl
they also carry $\tau_{r\theta}$ and $\tau_{x\theta}$. Near the axis
$\tau_{\theta\theta}/r$ is as large as the radial stress terms it balances: in
Homann's stagnation-point flow [8] $\tau_{rr} = \tau_{\theta\theta}$, and
$\tau_{rr}/r$ and $\tau_{\theta\theta}/r$ cancel exactly while each grows like
$1/r$ towards the axis.

The azimuthal derivatives of the $i$ and $j$ faces are differences across the
rotated ghost cells divided by the $k$ edge of the cell, which is a chord: they
carry a relative error of order $\delta^2$ that the mesh does not reduce
($\delta^2/12$ of the hoop stress $2\mu u_r/r^2$ in the radial balance,
$2.5\times10^{-5}$ for $\delta = 1^\circ$). Keep the wedge thin.

Put the axis nodes at a small radius ($10^{-6}$ of the radial size of the first
cell), turned with the other nodes of their $k$ plane: a viscous run on a wedge
whose boundary faces on the axis have zero area stops at setup, whether its
axis nodes are exactly at $r = 0$ or those of the two $k$ planes coincide, with
every `diffusive-metric`. On such an axis the faces have zero area and no
normal, so their ghost cells copy the cells. With `mean` the wedge faces and the
faces $i$ of the first row of cells take the derivative in $r$ across the axis
with the ghost cell one cell away and would see half of
$\partial u_r/\partial r$: an error of the hoop stress of those cells that
refinement does not reduce ($1/6$ of $2\mu u_r/r^2$ in Homann's flow, on 32 and
64 cells). With `stencil` the copy sits at the cell and the faces take the
difference on one side of it: in Homann's flow, linear in $r$, the first row
gets the hoop stress of a small axis radius, so the error of $1/6$ is gone; a
field that is not linear in $r$ stays farther from the result of the small
radius. The stop stays with every metric, a declared limit: the axis at $r = 0$
is not treated completely (ghost cells beyond the axis, not copies).

The value of a cell is the mean of the solution over the cell, and the mean of
a linear field over a cell is its value at the centroid of the cell volume
[13, Eq. (3.1) and (3.12)]. ARES places the value at the cell centre, the mean
of the eight nodes. The two points coincide on a parallelepiped, but not next
to the axis of a wedge: in the first row of cells, whose inner nodes are at a
small radius, the centroid lies at two thirds of the height of the cell from the
axis and the cell centre at one half, at every refinement. There the operators
that take the cell centres are not exact for a linear field given by its cell
means: with uniform radial spacing and `stencil`, the radial gradient at the
faces between the first and the second row is $1/9$ too small ($1/45$ between
the second and the third), the cell gradient of the turbulence sources of the
first row is off by $1/9$, and the term $u_r/r$ of the hoop stress is $1/3$ too
large in the first row ($1/27$ in the second). These errors do not decrease as
the mesh is refined. In a one-dimensional model of the viscous terms of the
radial momentum next to the axis, with these positions, the error of the
solution decreases at first order in the first row and at second order in its
mean weighted by the volume.

Next to the axis the hoop stress also stiffens the explicit step. In the radial
momentum of the first row of cells, at $r \approx \Delta r/2$, it is a damping
$(2\mu+\lambda)\,u_r/r^2 \approx \tfrac{16}{3}\,\mu\,u_r/\Delta r^2$ ($\mu$ with
the eddy viscosity), with couplings to the neighbours through $\lambda\,\nabla\cdot\mathbf{u}$.
The local time step includes them ([Time Integration](time-integration.md#local-time-step)),
so the `vnn` limits of the rule hold next to the axis too.

#### 3-D sector

On a 3-D sector (type `200` with more than one cell in $k$) faces 5 and 6
make one interior face of the whole annulus, computed once, by the record of
face 6, for the inviscid and for the viscous flux. Beyond face 6 lies the cell
next to face 5 turned forward by the sector angle, a turn $R$ about $x$: the
$i$ and $j$ differences of the image row are those of that row turned, the
metric of the image is that of its cell turned, $M R^{\mathsf T}$ (each row of
$M = \partial\xi/\partial\mathbf{x}$ is a vector), and the wall distance and
roughness of the face are the mean of the two cells, as on an interior face.
The cell next to face 6 takes the flux and the cell next to face 5 the same
flux with its momentum turned back, $R^{\mathsf T}F$: mass, energy, turbulence
variables, axial momentum and angular momentum about $x$ are conserved. A
translated sector (faces 5 and 6 parallel) is the same with $R = I$.

The viscous face of the wedge and of the sector takes, with
`diffusive-metric = stencil`, the metric of the same differences of the
positions as the gradient takes of the values: across the face the image, at
the position of its ghost cell, minus the cell; along the face the differences
of the row of the cell and of the row of the image, which are those of the row
of the source turned by $R$, positions and values (a ghost cell of the row of
the source, behind a wall too, turned with it). A linear field that the
periodicity allows, $\mathbf{u}(R\mathbf{x}) = R\,\mathbf{u}(\mathbf{x})$, is so
exact at these faces on any mesh; where a wall ghost cell enters the row of the
cell or of the source, the face takes one metric per group of quantities, as an
interior face. With `mean` the face takes the inverse of the mean of the inverse
metrics of the cell and of its image, whose columns along the face are the mean
edges of the two cells, not the differences of the centres: where the cell size
grows along the face by a ratio $q$ the derivatives along the face carry the
factor $(1+q)^2/(4q)$ of the note *Tangential derivatives*
([Face metric `mean`](#face-metric-mean)), a relative error of $(q-1)^2/(4q)$
($8.3\times10^{-3}$ at $q = 1.2$), and in the first row next to an axis with
symmetry (`300`), whose ghost cell is the mirror image of the cell, $(3+q)/4$, a
relative error of $(q-1)/4$ ($5\times10^{-2}$ at $q = 1.2$); refinement at a
fixed ratio does not reduce them. For $u_r = U r$, whose viscous force is zero,
the radial viscous force of a cell of the wedge is, with the inverse-mean metric
on every face, that relative error times the hoop term
$\tau_{\theta\theta}/r$ on top of the error of order $\delta^2$ above, and with
`stencil` that error alone, $\sin^2(\delta/2)$ of the hoop term in every cell. A
viscous run on a sector whose boundary faces on the axis have zero area stops at
setup, as on the wedge, with every `diffusive-metric`: behind faces of zero area
the ghost cells copy the cells, and with `mean` the periodic faces of the first
row of cells would see half of $\partial u_r/\partial r$ ($1/6$ of the hoop term
in Homann's flow; with `stencil` the error of a small axis radius, next to the
periodic faces too). Put the axis nodes at a small radius, turned with the other
nodes of their $k$ plane.

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

Velocity, temperature and turbulence-variable gradients for the viscous flux are computed at each interface from a **10-point stencil**:

- **2 points** in the face-normal direction — the two cells sharing the face, giving the normal derivative directly;
- **8 points** in the two tangential directions — four cells per direction, each tangential derivative being the mean of the central differences of the two cells.

With $D\phi = (D_a\phi, D_b\phi, D_c\phi)$ the difference across the face and the
two tangential differences, the face gradient is $\nabla\phi_f = M_f^{T} D\phi$,
with $M_f$ the face metric tensor. The differenced quantities are the viscous
variables of the cells: the velocity, the temperature $T(p, h)$ of the fluid
tables and the turbulence variables per unit mass (`Visc_State`), computed once
per stage, after the ghost cells are filled, for every cell and for the first
layer of ghost cells (`Viscous_Variables`). Each cell stores its metric
$M = A^{-T}$ ([Metrics of the ghost cells](#metrics-of-the-ghost-cells)).
`diffusive-metric` in `[ARES-Numerics]` selects the face metric:

| `diffusive-metric` | Face metric | |
|---|---|---|
| `stencil` (default) | $M_f = \left(J_f^{T}\right)^{-1}$, $J_f^{T} = \left[D_a\mathbf{x},\ D_b\mathbf{x},\ D_c\mathbf{x}\right]$ | the same differences of the cell centres as of the field: exact for a linear field on any mesh (at the symmetry faces and the faces of type `200` for the linear fields that they allow), except at the chimera faces (`102`) and at the symmetry faces with one cell across them, which take the metric of the cell |
| (inverse-mean, not a value of the key) | $M_f = \left[\tfrac12\left(A_L + A_R\right)\right]^{-T}$ | the face-normal edge is the centre-to-centre vector: exact for a linear field on any stretching normal to the face; taken by the chimera faces (`102`) whatever the key, and with `mean` by the faces between blocks and of type `200` |
| `mean` | $M_f = \tfrac12\left(M_L + M_R\right)$ | faster; biased on stretched meshes |

The inverse-mean metric is $\left[\tfrac12\left(M_L^{-1} + M_R^{-1}\right)\right]^{-1}$ (`Face_Metric`):
the one-index difference across the face is divided by the centre-to-centre distance.

### Face metric `stencil`

With `stencil`, the columns of $J_f^{T}$ are the differences of the cell
centres (the mean of the eight nodes) that $D_a\phi$, $D_b\phi$ and $D_c\phi$
take of the field. For a linear field $\phi = \phi_0 + \mathbf{g}\cdot\mathbf{x}$
every difference is $D\phi = J_f\,\mathbf{g}$, so the face gradient is
$\mathbf{g}$ on skewed, stretched and non-smooth meshes and across folds at the
faces between blocks: the difference approximations of the field and of the
coordinates must be the same [5]. Vinokur states the rule for the gradient of a
potential (p. 26) and writes the face gradient of the viscous terms on the same
values as here, Eq. (90)–(91), with metrics from averaged areas and volumes,
Eq. (93)–(94), as in the cross-derivative terms of CFL3D [12]; `stencil`
applies the rule to that stencil. Its component along the centre-to-centre
vector is the difference of the two cell values over their distance, as in the
face-tangent and node-averaging schemes [9, 10], with the derivatives along the
stencil's tangential directions taken between means of four cells, at the mean
of their centres (near, not in, the plane of the face); on a 2-D mesh of
parallelograms it is the diamond scheme with least-squares weights [11], which
are exactly 1/4 there. Where the tangential differences of the stencil lie in
the plane of the face, its split of the face area vector is the over-relaxed
one of [13]. Exactness for linear fields does not certify accuracy on curved,
highly stretched cells [14]. One effect of the off-face position of the
tangential derivatives: on a wedge whose cells grow along $r$ by a ratio $q$,
the axial viscous force of Homann's flow converges at first order on the first
rows next to the axis, as the derivative of $u_r$ along $x$ at the faces $j$ is
that of the mean of the two rows of cells, off the face by about
$\Delta r\,(q-1)/4$. The inverse-mean metric takes the tangential columns from
the mean edges of the two cells instead. On smooth meshes they approach the
centre differences as the cells shrink, but where the mesh is not smooth (cells
of alternating size, nodes moved at random, folds at the faces between blocks)
the face gradient of a linear field stays wrong by several per cent or more at
any resolution. The column across the face is kept as $\tfrac12(a_L + a_R)$,
which is $\mathbf{x}_R - \mathbf{x}_L$, and so are the columns of the index
directions with one cell (1-D runs, the $k$ direction of 2-D meshes and
wedges). On 1-D runs and plane 2-D meshes the field does not vary along them,
and 1-D runs give the bits of the inverse-mean metric. On a wedge the $k$
neighbours are the cell turned by $\mp\delta$, so the velocity differs across
them, and that column, the mean edge of the cells, is the difference of the
turned positions divided by $\cos^2(\delta/2)$: the error of order $\delta^2$ of
the azimuthal derivatives ([Axisymmetric wedge](#axisymmetric-wedge)). The ghost
cells in the stencil take the position at which their value is right for a
linear field:

- faces between blocks (`101`, `201`): the centre of the cell across, plus the
  translation of a periodic pair;
- walls (`301`, `302`, and `103`, `104`, coupled to another solver) and `400`
  faces: the position of the class of each quantity,
  $2\mathbf{x}_f - \mathbf{x}_1$ for the quantities with a wall value and
  $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ for the extrapolated ones
  ([Ghost Cells of Walls and 400 Faces](#ghost-cells-of-walls-and-400-faces));
  a face with such a ghost cell in its stencil takes one metric per group of
  quantities (temperature, velocity, turbulence variables);
- chimera faces (`102`): $3\mathbf{x}_1 - 3\mathbf{x}_2 + \mathbf{x}_3$
  ($2\mathbf{x}_1 - \mathbf{x}_2$, or $\mathbf{x}_1$, in a block with two
  cells, or one, across the face), an approximation, since their values come
  from donor cells; the tangential differences of their faces take as
  neighbours of a ghost cell along the face the ghost cells of the neighbouring
  faces, and beyond the edges of the block face the row of its ghost cells
  extrapolated linearly (in the $k$ direction of a 2-D mesh, the ghost cell
  turned by the wedge angle); the MUSCL reconstruction of the first interior
  face takes the distance $l_0 = \Delta l_1$ of the ghost cell instead, so on
  cells stretched across these faces the slope behind the face is scaled by a
  factor of order one, which the limiter bounds;
- zero-gradient faces (`0`, `404`–`406`): the cell itself;
- symmetry (`300`): the mirror image of the cell across the face; on a face of
  zero area (an axis at $r = 0$) the cell itself, since the ghost value there is
  a copy of the cell, and the viscous flux of that face is zero;
- the $k$ faces of a wedge or of a 3-D sector (`200`): the source of the ghost
  cell, the cell at the other end of the block (on a wedge the cell itself),
  rotated by the angle of the wedge or sector; on a translated sector (no angle)
  the cell at the other end translated by the period, as across a periodic
  pair: plus the vector from the centre of its outer face to the centre of this
  face. The mesh check at setup has found every node of the last $k$ plane of
  the block to be the node of the first plane translated by one vector, within
  1 % of the move of the nodes per cell in $k$ and better than turned (the
  round-off of the coordinates, in a mesh written with enough digits), so that
  vector is the period of the block to the same tolerance.

Behind a symmetry face, a periodic pair and a rotational periodic face the
gradient is exact for the linear fields that the condition allows:
$\mathbf{g}\cdot\mathbf{n} = 0$ for a scalar and $G = RGR$, $R\mathbf{u}_0 =
\mathbf{u}_0$ for the velocity at a symmetry plane ($R = I - 2\mathbf{n}\mathbf{n}^{T}$),
$\mathbf{g}\cdot\mathbf{T} = 0$ across a periodic pair of translation $\mathbf{T}$
and across the $k$ faces of a translated sector of period $\mathbf{T}$,
$\mathbf{u}(R_\delta\mathbf{x}) = R_\delta\mathbf{u}(\mathbf{x})$ across a rotational
periodic face. For a field that breaks the condition the value of the ghost cell
(reflected, copied or rotated) is not the value of the field at any point, so no
position makes the stencil exact: the error is proportional to the part that the
condition forbids ($\mathbf{g}\cdot\mathbf{n}$, $G - RGR$, $\mathbf{g}\cdot\mathbf{T}$),
it stays on the faces next to the boundary, and it vanishes as the solution
satisfies the condition. Behind the other faces between blocks and the
zero-gradient faces the gradient is exact for every linear field; behind walls
and `400` faces for the linear fields of
[Ghost Cells of Walls and 400 Faces](#ghost-cells-of-walls-and-400-faces).

$M_f$ is computed once at setup and stored on every face of the blocks of each
rank, 27 values per cell (the mesh does not move); the cell centres come from
the nodes, which every rank holds at setup, so no message is needed. Where a
face or a cell is degenerate, the differences of its stencil can stop spanning
the space and the metric is not defined there: the run stops at setup with a
message that names the block, the cell and the face, a declared limit of
`stencil` (`mean` does not need this metric, but a viscous run with it stops on
a 3-D mesh at a boundary face of zero area, below). A wall on faces of zero
area, as on the axis of a 3-D mesh around it (a polar line), stops too: its
gradient is not defined there, and the axis takes the symmetry condition
(`300`). With `300` on the axis the polar line runs, and the face gradients next
to it are exact for linear fields. With `mean` a viscous run on a 3-D mesh with
a boundary face of zero area, as on a polar line with the axis exactly at
$r = 0$, stops at setup and names `stencil`: next to the axis the viscous force
of the first row of cells grows as the mesh is refined (a linear field, whose
force is zero).

With the axis nodes exactly at $r = 0$ the polar line runs with `stencil`, but
the solution next to the axis is of first order. On smooth fields the residual
of the two rows of cells next to the axis does not converge as the mesh is
refined: its viscous part grows as $1/h$ and its convective part does not
decrease. With `mean` the viscous residual of the cells next to a symmetry plane
(`300`) that bounds a sector of a polar mesh does not converge either: on smooth
fields that the plane allows it stays of order one as the mesh is refined, while
with `stencil` those cells converge as the others.

### Face metric `mean`

In one dimension $M = 1/h$, so `mean` divides the difference across the face by
the mean of $1/h$ instead of the centre-to-centre distance. With a cell-size
ratio $r = h_R/h_L$ it overstates the face-normal gradient, and so the viscous
and conductive fluxes, by

$$
\frac{(1+r)^2}{4r}:\qquad
+0.2\,\% \;(r=1.1),\quad +0.8\,\% \;(r=1.2),\quad +4.2\,\% \;(r=1.5),\quad +8.9\,\% \;(r=1.8)
$$

This bias does not decrease with grid refinement at a fixed growth ratio, and it
is present in the converged steady solution. `mean` is adequate where adjacent
cells grow by about 1.2 or less; boundary-layer meshes often grow faster near
the wall.

With `stencil`, a face between blocks takes one face gradient with its stencil
metric, as an interior face, and a symmetry face the gradient between the cell
and its mirror image. With `mean`, block-connection faces use the inverse-mean
metric, as chimera faces do with every option: across a large cell-size jump at
a block interface (86:1 makes it 22 times too steep) the mean makes the normal
gradient so steep that the explicit diffusion on the large cell becomes
unstable, although its own time step is fine.

!!! note "Tangential derivatives"
    With `mean` the tangential derivatives on a face are
    central differences over two cells, divided by the cell's own tangential
    edge. With stretching in a tangential direction they carry the same
    $(1+r)^2/(4r)$ factor, which `mean` does not correct. In a thin boundary
    layer this enters only through the streamwise variation of the shear on
    faces normal to the flow, so its effect is much smaller than that of the
    face-normal gradient. `stencil` divides them by the difference of the
    centres they span, and has no such factor.

### Wall faces

At a wall face (`301`, `302`) the gradient takes the difference between the
cell and the wall value across the face. With `stencil` it also takes the
differences of the wall values along the wall, between the neighbouring wall
faces of the same block face and type (central, or one-sided where the wall ends
or goes on in another block face, whose faces it does not take: exact there for
a linear wall value, first order for a curved one), with the metric of the same
differences of the positions: the vector from the face centre (the mean of its
four nodes) to the cell centre across the face, the differences of the face
centres along it. For a linear field whose wall values are its values at the
face centres the wall gradient is exact: the temperature of a `302` wall that
varies along the wall, the velocity of a no-slip wall on cells of any shape. The
wall temperatures of two neighbouring records are differenced as they are given:
a jump between them enters the derivative along the wall at the faces next to
it, of the order of the jump over the spacing of the faces, and on cells skewed
by θ their heat flux through $\tan\theta$. The difference across the face is
that of a wall of a finite-volume grid [5, Eq. (153)]; the differences along the
wall and the metric of the same positions are not taken from the references. On
a 2-D mesh the wall stencil has the structure of the boundary diamond of [11] (a
value at the midpoint of the boundary edge and values along the edge), there at
the vertices of the edge and here at the neighbouring face centres. A boundary
face that keeps only the normal distance, as in [13] (Eq. 3.57–3.59) and in
OpenFOAM (OpenFOAM-dev, commit dfa6c7e, `fvPatch.C`), has the error of `mean`
([Face metric `mean`](#face-metric-mean)). With `mean`, and for the quantities whose wall value is not known
along the wall (the temperature of `301` walls), the derivatives along the wall
are zero and the metric is that of the wall cell: on a cell skewed by θ the
normal gradient is then wrong by $(\partial\phi_w/\partial s)\tan\theta$, which
refinement does not reduce. The turbulence variables take no difference of their
wall values along the wall, with every option: their wall values are uniform
along a smooth wall (zero for $\tilde\nu$ of Spalart–Allmaras, for $k$ and for
the Reynolds stresses), but not those of $\omega$ of the $k$–$\omega$ models and
of SSG/LRR, $C\mu_l/(\rho\,d^2)$ with $d$ the distance of the cell centre from
the wall, which follow the height of the first cell, nor that of $\tilde\nu$ on a
rough wall, taken from the value of the cell. Their derivative along the wall is
left out, with the error above on skewed cells. The wall temperature of a `301`
wall and the wall values of the turbulence variables are estimated with the
metric of the cell.

### Metrics of the ghost cells

Each cell stores its metric $M = A^{-T}$, where the rows of $A$ are the cell edge
vectors $\partial\mathbf{x}/\partial\xi$, $\partial\mathbf{x}/\partial\eta$,
$\partial\mathbf{x}/\partial\zeta$ (each the mean of the four edges of the cell
along that index); the edge lengths give the physical spacings $\Delta l$ of the
MUSCL reconstruction. The ghost cells behind a boundary face take their metric
from nodes extrapolated linearly along the index lines,
$\mathbf{N}_{-g} = 2\mathbf{N}_{-g+1} - \mathbf{N}_{-g+2}$ (faces between
blocks: the cell across; symmetry: the mirror image; the wedge: the rotated
cell; a direction with one cell: a copy): both layers have the edge across the
face of the first cell, so their size across the face never vanishes nor turns
over however stretched the mesh, and the first is centred at
$2\mathbf{x}_f - \mathbf{x}_1$ ($\mathbf{x}_f$ the mean of the nodes of the face,
$\mathbf{x}_1$ that of the cell), the position of the ghost values that put a
wall value on the face. The MUSCL reconstruction of the first interior face
takes $l_0 = \Delta l_1$ there. Their edges along the face are extrapolated as
well, and turn a layer over where they change sharply across the first cell, as
where the index lines across the face spread apart away from it: on a plane 2-D
mesh where they grow 3 times or more the first layer, from 5/3 times the second
(their determinants are $(3 - r)/(1 + r)$ and $(5 - 3r)/(1 + r)$ of that of the
first cell, $r$ the ratio of the edges along the face of the node planes 1 and 0).
Behind a chimera face the first layer enters the viscous flux (below), and it
must keep the sign and, to $10^{-12}$, the size of the determinant of the first
cell, or the run stops at setup with its block, cell and face. Behind the other
faces no flux reads the metric of either layer: walls take the metric of the
cell, or the stencil of positions, the other faces carry no viscous flux that
needs it, and the wedge rotates the cell. The solver reads only the length
across the face, in the MUSCL reconstruction, and that is the length of the
first cell whatever the edges along the face: their sign is not checked, only a
layer of zero size stops the run. The symmetry faces fill their ghost layers one
at a time, each after the previous layer of every face: in a direction with
fewer cells than the layer the cell copied is a ghost cell of the opposite face,
so no face reads a ghost cell that another is writing, whatever the number of
threads.

### Chimera faces (`102`)

The ghost values of a chimera face come from the donor cells and weights that
the tool which wrote `bc.txt` computed for its own ghost cells; the face metric
of their viscous flux takes these linear ones, so where the tool builds its ghost
cells otherwise (on cells stretched across the face) the normal gradient of the
viscous flux is taken over a distance that differs by as much. The viscous flux
of a chimera face is the mean of two fluxes: one with the central differences
along the face of the row of the cell, one with those of the row of the ghost
cells, filled from the donors. With a stress linear in the gradient (laminar,
Spalart–Allmaras) the mean is the flux of the mean gradient, whose differences
along the face are the mean of those of the two rows, as on an interior face
[5, Eq. (90)–(91)]. Beyond the edges of the block face the row of the ghost cells
is extrapolated linearly, in the variables whose differences the stencil takes
(velocity, temperature and the turbulence variables per unit mass, at the
pressure of the ghost cell): the difference there is the one-sided one, exact
for a linear field (a copy of the ghost cell would halve it, and the face would
take 3/4 of the derivative of a linear field along the edge). In the $k$
direction of a 2-D mesh the neighbours of a ghost cell are the ghost cell turned
by the wedge angle, as those of a cell are; across a block one cell wide in
another direction they are the ghost cell itself. The ghost cells of a chimera
face take the wall distance and the roughness of their receiver cell (zeroth
order on purpose: a linear extrapolation can reach $y_n \le 0$ next to walls).

---

## Ghost Cells of Walls and 400 Faces

Behind a viscous wall (`301`, `302`, and the walls `103` and `104` coupled to
another solver) and an extrapolation face (`400`) the two layers of ghost cells
are filled once per stage, before the fluxes (`Ghost_Wall`), and the viscous
fluxes, the MUSCL reconstruction and the source terms of the turbulence models
read the same values. Extrapolated are the quantities whose differences the
viscous fluxes and the sources take, the temperature $T$, the velocity
$\mathbf{u}$ and the turbulence variables per unit mass $\phi = (\rho\phi)/\rho$,
and the pressure $p$; the ghost state follows from them through the fluid
tables, $h = h(p, T)$, $\rho = \rho(p, h)$, $(\rho\phi) = \rho\phi$. With $q_1$,
$q_2$ the values in the wall cell and in the next one across the face:

| Class | Quantities | First ghost cell | Position of its value |
|---|---|---|---|
| wall value | velocity at `103`, `104`, `301` and `302` (zero); the turbulence variables at the walls ($k$, $R_{ij}$, $\tilde\nu$ on a smooth wall: zero; $\omega$, $\tilde\nu$ on a rough wall: the wall value of the model) | $2q_w - q_1$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| wall temperature | $T$ at `302` | $\max(2T_w - T_1,\ T_w/2)$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| extrapolated | $p$; $T$ at `103`, `104`, `301` and `400`; the velocity and the turbulence variables at `400` | $q_1 + \psi(q_1 - q_2)$ | $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ |

$\mathbf{x}_f$ is the centre of the face (the mean of its four nodes) and
$s = l_0/l_1$ the ratio of the distances of the MUSCL reconstruction at the
first interior face, $l_0 = (\Delta l_g + \Delta l_1)/2 = \Delta l_1$ (the ghost
cell has the edge across the face of the first cell, see
[Metrics of the ghost cells](#metrics-of-the-ghost-cells)),
$l_1 = (\Delta l_1 + \Delta l_2)/2$: $s = 2\Delta l_1/(\Delta l_1 + \Delta l_2)$,
1 on uniform cells. $\psi \le s$ is the largest value that keeps the pressure,
the temperature and each turbulence variable within $[q_1/2,\ 2q_1]$ (with the
sign of $q_1$); the velocity is not limited ($\psi = s$). The wall values of the
turbulence variables are those of the wall fluxes, from the same wall state
(`302`: the state at $T_w$ and the pressure of the cell; `301`: the state of the
wall cell, as the wall flux takes it; the coupled walls `103` and `104`, whose
viscous flux the other solver gives: the state of the wall cell, as an
adiabatic wall, smooth or rough as their `ks`). Every wall takes the roughness
of the nearest wall, as its wall flux and the source terms of the cell do
(`SA-rough`; $k_s = 0$ gives the value of a smooth wall). With one cell across
the face the extrapolated quantities are copied (position $\mathbf{x}_1$). In an inviscid run the walls carry no wall
values and every quantity is extrapolated, as behind a `400` face. The second
layer takes $2q_w - q_2$, $\max(2T_w - T_2,\ T_w/2)$ and
$q_1 + 2\psi_2(q_1 - q_2)$; no flux reads it. The pressure and the temperature
of each layer are then kept within the fluid tables: $p$ within the pressures of
the tables and $T$ within the temperatures that the two rows of the tables
bracketing $p$ reach from $h_\text{min}$ to $h_\text{max}$, where $h(p, T)$ and
$\rho(p, h)$ read the tables without leaving them (where the two rows have no
temperature in common, the ghost layer takes the pressure and the enthalpy of the
wall cell; nothing else checks the tables: the state of the cells must stay
within them).

Properties, for every mesh and in 3-D, with the hypotheses written:

- **Positivity.** The ghost state has $p > 0$ and $T > 0$ (the extrapolated
  ones within a factor two of the cell, the temperature of a `302` wall at least
  $T_w/2$) and is a state of the fluid tables; the extrapolation amplifies $q_1$
  at most by $1 + s$. A quadratic
  extrapolation of the primitives, $3q_1 - 3q_2 + q_3$, amplifies by 3 and 6 on
  the two layers and can leave the tables.
- **Wall value.** For the first two classes without the floor, the value
  interpolated at the face centre, $(G_1 + q_1)/2$, is the wall value, for any
  profile (zero velocity, zero $k$, the model value of $\omega$). The floor
  $T_w/2$ acts only where $T_1 > 1.5\,T_w$, a first cell much hotter than a
  cooled wall that the mesh does not resolve; there the wall flux still takes
  $T_w$, and the ghost cell stays on the cold side of the wall.
- **Linear fields.** The rules are affine combinations of the values whose
  coefficients sum to one, so for a linear field the ghost value is the field at
  the position of the table: for the extrapolated class for every linear field,
  for the other two for the linear fields that take the wall value at the face
  centres (zero velocity on the wall, $T(\mathbf{x}_f) = T_w$ of each face).
  With these positions in the face metrics (`stencil`) the face gradients next
  to walls and `400` faces are exact for those fields, whatever the density:
  the quantities differenced are extrapolated themselves. So is the MUSCL
  reconstruction at the first interior face along a straight grid line: both
  positions are $l_0 = \Delta l_1$ from $\mathbf{x}_1$, the distance it assumes.
  The density gradient of the source terms takes the position of the
  temperature. At `103`, `104`, `301` and `400` the pressure and the temperature
  both sit at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$, so the density of
  the ghost cell is that of the fluid at the extrapolated $p$ and $T$ there. At
  `302` the temperature sits at $2\mathbf{x}_f - \mathbf{x}_1$, $\Delta l_1$ from
  $\mathbf{x}_1$ along the edge across the face of the wall cell, and the
  pressure at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$, along the line from
  $\mathbf{x}_2$ to $\mathbf{x}_1$, $\Delta l_1$ from $\mathbf{x}_1$ on a straight
  grid line: there the two positions are the same point; elsewhere they differ
  by about $\Delta l_1$ times the angle $\theta$ between the two directions, and
  the density is that of a single point where the pressure does not vary across
  the wall cell (a boundary layer).
- **Fluxes and sources.** One set of ghost cells: the wall flux routines write
  no ghost values, and the source terms take the positions of the table.

Where a limit acts the ghost value is no longer the value of the field at its
position, and the gradients next to the face are first order there: it happens
only where a quantity changes by more than a factor two between the first two
cells, where the first cell is 1.5 times hotter than a `302` wall, or where the
ghost state would leave the fluid tables. At the end of a run in which any of
them acted, one line reports the ghost fills in which each did (some of the
counts can be zero); a run in which none acted prints no such line:

```
 Ghost cells of walls and 400 faces, fills with the extrapolation cut: <n>, with the 302 temperature at T_w/2: <m>, brought within the fluid tables: <t>
```

Declared limits:

- the coupled walls `103` and `104` take the wall values of the turbulence
  variables of an adiabatic wall at the state of the cell: exact for $k$,
  $R_{ij}$ and $\tilde\nu$ on a smooth wall (zero); on a rough wall
  (`SA-rough`, $k_s > 0$) $\tilde\nu$ takes the value of the rough-wall
  condition from the state of the cell, as the wall flux of the other solver
  does, per unit mass with the density of the cell where that flux takes the
  density at the interface temperature; that of $\omega$, $C\mu/(\rho\,d^2)$,
  differs from the value at the state of the coupled surface by the ratio
  $\nu(T_1)/\nu(T_s)$ of the kinematic viscosities at the two temperatures: of
  order one but bounded, it enters only the source terms (the gradient of
  $\omega$ in the cross-diffusion of SST, $\nabla k\cdot\nabla\omega$) and the
  limited MUSCL reconstruction of the first face; their viscous flux comes from
  the other solver;
- the temperature of a `301` wall is extrapolated (consistent, exact for every
  linear field) rather than set from the heat flux.

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
5. M. Vinokur, "An analysis of finite-difference and finite-volume formulations of conservation laws," NASA CR-177416, 1986: pp. 15–16, Eq. (60a) (the area vector of a face of a hexahedron as half the cross product of its diagonals; the shape of each face "precisely defined and consistently used by the two neighboring cells" so that the volumes sum to the total volume; the doubly ruled face); p. 23, Eq. (90)–(91) (the transport terms on an auxiliary cell, with means of four cells along the face); p. 26, Eq. (101) ("the two difference approximations must be the same"); p. 35, Eq. (137)–(138) (first-order one-sided expressions for the geometry at a boundary, as second-order ones can give negative volumes); p. 40, Eq. (153) (the difference across a wall). Journal version: *J. Comput. Phys.*, 81, 1989.
6. R. Vandenhoeck, A. Lani, "Implicit high-order flux reconstruction solver for high-speed compressible flows," arXiv:1811.08666v1, 2018, Appendix A.2, eqs. (A.6)–(A.7) (axisymmetric source $(0, 0, p - \tau_{\theta\theta}, 0)^T$).
7. Eilmer 4, Gas Dynamics Toolkit, `src/eilmer/fvcell.d`, lines 1213–1228 (viscous source of axisymmetric flow), https://github.com/gdtk-uq/gdtk, commit 826496a; the code is described in N. N. Gibbons, K. A. Damm, P. A. Jacobs, R. J. Gollan, "Eilmer: an open-source multi-physics hypersonic flow solver," *Comput. Phys. Commun.*, 282, 108551, 2023, doi:10.1016/j.cpc.2022.108551.
8. F. Homann, "Der Einfluss grosser Zähigkeit bei der Strömung um den Zylinder und um die Kugel," *Z. Angew. Math. Mech.*, 16, 153–164, 1936 (not read; as cited by J. Magnaudet and M. Abbas, arXiv:2011.10755v1, 2020). The equation of the flow, f''' + 2ff'' − f'² + 1 = 0 with f(0) = f'(0) = 0, f'(∞) = 1, u_r = a r f'(η), u_x = −2(aν)^{1/2} f(η), as written by P. Rajamanickam, arXiv:2508.17982v4, Eq. (12), and by H. Guo, X. Zhuang, X. Meng, T. Rabczuk, arXiv:2010.05682v1, Eq. (7) with β₀ = 2, β = 1, whose Table 8 gives f''(0) = 1.31194.
9. B. Diskin, J. L. Thomas, E. J. Nielsen, H. Nishikawa, J. A. White, "Comparison of Node-Centered and Cell-Centered Unstructured Finite-Volume Discretizations: Viscous Fluxes," *AIAA J.*, 48(7), 1326–1338, 2010, doi:10.2514/1.44940, Sec. III.A.
10. H. Nishikawa, "Beyond Interface Gradient: A General Principle for Constructing Diffusion Schemes," AIAA Paper 2010-5093, 2010, App. C.
11. Y. Coudière, J.-P. Vila, P. Villedieu, "Convergence rate of a finite volume scheme for a two dimensional convection-diffusion problem," *ESAIM: M2AN*, 33(3), 493–516, 1999, Eq. (9) and Lemma 6.2.
12. CFL3D, `source/cfl3d/libs/gfluxv1.F` (cross-derivative terms of the viscous fluxes), https://github.com/nasa/CFL3D, commit 776d516; S. L. Krist, R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual (Version 5.0)*, NASA/TM-1998-208444, 1998, App. C.2, p. 249 (the volume at a cell interface as the mean of its neighbours, for the thin-layer terms; the manual, p. 27, leaves out the cross-derivative terms, which gfluxv1.F adds).
13. H. Jasak, *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*, PhD thesis, Imperial College, University of London, 1996, Eq. (3.1), (3.12) and (3.32).
14. D. J. Mavriplis, "Revisiting the Least-squares Procedure for Gradient Reconstruction on Unstructured Meshes," NASA/CR-2003-212683, 2003, p. 14.
