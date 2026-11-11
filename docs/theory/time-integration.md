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

### Unified coefficient form

The three schemes share a single kernel. At stage $k$ the update is

$$
\mathbf U^{(k)} = \mathbf U^n + c_{s,k}\bigl(\mathbf U^{(k-1)} - \mathbf U^n + \mathbf R^\ast\bigr),
$$

where $\mathbf R^\ast$ is the scaled residual of stage $k$ and the stage coefficients $c_{s,k}$ are stored in the `RKcoeff` matrix

$$
\mathbf C =
\begin{pmatrix}
1 & 0 & 0 \\[2pt]
1 & \tfrac12 & 0 \\[2pt]
1 & \tfrac14 & \tfrac23
\end{pmatrix}.
$$

The row index $s$ is the scheme order (1 = `euler`, 2 = `RK2`, 3 = `RK3`) and the column index $k$ is the stage, so a single coded loop realises all three integrators.

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

### Local time step

Each cell $i$ gets a local time step that combines convection and diffusion in one additive rule
([5], §6 and §9.3; [4], eq. 26–28):

$$
\frac{1}{\Delta t_i} =
\frac{1}{\text{CFL}} \sum_{d=1}^{n_d} \Bigl( \bigl|\mathbf{v}\cdot\nabla\xi_d\bigr| + a\,\bigl|\nabla\xi_d\bigr| \Bigr)
\;+\; \frac{D_{\max}}{\text{VNN}} \sum_{d=1}^{n_d} \bigl|\nabla\xi_d\bigr|^2 ,
$$

where $\nabla\xi_d$ is the row $d$ of the cell metric (its modulus is the inverse of the width of the
cell across its faces of direction $d$), $n_d$ the number of computational directions (one cell in
$k$ gives $n_d = 2$), $a = a(p,h)$ the real-fluid speed of sound and $D_{\max}$ the largest diffusivity of
the equations:

$$
D_{\max} = \max\!\left(
\frac{4}{3}\,\frac{\mu_\ell + \mu_t}{\rho},\;
\frac{k_\ell + \mu_t\,c_p/\text{Pr}_t}{\rho\,c_v},\;
\frac{\mu_\ell + \rho\tilde\nu}{\sigma\,\rho}
\right).
$$

The conductivity acts on the internal energy $\rho c_v T$, so the energy diffusivity is
$\gamma\nu/\text{Pr} \approx 1.94\,\nu$ in air, larger than the momentum one. The specific heat at
constant volume is that of the real fluid, $c_v = c_p/\gamma$ with
$\gamma = a^2\,(\partial\rho/\partial p)_T$ and
$(\partial\rho/\partial p)_T = (\partial\rho/\partial p)_h + (\partial\rho/\partial h)_p\,
\bigl[1/\rho + T\,(\partial\rho/\partial T)_p/\rho^2\bigr]$, all from the $(p,h)$ table
([Thermodynamics](thermo.md)); for a perfect gas it is $c_p - R$. The last term enters with
Spalart–Allmaras ($\sigma = 2/3$). The $k$ and $\omega$ diffusivities of SST and Wilcox 2006,
$(\mu_\ell + \sigma\mu_t)/\rho$ with $\sigma \le 1$, stay below the momentum one. The eddy viscosity
is that of the state of the cell at zero strain. Euler runs have no diffusive term.

In 3-D, on a cell whose three directions are skewed against each other, the viscous rate is
multiplied by $F = (2 + m + 1/m)/4$: with $a$, $b$, $c$ the cosines of the angles between the three
$\nabla\xi_d$, where $a^2 + b^2 + c^2 + 2abc > 1$, $m + 1$ is the largest eigenvalue of the
$3\times3$ matrix of the cosines ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$).
There the spectral radius of the viscous operator exceeds $4D_{\max}\sum_d|\nabla\xi_d|^2$ by a
factor below $9/8$, which $F$ bounds (exactly on rhombohedra); elsewhere, orthogonal cells and cells
with one direction orthogonal to the other two included, $F = 1$, and in 2-D the rate holds as it
is.

Without dimensional splitting the convective limits of the directions add up: an explicit scheme
on square 2-D cells needs a condition of the type $\tau(|a_x|/\Delta x + |a_y|/\Delta y) \le C_0$, not
the largest of the two terms ([5], §9.3, which states the form without a proof). For forward
Euler the form follows, with the diffusive rates added as well, from the positivity condition (6.1)
of [5]: with first-order upwind convection and central diffusion on a Cartesian mesh the
diagonal entry of a row is $-\sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr)$ and the others are not negative, so
$\tau \sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr) \le 1$: the rule above with CFL = 1 and VNN = 1/2. The
constants of Heun and RK3 in the table do not come from that source: they are the von Neumann limits
of the linearised MUSCL–upwind advection with central diffusion, proven for this stencil in 1-D and
carried to 2-D and 3-D, for any split of the Courant number, by the convexity of the symbol; they hold
for the van Leer slope, van Albada and MC, not for minmod and superbee. The source warns that such an
eigenvalue criterion can mislead for the non-normal operators that boundaries give, as an inflow
boundary with upwind differences ([5], §3.3, pp. 26–27): the limits are those of the interior.

??? note "Proof of the limits"
    The
    van Leer slope (van Albada and MC linearise the same way) linearised in smooth flow is
    the mean of the two one-sided slopes (Fromm's scheme), with the symbol
    $g(\theta) = \tfrac12[(1-\cos\theta)^2 + i\sin\theta(3-\cos\theta)]$, $z = -\nu g$; the constants
    hold for these three limiters, not for minmod and superbee, which take the smaller or the larger of
    the two slopes and have no linearisation where they are equal. Heun is stable if and only if
    $\nu \le 1$: $|R(z)|^2 - 1 = \tfrac{\nu}{16}(1-c)^2[\nu(\nu a - 2b)^2 + 4\nu b^2 - 16]$ with
    $c = \cos\theta$, $a = 5 - c^2$, $b = 1 - c$, increasing in $\nu$ and equal to
    $\tfrac{\nu}{16}(1-c)^4(1+c)(c-3)$ at $\nu = 1$. RK3 is stable if and only if
    $\nu \le \nu^\ast = 1.1757678815978$, a root of a polynomial of degree 24 (a cylindrical
    decomposition in exact arithmetic), and with the first-order flux if and only if
    $\nu \le 1.2563726633$, half the real stability interval $2.5127453$; diffusion gives VNN $\le 1/2$
    and $0.6282$ (the table rounds them down). The curve of $-g$ bounds a convex set $K$
    ($x'y'' - y'x'' = -\tfrac12(1-\cos\theta)(\cos^2\theta + \cos\theta + 4) \le 0$ along it), so a sum
    over directions is a point of $\nu K$, inside the stability region by the maximum principle: the
    limits of 1-D hold in 2-D and 3-D for any split of the Courant number, with no need for a convex
    stability region. With diffusion the additive rule is stable for Heun with CFL = 1 and VNN = 1/2,
    and exact at every Péclet number (the checkerboard mode reaches $z = -2$); for RK3 the region is not
    convex, and the convex hull of $1.17K$ and $[-2.48, 0]$ (and of $1.1757K$ and $[-2.512, 0]$) is
    shown inside it with a certificate in exact rational arithmetic.

| Scheme | real stability interval | CFL ≤, MUSCL, smooth flow | CFL ≤, first order (limiter at extrema and shocks) | VNN ≤ | CFL and VNN ≤ for positivity with the limiter |
|---|---|---|---|---|---|
| Heun (`RK2`) | $[-2, 0]$ | 1 | 1 | 1/2 | 1/2 |
| `RK3` | $[-2.5127, 0]$ | 1.17 | 1.25 | 0.62 | 1/2 |
| explicit Euler (`euler`) | $[-2, 0]$ | unstable | 1 | 1/2 | 1/2 |

On the 2-D axisymmetric wedge (type `200` or `300` on the $k$ faces) the stress on the wedge faces adds to the
rows of $u_r$ and $u_\theta$ of every cell terms that grow as $1/r^2$ next to the axis: in the row of
$u_r$ the diagonal $-\tfrac{4}{3}\nu_e/r^2$ and couplings $\tfrac{2}{3}\nu_e|\nabla\xi_d\cdot\mathbf{e}_r|/(2r)$
to the $u_r$ and $\tfrac{2}{3}\nu_e|\nabla\xi_d\cdot\mathbf{e}_x|/(2r)$ to the $u_x$ of the neighbours $d\pm1$, in the
row of $u_\theta$ the diagonal $-\nu_e/r^2$ and couplings $\nu_e|\nabla\xi_d\cdot\mathbf{e}_r|/(2r)$, with
$\nu_e = (\mu_\ell + \mu_t)/\rho$ and $r$ the distance of the cell centre from the axis (no coupling in the row of
$u_x$). The diffusive rate takes the increase of the row sum they give,

$$
s_h = \nu_e \max\!\left( \frac{4/3}{r^2} + \frac{2}{3r}\sum_d \bigl(|\nabla\xi_d\cdot\mathbf{e}_x| + |\nabla\xi_d\cdot\mathbf{e}_r|\bigr),\;
\frac{1}{r^2} + \frac{1}{r}\sum_d |\nabla\xi_d\cdot\mathbf{e}_r| \right),
\qquad
\frac{1}{\Delta t_i} \mathrel{+}= \frac{s_h}{4\,\text{VNN}} ,
$$

so that $\Delta t\,(4D_{\max}\sum_d|\nabla\xi_d|^2 + s_h) = 4\,\text{VNN}$. The limits of the table are proven away
from the axis on uniform meshes (in 2-D for any skew and aspect ratio, in 3-D with the factor $F$ above) and, next to the
axis, for the scalars; for the momentum next to the axis they rest on evidence.

??? note "What supports the limits of the table"
    - Away from the axis, on a uniform mesh of equal cells (parallelograms in 2-D, parallelepipeds in 3-D) with frozen
      coefficients and no boundaries, the viscous operator has a real symmetric Fourier symbol. With the face gradient (the
      difference across the face, the mean of the central differences of the two cells along it, mapped with the face
      metric, which on such a mesh is the cell metric), a mode $e^{i\boldsymbol\theta\cdot\mathbf c}$
      gives $\rho\,\dot{\hat{\mathbf u}} = -S\hat{\mathbf u}$ with $S = \mu\,\mathrm{tr}(B)\,I + (\mu+\lambda)B$, $B = G\Theta G^T$,
      the columns of $G$ the rows $\nabla\xi_d$ of the cell metric, $\Theta_{dd} = 4\sin^2(\theta_d/2)$,
      $\Theta_{de} = \sin\theta_d\sin\theta_e$; a scalar of diffusivity $D$ has the symbol $D\,\mathrm{tr}(\Theta\Gamma)$,
      $\Gamma_{de} = \nabla\xi_d\cdot\nabla\xi_e$ (the cross metric terms of skewed cells included).
      $\Theta = \mathbf c\mathbf c^T + \mathrm{diag}\bigl((1-\cos\theta_d)^2\bigr)$, $c_d = \sin\theta_d$, is positive
      semidefinite: the eigenvalues are real and not positive, the operator is normal, and the eigenvalue criterion is von
      Neumann's ([5], §2.4, eq. (2.10), p. 18, and §3.3, pp. 25–28).
    - In 2-D $\Theta \le 4I$: for any skew angle and aspect ratio the spectral radius is exactly $4D\sum_d|\nabla\xi_d|^2$ for a
      scalar and $4\nu\bigl(\sum_d|\nabla\xi_d|^2 + \lambda_{\max}(\Gamma)/3\bigr) \le \tfrac{16}{3}\nu\sum_d|\nabla\xi_d|^2$ for the
      momentum, so the rate $4D_{\max}\sum_d|\nabla\xi_d|^2$ bounds it, with equality for the scalar of largest diffusivity (the
      energy in air). In 3-D the sharp constant is $9/2$: with $a$, $b$, $c$ the cosines of the angles between the three
      $\nabla\xi_d$, the rate bounds the spectral radius if and only if $a^2 + b^2 + c^2 + 2abc \le 1$ (true for orthogonal
      cells, for cells with one direction orthogonal to the other two, and for a rhombohedron up to the edge angle
      $\arccos(-1/3) = 109.47^\circ$); elsewhere the radius exceeds it by a factor below $9/8$, and the local step multiplies the
      rate by $F = (2 + m + 1/m)/4$, where $m + 1 \in (2, 3)$ is the largest eigenvalue of the $3\times3$ matrix of the cosines
      ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$): a proven bound of the factor, exact on rhombohedra
      (1.04 at $115^\circ$, 1.11 at $119.5^\circ$), below $9/8$ on every cell, and 1 where the rate already holds. The real
      stability intervals of Heun, $[-2, 0]$, and of RK3, $[-2.5127, 0]$, then give VNN = 1/2 and 0.62 ($2.5127/4 = 0.628$).
    - Next to the axis the operator is not translation invariant. For the scalars (energy, turbulence variables) it stays
      self-adjoint in the inner product weighted by the cell volumes: the ghost cells of the wedge copy the scalars and the wedge
      faces do not update them, so on a mesh of equal rectangular cells in $(x, r)$ each face couples two cells symmetrically, and
      the row sums are at most $4D\sum_d|\nabla\xi_d|^2$ (up to a relative $\varepsilon/6$ for an axis at $\varepsilon\,\Delta r$,
      covered by $s_h$ for diffusivities below $0.6\,\nu_e/\varepsilon$). A self-adjoint operator whose row sums are within the
      rate of each cell is stable with the local step, in the norm weighted by $V/\Delta t$, and the limits of the table hold for
      the scalars on every row. The momentum operator is not self-adjoint in any inner product with cell weights: the cross terms
      $(\lambda+\mu)\,\partial^2/\partial x\,\partial r$ of the faces in $x$ and in $r$ are weighted by the radius of the cell,
      $(j-\tfrac12)\Delta r$, and of the face, $j\Delta r$. In the first row the rows of $u_r$ and $u_\theta$ are diagonally
      dominant (their diagonals carry the wedge terms $-\tfrac43\nu_e/r^2$ and $-\nu_e/r^2$); the row of $u_x$ is not: next to the
      axis (doubled weight of the upper face, mirror image across the axis) it couples to $u_r$ of four neighbours with
      $\nu/(3\Delta x\,\Delta r)$ each, four times the cross coupling of the plane operator, a deficit of
      $\tfrac43\nu/(\Delta x\,\Delta r)$; its row sum stays below the rate of the cell. No proof covers the momentum next to the
      axis. The support is evidence: on the discrete operator for cells with $dx/dr$ from 0.01 to 128 and wedges of 0.25° to 5°,
      its symmetric part is negative semidefinite, its eigenvalues are real to $6\cdot10^{-4}$, and the largest stable VNN with the
      local step is 0.500 to 0.577 with Heun and 0.628 to 0.725 with RK3, approaching the limits of the table, 1/2 and 0.62
      (0.628), from above for very long or very short cells; on a diffusion-dominated 1° wedge of 32 × 32 cells the largest stable
      VNN over 2000 local steps is about 0.5 with Heun and 0.63 to 0.65 with RK3 for $dx/dr$ = 1, 4 and 8. The eigenvalue
      criterion can mislead for non-normal operators ([5], pp. 26–27).

Next to the axis of square cells $s_h = 8\,\nu/\Delta r^2$; far from it $s_h$ vanishes as $1/r$. The wedge stress, and with
it $s_h$, applies with either sign of the wedge angle (a right-handed mesh with the axis on the last $j$ face has a negative
one). The step takes $s_h$ on every 2-D mesh with a wedge angle, whatever the type of its $k$ faces: the viscous flux of
faces of type `200` and that of symmetry faces (`300`), whose mirror image of a cell is, without swirl, the cell turned by
the wedge angle, both carry this stress; null (`0`) $k$ faces carry none, and there $s_h$ only makes the step smaller than
it needs to be.

These limits are linear (von Neumann, smooth solutions). With the limiter, positivity (a maximum principle:
in more than one dimension a second-order scheme is not TVD, Goodman et al. 1985 as cited in [6], p. 15 of the
report) is guaranteed only for smaller steps, for the
scalar advection of each direction with the limited flux: the van Leer slope is the flux limiter $\psi(\theta) =
\theta/(1+\theta)$, which satisfies $0 \le \psi \le 1$ and $0 \le \psi/\theta \le \mu$ with $\mu = 1$, and with
Heun and RK3 (SSP coefficient 1) positivity holds for $(1+\mu)\,\nu \le 1$, that is
$\sum_d \nu_d \le 1/2$; central diffusion is positive for $\text{VNN} \le 1/2$; with the additive rule
both hold for CFL ≤ 1/2 and VNN ≤ 1/2 ([5], §5.2, §6.2 and §6.4; experiments reported there give
positive results up to $\nu = 1$ in 1D and about $2/3$ in 2D). With moving shocks and Heun keep `cfl` at
0.5–0.7. With preconditioning,
$|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$ is replaced by the preconditioned eigenvalue
$\bigl(u_n (1-\alpha) + c^\ast\bigr)|\nabla\xi_d|$, with $c^\ast = \sqrt{\alpha^2 u_n^2 + U_r^2}$ and
$\alpha = (1 - \beta U_r^2)/2$ ([Low-Mach Preconditioning](preconditioning.md)).

The local time step takes the wave speeds of the cell alone. With a stationary shock captured on one face, the waves
that the HLLC flux lets into a cell on the slow side of the shock take the speeds of the Roe state of the face, which
exceed the cell's own $|u_n| + a$: there the Courant number of the face exceeds the CFL number. With strong shocks keep
CFL ≤ 0.8 for RK3. The local step changes the path of a steady run, not the steady solution it converges to.

**Recommended values.** For time-accurate runs with moving shocks and RANS models: with RK3, `cfl = 0.8` and
`vnn = 0.5`, 20 % below the limit 0.628 of the scheme (up to that limit, for example at `vnn = 0.6`, 4.5 % below it, the
linearised scheme is stable); with Heun, `cfl` from 0.5 to 0.7 and `vnn = 0.4`, below its limit 0.5. These are values
of the scheme, not of a case: a case that runs stable beyond them owes it to its mesh (a non-uniform spacing around the
cell that sets the step lowers the spectral radius of its diffusion), and its threshold is not a limit to carry
elsewhere. The viscous rate takes the largest diffusivity, that of the energy, and adds up the directions: a `vnn` tuned
on a rule that took the momentum diffusivity and the largest direction alone gives here a smaller step (by
$\gamma/\text{Pr}_t$, 1.56 for air with the default $\text{Pr}_t = 0.9$, in a stretched cell where the turbulent
viscosity sets it), and can be raised within the limits of the table. The CFL number of a run with local time steps,
as a steady run, is a parameter of the pseudo-time, not a Courant number of the physical time: it is not a value to
carry to a time-accurate run.

### Stiffness of the turbulence sources (time-accurate RANS)

In time-accurate RANS runs the turbulence sources are explicit, and each local step also includes
a bound $\hat d_i$ of their stiffness, the derivative $d_i$ of the destruction terms, from the state of the cell at
the start of the step alone: $1/\Delta t_i$ of the rule above gains the term $\hat d_i$,
so that $\Delta t\,d_i \le \Delta t\,\hat d_i \le 1$ and no stage makes $\tilde\nu$, $k$ or $\omega$ negative through
their destruction ([Turbulence](turbulence.md#time-accurate-runs-explicit-sources)). The Reynolds-stress model SSG/LRR
has no such bound: its $\hat d_i$ is zero. Steady runs take no such term: there the stiffness is handled, when
`point-implicit = .true.`, by the point-implicit factor
([Turbulence](turbulence.md#steady-runs-point-implicit-factor)).

### Global time step

In a time-accurate run all cells advance with the smallest local step:

$$
\Delta t = \min_i \Delta t_i .
$$

!!! tip "CFL ramp-up"
    `cfl-rise-threshold` linearly ramps the CFL from a low starting value to the target over the given number of iterations, improving robustness when the initial condition is far from steady state. With SST or Wilcox 2006, a time-accurate start far from the wall value of $\omega$ (100 times or more below it, in wall cells where diffusion sets the step) needs it: without the ramp the first step can make $\omega$ negative ([Turbulence](turbulence.md#time-accurate-runs-explicit-sources)).

---

## Implicit Residual Smoothing (IRS)

IRS raises the effective CFL limit by smoothing the residual with a Laplacian operator before the update. In each direction $d$ the smoothed residual satisfies

$$
\mathbf R^\ast_i = \frac{\mathbf R_i + \varepsilon(\mathbf R^\ast_{i-1} + \mathbf R^\ast_{i+1})}{1 + 2\varepsilon},
$$

approximated with a couple of Jacobi sweeps. The stable CFL limit is amplified by roughly $1+2\varepsilon$; typical $\varepsilon$ (`irs-beta`) is 0.1–0.5. Enable it with `irs = .true.`.

### Steady runs only

The smoothed residual is not the residual of the space discretization, so a run with IRS does not follow the
transient of the discretized equations: IRS accelerates the convergence of steady (pseudo-time) runs. It is switched
on by `irs = true` or by `irs-beta > 0` alone, and ARES refuses it at set-up when `time-accurate = true`. The change
does not vanish as $\Delta t \to 0$: the smoothing turns the semi-discrete system into
$(I - \varepsilon\,\delta^2)\,d\mathbf{U}/dt = -\mathbf{R}/V$ in each direction, which slows a Fourier mode of
wavenumber $k$ by $1 + 4\varepsilon\sin^2(k\Delta x/2)$: a change of order $\varepsilon (k\Delta x)^2$ for resolved
waves and of order one at discontinuities. IRS can only be switched on at set-up, where its work arrays are allocated:
a reading of `input.ini` during the run that turns it on is ignored, with a warning.

### Smoothed residuals

`irs-variables` chooses the residuals that are smoothed: `all`, the default, smooths
every residual, those of the turbulence equations too; `flow` smooths those of the
mean flow only (mass, momentum, energy), and the turbulence equations then advance
with the same local time step without smoothing: above their explicit limit they
become the unstable part. `all` is meant for Spalart–Allmaras: the k-omega models
can diverge with it, and with them `flow` is the choice. Smoothing changes the
path to the steady state, not the
steady state itself (a zero residual stays zero). On a run started with IRS,
`irs-beta` and `irs-variables` can change at runtime, and `irs = false` switches it
off.

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

**Interpolation weights.** In 3-D each fine cell is built from its parent coarse cell plus the three face-, three edge- and one corner-neighbour selected toward its octant:

$$
a_1 = \tfrac{27}{64},\quad
a_2 = \tfrac{9}{64},\quad
a_3 = \tfrac{3}{64},\quad
a_4 = \tfrac{1}{64},
$$

applied as $1\cdot a_1 + 3\cdot a_2 + 3\cdot a_3 + 1\cdot a_4 = 1$. In 2-D ($n_k=1$) the stencil collapses to the parent, two face-neighbours and one corner-neighbour:

$$
a_1 = \tfrac{9}{16},\quad
a_2 = \tfrac{3}{16},\quad
a_3 = \tfrac{1}{16},
$$

with $a_1 + 2\,a_2 + a_3 = 1$. The weights are normalised by construction, so the prolongation is conservative in the mean.

### Cycle

A V-cycle (pre-smooth → restrict → coarse solve → prolong → post-smooth) is the default; the driver can recurse for W-cycles. Multigrid can also accelerate each physical step of a transient run.

---

## References

1. C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.* 77 (1988).
2. S. Gottlieb, C.-W. Shu, E. Tadmor, "Strong stability-preserving high-order time discretization methods," *SIAM Rev.* 43 (2001).
3. A. Jameson, "Solution of the Euler equations for two-dimensional transonic flow by a multigrid method," *Appl. Math. Comput.* 13 (1983).
4. A. Arnone, M.-S. Liou, L. A. Povinelli, "Multigrid time-accurate integration of Navier–Stokes equations," NASA TM-106373 (ICOMP-93-37), 1993.
5. W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam, 2000.
6. E. M. Constantinescu, A. Sandu, "Multirate timestepping methods for hyperbolic conservation laws," *J. Sci. Comput.*, 33(3), 239–278, 2007, doi:10.1007/s10915-007-9151-y; the page cited is that of the authors' report, Virginia Tech, https://hdl.handle.net/10919/19488.
