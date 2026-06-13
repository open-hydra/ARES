# Governing Equations

ARES solves the compressible Navier–Stokes equations for a **single real fluid** in conservative form on structured multi-block grids. Both the inviscid (Euler) and viscous (Navier–Stokes) limits are supported, selected by `simulation-type`.

The governing system is

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla\!\cdot\!\mathbf{F}^{c}(\mathbf{U})
= \nabla\!\cdot\!\mathbf{F}^{v}(\mathbf{U})
+ \mathbf{S}
$$

where $\mathbf{U}$ is the vector of conservative variables, $\mathbf{F}^{c}$ the convective (inviscid) flux, $\mathbf{F}^{v}$ the viscous (diffusive) flux, and $\mathbf{S}$ the source-term vector (turbulence model, optional rotating frame).

Unlike an ideal-gas solver, ARES does **not** close the system with $p=\rho R T$. The thermodynamic state is closed by a **tabulated real-fluid equation of state** expressed on a pressure–enthalpy grid (see [Real-Fluid Thermodynamics](thermo.md)).

---

## State Vectors

ARES uses two representations: a **primitive** set built on pressure and enthalpy (used for reconstruction, boundary conditions, and I/O) and a **conservative** set used for the time update and flux balance.

### Primitive variables

$$
\mathbf{P} = [\,p,\; u,\; v,\; w,\; h,\; \mathbf{q}_\text{rans}\,]
$$

| Component | Symbol | Description |
|:---------:|:------:|-------------|
| Pressure | $p$ | Static pressure |
| Velocity | $u, v, w$ | Cartesian velocity components |
| Enthalpy | $h$ | Specific static enthalpy |
| RANS variables | $\mathbf{q}_\text{rans}$ | Turbulence-model variables (optional) |

The pair $(p, h)$ is the thermodynamic state: every other property (density $\rho$, temperature $T$, speed of sound $a$, specific heats) is obtained by interpolating the real-fluid table at $(p,h)$.

**RANS variables (if active):**

- **Spalart–Allmaras**: $\tilde\nu$
- **$k$–$\omega$ (SST / Wilcox 2006)**: $k$, $\omega$
- **SSG-LRR**: the Reynolds-stress components $R_{ij}$ and $\omega$

### Conservative variables

$$
\mathbf{U} = [\,\rho,\; \rho u,\; \rho v,\; \rho w,\; \rho E_0,\; \mathbf{q}_\text{rans}\,]
$$

with

$$
E_0 = e + \tfrac{1}{2}\lvert\mathbf{u}\rvert^2,
\qquad
e = h - \frac{p}{\rho}
$$

where $\rho = \rho(p,h)$ is read from the table.

### Transformation between primitive and conservative variables

**Forward ($\mathbf{P} \to \mathbf{U}$).** With $(p,h)$ known, $\rho(p,h)$ comes from the table, momentum is $\rho\mathbf{u}$, and the total energy follows from $e = h - p/\rho$ plus the kinetic energy.

**Inverse ($\mathbf{U} \to \mathbf{P}$).** Recovering $(p,h)$ from $(\rho, \rho\mathbf{u}, \rho E_0)$ is the *thermo inversion* reported during loading. The velocity is $\mathbf{u} = (\rho\mathbf{u})/\rho$, the static specific energy is $e = E_0 - \tfrac12|\mathbf{u}|^2$, and $(p,h)$ are recovered from the pair $(\rho, e)$ by inverting the real-fluid table (a local two-variable lookup/iteration). Because the table stores $(\rho, T, a, \dots)$ as functions of $(p,h)$, this inversion replaces the Newton iteration on temperature used by an ideal-gas solver.

!!! note "Why pressure-based primitives?"
    Carrying $p$ (rather than $\rho$) as the primitive thermodynamic variable is what makes the [low-Mach preconditioned](preconditioning.md) update natural: the acoustic stiffness as $M\to0$ is handled cleanly in pressure–velocity form.

---

## Convective Fluxes

The convective flux through a face of outward unit normal $\hat{\mathbf{n}}$ is

$$
\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}} =
\begin{bmatrix}
  \rho\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}}) \\
  \rho\,u\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_x \\
  \rho\,v\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_y \\
  \rho\,w\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_z \\
  \rho\,H_0\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}})
\end{bmatrix}
$$

where $H_0 = h + \tfrac12\lvert\mathbf{u}\rvert^2$ is the total specific enthalpy and $\dot m = \rho\,(\mathbf{u}\!\cdot\!\hat{\mathbf{n}})$ is the normal mass flux. The numerical convective flux is evaluated by a [Riemann solver](riemann-solvers.md) from the reconstructed left/right states.

---

## Viscous Fluxes

### Stress tensor

The Newtonian stress tensor with the Stokes hypothesis $\lambda = -\tfrac23\mu$:

$$
\tau_{ij} =
\mu\!\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)
- \frac{2}{3}\,\mu\,(\nabla\!\cdot\!\mathbf{u})\,\delta_{ij}
$$

where $\mu = \mu_\ell + \mu_t$ is the sum of the molecular and eddy viscosities. For two-equation models the isotropic Reynolds-stress part $-\tfrac23\rho k\,\delta_{ij}$ is added.

### Heat conduction

$$
\mathbf{q} = -\kappa\,\nabla T,
\qquad
\kappa = k_\ell + \frac{\mu_t\,c_p}{\mathrm{Pr}_t}
$$

The molecular conductivity $k_\ell(p,h)$ and $c_p(p,h)$ come from the real-fluid table; $\mathrm{Pr}_t$ is the turbulent Prandtl number from `[ARES-RANS]`.

### Diffusive flux vector

$$
\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}} =
\begin{bmatrix}
  0 \\
  \boldsymbol{\tau}\!\cdot\!\hat{\mathbf{n}} \\
  (\boldsymbol{\tau}\!\cdot\!\mathbf{u})\!\cdot\!\hat{\mathbf{n}} + \kappa\,\nabla T\!\cdot\!\hat{\mathbf{n}}
\end{bmatrix}
$$

i.e. viscous work plus heat conduction in the energy equation. Gradients are evaluated at faces with a wide stencil mapped to Cartesian coordinates by the face-metric tensor (see [Spatial Discretization](numerics.md)).

---

## Source Terms

| Physical process | Affected variables | Treatment |
|------------------|--------------------|-----------|
| Turbulence | RANS variables ($\tilde\nu$; $k,\omega$; $R_{ij},\omega$) | Production/destruction added to the spatial residual; wall corrections |
| Rotating frame *(optional)* | Momentum, energy | Coriolis and centrifugal source terms on the relative velocity |

ARES has no chemical-source term: it solves a single fluid whose properties are fixed by the real-fluid table, so there is no species transport or finite-rate kinetics.
