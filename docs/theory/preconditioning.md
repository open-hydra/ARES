# Low-Mach Preconditioning

## Why preconditioning?

ARES is a **density-based** (compressible) solver: it marches the coupled mass–momentum–energy system with explicit time stepping bounded by the acoustic CFL limit. This works well at transonic and supersonic speeds, but it degrades badly as the Mach number drops toward zero, for two coupled reasons:

1. **Stiffness.** The eigenvalues of the convective system are $u_n$ and $u_n\pm a$. As $M=u/a\to0$ the acoustic speeds $u_n\pm a$ dwarf the convective speed $u_n$; their ratio is $\sim 1/M$. The explicit step is set by the fast acoustic waves, so the slow convective field needs $\mathcal O(1/M)$ steps to converge.
2. **Excessive dissipation / accuracy loss.** Upwind Riemann solvers add dissipation proportional to the largest eigenvalue ($\sim a$). At low Mach this swamps the physical pressure–velocity coupling, and the discrete solution loses accuracy (the pressure field develops spurious $\mathcal O(M)$ checkerboard noise).

**Preconditioning** cures both by multiplying the time derivative with a matrix $\boldsymbol\Gamma$ that rescales the acoustic eigenvalues down toward the convective scale, **without changing the steady-state solution** (the preconditioner multiplies only the pseudo-time term, which vanishes at convergence):

$$
\boldsymbol\Gamma\,\frac{\partial \mathbf P}{\partial\tau}
+ \nabla\!\cdot\!\mathbf F^c(\mathbf P)
= \nabla\!\cdot\!\mathbf F^v(\mathbf P) + \mathbf S .
$$

ARES uses a **Weiss–Smith / Choi–Merkle-type** preconditioner formulated in the pressure-based primitive variables $(p,\mathbf u,h)$ — the natural set for this solver, since the thermodynamic state already carries pressure explicitly.

---

## Activation

Preconditioning is enabled through the integration-variable switch and tuned by two optional keys in `[ARES-Numerics]`:

```ini
[ARES-Numerics]
integration-variables    = prec     ; enable preconditioned update
preconditioning-Uref     = -1.0     ; reference velocity ( < 0 ⇒ use local sound speed )
preconditioning-eps-min  = -1.0     ; minimum cut-off  ( < 0 ⇒ default 0.10 )
riemann-solver           = HLLC Prec ; pair with a preconditioned flux
```

| Key | Default | Meaning |
|-----|:-------:|---------|
| `integration-variables = prec` | `prim` | Switch the update (and time step) to the preconditioned formulation |
| `preconditioning-Uref` | `-1.0` | Global reference velocity $U_\text{ref}$. If negative, the **local sound speed** is used (fully local preconditioning) |
| `preconditioning-eps-min` | `-1.0` | Lower cut-off $\varepsilon_\min$ on the preconditioning parameter. If negative, the built-in default **0.10** is used |

!!! warning "Use a matching Riemann solver"
    Preconditioning rescales the wave speeds, so the interface dissipation must be rescaled consistently. Pair `integration-variables = prec` with **`HLLC Prec`** or **`PLLF`**; using a non-preconditioned solver re-introduces the $\sim a$ dissipation and defeats the purpose.

---

## The reference velocity $U_r$

The heart of the method is the **local reference velocity** $U_r$, which sets the scale to which the acoustic eigenvalues are clipped. ARES computes it per cell from the local velocity magnitude and sound speed (routine `comp_Ur`):

$$
U_\text{ref} =
\begin{cases}
\texttt{preconditioning-Uref} & \text{if } > 0\\
a & \text{otherwise (local sound speed)}
\end{cases}
\qquad
\varepsilon_\min =
\begin{cases}
\texttt{preconditioning-eps-min} & \text{if } \ge 0\\
0.10 & \text{otherwise}
\end{cases}
$$

A local Mach-squared parameter is formed and clipped from below by $\varepsilon_\min$ and from above by 1:

$$
\varepsilon = \min\!\Bigl(1,\ \max\bigl(\varepsilon_\min,\ M_\text{loc}^2\bigr)\Bigr),
\qquad
M_\text{loc} = \frac{\lVert\mathbf u\rVert}{U_\text{ref}} .
$$

The reference velocity is then the local velocity, bounded into the band $[\varepsilon\,U_\text{ref},\ U_\text{ref}]$:

$$
U_r =
\begin{cases}
\varepsilon\,U_\text{ref} & \lVert\mathbf u\rVert < \varepsilon\,U_\text{ref}\quad(\text{stagnation floor})\\[4pt]
\lVert\mathbf u\rVert & \varepsilon\,U_\text{ref} \le \lVert\mathbf u\rVert \le U_\text{ref}\\[4pt]
U_\text{ref} & \lVert\mathbf u\rVert > U_\text{ref}\quad(\text{transonic/supersonic cap})
\end{cases}
$$

The three branches encode the standard cut-offs:

- **Stagnation floor** — near stagnation points $\lVert\mathbf u\rVert\to0$, so $U_r$ is held at $\varepsilon\,U_\text{ref}$ to prevent the preconditioner from becoming singular (and to keep some dissipation for stability).
- **Convective band** — in the bulk low-Mach flow, $U_r$ tracks the local velocity, exactly matching the dissipation to the convective scale.
- **Transonic cap** — once $\lVert\mathbf u\rVert$ reaches $U_\text{ref}$ (the sound speed, by default), $U_r$ is capped so the scheme smoothly reverts to the standard, un-preconditioned compressible solver. **Preconditioning therefore turns itself off automatically where it is not needed.**

---

## Effect on the eigenvalues and time step

With $U_r$ in hand, the preconditioned system has modified acoustic eigenvalues. The convective eigenvalue $u_n$ is unchanged, while the acoustic pair becomes

$$
\lambda^\pm = \tfrac12\,u_n\,(1+\alpha)\ \pm\ \sqrt{\tfrac14\,u_n^2(1-\alpha)^2 + U_r^2},
\qquad \alpha = 1 - \frac{U_r^2}{a^2},
$$

so a **preconditioned signal speed** $a' = \sqrt{\alpha^2 u_n^2 + U_r^2}$ replaces the physical sound speed $a$ in:

- the **CFL time step** — the much smaller $a'$ (instead of $a$) raises the stable $\Delta t$ at low Mach, removing the $\mathcal O(1/M)$ slowdown;
- the **Riemann dissipation** — `PLLF` and `HLLC Prec` use $a'$ / the preconditioned wave speeds, so the interface dissipation scales with the convective field rather than the acoustic field, restoring accuracy.

As $U_r\to a$ (high Mach) one has $\alpha\to0$, $a'\to a$, and the standard compressible scheme is recovered exactly.

---

## Practical guidance

- For low-Mach internal flows — the bundled [HTD](../vv/htd.md) and [Prt-correction](../vv/prt-correction.md) cases both run preconditioned — enabling `prec` with `HLLC Prec` both speeds convergence and sharpens the pressure field.
- Leave `preconditioning-Uref` negative (local sound speed) for general use. Set a positive global $U_\text{ref}$ when a single throughput velocity scale is meaningful for the whole domain and you want a uniform preconditioning level — this is what the bundled internal-flow cases do (e.g. `preconditioning-Uref = 600` for HTD, `120` for Prt-correction).
- Increase `preconditioning-eps-min` above 0.10 if the solver is noisy near stagnation regions (more dissipation there); decrease it for maximum low-Mach accuracy at some robustness cost.

---

## References

1. J. M. Weiss, W. A. Smith, "Preconditioning applied to variable and constant density flows," *AIAA J.* 33 (1995).
2. Y.-H. Choi, C. L. Merkle, "The application of preconditioning in viscous flows," *J. Comput. Phys.* 105 (1993).
3. E. Turkel, "Preconditioned methods for solving the incompressible and low speed compressible equations," *J. Comput. Phys.* 72 (1987).
