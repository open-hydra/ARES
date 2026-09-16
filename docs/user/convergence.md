# Convergence & Numerical Tuning

[Input Parameters](registry.md) lists *what* every `[ARES-Numerics]` and `[ARES-IO]` key is. This page is about *when to reach for them* — the practical levers that get a stubborn case to converge, and the ones that let you watch and steer a run while it is going.

---

## Steady state or time-accurate?

The single switch that changes the meaning of everything else is `time-accurate`.

| | `time-accurate = .false.` (default) | `time-accurate = .true.` |
|---|---|---|
| Time step | **local** — every cell marches at its own stability limit | **global** — every cell is overwritten with the domain-wide minimum |
| Physical time | set to the sentinel $-1$ and never advanced | accumulated as $t \mathrel{+}= \Delta t_\text{global}$ |
| Use for | steady RANS, fastest route to a converged field | unsteady/URANS, anything where the transient is the answer |

Local time stepping is a convergence accelerator, not a physical integration: cells with small volumes take small steps and large cells take large ones, so the field reaches steady state in far fewer iterations. The price is that intermediate states are meaningless. In steady mode the `time` column printed to the shell is the $-1$ sentinel — **watch the residual, not the time**.

---

## The two stability numbers

`cfl` (required) and `vnn` set the convective and viscous limits; the local step is the minimum over directions of both:

$$
\Delta t_\text{conv} = \mathrm{CFL}\,\frac{\Delta \ell_d}{|u_d| + a},
\qquad
\Delta t_\text{visc} = \mathrm{VNN}\,\frac{\rho\,\Delta \ell_d^2}{\mu_\ell + \mu_t}
$$

With `integration-variables = prec` the convective wave speed in the first expression is itself replaced by its preconditioned counterpart, so the step no longer collapses as the Mach number falls.

On a wall-resolved viscous mesh the near-wall cells are thin, $\Delta\ell^2$ collapses and **the VNN limit dominates** — raising `cfl` alone then buys nothing. If the step seems pinned regardless of `cfl`, `vnn` is what is binding.

---

## Getting a run off the ground: `cfl-rise-threshold`

The most common failure is not a bad scheme, it is the first fifty iterations. A uniform initial field against a strong boundary condition — an imposed heat flux, a large pressure ratio, an impulsive inflow — produces a transient far more violent than the solution it converges to.

`cfl-rise-threshold = N` ramps the step linearly over the first `N` iterations:

```fortran
if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter
```

At iteration 1 the step is $1/N$ of nominal, at iteration `N` it is full. `0` (the default) disables it.

Note that it scales $\Delta t$ itself, so it ramps the convective **and** viscous limits together. Typical values are a few hundred to a few thousand iterations — enough to cover the initial transient, short enough not to waste the run. If a case diverges in the first handful of iterations and is stable afterwards, this is the first key to try, before lowering `cfl` permanently.

---

## Buying a larger CFL: implicit residual smoothing

`irs = .true.` with `irs-beta = β` smooths the residual before the update,

$$
R^\ast = R + \beta\,\nabla^2 R^\ast
$$

solved with Jacobi sweeps applied one direction at a time — in each direction

$$
R^\ast_i = \frac{R_i + \beta\,(R^\ast_{i-1} + R^\ast_{i+1})}{1 + 2\beta}
$$

Smoothing the residual widens the stability envelope of the explicit scheme, so IRS is how you run at a CFL the bare Runge–Kutta would not tolerate. `irs-beta = 0` makes it a no-op; the useful range is a few tenths. It costs a sweep per direction per stage, so it pays only if it lets you raise `cfl` by more than that overhead. See [Time Integration](../theory/time-integration.md) for the derivation.

---

## Low-Mach cases: preconditioning

When the Mach number is low the acoustic and convective scales separate and a density-based scheme stalls. `integration-variables = prec` switches the update to preconditioned variables; `preconditioning-Uref` sets the reference velocity (negative ⇒ local sound speed) and `preconditioning-eps-min` the lower cut-off (negative ⇒ 0.10).

Pair it with a preconditioned flux — `riemann-solver = HLLC Prec` or `PLLF` — otherwise the dissipation of the Riemann solver still scales on the acoustic speed and you lose most of the benefit. Full treatment in [Low-Mach Preconditioning](../theory/preconditioning.md).

---

## Multigrid

`[ARES-Multigrid]` takes one `levelN-iter` per grid level, each becoming that level's iteration budget. Iterating on coarse levels first removes the low-frequency error cheaply; see [Time Integration](../theory/time-integration.md).

---

## Trading accuracy for robustness in space

Two keys, both in `[ARES-Numerics]`, ordered here from most robust to most accurate:

- `space-reconstruction` — `first-order`/`none` is unconditionally the most forgiving; `MUSCL` is second order and what you want for the final answer.
- `flux-limiter` — only read when `space-reconstruction = MUSCL`: specify one without the other and it is ignored with a warning, while MUSCL without a limiter silently falls back to Van Leer (also warned). `minmod` is the most dissipative and most robust; `vanalbada`, `vanleer` and `MC` sit in the middle; `superbee` is the sharpest and the most likely to ring. `LIMO3` is third-order.
- `riemann-solver` — `Rusanov` is the most dissipative, `HLLE` and `PLLF` intermediate, `HLLC` the most accurate. `HLLC Rotated` helps on grid-aligned shocks.

A reliable recipe for a case that will not start: run first-order until the field is established, then switch to `MUSCL` with `minmod`, and only then sharpen the limiter. Details in [Spatial Discretization](../theory/numerics.md) and [Riemann Solvers](../theory/riemann-solvers.md).

---

## Watching a run: the `[ARES-IO]` keys

Convergence you cannot see is convergence you cannot diagnose.

- `res-diter` — how often the residual history is appended. Keep it at `1` while tuning; this is the curve that tells you whether you are converging, stalled, or diverging.
- `shell-diter` — console refresh. Large values on long batch runs, small while watching.
- `sol-diter` / `sol-dtime` — field output by iteration or by physical time. Both default to effectively never, so a run left alone writes only at the end.
- `sol-overwrite = false` — write **numbered snapshots** instead of overwriting. This is the key to use when a run diverges and you need to see *where* it went wrong: without it the only surviving field is the corrupted one.
- `sol-variables`, `wall-variables` — which variable groups end up in the output; `ini-format` / `sol-format` — ASCII or binary, Tecplot or VTK. ASCII is portable and slow, binary the opposite.
- `[ARES-Probes]` — point histories by index or coordinates, with their own `diter`/`dtime`. The cheapest way to watch a single location converge without dumping fields.

---

## Steering a run without restarting it: `ini-diter`

Every `ini-diter` iterations the solver re-reads `input.ini` and re-validates the whole registry. Any key that the solver consults **every iteration** therefore takes effect live, with no restart:

- `cfl`, `vnn` — start conservative, watch the residual, raise them once the transient has passed
- `sol-diter`, `res-diter`, `shell-diter`, `sol-overwrite` — start dumping snapshots the moment a run starts misbehaving

!!! warning "What hot-reload does *not* change"
    Keys consumed **once at setup** are copied into the running state and will not change mid-run even though the registry value updates. The important one is `iter-threshold`: it is read into `domain%itermax` during setup, so editing it while the solver runs will not extend or shorten the run. The same applies to anything that binds the solver's structure — turbulence model, mesh, multigrid level count, output format.

Set `ini-diter` small enough to be responsive and large enough not to hammer the filesystem; the default is 10000. Bear in mind that the file is re-validated on every reload, so a syntax error introduced while editing will be caught mid-run.

---

## Quick triage

| Symptom | First thing to try |
|---|---|
| Diverges within the first iterations, from a uniform field | `cfl-rise-threshold` |
| Step pinned however high `cfl` goes | the VNN limit is binding — `vnn`, or the near-wall spacing |
| Stalls at low Mach, residual flat | `integration-variables = prec` + `HLLC Prec` |
| Residual oscillates without dropping | more dissipative `flux-limiter`, or `first-order` to re-establish |
| Stable but too slow | `irs` + higher `cfl`, or multigrid levels |
| Diverged and no usable output | `sol-overwrite = false`, rerun |
