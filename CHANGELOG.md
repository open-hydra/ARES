# Changelog

All notable changes to ARES are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.0.0] - 2026-11-17

### Changed (breaking)

- Input layout aligned with MOSE: `equations` + `turbulence` replace `simulation-type`; `[ARES-Turbulence]`, `schmidt`, `ic-format`.
  Old keys stop the run and say where they went (see `docs/user/input.md`).
- `time-scheme` is required.
- New names of the turbulence variables in input and output.

### Added

- Chimera boundary condition `102`, with per-level MPI and multigrid exchange.
- Boundary conditions `103` (multidomain coupling) and `104`, both with wall roughness.
- Wall `301` and wall output; probes placed by position, MPI-safe.
- SSG/LRR variants; SST asymptotic, `sst-production`, SST-RC; `omega-wall-bc` (`practical` or `asymptotic`).
- Point-implicit turbulence sources.
- `diffusive-metric`, `irs-variables`, `decomposition-map`; additive time step; reports of inflow through `400` faces.
- Mesh checks, wall distance, quasi-2D meshes.
- Restart counters; a warning when a newer numbered solution exists.
- Checks of `input.ini` and `bc.txt`: unknown keys, bad records, out-of-range cells and faces.
- Real-fluid states outside the tables are counted and reported.
- Wall-clock instrumentation, ported from MOSE.
- IPO in RELEASE builds (`USE_IPO`); ifx 2024.1 refused at `-O2` and above without IPO.
- A [Versioning](docs/development/versioning.md) page: the releases ARES requires of its libraries.
- This changelog.

### Changed

- ORION, OSLO and FLINT required by release range (ORION 2.1, OSLO 1.1, FLINT 2.1), no longer submodules.
- FiNeR moved to `lib/FiNeR`, pinned at v2.3.0; its libraries are fetched by the configure.
- A library at another release than ARES requires stops the configure (`-DARES_PIN_CHECK=WARNING` to override).
- MPI: exchange on the primary thread, short tags, output by the owner rank; stops are collective.
- `flux-limiter` accepts `mc` as well as `MC`.
- Documentation ported to Zensical; theory and user pages aligned with the code.
- `extract1d` of the tests is built from its source, no longer tracked.

### Fixed

- Prandtl correction of rough pipes; wall eddy viscosity in the boundary fluxes.
- SA, SA-rough, Wilcox and SA-comp fixes; RSM symmetry; SSG/LRR floor of omega.
- Multigrid: restriction, level checks, file names on coarse levels.
- HLLC Prec and time step with the true speed of sound; LIMO3, van Albada and QCR fixes.
- `CITATION.cff` is valid CFF 1.2.0 (GitHub's "Cite this repository" works).

## [0.1.0-alpha] - 2026-06-13

First tagged version.

[Unreleased]: https://github.com/open-hydra/ARES/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/open-hydra/ARES/compare/v0.1.0-alpha...v1.0.0
[0.1.0-alpha]: https://github.com/open-hydra/ARES/releases/tag/v0.1.0-alpha
