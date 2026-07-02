# Acknowledgements

ARES is developed and maintained by a team of contributors within the **Hydra** project. The author list is kept in [`AUTHORS.md`](https://github.com/open-hydra/ARES/blob/main/AUTHORS.md) at the repository root.

## Original Authors

- Andrea Giacomi

## Present Maintainers

- Andrea Giacomi
- Marco Grossi
- Alessio Sereno
- Vincenzo Barbato
- Alessandro Montanari
- Alex Falco
- Giacomo Passarani
- Marco Fabiani
- Gianluca Cocirla

---

## Built On

ARES stands on several open-source libraries:

| Library | Role | Author |
|---------|------|--------|
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Real-fluid thermodynamic & transport tables | M. Grossi |
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format structured-grid I/O | M. Grossi |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration parser | S. Zaghi |

- **FLINT** provides the tabulated $(p,h)$ real-fluid equation of state and transport properties that make ARES a real-fluid solver, with reference data from CoolProp / NASA correlations.
- **ORION** handles all structured-grid input and output across the Tecplot (ASCII / binary via TecIO) and VTK formats.
- **FiNeR** parses the `input.ini` configuration that drives both ARES and the registry-based parameter documentation.

It is also indebted to the wider CFD community whose methods it implements — the authors of the Riemann solvers, turbulence closures, and preconditioning techniques cited throughout the [Theory Guide](../theory/index.md).

---

## Documentation Tools

This documentation is built with:

- **[MkDocs](https://www.mkdocs.org/)** — static-site generator.
- **[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)** — theme and UI components.
- **MathJax** (via `pymdownx.arithmatex`) for the equations and **Mermaid** (via `pymdownx.superfences`) for the diagrams.

---

## License Compliance

ARES is distributed under the [GNU General Public License v3.0](license.md). The bundled dependencies carry licenses compatible with redistribution under GPL-3.0 (FLINT, ORION, FiNeR are open-source; TecIO is optional and only linked when explicitly enabled). When redistributing ARES or derivative works, the GPL-3.0 source-disclosure and share-alike terms apply.
