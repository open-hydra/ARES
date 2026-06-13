# Getting Started

Welcome to the ARES getting-started guide! This section will help you install ARES and run your first example.

!!! info "ARES within Hydra"
    ARES is the real-fluid compressible-flow solver of **Hydra**. The pre-processor **ATLAS** generates the data ARES requires — mesh, initial conditions, boundary conditions, and real-fluid thermodynamic tables. In the Hydra workflow these are driven by the same `input.ini`. See the [overview page](../overview.md#hydra-cfd-suite) for the ecosystem.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Build the library and the executable from source

    [:octicons-arrow-right-24: Install ARES](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Get up and running in minutes

    [:octicons-arrow-right-24: Quick start tutorial](quick-start.md)

</div>

## Scope of This Section

This getting started guide covers:

1. **Installation**
   - Step-by-step build instructions
   - Build options and customization

2. **Quick Start**
   - Run your first case — a turbulent flat plate validated against the NASA Turbulence Modeling Resource
