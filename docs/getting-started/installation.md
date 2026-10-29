# Installation

This document describes how to obtain and build **ARES**. The instructions describe the current `install.sh` script, the CMake configuration, and the Git submodule layout.

!!! note
    ARES has a dual nature: it is both a library and an executable. The installation process produces both the static library `libARES.a` and the main executable `bin/ARES`. If you only want to use ARES as a library, you can link against `libARES.a` without caring about the executable.

## Prerequisites

Before building ARES make sure your system provides the following tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`).
- **C / C++ compiler** – required only for optional components (TecIO, Cantera).
- **OpenMP / MPI** – needed for optional parallelization support.

**ifx 2024.1 and IPO.** ifx 2024.1 (tested: 2024.1.0, build 20240308) compiles a loop of the ghost metrics wrongly at `-O2` and `-O3` without interprocedural optimization: a minimal program in standard Fortran gives wrong values at `-O2` and `-O3`, also at the default level with no `-O` option, and the right ones at `-O0`, `-O1` and with `-ipo`. `BC_Extrapolate_Metrics` of ARES has that loop: ARES built so can stop at set-up on a negative volume, or run with wrong ghost metrics. CMake therefore refuses ifx 2024.1 at `-O2` or above without IPO, with a message: a RELEASE build with `-DUSE_IPO=OFF` (`--no-ipo`), or whose toolchain cannot do IPO, and the TESTING build. The RELEASE build with IPO (the default) and the DEBUG build (`-O0`) are accepted. Other versions of ifx have not been tested. The check runs only when ARES is built on its own: a project that embeds ARES (Hydra) sets its own flags.

### Git submodules

ARES depends on several repositories included as Git submodules.

| Path | Repository | Purpose |
|------|------------|---------|
| `lib/FLINT` | `github.com/MarcoGrossi92/FLINT` | Real-fluid thermodynamic & transport tables |
| `lib/ORION` | `github.com/MarcoGrossi92/ORION` | I/O routines (Tecplot, VTK, Plot3D, …) |
| `lib/third_party/FiNeR` | `github.com/szaghi/FiNeR` | INI file parser |

## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/open-hydra/ARES.git
cd ARES
# initialise submodules
./install.sh update            # or: git submodule update --init --recursive
```

To install ARES you may either use the bundled install script or invoke CMake manually. The script is the preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`. It also maintains a `CMakePresets.json` file that records the configuration used for the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global options**

* `-h`, `--help` – show the help message and exit.
* `-v`, `--verbose` – enable verbose logging.

**`build` command** — performs a clean configure + build cycle.

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compilers=gnu --use-openmp

# Intel compilers with MPI and TecIO
./install.sh build --compilers=intel --use-openmp --use-mpi --use-tecio
```

Options accepted by `build`:

* `--compilers=<gnu|intel>` – select the compiler family.
* `--use-openmp` – enable OpenMP parallelization.
* `--use-mpi` – enable MPI parallelization.
* `--use-tecio` – enable TecIO support (requires a C++ compiler).
* `--use-sundials` – build and link the bundled SUNDIALS library (requires a C compiler).
* `--use-cantera` – build and link Cantera (requires a C++ compiler).
* `--no-ipo` – disable interprocedural (link-time) optimization in the RELEASE build.
* `--include-orion=PATH` – use an external ORION tree instead of the submodule.
* `--include-flint=PATH` – same for FLINT.
* `--include-oslo=PATH` – same for OSLO.
* `--include-finer=PATH` – same for FiNeR.

After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

!!! note "Interprocedural optimization (IPO)"
    `USE_IPO` (on by default) adds `-ipo` to the RELEASE build of the Intel compilers. The flag is probed with a link step, since `-ipo` objects hold IR, not machine code; the link line gets `-ipo` too, and the static libraries are built with the `llvm-ar` / `llvm-ranlib` of oneAPI: in `<bin>/compiler` from oneAPI 2024 on, in `<bin>/../bin-llvm` up to oneAPI 2023, `<bin>` being the directory of `ifx` with its symbolic links resolved. `-DIPO_AR=... -DIPO_RANLIB=...` give the tools of any other layout. If the probe or the lookup fails, the build falls back to a plain RELEASE with a warning: `-ipo` leaves the compile and the link flags, and every other RELEASE flag stays. GNU's `-flto` is not offered. A build directory configured before keeps the flags it cached then: configure a new one (`./install.sh build` does).

**`compile` command** — re-runs CMake using the previously generated preset and rebuilds without clearing the build directory. Useful during development when only the source has changed.

```bash
./install.sh compile
```

**`update` command** — synchronises the Git submodules. By default it checks out the commit recorded in `.gitmodules`; `--remote` fetches the latest commit from each remote branch.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine-grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DUSE_MPI=OFF \
    -DORION_PATH=/path/to/ORION \        # optional
    -DFLINT_PATH=/path/to/FLINT \        # optional
    -DFINER_PATH=/path/to/FiNeR          # optional
cmake --build . --parallel
```

The resulting artifacts are placed in `build/` by default. The static library is `lib/libARES.a` and the main executable is `bin/ARES`.

## CMake presets

The `CMakePresets.json` produced by the install script records the cache variables used during configuration. You can then build later with:

```bash
cmake --preset default
cmake --build build
```

or with the `compile` command of the install script.

## Optional components

### MPI

`-DUSE_MPI=ON` (`--use-mpi`) builds ARES with MPI through the Fortran wrapper of the MPI library that CMake finds (`mpif90`, `mpiifx`). With GNU Fortran 10 or later the build compiles `Mod_MPI.f90` with `-fallow-argument-mismatch`: it calls the same MPI routine with buffers of different types and ranks.

**gfortran with Intel MPI.** The `mpif90` of Intel MPI 2021.12 puts on the include path the `mpi` module built for ifx, which gfortran cannot read ("Reading module .../include/mpi/mpi.mod ... Unexpected EOF"). Give CMake the modules built for gfortran:

```bash
cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DUSE_MPI=ON \
    "-DMPI_Fortran_COMPILER_INCLUDE_DIRS=$I_MPI_ROOT/include/mpi/gfortran/11.1.0;$I_MPI_ROOT/include"
```

The recipe was tested with gfortran 13.2.1: where `gfortran` is an older version, name the newer compiler, for example
`-DCMAKE_Fortran_COMPILER=gfortran-13`.

### TecIO

ARES can be built with support for TecIO, the library for writing Tecplot binary files (shipped with ORION). If enabled, ARES can read and write output in Tecplot binary (`.szplt`) format, which is useful for large datasets. Enabling TecIO requires a working C++ compiler.

## Library linking (advanced)

To use ARES from an external Fortran program:

```bash
gfortran -I/path/to/ARES/include \
         -L/path/to/ARES/lib \
         -lARES \
         your_program.f90 -o your_program
```

or, from a CMake project, add ARES as a sub-directory — the top-level `CMakeLists.txt` is written so that ARES can be embedded into a larger build (this is how the **Hydra** coupled solver consumes it).

## Next steps

* **[Quick Start Tutorial](quick-start.md)** – build, run, and verify your first case.

---
