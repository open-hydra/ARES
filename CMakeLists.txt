# CMake project file for 'ARES'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver ARES)

project (${solver} LANGUAGES Fortran CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # Uncomment the below if you want the user to choose a parallelization library
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
  # Interprocedural (link-time) optimization of the RELEASE build. On by default,
  # switchable off for toolchains where it misbehaves. It also self-disables if
  # unsupported.
  OPTION(USE_IPO "Use interprocedural optimization in RELEASE builds" ON)

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

  # ifx 2024.1 compiles a loop of the ghost metrics of ARES wrongly at -O2 and -O3 without interprocedural optimization:
  # a minimal program in standard Fortran (an automatic array of a derived type, bounds from a module variable, filled
  # element by element and passed to a procedure) gives wrong values at -O2 and -O3, also at the default level with no -O
  # option, and the right ones at -O0, -O1 and with -ipo. BC_Extrapolate_Metrics has that loop: ARES built so can stop
  # at set-up on a negative volume, or run with wrong ghost metrics. Such a build is refused. The level is the last -O
  # option of the flags of the build type (-O2, the default of ifx, when there is none); IPO is -ipo in those flags (the
  # RELEASE build adds it with USE_IPO, the default, and drops it when the toolchain cannot do it) or
  # CMAKE_INTERPROCEDURAL_OPTIMIZATION. This only reads the flags: it changes none.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 2024.1
     AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2024.2)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _ifx_bt)
    separate_arguments(_ifx_flags UNIX_COMMAND "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_ifx_bt}}")
    set(_ifx_level "no -O option (-O2)")
    set(_ifx_high TRUE)
    set(_ifx_ipo FALSE)
    if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_ifx_bt})
      set(_ifx_ipo TRUE)
    endif()
    foreach(_ifx_f IN LISTS _ifx_flags)
      if(_ifx_f MATCHES "^[-/]O[01]$")
        set(_ifx_level "${_ifx_f}")
        set(_ifx_high FALSE)
      elseif(_ifx_f MATCHES "^[-/]O([23s]|fast)?$")
        set(_ifx_level "${_ifx_f}")
        set(_ifx_high TRUE)
      elseif(_ifx_f MATCHES "^[-/]Q?ipo$")
        set(_ifx_ipo TRUE)
      endif()
    endforeach()
    if(_ifx_high AND NOT _ifx_ipo)
      message(FATAL_ERROR "ifx ${CMAKE_Fortran_COMPILER_VERSION} with ${_ifx_level} and without IPO compiles a loop of the "
                          "ghost metrics of ARES wrongly (an error of the compiler: a minimal program in standard Fortran "
                          "gives wrong values at -O2 and -O3 and the right ones at -O0, -O1 and with -ipo); ARES built so "
                          "can stop at set-up on a negative volume, or run with wrong ghost metrics. Build with IPO "
                          "(RELEASE with -DUSE_IPO=ON, the default), with -O0 or -O1 (DEBUG), or with gfortran. See "
                          "docs/getting-started/installation.md, \"ifx 2024.1 and IPO\".")
    endif()
  endif()

  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # ORION, OSLO and FLINT by release range (cmake/SemVerDeps.cmake), checked out in lib/, releases kept in deps.lock: all
  # three required before any is made available, as FLINT requires ORION and OSLO too. FLINT sets FLINT_SRC, the
  # sources src/lib compiles, in this directory
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
  semver_dependency(OSLO GIT_REPOSITORY https://github.com/MarcoGrossi92/OSLO.git VERSION 1.1 GIT_SUBMODULES ""
                    SOURCE_DIR lib/OSLO)
  semver_dependency(ORION GIT_REPOSITORY https://github.com/MarcoGrossi92/ORION.git VERSION 2.1 SOURCE_DIR lib/ORION)
  semver_dependency(FLINT GIT_REPOSITORY https://github.com/MarcoGrossi92/FLINT.git VERSION 2.0 GIT_SUBMODULES ""
                    SOURCE_DIR lib/FLINT)
  FetchContent_MakeAvailable(OSLO ORION FLINT)

  # The libraries of FiNeR (PENF, FACE, FLAP, BeFoR64, StringiFor), submodules of ARES under FINER_DEPS_PATH, at their
  # release tags, added before FiNeR in the order of their dependencies: FiNeR adds each of them from its own
  # src/third_party/ only when its target is not defined yet, so a build needs neither FoBiS nor the network
  if (NOT DEFINED FINER_DEPS_PATH)
    set(FINER_DEPS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/")
  endif()
  foreach (finer_dep PENF FACE FLAP BeFoR64 StringiFor)
    if (NOT TARGET ${finer_dep})
      add_subdirectory("${FINER_DEPS_PATH}/${finer_dep}" "${CMAKE_BINARY_DIR}/lib/third_party/${finer_dep}")
    endif()
  endforeach()

  # FiNeR
  if (NOT TARGET FiNeR)
    # ARES needs FiNeR v2.1.0 or later: before it, a comment line that holds "=" or a bracket was taken for an option or
    # a section header, and the reading of its section of input.ini ended there. A clone that pulls without
    # "git submodule update" keeps the FiNeR it had checked out, and it would build in silence. FiNeR names its version
    # in a VERSION file ("vX.Y.Z") since v2.0.7; a FiNeR without it is older
    set(finer_version "")
    if (EXISTS "${FINER_PATH}/VERSION")
      file(STRINGS "${FINER_PATH}/VERSION" finer_version LIMIT_COUNT 1)
      string(REGEX REPLACE "^[ \t]*v?([0-9][0-9.]*).*$" "\\1" finer_version "${finer_version}")
    endif()
    if (NOT "${finer_version}" VERSION_GREATER_EQUAL "2.1.0")
      if ("${finer_version}" STREQUAL "")
        set(finer_found "a FiNeR without a VERSION file, older than v2.0.7")
      else()
        set(finer_found "FiNeR v${finer_version}")
      endif()
      message(FATAL_ERROR "FiNeR at ${FINER_PATH} is ${finer_found}: ARES needs FiNeR v2.1.0 or later (it pins v2.3.0), "
                          "as before v2.1.0 a comment line of input.ini that holds \"=\" ends the reading of its section. "
                          "Run: git submodule update --init --recursive")
    endif()
    set ( finerD "${FINER_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/third_party/FiNeR/")
    add_subdirectory(${finerD})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

OPTION(BASIC_TEST "Build basic tests" OFF)
if(BASIC_TEST)
  add_subdirectory(src/test)
endif()

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)
endif()
