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

  # ORION
  if (NOT TARGET ORION)
    set ( orionD "${ORION_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
    add_subdirectory(${orionD})
  endif()

  # FLINT
  if (NOT TARGET FLINT)
    set ( flintD "${FLINT_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/FLINT/")
    add_subdirectory(${flintD})
  endif()

  # FiNeR
  if (NOT TARGET FiNeR)
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
