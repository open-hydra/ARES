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

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags 
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

  # When GCC is used with -flto (slim LTO), plain ar cannot index LTO IR symbols
  # in static archives, causing undefined-reference linker errors. Override the
  # archiver and ranlib with the LTO-aware gcc-ar / gcc-ranlib wrappers so that
  # every static library built by this project is correctly indexed.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    find_program(GCC_AR     NAMES gcc-ar-12 gcc-ar-11 gcc-ar-10 gcc-ar
                 HINTS /usr/bin DOC "LTO-aware ar for GCC")
    find_program(GCC_RANLIB NAMES gcc-ranlib-12 gcc-ranlib-11 gcc-ranlib-10 gcc-ranlib
                 HINTS /usr/bin DOC "LTO-aware ranlib for GCC")
    if(GCC_AR)
      set(CMAKE_AR     "${GCC_AR}"     CACHE FILEPATH "LTO-aware ar (GCC)"     FORCE)
    endif()
    if(GCC_RANLIB)
      set(CMAKE_RANLIB "${GCC_RANLIB}" CACHE FILEPATH "LTO-aware ranlib (GCC)" FORCE)
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
