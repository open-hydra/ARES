# Relative to this directory, NOT ${CMAKE_SOURCE_DIR}/src/lib: when ARES is
# built as a subproject (hydra) CMAKE_SOURCE_DIR is the superproject root, so
# the absolute form swept in every sibling under its src/lib -- MOSE, FUSS,
# GROOT, ARES's own src/app programs and FUSS's vendored TecIO examples.
file(GLOB_RECURSE ARES_SRC
    CONFIGURE_DEPENDS
    *.f90
)

add_library(ARESL ${ARES_SRC} ${FLINT_SRC})
add_library(ARESL::ARESL ALIAS ARESL)
target_link_libraries(ARESL PRIVATE FiNeR::FiNeR ORION)
link_FLINT_dependencies(ARESL)

target_include_directories(ARESL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    # Use the imported target so compile *and* link flags propagate
    # transitively to every consumer of ARESL (PUBLIC).
    # The old SET_TARGET_PROPERTIES(LINK_FLAGS) approach is deprecated and
    # its LINK_FLAGS are silently ignored for static libraries by CMake.
    target_link_libraries(ARESL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(ARESL PUBLIC USE_MPI)
    target_link_libraries(ARESL PUBLIC MPI::MPI_Fortran)
    # An mpi module paired with a compiler that cannot ignore the types of choice buffers declares no
    # interfaces for the routines with choice buffers (MPI 4.1, Section 19.1.3), as the gfortran modules of
    # Intel MPI do: gfortran >= 10 then stops on calls of one such routine with buffers of different types
    # or ranks. -fallow-argument-mismatch turns those errors into warnings; it is set in MPI builds only,
    # and on Mod_MPI.f90 only, the one file that makes such calls: the builds without MPI, and the other
    # files of ARES and FLINT in the MPI builds, keep the full check.
    IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        set_source_files_properties(parallel/Mod_MPI.f90 PROPERTIES COMPILE_OPTIONS -fallow-argument-mismatch)
    ENDIF()
ENDIF(USE_MPI)
