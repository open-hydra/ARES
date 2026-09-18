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
ENDIF(USE_MPI)
