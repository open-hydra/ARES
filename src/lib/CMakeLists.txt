file(GLOB_RECURSE ARES_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/src/lib/*.f90
)

add_library(ARESL ${ARES_SRC} ${FLINT_SRC})
add_library(ARESL::ARESL ALIAS ARESL)
target_link_libraries(ARESL PRIVATE FiNeR::FiNeR ORION)
link_FLINT_dependencies(ARESL)

target_include_directories(ARESL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(ARESL PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(ARESL PUBLIC USE_MPI)
    target_link_libraries(ARESL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
