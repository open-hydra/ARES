function(add_fortran_executable target source)
  
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE ARESL FiNeR::FiNeR)

  if (USE_TECIO)
    # tecio/teciompi is a C++ library, so CMake would otherwise pick the C++
    # driver to link this executable. The program's entry point is a Fortran
    # main, whose runtime `main` is provided by the Fortran driver, so force the
    # Fortran linker and pull in the C++ runtime (libstdc++) for tecio's symbols.
    set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${target} PRIVATE -pthread -lstdc++)
  endif()

  target_include_directories(${target} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if (USE_OPENMP)
    target_compile_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()

  if (USE_MPI)
    target_compile_definitions(${target} PRIVATE USE_MPI)
    target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
  endif()

endfunction()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_fortran_executable(${solver} "main.f90")
add_fortran_executable(DocGen "docgen.f90")
