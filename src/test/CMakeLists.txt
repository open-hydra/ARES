file(GLOB_RECURSE TESTS *.f90 *.F90)
set(TESTS ${TESTS} PARENT_SCOPE)

foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  set(EXE_DIR ${CMAKE_SOURCE_DIR}/test/basic/${EXE_NAME})

  add_executable(${EXE_NAME} ${TEST})
  set_target_properties(${EXE_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${EXE_DIR}
      OUTPUT_NAME ${EXE_NAME}.exe)
  target_link_libraries(${EXE_NAME} PRIVATE ${lib} FiNeR::FiNeR OSLO ORION)
  target_include_directories(${EXE_NAME} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if(USE_OPENMP)
    target_compile_options(${EXE_NAME} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${EXE_NAME} PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()
endforeach()