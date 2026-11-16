# Helpers for the libraries ARES builds from: initialize a git submodule (lib/FiNeR), and
# check that a library is at the commit or the release this ARES commit was tested with.

find_package(Git QUIET)

# A library that differs from what ARES requires (a submodule at another commit, an external clone of
# another release) is fatal; to build against a work-in-progress checkout on purpose (an ORION or FiNeR branch
# under test) configure with -DARES_PIN_CHECK=WARNING.
set(ARES_PIN_CHECK "FATAL_ERROR" CACHE STRING "Severity of a library that differs from the commit or release ARES requires: FATAL_ERROR (default) or WARNING")
if (NOT ARES_PIN_CHECK STREQUAL "FATAL_ERROR" AND NOT ARES_PIN_CHECK STREQUAL "WARNING")
  message(FATAL_ERROR "ARES_PIN_CHECK must be FATAL_ERROR or WARNING (got '${ARES_PIN_CHECK}')")
endif()

# -------- Initialize a git submodule from the correct git root --------
macro(ensure_submodule_initialized SUBMOD_PATH)
  if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
    message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
    get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
    execute_process(
      COMMAND git rev-parse --show-toplevel
      WORKING_DIRECTORY "${_sm_parent}"
      OUTPUT_VARIABLE _sm_git_root
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _sm_git_result)
    if(_sm_git_result EQUAL 0)
      file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
      execute_process(
        COMMAND git submodule update --init "${_sm_rel}"
        WORKING_DIRECTORY "${_sm_git_root}"
        RESULT_VARIABLE _sm_init_result)
      if(NOT _sm_init_result EQUAL 0)
        message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
      endif()
    else()
      message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
    endif()
  endif()
endmacro()

# -------- Check that a bundled submodule is at the commit this ARES commit records --------
# A mismatch typically follows `git checkout <branch>` without `git submodule update`: ARES would
# compile against a library version it was not tested with. Read-only git queries; nothing is checked for
# tarballs, plain copies or paths outside this tree (ares_check_version covers those for FiNeR).
# The hint syncs before updating: a clone made before a submodule changed URL or path (FiNeR, from
# lib/third_party/FiNeR to lib/FiNeR) may still point at the old one.
function(ares_check_submodule_pin SUBMOD_PATH)
  if (NOT GIT_FOUND OR "${SUBMOD_PATH}" STREQUAL "" OR NOT EXISTS "${SUBMOD_PATH}/.git")
    return()
  endif()
  file(REAL_PATH "${SUBMOD_PATH}" _abs)
  file(REAL_PATH "${CMAKE_SOURCE_DIR}" _top)
  file(RELATIVE_PATH _rel "${_top}" "${_abs}")
  if ("${_rel}" STREQUAL "" OR "${_rel}" MATCHES "^\\.\\.")
    return()
  endif()
  set(_git ${CMAKE_COMMAND} -E env GIT_OPTIONAL_LOCKS=0 ${GIT_EXECUTABLE})
  execute_process(COMMAND ${_git} ls-files -s -- "${_rel}" WORKING_DIRECTORY "${_top}"
                  OUTPUT_VARIABLE _entry RESULT_VARIABLE _r1 ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${_git} rev-parse --show-toplevel HEAD WORKING_DIRECTORY "${_abs}"
                  OUTPUT_VARIABLE _head RESULT_VARIABLE _r2 ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT _r1 EQUAL 0 OR NOT _r2 EQUAL 0 OR NOT "${_entry}" MATCHES "^160000 ([0-9a-f]+) ")
    return()   # not a registered submodule (gitlink) of this checkout
  endif()
  set(_pin "${CMAKE_MATCH_1}")
  string(REPLACE "\n" ";" _head "${_head}")
  list(LENGTH _head _n)
  if (NOT _n EQUAL 2)
    return()
  endif()
  list(GET _head 0 _root)
  list(GET _head 1 _sha)
  file(REAL_PATH "${_root}" _root)
  if (NOT "${_root}" STREQUAL "${_abs}")
    return()   # the directory has no repository of its own
  endif()
  if (NOT "${_sha}" STREQUAL "${_pin}")
    # Name each commit by its nearest tag (v2.3.0, v2.3.0-6-g9ecd696); a pin not fetched yet stays a sha
    string(SUBSTRING "${_pin}" 0 7 _pin_name)
    string(SUBSTRING "${_sha}" 0 7 _sha_name)
    execute_process(COMMAND ${_git} describe --tags "${_pin}" WORKING_DIRECTORY "${_abs}"
                    OUTPUT_VARIABLE _pin_tag RESULT_VARIABLE _r3 ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${_git} describe --tags "${_sha}" WORKING_DIRECTORY "${_abs}"
                    OUTPUT_VARIABLE _sha_tag RESULT_VARIABLE _r4 ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (_r3 EQUAL 0)
      set(_pin_name "${_pin_tag}")
    endif()
    if (_r4 EQUAL 0)
      set(_sha_name "${_sha_tag}")
    endif()
    message(${ARES_PIN_CHECK} "${_rel} is checked out at ${_sha_name}, but this ARES commit records ${_pin_name}. "
                    "ARES would be built against a library version it was not tested with. "
                    "Run `git submodule sync --recursive && git submodule update --init --recursive` and configure again "
                    "(or configure with -DARES_PIN_CHECK=WARNING to build against this checkout on purpose).")
  endif()
endfunction()

# -------- Check the release of a library against the one ARES pins --------
# A library from a path outside this tree (FETCHCONTENT_SOURCE_DIR_FINER) is not this tree's
# submodule, so ares_check_submodule_pin skips it. Called for the bundled library too, so a gitlink or tag
# bumped without its ARES_<LIB>_TAG is caught. HINT says how to get the pinned release.
function(ares_check_version NAME SOURCE_DIR FOUND PINNED HINT)
  if (NOT "${FOUND}" VERSION_EQUAL "${PINNED}")
    message(${ARES_PIN_CHECK} "${NAME} at ${SOURCE_DIR} is version '${FOUND}', but ARES pins ${NAME} ${PINNED}. "
                    "ARES would be built against a library version it was not tested with. ${HINT} "
                    "(or configure with -DARES_PIN_CHECK=WARNING to build against this ${NAME} on purpose).")
  endif()
endfunction()
