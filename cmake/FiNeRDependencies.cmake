# FiNeR's dependencies. FiNeR tracks none of the libraries its CMakeLists adds from src/third_party/: its
# author fetches them with `fobis fetch`, which clones each [dependencies] entry of FiNeR's fobos that is
# absent (a plain URL: its default branch) and leaves an existing clone as it is.
#
#   ensure_finer_dependencies(<FiNeR folder> <manage>)
#
# With <manage> TRUE (the submodule lib/FiNeR) it does the same without FoBiS, as
# ARES_FINER_DEPS (set in CMakeLists.txt) chooses:
#   FREE    as `fobis fetch`, reading the list from FiNeR's fobos; a clone left detached (by BACKUP) goes back
#           to its branch, wherever that stands (nothing is fetched). Delete src/third_party/<dep> to take its
#           latest commit.
#   BACKUP  each dependency at its release in ARES_FINER_BACKUP_PINS, as `fobis fetch --frozen` with a lock;
#           a release the clone lacks is fetched from the URL in the fobos.
# A directory that holds files but no repository of its own is an error.
# With <manage> FALSE (an external clone, FETCHCONTENT_SOURCE_DIR_FINER) nothing is cloned or moved, as for the
# clone itself: every dependency must be there, and with BACKUP at its release (ARES_PIN_CHECK otherwise).

if (NOT ARES_FINER_DEPS STREQUAL "FREE" AND NOT ARES_FINER_DEPS STREQUAL "BACKUP")
  message(FATAL_ERROR "ARES_FINER_DEPS must be FREE or BACKUP (got '${ARES_FINER_DEPS}')")
endif()

function(ensure_finer_dependencies FINER_ROOT MANAGE)
  file(REAL_PATH "${FINER_ROOT}" _root)

  # [dependencies] of FiNeR's fobos: deps_dir and one `NAME = URL` line per dependency
  set(_deps_dir "src/third_party")
  set(_deps "")
  set(_in_deps FALSE)
  file(STRINGS "${_root}/fobos" _lines)
  foreach(_line IN LISTS _lines)
    if ("${_line}" MATCHES "^\\[([A-Za-z0-9_-]+)\\]")
      set(_in_deps FALSE)
      if ("${CMAKE_MATCH_1}" STREQUAL "dependencies")
        set(_in_deps TRUE)
      endif()
    elseif (_in_deps AND "${_line}" MATCHES "^([A-Za-z0-9_.+-]+)[ \t]*=[ \t]*(.*[^ \t])[ \t]*$")
      set(_name "${CMAKE_MATCH_1}")
      set(_spec "${CMAKE_MATCH_2}")
      if (_name STREQUAL "deps_dir")
        set(_deps_dir "${_spec}")
      elseif ("${_spec}" MATCHES "::")   # tag=, branch=, rev=, semver=, ... : a FiNeR newer than 2.3.0
        message(FATAL_ERROR "FiNeR's fobos fetches ${_name} as '${_spec}': ARES handles a plain URL only "
                            "(its default branch). Extend ensure_finer_dependencies for the options.")
      else()
        list(APPEND _deps "${_name}")
        set(_url_${_name} "${_spec}")
      endif()
    endif()
  endforeach()
  if (NOT _deps)
    message(FATAL_ERROR "No dependency found in the [dependencies] of ${_root}/fobos")
  endif()

  set(_report "")
  foreach(_dep IN LISTS _deps)
    set(_dir "${_root}/${_deps_dir}/${_dep}")
    set(_url "${_url_${_dep}}")
    set(_release "")
    if (ARES_FINER_DEPS STREQUAL "BACKUP")
      foreach(_entry IN LISTS ARES_FINER_BACKUP_PINS)
        if ("${_entry}" MATCHES "^${_dep}=(.+)$")
          set(_release "${CMAKE_MATCH_1}")
        endif()
      endforeach()
      if (_release STREQUAL "")
        message(FATAL_ERROR "FiNeR dependency ${_dep} has no backup release: add one to ARES_FINER_BACKUP_PINS")
      endif()
    endif()

    set(_files "")
    if (EXISTS "${_dir}")
      file(GLOB _files "${_dir}/*" "${_dir}/.*")
    endif()
    if (NOT _files AND NOT MANAGE)
      message(FATAL_ERROR "FiNeR at ${_root} has no dependency ${_dep} in ${_deps_dir}: fetch its dependencies "
                          "there (fobis fetch --no-build), or build the FiNeR that ARES checks out "
                          "(without FETCHCONTENT_SOURCE_DIR_FINER).")
    elseif (NOT _files)
      message(STATUS "FiNeR dependency ${_dep} missing. Cloning ${_url}...")
      execute_process(COMMAND git clone "${_url}" "${_dir}" RESULT_VARIABLE _result)
      if (NOT _result EQUAL 0)
        message(FATAL_ERROR "Failed to clone FiNeR dependency ${_dep} into ${_dir}")
      endif()
    endif()
    # A repository of its own: git commands in a plain copy would answer for the repository around it (FiNeR's)
    set(_own_repo FALSE)
    execute_process(COMMAND git rev-parse --show-toplevel WORKING_DIRECTORY "${_dir}"
                    OUTPUT_VARIABLE _git_root RESULT_VARIABLE _result
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if (_result EQUAL 0)
      file(REAL_PATH "${_git_root}" _git_root)
      if ("${_git_root}" STREQUAL "${_dir}")
        set(_own_repo TRUE)
      endif()
    endif()
    if (MANAGE AND NOT _own_repo)
      message(FATAL_ERROR "FiNeR dependency ${_dir} exists but is not a git checkout: "
                          "remove it to let ARES clone ${_dep} from ${_url}.")
    endif()

    if (NOT _release STREQUAL "")
      # BACKUP: detached at the commit of the release tag (a clone of ARES fetches a tag it lacks)
      set(_tag "refs/tags/${_release}")
      set(_commit "")
      set(_head "")
      if (_own_repo)
        execute_process(COMMAND git cat-file -e "${_tag}" WORKING_DIRECTORY "${_dir}" RESULT_VARIABLE _result ERROR_QUIET)
        if (NOT _result EQUAL 0 AND MANAGE)
          execute_process(COMMAND git fetch "${_url}" "${_tag}:${_tag}" WORKING_DIRECTORY "${_dir}"
                          RESULT_VARIABLE _result)
          if (NOT _result EQUAL 0)
            message(FATAL_ERROR "Failed to fetch release ${_release} of FiNeR dependency ${_dep} from ${_url}")
          endif()
        endif()
        execute_process(COMMAND git rev-parse "${_tag}^{commit}" WORKING_DIRECTORY "${_dir}"
                        OUTPUT_VARIABLE _commit OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${_dir}"
                        OUTPUT_VARIABLE _head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      endif()
      if ("${_commit}" STREQUAL "" OR NOT "${_head}" STREQUAL "${_commit}")
        if (NOT MANAGE)
          message(${ARES_PIN_CHECK} "FiNeR dependency ${_dep} in ${_dir} is not at ${_release}, its backup release "
                  "(ARES_FINER_DEPS=BACKUP). Check the release out there (git -C ${_dir} checkout ${_release}), or "
                  "configure with -DARES_PIN_CHECK=WARNING to build it as it is.")
        else()
          message(STATUS "FiNeR dependency ${_dep} at '${_head}', moving it to its backup release ${_release}...")
          execute_process(COMMAND git -c advice.detachedHead=false checkout --detach "${_tag}"
                          WORKING_DIRECTORY "${_dir}" RESULT_VARIABLE _result)
          execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${_dir}"
                          OUTPUT_VARIABLE _head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
          if (NOT _result EQUAL 0 OR "${_commit}" STREQUAL "" OR NOT "${_head}" STREQUAL "${_commit}")
            message(FATAL_ERROR "Failed to check out FiNeR dependency ${_dep} at ${_release} in ${_dir}")
          endif()
        endif()
      endif()
    elseif (MANAGE)
      # FREE: on the branch `git clone` gave it, wherever that branch stands (no fetch, as fobis without --update)
      execute_process(COMMAND git symbolic-ref -q HEAD WORKING_DIRECTORY "${_dir}"
                      RESULT_VARIABLE _attached OUTPUT_QUIET ERROR_QUIET)
      if (NOT _attached EQUAL 0)
        execute_process(COMMAND git symbolic-ref --short refs/remotes/origin/HEAD WORKING_DIRECTORY "${_dir}"
                        OUTPUT_VARIABLE _branch RESULT_VARIABLE _result OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        string(REGEX REPLACE "^origin/" "" _branch "${_branch}")
        if (NOT _result EQUAL 0 OR _branch STREQUAL "")
          message(FATAL_ERROR "FiNeR dependency ${_dep} is detached and its default branch is unknown: "
                              "check out a branch in ${_dir}, or remove it to let ARES clone it again.")
        endif()
        message(STATUS "FiNeR dependency ${_dep} detached, moving it back to its branch ${_branch}...")
        execute_process(COMMAND git checkout -q "${_branch}" WORKING_DIRECTORY "${_dir}" RESULT_VARIABLE _result)
        if (NOT _result EQUAL 0)
          message(FATAL_ERROR "Failed to check out branch ${_branch} of FiNeR dependency ${_dep} in ${_dir}")
        endif()
      endif()
    endif()

    set(_version "(not a git checkout)")
    if (_own_repo)
      execute_process(COMMAND git describe --tags --always HEAD WORKING_DIRECTORY "${_dir}"
                      OUTPUT_VARIABLE _version OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    list(APPEND _report "${_dep} ${_version}")
  endforeach()
  list(JOIN _report ", " _report)
  if (MANAGE)
    message(STATUS "FiNeR dependencies (ARES_FINER_DEPS=${ARES_FINER_DEPS}): ${_report}")
  else()
    message(STATUS "FiNeR dependencies of ${_root}, as they are (ARES_FINER_DEPS=${ARES_FINER_DEPS}): ${_report}")
  endif()
endfunction()
