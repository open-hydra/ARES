# SemVerDeps: the in-house dependencies of a project, required by version range and fetched with FetchContent, one
# copy of each per build tree
#
#   include(SemVerDeps)
#   semver_dependency(<name> GIT_REPOSITORY <url> VERSION <requirement>... [SOURCE_DIR <dir>]
#                     [<FetchContent_Declare option>...])
#   semver_dependency(<name> GIT_REPOSITORY <url> GIT_TAG <tag> [SOURCE_DIR <dir>] [<FetchContent_Declare option>...])
#   FetchContent_MakeAvailable(<name>...)
#
# A requirement is a caret range on the vX.Y.Z tags of <url>: 1.7 accepts 1.7.0 and later up to 2.0.0 excluded, 1
# accepts 1.0.0 and later, 0.3 accepts 0.3.0 up to 0.4.0 excluded. Several requirements are alternatives: VERSION 1.7 2
# accepts 1.x from 1.7.0 on and 2.x. GIT_TAG pins one tag instead (third-party libraries).
#
# SOURCE_DIR <dir> (relative to the calling directory) keeps the sources out of the build tree: SemVerDeps checks the
# release out there itself, since FetchContent would delete the folder and clone it again in each new build folder.
# The folder is cloned if missing; at the release, it is used as it is, local changes included; at another release,
# it is moved to the one selected if it is as SemVerDeps left it (a detached HEAD without changes), otherwise the
# configure stops and leaves it alone. GIT_SUBMODULES selects the submodules to initialize, as for FetchContent: those
# not populated are initialized at each configure, and those populated keep their commit and their changes, except
# when the checkout is moved. As for FetchContent, FETCHCONTENT_FULLY_DISCONNECTED forbids any clone or fetch, and
# FETCHCONTENT_UPDATES_DISCONNECTED(_<NAME>) any in an existing checkout: the configure stops when one is needed.
#
# The first call for <name> in the build tree selects the release; a later call (from a nested project) checks it
# against its own requirement and stops the configure when it is not satisfied. The release is, in order:
#   1. the checkout given with -DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>, as it is (its version is the newest vX.Y.Z tag
#      in its history; pre-releases are skipped)
#   2. the release of the lock file deps.lock, in the top-level source directory, if it satisfies the requirement
#   3. the newest release on <url> that satisfies the requirement (git ls-remote); pre-releases are skipped
# At the end of the configure deps.lock lists the releases selected: the next configure uses them without the
# network. -DSEMVER_UPDATE=ALL, or -DSEMVER_UPDATE=<name>;<name>, selects again from the remotes; it is then OFF again.
#
# Every project carries a copy of this file: the newest copy in a build tree is the one used.

set(_semver_deps_version 1.1.1)
get_property(_semver_deps_loaded GLOBAL PROPERTY SEMVER_DEPS_VERSION)
if(_semver_deps_loaded AND _semver_deps_loaded VERSION_GREATER_EQUAL _semver_deps_version)
  return()
endif()
if(_semver_deps_loaded)
  message(STATUS "SemVerDeps ${_semver_deps_version} of ${PROJECT_NAME} replaces ${_semver_deps_loaded}")
endif()
set_property(GLOBAL PROPERTY SEMVER_DEPS_VERSION ${_semver_deps_version})

if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "SemVerDeps needs CMake 3.19 or later")
endif()
include(FetchContent)
find_package(Git QUIET)

set(SEMVER_UPDATE OFF CACHE STRING
    "Dependencies to select again from their remotes at the next configure: ALL, or a list of names")
set(SEMVER_LOCK_FILE "${CMAKE_SOURCE_DIR}/deps.lock" CACHE FILEPATH "Releases of the dependencies selected")

# _semver_range(<requirement> <out_min> <out_max>): the versions a caret requirement accepts, from min up to max
# excluded
function(_semver_range req out_min out_max)
  if(NOT req MATCHES "^v?([0-9]+)(\\.([0-9]+))?(\\.([0-9]+))?$")
    message(FATAL_ERROR "SemVerDeps: '${req}' is not a version requirement (X, X.Y or X.Y.Z)")
  endif()
  set(_maj "${CMAKE_MATCH_1}")
  set(_min "${CMAKE_MATCH_3}")
  set(_pat "${CMAKE_MATCH_5}")
  # The leftmost non-zero number is the one that cannot change (a 0.y release may break 0.y-1)
  if(_maj GREATER 0 OR _min STREQUAL "")
    math(EXPR _next "${_maj} + 1")
    set(_max "${_next}.0.0")
  elseif(_min GREATER 0 OR _pat STREQUAL "")
    math(EXPR _next "${_min} + 1")
    set(_max "0.${_next}.0")
  else()
    math(EXPR _next "${_pat} + 1")
    set(_max "0.0.${_next}")
  endif()
  if(_min STREQUAL "")
    set(_min 0)
  endif()
  if(_pat STREQUAL "")
    set(_pat 0)
  endif()
  set(${out_min} "${_maj}.${_min}.${_pat}" PARENT_SCOPE)
  set(${out_max} "${_max}" PARENT_SCOPE)
endfunction()

# _semver_satisfies(<version> <out> <requirement>...): whether <version> satisfies one of the requirements
function(_semver_satisfies version out)
  set(_ok FALSE)
  foreach(_req IN LISTS ARGN)
    _semver_range("${_req}" _min _max)
    if(version VERSION_GREATER_EQUAL _min AND version VERSION_LESS _max)
      set(_ok TRUE)
    endif()
  endforeach()
  set(${out} ${_ok} PARENT_SCOPE)
endfunction()

# _semver_tag_version(<tag> <out>): X.Y.Z of a vX.Y.Z tag, empty for any other tag (pre-releases included)
function(_semver_tag_version tag out)
  if(tag MATCHES "^v?([0-9]+\\.[0-9]+\\.[0-9]+)$")
    set(${out} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  else()
    set(${out} "" PARENT_SCOPE)
  endif()
endfunction()

# _semver_newest(<url> <out_tag> <out_version> <requirement>...): the newest vX.Y.Z tag on <url> that satisfies one
# of the requirements
function(_semver_newest url out_tag out_version)
  if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "SemVerDeps: git is needed to list the releases of ${url}")
  endif()
  # GIT_TERMINAL_PROMPT=0: fail instead of waiting for a password nobody types
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env GIT_TERMINAL_PROMPT=0 "${GIT_EXECUTABLE}" ls-remote --tags --refs "${url}"
    OUTPUT_VARIABLE _refs ERROR_VARIABLE _err RESULT_VARIABLE _rc TIMEOUT 60)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "SemVerDeps: cannot list the tags of ${url} (${_rc}):\n${_err}")
  endif()
  set(_best_tag "")
  set(_best_version "")
  string(REPLACE "\n" ";" _lines "${_refs}")
  foreach(_line IN LISTS _lines)
    if(_line MATCHES "refs/tags/(.+)$")
      set(_tag "${CMAKE_MATCH_1}")
      _semver_tag_version("${_tag}" _version)
      if(_version)
        _semver_satisfies("${_version}" _ok ${ARGN})
        if(_ok AND (_best_version STREQUAL "" OR _version VERSION_GREATER _best_version))
          set(_best_tag "${_tag}")
          set(_best_version "${_version}")
        endif()
      endif()
    endif()
  endforeach()
  set(${out_tag} "${_best_tag}" PARENT_SCOPE)
  set(${out_version} "${_best_version}" PARENT_SCOPE)
endfunction()

# _semver_read_lock(): the entries of the lock file, once per configure, as the global properties
# SEMVER_LOCKED_<NAME> (the tag) and SEMVER_LOCKED_NAMES
function(_semver_read_lock)
  get_property(_read GLOBAL PROPERTY SEMVER_LOCK_READ)
  if(_read)
    return()
  endif()
  set_property(GLOBAL PROPERTY SEMVER_LOCK_READ TRUE)
  if(NOT EXISTS "${SEMVER_LOCK_FILE}")
    return()
  endif()
  # Editing or pulling the lock file configures again
  set_property(DIRECTORY "${CMAKE_SOURCE_DIR}" APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${SEMVER_LOCK_FILE}")
  file(STRINGS "${SEMVER_LOCK_FILE}" _lines REGEX "^[^#]")
  foreach(_line IN LISTS _lines)
    if(_line MATCHES "^[ \t]*([A-Za-z0-9_.+-]+)[ \t]+([^ \t]+)[ \t]*$")
      string(TOUPPER "${CMAKE_MATCH_1}" _N)
      set_property(GLOBAL APPEND PROPERTY SEMVER_LOCKED_NAMES "${CMAKE_MATCH_1}")
      set_property(GLOBAL PROPERTY SEMVER_LOCKED_${_N} "${CMAKE_MATCH_2}")
    else()
      message(WARNING "SemVerDeps: ${SEMVER_LOCK_FILE}: line '${_line}' is not '<name> <tag>', ignored")
    endif()
  endforeach()
endfunction()

# _semver_submodules(<dir> <out_all> <out_missing> <path>...): the submodules of the checkout <dir> among <path>... (all
# of them without paths), and those of them not populated yet
function(_semver_submodules dir out_all out_missing)
  execute_process(COMMAND "${GIT_EXECUTABLE}" submodule status -- ${ARGN} WORKING_DIRECTORY "${dir}"
                  OUTPUT_VARIABLE _status ERROR_QUIET)
  set(_all "")
  set(_missing "")
  string(REPLACE "\n" ";" _lines "${_status}")
  foreach(_line IN LISTS _lines)
    # <state><commit> <path>[ (<describe>)]: state - for a submodule not populated
    if(_line MATCHES "^(.)[0-9a-f]+ ([^ ]+)")
      list(APPEND _all "${CMAKE_MATCH_2}")
      if(CMAKE_MATCH_1 STREQUAL "-")
        list(APPEND _missing "${CMAKE_MATCH_2}")
      endif()
    endif()
  endforeach()
  set(${out_all} "${_all}" PARENT_SCOPE)
  set(${out_missing} "${_missing}" PARENT_SCOPE)
endfunction()

# _semver_checkout(<name> <dir> <url> <tag> <submodules_given> <submodule>...): the release <tag> of <url> in <dir>
# (see SOURCE_DIR at the top)
function(_semver_checkout name dir url tag submodules_given)
  if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "SemVerDeps: git is needed to check out ${name} in ${dir}")
  endif()
  string(TOUPPER "${name}" _N)
  # The network, as FetchContent: no fetch in an existing checkout, of a tag or of submodules, with
  # FETCHCONTENT_FULLY_DISCONNECTED or FETCHCONTENT_UPDATES_DISCONNECTED(_<NAME>), and no clone with the first one. A
  # move to a tag that the checkout has needs none
  set(_no_fetch "")
  foreach(_switch FETCHCONTENT_FULLY_DISCONNECTED FETCHCONTENT_UPDATES_DISCONNECTED
                  FETCHCONTENT_UPDATES_DISCONNECTED_${_N})
    if(${_switch} AND NOT _no_fetch)
      set(_no_fetch ${_switch})
    endif()
  endforeach()
  # The submodules, as FetchContent: all of them, or those of GIT_SUBMODULES (none for GIT_SUBMODULES "")
  set(_with_submodules TRUE)
  if(submodules_given AND NOT ARGN)
    set(_with_submodules FALSE)
  endif()
  set(_git ${CMAKE_COMMAND} -E env GIT_TERMINAL_PROMPT=0 "${GIT_EXECUTABLE}" -c advice.detachedHead=false)
  set(_moved FALSE)
  if(EXISTS "${dir}/.git")
    execute_process(COMMAND ${_git} rev-parse HEAD WORKING_DIRECTORY "${dir}"
                    OUTPUT_VARIABLE _head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    execute_process(COMMAND ${_git} rev-parse --verify --quiet "refs/tags/${tag}^{commit}" WORKING_DIRECTORY "${dir}"
                    OUTPUT_VARIABLE _wanted OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT (_head AND _head STREQUAL _wanted))
      # At another commit: moved only from a detached HEAD without changes, those of its submodules included
      # (symbolic-ref succeeds on a branch)
      execute_process(COMMAND ${_git} symbolic-ref --quiet HEAD WORKING_DIRECTORY "${dir}"
                      OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _branch_rc)
      execute_process(COMMAND ${_git} status --porcelain --untracked-files=no WORKING_DIRECTORY "${dir}"
                      OUTPUT_VARIABLE _changes ERROR_QUIET)
      if(_branch_rc EQUAL 0 OR _changes)
        message(FATAL_ERROR "SemVerDeps: ${dir} is not at ${name} ${tag}, the release selected, and it is on a branch "
                            "or has changes, which the configure does not touch. Commit or stash them and check the "
                            "tag out (git -C ${dir} checkout ${tag}), or build the folder as it is with\n"
                            "  -DFETCHCONTENT_SOURCE_DIR_${_N}=${dir}")
      endif()
      set(_all "")
      if(_with_submodules)
        _semver_submodules("${dir}" _all _missing ${ARGN})
      endif()
      if(_no_fetch AND (NOT _wanted OR _all))
        message(FATAL_ERROR "SemVerDeps: ${dir} is not at ${name} ${tag}, the release selected, and moving it needs "
                            "the remote (for the tag or the submodules), which ${_no_fetch} forbids")
      endif()
      message(STATUS "${name}: moving ${dir} to ${tag}")
      if(NOT _wanted)
        execute_process(COMMAND ${_git} fetch --quiet --depth 1 "${url}" tag "${tag}" WORKING_DIRECTORY "${dir}"
                        RESULT_VARIABLE _rc ERROR_VARIABLE _err)
        if(NOT _rc EQUAL 0)
          message(FATAL_ERROR "SemVerDeps: cannot fetch ${name} ${tag} from ${url}:\n${_err}")
        endif()
      endif()
      execute_process(COMMAND ${_git} checkout --quiet --detach "${tag}" WORKING_DIRECTORY "${dir}"
                      RESULT_VARIABLE _rc ERROR_VARIABLE _err)
      if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "SemVerDeps: cannot check ${name} ${tag} out in ${dir}:\n${_err}")
      endif()
      set(_moved TRUE)
    endif()
  else()
    file(GLOB _content "${dir}/*")
    if(_content)
      message(FATAL_ERROR "SemVerDeps: ${dir}, the folder of ${name}, is not empty and is not a git checkout")
    endif()
    if(FETCHCONTENT_FULLY_DISCONNECTED)
      message(FATAL_ERROR "SemVerDeps: ${name} is not in ${dir}, and FETCHCONTENT_FULLY_DISCONNECTED forbids cloning it")
    endif()
    message(STATUS "${name}: cloning ${tag} in ${dir}")
    execute_process(COMMAND ${_git} clone --quiet --depth 1 --branch "${tag}" "${url}" "${dir}"
                    RESULT_VARIABLE _rc ERROR_VARIABLE _err)
    if(NOT _rc EQUAL 0)
      message(FATAL_ERROR "SemVerDeps: cannot clone ${name} ${tag} from ${url}:\n${_err}")
    endif()
    # The first clone takes its submodules too, as FetchContent's download
    set(_moved TRUE)
    set(_no_fetch "")
  endif()
  if(NOT _with_submodules)
    return()
  endif()
  # After a clone or a move, all the submodules at the commits of the release (those of a moved checkout had no
  # changes); otherwise those not populated yet, at every configure, and the others keep their commit and changes
  _semver_submodules("${dir}" _all _missing ${ARGN})
  if(_moved)
    set(_update "${_all}")
  else()
    set(_update "${_missing}")
  endif()
  if(NOT _update)
    return()
  endif()
  if(_no_fetch)
    list(JOIN _update ", " _shown)
    message(FATAL_ERROR "SemVerDeps: the submodules ${_shown} of ${dir} are not populated, and ${_no_fetch} forbids "
                        "cloning them")
  endif()
  execute_process(COMMAND ${_git} submodule update --quiet --init --recursive --depth 1 -- ${_update}
                  WORKING_DIRECTORY "${dir}" RESULT_VARIABLE _rc ERROR_VARIABLE _err)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "SemVerDeps: cannot initialize the submodules of ${name} in ${dir}:\n${_err}")
  endif()
endfunction()

function(semver_dependency name)
  # The arguments by index: an option of FetchContent_Declare after VERSION is not a requirement, and an empty value
  # (GIT_SUBMODULES "") reaches FetchContent_Declare
  set(_arg_GIT_REPOSITORY "")
  set(_arg_GIT_TAG "")
  set(_arg_SOURCE_DIR "")
  set(_arg_VERSION "")
  set(_forward "")
  # GIT_SUBMODULES apart: forwarded to FetchContent, or applied by _semver_checkout with SOURCE_DIR
  set(_submodules_given FALSE)
  set(_submodules "")
  set(_submodules_forward "")
  set(_i 1)
  while(_i LESS ARGC)
    set(_a "${ARGV${_i}}")
    math(EXPR _i "${_i} + 1")
    if(_a STREQUAL "GIT_REPOSITORY" OR _a STREQUAL "GIT_TAG" OR _a STREQUAL "SOURCE_DIR")
      if(NOT _i LESS ARGC)
        message(FATAL_ERROR "semver_dependency(${name}): ${_a} has no value")
      endif()
      set(_arg_${_a} "${ARGV${_i}}")
      math(EXPR _i "${_i} + 1")
    elseif(_a STREQUAL "VERSION")
      while(_i LESS ARGC AND "${ARGV${_i}}" MATCHES "^v?[0-9]+(\\.[0-9]+)?(\\.[0-9]+)?$")
        list(APPEND _arg_VERSION "${ARGV${_i}}")
        math(EXPR _i "${_i} + 1")
      endwhile()
    elseif(_a STREQUAL "GIT_SUBMODULES")
      # Its values: up to the next option, an upper-case word
      set(_submodules_given TRUE)
      string(APPEND _submodules_forward " [==[${_a}]==]")
      while(_i LESS ARGC AND NOT "${ARGV${_i}}" MATCHES "^[A-Z][A-Z0-9_]*$")
        list(APPEND _submodules "${ARGV${_i}}")
        string(APPEND _submodules_forward " [==[${ARGV${_i}}]==]")
        math(EXPR _i "${_i} + 1")
      endwhile()
    else()
      string(APPEND _forward " [==[${_a}]==]")
    endif()
  endwhile()
  if(NOT _arg_GIT_REPOSITORY)
    message(FATAL_ERROR "semver_dependency(${name}): GIT_REPOSITORY is missing")
  endif()
  if((_arg_VERSION AND _arg_GIT_TAG) OR NOT (_arg_VERSION OR _arg_GIT_TAG))
    message(FATAL_ERROR "semver_dependency(${name}): give either VERSION or GIT_TAG")
  endif()
  string(TOUPPER "${name}" _N)
  if(_arg_VERSION)
    list(TRANSFORM _arg_VERSION PREPEND "^" OUTPUT_VARIABLE _required)
    list(JOIN _required " || " _required)
  else()
    set(_required "${_arg_GIT_TAG}")
  endif()
  set_property(GLOBAL APPEND PROPERTY SEMVER_${_N}_REQUIRED_BY "${PROJECT_NAME} ${_required}")

  # -------- Already selected in this build tree: check it --------
  get_property(_selected GLOBAL PROPERTY SEMVER_${_N}_SOURCE SET)
  if(_selected)
    if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
      message(WARNING "${name} was selected by a dependency before the top-level project ${PROJECT_NAME} required it: "
                      "declare every dependency of the top-level project before making any available")
    endif()
    get_property(_tag GLOBAL PROPERTY SEMVER_${_N}_TAG)
    get_property(_version GLOBAL PROPERTY SEMVER_${_N}_VERSION)
    get_property(_source GLOBAL PROPERTY SEMVER_${_N}_SOURCE)
    get_property(_by GLOBAL PROPERTY SEMVER_${_N}_BY)
    if(_arg_GIT_TAG)
      if(NOT _tag STREQUAL _arg_GIT_TAG)
        message(WARNING "${PROJECT_NAME} pins ${name} ${_arg_GIT_TAG}; the build uses ${name} ${_tag}, chosen by ${_by}")
      endif()
      return()
    endif()
    if(NOT _version)
      message(STATUS "${name}: cannot check ${_tag} against ${PROJECT_NAME}'s requirement ${_required}")
      return()
    endif()
    _semver_satisfies("${_version}" _ok ${_arg_VERSION})
    if(_ok)
      return()
    endif()
    set(_what "${name} ${_tag} does not satisfy ${PROJECT_NAME}, which requires ${_required}.\n"
              "It was chosen by ${_by}")
    if(_source STREQUAL "local")
      message(WARNING ${_what} " (the checkout FETCHCONTENT_SOURCE_DIR_${_N}).")
    elseif(_source STREQUAL "pinned")
      message(FATAL_ERROR ${_what} ", which pins it.")
    elseif(_source STREQUAL "lock")
      message(FATAL_ERROR ${_what} " from ${SEMVER_LOCK_FILE}. If a newer release satisfies both, select it with\n"
                          "  cmake -DSEMVER_UPDATE=${name} <build-dir>")
    else()
      message(FATAL_ERROR ${_what} ", and it is the newest release that ${_by} accepts: the two projects need "
                          "different releases of ${name}. Build with versions of ${PROJECT_NAME} and ${_by} that "
                          "require compatible releases of ${name}.")
    endif()
    return()
  endif()

  # -------- First in this build tree: select the release --------
  _semver_read_lock()
  get_property(_locked GLOBAL PROPERTY SEMVER_LOCKED_${_N})
  set(_update FALSE)
  string(TOUPPER "${SEMVER_UPDATE}" _update_list)
  if(_update_list MATCHES "^(ALL|ON|TRUE|YES|1)$" OR _N IN_LIST _update_list)
    set(_update TRUE)
  endif()

  set(_tag "")
  set(_version "")
  if(FETCHCONTENT_SOURCE_DIR_${_N})
    set(_source "local")
    set(_tag "${FETCHCONTENT_SOURCE_DIR_${_N}}")
    if(GIT_EXECUTABLE)
      # The newest vX.Y.Z tag merged into HEAD: the nearest tag (git describe) may be a pre-release, such as
      # v2.2.0-rc1, with a release before it
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" tag --merged HEAD --list "v[0-9]*"
        WORKING_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_${_N}}"
        OUTPUT_VARIABLE _merged OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _rc)
      if(_rc EQUAL 0)
        string(REPLACE "\n" ";" _merged "${_merged}")
        foreach(_t IN LISTS _merged)
          _semver_tag_version("${_t}" _v)
          if(_v AND (_version STREQUAL "" OR _v VERSION_GREATER _version))
            set(_version "${_v}")
            set(_tag "${_t}+ (${FETCHCONTENT_SOURCE_DIR_${_N}})")
          endif()
        endforeach()
      endif()
    endif()
    if(_version AND _arg_VERSION)
      _semver_satisfies("${_version}" _ok ${_arg_VERSION})
      if(NOT _ok)
        message(WARNING "${name} ${_tag} does not satisfy ${PROJECT_NAME}, which requires ${_required}")
      endif()
    endif()
  elseif(_arg_GIT_TAG)
    set(_source "pinned")
    set(_tag "${_arg_GIT_TAG}")
    _semver_tag_version("${_tag}" _version)
  else()
    if(_locked AND NOT _update)
      _semver_tag_version("${_locked}" _version)
      if(_version)
        _semver_satisfies("${_version}" _ok ${_arg_VERSION})
      endif()
      if(_version AND _ok)
        set(_source "lock")
        set(_tag "${_locked}")
      else()
        set(_version "")
        message(STATUS "${name}: ${_locked} of the lock file does not satisfy ${PROJECT_NAME}'s ${_required}")
      endif()
    endif()
    if(NOT _tag)
      if(FETCHCONTENT_FULLY_DISCONNECTED)
        message(FATAL_ERROR "${name}: no release in the lock file for ${_required}, and FETCHCONTENT_FULLY_DISCONNECTED "
                            "forbids looking for one on ${_arg_GIT_REPOSITORY}")
      endif()
      _semver_newest("${_arg_GIT_REPOSITORY}" _tag _version ${_arg_VERSION})
      if(NOT _tag)
        message(FATAL_ERROR "${name}: no release on ${_arg_GIT_REPOSITORY} satisfies ${PROJECT_NAME}'s ${_required}")
      endif()
      set(_source "remote")
    endif()
  endif()

  message(STATUS "${name} ${_tag} (from ${_source}; ${PROJECT_NAME} requires ${_required})")
  set_property(GLOBAL APPEND PROPERTY SEMVER_NAMES "${name}")
  set_property(GLOBAL PROPERTY SEMVER_${_N}_TAG "${_tag}")
  set_property(GLOBAL PROPERTY SEMVER_${_N}_VERSION "${_version}")
  set_property(GLOBAL PROPERTY SEMVER_${_N}_SOURCE "${_source}")
  set_property(GLOBAL PROPERTY SEMVER_${_N}_BY "${PROJECT_NAME}")
  if(_arg_SOURCE_DIR AND NOT _source STREQUAL "local")
    # Checked out by SemVerDeps: FetchContent only adds the folder
    get_filename_component(_dir "${_arg_SOURCE_DIR}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    _semver_checkout("${name}" "${_dir}" "${_arg_GIT_REPOSITORY}" "${_tag}" ${_submodules_given} ${_submodules})
    cmake_language(EVAL CODE "FetchContent_Declare([==[${name}]==] SOURCE_DIR [==[${_dir}]==]${_forward})")
  else()
    # With a local checkout (FETCHCONTENT_SOURCE_DIR_<NAME>) FetchContent does not use GIT_TAG
    if(NOT _source STREQUAL "local")
      set(_forward " GIT_TAG [==[${_tag}]==] GIT_SHALLOW TRUE${_forward}")
    endif()
    cmake_language(EVAL CODE "FetchContent_Declare([==[${name}]==] GIT_REPOSITORY [==[${_arg_GIT_REPOSITORY}]==]${_forward}${_submodules_forward})")
  endif()

  get_property(_scheduled GLOBAL PROPERTY SEMVER_FINALIZE_SCHEDULED)
  if(NOT _scheduled)
    set_property(GLOBAL PROPERTY SEMVER_FINALIZE_SCHEDULED TRUE)
    cmake_language(DEFER DIRECTORY "${CMAKE_SOURCE_DIR}" CALL _semver_finalize)
  endif()
endfunction()

# _semver_finalize(): at the end of the configure, the summary, and the lock file with the releases selected
function(_semver_finalize)
  get_property(_names GLOBAL PROPERTY SEMVER_NAMES)
  get_property(_locked_names GLOBAL PROPERTY SEMVER_LOCKED_NAMES)
  message(STATUS "Dependencies (name, release, from, required by):")
  set(_entries "")
  foreach(_name IN LISTS _names)
    string(TOUPPER "${_name}" _N)
    get_property(_tag GLOBAL PROPERTY SEMVER_${_N}_TAG)
    get_property(_source GLOBAL PROPERTY SEMVER_${_N}_SOURCE)
    get_property(_required_by GLOBAL PROPERTY SEMVER_${_N}_REQUIRED_BY)
    list(JOIN _required_by ", " _required_by)
    message(STATUS "  ${_name} ${_tag}, ${_source}, ${_required_by}")
    if(_source STREQUAL "lock" OR _source STREQUAL "remote")
      list(APPEND _entries "${_name} ${_tag}")
    endif()
  endforeach()
  # Entries of releases not selected in this configure stay (a dependency switched off by an option keeps its
  # release), and so does the entry of a dependency taken from a local checkout
  foreach(_name IN LISTS _locked_names)
    string(TOUPPER "${_name}" _N)
    get_property(_source GLOBAL PROPERTY SEMVER_${_N}_SOURCE)
    if(NOT (_source STREQUAL "lock" OR _source STREQUAL "remote"))
      get_property(_tag GLOBAL PROPERTY SEMVER_LOCKED_${_N})
      list(APPEND _entries "${_name} ${_tag}")
    endif()
  endforeach()
  if(_entries)
    list(SORT _entries CASE INSENSITIVE)
    list(JOIN _entries "\n" _entries)
    string(CONCAT _content "# Releases of the dependencies, selected by cmake/SemVerDeps.cmake: <name> <tag>\n"
                           "# cmake -DSEMVER_UPDATE=ALL (or =<name>) <build-dir> selects the newest ones again\n"
                           "${_entries}\n")
    set(_old "")
    if(EXISTS "${SEMVER_LOCK_FILE}")
      file(READ "${SEMVER_LOCK_FILE}" _old)
    endif()
    if(NOT _old STREQUAL _content)
      file(WRITE "${SEMVER_LOCK_FILE}" "${_content}")
      message(STATUS "Dependencies: ${SEMVER_LOCK_FILE} written")
    endif()
  endif()
  if(NOT SEMVER_UPDATE STREQUAL "OFF")
    set(SEMVER_UPDATE OFF CACHE STRING
        "Dependencies to select again from their remotes at the next configure: ALL, or a list of names" FORCE)
  endif()
endfunction()
