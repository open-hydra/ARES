######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
ENDIF(BT STREQUAL "RELEASE")

# Intel IPO needs two more things that CMake will not arrange on its own, and
# both have to be re-established on *every* configure: the archiver is carried
# in plain (non-cache) variables that CMakeFortranCompiler.cmake resets to GNU
# ar each time project() runs.  Wrapped in a macro so it can also be applied on
# the early-return path below, where the flags are already cached and the rest
# of this file is skipped.  A macro, not a function: it must set variables in
# the caller's scope.  The flags change in the cache only, through a variable of
# another name: a plain CMAKE_Fortran_FLAGS_RELEASE would hide the cache entry
# (policy CMP0126 is NEW with cmake_minimum_required 3.23), so the flags that
# SET_COMPILE_FLAG adds below would reach the cache but not the targets.  Every
# configure gives the same flags.
MACRO(SET_IPO_ARCHIVER)
IF(CMAKE_Fortran_FLAGS_RELEASE MATCHES "(^| )[-/]Q?ipo( |$)")

    # 1. An IPO-aware archiver.  GNU ar/ranlib cannot index IR objects ("archive
    #    has no index" / "file format not recognized") and xiar merely forwards to
    #    ar, so both fail.  oneAPI ships llvm-ar/llvm-ranlib with the compiler: in
    #    <bin>/compiler from oneAPI 2024 on, in <bin>/../bin-llvm up to oneAPI 2023
    #    (where CMake's own IntelLLVM module looks), <bin> the directory of the
    #    compiler with its symbolic links resolved.  -DIPO_AR=... -DIPO_RANLIB=...
    #    give them by hand, for any other layout.
    GET_FILENAME_COMPONENT(_ipo_bindir "${CMAKE_Fortran_COMPILER}" REALPATH)
    GET_FILENAME_COMPONENT(_ipo_bindir "${_ipo_bindir}" DIRECTORY)
    FIND_PROGRAM(IPO_AR     NAMES llvm-ar
                 HINTS "${_ipo_bindir}/compiler" "${_ipo_bindir}/../bin-llvm" "${_ipo_bindir}")
    FIND_PROGRAM(IPO_RANLIB NAMES llvm-ranlib
                 HINTS "${_ipo_bindir}/compiler" "${_ipo_bindir}/../bin-llvm" "${_ipo_bindir}")

    IF(IPO_AR AND IPO_RANLIB)
        # 2. The link line must carry -ipo too: IPO objects hold IR, not ELF, so the
        #    real optimisation and code generation happen at link time.  Once: the
        #    macro runs at every configure.
        IF(NOT CMAKE_EXE_LINKER_FLAGS_RELEASE MATCHES "(^| )[-/]Q?ipo( |$)")
            SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -ipo"
                CACHE STRING "Set the CMAKE_EXE_LINKER_FLAGS_RELEASE flags" FORCE)
        ENDIF()
        SET(CMAKE_AR     "${IPO_AR}"     CACHE FILEPATH "Archiver used for IPO objects" FORCE)
        SET(CMAKE_RANLIB "${IPO_RANLIB}" CACHE FILEPATH "Ranlib used for IPO objects"   FORCE)
        # CMakeFindBinUtils has already left plain variables of the same names in
        # this scope, and a plain variable shadows the cache entry when CMake
        # expands <CMAKE_AR>/<CMAKE_RANLIB> in the archive rules.  Setting only
        # the cache leaves the generated link.txt still calling /usr/bin/ar, which
        # produces an unindexed archive ("archive has no index").  Set both; these
        # propagate to every subdirectory added after this point.
        SET(CMAKE_AR     "${IPO_AR}")
        SET(CMAKE_RANLIB "${IPO_RANLIB}")
        MESSAGE(STATUS "IPO                : enabled (archiver ${IPO_AR})")
    ELSE()
        MESSAGE(WARNING "IPO was requested and -ipo compiles, but no llvm-ar/llvm-ranlib "
                        "was found with ${CMAKE_Fortran_COMPILER} (in <bin>/compiler, "
                        "<bin>/../bin-llvm or <bin>, <bin> its directory). Disabling IPO: "
                        "a static library of IPO objects cannot be indexed by GNU ar. "
                        "Give the tools with -DIPO_AR=... -DIPO_RANLIB=..., or re-run "
                        "with -DUSE_IPO=OFF to silence this.")
        # -ipo out of the compile and the link flags, every occurrence, in the cache only
        SET(_ipo_flags "${CMAKE_Fortran_FLAGS_RELEASE}")
        WHILE(_ipo_flags MATCHES "(^| )[-/]Q?ipo( |$)")
            STRING(REGEX REPLACE "(^| )([-/]Q?ipo)( |$)" " " _ipo_flags "${_ipo_flags}")
        ENDWHILE()
        SET(CMAKE_Fortran_FLAGS_RELEASE "${_ipo_flags}"
            CACHE STRING "Set the CMAKE_Fortran_FLAGS_RELEASE flags" FORCE)
        SET(_ipo_flags "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
        WHILE(_ipo_flags MATCHES "(^| )[-/]Q?ipo( |$)")
            STRING(REGEX REPLACE "(^| )([-/]Q?ipo)( |$)" " " _ipo_flags "${_ipo_flags}")
        ENDWHILE()
        SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${_ipo_flags}"
            CACHE STRING "Set the CMAKE_EXE_LINKER_FLAGS_RELEASE flags" FORCE)
        UNSET(_ipo_flags)
    ENDIF()
ENDIF()
ENDMACRO(SET_IPO_ARCHIVER)

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    # The flags are cached, but the archiver is not -- redo it before bailing out.
    SET_IPO_ARCHIVER()
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

# Don't add underscores in symbols for C-compatability.
# With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
# MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while this
# flag is set: the MPI library exports the blocks with an underscore (mpipriv1_), ARES refers to them
# without one. Do not use them (src/lib/parallel/Mod_MPI.f90).
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-underscoring")

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         "/QxHost"       # Intel Windows
                         ${GNUNATIVE}    # GNU
                         "-ta=host"      # Portland Group
                )
# Add preprocessor flag
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-cpp"        # Intel or GNU
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                  Fortran  "-ffree-line-length-none"     # Intel o GNU
                           "-ffixed-line-length-none"    # Intel
                           "-extend-source"
                )        
###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g or /debug:full) are already on by default

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran REQUIRED "-O0" # All compilers not on Windows
                                  "/Od" # Intel Windows
                )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all" # Intel
                         "/warn:all" # Intel Windows
                         "-Wall"     # GNU
                                     # Portland Group (on by default)
                )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "/traceback"   # Intel Windows
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                )

# Check array bounds
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-check bounds"  # Intel
                         "/check:bounds"  # Intel Windows
                         "-fcheck=bounds" # GNU (New style)
                         "-fbounds-check" # GNU (Old style)
                         "-Mbounds"       # Portland Group
                )

#####################
### TESTING FLAGS ###
#####################

# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default

# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-funroll-loops" # GNU
                         "-unroll"        # Intel
                         "/unroll"        # Intel Windows
                         "-Munroll"       # Portland Group
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-inline"            # Intel
                         "/Qinline"           # Intel Windows
                         "-finline-functions" # GNU
                         "-Minline"           # Portland Group
                )

# Interprocedural optimization (Intel/PGI only).
#
# Turn it off with -DUSE_IPO=OFF (install.sh --no-ipo) if a machine's toolchain
# misbehaves.  Nothing here is assumed: the flag is probed with a link step, and
# if either the probe or the archiver lookup below fails the build quietly falls
# back to a plain RELEASE rather than producing something that dies at link time.
#
# GNU's -flto is deliberately not offered: it is fragile with mixed-language
# static libraries, MPI wrappers, and certain linker configurations.
IF(USE_IPO)
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                     Fortran LINK
                             "-ipo"              # Intel
                             "/Qipo"             # Intel Windows
                             "-Mipa=fast,inline" # Portland Group
                    )
ENDIF(USE_IPO)

SET_IPO_ARCHIVER()

# Single-file optimizations
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran "-ip"  # Intel
#                         "/Qip" # Intel Windows
#                )

# Vectorize code
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran "-vec-report0"  # Intel
#                         "/Qvec-report0" # Intel Windows
#                         "-Mvect"        # Portland Group
#                )
