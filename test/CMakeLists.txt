# The CI tests of ARES: fast cases with their checks, each a test/<case>/run.sh that exits 0 when every check
# passes. They need python3 (standard library only), not ATLAS nor CoolProp: the real-fluid tables are committed.
#   ctest --test-dir build -L ci --output-on-failure
# The validation cases (flat plates, HTD, Prt-correction, compressible-nozzle) take minutes to hours and are run by
# hand (docs/development/testing.md).

# ares_test(<name> <case folder> <processors> <description> [run.sh arguments])
function(ares_test name dir procs desc)
  add_test(NAME ${name} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/run.sh ${ARGN})
  set_tests_properties(${name} PROPERTIES
    ENVIRONMENT "ARES_BIN=$<TARGET_FILE:${solver}>"
    LABELS "ci"
    PROCESSORS ${procs}
    TIMEOUT 600)
endfunction()

ares_test(RealFluidShockTube real-fluid-shock-tube 1
  "CO2 shock tube vs the exact solution of the CoolProp equation of state")
ares_test(LowMachCylinder low-mach-cylinder 3
  "Preconditioning: Mach-independent Cp of CO2 past a cylinder, prim loses it")

# One test per turbulence model, each in its own runs/<model>: short runs from a turbulent boundary layer (smoke tests,
# not accuracy: the flat-plate validations are test/flat-plate-*). SAcomp is refused by ARES at set-up.
foreach(model SA SA-R SA-RC SA-QCR2000 SA-rough SA-rough-QCR2000 SST SST-RC Wilcox2006
              SSGLRR SSGLRR-SD SSGLRR-w2019 SSGLRR-w2019-SD)
  ares_test(TurbulenceSmoke.${model} turbulence-smoke 1 "300 iterations of ${model} on a coarse flat plate" ${model})
endforeach()
