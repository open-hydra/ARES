import numpy as np

# =======================================================================
# Build ARES / ATLAS-compatible real-fluid lookup tables
# (thermo.dat, transport.dat, phase.txt) from a simple constant-property
# fluid model.
#
# This mirrors the file format produced by ATLAS' real_fluid module
# (src/GPB/real_fluid/io.py), i.e. a Tecplot DATAPACKING=BLOCK table on a
# (pressure, enthalpy) grid, so that the resulting files can be read
# directly by ARES.
# =======================================================================

# -----------------------------------------------------------------------
# Fluid model (constant properties) -- values from make-fluid-table.py
# -----------------------------------------------------------------------
name  = ""              # filename prefix (e.g. "gas-"); "" -> thermo.dat ...
fluid = "Water"         # phase descriptor written in phase.txt

rho = 1002.400000000     # density            [kg/m^3]
cp  = 0.4148300000E+04  # specific heat      [J/(kg K)]
mu  = 0.8906400000E-03  # dynamic viscosity  [Pa s]
Pr  = 6.033               # Prandtl number     [-]
k   = mu * cp / Pr      # thermal conductivity [W/(m K)]
#k   = 5.278             # thermal conductivity [W/(m K)]
a   = 1500.0            # speed of sound     [m/s] (constant; incompressible)

# -----------------------------------------------------------------------
# Table bounds and resolution
# -----------------------------------------------------------------------
pmin = 1e5             # min pressure       [Pa]
pmax = 100e5            # max pressure       [Pa]
Tmin = 300.0            # min temperature    [K]
Tmax = 600.0            # max temperature    [K]

NP = 600                 # number of pressure points
NH = 600                 # number of enthalpy points

T_ref = 300.0          # enthalpy reference temperature [K]

# -----------------------------------------------------------------------
# Convert (T, p) bounds to enthalpy bounds, build (p, h) grid
# -----------------------------------------------------------------------
# Constant-cp model:  h = cp * (T - T_ref)  ->  T = T_ref + h / cp
hmin = cp * (Tmin - T_ref)
hmax = cp * (Tmax - T_ref)

# Reference state values to subtract (zero h and s at Tmin, pmin)
h_ref = hmin
s_ref = cp * np.log(Tmin)

p = np.linspace(pmin, pmax, NP)
h = np.linspace(hmin, hmax, NH)

# -----------------------------------------------------------------------
# Compute properties on the grid (p fastest -> matches I=NP, J=NH)
# -----------------------------------------------------------------------
thermo = []
transport = []

for j in range(NH):
    for i in range(NP):
        T = T_ref + h[j] / cp

        rho_  = rho
        drdT  = 0.0                      # incompressible / constant rho
        drdh  = 0.0
        drdp  = 0.0
        cp_   = cp
        s     = cp * np.log(T) - s_ref   # shifted to zero at Tmin
        ss    = a
        hh    = h[j] - h_ref             # shifted to zero at hmin

        thermo.append([p[i], hh, rho_, T, drdT, drdh, cp_, s, drdp, ss])
        transport.append([p[i], hh, mu, k])


# -----------------------------------------------------------------------
# Writers (identical layout to ATLAS src/GPB/real_fluid/io.py)
# -----------------------------------------------------------------------
def write_phase(name, fluid):
    with open(name + "phase.txt", 'w') as f:
        f.write("real-fluid phase\n")
        f.write(f"{fluid}\n")


def write_block_table(filename, title, variables, varloc, NP, NH, rows):
    with open(filename, 'w') as f:
        f.write(f'TITLE = "{title}"\n')
        f.write("VARIABLES = " + ", ".join(f'"{v}"' for v in variables) + "\n")
        f.write(f"ZONE T=real-gas, I={NP}, J={NH}, K=1, "
                f"DATAPACKING=BLOCK, VARLOCATION=({varloc}=NODAL)\n")
        n_vars = len(rows[0])
        for col in range(n_vars):
            for row in rows:
                f.write(f"{row[col]}\n")


def write_thermo(name, NP, NH, rows):
    write_block_table(
        name + "thermo.dat",
        "Mass Thermodynamic Properties",
        ["Pressure", "Enthalpy", "Density", "Temperature", "dRho/dT",
         "dRho/dh", "Cp", "Entropy", "dRho/dp", "SoundSpeed"],
        "[1-10]", NP, NH, rows)


def write_transport(name, NP, NH, rows):
    write_block_table(
        name + "transport.dat",
        "Transport Properties",
        ["Pressure", "Enthalpy", "Viscosity", "Conductivity"],
        "[1-4]", NP, NH, rows)


# -----------------------------------------------------------------------
# Write output
# -----------------------------------------------------------------------
write_phase(name, fluid)
write_thermo(name, NP, NH, thermo)
write_transport(name, NP, NH, transport)

print(f"Wrote {name}phase.txt, {name}thermo.dat and {name}transport.dat "
      f"({NP} x {NH} grid)")
