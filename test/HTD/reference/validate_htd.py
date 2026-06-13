import numpy as np
import matplotlib.pyplot as pl
import re
import pandas as pd

# DATI NASA ----- TEST 24-1027 --------------
"""EXPERIMENTAL HEAT-TRANSFER RESULTS FOR CRYOGENIC
HYDROGEN FLOWING IN TUBES AT SUBCRITICAL AND
SUPERCRITICAL PRESSURES TO 800 POUNDS
PER SQUARE INCH ABSOLUTE"""


def read_tecplot_dat(filepath: str, zone: str = None) -> pd.DataFrame:
    """
    Legge un file .dat in formato Tecplot ASCII e restituisce
    un DataFrame pandas con tutte le variabili come colonne.

    Parametri
    ----------
    filepath : str
        Percorso del file .dat da leggere.
    zone : str, optional
        Nome della zona da leggere (es. "B1" o "B2").
        Se None, legge e concatena tutte le zone.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    var_names = []
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("variables"):
            raw = re.sub(r"(?i)variables\s*=\s*", "", line.strip())
            var_names = [v.strip().strip('"') for v in raw.split(",")]
            data_start = i + 1
            break

    if not var_names:
        raise ValueError(f"Nessuna riga 'variables' trovata in: {filepath}")

    zones = {}
    current_zone = "__default__"
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.lower().startswith("zone"):
            match = re.search(r'T\s*=\s*"([^"]+)"', stripped, re.IGNORECASE)
            current_zone = match.group(1) if match else stripped
            zones.setdefault(current_zone, [])
        elif stripped.lower().startswith(("title", "text", "geometry")):
            continue
        else:
            try:
                values = [float(v) for v in stripped.split()]
                if len(values) == len(var_names):
                    zones.setdefault(current_zone, []).append(values)
            except ValueError:
                continue

    if zone is not None:
        if zone not in zones:
            raise KeyError(f"Zona '{zone}' non trovata. Zone disponibili: {list(zones.keys())}")
        rows = zones[zone]
    else:
        rows = []
        for z_rows in zones.values():
            rows.extend(z_rows)

    return pd.DataFrame(rows, columns=var_names)


def read_thermo_block(filepath: str):
    """
    Legge thermo.dat in formato Tecplot BLOCK (DATAPACKING=BLOCK).

    Layout: I-axis = pressione (varia piu' veloce), J-axis = entalpia (piu' lenta).
    Variabili nell'ordine: Pressure, Enthalpy, Density, Temperature, dRho/dT,
                           dRho/dh, Cp, Entropy, dRho/dp, SoundSpeed.

    Parametri
    ----------
    filepath : str
        Percorso del file thermo.dat.

    Ritorna
    -------
    pmin, dp, NP, hmin, dh, NH : float/int
        Estremi e passi della griglia (p,h); NP e NH sono il numero di nodi.
    Ttab : ndarray shape (NP, NH)
        Temperatura tabulata, Ttab[i,j] = T(p_i, h_j).
    """
    # Trova la riga ZONE e leggi I= (NP), J= (NH)
    n_header = 0
    NP = NH = 0
    with open(filepath, 'r') as f:
        for line in f:
            n_header += 1
            if re.search(r'\bZONE\b', line, re.IGNORECASE):
                m_i = re.search(r'\bI\s*=\s*(\d+)', line, re.IGNORECASE)
                m_j = re.search(r'\bJ\s*=\s*(\d+)', line, re.IGNORECASE)
                if m_i:
                    NP = int(m_i.group(1))
                if m_j:
                    NH = int(m_j.group(1))
                break

    if NP == 0 or NH == 0:
        raise ValueError(f"Impossibile leggere I= / J= dall'header ZONE di {filepath}")

    total = NP * NH

    # Legge i blocchi: Pressure(0), Enthalpy(1), Density(2-skip), Temperature(3)
    # max_rows limita la lettura ai primi 4 blocchi (evita di caricare tutti i 10 blocchi)
    all_4 = np.loadtxt(filepath, skiprows=n_header, max_rows=4 * total, dtype=np.float64)

    p_block = all_4[0 * total : 1 * total]
    h_block = all_4[1 * total : 2 * total]
    T_block = all_4[3 * total : 4 * total]

    pmin = p_block[0]
    dp   = p_block[1] - p_block[0]   # p varia piu' veloce (I-axis): p_block[1] = p(i=1,j=0)
    hmin = h_block[0]
    dh   = h_block[NP] - h_block[0]  # h cambia dopo NP valori: h_block[NP] = h(i=0,j=1)

    # Reshape: dati scritti j-esterno, i-interno → (NH, NP).T = (NP, NH)
    Ttab = T_block.reshape(NH, NP).T

    return pmin, dp, NP, hmin, dh, NH, Ttab


def T_val(p, pmin, dp, h, hmin, dh, Ttab):
    """Interpolazione bilineare di Ttab(NP,NH) al punto (p, h)."""
    i = int((p - pmin) / dp)
    j = int((h - hmin) / dh)

    Vij   = Ttab[i,   j  ]
    Viij  = Ttab[i+1, j  ]
    Vijj  = Ttab[i,   j+1]
    Viijj = Ttab[i+1, j+1]

    A = Vij
    B = (Viij  - Vij) / dp
    C = (Vijj  - Vij) / dh
    E = (Vij + Viijj - Viij - Vijj) / (dp * dh)

    dist_p = p - pmin - i * dp
    dist_h = h - hmin - j * dh

    return A + B * dist_p + C * dist_h + E * dist_p * dist_h


# ----------------------------------------------------------------
# Legge tabelle termodinamiche dalla cartella INPUT (BLOCK format)
# ----------------------------------------------------------------
print("Lettura thermo.dat ...")
pmin, dp, NP, hmin, dh, NH, Ttab = read_thermo_block('../INPUT/thermo.dat')
pmax = pmin + (NP - 1) * dp
hmax = hmin + (NH - 1) * dh
print(f"  p: {pmin:.4e} .. {pmax:.4e}  dp={dp:.4e}  NP={NP}")
print(f"  h: {hmin:.4e} .. {hmax:.4e}  dh={dh:.4e}  NH={NH}")

# ----------------------------------------------------------------
# Dati sperimentali NASA
# ----------------------------------------------------------------
x = np.array([2.00, 3.50, 5.00, 6.50, 7.25, 8.00, 8.75, 10.00, 11.50, 13.00, 14.50, 16.00])  # inches
xm = x * 0.0254  # metri

Tb_Rank = np.array([64.3, 69.2, 73.7, 77.9, 79.9, 81.8, 83.7, 87.0, 91.0, 95.2, 99.7, 104.6])
Tb_NASA = Tb_Rank * 5 / 9
Tw_Rank = np.array([495.0, 623.0, 684.0, 678.0, 646.0, 648.0, 589.0, 561.0, 522.0, 498.0, 473.0, 461.0])
Tw_NASA = Tw_Rank * 5 / 9


# ----------------------------------------------------------------
# Legge dati CFD (ARES) da 1d.dat — B1 (ingresso) + B2 (tubo riscaldato)
# ----------------------------------------------------------------
df_cfd = read_tecplot_dat("../OUTPUT/1d.dat", zone=None)

x_CFD  = df_cfd["x[m]"].to_numpy()
Tw_CFD = df_cfd["Tw[K]"].to_numpy()
Tb_CFD = df_cfd["T[K]"].to_numpy()

# Temperatura bulk alternativa dalla tabella (decommentare se necessario)
# hb_CFD = df_cfd["h[J/kg]"].to_numpy()
# Pave   = df_cfd["p[Pa]"].to_numpy()
# Tb_CFD = np.array([T_val(Pave[k], pmin, dp, hb_CFD[k], hmin, dh, Ttab)
#                    for k in range(len(hb_CFD))])

# ----------------------------------------------------------------
# Plot: Temperatura Bulk
# ----------------------------------------------------------------
pl.figure()
pl.minorticks_on()
pl.plot(x_CFD, Tb_CFD, '-', color='black', label='ARES')
pl.plot(xm, Tb_NASA, 'o', label='Experimental', markerfacecolor='white', markeredgecolor='red', markersize=8)
pl.grid(True, which='both', linestyle='-', alpha=0.3)
pl.xlim(0, 0.45)
pl.ylim(30, 60)
pl.xlabel('Tube Length [m]')
pl.ylabel('Bulk Temperature [K]')
pl.title('Bulk Temperature along the Tube')
pl.legend()

# ----------------------------------------------------------------
# Plot: Temperatura a Parete
# ----------------------------------------------------------------
pl.figure()
ax = pl.gca()
ax.set_box_aspect(1)
pl.minorticks_on()
pl.plot(x_CFD, Tw_CFD, '-', color='black', label='ARES')
pl.plot(xm, Tw_NASA, 'o', label='Experimental', markerfacecolor='white', markeredgecolor='red', markersize=8)
pl.grid(True, which='both', linestyle='-', alpha=0.3)
pl.xlim(0, 0.45)
pl.ylim(50, 450)
pl.xlabel('Tube Length [m]')
pl.ylabel('Wall Temperature [K]')
pl.title('Wall Temperature along the Tube')
pl.legend()


pl.show()
