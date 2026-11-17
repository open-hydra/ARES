#!/usr/bin/env python3
"""Quasi-1D isentropic MDM solution of the planar nozzle (CoolProp HEOS, gas phase imposed), from the total state
of test M1.5H (p0 = 904388 Pa, T0 = 542.13 K): subsonic up to the throat, supersonic past it. Writes the
two-column profiles along x that [ICB-Block1] reads (IC/p.dat, IC/T.dat, IC/u.dat) and prints the inlet mass flux
and static temperature of the [inflow] section (BC 404: g, T), the throat and exit states.

Run in ct-env:  python3 IC/make_ic.py"""
import os
import numpy as np
import CoolProp.CoolProp as CP

p0, T0 = 904388.0, 542.13
here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the case directory

# wall profile (half height = flow area per unit depth and per half) from the ARES mesh
with open(os.path.join(here, 'MESH', 'mesh.tec')) as fh:
    fh.readline(); z = fh.readline()
    ni = int(z.split('I=')[1].split(',')[0]); nj = int(z.split('J=')[1].split(',')[0])
    nk = int(z.split('K=')[1].split(',')[0])
    v = np.loadtxt(fh).reshape(3, nk, nj, ni)
xw, yw = v[0, 0, -1, :], v[1, 0, -1, :]
it = int(np.argmin(yw))

gas = CP.AbstractState('HEOS', 'MDM'); gas.specify_phase(CP.iphase_gas)
gas.update(CP.PT_INPUTS, p0, T0); h0, s0 = gas.hmass(), gas.smass()
ps = np.linspace(0.9999 * p0, 1.0e5, 20000)
T, rho, u, a = (np.empty_like(ps) for _ in range(4))
for n, p in enumerate(ps):
    gas.update(CP.PSmass_INPUTS, p, s0)
    T[n], rho[n], a[n] = gas.T(), gas.rhomass(), gas.speed_sound()
    u[n] = np.sqrt(max(2.0 * (h0 - gas.hmass()), 0.0))
G = rho * u; k = int(np.argmax(G)); Gs = G[k]

def state(A_ratio, supersonic):
    g = Gs / A_ratio
    if supersonic:
        n = k + int(np.argmin(np.abs(G[k:] - g)))
    else:
        n = int(np.argmin(np.abs(G[:k + 1] - g)))
    return ps[n], T[n], u[n], u[n] / a[n]

prof = np.array([state(yw[i] / yw[it], i > it) for i in range(ni)])
os.makedirs(os.path.join(here, 'IC'), exist_ok=True)
for c, name in enumerate(('p', 'T', 'u')):
    np.savetxt(os.path.join(here, 'IC', f'{name}.dat'), np.column_stack([xw, prof[:, c]]), fmt='%.10e')

pin, Tin, uin, Min = prof[0]
print(f'throat  x = {xw[it]:.4f} m  p*/p0 = {ps[k] / p0:.4f}  T* = {T[k]:.2f} K  G* = {Gs:.2f} kg/m2/s')
print(f'inlet   g = {Gs * yw[it] / yw[0]:.3f} kg/m2/s  T = {Tin:.3f} K  p = {pin:.1f} Pa  u = {uin:.3f} m/s  M = {Min:.4f}')
print(f'exit    p = {prof[-1, 0]:.0f} Pa  T = {prof[-1, 1]:.2f} K  M = {prof[-1, 3]:.3f}  p0/p = {p0 / prof[-1, 0]:.3f}')
print(f'inlet   rho = {Gs * yw[it] / yw[0] / uin:.3f} kg/m3')
