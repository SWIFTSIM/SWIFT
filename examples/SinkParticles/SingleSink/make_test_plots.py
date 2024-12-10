################################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

'''
Make some very simple plots showing the evolution of the single sink.

Run with python make_test_plots.py <snapshot-directory-name>

'''


import numpy as np
from sys import argv
from glob import glob
import matplotlib.pyplot as plt
import unyt
import h5py

params = {"text.usetex": True, 'axes.labelsize': 16, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'lines.linewidth' : 2, 'axes.titlesize' : 16, 'font.family' : 'serif'}
plt.rcParams.update(params)

# Units
G_cgs = 6.6743e-8 * unyt.cm**3 * unyt.g**-1 * unyt.s**-2
unit_mass_cgs = 1.988409870698051e+43 * unyt.g
unit_density_cgs = 6.767905323247329e-22 * unyt.Unit('g/cm**3')
unit_velocity_cgs = (1.*unyt.Unit("km/s")).to('cm/s')


# Basic Bondi-Hoyle prediction from the sink's starting mass in a constant medium
def simple_bondi_hoyle(t, m, rho, v, cs):

    m_sink = np.zeros(len(t),dtype=np.float64) * unyt.g

    m_sink[0] = m[0] * unit_mass_cgs

    rho_0 = rho[0] * unit_density_cgs
    v_0 = v[0] * unit_velocity_cgs
    cs_0 = cs[0] * unit_velocity_cgs

    timestep = ((t[1] - t[0]) * unyt.Gyr).to('s')

    for i in range(len(times)-1):

        numerator = 4. * np.pi * G_cgs**2 * m_sink[i]**2 * rho_0
        denominator = np.power(v_0**2 + cs_0**2, 3./2.)

        m_sink[i+1] = m_sink[i] + (numerator/denominator) * timestep

    return m_sink.to('Msun')


snapshots = glob(f'./{argv[1]}/*.hdf5')

times = np.empty(len(snapshots),dtype=np.float32)
mass_evo = np.empty(len(snapshots),dtype=np.float32)
subgrid_mass_evo = np.empty(len(snapshots),dtype=np.float32)
rho_evo = np.empty(len(snapshots),dtype=np.float32)
v_evo = np.empty(len(snapshots),dtype=np.float32)
cs_evo = np.empty(len(snapshots),dtype=np.float32)
hsml_evo = np.empty(len(snapshots),dtype=np.float32)


for s, snap in enumerate(snapshots):

    with h5py.File(snap) as f:
        times[s] = f['Header'].attrs['Time'][0]
        mass_evo[s] = f['PartType3/Masses'][0]
        subgrid_mass_evo[s] = f['PartType3/SubgridMasses'][0]
        rho_evo[s] = f['PartType3/GasDensities'][0]
        cs_evo[s] = f['PartType3/GasSoundSpeeds'][0]
        hsml_evo[s] = f['PartType3/SmoothingLengths'][0]

        v = f['PartType3/GasVelocities'][0] - f['PartType3/Velocities'][0]
        v_evo[s] = np.sqrt(v[0]**2+v[1]**2+v[2]**2)


# Run the Bondi-Hoyle prediction
m_sink_bondi_prediction = simple_bondi_hoyle(times,mass_evo,rho_evo,v_evo,cs_evo)

# Normalise time to go up to 1
t = times/times[-1]

# Time evolution of the sink mass, target mass, and Bondi-Hoyle prediction
fig, ax = plt.subplots(1,figsize=(8,6))
ax.plot(t,np.log10(m_sink_bondi_prediction),label=r'Simple Bondi-Hoyle, constant $\rho$, $v$, $c_{\rm s}$')
ax.plot(t,np.log10(subgrid_mass_evo * 1e10),label=r'Subgrid mass')
ax.plot(t,np.log10(mass_evo * 1e10),label=r'Dynamical mass')

ax.set_xlabel(r"$t/t_{\rm final}$")
ax.set_ylabel(r"$\log_{10}(M_{\rm sink})$")
ax.legend(fontsize=14)
ax.set_ylim(6.68,7.)
fig.savefig('mass_evolution.png',bbox_inches='tight')

# Time evolution of the total mass accreted relative to the target mass
fig, ax = plt.subplots(1,figsize=(8,6))
ax.plot(t,mass_evo/subgrid_mass_evo)
ax.set_xlabel(r"$t/t_{\rm final}$")
ax.set_ylabel(r"$m_{\rm accr}^{\rm gas}/m_{\rm target}^{\rm gas}$")
fig.savefig('target_mass_ratio.png',bbox_inches='tight')

# Time evolution of the smoothing length
fig, ax = plt.subplots(1,figsize=(8,6))
ax.plot(t,hsml_evo)
ax.set_xlabel(r"$t/t_{\rm final}$")
ax.set_ylabel(r"$h_{\rm sink}$ [kpc]")
fig.savefig('smoothing_length.png',bbox_inches='tight')
