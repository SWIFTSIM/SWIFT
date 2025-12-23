import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from swiftsimio import load

import unyt
import glob
import argparse

from tqdm import tqdm

parser = argparse.ArgumentParser(description='plot the energy evolution')
parser.add_argument('snapdir', type=str,
                    help='the directory containing the snapshots')
parser.add_argument('savefile', type=str,
                    help='the basename of the saved image')
parser.add_argument('--snapbasename', type=str, nargs='?', default='BK',
                    help='the basename of the snapshots')
args = parser.parse_args()

plt.rcParams.update({
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.titlesize": 20,
    "savefig.dpi": 100
})

def get_magnetic_energy (data, unit):
    try:
        B_vec = data.gas.magnetic_flux_densities.to_physical()

        B = np.sqrt(np.sum(np.square(B_vec), axis=1))

        E_B_part = B**2 / (2 * unyt.mu_0)

        n_part = data.gas.densities / data.gas.masses

        E_B = np.sum( E_B_part / n_part )

    except AttributeError:
        E_B = np.nan * unyt.erg

    E_B.convert_to_units(unit)

    return E_B

def get_thermal_energy (data, unit):
    u = data.gas.internal_energies
    m = data.gas.masses

    E = np.sum(u*m)

    E.convert_to_units(unit)

    return E

def get_radiated_energy (data, unit):
    try:
        Er = np.sum(data.gas.radiated_energies)

    except AttributeError:
        Er = 0 * unyt.erg

    Er.convert_to_units(unit)

    return Er

def get_kinetic_energy (data, unit):
    v_vec = data.gas.velocities.to_physical()
    v = np.sqrt(np.sum(np.square(v_vec), axis=1))

    m = data.gas.masses.to_physical()

    E_kin_part = m * np.square(v) / 2
    E_kin = np.sum(E_kin_part)

    E_kin.convert_to_units(unit)

    return E_kin

def get_dedner_energy (data, unit):
    try:
        psi_ch = data.gas.dedner_scalars_over_cleaning_speeds.to_physical()
        rho = data.gas.densities.to_physical()
        m   = data.gas.masses

        #inserting mu_0 in the expression for some reason gives errors
        u_psi = psi_ch**2 / (2 * rho)
        u_psi /= unyt.mu_0 

        E_psi = np.sum(u_psi * m)
    except AttributeError:
        E_psi = np.nan * unyt.erg

    E_psi.convert_to_units(unit)

    return E_psi

def get_injected_energy (data, unit):
    E_inj_part = data.gas.forcing_injected_energies
    
    E_inj = np.sum(E_inj_part)

    E_inj.convert_to_units(unit)

    return E_inj

def get_time (data, unit):
    t = data.metadata.time
    t.convert_to_units(unit)
    return t

def get_energy_evolution (path):
    """path contains the snapshots"""
    ut = unyt.Myr
    uE = 1e53 * unyt.erg

    fnames = sorted(glob.glob(path + args.snapbasename + '_????.hdf5'))
    data_dict = {'times': np.zeros(len(fnames)),
                 'E_B':   np.zeros(len(fnames)),
                 'E_psi': np.zeros(len(fnames)),
                 'E_th':  np.zeros(len(fnames)),
                 'E_k':   np.zeros(len(fnames)),
                 'E_r':   np.zeros(len(fnames)),
                 'E_inj': np.zeros(len(fnames)),
                 'u_t':   ut,
                 'u_E':   uE}

    for i in tqdm(range(len(fnames)), desc='looping over snapshots'):
        snapshot = fnames[i]
        data = load(snapshot)

        t = get_time(data, ut)

        E_B   = get_magnetic_energy(data, uE)
        E_psi = get_dedner_energy(data, uE)
        E_th  = get_thermal_energy(data, uE)
        E_k   = get_kinetic_energy(data, uE)
        E_r   = get_radiated_energy(data, uE)
        E_inj = get_injected_energy(data, uE)

        data_dict['times'][i] = t.value
        data_dict['E_B'][i]   = E_B.value
        data_dict['E_psi'][i] = E_psi.value
        data_dict['E_th'][i]  = E_th.value
        data_dict['E_k'][i]   = E_k.value
        data_dict['E_r'][i]   = E_r.value
        data_dict['E_inj'][i] = E_inj.value

    data_dict['E_tot'] = data_dict['E_B'] + data_dict['E_th'] + data_dict['E_k'] + \
                         data_dict['E_r'] + data_dict['E_psi']

    return data_dict

data_dict = get_energy_evolution(args.snapdir)

ratio = (data_dict['E_tot'] - data_dict['E_inj']) / data_dict['E_tot'][0]
for i in range(len(ratio)):
    print(ratio[i])

SN_rate = 1 / (0.00125 * 2.5577e1)  # in Myr^-1

SN_injections = data_dict['times'] * SN_rate

### Plotting ###

fig = plt.figure()

ax = fig.add_subplot(111)
topax = ax.twiny()

ax.plot(data_dict['times'], data_dict['E_tot']-data_dict['E_inj'], color='black',
        label=r'$E_\mathrm{tot}-E_\mathrm{inj}$')

ax.plot(data_dict['times'], data_dict['E_th'], label=r'$E_\mathrm{th}$')
ax.plot(data_dict['times'], data_dict['E_k'], label=r'$E_\mathrm{k}$')
ax.plot(data_dict['times'], data_dict['E_B'], label=r'$E_B$')
ax.plot(data_dict['times'], data_dict['E_r'], label=r'$E_r$')
ax.plot(data_dict['times'], data_dict['E_psi'], label=r'$E_\psi$')
ax.plot(data_dict['times'], data_dict['E_inj'], label=r'$E_\mathrm{inj}$')

ax.set_xlim(0, data_dict['times'].max())

ax.set_xlabel(r'$t \ \mathrm{[Myr]}$')
ax.set_ylabel(r'$E_\mathrm{X} \ \mathrm{[10^{53} \ erg]}$')

topax.set_xlim(0, SN_injections.max())
topax.set_xticks([0,5,10,15,20,25,30])
topax.set_xlabel('SN injection')

topax.xaxis.set_minor_locator(AutoMinorLocator())
topax.grid(visible=True, which='major', axis='x', alpha=0.5)
topax.grid(visible=True, which='minor', axis='x', alpha=0.2)

ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.grid(visible=True, which='major', axis='y', alpha=0.5)
ax.grid(visible=True, which='minor', axis='y', alpha=0.2)

ax.set_zorder(1)
ax.set_frame_on(False)
topax.set_frame_on(True)

ax.legend(facecolor='white', framealpha=1)

ax.set_yscale('log')
ax.set_ylim(1e-7)

fig.savefig(args.savefile)

