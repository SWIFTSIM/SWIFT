import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from swiftsimio import load

import unyt
import glob

from tqdm import tqdm

plt.rcParams.update({
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.titlesize": 20,
    "savefig.dpi": 100
})

snapdir  = 'snapshots/'          # directory containing the snapshots
savefile = 'EnergyEvolution.png' # name of the saved image

sn_time =  0.00125               # time in internal units between injections

def get_thermal_energy (data, unit):
    """returns the total thermal energy of the particles
       E_th = sum (u_i * m_i)
         
      data:  the data snapshot
      unit:  the energy unit to convert to"""
    u = data.gas.internal_energies.to_physical()
    m = data.gas.masses.to_physical()

    E = np.sum(u*m)

    E.convert_to_units(unit)

    return E

def get_radiated_energy (data, unit):
    """returns the total radiated energy
       cumulated sum over time on particle basis

       data:  the data snapshot
       unit:  the energy unit to convert to"""
    try:
        Er = np.sum(data.gas.radiated_energies)

    except AttributeError:
        Er = 0 * unyt.erg

    Er.convert_to_units(unit)

    return Er

def get_kinetic_energy (data, unit):
    """returns the total kinetic energy of the particles
       E_k = sum (m_i * (v_i)**2 / 2)

       data:  the data snapshot
       unit:  the energy unit to convert to"""
    v_vec = data.gas.velocities.to_physical()
    v = np.sqrt(np.sum(np.square(v_vec), axis=1))

    m = data.gas.masses.to_physical()

    E_kin_part = m * np.square(v) / 2
    E_kin = np.sum(E_kin_part)

    E_kin.convert_to_units(unit)

    return E_kin

def get_injected_energy (data, unit):
    """returns the total injected energy
       cumulated sum over time on particle basis

       data:  the data snapshot
       unit:  the energy unit to convert to"""
    E_inj_part = data.gas.forcing_injected_energies
    
    E_inj = np.sum(E_inj_part)

    E_inj.convert_to_units(unit)

    return E_inj

def get_time (data, unit):
    """returns the time of the snapshot

       data: the data snapshot
       unit: the time unit to convert to"""
    t = data.metadata.time
    t.convert_to_units(unit)
    return t

def get_energy_evolution (path):
    """Loops over all the snapshots in the directory to obtain the
       energy in the different components. Returns a dictionary
       with arrays of the time evolution of the different components

       path:  the directory containing the snapshots"""
    # specify the units we are using
    ut = unyt.Myr
    uE = 1e53 * unyt.erg

    # get the sorted filenames 
    fnames = sorted(glob.glob(path + 'BK_????.hdf5'))

    # initialise the dictionary for storage
    data_dict = {'times': np.zeros(len(fnames)),
                 'E_th':  np.zeros(len(fnames)),
                 'E_k':   np.zeros(len(fnames)),
                 'E_r':   np.zeros(len(fnames)),
                 'E_inj': np.zeros(len(fnames)),
                 'u_t':   ut,
                 'u_E':   uE}

    for i in tqdm(range(len(fnames)), desc='looping over snapshots'):
        snapshot = fnames[i]
        data = load(snapshot)

        # get the different components
        t = get_time(data, ut)

        E_th  = get_thermal_energy(data, uE)
        E_k   = get_kinetic_energy(data, uE)
        E_r   = get_radiated_energy(data, uE)
        E_inj = get_injected_energy(data, uE)
        
        # store the components in the dictionary
        data_dict['times'][i] = t.value
        data_dict['E_th'][i]  = E_th.value
        data_dict['E_k'][i]   = E_k.value
        data_dict['E_r'][i]   = E_r.value
        data_dict['E_inj'][i] = E_inj.value

    # compute some auxiliary data
    sn_rate = 1 / (sn_time * data.metadata.units.time)  # injection rate
    sn_rate.convert_to_units(1/unyt.Myr)

    # amount of injections happened at specified times
    data_dict['sn_injections'] = sn_rate.value * data_dict['times'] 

    # total energy of the system plus the radiated energy
    # this minus the injected energy should be constant
    data_dict['E_tot'] = data_dict['E_th'] + data_dict['E_k'] + \
                         data_dict['E_r']

    return data_dict

def plot_energies (data_dict):
    """Plots the energies nicely and saves the image"""
    fig = plt.figure()

    ax = fig.add_subplot(111)

    # we want a top axis showing the amount of injections there have been
    topax = ax.twiny()  

    # plot the lines
    ax.plot(data_dict['times'], data_dict['E_tot']-data_dict['E_inj'], color='black',
            label=r'$E_\mathrm{tot}-E_\mathrm{inj}$')

    ax.plot(data_dict['times'], data_dict['E_th'], label=r'$E_\mathrm{th}$')
    ax.plot(data_dict['times'], data_dict['E_k'], label=r'$E_\mathrm{k}$')
    ax.plot(data_dict['times'], data_dict['E_r'], label=r'$E_r$')
    ax.plot(data_dict['times'], data_dict['E_inj'], label=r'$E_\mathrm{inj}$')

    ax.legend(facecolor='white', framealpha=1)

    # limits
    ax.set_xlim(0, data_dict['times'].max())
    ax.set_ylim(1e-2, 6e1)
    topax.set_xlim(0, data_dict['sn_injections'].max())
    
    # labels
    ax.set_xlabel(r'$t \ \mathrm{[Myr]}$')
    ax.set_ylabel(r'$E_\mathrm{X} \ \mathrm{[10^{53} \ erg]}$')
    topax.set_xlabel('injections')
    
    # ticks and grids
    topax.set_xticks([0,5,10,15,20,25,30])

    topax.xaxis.set_minor_locator(AutoMinorLocator())
    topax.grid(visible=True, which='major', axis='x', alpha=0.5)
    topax.grid(visible=True, which='minor', axis='x', alpha=0.2)

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(visible=True, which='major', axis='y', alpha=0.5)
    ax.grid(visible=True, which='minor', axis='y', alpha=0.2)

    # ordering
    ax.set_zorder(1)
    ax.set_frame_on(False)
    topax.set_frame_on(True)

    # scales  
    ax.set_yscale('log')

    # saving and done
    fig.savefig(savefile)
    
# obtain all the data
data_dict = get_energy_evolution(snapdir)

# and plot it nicely
plot_energies(data_dict)




