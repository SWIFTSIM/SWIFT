import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

n_snaps = 41

#some constants
OMEGA = 0.3 # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
h = 0.67777 # hubble parameter
gamma = 5./3.
eta = 1.2349
H_0_cgs = 100. * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

#read some header/parameter information from the first snapshot

filename = "CoolingHalo_0000.hdf5"
f = h5.File(filename,'r')
params = f["Parameters"]
unit_mass_cgs = float(params.attrs["InternalUnitSystem:UnitMass_in_cgs"])
unit_length_cgs = float(params.attrs["InternalUnitSystem:UnitLength_in_cgs"])
unit_velocity_cgs = float(params.attrs["InternalUnitSystem:UnitVelocity_in_cgs"])
unit_time_cgs = unit_length_cgs / unit_velocity_cgs
v_c = float(params.attrs["SoftenedIsothermalPotential:vrot"])
v_c_cgs = v_c * unit_velocity_cgs
header = f["Header"]
N = header.attrs["NumPart_Total"][0]
box_centre = np.array(header.attrs["BoxSize"])

#calculate r_vir and M_vir from v_c
r_vir_cgs = v_c_cgs / (10. * H_0_cgs * np.sqrt(OMEGA))
M_vir_cgs = r_vir_cgs * v_c_cgs**2 / CONST_G_CGS

potential_energy_array = []
internal_energy_array = []
kinetic_energy_array = []
time_array_cgs = []

for i in range(n_snaps):

    filename = "CoolingHalo_%04d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType0/Coordinates"]
    coords = np.array(coords_dset)
#translate coords by centre of box
    header = f["Header"]
    snap_time = header.attrs["Time"]
    snap_time_cgs = snap_time * unit_time_cgs
    time_array_cgs = np.append(time_array_cgs,snap_time_cgs)
    coords[:,0] -= box_centre[0]/2.
    coords[:,1] -= box_centre[1]/2.
    coords[:,2] -= box_centre[2]/2.
    radius = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    radius_cgs = radius*unit_length_cgs
    radius_over_virial_radius = radius_cgs / r_vir_cgs

    r = radius_over_virial_radius
    total_potential_energy = np.sum(v_c**2*np.log(r))
    potential_energy_array = np.append(potential_energy_array,total_potential_energy)

    vels_dset = f["PartType0/Velocities"]
    vels = np.array(vels_dset)
    speed_squared = vels[:,0]**2 + vels[:,1]**2 + vels[:,2]**2
    total_kinetic_energy = 0.5 * np.sum(speed_squared)
    kinetic_energy_array = np.append(kinetic_energy_array,total_kinetic_energy)

    u_dset = f["PartType0/InternalEnergy"]
    u = np.array(u_dset)
    total_internal_energy = np.sum(u)
    internal_energy_array = np.append(internal_energy_array,total_internal_energy)

#put energies in units of v_c^2 and rescale by number of particles

pe = potential_energy_array / (N*v_c**2)
ke = kinetic_energy_array / (N*v_c**2)
ie = internal_energy_array / (N*v_c**2)
te = pe + ke + ie

dyn_time_cgs = r_vir_cgs / v_c_cgs
time_array = time_array_cgs / dyn_time_cgs

plt.plot(time_array,ke,label = "Kinetic Energy")
plt.plot(time_array,pe,label = "Potential Energy")
plt.plot(time_array,ie,label = "Internal Energy")
plt.plot(time_array,te,label = "Total Energy")
plt.legend(loc = "upper right")
plt.xlabel(r"$t / t_{dyn}$")
plt.ylabel(r"$E / v_c^2$")
plt.title(r"$%d \, \, \mathrm{particles} \,,\, v_c = %.1f \, \mathrm{km / s}$" %(N,v_c))
plt.ylim((-4,2))
#plot_filename = "density_profile_%03d.png" %i
plt.show()

