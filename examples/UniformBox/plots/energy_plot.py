import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

n_snaps = 102
directory = '../'
snap_name_base = "Isothermal_"
plot_dir = "random_orbits_1d"
plot_filename_base = "snap"

v_rot = 200.
r_0 = 100.

e_kin_array = []
e_pot_array = []
e_total_array = []
snap_numbers = np.linspace(0,n_snaps)
for i in range(n_snaps):
    snap_number = "%03d" %i
    filename = directory + snap_name_base + snap_number + ".hdf5"
    f = h5.File(filename)

    header = f["Header"]
    n_parts = header.attrs["NumPart_ThisFile"][1]
    pos = np.array(f["PartType1/Coordinates"])
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    r = np.sqrt(x**2 + y**2 + z**2)
    e_pot_particles = v_rot**2 * np.log(r/r_0)
    e_pot_total = np.sum(e_pot_particles)
    e_pot_array = np.append(e_pot_array,e_pot_total)
    v = np.array(f["PartType1/Velocities"])
    e_kin_particles = 0.5*v**2
    e_kin_total = np.sum(e_kin_particles)
    e_kin_array = np.append(e_kin_array,e_kin_total)
    e_total_array = np.append(e_total_array,e_kin_total+e_pot_total)

print snap_numbers.shape
print e_kin_array.shape
plt.plot(e_kin_array,label = "Kinetic Energy")
plt.plot(e_pot_array,label = "Potential Energy")
plt.plot(e_total_array,label = "Total Energy")
plt.legend(loc = "upper right")    
plot_filename = "%s/energy_plot.png" %plot_dir
plt.savefig(plot_filename,format = "png")
plt.close()
