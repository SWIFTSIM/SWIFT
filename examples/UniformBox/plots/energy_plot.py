import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

n_snaps = 101
directory = '../'
snap_name_base = "uniformBox_"
plot_dir = "./"
plot_title = r"$\Lambda \, = \, 1.0$"
plot_filename = "energy_plot_lambda_1.png"

e_kin_array = []
e_therm_array = []
e_total_array = []
snap_numbers = np.linspace(0,n_snaps)
for i in range(n_snaps):
    snap_number = "%03d" %i
    filename = directory + snap_name_base + snap_number + ".hdf5"
    f = h5.File(filename)
    print "Reading snap number %d" %i
    header = f["Header"]
    n_parts = header.attrs["NumPart_ThisFile"][1]
    u = np.array(f["PartType0/InternalEnergy"])
    u_total = np.sum(u)
    e_therm_array = np.append(e_therm_array,u_total)
    v = np.array(f["PartType0/Velocities"])
    e_kin_particles = 0.5*v**2
    e_kin_total = np.sum(e_kin_particles)
    e_kin_array = np.append(e_kin_array,e_kin_total)
    e_total_array = np.append(e_total_array,e_kin_total+u_total)

plt.plot(e_kin_array,label = "Kinetic Energy")
plt.plot(e_therm_array,label = "Internal Energy")
plt.plot(e_total_array,label = "Total Energy")
plt.xlabel("Snapshot number")
plt.ylabel("Energy (code units)")
plt.legend(loc = "lower right")    
plt.title(plot_title)
full_plot_filename = "%s/%s" %(plot_dir,plot_filename)
plt.savefig(full_plot_filename,format = "png")
plt.close()
