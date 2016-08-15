import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

n_snaps = 101
directory = '../'
snap_name_base = "uniformBox_"
plot_dir = "./"
cooling_lambda = 1.0
u_floor = 0.0
delta_time = 0.1
plot_title = r"$\Lambda \, = \, %1.3f \, \, \, u_{floor} = %1.3f$" %(cooling_lambda,u_floor)
plot_filename = "energy_plot_creasey_lambda_1p0eminus50.png"
e_kin_array=[]
e_therm_array=[]
e_total_array = []
time = np.zeros(n_snaps-1)
analytic_solution = np.zeros(n_snaps-1)
for i in range(1,n_snaps):
    snap_number = "%03d" %i
    filename = directory + snap_name_base + snap_number + ".hdf5"
    f = h5.File(filename)
    print "Reading snap number %d" %i
    header = f["Header"]
    n_parts = header.attrs["NumPart_ThisFile"][0]
    time[i-1] = i*delta_time
    u = np.array(f["PartType0/InternalEnergy"])
    u_total = np.sum(u)
    e_therm_array = np.append(e_therm_array,u_total)
    v = np.array(f["PartType0/Velocities"])
    e_kin_particles = 0.5*v**2
    e_kin_total = np.sum(e_kin_particles)
    e_kin_array = np.append(e_kin_array,e_kin_total)
    e_total_array = np.append(e_total_array,e_kin_total+u_total)
    analytic_solution[i-1] = e_total_array[0] - n_parts*(time[i-1] - time[0])*cooling_lambda


# plt.plot(e_kin_array,label = "Kinetic Energy")
# plt.plot(e_therm_array,label = "Internal Energy")
plt.plot(time,e_total_array,'r',label = "Total Energy")
plt.plot(time,analytic_solution,'--k',label = "Analytic Solution")
plt.xlabel("Time (seconds)")
plt.ylabel("Energy (erg/g)")
plt.ylim((0,1000))
plt.legend(loc = "upper right")    
plt.title(plot_title)
full_plot_filename = "%s/%s" %(plot_dir,plot_filename)
plt.savefig(full_plot_filename,format = "png")
plt.close()
