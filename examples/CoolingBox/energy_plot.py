import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import sys

filename = "./energy.txt"
plot_dir = "./"
cooling_lambda = 1.0e-22
n_H = 0.1
T_init = 1.0e5
T_floor = 1.0e4
plot_title = r"$\Lambda \, = \, %1.1g \mathrm{erg}\mathrm{cm^3}\mathrm{s^{-1}} \, \, T_{init} = %1.1g\mathrm{K} \, \, T_{floor} = %1.1g\mathrm{K} \, \, n_H = %1.1g\mathrm{cm^{-3}}$" %(cooling_lambda,T_init,T_floor,n_H)
plot_filename = "energy_plot_creasey_no_cooling_T_init_1p0e5_n_H_0p1.png"
#analytic_solution = np.zeros(n_snaps-1)
array = np.genfromtxt(filename,skip_header = 1)
time = array[:,0]
total_energy = array[:,2]


# plt.plot(e_kin_array,label = "Kinetic Energy")
# plt.plot(e_therm_array,label = "Internal Energy")
plt.plot(time,total_energy,'r',label = "Total Energy")
#plt.plot(time,analytic_solution,'--k',label = "Analytic Solution")
plt.xlabel("Time (code units)")
plt.ylabel("Energy (code units)")
#plt.ylim((0,1.0e8))
plt.legend(loc = "upper right")    
#plt.title(plot_title)
full_plot_filename = "%s/%s" %(plot_dir,plot_filename)
if (int(sys.argv[1])==0):
    plt.show()
else:
    plt.savefig(full_plot_filename,format = "png")
    plt.close()
