# Plots contribution to cooling rates from each of the different metals
# based on cooling_output.dat and cooling_element_*.dat files produced
# by testCooling. 

import matplotlib.pyplot as plt
import numpy as np

k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27
erg_to_J = 1.0e-7
gas_gamma=5.0/3.0

# Primoridal ean molecular weight as a function of temperature
def mu(T):
    H_frac=0.752
    T_trans=10000.0
    if T > T_trans:
        return 4. / (8. - 5. * (1. - H_frac))
    else:
        return 4. / (1. + 3. * H_frac)

# Temperature of some primoridal gas with a given internal energy
def T(u):
    H_frac=0.752
    T_trans=10000.0
    T_over_mu = (gas_gamma - 1.) * u * mH_in_kg / k_in_J_K
    ret = np.ones(np.size(u)) * T_trans

    # Enough energy to be ionized?
    mask_ionized = (T_over_mu > (T_trans+1) / mu(T_trans+1))
    if np.sum(mask_ionized)  > 0:
        ret[mask_ionized] = T_over_mu[mask_ionized] * mu(T_trans*10)

    # Neutral gas?
    mask_neutral = (T_over_mu < (T_trans-1) / mu((T_trans-1)))
    if np.sum(mask_neutral)  > 0:
        ret[mask_neutral] = T_over_mu[mask_neutral] * mu(0)

    return ret

# Number of metals tracked by EAGLE cooling
elements = 11

# Declare arrays of internal energy and cooling rate
u = [[] for i in range(elements+1)]
cooling_rate = [[] for i in range(elements+1)]
Temperature = [[] for i in range(elements+1)]

# Read in total cooling rate
#file_in = open('cooling_output.dat', 'r')
#for line in file_in:
#	data = line.split()
#	u.append(float(data[0]))
#	cooling_rate[0].append(-float(data[1]))
#
#file_in.close()

n_iz = 10
for iz in range(n_iz):
	file_in = open('cooling_element_'+str(iz)+'.dat','r')
	z = float(file_in.readline())
	for line in file_in:
		data = line.split() 
		u[iz+1].append(float(data[0]))
		cooling_rate[iz+1].append(-float(data[1]))
	file_in.close()
	a = 1.0/(1.0 + z)
	#u[iz+1] = np.asarray(u[iz+1]) / a**(3 * (gas_gamma - 1.))
	#Temperature[iz+1] = T(u[iz+1] * erg_to_J)

# Plot
ax = plt.subplot(111)
#p0, = plt.loglog(u, cooling_rate[0], linewidth = 0.5, color = 'k', label = 'Total')
p1,  = plt.loglog(u[1], cooling_rate[1], linewidth = 0.5, color = 'k',             label = 'z = 10')
p2,  = plt.loglog(u[2], cooling_rate[2], linewidth = 0.5, color = 'b',             label = 'z = 9')
p3,  = plt.loglog(u[3], cooling_rate[3], linewidth = 0.5, color = 'g',             label = 'z = 8')
p4,  = plt.loglog(u[4], cooling_rate[4], linewidth = 0.5, color = 'r',             label = 'z = 7')
p5,  = plt.loglog(u[5], cooling_rate[5], linewidth = 0.5, color = 'c',             label = 'z = 6')
p6,  = plt.loglog(u[6], cooling_rate[6], linewidth = 0.5, color = 'm',             label = 'z = 5')
p7,  = plt.loglog(u[7], cooling_rate[7], linewidth = 0.5, color = 'y',             label = 'z = 4')
p8,  = plt.loglog(u[8], cooling_rate[8], linewidth = 0.5, color = 'lightgray',     label = 'z = 3')
p9,  = plt.loglog(u[9], cooling_rate[9], linewidth = 0.5, color = 'olive',         label = 'z = 2')
p10, = plt.loglog(u[10], cooling_rate[10], linewidth = 0.5, color = 'saddlebrown', label = 'z = 1')
p11, = plt.loglog(u[1],  -np.asarray(cooling_rate[1]),  linewidth = 0.5, linestyle = '--', color = 'k')
p12, = plt.loglog(u[2],  -np.asarray(cooling_rate[2]),  linewidth = 0.5, linestyle = '--', color = 'b')
p13, = plt.loglog(u[3],  -np.asarray(cooling_rate[3]),  linewidth = 0.5, linestyle = '--', color = 'g')
p14, = plt.loglog(u[4],  -np.asarray(cooling_rate[4]),  linewidth = 0.5, linestyle = '--', color = 'r')
p15, = plt.loglog(u[5],  -np.asarray(cooling_rate[5]),  linewidth = 0.5, linestyle = '--', color = 'c')
p16, = plt.loglog(u[6],  -np.asarray(cooling_rate[6]),  linewidth = 0.5, linestyle = '--', color = 'm')
p17, = plt.loglog(u[7],  -np.asarray(cooling_rate[7]),  linewidth = 0.5, linestyle = '--', color = 'y')
p18, = plt.loglog(u[8],  -np.asarray(cooling_rate[8]),  linewidth = 0.5, linestyle = '--', color = 'lightgray')
p19, = plt.loglog(u[9],  -np.asarray(cooling_rate[9]),  linewidth = 0.5, linestyle = '--', color = 'olive')
p20, = plt.loglog(u[10], -np.asarray(cooling_rate[10]), linewidth = 0.5, linestyle = '--', color = 'saddlebrown')
ax.set_position([0.15,0.15,0.75,0.75])
#plt.xlim([1e10,1e17])
#plt.ylim([1e-24,1e-21])
plt.xlabel("Internal energy ${\\rm{[erg \cdot g^{-1}]}}$", fontsize = 14)
plt.ylabel("${\Lambda/n_h^2 }$ ${\\rm{[erg \cdot cm^3 \cdot s^{-1}]}}$", fontsize = 14)
plt.legend(handles = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10])
plt.show()

# Read in contributions to cooling rates from each of the elements
#for elem in range(elements):
#	file_in = open('cooling_element_'+str(elem)+'.dat','r')
#	for line in file_in:
#        	data = line.split()
#        	cooling_rate[elem+1].append(-float(data[0]))
#	file_in.close()

## Plot
#ax = plt.subplot(111)
#p0, = plt.loglog(u, cooling_rate[0], linewidth = 0.5, color = 'k', label = 'Total')
#p1, = plt.loglog(u, cooling_rate[1], linewidth = 0.5, color = 'k', linestyle = '--', label = 'H + He')
#p2, = plt.loglog(u, cooling_rate[3], linewidth = 0.5, color = 'b', label = 'Carbon')
#p3, = plt.loglog(u, cooling_rate[4], linewidth = 0.5, color = 'g', label = 'Nitrogen')
#p4, = plt.loglog(u, cooling_rate[5], linewidth = 0.5, color = 'r', label = 'Oxygen')
#p5, = plt.loglog(u, cooling_rate[6], linewidth = 0.5, color = 'c', label = 'Neon')
#p6, = plt.loglog(u, cooling_rate[7], linewidth = 0.5, color = 'm', label = 'Magnesium')
#p7, = plt.loglog(u, cooling_rate[8], linewidth = 0.5, color = 'y', label = 'Silicon')
#p8, = plt.loglog(u, cooling_rate[9], linewidth = 0.5, color = 'lightgray', label = 'Sulphur')
#p9, = plt.loglog(u, cooling_rate[10], linewidth = 0.5, color = 'olive', label = 'Calcium')
#p10, = plt.loglog(u, cooling_rate[11], linewidth = 0.5, color = 'saddlebrown', label = 'Iron')
#ax.set_position([0.15,0.15,0.75,0.75])
#plt.xlim([1e12,1e17])
#plt.ylim([1e-24,1e-21])
#plt.xlabel("Internal energy ${\\rm{[erg \cdot g^{-1}]}}$", fontsize = 14)
#plt.ylabel("${\Lambda/n_h^2 }$ ${\\rm{[erg \cdot cm^3 \cdot s^{-1}]}}$", fontsize = 14)
##plt.legend(handles = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10])
#plt.show()
