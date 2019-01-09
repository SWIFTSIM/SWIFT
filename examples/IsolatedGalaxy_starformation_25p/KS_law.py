#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt

# for the moment only use a single snapshot
i = 29 

# Load the snapshot
#f = h5py.File('output_%04dcopy.hdf5'%i, 'r')
f = h5py.File('output_%04d.hdf5'%i, 'r')

# Load the stars
starparticle = f['PartType4']

# Load the properties of the stars
coordinates = starparticle['Coordinates'][:,:]
newstarflag = starparticle['NewStarFlag'][:]
boxsize = f['Header'].attrs['BoxSize']
masses_star = starparticle['Masses'][:]

# Load the positions of the star particles and only take
# the new star particles
x = coordinates[:,0] - boxsize[0]/2.
y = coordinates[:,1] - boxsize[1]/2.
r = (x**2 + y**2)**.5
birthdensity = starparticle['BirthDensity'][:]
birthdensity=birthdensity[newstarflag==1]
r_new = r[newstarflag==1]
masses_new = masses_star[newstarflag==1]

# Load the gas particles
gasparticle = f['PartType0']

# Load the properties of the gas particles
coordinates_gas = gasparticle['Coordinates'][:,:]
masses_gas = gasparticle['Masses'][:]

# Load the positions of the gas particles
x_gas = coordinates_gas[:,0] - boxsize[0]/2.
y_gas = coordinates_gas[:,0] - boxsize[1]/2.
r_gas = (x_gas**2 + y_gas**2)**.5

r_array = np.linspace(0,15,100)
area = np.zeros(len(r_array)-1)
SigmaSFR = np.zeros(len(r_array)-1)
Sigmagas = np.zeros(len(r_array)-1)
massSFR_array = np.zeros(len(r_array)-1)
massGAS_array = np.zeros(len(r_array)-1)
for i in range(1,len(r_array)):
    area_loop = np.pi*(r_array[i]**2 - r_array[i-1]**2)
    mask = (r_new > r_array[i-1]) & (r_new < r_array[i])
    mask2 = (r_gas > r_array[i-1]) & (r_gas < r_array[i])
    massSFR = np.sum(masses_new[mask])
    massSFR_array[i-1] = massSFR
    massGAS = np.sum(masses_gas[mask2])
    massGAS_array[i-1] = massGAS
    #print('%d  area=%1.4f,  M_SFR=%1.4f,  M_gas=%1.4f'%(i,area_loop,massSFR,massGAS))
    area[i-1] = area_loop
    SigmaSFR[i-1] = massSFR/area_loop
    Sigmagas[i-1] = massGAS/area_loop

# Put SF law in correct units
SigmaSFR *= 1e10/0.29e8
Sigmagas *= 1e10/(1e3**2)
SigmaSFR[SigmaSFR==0] = 1e-6
print(Sigmagas)

Sgas_array =np.logspace(-1,3,40) 
Delta_Sgas = np.log10(Sgas_array[1]) - np.log10(Sgas_array[0])
SigmaSFR_avg = np.zeros(len(Sgas_array))
Sigmagas_avg = np.zeros(len(Sgas_array))

for i in range(0,len(Sgas_array)):
    totgasmass = 0.0
    totnewstarmass = 0.0
    area_loop = 0.0
    for j in range(0,len(Sigmagas)):
        #print('%d %1.4e %1.4e %1.4e %1.4e'%(j,Sigmagas[j],Sgas_array[i],10**(-Delta_Sgas),10**(Delta_Sgas)))
        if (Sigmagas[j]>Sgas_array[i]*10**(-Delta_Sgas)) and (Sigmagas[j]<Sgas_array[i]*10**(Delta_Sgas)):
            totgasmass += massGAS_array[j]
            totnewstarmass += massSFR_array[j]
            area_loop += area[j]
    if area_loop==0.0:
        SigmaSFR_avg[i] = 0.0
        Sigmagas_avg[i] = 0.0
    else:
        SigmaSFR_avg[i] = totnewstarmass/area_loop
        Sigmagas_avg[i] = totgasmass/area_loop

print(SigmaSFR_avg)
print(Sigmagas_avg)

# Put SF law in correct units
SigmaSFR_avg *= 1e10/0.29e8
Sigmagas_avg *= 1e10/(1e6)
SigmaSFR_avg[SigmaSFR_avg==0]=1e-6

def KS_law(Sgas,n=1.4,A=1.515):
    return A*1e-4*(Sgas)**n

Sigmagaslaw = np.logspace(-1,3,100) 
SigmaSFRlaw = KS_law(Sigmagaslaw)

Z = 0
ncrit = 10 # g/cm^3
ncrit = 10 # M_sun /pc^2

# plot SF law
plt.loglog(Sigmagas,SigmaSFR,'o',label='Data')
#plt.loglog(Sigmagas_avg,SigmaSFR_avg,'o',label='binned Data')
plt.loglog(Sigmagaslaw, SigmaSFRlaw,'r-.',label='KS law')
#plt.loglog(np.logspace(-1,3,40),1e-5*np.ones(40),'o')
plt.axvline(x=10,linestyle='--', color='r',label='$\Sigma_c = 10$')
plt.xlabel('$\Sigma_{gas}$')
plt.ylabel('$\Sigma_{SFR}$')
plt.legend()
plt.show()
   


'''
plt.hist(birthdensity,bins=100,range=[0,2.5])
plt.show()

plt.hist(x_gas,bins=100)
plt.show()

plt.hist(r,bins=100)
plt.show()

plt.hist(r_gas,bins=100)
plt.show()
'''
