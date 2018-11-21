#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import scipy.optimize as sco



Nmax = 2001
steps = 10

iterarray = np.arange(0,Nmax+1,steps)
Lxtot = np.zeros(len(iterarray))
Lytot = np.zeros(len(iterarray))
Lztot = np.zeros(len(iterarray))
Ltot = np.zeros(len(iterarray))
time_array = np.zeros(len(iterarray))


for i in np.arange(0,Nmax+1,steps):
    print(i)
    f = h5py.File('output_%04d.hdf5'%i, 'r')

    boxsize = f['Header'].attrs['BoxSize']/2.
    
    time_array[int(i/steps)] = f['Header'].attrs['Time']

    particles = f['PartType4']
    coordinates = particles['Coordinates'][:,:]
    velocities = particles['Velocities'][:,:]
    masses = particles['Masses'][:]

    R = ((coordinates[:,0]-boxsize[0])**2 + (coordinates[:,1]-boxsize[1])**2)**.5
    X = np.abs(coordinates[:,0]-boxsize[0])
    Y = np.abs(coordinates[:,1]-boxsize[1])
    Z = np.abs(coordinates[:,2]-boxsize[2])

    vx = velocities[:,0]
    vy = velocities[:,1]
    vz = velocities[:,2]

    Lx = (Y*vz - Z*vy)*masses
    Ly = (Z*vx - X*vz)*masses
    Lz = (X*vy - Y*vx)*masses
    
    L = (Lx**2 + Ly**2 + Lz**2)**.5

    Lxtot[int(i/steps)] = np.sum(Lx)
    Lytot[int(i/steps)] = np.sum(Ly)
    Lztot[int(i/steps)] = np.sum(Lz)
    Ltot[int(i/steps)] = np.sum(L)

time_array[-1]=2.
plt.plot(time_array,Lxtot/Lxtot[0]-1,label='Lx total')    
plt.plot(time_array,Lytot/Lytot[0]-1,label='Ly total')    
plt.plot(time_array,Lztot/Lztot[0]-1,label='Lz total')    
plt.plot(time_array,Ltot/Ltot[0]-1,label='L total')    
plt.xlabel('Time')
plt.ylabel('ratio between current and zero angular momentum')
plt.legend()
plt.show() 

plt.semilogy(time_array, np.absolute(Ltot/Ltot[0]-1))
plt.xlabel('Time (Gyr)')
plt.ylabel('Fractional change of total angular momentum')
plt.savefig('Total_angular_momentum.png')
plt.show()
plt.close()







