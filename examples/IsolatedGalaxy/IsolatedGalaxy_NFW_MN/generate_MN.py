###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2019 Alejandro Benitez-Llambay (alejandro.b.llambay@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import h5py

G = 4.299e-6

def mcum(r, c):
    a = np.log(1+c*r)-c*r/(1.0+c*r)
    b = np.log(1+c)-c/(1.0+c)
    return a/b

def MN_Vcirc(R, z, Rd=4.0, Md=3e10, Zd=0.4):
    return np.sqrt(G*Md*R**2/(R**2+(Rd+np.sqrt(z**2+Zd**2))**2)**(3./2.))

def get_vcirc_expo(R, Mgas=3e10, Rd=4.0):
    sigma0 = Mgas/(2.0*np.pi*Rd**2)
    sigma  = sigma0*np.exp(-R/Rd)
    y = R/(2.0*Rd)
    I0 = special.i0(y)
    K0 = special.k0(y)
    I1 = special.i1(y)
    K1 = special.k1(y)
    return np.sqrt(4.0*np.pi*G*sigma0*Rd*y**2*(I0*K0-I1*K1))

def dP1(x):
    return x*np.exp(-x)

def dP2(x):
    return 1.0/np.cosh(x)**2

hlittle = 0.704
rhoc = 2.775e2 * hlittle**2

#halo mass and concentration
#-------------------
m200 = 1.5e12 #Msun
c    = 10.0  
r200 = (m200/(200*rhoc*4./3.*np.pi))**(1./3.)
#-------------------

#Gaseous disk
#-----------------
npart = int(1e5)    #This controls the resolution
md    = 7e8
rd    = 0.02*r200
zd    = 0.1*rd
#-----------------

#MN parameters
#-----------------
MN_Md = 3e10 #Msun
MN_Rd = 4.0  #kpc
MN_Zd = 0.4  #kpc
#-----------------

#Outputfile
#------------------------
outputfile = 'MN_ICs.hdf5'  #Output file for SWIFT
#------------------------

x = []
y = []
z = []
phi = []

while(len(x) <= npart):
    n = int(4.5e5)
    rmin = 0
    rmax = 8
    
    zmin = -6
    zmax = 6
    
    r0  = rmin+(rmax-rmin)*np.random.rand(n)
    z0 = zmin+(zmax-zmin)*np.random.rand(n)
    
    ymin = 0
    ymax1 = np.max(dP1(r0))
    ymax2 = np.max(dP2(z0))
    
    y1 = ymin+(ymax1-ymin)*np.random.rand(n)
    y2 = ymin+(ymax2-ymin)*np.random.rand(n)
    
    k, = np.where( (y1 <= dP1(r0)) & (y2 <= dP2(z0) ) )

    phi0 = 2*np.pi*np.random.rand(len(k))
    phi.extend(phi0)
    x.extend(r0[k]*np.cos(phi0))
    y.extend(r0[k]*np.sin(phi0))
    z.extend(z0[k])

sample = np.random.randint(0,len(x), npart)
x = np.asarray(x)[:npart]
y = np.asarray(y)[:npart]
z = np.asarray(z)[:npart]
phi = np.asarray(phi)[:npart]

r = np.sqrt(x**2+y**2+z**2)

#Set positions
pgas = np.ones([npart,3])
pgas[:,0] = x*rd 
pgas[:,1] = y*rd
pgas[:,2] = z*zd
vgas = np.zeros([npart,3])
rhogas = np.zeros(npart)
ugas = np.ones(npart)*0
mgas = np.ones(npart)*md/npart*1e-10

#Set velocities
Rgas = np.sqrt(pgas[:,0]**2+pgas[:,1]**2)
ksort = np.argsort(Rgas)

vcirc_gas = get_vcirc_expo(Rgas, Mgas=md, Rd=rd)   #Approximate the velocity of the disk (better than spherical)
vcirc_dm  = np.sqrt(G*mcum(Rgas/r200,c)*m200/(Rgas))
vcirc_MN  = MN_Vcirc(Rgas, pgas[:,2], Md=MN_Md, Rd=MN_Rd, Zd=MN_Zd)
vcirc_tot = np.sqrt(vcirc_dm**2+vcirc_gas**2+vcirc_MN**2)

vgas[:,0] = vcirc_tot*np.sin(phi)
vgas[:,1] = -vcirc_tot*np.cos(phi)


boxsize = 400.0  #Boxsize for SWIFT

for i in range(3):
    pgas[:,i] += 200

    
with h5py.File(outputfile,'w') as snap:
    header = snap.create_group('Header')
    
    particles = snap.create_group('PartType0')
    particles.create_dataset('Coordinates',data=pgas,dtype=np.float32)
    particles.create_dataset('Velocities',data=vgas,dtype=np.float32)
    particles.create_dataset('InternalEnergy',data=ugas,dtype=np.float32)
    particles.create_dataset('Densities',data=rhogas,dtype=np.float32)
    particles.create_dataset('Masses',data=mgas,dtype=np.float32)
    particles.create_dataset('SmoothingLength',data=np.ones(len(mgas)),dtype=np.float32)
    particles.create_dataset('ParticleIDs',data=np.arange(len(mgas)),dtype='uint')

    header.attrs['NumPart_ThisFile'] = np.array([len(mgas),0,0,0,0,0], dtype='uint')
    header.attrs['NumPart_Total'] = np.array([len(mgas),0,0,0,0,0], dtype='uint')
    header.attrs['MassTable'] = np.array([0,0,0,0,0,0], dtype='double')
    header.attrs['Time'] = 0
    header.attrs['Redshift'] = 0
    header.attrs['Flag_Entropy_ICs'] = 0
    header.attrs['NumFilesPerSnapshot'] = 1
    header.attrs['NumPart_Total_HighWord'] = np.array([0,0,0,0,0,0], dtype='int')
    header.attrs['BoxSize'] = boxsize

    snap.close()







