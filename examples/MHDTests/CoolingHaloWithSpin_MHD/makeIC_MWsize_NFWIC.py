###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durham.ac.uk)
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

# Halo parameters similar to R. Pakmor, V. Springel, 1212.1452

import h5py
import sys
import numpy as np
import math
import random

# Generates N particles in a spherically symmetric distribution with density profile ~r^(-2)
# usage: python3 makeIC.py 1000: generate 1000 particles

# Some constants cgs
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
CONST_MU0_CGS = 4 * np.pi * 1e-2
MSOL_IN_CGS = 1.9891e33 # Solar mass 
kb_cgs = 1.38e-16 # boltzmann constant
m_H_cgs = 1.68e-24 # atomic hydrogen mass 
# First set unit velocity and then the circular velocity parameter for the isothermal potential
const_unit_velocity_in_cgs = 1.0e5  # kms^-1

# Cosmological parameters
OMEGA = 0.3  # Cosmological matter fraction at z = 0
h = 0.67777  # hubble parameter
# Find H_0, the inverse Hubble time, in cgs
H_0_cgs = 100.0 * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# Gas parameters
gamma = 5.0 / 3.0
T0_cgs = 1e5 # gas temperature on the edge of the box (if we want to set this manually)
nH_max_cgs = 1e0 # maximum hydrogen number density

# DM halo parameters
spin_lambda = 0.05  # spin parameter
f_b = 0.17  # baryon fraction
c_200 = 7.2 # concentration parameter

# Set the magnitude of the uniform seed magnetic field
B0_Gaussian_Units = 1e-9 #1e-6  # 1 micro Gauss
B0_cgs = np.sqrt(CONST_MU0_CGS / (4.0 * np.pi)) * B0_Gaussian_Units

# SPH
eta = 1.3663 # kernel smoothing

# From this we can find the virial radius, the radius within which the average density of the halo is
# 200. * the mean matter density

# Set M200 and get R200 and V200
M_200_cgs = 1e12 * MSOL_IN_CGS 
rhoc_cgs = 3*H_0_cgs**2/(8*np.pi*CONST_G_CGS)
r_200_cgs = (3*M_200_cgs/(4*np.pi*rhoc_cgs))**(1/3)
v_200_cgs = np.sqrt(CONST_G_CGS*M_200_cgs/r_200_cgs)
v_200 = v_200_cgs / const_unit_velocity_in_cgs 
T_200_cgs = m_H_cgs*v_200_cgs**2/(2*kb_cgs)

# Gas parameters
gamma = 5.0 / 3.0
T0_cgs = T_200_cgs #1e5 # gas temperature on the edge of the box (if we want to set this manually)
nH_max_cgs = 1e0 # maximum hydrogen number density
print("T_200 = %E" % T_200_cgs)

# Now set the unit length and mass
const_unit_mass_in_cgs = M_200_cgs
const_unit_length_in_cgs = r_200_cgs
print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs)
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)

# Now set the magnetic field unit
const_unit_magnetic_field_in_cgs = 1e-7  # 1muG

# Derived quantities
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_unit_current_in_cgs = const_unit_mass_in_cgs / (
    const_unit_magnetic_field_in_cgs * const_unit_time_in_cgs ** 2
)
print("UnitTime_in_cgs:     ", const_unit_time_in_cgs)
print("UnitCurrent_in_cgs:  ", const_unit_current_in_cgs)
const_G = (
    CONST_G_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)
print("G=", const_G)

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 4.0
G = const_G
N = int(sys.argv[1])  # Number of particles

# Create the file
filename = "CoolingHalo.hdf5"
file = h5py.File(filename, "w")

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = (
    const_unit_length_in_cgs / const_unit_velocity_in_cgs
)
grp.attrs["Unit current in cgs (U_I)"] = const_unit_current_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# NFW-like gas density profile
def rho_r(r_value,f_b,M_200_cgs,r_200_cgs,c_200):
    rho_0 = M_200_cgs/(np.log(1+c_200)-c_200/(1+c_200))/(4*np.pi*r_200_cgs**3/c_200**3)
    result_cgs = rho_0*f_b/(c_200*r_value*(1+c_200*r_value)**2)
    # Apply density cut
    rho_max_cgs = nH_max_cgs*m_H_cgs
    result_cgs = np.array(result_cgs)
    result_cgs[result_cgs>rho_max_cgs]=rho_max_cgs
    return result_cgs 

# NFW-like gas mass inside a sphere with radius R 
def Mgas_r(r_value,f_b,M_200_cgs,r_200_cgs,c_200):
    M_0 = M_200_cgs/(np.log(1+c_200)-c_200/(1+c_200))
    return M_0 * f_b * (np.log(1+c_200*r_value)-r_value/(1+c_200*r_value)) 

# NFW Gravitational acceleration
def a_NFW(r_value, M_200_cgs,r_200_cgs, c_200):
    a_pref = CONST_G_CGS*M_200_cgs/(np.log(1+c_200)-c_200/(1+c_200))/r_200_cgs**2
    return a_pref*((r_value/(r_value+1/c_200))-np.log(1+c_200*r_value))/r_value**2

# Integrate rho_gas*a_NFW
def integrate(r_min, r_max, f_b, M_200_cgs, r_200_cgs, c_200, Nsteps = 1000):
    # Perform the integration
    r_range = np.linspace(r_min, r_max, Nsteps)
    dr = np.abs((r_max-r_min)/Nsteps)
    integrands = rho_r(r_range, f_b, M_200_cgs, r_200_cgs, c_200) * a_NFW(r_range, M_200_cgs,r_200_cgs, c_200)
    result_cgs = np.sum(integrands*dr)*r_200_cgs
    return result_cgs

# NFW-like gas hydrostatic equilibrium internal energy profile
def u_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200):
    result_cgs = (P0_cgs-integrate(r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200))/(gamma-1)/rho_r(r_value,f_b,M_200_cgs,r_200_cgs,c_200)
    return result_cgs


# NFW-like gas hydrostatic equilibrium temperature profile
def T_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200):    
    result_cgs = u_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200) * (gamma-1) * m_H_cgs/kb_cgs
    return result_cgs

# Masses
gas_mass = (
    Mgas_r(boxSize*np.sqrt(3)/2,f_b,M_200_cgs,r_200_cgs,c_200)/M_200_cgs
)  # get total gas mass within a sphere of radius sqrt(3)/2*Lbox
gas_particle_mass = gas_mass / float(N)
print('Gas particle mass is %E ' % (gas_particle_mass*1e12))

# Unnormalized mass shell distribution
r = np.logspace(np.log10(1e-6*boxSize),np.log10(boxSize * np.sqrt(3.0) / 2.0),round(10*N**(1/3)))
#r = np.linspace(1e-6*boxSize,boxSize * np.sqrt(3.0) / 2.0,round(10*N**(1/3)))
f_r_distr_func = 4*np.pi*r**2 * rho_r(r, f_b, M_200_cgs, r_200_cgs, c_200)

# Plot
#rho0_cgs = rho_r(boxSize*np.sqrt(3)/2,f_b,M_200_cgs,r_200_cgs,c_200) #gas density on the edge
#P0_cgs = rho0_cgs*kb_cgs*T0_cgs/m_H_cgs # gas pressure on the edge of the box
#densities_vs_r = rho_r(r,f_b,M_200_cgs,r_200_cgs,c_200)
#u_plot = [u_vs_r(0, r[i], boxSize*np.sqrt(3)/2, f_b, M_200_cgs, r_200_cgs, c_200)/const_unit_velocity_in_cgs**2 for i in range(len(r))] # gas particle internal energies
#u_plot = [u_vs_r(P0_cgs, r[i], boxSize*np.sqrt(3)/2, f_b, M_200_cgs, r_200_cgs, c_200)/const_unit_velocity_in_cgs**2 for i in range(len(r))] # gas particle internal energies
#T_plot = [T_vs_r(P0_cgs, r[i], boxSize*np.sqrt(3)/2, f_b, M_200_cgs, r_200_cgs, c_200) for i in range(len(r))]

# Mass shell distribution function
dM_dr = 4*np.pi*r**2 * rho_r(r, f_b, M_200_cgs, r_200_cgs, c_200)

# Normalize f(r) to get a probability density function (PDF)
pdf = dM_dr / np.trapz(dM_dr, r)  # Normalize with trapezoidal integration
# Compute CDF
cdf = np.cumsum(pdf * np.diff(np.concatenate(([0], r))))  # Cumulative sum approximation
# set seed for random number
np.random.seed(1234)
# Generate N uniform random numbers in [0,1]
uniform_samples = np.random.rand(N)
# Use inverse transform sampling (interpolate CDF to find r values)
sampled_r_values = np.interp(uniform_samples, cdf, r)

# fibonacci-based grid spherical distribution
def fibonacci_sphere(Np):
    """Returns `Np` points distributed on a unit sphere using a Fibonacci lattice."""
    indicesp = np.arange(0, Np, dtype=float) + 0.5
    phip = (1 + np.sqrt(5)) / 2  # Golden ratio
    thetap = 2 * np.pi * indicesp / phip
    zp = 1 - (2 * indicesp / Np)
    mask_badzp = np.abs(zp)>1
    zp[mask_badzp] = np.sign(zp[mask_badzp])*1
    radiusp = np.sqrt(1 - zp**2)
    xp, yp = radiusp * np.cos(thetap), radiusp * np.sin(thetap)

    coordsp = np.vstack((xp, yp, zp)).T
    return coordsp

def generate_coordinates(r_valuesp, N_vs_r_dr):
    """Generates 3D coordinates based on a given radial distribution and particle counts."""
    all_coords = []
    
    for i, (rp, Np) in enumerate(zip(r_valuesp, N_vs_r_dr)):
        if Np > 2:  # Avoid empty shells
            sphere_points = fibonacci_sphere(Np)  # Generate points on unit sphere
            scaled_points = sphere_points * rp  # Scale by radius
            all_coords.append(scaled_points)
 
    return np.vstack(all_coords)  # Stack all coordinates into a single array




# random-based spherical distribution
def generate_coordinates_rand(sampled_r_values):
    radius = (
          sampled_r_values # np.random.rand(N)
    )  # the diagonal extent of the cube
    ctheta = -1.0 + 2 * np.random.rand(N)
    stheta = np.sqrt(1.0 - ctheta ** 2)
    phi = 2 * math.pi * np.random.rand(N)
    coords = np.zeros((N, 3))
    coords[:, 0] = radius * stheta * np.cos(phi)
    coords[:, 1] = radius * stheta * np.sin(phi)
    coords[:, 2] = radius * ctheta

# Example usage:
# intervals = np.array([[1, 3], [2, 6], [5, 10]])
# x_min = 1
# x_max = 10
# print(minimal_covering_intervals(intervals, x_min, x_max))


def minimal_covering_intervals(intervals, x_min, x_max):

    # Sort intervals based on their starting points (and ending points in case of ties)
    
    i = 0
    indeces = np.array([0])
    counter = 0

    while i<len(intervals)-1:
        distances_to_left = intervals[:,0]-intervals[i,1] # select ith interval and find distances to left edges
        indeces_non_intersect = np.argwhere(distances_to_left>0) # find intervals that don't intersect
        if len(indeces_non_intersect)>0: # if there are still non-intersecting intervals to the left
            ilast = np.min(indeces_non_intersect)-1 # find index of last intersecting interval
            #print('way 1.1 i=',i)
        else:
            #print('way 1.2 i=',i)
            break
        
        if i<ilast: 
            indeces = np.append(indeces,ilast)
            i=ilast
            print('way 2.1 i=',i)
        else:
            print(ilast)
            i=ilast+1
            print('way 2.2 i=',i)
        
        # Loop safety
        counter+=1
        if counter>len(intervals):
            print('program loop')
            break

    return indeces

# Generate sphere
d_bin = (gas_particle_mass*M_200_cgs/rho_r(r, f_b, M_200_cgs, r_200_cgs, c_200))**(1/3)/const_unit_length_in_cgs
dN_bin = np.round(4*np.pi*r**2/(d_bin)**2,0)
bins = np.vstack([r-d_bin/2,r+d_bin/2]).T
indeces_selected = minimal_covering_intervals(bins, d_bin[0]/2, r[-1])

#dN_bin *= N/np.sum(dN_bin)
#print(np.sum(dN_bin[indeces_selected]),N)

import matplotlib.pyplot as plt
plt.scatter(r[indeces_selected],dN_bin[indeces_selected])
#plt.yscale('log')
plt.savefig('test_N_vs_r.png')

coords = generate_coordinates(r[indeces_selected],dN_bin[indeces_selected])

# Masses
#gas_mass = (
#    Mgas_r(boxSize*np.sqrt(3)/2,f_b,M_200_cgs,r_200_cgs,c_200)/M_200_cgs
#)  # get total gas mass within a sphere of radius sqrt(3)/2*Lbox
#gas_particle_mass = gas_mass / float(N)
print('Gas particle mass is %E ' % (gas_particle_mass*1e12))

# shift to centre of box
coords += np.full(coords.shape, boxSize / 2.0)
print("x range = (%f,%f)" % (np.min(coords[:, 0]), np.max(coords[:, 0])))
print("y range = (%f,%f)" % (np.min(coords[:, 1]), np.max(coords[:, 1])))
print("z range = (%f,%f)" % (np.min(coords[:, 2]), np.max(coords[:, 2])))

print(np.mean(coords[:, 0]))
print(np.mean(coords[:, 1]))
print(np.mean(coords[:, 2]))

# now find the particles which are within the box

x_coords = coords[:, 0]
y_coords = coords[:, 1]
z_coords = coords[:, 2]

ind = np.where(x_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(x_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(y_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(y_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(z_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(z_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

# count number of particles

N = x_coords.size

print("Number of particles in the box = ", N)

# make the coords and radius arrays again
coords = np.zeros((N, 3))
coords[:, 0] = x_coords
coords[:, 1] = y_coords
coords[:, 2] = z_coords


# Save particle arrangement to .vtk file to view
import vtk
vtk_points = vtk.vtkPoints()
for x, y, z in coords:
    vtk_points.InsertNextPoint(x, y, z)

poly_data = vtk.vtkPolyData()
poly_data.SetPoints(vtk_points)

writer = vtk.vtkPolyDataWriter()
writer.SetFileName("output.vtk")
writer.SetInputData(poly_data)
writer.Write() # you can use paraview to watch how particles are arranged


radius = np.sqrt(
    (coords[:, 0] - boxSize / 2.0) ** 2
    + (coords[:, 1] - boxSize / 2.0) ** 2
    + (coords[:, 2] - boxSize / 2.0) ** 2
)

# now give particle's velocities and magnetic fields
v = np.zeros((N, 3))
B = np.zeros((N, 3))

# first work out total angular momentum of the halo within the virial radius
# we work in units where r_vir = 1 and M_vir = 1
Total_E = v_200 ** 2 / 2.0
J = spin_lambda * const_G / np.sqrt(Total_E)
print("J =", J)
# all particles within the virial radius have omega parallel to the z-axis, magnitude
# is proportional to 1 over the radius
omega = np.zeros((N, 3))
for i in range(N):
    omega[i, 2] = 3.0 * J / radius[i]
    v[i, :] = np.cross(omega[i, :], (coords[i, :] - boxSize / 2.0))
    #B[i, 2] = B0_cgs / const_unit_magnetic_field_in_cgs
    B[i, 0] = B0_cgs / const_unit_magnetic_field_in_cgs

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Particle group
grp = file.create_group("/PartType0")

ds = grp.create_dataset("Coordinates", (N, 3), "d")
ds[()] = coords
coords = np.zeros(1)

ds = grp.create_dataset("Velocities", (N, 3), "f")
ds[()] = v
v = np.zeros(1)

# Magnetic field
ds = grp.create_dataset("MagneticFluxDensities", (N, 3), "f")
ds[()] = B
B = np.zeros(1)

# All particles of equal mass
m = np.full((N,), gas_particle_mass)
ds = grp.create_dataset("Masses", (N,), "f")
ds[()] = m
m = np.zeros(1)

# Smoothing lengths
l = (4.0 * np.pi * radius ** 2 / N) ** (
    1.0 / 3.0
)  # local mean inter-particle separation
h = np.full((N,), eta * l)
ds = grp.create_dataset("SmoothingLength", (N,), "f")
ds[()] = h
h = np.zeros(1)

# Internal energies
rho0_cgs = rho_r(boxSize*np.sqrt(3)/2,f_b,M_200_cgs,r_200_cgs,c_200) #gas density on the edge
P0_cgs = rho0_cgs*kb_cgs*T0_cgs/m_H_cgs # gas pressure on the edge of the box
u = [u_vs_r(P0_cgs, radius[i], boxSize*np.sqrt(3)/2, f_b, M_200_cgs, r_200_cgs, c_200)/const_unit_velocity_in_cgs**2 for i in range(N)] # gas particle internal energies
#u_vs_r(radius)
u = np.full((N,), u)
ds = grp.create_dataset("InternalEnergy", (N,), "f")
ds[()] = u
u = np.zeros(1)

# Particle IDs
ids = 1 + np.linspace(0, N, N, endpoint=False)
ds = grp.create_dataset("ParticleIDs", (N,), "L")
ds[()] = ids

file.close()
