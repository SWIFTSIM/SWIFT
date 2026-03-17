'''
IC making file written by Filip Husko,
Last edits - Lily Magnus March 2026
'''

#!/usr/bin/env python3
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.interpolate as sci
from tqdm import tqdm
import sys
import time
from shutil import copyfile
from unyt import unyt_array
import unyt
np.random.seed(420)
import swiftsimio.visualisation as swvis
from swiftsimio.visualisation.smoothing_length.generate import generate_smoothing_lengths
import scipy.optimize as scop
from scipy.optimize import root_scalar
from swiftsimio.objects import cosmo_array, cosmo_factor
from swiftsimio.units import cosmo_units
from swiftsimio import Writer
import swiftsimio as sw


# definition of constants
G=4.3e4 #cm^2 s^-2
mu = 0.62
proton_mass = 1.6726e-24 # gram
kB = 1.38e-16 # erg / K
gamma = 5./3.

def V200_func(M200):
    return 311.7*(M200/1e13)**(1/3)

def R200_func(M200):
    return 442.7*(M200/1e13)**(1/3)

def vcirc_dmhalo(r,a,Mtotal):
    return G * Mtotal * r / (r+a)**2

def vcirc_dmhalo_NFW(r,Rs,M200,c):
    return G * Mp_NFW(r, Rs, M200, c) / r

def Mp_NFW(r,Rs,M200, c):
    return M200 / (np.log(1+c) - c/(1+c)) * (np.log(1+r/Rs) - r/(r+Rs))

def Mp_hern(r,a,Mtotal):
    return Mtotal * r**2 / (r+a)**2

def rho_NFW(r, Rs, M200, c):
    return M200/(4*np.pi*Rs**3*(np.log(1+c) - c/(1+c))) / (r/Rs) / (1+r/Rs)**2

def rho_NFW_minus_rhox(r, Rs, M200, c, rhox):
    return rho_NFW(r, Rs, M200, c) - rhox

def average_rho_NFW_minus_rhox(r, Rs, M200, c, rhox):
    return 3*Mp_NFW(r, Rs, M200, c)/(4*np.pi*r**3) - rhox

def rho_hern(r,a,Mtotal):
    return Mtotal/(2*np.pi * a**3) / (r/Rs) / (1+r/Rs)**3

def rho_hern_minus_rhox(r,a, Mtotal, rhox):
    return rho_hern(r,a,Mtotal) - rhox

def average_rho_hern_minus_rhox(r,a, Mtotal, rhox):
    return 3*Mp_hern(r,a, Mtotal)/(4*np.pi*r**3) - rhox

def M_stellar_bulge(coordinates, masses,radius_array):
    r2 = coordinates[:,0]**2 + coordinates[:,1]**2 + coordinates[:,2]**2
    r = r2**.5
    
    mass_encl = np.zeros(len(radius_array))
    for i in range(0,len(radius_array)):
        mask = r<radius_array[i]
        mass_encl[i] = np.sum(masses[mask])

    return mass_encl

def gas_fraction_rho500(m500,alpha,mt,omegab=0.0463,omegam=0.2793):
    return 0.5*omegab/omegam * (1+ np.tanh(np.log10(m500/mt)/alpha))

def mean_density(r_array, rho_avg, rho_500,rho_200,rho_c):
    plt.plot(r_array,rho_avg)
    plt.plot(r_array,np.ones(len(r_array))*rho_500,linestyle="--")
    plt.plot(r_array,np.ones(len(r_array))*rho_200,linestyle="--")
    plt.plot(r_array,np.ones(len(r_array))*rho_c,linestyle="--")
    plt.xlabel("Radius [$\\rm kpc$]")
    plt.ylabel("Average density [$\\rm M_\odot \\rm kpc^{-3}$]")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("Mean_density_guess.png")
    plt.close()

def find_gas_mass(Mdm_encl, Mstar_encl, r_array):
    # Calculate the point at which we are close to r500
    rho_avg = (Mdm_encl + Mstar_encl) / (4*np.pi/3 * r_array**3)
    rho_c = 130
    rho_500 = rho_c*500 # solar mass / kpc^3
    rho_200 = rho_c*200
    mean_density(r_array, rho_avg, rho_500,rho_200,rho_c)

    idx_min_mass_500 = np.argmin(np.abs(rho_avg-rho_500))

    r500_guess = r_array[idx_min_mass_500]
    m500_guess = Mdm_stars_encl[idx_min_mass_500]
    m500_gas_fraction = gas_fraction_rho500(m500_guess, alpha, mt)
    mgas_mass_guess = m500_gas_fraction * m500_guess
    return mgas_mass_guess

def get_gas_T(Tmin, circ_v_tot):
    return Tmin + mu * proton_mass / (kB * gamma) * circ_v_tot

def get_gas_T_center(Tmin, rmin, r,a):
    return np.log10( Tmin / (1 + np.exp((r-a*rmin)/rmin)) )

def get_gas_T_center2(Tmin, rmin, r):
    result = np.zeros(len(r))
    result[r<rmin] = Tmin
    result[r>=rmin] = Tmin*(1- (r[r>=rmin]-rmin) / (5*rmin))
    result[result<0] = 0
    print(result)
    return np.log10(result)

def particle_mass_function(r_array,rtyp):
    mass_factor = np.ones(len(r_array))
    if mass_scheme != "default":
        mask = r_array > rtyp
        mass_factor[mask] = ms_a * (r_array[mask]/rtyp)**ms_n - ms_b
    return mass_factor

def mass_scheme_used(mass_scheme):
    if mass_scheme=="default":
        ms_n = 1
        ms_a = 1
        print("Generate ICs with constant mass")
    elif mass_scheme=="linear":
        ms_n = 1
        ms_a = 1
        print("Generate ICs with linear1")
    elif mass_scheme=="quadratic":
        ms_n = 2
        ms_a = 1
        print("Generate ICs with quadratic")
    elif mass_scheme=="cubic":
        ms_n = 3
        ms_a = 1
        print("Generate ICs with cubic")
    else: 
        raise TypeError("Invalid choice of mass scheme!") 
    
    ms_b = ms_a - 1
    return ms_n, ms_a, ms_b

def plots_temperature_profiles(r, T, circ):
    plt.plot(np.log10(r),np.log10(T))
    plt.xlabel("log radius [kpc]")
    plt.ylabel("log temperature [K]")
    plt.savefig(f"Updated_temperature_profiles_{mass_scheme}.png")
    plt.close()

    plt.plot(np.log10(r_array), circ/(kB * gamma * (T))*mu*proton_mass)
    plt.xlabel("log radius [kpc]")
    plt.ylabel("circular velocity / sound speed")
    plt.ylim(0,1.01)
    plt.savefig(f"Ratio_circ_velocity_sound_speed_{mass_scheme}.png")
    plt.close() 

def calculate_pressure(r, T, circ):
    log_pressure_array = np.zeros(len(r))
    f_circ_v = sci.interp1d(r, circ)
    f_T = sci.interp1d(r, T)
    # use an estimate for the density at Rmin
    rho0 = 10.
    P0 = rho0 * kB * T[0]
    log_pressure_array[0] = np.log10(P0)

    # loop over all radii
    for i in range(1,len(r)):
        r_mean = 10**( (np.log10(r[i]) + np.log10(r[i-1]))/2.)
        delta_logr = np.log10(r[i]) - np.log10(r[i-1])
        log_pressure_array[i] = log_pressure_array[i-1] - f_circ_v(r_mean) * mu * proton_mass / (kB * f_T(r_mean)) * delta_logr

    plt.plot(np.log10(r), log_pressure_array)
    plt.xlabel("log radius [kpc]")
    plt.ylabel("log pressure [Barye]")
    plt.savefig(f"Pressure_profile_guess_{mass_scheme}.png")
    plt.close()

    plt.plot(np.log10(r), log_pressure_array - np.log10(kB) - np.log10(T))
    plt.xlabel("log radius [kpc]")
    plt.ylabel("log density [$\\rm cm^{-3}$]")

    plt.savefig(f"density_profile_guess_{mass_scheme}.png")
    plt.close()
    return log_pressure_array

# Debackere model
alpha = 1.35
mt = 10**13.94
amount_of_sigma_used = 0.

# Parameters of the model
mu_angular = 1.25
lambda_angular = 0.05 # 0.05
s_angular = 1.0
adaptive_mass_radius = float(sys.argv[7]) #kpc
print(f"adaptive mass radius {adaptive_mass_radius} kpc")
minimum_radius = 0.1 # kpc
NFW = False

# get the arguments of the code
file_name = str(sys.argv[1])
file_name_save = str(sys.argv[2])
mass_scheme = str(sys.argv[3])

# parameters of the galaxy:
try:
    galaxy = str(sys.argv[4])
except:
    galaxy = "H13"

try: 
    cc = float(sys.argv[5])
except:
    cc = 0

gas_mass_used = 0.

if galaxy=="H12":
    M200 = 137e10
    bulge_fraction = 0.00
    mBH = 1e6
    # set initial value in model
    Tmin_model = float(10**6)    # K
    minimum_temperature = True 
    Tscale_length = 3.5 # kpc   
    gas_mass_used = 5e10 # this is within R500
    if cc ==0:
        cc = 9
elif galaxy=="H13":
    M200 = 1e13
    bulge_fraction = 0.01
    mBH = 10**8.4
    # set initial value in model
    Tmin_model = float(10**5.75)    # K
    minimum_temperature = False 
    Tscale_length = 30. # kpc   
    gas_mass_used = 1.6218e+11 #-8e12*0.015
    if cc ==0:
        cc = 7.2
elif galaxy=="H13.5":
    M200 = 3e13
    bulge_fraction = 0.005
    mBH = 10**8.6
    # set initial value in model
    Tmin_model = float(10**6.75)    # K
    minimum_temperature = True 
    Tscale_length = 35. # kpc
    gas_mass_used = 1.1445e12 #- 8e13*0.015
    if cc==0:
        cc=6.4
elif galaxy=="H14":
    M200 = 1e14
    bulge_fraction = 0.0025 
    mBH = 10**(8.7)
    # set initial value in model
    Tmin_model = float(10**6.5)    # K
    minimum_temperature = True # True #True 
    Tscale_length = 50. # kpc
    gas_mass_used = 10**(12.637) #- 8e13*0.015
    if cc==0:
        cc=5.6
elif galaxy=="H14.5":
    M200 = 3e14
    bulge_fraction = 0.00125
    mBH = 10**9.0
    # set initial value in model
    Tmin_model = float(10**7.00)    # K
    minimum_temperature = True 
    Tscale_length = 70. # kpc
    gas_mass_used = 3e13 #- 8e13*0.015
    if cc==0:
        cc=4.9
elif galaxy=="H15":
    M200 = 1e15
    bulge_fraction = 0.0006
    mBH = 6.5*10**9
    # set initial value in model
    Tmin_model = float(10**7.4)    # K
    minimum_temperature = True 
    Tscale_length = 100. # kpc
    gas_mass_used = 1.2e14 #- 8e13*0.015
    if cc==0:
        cc=4.0

try:
    Tmin_model = 10**float(sys.argv[6])
except:
    print("We use the default Tmin forthe model")

print(f"Used parameters c={cc:1.4f} Tmin={Tmin_model:1.4e}")

V200 = V200_func(M200)
R200 = R200_func(M200)
Mdm = M200 * (1-bulge_fraction)
Mstar = M200 * bulge_fraction
Rs = R200/cc
bb = 2./cc**2 * (np.log(1+cc) - cc/(1+cc))
a_hernquist = (bb+bb**.5)/(1-bb) * R200
Mdm_hern = Mdm * (R200+a_hernquist)**2 / R200**2

# find the correct coefficients of the mass scheme
ms_n, ms_a, ms_b = mass_scheme_used(mass_scheme)

print(f"mass_scheme: {mass_scheme}")
print(f"ms_n={ms_n}, ms_a={ms_a}, ms_b={ms_b}")
print(f"adaptive_mass_radius={adaptive_mass_radius}")
print(f"particle_mass at R200: {particle_mass_function(np.array([R200]), adaptive_mass_radius)}")
print('-----------')
print(f"Model used:")
print(f"We are running with model: {galaxy}")
print(f"The used mass profile is {mass_scheme}")
print(f"n: {ms_n:1.4f}")
print(f"a: {ms_a:1.4f}")
print(f"b: {ms_b:1.4f}")
print("#"*60)
print("Halo properties")
print(f"Properties:")
print(f"M200: {M200:1.4e} Msun")
print(f"R200: {R200:1.4f} kpc")
print(f"V200: {V200:1.4f} km/s")
print(f"concentration: {cc:1.3f}")
print(f"Rs_NFW: {Rs:1.4f} kpc")
print(f"A_Hernquist: {a_hernquist:1.4f} kpc (Hernquist equivalent scale length)")
print(f"M_dm: {Mdm_hern:1.4e} Msun (total Hernquist mass)")
print("#"*60)
print(f"Stellar properties:")
print(f"Stellar bulge fraction: {bulge_fraction:1.6f}")
print(f"Mstar: {Mstar:1.4e} Msun")
print(f"Reff: Implicit in IC file")
print("#"*60)
print(f"Black hole properties")
print(f"M_BH: {mBH:1.4e} Msun")
print("#"*60)
if minimum_temperature:
    print(f"Constrain on minimum temperature at r=0 are:")
    print(f"Tmin: {Tmin_model:1.4e} K")
    print(f"Tscale length: {Tscale_length:1.4f} kpc")
    print("#"*60)
print(f"We will use the fit in Debackere et al. (2020) to find the gas fraction:")
print(f"alpha: {alpha:1.4f}")
print(f"m_t: {mt:1.4e} Msun")
print(f"Amount of sigma used: {amount_of_sigma_used:1.2f}")
print("#"*60)

# Calculate M500, R500
rho_c = 130
rho_500 = rho_c*500 # solar mass / kpc^3
rho_200 = rho_c*200

results_R500_old = root_scalar(rho_NFW_minus_rhox, args=(Rs, M200, cc, rho_500), method="bisect", bracket=[1, R200] )
results_R1000_old = root_scalar(rho_NFW_minus_rhox, args=(Rs, M200, cc, 2*rho_500), method="bisect", bracket=[1, R200] )

results_R500 = root_scalar(average_rho_NFW_minus_rhox, args=(Rs, M200, cc, rho_500), method="bisect", bracket=[1, R200] )
results_R1000 = root_scalar(average_rho_NFW_minus_rhox, args=(Rs, M200, cc, 2*rho_500), method="bisect", bracket=[1, R200] )

R500_old = results_R500_old.root
R1000_old = results_R1000_old.root
R500 = results_R500.root
R1000 = results_R1000.root
M500 = Mp_NFW(R500, Rs, M200, cc) 
M1000 = Mp_NFW(R1000, Rs, M200, cc) 

print(f"Other halo information (NFW):")
print(f"R500: {R500:1.4f} kpc")
print(f"M500: {M500:1.4e} Msun")
print(f"R500_old: {R500_old:1.4f} kpc")
print(f"R1000: {R1000:1.4f} kpc")
print(f"M1000: {M1000:1.4e} Msun")
print(f"R1000_old: {R1000_old:1.4f} kpc")
print("#"*60)
print(f"Hernquist equivalent halo information:")

results_R500 = root_scalar(average_rho_hern_minus_rhox, args=(a_hernquist, Mdm_hern, rho_500), method="bisect", bracket=[1, R200] )
results_R1000 = root_scalar(average_rho_hern_minus_rhox, args=(a_hernquist, Mdm_hern, 2*rho_500), method="bisect", bracket=[1, R200] )

R500_hern = results_R500.root
R1000_hern = results_R1000.root
M500_hern = Mp_hern(R500, a_hernquist, Mdm_hern) 
M1000_hern = Mp_hern(R1000, a_hernquist, Mdm_hern) 

print(f"R500: {R500_hern:1.4f} kpc (Hernquist value)")
print(f"M500: {M500_hern:1.4e} Msun (Hernquist value)")
print(f"R1000: {R1000_hern:1.4f} kpc (Hernquist value)")
print(f"M1000: {M1000_hern:1.4e} Msun (Hernquist value)")
print("#"*60)

m500_gas_fraction = gas_fraction_rho500(M500, alpha, mt)
m500_gas_fraction_1std = 0.015
hydrostatic_gas_fraction = m500_gas_fraction + amount_of_sigma_used*m500_gas_fraction_1std
mgas_mass_guess = (m500_gas_fraction + amount_of_sigma_used*m500_gas_fraction_1std) * M500

m500_gas_fraction_hern = gas_fraction_rho500(M500_hern, alpha, mt)
hydrostatic_gas_fraction_hern = m500_gas_fraction_hern + amount_of_sigma_used*m500_gas_fraction_1std
mgas_mass_guess_hern = (m500_gas_fraction_hern + amount_of_sigma_used*m500_gas_fraction_1std) * M500_hern

print(f"The hydrostatic gas fraction at R500: {hydrostatic_gas_fraction:1.4f}")
print(f"The hydrostatic gas fraction at R500 (Hernquist case): {hydrostatic_gas_fraction_hern:1.4f}")
print(f"Total gass mass within R500: {mgas_mass_guess:1.4e} Msun")
print(f"Total gass mass within R500 (Hernquist case): {mgas_mass_guess_hern:1.4e} Msun")
if NFW:
    print(f"We will use a NFW profile to determine the properties")
else:
    print(f"We will use a Hernquist profile to determine the properties")
    R500 = R500_hern
    M500 = M500_hern
    m500_gas_fraction = m500_gas_fraction_hern
    mgas_mass_guess = mgas_mass_guess_hern

if gas_mass_used > 0:
    print(f"we overwrite the gas fraction to match the BAHAMAS sample")
    mgas_mass_guess = gas_mass_used
    print(f"The used gas mass is {mgas_mass_guess:1.4e}")

print("#"*60)

rmax = 3*R200
print("Profile properties:")
print(f"Minimum radius used: {minimum_radius:1.4f} kpc")
print(f"End of high resolution region: {adaptive_mass_radius:1.4f} kpc")
print(f"Maximum radius used (3xR200): {rmax:1.4f} kpc")

# copy the existing file to a new file that will get the gas
copyfile(file_name, file_name_save)

# load the stellar data from the IC file that will get the gas 
# particles
f = h5.File(file_name_save, "r+")

Coordinates_gas_disc = f["/PartType0/Coordinates"][:,:]
Masses_gas_disc = f["/PartType0/Masses"][:]
Velocities_gas_disc = f["/PartType0/Velocities"][:,:]
ids_gas_disc = f["/PartType0/ParticleIDs"][:]
smoothing_lengths_gas_disc = f["/PartType0/SmoothingLength"][:]
internal_energies_gas_disc = f["/PartType0/InternalEnergy"][:]

Coordinates = f["/PartType4/Coordinates"][:,:]
Masses = f["/PartType4/Masses"][:] * 1e10
Velocities = f["/PartType4/Velocities"][:,:]
ids = f["/PartType4/ParticleIDs"][:]
# update the IDs by adding 1, such that we can give the BH the 
# lowest ID
offset_stars = int(1e9)
ids_update = ids + offset_stars
ids_data = f["/PartType4/ParticleIDs"]
ids_data[...] = ids_update

# get the boxsize and the coordinates of the star particles
boxsize = f["/Header"].attrs["BoxSize"]

# --- Box size ---                                                                                                                                                                             
boxsize_array = np.array([boxsize[0], boxsize[1], boxsize[2]])
boxsize_writer = boxsize_array * unyt.kpc  # or Mpc if your arrays are in Mpc                                                                                                                                      
Coordinates[:,0] -= boxsize[0]/2.
Coordinates[:,1] -= boxsize[1]/2.
Coordinates[:,2] -= boxsize[2]/2.
rs_stars = np.sqrt(Coordinates[:,0]**2+Coordinates[:,1]**2+Coordinates[:,2]**2)
choice_stars = rs_stars<50

Masses_int=f["/PartType4/Masses"][:][choice_stars]
N_stars = np.size(f["/PartType4/Masses"][:][choice_stars])
Velocities_int=f["/PartType4/Velocities"][:,:][choice_stars]
ids_int = f["/PartType4/ParticleIDs"][:][choice_stars]
coords_stars_int = f["/PartType4/Coordinates"][:,:]
coords_stars_int[:,0] = Coordinates[:,0] + boxsize[0]/2.
coords_stars_int[:,1] = Coordinates[:,1] + boxsize[1]/2.
coords_stars_int[:,2] = Coordinates[:,2] + boxsize[2]/2.
coords_stars_int = coords_stars_int[choice_stars]

del f["/PartType4/Coordinates"]
del f["/PartType4/Velocities"]
del f["/PartType4/Masses"]
del f["/PartType4/ParticleIDs"]
    
g = f["/PartType4"]
ds = g.create_dataset("Coordinates", (N_stars,3), "d")    
ds[()] = coords_stars_int
ds = g.create_dataset("Velocities", (N_stars,3), "f")    
ds[()] = Velocities_int
ds = g.create_dataset("Masses", (N_stars,), "f")    
ds[()] = Masses_int
ds = g.create_dataset("ParticleIDs", (N_stars,), "L")    
ds[()] = ids_int

Coordinates=np.array(coords_stars_int)
del coords_stars_int

N_radii = int(1e4)

r_interest = np.array([R200, 2*R200, 3*R200])
Mass_maximal = particle_mass_function(r_interest,adaptive_mass_radius)

print(f"Resolution high resolution regions: {Masses[0]:1.4e} Msun")
print(f"Resolution of stellar particles: {Masses[0]:1.4e} Msun")
print(f"IDs of the stellar particles in the IC are offset by {offset_stars:1.4e}")
print(f"Particle mass at R200: {Mass_maximal[0]*Masses[0]:1.4e} Msun")
print(f"Particle mass at 2xR200: {Mass_maximal[1]*Masses[0]:1.4e} Msun")
print(f"Maximal particle mass (r=3xR200): {Mass_maximal[2]*Masses[0]:1.4e} Msun")
print(f"Number of points used to determine radial profile: {N_radii:d}")


# create two arrays that go from the minimum radius to the adaptive radius 
# from the adaptive radius to the maximum radius. 
r_inner = np.linspace(minimum_radius, adaptive_mass_radius, int(0.1*N_radii))
r_outer = np.linspace(adaptive_mass_radius, rmax, int(0.9*N_radii)+1, endpoint=True)


print(f"high resolution region: {len(r_inner):d}")
print(f"lower resolution region: {len(r_outer):d}")
print("#"*60)

# merge both arrays
r_array = np.append(r_inner,r_outer[1:])

# calculate the enclosed dark matter mass 
if NFW:
    Mdm_encl = Mp_NFW(r_array, Rs, M200, cc)
else:
    Mdm_encl = Mp_hern(r_array,a_hernquist, Mdm_hern)

# calculate the enclosed stellar mass
Mstar_encl = M_stellar_bulge(Coordinates, Masses_int, r_array)

# calculate the total enclosed mass (stars + DM)
Mdm_stars_encl = Mdm_encl+Mstar_encl

plt.plot(np.log10(r_array), np.log10(Mstar_encl),label="stars")
plt.plot(np.log10(r_array), np.log10(Mdm_encl), label="DM")
plt.xlabel("log r [kpc]")
plt.ylabel("log Mass")
plt.savefig("enclosed_mass_hern.png")
plt.close()

# Calculate the circular velocity assuming gas has no contribution
# stellar part:
circular_velocity = G * Mstar_encl / r_array

# Dark matter part:
if NFW:
    circular_vel_dm = vcirc_dmhalo_NFW(r_array,Rs,M200,cc)
else:
    circular_vel_dm = vcirc_dmhalo(r_array, a_hernquist, Mdm_hern)

# total circular velocity profile
circ_v_tot = circular_velocity + circular_vel_dm

# get the temperature profile
# get the temperature profile
log_T_profile = np.log10(get_gas_T(0, circ_v_tot))
log_T_center = get_gas_T_center(Tmin_model, Tscale_length, r_array,2.5)

# include or not include the minimum temperature at the center
if minimum_temperature:
    T_profile = 10**log_T_profile + 10**log_T_center
else:
    T_profile = 10**log_T_profile

# Make some plots of the temperature profile
plots_temperature_profiles(r_array, T_profile, circ_v_tot)

print(Tmin_model)
print(np.shape(T_profile))
print(np.shape(circ_v_tot))
np.savetxt(f"temp_vs_circ_v_{Tmin_model:1.4e}.txt", (T_profile, circ_v_tot), fmt="%1.6e")

# calculate the pressure from the circular velocity and the temperature
log_pressure_array = calculate_pressure(r_array, T_profile, circ_v_tot)

# calculate the bracketting values of the radius
idx_r = np.argmin(np.abs(r_array-R500))
if r_array[idx_r] > R500:
    R2 = r_array[idx_r]
    R1 = r_array[idx_r-1]
else:
    R2 = r_array[idx_r+1]
    R1 = r_array[idx_r]

delta_R = R2 - R1

# calculate the density profile arbitrary normalization
log_rho = log_pressure_array - np.log10(kB) - np.log10(T_profile)
rho = 10**log_rho

dr = np.diff(r_array)
r_mean = (r_array[1:] + r_array[:-1])/2.
rho_mean = (rho[1:] + rho[:-1])/2.

# calculate the enclosed mass
mass_gas = rho_mean * mu * proton_mass * 4 * np.pi * r_mean**2 * dr * 3.0857e21**3
mass_gas_encl = np.cumsum(mass_gas) / 1.9891e33

if r_array[idx_r] > R500:
    M2 = mass_gas_encl[idx_r]
    M1 = mass_gas_encl[idx_r-1]
else:
    M2 = mass_gas_encl[idx_r+1]
    M1 = mass_gas_encl[idx_r]

M500_estimate = M1 * (R2 - R500)/delta_R + M2 * (R500-R1)/delta_R

correction_factor = mgas_mass_guess / M500_estimate

print(f"Match the correct M500 value and solve the hydrostatic equilibrium ODE")
print(f"Current initial estimate of M500_gas: {M500_estimate:1.4e} Msun")
print(f"Desired M500_gas: {mgas_mass_guess:1.4e} Msun")
print(f"Estimated density at minimum distance: {rho[0]:1.4e} cm^-3")
print(f"Correction factor becomes: {correction_factor:1.4e}")

rho *= correction_factor
print(f"Applied correction factor.")
print(f"corrected density at minimum distance: {rho[0]:1.4e} cm^-3")
print("#"*60)

if NFW:
    np.savetxt("NFW_temp.txt",np.transpose([r_array,T_profile,rho]))
else:
    np.savetxt("hernquist_temp.txt",np.transpose([r_array,T_profile,rho]))

# Calculate the angular momentum stuff convert to cgs
j0_angular = 2**.5 * V200 * R200 * lambda_angular / (- mu_angular * np.log(1.-1./mu_angular) - 1.) * 3.086e26

jmax_angular = j0_angular / (mu_angular - 1.0)
    
print(f"Setting up the angular momentum profile, e.g. velocity profile")
print(f"\lambda: {lambda_angular:1.4f}")
print(f"\mu: {mu_angular:1.4f}")
print(f"j0: {j0_angular/3.086e26:1.4e} km/s kpc")
print(f"jmax: {jmax_angular/3.086e26:1.4e} km/s kpc")


# generate particles for the inner part
dr_inner = r_inner[1]-r_inner[0]

mask_inner = r_array<=adaptive_mass_radius

mass_gas_inner = rho[mask_inner] * mu * proton_mass * 4 * np.pi * r_inner**2 * dr_inner * 3.0857e21**3
mass_gas_encl_inner = np.cumsum(mass_gas_inner) / 1.9891e33

Mres = Masses[0]
N_gas_particles_per_shell = mass_gas_inner / 1.9891e33 / Mres
#print(N_gas_particles_per_shell)
N_particles_inner = int(np.sum(N_gas_particles_per_shell))

print(f"Number of inner gas particles = {N_particles_inner:d}")

# construct interpolation functions that we need:
max_particles_per_shell = np.max(N_gas_particles_per_shell)
f_probability = sci.interp1d(r_inner, N_gas_particles_per_shell, kind="linear")

# functions that are also applicable for all distances:
f_circ = sci.interp1d(r_array,circ_v_tot,kind="linear")

f_Mencl = sci.interp1d(r_array,Mdm_stars_encl,kind="linear")
f_T = sci.interp1d(r_array,T_profile,kind="linear")

# calculate the angular momentum of the gas
j_result_update = j0_angular * (f_Mencl(r_array)/M200) * (1./(mu_angular-f_Mencl(r_array)/M200))
v_found = np.log10(j_result_update/(r_array*3.086e21))
last_v_value = v_found[r_array<R200][-1]
v_found[r_array>R200] = last_v_value  + np.log10(R200/r_array[r_array>R200])
f_rotation_v = sci.interp1d(r_array,10**v_found,kind="cubic")


# start of random process!!!
# only for the inner particles

phi_random_inner = np.random.uniform(0,2*np.pi,size=N_particles_inner)
theta_random_inner = np.arccos(np.random.uniform(-1,1,size=N_particles_inner))

# generate a much larger array to find enough particles to do the try
array_length_inner = N_particles_inner * 50

# get random radii
r_values_guess_inner = np.random.uniform(r_array[0], adaptive_mass_radius, array_length_inner)

# get random y values
y_value_inner = np.random.uniform(0,max_particles_per_shell, array_length_inner)

# get maximum y value expected at radii
y_max_value_inner = f_probability(r_values_guess_inner)

# allowed to use this numer?
mask_inner = y_max_value_inner >= y_value_inner
# make the array of the radii that we can use according to rejection 
# sampling 
r_values_using_inner = r_values_guess_inner[mask_inner][:N_particles_inner]

inner_particle_mass = np.ones(N_particles_inner) * Mres


accepted_inner = np.sum(mask_inner)

print("Accepted inner:", accepted_inner)
print("Required inner:", N_particles_inner)

print("Length after slicing:", len(r_values_using_inner))

# prepare to generate the random numbers in the outer part

dr_outer = r_outer[1]-r_outer[0]

mask_outer = r_array>=adaptive_mass_radius

mass_gas_outer = rho[mask_outer] * mu * proton_mass * 4 * np.pi * r_outer**2 * dr_outer * 3.0857e21**3
mass_gas_encl_outer = np.cumsum(mass_gas_outer) / 1.9891e33

# calculate the particle mass at different radii
particle_mass = particle_mass_function(r_outer,adaptive_mass_radius)

N_gas_particles_per_shell_outer = mass_gas_outer / 1.9891e33 / (Mres*particle_mass)
N_particles_outer = int(np.sum(N_gas_particles_per_shell_outer))

print(f"Number of outer gas particles = {N_particles_outer:d}")


print(f"particle_mass min: {particle_mass.min():.4f}, max: {particle_mass.max():.4f}")
print(f"N_particles_outer: {N_particles_outer}")
print(f"N_particles_inner: {N_particles_inner}")

# construct interpolation functions that we need:
max_particles_per_shell = np.max(N_gas_particles_per_shell_outer)
f_probability2 = sci.interp1d(r_outer, N_gas_particles_per_shell_outer, kind="linear")

phi_random_outer = np.random.uniform(0,2*np.pi,size=N_particles_outer)
theta_random_outer = np.arccos(np.random.uniform(-1,1,size=N_particles_outer))

# generate a much larger array to find enough particles to do the try
array_length_outer = N_particles_outer * 70

# get random radii
r_values_guess_outer = np.random.uniform(adaptive_mass_radius, r_array[-1], array_length_outer)

# get random y values
y_value_outer = np.random.uniform(0,max_particles_per_shell, array_length_outer)

# get maximum y value expected at radii
y_max_value_outer = f_probability2(r_values_guess_outer)

# allowed to use this numer?
mask_outer = y_max_value_outer >= y_value_outer
# make the array of the radii that we can use according to rejection 
# sampling 
r_values_using_outer = r_values_guess_outer[mask_outer][:N_particles_outer]


outer_particle_mass = particle_mass_function(r_values_using_outer, adaptive_mass_radius)*Mres

final_radius = np.append(r_values_using_inner, r_values_using_outer)
final_mass_particles = np.append(inner_particle_mass, outer_particle_mass)
final_phi = np.append(phi_random_inner, phi_random_outer)
final_theta = np.append(theta_random_inner, theta_random_outer)

print(f"Realized gas mass within R500: {np.sum(final_mass_particles[final_radius<R500]):1.4e} Msun")

# get the particle temperature
T_particles = f_T(final_radius)

R_using = final_radius * np.sin(final_theta)

v_rot = f_rotation_v(np.maximum(final_radius,r_array[0]))

# get the internal energy
internal_energy = kB * T_particles / (0.62 * proton_mass * (gamma-1))

# get the position + velocity
x_gas = final_radius * np.sin(final_theta) * np.cos(final_phi)
y_gas = final_radius * np.sin(final_theta) * np.sin(final_phi)
z_gas = final_radius * np.cos(final_theta)

#box_size = 2200


x_gas += boxsize[0] / 2. #boxsize[0]/2.
y_gas += boxsize[1] / 2. #boxsize[1]/2.
z_gas += boxsize[2] / 2. #boxsize[2]/2.

#v_rot = 0.

velocity_x = - v_rot * np.sin(final_phi)
velocity_y = v_rot * np.cos(final_phi)
velocity_z = 0* v_rot

coords = np.zeros((len(final_radius),3))
coords[:,0] = x_gas
coords[:,1] = y_gas
coords[:,2] = z_gas
selection_gas = (coords[:,0]>0) & (coords[:,0]<boxsize[0]) & (coords[:,1]>0) & (coords[:,1]<boxsize[1]) & (coords[:,2]>0) & (coords[:,2]<boxsize[2])

coords = unyt_array(coords, "kpc")
boxsize = unyt_array(np.array([boxsize[0],boxsize[1],boxsize[2]]), "kpc")
#print(coords)
coords = cosmo_array(coords, unyt.kpc, comoving=True, scale_factor = 1, scale_exponent=1)

number_of_gas_particles = len(x_gas)
print(f"Number of gas particles = {number_of_gas_particles:d}")
print("Calculate the smoothing lengths:")
#smoothing_lengths = np.ones(len(x_gas))

smoothing_lengths = sw.visualisation.smoothing_length.generate.generate_smoothing_lengths(
        coords[selection_gas], 
        boxsize=boxsize, 
        kernel_gamma=1.936492, 
        neighbours=57,
        dimension=3,
        speedup_fac = 2)

print("Finished calculating the smoothing lengths")
#print(smoothing_lengths)
print(np.max(smoothing_lengths))

# store everything in the IC file:
N = len(x_gas[selection_gas])+len(Masses_gas_disc)

coord_dtype = f["PartType0/Coordinates"].dtype
vels_dtype = f["PartType0/Velocities"].dtype
m_dtype = f["PartType0/Masses"].dtype
h_dtype = f["PartType0/SmoothingLength"].dtype
u_dtype = f["PartType0/InternalEnergy"].dtype
id_dtype = f["PartType0/ParticleIDs"].dtype
 
del f["/PartType0"]
grp = f.create_group("/PartType0")

coords = np.zeros((N,3))
coords[:,0] = np.concatenate((x_gas[selection_gas],Coordinates_gas_disc[:,0]),axis=None)
coords[:,1] = np.concatenate((y_gas[selection_gas],Coordinates_gas_disc[:,1]),axis=None)
coords[:,2] = np.concatenate((z_gas[selection_gas],Coordinates_gas_disc[:,2]),axis=None)

ds = grp.create_dataset("Coordinates", (N,3), coord_dtype)
ds[()]=coords

v = np.zeros((N,3))
v[:,0] = np.concatenate((velocity_x[selection_gas]*1e-5,Velocities_gas_disc[:,0]),axis=None)
v[:,1] = np.concatenate((velocity_y[selection_gas]*1e-5,Velocities_gas_disc[:,1]),axis=None)
v[:,2] = np.concatenate((velocity_z[selection_gas]*1e-5,Velocities_gas_disc[:,2]),axis=None)

ds = grp.create_dataset("Velocities", (N,3), vels_dtype)
ds[()]=v

ds = grp.create_dataset("Masses", (N,), m_dtype)
ds[()]=np.concatenate((final_mass_particles[selection_gas]/1e10,Masses_gas_disc),axis=None)

ds = grp.create_dataset("SmoothingLength", (N,), h_dtype)
ds[()]=np.concatenate((np.array(smoothing_lengths.value),smoothing_lengths_gas_disc),axis=None)

ds = grp.create_dataset("InternalEnergy", (N,), u_dtype)
ds[()]=np.concatenate((internal_energy[selection_gas]*1e-10,internal_energies_gas_disc),axis=None)
N_particles_outer
N_particles_inner


ids_inner = int(1e12) + np.arange(1,N_particles_inner+1,1)
ids_outer = int(1e14) + np.arange(1,N_particles_outer+1,1)
new_ids = np.append(ids_inner, ids_outer)
ds = grp.create_dataset("ParticleIDs", (N,), id_dtype)
#ids_special = final_radius[selection_gas] > 1000
new_ids = np.concatenate((new_ids[selection_gas],ids_gas_disc+np.amax(new_ids)),axis=None)
#new_ids[ids_special] = 0
ds[()]=new_ids


if '/PartType5' not in f:
    print("BH data not available, generating now...")
    # Add the black hole
    bh_coordinates = boxsize/2.+0.001*unyt.kpc 
    spin=0.5 
    del f["/PartType5"]
    grp = f.create_group("/PartType5")
    ds = grp.create_dataset("Coordinates", (1,3), "d")
    
    #bh_coordinates[0] += 1 *0.3 * unyt.kpc
    bh_coordinates = np.array([boxsize[0]/2., boxsize[1]/2., boxsize[2]/2.]) + 0.001
    ds[()] = bh_coordinates 
    
    v = np.zeros((1,3))
    ds = grp.create_dataset("Velocities", (1,3), "f")
    #v[0,1] = 250
    ds[()]=v
    
    ds = grp.create_dataset("SmoothingLength", (1,), "f")
    ds[()]=np.ones(1)
    
    ds = grp.create_dataset("ParticleIDs", (1,), "L")
    ds[()]=np.ones(1)*int(1)
    
    ds = grp.create_dataset("Masses", (1,), "f")
    ds[()]=np.ones(1)*mBH/1e10
    
    ds = grp.create_dataset("SubgridMasses", (1,), "f")
    ds[()]=np.ones(1)*mBH/1e10
    
    ds = grp.create_dataset("EnergyReservoir", (1,), "f")
    ds[()]=np.zeros(1)
    
    spin=0.5
    ds = grp.create_dataset("Spins", (1,), "f")
    ds[()]=np.ones(1)*spin
    
    ds = grp.create_dataset("AngularMomentumDirections", (1,3), "f")
    ds[()] = np.array([[0],[0],[1]]).T

    

# Also update the header
header = f["Header"]
ThisFile = header.attrs["NumPart_ThisFile"]
new_values = np.zeros(6)
new_values[4] = N_stars
new_values[0] = N
new_values[5] = 1
Total = header.attrs["NumPart_Total"]
new_values2 = np.zeros(6)
new_values2[4] = N_stars
new_values2[0] = N
new_values2[5] = 1
boxsize = np.array([boxsize[0],boxsize[1],boxsize[2]]) #header.attrs["BoxSize"]

del f["Header"]
grp = f.create_group("/Header")
grp.attrs["BoxSize"] = boxsize
grp.attrs["NumPart_Total"] = new_values2
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = new_values
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

f.close()
