import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas

# Some constants cgs
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
MSOL_IN_CGS = 1.9891e33 # Solar mass 
kb_cgs = 1.38e-16 # boltzmann constant
m_H_cgs = 1.68e-24 # atomic hydrogen mass 
# First set unit velocity and then the circular velocity parameter for the isothermal potential
const_unit_velocity_in_cgs = 1.0e5  # kms^-1

h = 0.681 #0.67777  # hubble parameter
H_0_cgs = 100.0 * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# From this we can find the virial radius, the radius within which the average density of the halo is
# 200. * the mean matter density

# Set M200 and get R200 and V200
f_b = 0.17
c_200 = 7.2
nH_max_cgs = 1e2
M_200_cgs = 1e12 * MSOL_IN_CGS 
rhoc_cgs = 3*H_0_cgs**2/(8*np.pi*CONST_G_CGS)
r_200_cgs = (3*M_200_cgs/(4*np.pi*rhoc_cgs*200))**(1/3)
v_200_cgs = np.sqrt(CONST_G_CGS*M_200_cgs/r_200_cgs)
v_200 = v_200_cgs / const_unit_velocity_in_cgs 
T_200_cgs = m_H_cgs*v_200_cgs**2/(2*kb_cgs)

# load reference
with_reference = False

if with_reference:
    import pandas as pd
    rho_vs_x_data = pd.read_csv('./1Dreference/rho_vs_x.csv',header=None)
    Bx_vs_x_data = pd.read_csv('./1Dreference/Bx_vs_x.csv',header=None)
    By_vs_x_data = pd.read_csv('./1Dreference/By_vs_x.csv',header=None)
    Berr_vs_x_data = pd.read_csv('./1Dreference/Berr_vs_x.csv',header=None)

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# Where we want to slice the slice
y0 = 0.5 #0.5

# Load snapshot
filename = args.input
data = load(filename)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]
n_gas = data.metadata.n_gas

# Retrieve particle attributes of interest
rho = data.gas.densities
rho.convert_to_units(unyt.g*unyt.cm**(-3))

P = data.gas.pressures
T = data.gas.temperatures

B = data.gas.magnetic_flux_densities

Bx, By = B[:, 0], B[:, 1]

normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)

divB = data.gas.magnetic_divergences

h = data.gas.smoothing_lengths

errB = np.log10(h * abs(divB) / normB)

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_densities = data.gas.masses * rho

data.gas.mass_weighted_pressures = data.gas.masses * P

data.gas.mass_weighted_Bx = data.gas.masses * Bx

data.gas.mass_weighted_By = data.gas.masses * By

data.gas.mass_weighted_errB = data.gas.masses * errB

data.gas.mass_weighted_T = data.gas.masses * T

common_arguments = dict(
    data=data, z_slice=0.5 * data.metadata.boxsize[2], resolution=512, parallel=True
)

mass_map = slice_gas(**common_arguments, project="masses")

mass_weighted_density_map = slice_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_pressure_map = slice_gas(
    **common_arguments, project="mass_weighted_pressures"
)

mass_weighted_Bx_map = slice_gas(**common_arguments, project="mass_weighted_Bx")

mass_weighted_By_map = slice_gas(**common_arguments, project="mass_weighted_By")

mass_weighted_errB_map = slice_gas(**common_arguments, project="mass_weighted_errB")
mass_weighted_T_map = slice_gas(**common_arguments, project="mass_weighted_T")

# Take out mass dependence
density_map = mass_weighted_density_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map
Bx_map = mass_weighted_Bx_map / mass_map
By_map = mass_weighted_By_map / mass_map
errB_map = mass_weighted_errB_map / mass_map
T_map = mass_weighted_T_map / mass_map

map_pixel_length = len(mass_map)

n_H = density_map/m_H_cgs

x = np.linspace(0.0, 1.0, map_pixel_length)
slice_ind = int(np.floor(y0 * map_pixel_length))


# Plot maps
plt.rcParams.update({"font.size": 16})

nx = 1
ny = 3
fig, axs = plt.subplots(ny, nx, figsize=((10*nx, 5*ny)), sharex=True)
fig.subplots_adjust(hspace=0.1)

Lbox_kpc = data.metadata.boxsize.to(PARSEC_IN_CGS*1e3*unyt.cm).value[0]
#locs = [map_pixel_length / 4, map_pixel_length / 2, 3 * map_pixel_length / 4]
#labels = [-Lbox_kpc / 4, 0, Lbox_kpc / 4]
x = np.linspace(-Lbox_kpc/2, Lbox_kpc/2, map_pixel_length)
slice_ind = int(np.floor(y0 * map_pixel_length))

axs[0].plot(x, n_H[:, slice_ind], "k-",color='black')
axs[0].set_yscale('log')
axs[0].set_yticks(np.logspace(-7, 2, 10))
axs[0].set_ylim(1e-7, 1e2)


axs[1].plot(x, T_map[:, slice_ind], "k-",color='black')
axs[1].set_yscale('log')

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
def integrate(r_min, r_max, f_b, M_200_cgs, r_200_cgs, c_200, Nsteps = 10000):
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

r_200_kpc = r_200_cgs / (PARSEC_IN_CGS*1e3)
n_H_analytic = rho_r(np.abs(x)/r_200_kpc,f_b,M_200_cgs,r_200_cgs,c_200)/m_H_cgs

# Gas parameters
gamma = 5.0 / 3.0
T0_cgs = T_200_cgs #1e5 # gas temperature on the edge of the box (if we want to set this manually)

rho0_cgs = rho_r(Lbox_kpc/r_200_kpc*np.sqrt(3)/2,f_b,M_200_cgs,r_200_cgs,c_200) #gas density on the edge
P0_cgs = rho0_cgs*kb_cgs*T0_cgs/m_H_cgs # gas pressure on the edge of the box

T_analytic = [T_vs_r(P0_cgs, np.abs(x[i]/r_200_kpc), Lbox_kpc/r_200_kpc*np.sqrt(3)/2, f_b, M_200_cgs, r_200_cgs, c_200) for i in range(map_pixel_length)] # gas particle internal energies


locs = [map_pixel_length / 4, map_pixel_length / 2, 3 * map_pixel_length / 4]
labels = [-Lbox_kpc / 2, 0, Lbox_kpc / 2]

axs[0].set_xlabel(r"$x$ [kpc]")
axs[0].set_xticks(locs, labels)
axs[0].set_xlim(0, map_pixel_length-1)
axs[0].plot(x, n_H_analytic, "k-",color='red')
axs[0].set_ylabel(r"$n_H(x,y_0)$ $[cm^{-3}]$")

axs[1].set_xlabel(r"$x$ [kpc]")
axs[1].set_xticks(locs, labels)
axs[1].set_xlim(0, map_pixel_length-1)
axs[1].plot(x, T_analytic, "k-",color='red')
axs[1].set_ylabel(r"$T(x,y_0)$ $[K]$")
axs[1].set_yscale('log')
axs[1].set_yticks(np.logspace(4, 8, 5))
axs[1].set_ylim(1e4, 1e6)

#axs[1].plot(x, pressure_map[:, slice_ind], "k-",color='black')
#axs[1].set_yticks(np.arange(-0.5, 0.55, 0.25))
#axs[1].set_ylabel(r"$B_x(x,y_0)$")
#axs[1].set_ylim(-0.5, 0.5)

#axs[2].plot(x, By_map[:, slice_ind], "k-", color='black',label='SWIFT')
#axs[2].set_yticks(np.arange(-0.5, 0.55, 0.25))
#axs[2].set_ylabel(r"$B_y(x,y_0)$")
#axs[2].set_ylim(-0.5, 0.5)

#axs[3].plot(x, errB_map[:, slice_ind], "k-",color='black', label = 'SWIFT')
#axs[3].set_yticks(np.arange(-6.0, 1.0, 1.0))
#axs[3].set_xlabel(r"$x$")
#axs[3].set_ylabel(r"$\mathrm{log}_{10} \left( h \quad \nabla \cdot B / |B| \right)$")
#axs[3].set_ylim(-6.5, 0.5)

if with_reference:
    axs[0].plot(rho_vs_x_data[0],rho_vs_x_data[1],label='MFM $256^2$',color='red')
    #axs[1].plot(Bx_vs_x_data[0],Bx_vs_x_data[1],label='MFM $256^2$',color='red')
    #axs[2].plot(By_vs_x_data[0],By_vs_x_data[1],label='MFM $256^2$',color='red')
    #axs[3].plot(Berr_vs_x_data[0],Berr_vs_x_data[1],label='MFM $256^2$',color='red')

#axs[2].legend()
# Add panel with infromation about the run
Ninfo = 2
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=axs[Ninfo].transAxes
)

axs[Ninfo].text(
    0.5,
    0.8,
    "Orszag Tang Vortex at $t=%.2f$, $y_0 = %.4f$" % (data.metadata.time, y0),
    **text_common_args,
)
axs[Ninfo].text(0.5, 0.7, "SWIFT %s" % git.decode("utf-8"), **text_common_args)
axs[Ninfo].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
axs[Ninfo].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
axs[Ninfo].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
axs[Ninfo].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
axs[Ninfo].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
axs[Ninfo].text(
    0.5, 0.1, r"Physical resistivity: $\eta$: $%.2f$ " % (eta), **text_common_args
)
axs[Ninfo].text(
    0.5, 0.0, r"Number of particles: $N_p$: $%.0f$ " % (n_gas), **text_common_args
)

axs[Ninfo].axis("off")



for ax in axs:
    ax.minorticks_on()
    ax.grid()

fig.tight_layout()
plt.savefig(args.output)
