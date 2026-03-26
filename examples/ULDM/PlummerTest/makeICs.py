#!/usr/bin/env python3


import numpy as np
from pNbody import ic
from pNbody import Nbody
from astropy import units as u
import argparse

def generate_plummer_multi_grid(N_target, a=1.0, scales=[1.0, 5.0, 20.0], grid_res=100):
    """
    Generates a Plummer profile using nested grids to handle high density contrast.
    
    N_target : Total desired particles across all grids
    a        : Plummer scale radius
    scales   : List of L_max values (e.g., [1, 10, 100]) from smallest to largest
    grid_res : Resolution per axis (constant for each level to save RAM)
    """
    all_pos = []
    
    # 1. Calculate the total theoretical mass (integral of rho) 
    # to normalize N_target correctly across the entire volume.
    # The integral of the Plummer profile to infinity is proportional to 1.
    # We use a large scale limit to find the global normalization.
    max_L = scales[-1]
    
    # Pre-calculate global normalization factor
    # We want: N_total = Factor * Integral(rho dV)
    # For simplicity, we sample based on the fraction of mass in each shell.
    
    def get_rho(r_sq):
        return (1 + r_sq / a**2)**(-2.5)

    last_L = 0
    for L_max in sorted(scales):
        # Create the grid for this scale
        coord = np.linspace(-L_max, L_max, grid_res)
        xv, yv, zv = np.meshgrid(coord, coord, coord)
        r_sq = xv**2 + yv**2 + zv**2
        r = np.sqrt(r_sq)
        
        # MASK: Only keep points in the current "shell" 
        # (Inside L_max but outside the previous smaller L_max)
        # For the first grid (center), last_L is 0.
        shell_mask = (np.abs(xv) <= L_max) & (np.abs(yv) <= L_max) & (np.abs(zv) <= L_max)
        if last_L > 0:
            inside_mask = (np.abs(xv) <= last_L) & (np.abs(yv) <= last_L) & (np.abs(zv) <= last_L)
            shell_mask = shell_mask & ~inside_mask
            
        rho = get_rho(r_sq)
        
        # Calculate local probability
        # Volume of one cell in this grid
        dv = (2 * L_max / (grid_res - 1))**3
        
        # Expected number of particles in this grid level
        # We normalize so that the total sum over all space would be N_target
        # The analytical 3D integral of (1+r^2)^-2.5 is 4/3 * PI
        # Normalization factor C = N_target / (Integral of rho)
        # Integral of rho over all space is (4/3) * pi * a^3
        total_vol_integral = (4/3) * np.pi * (a**3)
        p_accept = (N_target / total_vol_integral) * rho * dv
        
        # Apply shell mask
        p_accept[~shell_mask] = 0
        
        # Draw particles
        rand_vals = np.random.uniform(0, 1, size=rho.shape)
        mask = rand_vals < p_accept
        
        level_pos = np.stack((xv[mask], yv[mask], zv[mask]), axis=-1)
        all_pos.append(level_pos)
        
        last_L = L_max
        print(f"Level L={L_max}: Generated {len(level_pos)} particles")

    # Combine all levels
    final_pos = np.concatenate(all_pos, axis=0)
    n = len(final_pos)
    
    # Standard Nbody object creation (as per your previous snippet)
    vel = np.zeros([n, 3])
    mass = np.ones([n]) * 1. / n # Assuming total mass M=1 distributed equally

    nb = Nbody(
        status='new',
        p_name="plum_multiscale.hdf5",
        pos=final_pos,
        vel=vel,
        mass=mass,
        ftype='swift')

    return nb

def generate_plummer_grid(N_target, a=1.0, L_max=10.0, grid_res=100):
    """
    Samples a Plummer profile on a regular grid while respecting N_target.
    N_target : Approximate number of desired particles
    a        : Scale radius parameter
    L_max    : Grid extension (from -L_max to L_max)
    grid_res : Grid resolution along one axis
    """

    # 1. Creation of the regular 3D grid
    coord = np.linspace(-L_max, L_max, grid_res)
    xv, yv, zv = np.meshgrid(coord, coord, coord)
    r_sq = xv**2 + yv**2 + zv**2
    
    # 2. Calculation of the Plummer density (volume mass profile)
    # rho(r) = (3M / 4*pi*a^3) * (1 + (r/a)^2)^(-5/2)
    # Constants are ignored as we normalize relative to N_target
    rho = (1 + r_sq / a**2)**(-2.5)
    
    # 3. Calculation of the normalization factor
    # The probability p of occupying a node is p = rho * factor
    # We want the sum of probabilities to equal N_target
    sum_rho = np.sum(rho)
    
    # Acceptance probability for each point
    p_accept = (N_target / sum_rho) * rho
    
    # Safety: probability cannot exceed 1
    if np.max(p_accept) > 1:
        print(f"Warning: grid_res is too low for N_target={N_target}.")
        print("Certain areas are saturated (p > 1). Increase grid_res.")
        p_accept = np.clip(p_accept, 0, 1)

    # 4. Random draw on the grid
    rand_vals = np.random.uniform(0, 1, size=rho.shape)
    mask = rand_vals < p_accept
    
    # Extraction of coordinates
    pos = np.stack((xv[mask], yv[mask], zv[mask]), axis=-1)

    n = len(pos)

    nb = Nbody(
        status='new',
        p_name="plum.hdf5",
        pos=pos,
        ftype='swift')

    return nb



####################################################################
# option parser
####################################################################

description="""Generate a Plummer sphere"""
epilog     ="""
Examples:
--------
"""

parser = argparse.ArgumentParser(description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("-n",
                  action="store",
                  dest="n",
                  type=int,
                  default=100000,
                  help="number of particles",
                  metavar=" INT")   

parser.add_argument("--regular-grid",
                  action="store_true",
                  default=False,
                  help="Put particles on a regular grid") 

parser.add_argument("-N",
                   action="store",
                   type=int,
                   default=256,
                   help="Grid resolution: number of points per axis (total nodes = N^3)",
                   metavar=" INT") 

parser.add_argument("-o",
                    action="store",
                    type=str,
                    dest="outputfilename",
                    default='snap.hdf5',
                    help="Name of the output file")  


opt = parser.parse_args()

 
# define units
u_Length   = 1* u.kpc
u_Mass     = 10**10 * u.M_sun
u_Velocity = 1* u.km/u.s
u_Time     = u_Length/u_Velocity
toMsol     = u_Mass.to(u.M_sun).value
G          = 4.299581e+04

M = 1        # 1e10 Msol
a = 0.1      # scale radius
R = 50       # kpc
Rt = 5       # truncation radius
L = 20       # boxsize [kpc}

outputfile = opt.outputfilename


# create the pNbody object

if opt.regular_grid:
  #nb = generate_plummer_grid(N_target=10000, a=a, L_max=1, grid_res=opt.N)
  nb = generate_plummer_multi_grid(opt.n, a=a, scales=[0.01*Rt,0.03*Rt,0.1*Rt, 0.3*Rt, Rt], grid_res=opt.N)
else:
  nb = ic.plummer(opt.n,1,1,1,a,R,irand=1,ftype='swift')

nb.mass = np.ones(nb.nbody)*M/nb.nbody
nb = nb.selectc(nb.rxyz()<Rt)

hydro=True
if hydro:
  # set particles to gas
  nb.set_tpe(0)
  nb.u_init   = np.ones(nb.nbody)*0.05  
  nb.rsp_init = nb.get_rsp_approximation()
  #nb.rsp_init = np.clip(nb.rsp_init,0,5)
else:
  # set particles to dm
  nb.set_tpe(1)
  

# add units
nb.UnitLength_in_cm         = u_Length.to(u.cm).value
nb.UnitMass_in_g            = u_Mass.to(u.g).value
nb.UnitVelocity_in_cm_per_s = u_Velocity.to(u.cm/u.s).value
nb.Unit_time_in_cgs         = u_Time.to(u.s).value

nb.hubblefactorcorrection      = False
nb.comovingtoproperconversion  = False
nb.atime                       = 1

nb.boxsize = np.array([L,L,L])

nb.rename(outputfile)
nb.write()





