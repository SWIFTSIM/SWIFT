#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import os
from scipy import integrate, optimize
import numpy as np
from scipy.interpolate import interp1d
import vtk


# Files to read from and write to



def open_snapshot(fileName):
    
    """ Open the input file and read the header and particle positions 
    param: fileName - the name of the input file
    return: a dictionary containing the BoxSize, afact, N_in, and pos_in 
    """

    infile = h5py.File(fileName, "r")

    head = infile["/Header"]
    units = infile["/Units"]
    BoxSize = head.attrs["BoxSize"]
    NumPart_Total = head.attrs["NumPart_Total"]
    NumPart_ThisFile = head.attrs["NumPart_ThisFile"]

    UI = units.attrs['Unit current in cgs (U_I)']
    UL = units.attrs['Unit length in cgs (U_L)']
    UM = units.attrs['Unit mass in cgs (U_M)']
    UT = units.attrs['Unit temperature in cgs (U_T)']  
    UTM = units.attrs['Unit time in cgs (U_t)']  

    # Read PartType0 data
    pos = infile["/PartType0/Coordinates"][:, :]
    u = infile["/PartType0/InternalEnergy"][:]
    m = infile["/PartType0/Masses"][:]
    ids = infile["/PartType0/ParticleIDs"][:]
    h = infile["/PartType0/SmoothingLength"][:]
    vel = infile["/PartType0/Velocities"][:, :]

    infile.close()


    data = {
        'PartType0':
        {
        'pos': pos,
        'u': u,
        'm': m,
        'ids': ids,
        'h': h,
        'vel': vel
        },

        'Units':
        {
            'UI': UI,
            'UL': UL,
            'UM': UM,
            'UT': UT,
            'UTM': UTM
        },

        'BoxSize': BoxSize,
        'NumPart_Total':NumPart_Total,
        'NumPart_ThisFile':NumPart_ThisFile,
    }

    return data



def add_atmospere_particles(data,
                            pars,
                            R_max,  # Maximum radius of the atmosphere in BoxSize units
                            IAsource='IAfile',  # Source of atmosphere particles, e.g., 'grid' or 'random',
                            ):
    """ Add atmosphere particles to the data dictionary
    param: data - the dictionary containing BoxSize, afact, N_in, and pos_in
    param: pars - parameters for density profile
    param: IAsource - mode for atmosphere initial particle arrangement
    return: the updated data dictionary with atmosphere particles added
    """
   
    mp = np.mean(data['PartType0']['m'])  # Mean mass of a particle in the simulation
    #BoxSize = data['BoxSize'][0]  # Assuming BoxSize is a 1D array

    # Estimate required number of particles
    Natm = estimate_Natm(pars['density'], mp, R_max)
    print('Estimated number of atmosphere particles:', Natm)
    if Natm > 1e6:
        print('Warning: Number of atmosphere particles is very high')
        exit(1)

    # Calculate required number of particles for initial arrangement file
    Natm = int(Natm * 6 / np.pi) # Renormalize for spheical cut

    # Make initial particle arrangement
    pos= np.zeros((Natm, 3))  # Initialize positions of atmosphere particles
    if IAsource == "grid":
        Nside = int(Natm ** (1 / 3))
        grid_1d = np.linspace(0.0, 1.0, Nside)  # Cube side is 1, centered at (0.5,0.5,0.5)
        # Create a 3D grid of coordinates
        x, y, z = np.meshgrid(grid_1d, grid_1d, grid_1d, indexing="ij")
        # Reshape into a list of points
        pos = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T
        pos += 0.5 / Nside
    elif IAsource == "IAfile":
        # Loading initial arrangement file
        IAfile = "glassCube_128.hdf5"
        pos, _ = open_IAfile(IAfile)
        lcut = (Natm / len(pos)) ** (1 / 3)
        pos -= 0.5
        mask = (
            (np.abs(pos[:, 0]) <= lcut / 2)
            & (np.abs(pos[:, 1]) <= lcut / 2)
            & (np.abs(pos[:, 2]) <= lcut / 2)
        )
        pos = pos[mask]
        pos /= lcut
        pos += 0.5
    else:
        print('not impl yet')

    pos -= 0.5 # shift initial arrangement center at zero, prepare for cutting and rescaling

    # Cut a sphere
    r_uni = np.linalg.norm(pos, axis=1)
    pos = pos[r_uni <= 0.5]
    print('len par is',len(pos))

    # Rescale particle positions to match the required density profile
    pos *= 2*R_max # Scale positions to make sphere with radius = Rmax
    rtab, Ftab = tabulate_F_function(pars['density'], R_max) 
    F_rho = interp1d(rtab, Ftab, kind='linear', bounds_error=False, fill_value=(Ftab[0], Ftab[-1]))  # Build CDF function for the selected density
    r_uni = np.linalg.norm(pos, axis=1)
    F_uni = r_uni**3 / R_max**3 # data['BoxSize'][0]**3
    r_target = np.array([invert_density_cdf(F, F_rho, R_max,1e-6*R_max) for F in F_uni]) # Invert CDF
    pos *= (r_target / r_uni)[:, np.newaxis] # Rescale coordinates
    pos += data['BoxSize'] / 2  # Shift positions to be centered at BoxSize/2


    # cut atmosphere particles close to the galactic disk
    pos = cut_disk_region_from_atmosphere( data, pos)

    # Debug: printing the positions of the particles
    #print(pos[mask_cut])
    print(len(pos))

    # write particle positions
    data['PartType0']['pos'] = np.concatenate((data['PartType0']['pos'], pos), axis=0)

    # write particle smoothing lengths
    r_uni = np.linalg.norm(pos, axis=1)
    h = 1.595*(mp/rho(r_uni,pars['density']))**(1/3)
    data['PartType0']['h'] = np.concatenate((data['PartType0']['h'], h), axis=0)

    # write particle masses
    m = mp*np.ones(pos.shape[0])
    data['PartType0']['m'] = np.concatenate((data['PartType0']['m'],m), axis=0)

    # write particle indexes
    idsmax = np.max(data['PartType0']['ids'])
    ids = np.arange(idsmax+1,idsmax+1+pos.shape[0])
    data['PartType0']['ids'] = np.concatenate((data['PartType0']['ids'], ids), axis=0)
    
    # write particle energies
    umin = np.min( data['PartType0']['u'] )
    u = umin*np.ones(pos.shape[0])
    data['PartType0']['u'] = np.concatenate((data['PartType0']['u'], u), axis=0)

    # write particle velocities
    v = np.zeros(pos.shape)
    data['PartType0']['vel'] = np.concatenate((data['PartType0']['vel'], v), axis=0)

    # update particle count
    data['NumPart_Total'][0] = data['PartType0']['pos'].shape[0]
    data['NumPart_ThisFile'][0] = data['PartType0']['pos'].shape[0]

    return data



import numpy as np
from scipy.spatial import cKDTree

def cut_disk_region_from_atmosphere(data, pos_atm):

    """ 
    Cut atmosphere particles which lie inside the disk
    param: pos_disk - positions of the disk particles
    param: h_disk - smoothing lengths of the disk particles
    param: pos_atm - positions of the atmosphere particles
    return: positions outside the disk
    """

    # read atmosphere particles positions
    pos_disk = data['PartType0']['pos']
    h_disk  = data['PartType0']['h']

    # Build KDTree from disk positions
    tree_disk = cKDTree(pos_disk)

    # Query the closest disk particle for each atmosphere particle
    distances, indices = tree_disk.query(pos_atm, k=1)  # Nearest neighbor

    # Compare distance to h_disk of that disk particle
    inside_disk = distances < 3*h_disk[indices]  # Boolean mask

    # Mask atmosphere particles outside disk region
    pos_atm_outside_disk = pos_atm[~inside_disk]

    return pos_atm_outside_disk 



def open_IAfile(path_to_file):
    """ 
    Open initial particle arrangement file (used for atmosphere generation)
    param: path_to_file - path to .hdf5
    return: positions and smoothing lengths of particles
    """
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h



def estimate_Natm(pars, mp, R_max):
    """ 
    Estimate number of atmosphere particles given particle mass and atmosphere radius
    param: pars - parameters for density profile in snapshot units
    param: mp - particle mass in snapshot units
    param: R_max - atmosphere radius in snapshot units
    return: number of particles in the atmosphere
    """
    # Compute total mass in the profile
    def integrand(r):
        return rho(r, pars) * 4 * np.pi * r**2

    M_tot, _ = integrate.quad(integrand, 1e-6*R_max, R_max)

    Natm = int(M_tot / mp)
    return Natm



def rho(r, pars):
    """ 
    Calculate the density profile at given positions
    param: r - distance to the center
    param: pars - parameter dict for density profile.
        pars['profile'] - name of the profile
        other fields - parameters specific for the profile
    return: the density at the given positions
    """

    # Define density profiles
    if pars['profile']=='uniform':
        # Constant density profile
        rho0 = pars['rho0']
        density= rho0 * np.ones_like(r)
    elif pars['profile']=='isothermal':
        # Isothermal profile
        rho0 = pars['rho0']
        r0 = pars['r0']
        density = rho0 / (1 + (r/r0)**2)
    elif pars['profile']=='beta':
        # Beta profile
        rho0 = pars['rho0']
        r0 = pars['r0']
        beta = pars['beta']
        density = rho0 * (1 + (r/r0)**2)**( - 3 * beta / 2)
   
    return density


def tabulate_F_function(pars, R_max, N_points=1000):
    """ 
    Tabulate the mass function for a given density profile
    param: pars - parameters of the density profile
    param: R_max - maximum radius for the mass function
    param: N_points - number of points to sample in the mass function
    return: a function that computes the cumulative mass distribution
    """
    r = np.linspace(1e-6*R_max, R_max, N_points)
    M = np.array([integrate.quad(lambda s: rho(s, pars) * s**2, 0, ri)[0] for ri in r])
    M /= M[-1]  # Normalize to 1 at R_max

    return r, M



def invert_density_cdf(F_uni, F_rho, R_max, Rmin):
    """ Invert cdf function to get position rescaling
    param: F_uni - uniform density CDF of particles
    param: F_rho - chosen density profile CDF
    param: R_max - maximal radius of the atmosphere
    param: R_min - minimal radius, regularizer
    return: a function that computes the cumulative mass distribution
    """
    def objective(r):
        return F_rho(r) - F_uni
    return optimize.root_scalar(objective, bracket=[Rmin, R_max], method="bisect").root



def add_magnetic_fields(data,
                        B0_Gaussian_Units=1e-3,  # micro Gauss
                        type="uniform",
                        ):
    """ 
    Add magnetic fields to the data dictionary
    param: data - the dictionary containing BoxSize, afact, N_in, and pos_in
    param: B0_Gaussian_Units - the magnitude of the uniform seed magnetic field in Gaussian units
    param: type - the type of magnetic field to add (default is "uniform")
    return: the updated data dictionary with magnetic fields added
    """
  
    # initialize the magnetic field and vector potential arrays
    Ngas = data['PartType0']['m'].shape[0]
    B = np.zeros((Ngas, 3))
    A = np.zeros((Ngas, 3))

    # Set the magnitude of the uniform seed magnetic field
    CONST_MU0_CGS = 4 * np.pi
    B0_cgs = np.sqrt(CONST_MU0_CGS / (4.0 * np.pi)) * B0_Gaussian_Units

    # Set the magnetic field and vector potential based on the type
    if type == "uniform":
        B[:,0] = B0_cgs

    # Set the magnetic field unit
    const_unit_magnetic_field_in_cgs = 1e-7  # 1 micro Gauss in cgs
    const_unit_mass_in_cgs = 1.9891e43
    const_unit_length_in_cgs = 3.08567758e21
    const_unit_velocity_in_cgs = 1e5

    # Derived quantities
    const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
    const_unit_current_in_cgs = const_unit_mass_in_cgs / (
        const_unit_magnetic_field_in_cgs * const_unit_time_in_cgs ** 2
    )

    # Store the magnetic field and vector potential in the data dictionary
    data['PartType0']["B"] = B
    data['PartType0']["A"] = A
    data['Units']['UI'] = const_unit_current_in_cgs  # Current unit in cgs

    return data


def print_to_vtk(data, fileName):
    """ 
    Print the data dictionary to a VTK file
    param: data - the dictionary containing BoxSize, afact, N_in, and pos_in
    param: fileName - the name of the output VTK file
    """

    vtk_points = vtk.vtkPoints()
    for x, y, z in data['PartType0']['pos'] - data['BoxSize'] / 2:
        vtk_points.InsertNextPoint(x, y, z)

    print("Number of particles in VTK file:", vtk_points.GetNumberOfPoints())
    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(vtk_points)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fileName)
    writer.SetInputData(poly_data)
    writer.Write()  # one can use paraview to watch how particles are arranged

    return


def makeIC(fileInputName, fileOutputName, PrintToVTK = True):

    """
    Make a new file from an input snapshot file
    param: fileInputName - the name of the snapshot to be read
    param: fileOutputName - the name of the output snapshot file
    param: PrintToVTK - whether to print the particles to a VTK file
    return: None
    """

    # Load snapshot
    data = open_snapshot(fileInputName)   

    # Atmosphere parameter setup. Profile follows U. P. Steinwandel et al. 2019
    rho_units_snapshot = data['Units']['UM'] * data['Units']['UL']**(-3)
    rho0 = 5e-26 / rho_units_snapshot  
    R_max = 50 # maximal radius in kpc
    pars = {'density':{'profile': 'beta', 'rho0': rho0, 'r0': 0.33, 'beta': 2/3 },
            'temperature':{'profile': 'constant', 'T0': 1e4 }
            }  # form parameter dict

    # Add atmosphere
    data = add_atmospere_particles(data, pars, R_max=R_max)

    # Debug: save particle positions to .vtk
    if PrintToVTK:
        print_to_vtk(data, 'fid.vtk')

    # Add magnetic fields
    data = add_magnetic_fields(data,
                               B0_Gaussian_Units=1e-3,  # micro Gauss        
                               )


    # Generate and open output file
    os.system("cp " + fileInputName + " " + fileOutputName)
    fileOutput = h5py.File(fileOutputName, "a")

    # Particle group
    grp = fileOutput.require_group("/PartType0")


    Ntot = len(data['PartType0']['pos'])
    write_or_replace_dataset(grp, "Coordinates", data['PartType0']['pos'], size=(Ntot,3), dtype="d")
    write_or_replace_dataset(grp, "Velocities", data['PartType0']['vel'], size=(Ntot,3), dtype="f")
    write_or_replace_dataset(grp, "InternalEnergy", data['PartType0']['u'], size=(Ntot,), dtype="f")
    write_or_replace_dataset(grp, "SmoothingLength", data['PartType0']['h'], size=(Ntot,), dtype="f")
    write_or_replace_dataset(grp, "ParticleIDs", data['PartType0']['ids'], size=(Ntot,), dtype="L")
    write_or_replace_dataset(grp, "Masses", data['PartType0']['m'], size=(Ntot,), dtype="f")
    write_or_replace_dataset(grp, "MagneticFluxDensities", data['PartType0']['B'], size=(Ntot,3), dtype="f")
    write_or_replace_dataset(grp, "MagneticVectorPotentials", data['PartType0']['A'], size=(Ntot,3), dtype="f")

   # grp.create_dataset("MagneticFluxDensities", data=data['PartType0']['B'], dtype="f")
    #grp.create_dataset("MagneticVectorPotentials", data=data['PartType0']['A'], dtype="f")

    # Update gas particle properties

    # Change current unit to something more resonable
    unitSystem = fileOutput["/Units"]
    unitSystem.attrs.modify("Unit current in cgs (U_t)", data['Units']['UTM'])
    unitSystem.attrs.modify("Unit current in cgs (U_L)", data['Units']['UL'])
    unitSystem.attrs.modify("Unit current in cgs (U_M)", data['Units']['UM'])
    unitSystem.attrs.modify("Unit current in cgs (U_T)", data['Units']['UT'])
    unitSystem.attrs.modify("Unit current in cgs (U_I)", data['Units']['UI'])


    grp = fileOutput.require_group("/Header")
    grp.attrs.modify("NumPart_Total", data['NumPart_Total'])
    grp.attrs.modify("NumPart_ThisFile", data['NumPart_ThisFile'])


    fileOutput.close()

def write_or_replace_dataset(group, name, datai, size, dtype):

    """ Write or replace a field in .hdf5
    param: group - group to write dataset into
    param: name - name of dataset
    param: datai - data to write
    param: size - size of dataset
    param: dtype - data type of dataset
    return: None
    """

    if name in group:
        del group[name]
    ds = group.create_dataset(name, size, dtype=dtype)
    ds[()] = datai

#fileInputName = sys.argv[1]
#fileOutputName = sys.argv[2]

makeIC('fid.hdf5', 'fid_B.hdf5')
#makeIC(fileInputName,fileOutputName)



# References:
# 1. U P Steinwandel, M C Beck, A Arth, K Dolag, B P Moster, P Nielaba, 
#    Magnetic buoyancy in simulated galactic discs with a realistic circumgalactic medium, 
#    Monthly Notices of the Royal Astronomical Society, Volume 483, Issue 1, February 2019, Pages 1008–1028, 
#    https://doi.org/10.1093/mnras/sty3083