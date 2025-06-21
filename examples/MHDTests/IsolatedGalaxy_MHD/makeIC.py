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


# Files to read from and write to
#fileInputName = sys.argv[1]
#fileOutputName = sys.argv[2]



def open_snapshot(fileName):
    
    """ Open the input file and read the header and particle positions 
    param: fileName - the name of the input file
    return: a dictionary containing the BoxSize, afact, N_in, and pos_in 
    """

    infile = h5py.File(fileName, "r")

    head = infile["/Header"]
    units = infile["/Units"]
    BoxSize = head.attrs["BoxSize"]
    N_in = head.attrs["NumPart_Total"][0]

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
        'N_in': N_in,
        
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
    Natm = estimate_Natm(pars, mp, R_max)
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
    rtab, Ftab = tabulate_F_function(pars, R_max) 
    F_rho = interp1d(rtab, Ftab, kind='linear', bounds_error=False, fill_value=(Ftab[0], Ftab[-1]))  # Build CDF function for the selected density
    r_uni = np.linalg.norm(pos, axis=1)
    F_uni = r_uni**3 / R_max**3 # data['BoxSize'][0]**3
    r_target = np.array([invert_density_cdf(F, F_rho, R_max,1e-6*R_max) for F in F_uni]) # Invert CDF
    pos *= (r_target / r_uni)[:, np.newaxis] # Rescale coordinates
    pos += data['BoxSize'] / 2  # Shift positions to be centered at BoxSize/2

    # adding atmosphere particles to the data
    data['PartType0']['pos'] = np.concatenate((data['PartType0']['pos'], posatm), axis=0)

    return data



def open_IAfile(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h



def estimate_Natm(pars, mp, R_max):
    """ Estimate number of atmosphere particles given particle mass and atmosphere radius
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
    """ Calculate the density profile at given positions
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
        density = rho0 / (1 + (r/r0)**2)**(3 * beta / 2)
   
    return density


def tabulate_F_function(pars, R_max, N_points=1000):
    """ Tabulate the mass function for a given density profile
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
    """ Add magnetic fields to the data dictionary
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
    """ Print the data dictionary to a VTK file
    param: data - the dictionary containing BoxSize, afact, N_in, and pos_in
    param: fileName - the name of the output VTK file
    """
    import vtk

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


def makeIC(fileInputName, fileOutputName):

    data = open_snapshot(fileInputName)   

    # Atmosphere parameter setup
    m_H = 1.6726219e-24 / data['Units']['UM']  # Hydrogen mass in snapshot units
    n0 = 1e-4 / (data['Units']['UL']**(-3))  # Example number density in snapshot units
    rho0 = n0 * m_H  # Convert number density to mass density

    R_max = 60 # maximal radius
    pars = {'profile': 'isothermal', 'rho0': rho0, 'r0':10 }  # form parameter dict

    # Add atmosphere
    data = add_atmospere_particles(data, pars, R_max=R_max)

    # Debug: save atmosphere to .vtk
    print('Atmosphere particles added:', posatm)
    print_to_vtk(data, 'fid.vtk')

    # Add magnetic fields
    # data = add_magnetic_fields(data,
    #                            B0_Gaussian_Units=1e-3,  # micro Gauss        
    #                            )

    os.system("cp " + fileInputName + " " + fileOutputName)

    # File
    fileOutput = h5py.File(fileOutputName, "a")

    # Particle group
    grp = fileOutput.require_group("/PartType0")
    #grp.create_dataset("MagneticFluxDensities", data=data['PartType0']['B'], dtype="f")
    #grp.create_dataset("MagneticVectorPotentials", data=data['PartType0']['A'], dtype="f")

    # Change current unit to something more resonable
    unitSystem = fileOutput["/Units"]
    unitSystem.attrs.modify("Unit current in cgs (U_I)", data['Units']['UI'])
    fileOutput.close()

makeIC('fid.hdf5', 'fid_B.hdf5')