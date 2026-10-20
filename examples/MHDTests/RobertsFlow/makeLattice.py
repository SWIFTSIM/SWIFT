"""
Create a lattice file

Nikyta Shchutskyi 2026

based on makeBCC.py from Brio Wu test

"""

import numpy as np
import h5py

from unyt import cm, g, s, erg
from unyt.unit_systems import cgs_unit_system

from swiftsimio import Writer


def generate_sc(n, Lbox):

    """

    Generate SC lattice in [0, side_length)^3

    """
    a = Lbox/n

    # Primitive cubic grid
    gridX = np.linspace(0.0, Lbox[0], nx, endpoint=False)
    gridY = np.linspace(0.0, Lbox[1], ny, endpoint=False)
    gridZ = np.linspace(0.0, Lbox[2], nz, endpoint=False)

    X, Y, Z = np.meshgrid(gridX, gridY, gridZ, indexing='ij')

    corners = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
 
    positions = corners

    shift = a / 2

    positions = (positions + shift[None,:]) % Lbox[None,:] #side_length

    #positions %= side_length

    return positions


def generate_bcc(n, Lbox):

    """

    Generate BCC lattice in [0, side_length)^3

    """
    a = Lbox/n

    # Primitive cubic grid
    gridX = np.linspace(0.0, Lbox[0], nx, endpoint=False)
    gridY = np.linspace(0.0, Lbox[1], ny, endpoint=False)
    gridZ = np.linspace(0.0, Lbox[2], nz, endpoint=False)

    X, Y, Z = np.meshgrid(gridX, gridY, gridZ, indexing='ij')

    corners = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
    
    # Body-centered points (shifted by half a cell)

    centers = corners + a[None,:]/2

    # Combine both
    
    positions = np.vstack([corners, centers])
    
    shift = a / 4

    positions = (positions + shift[None,:]) % Lbox[None,:] #side_length

    return positions


def generate_fcc(num_on_side, side_length=1.0):

    """

    Generate BCC lattice in [0, side_length)^3

    """
    a = Lbox/n

    # Primitive cubic grid
    gridX = np.linspace(0.0, Lbox[0], nx, endpoint=False)
    gridY = np.linspace(0.0, Lbox[1], ny, endpoint=False)
    gridZ = np.linspace(0.0, Lbox[2], nz, endpoint=False)

    X, Y, Z = np.meshgrid(gridX, gridY, gridZ, indexing='ij')

    corners = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
        
    # Face-centered shifts

    shift_xy = np.array([a[0]/2, a[1]/2, 0.0])

    shift_xz = np.array([a[0]/2, 0.0, a[2]/2])

    shift_yz = np.array([0.0, a[1]/2, a[2]/2])

    face_xy = corners + shift_xy

    face_xz = corners + shift_xz

    face_yz = corners + shift_yz

    # Combine all

    positions = np.vstack([corners, face_xy, face_xz, face_yz])

    # Periodic wrap

    shift = a / 4

    positions = (positions + shift[None,:]) % Lbox[None,:]

    #positions %= side_length

    return positions


def generate_distortion(cube, noise_width, n, Lbox):

    a = Lbox/n
    dmin = np.min(Lbox/n)

    noise_amp = dmin
    noise_width *= dmin 

    noise = noise_amp * np.random.normal(loc=0.0, scale=noise_width, size=cube.shape)

    cube_noisy = cube + noise

    cube_noisy = cube_noisy % Lbox[None,:]

    return cube_noisy


def write_out_glass(filename, cube, Lbox):

    N = len(cube)
    h = (Lbox[0]*Lbox[1]*Lbox[2]/N)**(1/3) * np.ones(N) 

    print('Np total',N)
 
    fileOutput = h5py.File(filename, "w")

    # Header
    grp = fileOutput.create_group("/Header")
    grp.attrs["BoxSize"] = Lbox  #####
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3

    # Units
    grp = fileOutput.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = fileOutput.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=cube, dtype="d")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")

    fileOutput.close()


    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="Generate a BCC lattice")

    #parser.add_argument(
    #    "-n",
    #    "--numparts",
    #    help="Number of particles on a side. Default: 32",
    #    default=32,
    #    type=int,
    #)
    parser.add_argument(
        "-LT",
        "--lattice_type",
        help="Type of the lattice: SC, BCC, FCC, HCP",
        default='SC',
        type=str,
    )
    parser.add_argument(
        "-nx",
        "--npart_x",
        help="Number of particles in x direction",
        default=16,
        type=int,
    )
    parser.add_argument(
        "-ny",
        "--npart_y",
        help="Number of particles in x direction",
        default=16,
        type=int,
    )
    parser.add_argument(
        "-nz",
        "--npart_z",
        help="Number of particles in x direction",
        default=16,
        type=int,
    )
    parser.add_argument(
        "-noiseW",
        "--noise_width",
        help="relative width of the noise",
        default=0.0,
        type=float,
    )



    import_to_vtk=False

    args = parser.parse_args()

    # set box geometry
    nx = args.npart_x
    ny = args.npart_y
    nz = args.npart_z
    n = np.array([nx,ny,nz])
    if ( (nx<nz) | (ny<nz) ):
        maxNperp = max([nx,ny])
        Lbox = np.array([nx/maxNperp, ny/maxNperp, nz/maxNperp])
    else:
        maxN = max([nx,ny,nz])
        Lbox = np.array([nx/maxN, ny/maxN, nz/maxN])

    LT = args.lattice_type
    if LT=='SC':
        cube = generate_sc(n, Lbox)
        output = 'SCcube_'+str(nx)+'_'+str(ny)+'_'+str(nz)+'.hdf5'
    elif LT=='BCC':
        cube = generate_bcc(n, Lbox)
        output = 'BCCcube_'+str(nx)+'_'+str(ny)+'_'+str(nz)+'.hdf5'
    elif LT=='FCC':
        cube = generate_fcc(n, Lbox)
        output = 'FCCcube_'+str(nx)+'_'+str(ny)+'_'+str(nz)+'.hdf5'

    # add distortion
    cube = generate_distortion(cube,args.noise_width,n,Lbox) 

    if import_to_vtk:
        import vtk
        vtk_points = vtk.vtkPoints()
        for x, y, z in cube:
            vtk_points.InsertNextPoint(x, y, z)
        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(vtk_points)
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName("output.vtk")
        writer.SetInputData(poly_data)
        writer.Write()  # one can use paraview to watch how particles are arranged

    write_out_glass(output, cube, Lbox)
