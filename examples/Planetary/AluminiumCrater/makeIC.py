"""
Creates the BCC ICs for the blob test.
"""

import numpy as np
import h5py

import sys


def generate_cube(num_on_side, side_length=1.0):
    """
    Generates a cube
    """

    values = np.linspace(0.0, side_length, num_on_side + 1)[:-1]

    positions = np.empty((num_on_side ** 3, 3), dtype=float)

    for x in range(num_on_side):
        for y in range(num_on_side):
            for z in range(num_on_side):
                index = x * num_on_side + y * num_on_side ** 2 + z

                positions[index, 0] = values[x]
                positions[index, 1] = values[y]
                positions[index, 2] = values[z]

    return positions


def generate_bcc_lattice(num_on_side, side_length=1.0):
    cube = generate_cube(num_on_side // 2, side_length)

    mips = side_length / num_on_side

    positions = np.concatenate([cube, cube + mips])
        
    positions += 0.5 / num_on_side
    
    return positions


def generate_surface(
    num_on_side_highres, num_on_side_lowres, num_copies_highres, num_copies_lowres_eitherside, num_copies_y, side_length, mass_highres, mass_lowres, rho_surface,
):
    transition_regions = int(np.log2(num_on_side_highres) - np.log2(num_on_side_lowres) - 1)
    transition_resolutions = [int(2**(i + np.log2(num_on_side_lowres) + 1)) for i in range(transition_regions)]
    transition_lattices = []
    
    for i in range(transition_regions):
         num_on_side_i = transition_resolutions[i]    
         #print(num_on_side_i)   
         lattice_i = generate_bcc_lattice(num_on_side_i)                    
         transition_lattices.append(lattice_i)                   
    
    bcc_lattice_highres = generate_bcc_lattice(num_on_side_highres)
    bcc_lattice_lowres = generate_bcc_lattice(num_on_side_lowres)
    
    lower_boundary_highres = num_copies_lowres_eitherside
    upper_boundary_highres = num_copies_lowres_eitherside + num_copies_highres                         
    bcc_lattice_full = []
    masses = []
    
    
    for i in range(num_copies_highres + 2 * num_copies_lowres_eitherside):
        for j in range(num_copies_highres + 2 * num_copies_lowres_eitherside):
            if (lower_boundary_highres <= i < upper_boundary_highres) and (lower_boundary_highres <= j < upper_boundary_highres):
                
                bcc_lattice_full.append(list(bcc_lattice_highres + i * np.array([1, 0, 0])  + j * np.array([0, 0, 1])))
                masses.append(mass_highres * np.ones_like(bcc_lattice_highres[:,0]))
               
            elif (lower_boundary_highres - transition_regions <= i < upper_boundary_highres + transition_regions) and (lower_boundary_highres - transition_regions <= j < upper_boundary_highres + transition_regions):
                 #print("here")
                 
                 for i_transition in range(transition_regions):

                     lower_boundary_transitions = lower_boundary_highres - i_transition - 1
                     upper_boundary_transitions = upper_boundary_highres + i_transition + 1       
                     if (lower_boundary_transitions <= i < upper_boundary_transitions) and (lower_boundary_transitions <= j < upper_boundary_transitions):   

                             transition_idx = transition_regions - i_transition - 1
        
                             bcc_lattice_full.append(list(transition_lattices[transition_idx] + i * np.array([1, 0, 0])  + j * np.array([0, 0, 1])))
                             #print(transition_regions)
                            # print(transition_idx)
                            # print(transition_resolutions[transition_idx])
                             mass_transition = rho_surface * (side_length / (transition_resolutions[transition_idx]  / 2))**3 / 2
                             masses.append(mass_transition * np.ones_like(transition_lattices[transition_idx][:,0]))
                                
                             break
                 

            else:
                             
                bcc_lattice_full.append(list(bcc_lattice_lowres + i * np.array([1, 0, 0])  + j * np.array([0, 0, 1])))
                masses.append(mass_lowres * np.ones_like(bcc_lattice_lowres[:,0]))


    bcc_lattice_full = np.concatenate(bcc_lattice_full)
    masses = np.concatenate(masses)
    
    
    bcc_lattice_full = np.concatenate(
        [bcc_lattice_full + x * np.array([0.0, 1.0, 0.0]) for x in range(num_copies_y)]
    )
    masses = np.concatenate(
        [masses for x in range(num_copies_y)]
    )

    # Now size it appropriately
    bcc_lattice_full *= side_length


    positions = bcc_lattice_full
    
    

    # Now we can generate the velocities
    velocities = np.zeros_like(positions)
    
    return positions, velocities, masses


def generate_impactor(num_on_side, side_length, num_copies_highres, num_copies_lowres_eitherside, num_copies_y, impactor_radius, speed, angle, mass_highres):
    """
    Generate the positions and velocities for the blob.
    """

    factor = 2 * int(impactor_radius / side_length + 1)
    
    print(factor)
    
    bcc_lattice = generate_bcc_lattice(factor * num_on_side)

    # Update to respect side length
    bcc_lattice *= factor * side_length

    # Find the radii
    squared_radius = np.sum((bcc_lattice - 0.5 * factor * side_length) ** 2, axis=1)
    inside_sphere = squared_radius <= (impactor_radius)**2

    # Now select out particles
    bcc_lattice = bcc_lattice[inside_sphere]

    # Move to the correct x_position
    #bcc_lattice[:, 0] = bcc_lattice[:, 0] + blob_centre_x - 0.5 * num_copies_x_surface *side_length
    
    num_copies_surface = num_copies_highres + 2 * num_copies_lowres_eitherside
    centre = 0.5 * num_copies_surface * side_length
    bcc_lattice[:, 0] = bcc_lattice[:, 0] + centre - 0.5 * factor * side_length
    bcc_lattice[:, 1] = bcc_lattice[:, 1] + num_copies_y * side_length + impactor_radius - 0.5 * factor * side_length
    bcc_lattice[:, 2] = bcc_lattice[:, 2] + centre - 0.5 * factor * side_length
 #   buffer = 2 * impactor_radius
 #   bcc_lattice[:, 0] += buffer * np.sin(angle)
 #   bcc_lattice[:, 1] += buffer * np.cos(angle)
    
    
    masses = mass_highres * np.ones_like(bcc_lattice[:,0])
    
    positions = bcc_lattice

    # Generate velocities
    velocities = np.zeros_like(positions)
    velocities[:,0] = -speed * np.cos(angle)
    velocities[:,1] = -speed * np.sin(angle)

    return positions, velocities, masses


def write_out_ics(
    filename,
    num_on_side_highres,
    side_length=1.2e-3,#10e3,
    duplications_highres=10,
    duplications_lowres_eitherside=10,
    duplications_y=10,
    impactor_radius=0.5e-3,
):
    
    num_on_side_lowres = 4
    
  
    pppr = num_on_side_highres * (impactor_radius / side_length) / np.sqrt(2) # for bcc
    print("pppr = " + str(pppr))
    
    
    #See Ong et al. 2010 for these 
    
    mat_id_surface = 1000
    rho_surface = 2700 
    u_surface = 1e5#1#small
    
    mat_id_impactor = 1000
    rho_impactor = 2700 
    u_impactor = 1e5#1#small
    
    
    density_ratio = rho_impactor / rho_surface
    
    mass_highres = rho_surface * (side_length / (num_on_side_highres / 2))**3 / 2
    mass_lowres = rho_surface * (side_length / (num_on_side_lowres / 2))**3 / 2
    

    

    positions_surface, velocities_surface, masses_surface = generate_surface(
        num_on_side_highres,
        num_on_side_lowres,
        duplications_highres,
        duplications_lowres_eitherside,
        duplications_y,
        side_length,
        mass_highres,
        mass_lowres,
        rho_surface,
    )
    
 

    speed = 5e3
    angle = 50
    angle *= (np.pi / 180)
    
    positions_impactor, velocities_impactor, masses_impactor = generate_impactor(
        int(num_on_side_highres * np.cbrt(density_ratio)),
        side_length,
        duplications_highres,
        duplications_lowres_eitherside,
        duplications_y,
        impactor_radius,
        speed,
        angle,
        mass_highres,
    )
    
    # This is just an estimate to get the crater at centre of high res
    offset = 0.25 * side_length * duplications_highres * np.cos(angle)
    positions_impactor[:, 0] += offset
    
    coordinates = np.concatenate([positions_surface, positions_impactor])
    velocities = np.concatenate([velocities_surface, velocities_impactor])
    
    masses = np.concatenate([masses_surface, masses_impactor])#rho_surface * np.ones(coordinates.shape[0], dtype=float) * (side_length / (num_on_side / 2))**3 / 2
    


    internal_energies_surface = u_surface * np.ones(len(positions_surface), dtype=float)
    internal_energies_impactor = u_impactor * np.ones(len(positions_impactor), dtype=float)
    internal_energies = np.concatenate([internal_energies_surface, internal_energies_impactor])



    densities_surface = rho_surface * np.ones(len(positions_surface), dtype=float)
    densities_impactor = rho_impactor * np.ones(len(positions_impactor), dtype=float)
    densities = np.concatenate([densities_surface, densities_impactor])



    smoothinglengths = np.cbrt(masses / densities)
    
    
    
    mat_surface = mat_id_surface * np.ones(len(positions_surface), dtype=float)
    mat_impactor = mat_id_impactor * np.ones(len(positions_impactor), dtype=float)
    materials = np.concatenate([mat_surface, mat_impactor])

    boundary_y = 0.05 * duplications_y * side_length
    particle_ids_surface = np.empty(len(positions_surface))
    mask_boundary = positions_surface[:,1] <= np.min(positions_surface[:,1]) + boundary_y
    mask_boundary = np.logical_or(mask_boundary, positions_surface[:,0] <= np.min(positions_surface[:,0]) + boundary_y)
    mask_boundary = np.logical_or(mask_boundary, positions_surface[:,0] >= np.max(positions_surface[:,0]) - boundary_y)
    mask_boundary = np.logical_or(mask_boundary, positions_surface[:,2] <= np.min(positions_surface[:,2]) + boundary_y)
    mask_boundary = np.logical_or(mask_boundary, positions_surface[:,2] >= np.max(positions_surface[:,2]) - boundary_y)
    
    
    
    particle_ids_surface[mask_boundary] = np.linspace(1, len(particle_ids_surface[mask_boundary]), len(particle_ids_surface[mask_boundary]))
    particle_ids_surface[~mask_boundary] = np.linspace(len(particle_ids_surface[mask_boundary]) + 1, len(particle_ids_surface), len(particle_ids_surface[~mask_boundary]))
    
    particle_ids_impactor = np.linspace(len(positions_surface) + 1, len(positions_impactor)+len(positions_surface), len(positions_impactor))
    particle_ids = np.concatenate([particle_ids_surface, particle_ids_impactor])
    
    boundary_particles = len(particle_ids_surface[mask_boundary])
    print(
    "You need to compile the code with " "--enable-boundary-particles=%i" % boundary_particles
)
    
    numPart = coordinates.shape[0]
    print(numPart)
    pos = coordinates
    v = velocities
    m =  masses
    rho =  densities
    h = smoothinglengths
    u = internal_energies
    ids = particle_ids
    mat = materials
    file = h5py.File(filename, "w")

    # Header
    grp = file.create_group("/Header")
    grp.attrs["BoxSize"] = [(duplications_highres + 2 * duplications_lowres_eitherside) * side_length, 4 * duplications_y * side_length , (duplications_highres + 2 * duplications_lowres_eitherside) * side_length]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 3

    # Units
    grp = file.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 100.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1000.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = file.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=v, dtype="f")
    grp.create_dataset("Masses", data=m, dtype="f")
    grp.create_dataset("Density", data=rho, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=u, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
    grp.create_dataset("MaterialIDs", data=mat, dtype="i")

    return


if __name__ == "__main__":
    write_out_ics(filename="impact.hdf5", num_on_side_highres=32,)#32
