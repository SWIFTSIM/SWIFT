import numpy as np
import h5py
import sys

def generate_cubic_lattice(particle_sep, side_x, side_y, side_z):

    N_x = int(side_x / particle_sep)
    N_y = int(side_y / particle_sep)
    N_z = int(side_z / particle_sep)

    A2_positions = np.zeros((N_x * N_y * N_z, 3))

    for i in range(N_x):
        for j in range(N_y):
            for k in range(N_z):
                index = i * N_y * N_z + j * N_z + k

                x = i / float(N_x) + 1.0 / (2.0 * N_x)
                y = j / float(N_y) + 1.0 / (2.0 * N_y)
                z = k / float(N_z) + 1.0 / (2.0 * N_z)

                A2_positions[index, 0] = x * side_x
                A2_positions[index, 1] = y * side_y
                A2_positions[index, 2] = z * side_z

    # Adjust lattice so that 0 is in centre of x, y and at top of z
    # This corresponds to the point of impact
    A2_positions[:, 0] -= np.mean(A2_positions[:, 0])
    A2_positions[:, 1] -= np.mean(A2_positions[:, 1])
    A2_positions[:, 2] -= np.max(A2_positions[:, 2]) + 0.5 * particle_sep

    return A2_positions



def generate_highres(
    approx_highres_width,
    approx_highres_depth,
    particle_sep_highres,
    particle_sep_lowres,
    density
):

    # Width and depth of high resolution region. Must be multiple of 2*particle_sep_lowres
    # Add an extra 2*particle_sep_lowres so that it "rounds" up
    highres_width = approx_highres_width - approx_highres_width%(2*particle_sep_lowres) + 2*particle_sep_lowres
    highres_depth = approx_highres_depth - approx_highres_depth%(2*particle_sep_lowres) + 2*particle_sep_lowres
    
    A2_pos_highres = generate_cubic_lattice(particle_sep_highres, highres_width,
                                                highres_width, highres_depth)
    N_particles = len(A2_pos_highres)
    particle_volume = (highres_width * highres_width * highres_depth) / N_particles
    A1_m_highres = (density * particle_volume) * np.ones(N_particles)

    return A2_pos_highres, A1_m_highres, highres_width, highres_depth


def generate_lowres(
    approx_surface_width,
    approx_surface_depth,
    highres_width,
    highres_depth,
    transition,
    particle_sep_lowres,
    particle_sep_highres,
    density
):

    # Calculate the width of the surface. Must be whole number of partiles
    # Approx length of low resolution region on one side of the centre
    approx_lowres_one_side = 0.5 * (approx_surface_width - highres_width -
                                                                2 * transition)
    # Make multiple of particle_sep_lowres and add an extra particle_sep_lowres so that it "rounds" up
    lowres_one_side = approx_lowres_one_side - approx_lowres_one_side%particle_sep_lowres + particle_sep_lowres
    surface_width = 2 * lowres_one_side + 2 * transition + highres_width

    # Calculate the depth of the surface. Must be whole number of partiles
    # Approx depth of low resolution region
    approx_lowres_depth = approx_surface_depth - highres_depth - transition
    # Make multiple of particle_sep_lowres and add an extra particle_sep_lowres so that it "rounds" up
    lowres_depth = approx_lowres_depth - approx_lowres_depth%particle_sep_lowres + particle_sep_lowres
    surface_depth = lowres_depth + transition + highres_depth

    # Generate full surace made up only of a low resolution region
    A2_pos_surface = generate_cubic_lattice(particle_sep_lowres,
                                surface_width, surface_width, surface_depth)
    N_particles = len(A2_pos_surface)
    particle_volume = (surface_width * surface_width * surface_depth) / N_particles
    A1_m_surface = (density * particle_volume) * np.ones(N_particles)

    # Make masks to select the low resolution region
    mask_lowres_x = np.logical_or(
                             A2_pos_surface[:, 0] < -(0.5 * highres_width + transition),
                             A2_pos_surface[:, 0] > (0.5 * highres_width + transition),
                            )
    mask_lowres_y = np.logical_or(
                             A2_pos_surface[:, 1] < -(0.5 * highres_width + transition),
                             A2_pos_surface[:, 1] > (0.5 * highres_width + transition),
                            )
    mask_lowres_z = A2_pos_surface[:, 2] < -(highres_depth + transition)

    mask_lowres = np.logical_or(mask_lowres_x, mask_lowres_y)
    mask_lowres = np.logical_or(mask_lowres, mask_lowres_z)

    # Select low region from full surface
    A2_pos_lowres = A2_pos_surface[mask_lowres]
    A1_m_lowres = A1_m_surface[mask_lowres]

    """
    # Slightly adjust z positions since we're shifting grid up to align with surface
    mask = A2_pos_lowres[:, 2] < -highres_depth-transition
    first_layer_below = np.max(A2_pos_lowres[mask, 2])
    adjustment_needed = (-highres_depth - transition - 0.25 * (particle_sep_lowres + particle_sep_highres)) / first_layer_below
    A2_pos_lowres[:, 2] *= adjustment_needed
    surface_depth *= adjustment_needed
    """
    return A2_pos_lowres, A1_m_lowres, surface_width, surface_depth




def generate_transition(
    transition,
    transition_regions,
    highres_width,
    highres_depth,
    particle_sep_highres,
    density
):


    A2_pos_transition = np.empty((0,3))
    A1_m_transition = np.empty(0)


    last_transition_distance = 0
    for i in range(transition_regions):

        if i > 0:
            # For most regions, halve the size
            current_transition_distance = last_transition_distance + transition / 2**(transition_regions - i)
        else:
            # For last region, add an extra region of the same size
            current_transition_distance = last_transition_distance + 2 * transition / 2**(transition_regions - i)
        current_particle_sep =  2**(i + 1) * particle_sep_highres

        current_width = highres_width + 2 * current_transition_distance
        current_depth = highres_depth + current_transition_distance

        A2_pos_current = generate_cubic_lattice(current_particle_sep, current_width,
                                                    current_width, current_depth)

        N_particles = len(A2_pos_current)
        particle_volume = (current_width * current_width * current_depth) / N_particles
        A1_m_current = (density * particle_volume) * np.ones(N_particles)

        mask_current_x = np.logical_or(
                                    A2_pos_current[:, 0] < -(0.5 * highres_width + last_transition_distance),
                                    A2_pos_current[:, 0] > (0.5 * highres_width + last_transition_distance),
                                )
        mask_current_y = np.logical_or(
                                 A2_pos_current[:, 1] < -(0.5 * highres_width + last_transition_distance),
                                 A2_pos_current[:, 1] > (0.5 * highres_width + last_transition_distance),
                                )
        mask_current_z = A2_pos_current[:, 2] < -(highres_depth + last_transition_distance)

        mask_current = np.logical_or(mask_current_x, mask_current_y)
        mask_current = np.logical_or(mask_current, mask_current_z)

        A2_pos_transition = np.append(A2_pos_transition, A2_pos_current[mask_current], axis=0)
        A1_m_transition = np.append(A1_m_transition, A1_m_current[mask_current])

        last_transition_distance = current_transition_distance

    return A2_pos_transition, A1_m_transition

def generate_surface(
    particle_sep_highres,
    particle_sep_lowres,
    approx_surface_width,
    approx_surface_depth,
    approx_highres_width,
    approx_highres_depth,
    transition,
    mat_id_surface,
    density_surface,
    u_surface,
):

        (
            A2_pos_highres,
            A1_m_highres,
            highres_width,
            highres_depth
        ) = generate_highres(
            approx_highres_width,
            approx_highres_depth,
            particle_sep_highres,
            particle_sep_lowres,
            density_surface
        )

        # check whether a transition region is needed
        lowres_power_of_2 = int(np.log2(particle_sep_lowres / particle_sep_highres))
        transition_regions = lowres_power_of_2 - 1

        if transition_regions > 0:
            (
                A2_pos_lowres,
                A1_m_lowres,
                surface_width,
                surface_depth
            ) = generate_lowres(
                approx_surface_width,
                approx_surface_depth,
                highres_width,
                highres_depth,
                transition,
                particle_sep_lowres,
                particle_sep_highres,
                density_surface
            )

            (
                A2_pos_transition,
                A1_m_transition
            ) = generate_transition(
                transition,
                transition_regions,
                highres_width,
                highres_depth,
                particle_sep_highres,
                density_surface
            )

            A2_pos_surface = np.concatenate([A2_pos_highres, A2_pos_transition])
            A2_pos_surface = np.concatenate([A2_pos_surface, A2_pos_lowres])

            A1_m_surface = np.concatenate([A1_m_highres, A1_m_transition])
            A1_m_surface = np.concatenate([A1_m_surface, A1_m_lowres])

        else:
            (
                A2_pos_lowres,
                A1_m_lowres,
                surface_width,
                surface_depth
            ) = generate_lowres(
                approx_surface_width,
                approx_surface_depth,
                highres_width,
                highres_depth,
                0,
                particle_sep_lowres,
                particle_sep_highres,
                density_surface
            )

            A2_pos_surface = np.concatenate([A2_pos_highres, A2_pos_lowres])
            A1_m_surface = np.concatenate([A1_m_highres, A1_m_lowres])

        numPart_surface = len(A1_m_surface)
        A2_vel_surface = np.zeros_like(A2_pos_surface)
        A1_rho_surface = density_surface * np.ones(numPart_surface)
        A1_u_surface = u_surface * np.ones(numPart_surface)
        A1_mat_id_surface = mat_id_surface * np.ones(numPart_surface)

        return (
            A2_pos_surface,
            A2_vel_surface,
            A1_m_surface,
            A1_rho_surface,
            A1_u_surface,
            A1_mat_id_surface,
            surface_width,
            surface_depth,
            highres_width,
            highres_depth,
        )

def generate_impactor(
    impactor_radius,
    impactor_angle,
    impactor_speed,
    particle_sep_highres,
    mat_id_impactor,
    density_impactor,
    u_impactor,
):

    A2_pos_cube = generate_cubic_lattice(particle_sep_highres, 3 * impactor_radius,
                                                3 * impactor_radius, 3 * impactor_radius)
    A2_pos_cube[:, 2] += 1.5 * impactor_radius
    N_particles = len(A2_pos_cube)
    particle_volume = (3 * impactor_radius)**3 / N_particles
    A1_m_cube = (density_impactor * particle_volume) * np.ones(N_particles)


    # Now select particles in impactor
    A1_r_cube = np.sqrt(A2_pos_cube[:, 0]**2 + A2_pos_cube[:, 1]**2 + A2_pos_cube[:, 2]**2)
    mask_sphere = A1_r_cube <= impactor_radius
    A2_pos_impactor = A2_pos_cube[mask_sphere]
    A1_m_impactor = A1_m_cube[mask_sphere]

    # Move above surface
    A2_pos_impactor[:, 0] += 2 * impactor_radius  * np.cos(np.deg2rad(impactor_angle))
    A2_pos_impactor[:, 2] += 2 * impactor_radius  * np.sin(np.deg2rad(impactor_angle))

    # Generate velocities
    A2_vel_impactor = np.zeros_like(A2_pos_impactor)
    A2_vel_impactor[:,0] = -impactor_speed * np.cos(np.deg2rad(impactor_angle))
    A2_vel_impactor[:,2] = -impactor_speed * np.sin(np.deg2rad(impactor_angle))

    numPart_impactor = len(A1_m_impactor)
    A1_rho_impactor = density_impactor * np.ones(numPart_impactor)
    A1_u_impactor = u_impactor * np.ones(numPart_impactor)
    A1_mat_id_impactor = mat_id_impactor * np.ones(numPart_impactor)

    return (
        A2_pos_impactor,
        A2_vel_impactor,
        A1_m_impactor,
        A1_rho_impactor,
        A1_u_impactor,
        A1_mat_id_impactor,
    )


def write_out_ics(
    filename,
    impactor_radius,
    impactor_angle,
    impactor_speed,
    particle_sep_highres,
    particle_sep_lowres,
    approx_surface_width,
    approx_surface_depth,
    approx_highres_width,
    approx_highres_depth,
    transition,
    mat_id_surface,
    density_surface,
    u_surface,
    mat_id_impactor,
    density_impactor,
    u_impactor,
):

    (
        A2_pos_surface,
        A2_vel_surface,
        A1_m_surface,
        A1_rho_surface,
        A1_u_surface,
        A1_mat_id_surface,
        surface_width,
        surface_depth,
        highres_width,
        highres_depth,
    ) = generate_surface(
        particle_sep_highres,
        particle_sep_lowres,
        approx_surface_width,
        approx_surface_depth,
        approx_highres_width,
        approx_highres_depth,
        transition,
        mat_id_surface,
        density_surface,
        u_surface
    )

    (
        A2_pos_impactor,
        A2_vel_impactor,
        A1_m_impactor,
        A1_rho_impactor,
        A1_u_impactor,
        A1_mat_id_impactor,
    ) = generate_impactor(
        impactor_radius,
        impactor_angle,
        impactor_speed,
        particle_sep_highres,
        mat_id_impactor,
        density_impactor,
        u_impactor,
    )
    
    # impactor offset to crudely centre crater in high resolution region
    offset = 2 * impactor_radius * np.cos(np.deg2rad(impactor_angle))
    A2_pos_impactor[:, 0] += offset

    coordinates = np.concatenate([A2_pos_surface, A2_pos_impactor])
    coordinates[:, 0] += 0.5 * surface_width
    coordinates[:, 1] += 0.5 * surface_width
    coordinates[:, 2] += surface_depth

    
    velocities = np.concatenate([A2_vel_surface, A2_vel_impactor])
    masses = np.concatenate([A1_m_surface, A1_m_impactor])
    densities = np.concatenate([A1_rho_surface, A1_rho_impactor])
    smoothinglengths = np.cbrt(masses / densities)
    internal_energies = np.concatenate([A1_u_surface, A1_u_impactor])
    materials = np.concatenate([A1_mat_id_surface, A1_mat_id_impactor])

    # Set particle IDs so that boundary particles have the lowest IDs
    boundary = 0.05 * surface_depth
    A1_ids_surface = np.empty(len(A2_pos_surface))
    # Masks for boundary particles
    mask_boundary = A2_pos_surface[:,2] <= np.min(A2_pos_surface[:,2]) + boundary
    mask_boundary = np.logical_or(mask_boundary, A2_pos_surface[:,0] <= np.min(A2_pos_surface[:,0]) + boundary)
    mask_boundary = np.logical_or(mask_boundary, A2_pos_surface[:,0] >= np.max(A2_pos_surface[:,0]) - boundary)
    mask_boundary = np.logical_or(mask_boundary, A2_pos_surface[:,1] <= np.min(A2_pos_surface[:,1]) + boundary)
    mask_boundary = np.logical_or(mask_boundary, A2_pos_surface[:,1] >= np.max(A2_pos_surface[:,1]) - boundary)

    # Set IDs
    A1_ids_surface[mask_boundary] = np.linspace(1, len(A1_ids_surface[mask_boundary]), len(A1_ids_surface[mask_boundary]))
    A1_ids_surface[~mask_boundary] = np.linspace(len(A1_ids_surface[mask_boundary]) + 1, len(A1_ids_surface), len(A1_ids_surface[~mask_boundary]))
    A1_ids_impactor = np.linspace(len(A2_pos_surface) + 1, len(A2_pos_impactor)+len(A2_pos_surface), len(A2_pos_impactor))
    particle_ids = np.concatenate([A1_ids_surface, A1_ids_impactor])

    boundary_particles = len(A1_ids_surface[mask_boundary])
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
    grp.attrs["BoxSize"] = [surface_width, surface_width, 2 * surface_depth]
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

    # Parameters of impact
    impactor_radius = 0.5e-3
    pppr = 10
    impactor_angle = 50
    impactor_speed = 5e3
    lowres_power_of_2 = 3

    # Particle separation in high resolution region
    particle_sep_highres = impactor_radius / pppr
    # Particle separation in low resolution region
    particle_sep_lowres = 2**lowres_power_of_2 * particle_sep_highres


    # Surface dimensions
    approx_surface_width = 40e-3
    approx_surface_depth = 20e-3
    approx_highres_width = 10e-3
    approx_highres_depth = 2e-3
    transition = 4 * particle_sep_lowres#2e-3

    # Surface and impactor properties
    mat_id_surface = 1000
    density_surface = 2700
    u_surface = 1e5
    mat_id_impactor = 1000
    density_impactor = 2700
    u_impactor = 1e5

    write_out_ics(filename="impact.hdf5",
                  impactor_radius = impactor_radius,
                  impactor_angle = impactor_angle,
                  impactor_speed = impactor_speed,
                  particle_sep_highres = particle_sep_highres,
                  particle_sep_lowres = particle_sep_lowres,
                  approx_surface_width = approx_surface_width,
                  approx_surface_depth = approx_surface_depth,
                  approx_highres_width = approx_highres_width,
                  approx_highres_depth = approx_highres_depth,
                  transition = transition,
                  mat_id_surface = mat_id_surface,
                  density_surface = density_surface,
                  u_surface = u_surface,
                  mat_id_impactor = mat_id_impactor,
                  density_impactor = density_impactor,
                  u_impactor = u_impactor,
    )
