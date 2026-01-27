import h5py as h5
import numpy as np

# Define some functions that will be used at the end to write the actual IC data
def write_header(f, boxsize, flag_entropy, np_total, np_total_hw, other=False):
    # We'll first build a dictionary to iterate through.

    default_attributes = {
        "BoxSize": boxsize,
        "Flag_Entropy_ICs": flag_entropy,
        "NumPart_Total": np_total,
        "NumPart_Total_HighWord": np_total_hw,
        "NumFilesPerSnapshot": 1,  # required for backwards compatibility
        "NumPart_ThisFile": np_total,  # Also required for bw compatibility
    }

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    header = f.create_group("Header")

    # For some reason there isn't a direct dictionary interface (at least one
    # that is documented, so we are stuck doing this loop...

    for name, value in attributes.items():
        header.attrs[name] = value

    return


def write_runtime_pars(f, periodic_boundary, other=False):
    # First build the dictionary

    default_attributes = {"PeriodicBoundariesOn": periodic_boundary}

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    runtime = f.create_group("RuntimePars")

    for name, value in attributes.items():
        runtime.attrs[name] = value

    return


def write_units(f, current, length, mass, temperature, time, other=False):
    # First build the dictionary

    default_attributes = {
        "Unit current in cgs (U_I)": current,
        "Unit length in cgs (U_L)": length,
        "Unit mass in cgs (U_M)": mass,
        "Unit temperature in cgs (U_T)": temperature,
        "Unit time in cgs (U_t)": time,
    }

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    units = f.create_group("Units")

    for name, value in attributes.items():
        units.attrs[name] = value

    return


def write_block(f, part_type, pos, vel, ids, mass, int_energy, smoothing, other=False):

    # Build the dictionary

    if part_type == 0:
        default_data = {
            "Coordinates": pos,
            "Velocities": vel,
            "ParticleIDs": ids,
            "Masses": mass,
            "InternalEnergy": int_energy,
            "SmoothingLength": smoothing,
        }
    elif part_type == 4:
        default_data = {
            "Coordinates": pos,
            "Velocities": vel,
            "ParticleIDs": ids,
            "Masses": mass,
            "SmoothingLength": smoothing,
        }
    else:
        default_data = {
            "Coordinates": pos,
            "Velocities": vel,
            "ParticleIDs": ids,
            "Masses": mass,
            "SmoothingLength": smoothing,
        }

    if other:
        data = dict(default_data, **other)
    else:
        data = default_data

    particles = f.create_group("PartType" + str(part_type))

    for name, value in data.items():
        particles.create_dataset(name, data=value)

    return


### DEFINE CONSTANTS ###

# Box size and shape in internal units. We assume box shape of [box_aspect_ratio, box_aspect_ratio, 1] * box_size
box_size = 500.0  # Make sure this is consistent with parameter file
box_aspect_ratio = 0.5  # Make sure this is consistent with parameter file

# Particle mass (of the ambient medium and jet particles, separately) in internal units
m_part = 1.81e-4  # EAGLE res
m_part_jet = 1.81e-4  # EAGLE res
print("Resolution in the ambient medium:\n", m_part)
print("Resolution in the jet:\n", m_part_jet)

# Density of the ambient medium in internal units.
ambient_medium_density = (
    1.47e-5
)  # Corresponds to 0.01 cm^-3, typical of galaxy cluster centers (of halo mass 10**15 Msun)

# Temperature of the ambient medium in Kelvin.
ambient_medium_temperature = (
    10 ** 7
)  # Typical of cool core galaxy clusters (of halo mass 10**15 Msun)

## Some jet-related parameters ##
# Radius of the jet launching region. This should be a small fraction of the box size (and simulated jets/lobes), but not too small, as it can slow down the simulation.
R = 5.0  # kpc

# (Half-) opening angle of the jet
opening_angle = 10.0  # in degrees
cos_opening_angle = np.cos(opening_angle / 180.0 * np.pi)

# We need to compute how many particles we want to launch. Assume some jet parameters for this
# Jet power, duration that it's on, and total launched energy
P_jet = (
    1.56e8
)  # 10^45 erg/s in internal units. Make sure it is consistent with parameter file
T_jet = (
    0.10228
)  # 100 Myr in internal units. Make sure it is consistent with parameter file
E_jet = P_jet * T_jet

# Jet velocity
v_jet = 7.5e3  # km/s in internal units. Make sure it is consistent with parameter file

# Finally we can compute total number of particles to launch
N_jet_target = int(np.floor(E_jet / (0.5 * m_part_jet * v_jet ** 2)))
print("Target number of particles to launch:\n", N_jet_target)


### CREATE AMBIENT MEDIUM ###

# The number of particles we need to create. Obtained simply as N = density / V
N_ambient_medium = int(
    np.floor(
        ambient_medium_density * box_size * (box_size * box_aspect_ratio) ** 2 / m_part
    )
)
print("Number of particles in the ambient medium:\n", N_ambient_medium)

# Create a rectangular lattice from that number. Begin by defining the longest side of the rectangle
N_z = int(
    np.floor(N_ambient_medium ** (1.0 / 3.0) * 1 / box_aspect_ratio ** (2.0 / 3.0))
)  # This isn't a simple cube root of N_ambient_medium because we have a rectangular box

# Shorter two sides
N_xy = int(np.floor(N_z * box_aspect_ratio))

# Create the actual lattice
x_grid = np.linspace(
    -(1.0 - 1.0 / (1.0 + N_xy)) * box_size * 0.5 * box_aspect_ratio,
    (1.0 - 1.0 / (1.0 + N_xy)) * box_size * 0.5 * box_aspect_ratio,
    N_xy,
)
y_grid = np.linspace(
    -(1.0 - 1.0 / (1.0 + N_xy)) * box_size * 0.5 * box_aspect_ratio,
    (1.0 - 1.0 / (1.0 + N_xy)) * box_size * 0.5 * box_aspect_ratio,
    N_xy,
)
z_grid = np.linspace(
    -(1.0 - 1.0 / (1.0 + N_z)) * box_size * 0.5,
    (1.0 - 1.0 / (1.0 + N_z)) * box_size * 0.5,
    N_z,
)
matrix = np.meshgrid(x_grid, y_grid, z_grid)

# Move the particles to occupy rectangular box rather than be centred around [0, 0, 0]
x_ambient_medium = np.ndarray.flatten(matrix[0]) + box_aspect_ratio * box_size * 0.5
y_ambient_medium = np.ndarray.flatten(matrix[1]) + box_aspect_ratio * box_size * 0.5
z_ambient_medium = np.ndarray.flatten(matrix[2]) + box_size * 0.5

# Final ambient medium particle number
N_ambient_medium = np.size(x_ambient_medium)


### CREATE JET RESERVOIR ###
# We first create a cubical lattice and then cut out cones

# Number of particles that we need to occupy a cubical lattice such that once we cut out cones, the remaining number is roughly the number we want to launch.
# This won't give us the exact number, but is rough. The numerical multiplier is extra padding to ensure we get at least the number we want to launch
N_cone_box = 6 / np.pi * 1.0 / (1.0 - cos_opening_angle) * N_jet_target * 1.1

# Side length of the necessary cube
N_side = int(np.floor(N_cone_box ** (1.0 / 3.0)))

# Create the lattice
x_grid = np.linspace(-R, R, N_side)
y_grid = np.linspace(-R, R, N_side)
z_grid = np.linspace(-R, R, N_side)
matrix = np.meshgrid(x_grid, y_grid, z_grid)

# Jet particle positions
x_pos_jet = np.ndarray.flatten(matrix[0])
y_pos_jet = np.ndarray.flatten(matrix[1])
z_pos_jet = np.ndarray.flatten(matrix[2])

# Jet particle radii from origin
r_pos_jet = np.sqrt((x_pos_jet) ** 2 + (y_pos_jet) ** 2 + (z_pos_jet) ** 2)

# Cosines relative to z-axis (the launching axis) of all cubical lattice particles
cos_theta = np.absolute(z_pos_jet) / r_pos_jet

# Now select only the cone particles
selection_cone = (r_pos_jet < 1.0001 * R) & (cos_theta > cos_opening_angle)
x_pos_jet = x_pos_jet[selection_cone]
y_pos_jet = y_pos_jet[selection_cone]
z_pos_jet = z_pos_jet[selection_cone]
r_pos_jet = r_pos_jet[selection_cone]

# Number of jet particles in the cone we've ended up with
N_jet_final = np.size(x_pos_jet)
print("Final particle number in the jet cone reservoir:\n", N_jet_final)

# Check if we have enough
if N_jet_final < N_jet_target:
    print(
        "The final particle number in the cone reservoir is less than the target number. Increase the numerical prefactor in the N_cone_box variabe."
    )
    exit()

# Move particles from origin
x_pos_jet += box_aspect_ratio * 0.5 * box_size
y_pos_jet += box_aspect_ratio * 0.5 * box_size
z_pos_jet += 0.5 * box_size

# Sort particles so that first one in the array is the farthest one, and so forth.
sort = np.argsort(r_pos_jet)
x_pos_jet = x_pos_jet[sort][::-1]
y_pos_jet = y_pos_jet[sort][::-1]
z_pos_jet = z_pos_jet[sort][::-1]


### DO SOME FINAL PREPARATIONS FOR WRITING ###

# Total particle number
N_total = N_jet_final + N_ambient_medium

# Final gas positions
x_gas = np.concatenate((x_pos_jet, x_ambient_medium), axis=None)
y_gas = np.concatenate((y_pos_jet, y_ambient_medium), axis=None)
z_gas = np.concatenate((z_pos_jet, z_ambient_medium), axis=None)

# Final gas velocities
v_x_gas = np.zeros(N_total)
v_y_gas = np.zeros(N_total)
v_z_gas = np.zeros(N_total)

# Transpose the coordinate and velocity arrays for final writing
positions_gas = np.array([x_gas, y_gas, z_gas]).T
velocities_gas = np.array([v_x_gas, v_y_gas, v_z_gas]).T

# Compute internal energy instead of temperature, apply to all particles, not just ambient medium
internal_energy = 0.0207 * ambient_medium_temperature  # per unit mass (internal unit)
internal_energies_gas = internal_energy * np.ones(N_total)

# Assign particle IDs
IDs_gas = np.concatenate(
    (np.arange(1,N_jet_final+1), 10 ** 7 * 2 + np.arange(N_ambient_medium)), axis=None
)

# Assign masses
masses_gas = np.concatenate(
    (m_part_jet * np.ones(N_jet_final), m_part * np.ones(N_ambient_medium)), axis=None
)

# Calculate smoothing lengths, from the density of the ambient medium and jet cones. This is approximate! SWIFT recalculates this anyway.
smoothing_length_ambient_medium = 0.5 * (50 * m_part / (ambient_medium_density)) ** (
    1.0 / 3.0
)
volume_jet_cone = 4.0 * np.pi / 3.0 * R ** 3 * np.sin(opening_angle) ** 2
smoothing_length_jet_cone = 0.5 * (50 / (N_jet_final / volume_jet_cone)) ** (1.0 / 3.0)
smoothing_lengths_gas = np.concatenate(
    (
        smoothing_length_jet_cone * np.ones(N_jet_final),
        smoothing_length_ambient_medium * np.ones(N_ambient_medium),
    ),
    axis=None,
)


### WRITE ICs ####
print("Writing...")

# Write ICs
with h5.File("ICs.hdf5", "w") as f:
    write_header(
        f,
        boxsize=[box_size * box_aspect_ratio, box_size * box_aspect_ratio, box_size],
        flag_entropy=0,
        np_total=[N_total, 0, 0, 0, 0, 0],
        np_total_hw=[0, 0, 0, 0, 0, 0],
    )

    write_runtime_pars(f, periodic_boundary=1)

    write_units(
        f,
        current=1.0,
        length=3.08566 * 10 ** 21,  # kpc
        mass=1.98848 * 10 ** 43.0,  # 10**10 Msun
        temperature=1.0,  # Kelvin
        time=3.156 * 10 ** 16.0,  # Gyr
    )

    write_block(
        f,
        part_type=0,
        pos=positions_gas,
        vel=velocities_gas,
        ids=IDs_gas,
        mass=masses_gas,
        int_energy=internal_energies_gas,
        smoothing=smoothing_lengths_gas,
    )

print("Done!")
