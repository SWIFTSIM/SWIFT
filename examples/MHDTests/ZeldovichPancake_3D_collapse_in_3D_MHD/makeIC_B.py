### This file is part of SWIFT.
import h5py
from numpy import *

# Generates a swift IC file for the 3D collapse in the Zeldovich pancake model, with magnetic field

# Parameters
B_0 = 1.0e-9 # Initial magnetic field strength (in T)
T_i = 100.0  # Initial temperature of the gas (in K)
# Redshift of caustic formation, which can have 3 different values for each direction now as well
z_c_x = 5.0
z_c_y = 1.0
z_c_z = 0.5
z_i = 100.0  # Initial redshift
gamma = 5.0 / 3.0  # Gas adiabatic index
numPart_1D = 32  # Number of particles along each dimension
fileName = "zeldovichPancake.hdf5"

# Some units
Mpc_in_m = 3.08567758e22
Msol_in_kg = 1.98848e30
Gyr_in_s = 3.08567758e19  ### so this is not actually a Gigayear, but 1 Mpc / (1 km/s) !!!
mH_in_kg = 1.6737236e-27

# Some constants
kB_in_SI = 1.38064852e-23
G_in_SI = 6.67408e-11

# Some useful variables in h-full units
H_0 = 1.0 / Mpc_in_m * 10 ** 5  # h s^-1
rho_0 = 3.0 * H_0 ** 2 / (8 * math.pi * G_in_SI)  # h^2 kg m^-3
lambda_i = 64.0 / H_0 * 10 ** 5  # h^-1 m (= 64 h^-1 Mpc)
x_min = -0.5 * lambda_i
x_max = 0.5 * lambda_i

# SI system of units
unit_l_in_si = Mpc_in_m
unit_m_in_si = Msol_in_kg * 1.0e10
unit_t_in_si = Gyr_in_s
unit_v_in_si = unit_l_in_si / unit_t_in_si
unit_u_in_si = unit_v_in_si ** 2
unit_I_in_si = 1
unit_B_in_si = unit_m_in_si * unit_t_in_si**(-2) * unit_I_in_si**(-1)

# Total number of particles
numPart = numPart_1D ** 3

# ---------------------------------------------------

# Get the frequency of the initial perturbation
## We now set frequencies for all three directions; keep them the same for now
k_x = 2.0 * pi / lambda_i
k_y = 2.0 * pi / lambda_i
k_z = 2.0 * pi / lambda_i

# Get the redshift prefactor for the initial positions
zfac_x = (1.0 + z_c_x) / (1.0 + z_i)
zfac_y = (1.0 + z_c_y) / (1.0 + z_i)
zfac_z = (1.0 + z_c_z) / (1.0 + z_i)

# Set box size and interparticle distance
boxSize = x_max - x_min
delta_x = boxSize / numPart_1D  # = delta_y = delta_z

# Get the particle mass
a_i = 1.0 / (1.0 + z_i)
m_i = boxSize ** 3 * rho_0 / numPart

# Build the arrays
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
h = zeros(numPart)
u = zeros(numPart)

# Define the magnetic field
B = zeros((numPart, 3))
Bx = B_0
By = 0
Bz = 0
B[:, 0] = Bx
B[:, 1] = By
B[:, 2] = Bz

# Set the particles on the left
for i in range(numPart_1D):
    for j in range(numPart_1D):
        for k in range(numPart_1D):
            index = i * numPart_1D ** 2 + j * numPart_1D + k

            # in x
            q_x = x_min + (i + 0.5) * delta_x
            coords[index, 0] = q_x - zfac_x * sin(k_x * q_x) / k_x - x_min

            # in y
            q_y = x_min + (j + 0.5) * delta_x
            coords[index, 1] = q_y - zfac_y * sin(k_y * q_y) / k_y - x_min

            # in z
            q_z = x_min + (k + 0.5) * delta_x
            coords[index, 2] = q_z - zfac_z * sin(k_z * q_z) / k_z - x_min

            # The rest, T and v are also changed to 3D
            T = T_i * (1.0 / (1.0 - zfac_x * cos(k_x * q_x)) / (1.0 - zfac_y * cos(k_y * q_y)) / (1.0 - zfac_z * cos(k_z * q_z)) ) ** (2.0 / 3.0)
            u[index] = kB_in_SI * T / (gamma - 1.0) / mH_in_kg
            h[index] = 1.2348 * delta_x
            m[index] = m_i
            v[index, 0] = -H_0 * (1.0 + z_c_x) / sqrt(1.0 + z_i) * sin(k_x * q_x) / k_x
            v[index, 1] = -H_0 * (1.0 + z_c_y) / sqrt(1.0 + z_i) * sin(k_y * q_y) / k_y
            v[index, 2] = -H_0 * (1.0 + z_c_z) / sqrt(1.0 + z_i) * sin(k_z * q_z) / k_z

# Unit conversion
coords /= unit_l_in_si
v /= unit_v_in_si
m /= unit_m_in_si
h /= unit_l_in_si
u /= unit_u_in_si
B /= unit_B_in_si

boxSize /= unit_l_in_si

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize]
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
grp.attrs["Unit length in cgs (U_L)"] = 100.0 * unit_l_in_si
grp.attrs["Unit mass in cgs (U_M)"] = 1000.0 * unit_m_in_si
grp.attrs["Unit time in cgs (U_t)"] = 1.0 * unit_t_in_si
grp.attrs["Unit current in cgs (U_I)"] = 1.0 * unit_I_in_si
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=coords, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
