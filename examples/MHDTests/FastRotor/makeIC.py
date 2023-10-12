import h5py
import numpy as np
import matplotlib.pyplot as plt

###---------------------------###

# Parameters

r_in = 0.1
rho_in_0 = 10.0
rho_out_0 = 1.0
P_0 = 1.0
B_0 = 0.0  # 2.5 / np.sqrt(np.pi)
omega_0 = 20.0

fileOutputName = "FastRotor_no_B.hdf5"

###---------------------------###

N_side = 8
N_face = N_side ** 2
N_unit_cell = N_side ** 3

unit_cell = np.zeros((N_unit_cell, 3))
unit_cell_smoothing_lengths = np.ones(N_unit_cell)

for i in range(0, N_unit_cell):
    unit_cell[i, 0] = (i % N_side) / N_side
    # print(unit_cell[i,0])
    unit_cell[i, 1] = ((i // N_side) % N_side) / N_side
    unit_cell[i, 2] = ((i // N_face) % N_side) / N_side

# plt.scatter(unit_cell[:,0],unit_cell[:,1])
# plt.savefig("test.png")

# plt.scatter(unit_cell[:,0],unit_cell[:,2])
# plt.savefig("test2.png")

# plt.scatter(unit_cell[:,1],unit_cell[:,2])
# plt.savefig("test3.png")

N_UC_reps = 10
cx = N_UC_reps
cy = N_UC_reps
cz = 1

lx = 1.0
ly = 1.0
lz = 1.0 / float(N_UC_reps)
vol = lx * ly * lz

pos = np.zeros((int(cx * cy * cz * N_unit_cell), 3))
h = np.zeros(int(cx * cy * cz * N_unit_cell))

N_particles = N_unit_cell * cx * cy * cz

for i in range(0, cx):
    for j in range(0, cy):
        for k in range(0, cz):
            ind_low = i * N_unit_cell
            ind_upper = ind_low + N_unit_cell

            pos[ind_low:ind_upper, 0] = unit_cell[:, 0] + i

            pos[ind_low:ind_upper, 1] = unit_cell[:, 1] + j

            pos[ind_low:ind_upper, 2] = unit_cell[:, 2] + k

            h[ind_low:ind_upper] = unit_cell_smoothing_lengths

pos[:, 0] = pos[:, 0] * lx / cx
pos[:, 1] = pos[:, 1] * ly / cy
pos[:, 2] = pos[:, 2] * lz / cz


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c="b")

plt.savefig("test.png")

quit()

###---------------------------###

# Generate extra arrays

v = np.zeros((N, 3))
B = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))
