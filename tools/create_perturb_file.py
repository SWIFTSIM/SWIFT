#!/usr/bin/env python3

# Utility to create an HDF5 file with cosmological transfer functions using
# CLASS (python package classy)

import numpy as np
import h5py
from classy import Class

# Name of the perturbations file to be create
fname = "perturbations.hdf5"

# Open the file
f = h5py.File("perturbations.hdf5", mode="w")

# Cosmological parameters
h = 0.681
Omega_b = 0.0486
Omega_cdm = 0.2560110606
A_s = 2.0993736148e-09
n_s = 0.967

# Neutrino and radiation parameters
T_cmb = 2.7255
T_ncdm = 0.71611
N_ur = 2.0308
N_ncdm = 1
deg_ncdm = [1]
m_ncdm = [0.06]

# Maximum wavenumber and redshift
kmax = 30.0
zmax = 1e3
amin = 1.0 / (zmax + 1)

# CLASS output distance unit
Mpc_cgs = 3.085677581282e24

# CLASS parameters
params = {
    "h": h,
    "Omega_b": Omega_b,
    "Omega_cdm": Omega_cdm,
    "T_cmb": T_cmb,
    "N_ncdm": N_ncdm,
    "N_ur": N_ur,
    "T_ncdm": T_ncdm,
    "deg_ncdm": "".join(str(x) + "," for x in deg_ncdm)[:-1],
    "m_ncdm": "".join(str(x) + "," for x in m_ncdm)[:-1],
    "A_s": A_s,
    "n_s": n_s,
    "output": "dTk, vTk",
    "z_max_pk": zmax,
    "P_k_max_1/Mpc": kmax,
    "reio_parametrization": "reio_none",
    "YHe": "BBN",
    "k_output_values": kmax,
    "k_per_decade_for_pk": 100,
}

print("Running CLASS")

# Run CLASS
model = Class()
model.set(params)
model.compute()

# Extract wavenumbers and prepare redshifts
k = model.get_transfer(0)["k (h/Mpc)"] * h
a = np.exp(np.arange(0, np.log(amin), -0.01))[::-1]
z = 1.0 / a - 1.0
nk = len(k)
nz = len(z)

print("We have", nk, "wavenumbers and", nz, "redshifts")

keys = model.get_transfer(0).keys()

print("Available transfer functions:")
print(keys)

# Prepare dictionary
pt = {}
for key in keys:
    pt[key] = np.zeros((nz, nk))

# Extract transfer functions
for i in range(nz):
    pti = model.get_transfer(z[i])
    for key in pti:
        pt[key][i, :] = pti[key]

# Export the perturbations file
f.create_group("Functions")
f["Redshifts"] = z
f["Wavenumbers"] = k
f.create_group("Units")
f["Units"].attrs["Unit length in cgs (U_L)"] = Mpc_cgs

# Write the perturbations
for key in keys:
    f["Functions/" + key.replace("/", "\\")] = pt[key]

# Close the file
f.close()

print("Done.")
