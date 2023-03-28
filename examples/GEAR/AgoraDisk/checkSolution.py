import h5py
import numpy as np

expected_total_mass = 227413100.0

f = h5py.File("snapshot_0025.hdf5", "r")

# get the total stellar mass
mass = f["StarsParticles/Masses"]
total_mass = sum(mass) * 1e10

error = np.fabs((expected_total_mass - total_mass) / expected_total_mass)

if error < 0.01:
    print(48 * "!")
    print(r"Well done !")
    print(r"The stellar mass is %g Msol," % (total_mass))
    print(r"The expected stellar mass is %g Msol" % (expected_total_mass))
    print(r"This represents an error of %d %% !" % (100 * error))
    print(48 * "!")
else:
    print(48 * "!")
    print(r"Too bad !")
    print(r"The stellar mass is %g Msol," % (total_mass))
    print(r"While the expected stellar mass is %g Msol" % (expected_total_mass))
    print(r"This represents an error of %d %% !" % (100 * error))
    print(r"This is beyond the requirements.")
    print(48 * "!")
