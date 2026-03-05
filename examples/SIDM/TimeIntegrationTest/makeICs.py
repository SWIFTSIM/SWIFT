import h5py
import numpy as np

num_sidm = 100
nparttype = 8
L = 10

# define standard units
UnitMass_in_cgs = 1.988409870698051e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.0856775814913673e21  # kpc in centimeters
UnitVelocity_in_cgs = 1e5  # km/s in centimeters per second
UnitCurrent_in_cgs = 1  # Amperes
UnitTemp_in_cgs = 1  # Kelvin
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs
UnitLength = UnitLength_in_cgs
UnitTime = UnitTime_in_cgs 
UnitVelocity = UnitVelocity_in_cgs

with h5py.File("sidm_ics.hdf5", "w") as f:
    # Header attributes
    f.create_group("Header")
    f["Header"].attrs["NumPart_ThisFile"] = np.array([0, 0,0, 0, 0, 0, 0, num_sidm], dtype=np.uint32)
    f["Header"].attrs["NumPart_Total"] = np.array([0, 0, 0, 0, 0, 0, 0, num_sidm], dtype=np.uint32)
    f["Header"].attrs["NumPart_Total_HighWord"] = np.zeros(nparttype, dtype=np.uint32)
    f["Header"].attrs["BoxSize"] = [L, L, L]  # example box size
    f["Header"].attrs["Time"] = 0.0
    f["Header"].attrs["NumFileOutputsPerSnapshot"] = 1
    f["Header"].attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    f["Header"].attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0, 0, 0]
    f["Header"].attrs["Dimension"] = 3

    # Units
    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
    grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
    grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
    grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
    grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

    # Create groups for each particle type
    # SIDM particles in PartType7
    parttype7 = f.create_group("PartType7")
    parttype7.create_dataset("Coordinates", data=np.random.rand(num_sidm, 3))
    parttype7.create_dataset("Velocities", data=np.random.rand(num_sidm, 3))
    parttype7.create_dataset("ParticleIDs", data=np.arange(1, num_sidm+1))
    parttype7.create_dataset("Masses", data=np.ones(num_sidm))