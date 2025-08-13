import h5py as h5
import numpy as np

f = h5.File("M1_R0.01296_Z1_S0_A0_B0.1_I1_Res43_n2_sol0.5_42.hdf5", "a")

# Fix header
f["Header"].attrs.create("Flag_Entropy_ICs", 0)
f["Header"].attrs.create("BoxSize", 5.18e-02)
f["Header"].attrs.create("NumPart_Total_HighWord", [0, 0, 0, 0, 0, 0])

N = f["Header"].attrs.get("NumPart_ThisFile")[0]

B = f["/PartType0/MagneticField"][:, :] 
B *= 81 / 3.409e-11

f["/PartType0"].create_dataset("SmoothingLength", (N,), 'f', data=np.ones(N) * 0.001)
f["/PartType0"].create_dataset("MagneticFluxDensities", (N, 3), 'f', data=B)

f.close()
