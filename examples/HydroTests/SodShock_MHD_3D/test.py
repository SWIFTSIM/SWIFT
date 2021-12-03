from swiftsimio import load
import numpy as np
import matplotlib.pyplot as plt

data = load("sodShock_0005.hdf5")

#print(data.metadata.gas_properties.field_names)

attrs = data.metadata.gas_properties.field_names

ids = data.gas.particle_ids
idx = np.where(ids == 2089)[0][0]
print(idx)

#ii = 0

for attr in attrs:
	#print(ii)
	#ii += 1
	#print(attr + "Hello World")
	prop = getattr(data.gas, attr)
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print(attr + ": ", end="")
	print(prop[idx])

quit()
