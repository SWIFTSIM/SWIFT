from swiftsimio import load
import numpy as np
import matplotlib.pyplot as plt

attrnum = 15
snapnum = 10
file_prefix = "sodShock_000"

break_ptcle_id = 1984

attrev = np.zeros((attrnum,snapnum))

for ii in range(snapnum):

	filename = file_prefix + str(ii) + ".hdf5"
	data = load(filename)

	attrs = data.metadata.gas_properties.field_names

	ids = data.gas.particle_ids
	idx = np.where(ids == break_ptcle_id)[0][0]

	for jj in range(attrnum):
		#print(attrs[jj])
		prop = getattr(data.gas, attrs[jj])
		val = prop[idx].value
		#print(val)
		if val.ndim > 0: 
			attrev[jj,ii] = val[2]
		else:
			attrev[jj,ii] = val

fig, ax = plt.subplots(nrows=3, ncols=5)

count = 0

for ii in range(3):
	for jj in range(5):
		ax[ii,jj].plot(attrev[count,:], label=attrs[count])

		leg = ax[ii,jj].legend(handlelength=0, handletextpad=0, fancybox=True)
		for item in leg.legendHandles:
    			item.set_visible(False)

		count += 1

plt.tight_layout()
plt.show() 

quit()

