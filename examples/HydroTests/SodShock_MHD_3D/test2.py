from swiftsimio import load
import numpy as np
import matplotlib.pyplot as plt

attrnum = 16
snapnum = 9
file_prefix = "sodShock_000"

break_ptcle_id = 19823

attrev = np.zeros((attrnum,snapnum))

for ii in range(snapnum):

	filename = file_prefix + str(ii) + ".hdf5"
	data = load(filename)

	attrs = data.metadata.gas_properties.field_names

	print(attrs)
	quit()

	ids = data.gas.particle_ids
	idx = np.where(ids == break_ptcle_id)[0][0]

	for jj in range(attrnum):
		#print(attrs[jj])
		prop = getattr(data.gas, attrs[jj])
		val = prop[idx].value
		#print(val)
		if val.ndim > 0: 
			attrev[jj,ii] = val[1]
		else:
			attrev[jj,ii] = val

#count = 8
#plt.plot(attrev[count,:],label=attrs[count])
#plt.show()
#quit()

fig, ax = plt.subplots(nrows=int(4), ncols=int(4))

count = 0

for ii in range(int(4)):
	for jj in range(int(4)):
		ax[ii,jj].plot(attrev[count,:], label=attrs[count])

		leg = ax[ii,jj].legend(handlelength=0, handletextpad=0, fancybox=True)
		for item in leg.legendHandles:
    			item.set_visible(False)

		count += 1

#plt.tight_layout()
plt.show() 

quit()
