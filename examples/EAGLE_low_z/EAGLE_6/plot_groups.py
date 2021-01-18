import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

default_group_id = 2147483647

in_file = sys.argv[1]

if len(sys.argv) < 3:
    plot_id = -1
else:
    plot_id = int(sys.argv[2])

print("\nReading file: {}".format(in_file))

f = h5.File(in_file, "r")

group_ids = f["/PartType1/FOFGroupIDs"][:]
pos = f["/PartType1/Coordinates"][:, :]

if plot_id == -1:
    mask_groups = np.where(group_ids != default_group_id)
else:
    print("\nOnly plotting group: {}".format(plot_id))
    mask_groups = np.where(group_ids == plot_id)
    print(mask_groups)

group_pos = pos[mask_groups]
big_group_ids = group_ids[mask_groups]
pos_x = pos[:, 0]
pos_y = pos[:, 1]
pos_z = pos[:, 2]
group_x = group_pos[:, 0]
group_y = group_pos[:, 1]
group_z = group_pos[:, 2]

plt.figure(figsize=(20,12))
ax1 = plt.subplot(121)
plt.scatter(pos_x, pos_y, s=1, edgecolors='none', marker=',')
plt.xlim(0, np.max(pos_x))
plt.ylim(0, np.max(pos_y))

plt.subplot(122, sharex=ax1, sharey=ax1)
plt.scatter(group_x, group_y, c=big_group_ids, s=1, edgecolors='none', marker=',')
#plt.scatter(group_x, group_y, c=big_group_ids, s=20, marker='o')
plt.xlim(0, np.max(pos_x))
plt.ylim(0, np.max(pos_y))
plt.colorbar()

#fig = plt.figure()
#ax = plt.axes(projection='3d')
##ax.scatter3D(pos_x, pos_y, pos_z, c=group_ids, cmap='Greens', s=1, edgecolors='none', marker=',');
#ax.scatter3D(group_x, group_y, group_z, c=big_group_ids, s=1, edgecolors='none', marker=',');

plt.show()
