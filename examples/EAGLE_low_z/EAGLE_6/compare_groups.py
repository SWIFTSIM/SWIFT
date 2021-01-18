import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

default_group_id = 2147483647

swift_file = sys.argv[1]
vr_file = sys.argv[2]

if len(sys.argv) < 3:
    group_id = -1
else:
    group_id = int(sys.argv[3])

print("\nComparing files: {} and {}".format(swift_file, vr_file))

swift = h5.File(swift_file, "r")
vr_parts = h5.File(vr_file + ".catalog_particles.0", "r")
vr_parts_ub = h5.File(vr_file + ".catalog_particles.unbound.0", "r")
vr_props = h5.File(vr_file + ".properties.0", "r")
vr_groups = h5.File(vr_file + ".catalog_groups.0", "r")

sw_group_ids = swift["/PartType1/FOFGroupIDs"][:]
sw_pids = swift["/PartType1/ParticleIDs"][:]
sw_pos = swift["/PartType1/Coordinates"][:, :]

vr_pids = vr_parts["/Particle_IDs"]
vr_pids_ub = vr_parts_ub["/Particle_IDs"]
vr_offset = vr_groups["/Offset"]
vr_offset_ub = vr_groups["/Offset_unbound"]

if group_id == -1:
    mask_groups = np.where(sw_group_ids != default_group_id)
else:
    print("\nOnly comparing group: {}".format(group_id))
    mask_groups = np.where(sw_group_ids == group_id)

sw_group_pos = sw_pos[mask_groups]
sw_pids = sw_pids[mask_groups]

vr_pids = np.concatenate((vr_pids[vr_offset[group_id - 1]:vr_offset[group_id]], 
                          vr_pids_ub[vr_offset_ub[group_id - 1]:vr_offset_ub[group_id]]))

print("SWIFT group size: {}".format(len(sw_pids)))
print("VELOCIraptor group size: {}".format(len(vr_pids)))

id_list = set(sw_pids).symmetric_difference(set(vr_pids))

if len(id_list) == 0:
    print("Particles in group {} match in both SWIFT and VELOCIraptor!".format(group_id))
else:
    print("Particles in group {} don't match!".format(group_id))
    print("SWIFT group size: {}".format(len(sw_pids)))
    print("VELOCIraptor group size: {}".format(len(vr_pids)))
    #print("There are {} differences.".format(len(id_list))
