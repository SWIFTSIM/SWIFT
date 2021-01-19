import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

default_group_id = 2147483647
nr_groups = 820

swift_file = sys.argv[1]
vr_file = sys.argv[2]

if len(sys.argv) < 4:
    group_id = -1
else:
    group_id = int(sys.argv[3])

print("\nComparing files: {} and {}".format(swift_file, vr_file))

# Open SWIFT and VELOCIraptor files
swift = h5.File(swift_file, "r")
vr_parts = h5.File(vr_file + ".catalog_particles.0", "r")
vr_parts_ub = h5.File(vr_file + ".catalog_particles.unbound.0", "r")
vr_props = h5.File(vr_file + ".properties.0", "r")
vr_groups = h5.File(vr_file + ".catalog_groups.0", "r")

# Extract particle properties into arrays
sw_group_ids = swift["/PartType1/FOFGroupIDs"][:]
sw_pids_global = swift["/PartType1/ParticleIDs"][:]
sw_pos = swift["/PartType1/Coordinates"][:, :]

vr_pids_global = vr_parts["/Particle_IDs"]
vr_pids_ub = vr_parts_ub["/Particle_IDs"]
vr_offset = vr_groups["/Offset"]
vr_offset_ub = vr_groups["/Offset_unbound"]

# Find which groups to compare
if group_id == -1:
    groups = np.linspace(1, nr_groups, num=nr_groups, dtype=int)
else:
    print("\nOnly comparing group: {}".format(group_id))
    groups = [group_id]
    
# Loop over groups and compare particle IDs between SWIFT and VELOCIraptor
for i, gid in enumerate(groups):
    mask_groups = np.where(sw_group_ids == gid)
    
    sw_group_pos = sw_pos[mask_groups]
    sw_pids = sw_pids_global[mask_groups]

    # Index to end of array if at the end of the list
    if i == nr_groups - 1:
        vr_pids = np.concatenate((vr_pids_global[vr_offset[gid - 1]:], 
                                vr_pids_ub[vr_offset_ub[gid - 1]:]))
    else:
        vr_pids = np.concatenate((vr_pids_global[vr_offset[gid - 1]:vr_offset[gid]], 
                                vr_pids_ub[vr_offset_ub[gid - 1]:vr_offset_ub[gid]]))
    
    id_list = set(sw_pids).symmetric_difference(set(vr_pids))
    
    if len(id_list) == 0:
        print("Particles in group {} match! SWIFT group size: {}, VELOCIraptor group size: {}".format(gid, len(sw_pids), len(vr_pids)))
    else:
        print("Particles in group {} don't match! SWIFT group size: {}, VELOCIraptor group size: {}".format(gid, len(sw_pids), len(vr_pids)))
        print("There are {} differences.".format(len(id_list)))
