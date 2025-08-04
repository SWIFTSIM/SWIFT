import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

def load_and_extract(filename, type_='masses'):
    """       
    Load the data and extract relevant info.                                                                                                    
    """

    with h5.File(filename, "r") as f:
        try:
            pos = f["PartType0/Coordinates"][...]
        except:
            pos = f["/PartType0/Coordinates"][...]
        h = f["/PartType0/SmoothingLengths"][...]
        ids = f["PartType0/ParticleIDs"][...]
        masses = f["PartType0/Masses"][...]
        
        return h, ids, masses

h_list_large = []
h_list_norm = []
for i in range(3):
    filename = f"uniformBox_{i:04d}.hdf5"
    h,pos,m = load_and_extract(filename, type_='masses')
    h_list_large.append(h[0])
    h_list_norm.append(h[2])
    print('-------')
    print(f'large particle info h:{h[0]}, m:{m[0]}, pos:{pos[0]}')
    print('--------')
    print(f'normal particle info h:{h[2]}, m:{m[2]}, pos:{pos[2]}') 


snap_list = np.arange(0,len(h_list_large),1)
plt.plot(snap_list, h_list_large, label='small particle')
plt.plot(snap_list, h_list_norm, label='normal particle')
plt.title('mass-weighted')
plt.xlabel('snap')
plt.ylabel('h')
plt.legend()
plt.savefig('uniformbox.png')
