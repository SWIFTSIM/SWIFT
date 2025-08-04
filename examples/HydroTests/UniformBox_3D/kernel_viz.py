"""
SPH description plot. Somewhat complex, but just shows:
3D distribution with kernels.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import math
import h5py as h5
from scipy.optimize import newton
from matplotlib.colors import LogNorm
import swiftsimio as sw
from swiftsimio.visualisation.projection_backends.fast import scatter

from scipy.spatial.distance import pdist, squareform

kernel_gamma_3D = 2

def create_ax_right(ax):
    filename = "uniformBox_0001.hdf5"
    with h5.File(filename, "r") as f:
        positions = f["PartType0/Coordinates"][...]
        if positions.shape[0] == 3 and positions.shape[1] > 3:
            positions = positions.T
        smoothing_length = f["PartType0/SmoothingLengths"][...]
        densities = f['PartType0']['Densities'][...]
        masses = f["PartType0/Masses"][...]
        vels = f["PartType0/Velocities"][...]
        press = f["PartType0/Pressures"][...]
        
    idx = np.argmin(np.abs(positions[:,0] - 0.5))
    
    # Slice thin z slab
    z_center = positions[:,2][idx]
    print(vels[idx])
    z_tol = 1.0
    mask = np.abs(positions[:, 2] - z_center) < z_tol
    positions_3D = positions[mask]
    masses_3D = masses[mask]
    hsml_out = smoothing_length[mask]

    # Zoomed region
    zoom_box = [positions[:,0][idx]-0.2, positions[:,0][idx]+0.2, positions[:,1][idx]-0.2,positions[:,1][idx]+0.2]
    zoom_mask = ( 
    (positions_3D[:, 0] > zoom_box[0]) & (positions_3D[:, 0] < zoom_box[1]) &
    (positions_3D[:, 1] > zoom_box[2]) & (positions_3D[:, 1] < zoom_box[3])
    )
    
    positions_zoom_3D = positions_3D[zoom_mask]
    masses_zoom = masses_3D[zoom_mask]
    hsml_zoom_out = hsml_out[zoom_mask]
    # Select random particles in zoom region (full 3D positions)
    target = positions[idx]
    distances = np.linalg.norm(positions_zoom_3D - target, axis=1)
    selected_idx = np.argmin(distances)
    selected_pos_3D = positions_zoom_3D[selected_idx:selected_idx+1]
    #kernel_gamma = 2 in 3D
    hsml_out_sel = hsml_zoom_out[selected_idx:selected_idx+1] *kernel_gamma_3D

    colors = ['blue']
    
    # Plot zoomed-in region
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_title("Zoomed-In Kernel Comparison")
    
    # Background particles
    ax.scatter(positions_zoom_3D[:, 0], positions_zoom_3D[:, 1], s=4, color='lightgray')
    
    # Plot selected particles with kernel circles
    for i, pos in enumerate(selected_pos_3D):
        ax.plot(pos[0], pos[1], 'o', color=colors[i], label=f"Part {i}")
        
        # SWIFT HDF5 output
        circ_out = Circle(pos, radius=2*hsml_out_sel[i], edgecolor=colors[i], linestyle='dotted', fill=False, label=f"SWIFT h_{i}")
        ax.add_patch(circ_out)
    
    ax.set_xlim(zoom_box[0], zoom_box[1])
    ax.set_ylim(zoom_box[2], zoom_box[3])
    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(fontsize=8, loc='upper right')
    plt.tight_layout()
    # Compute distances to all particles in the slice
    distances = np.linalg.norm(positions_3D - selected_pos_3D[0], axis=1)
    
    # Count neighbors within h
    h = hsml_out_sel[0]
    within_h = distances < h
    num_within_h = np.sum(within_h)
    
    # Count neighbors within kernel support radius (γh)
    within_gamma_h = distances < kernel_gamma_3D * h
    
    # Print results
    print(f"Number of neighbors within h: {num_within_h}")
    
if __name__ == "__main__":
    # Actually make the plots

    
    # First we need to set up our image.
    fig,ax = plt.subplots(figsize=(8, 8))
    create_ax_right(ax)
    
    plt.savefig("nhb_kernel_test.png")
