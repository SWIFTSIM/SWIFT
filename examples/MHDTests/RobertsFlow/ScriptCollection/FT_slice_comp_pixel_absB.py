import numpy as np
import h5py
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",    # Computer Modern
})


def periodic_continue(Q):

    """
    Periodically continue the image
    Parameters: Q: quantity, with size (N_x, N_y)
    Returns: 
                Qp: quantity, with size (N_x+1, N_y+1)
    """

    M, N = Q.shape
    #Q = np.arange(M*N).reshape(M, N)  # example data

    # Periodic continuation
    Qp = np.zeros((M+1, N+1), dtype=Q.dtype)
    Qp[:-1, :-1] = Q         # copy original
    Qp[-1, :-1] = Q[0, :]     # wrap last row
    Qp[:-1, -1] = Q[:, 0]     # wrap last column
    Qp[-1, -1]  = Q[0, 0]     # wrap bottom-right corner
    return Qp

def compute_L2(in_data, ref_data, phase=False):
    #in_data=in_data[:-1,:-1] 
    #ref_data=ref_data[:-1,:-1] 
 
    diff = ref_data - in_data

    diff = np.abs(diff)

    #diff = (diff + np.pi/2) % (2*np.pi/2) - np.pi/2
    #diff = np.abs(diff)
 
    L2 = np.sqrt(np.mean(diff**2))
    L2_rel = L2 #/ np.sqrt(np.sum(np.abs(ref_data)**2))
    return L2_rel

def read_h5(file):
    with h5py.File(file, "r") as f:
        snap={
         #'r1':f["r1"][:],
         #'r2':f["r2"][:],
         #'r3':f["r3"][:],
         #'theta1':f["theta1"][:],
         #'theta2':f["theta2"][:],
         #'theta3':f["theta3"][:],
         #'B1':f["r1"][:] * np.exp(1j*f["theta1"][:]),
         #'B2':f["r2"][:] * np.exp(1j*f["theta2"][:]),
         #'B3':f["r3"][:] * np.exp(1j*f["theta3"][:]),
         'B1':f["ReB1"][:] +1j*f["ImB1"][:],
         'B2':f["ReB2"][:] +1j*f["ImB2"][:],
         'B3':f["ReB3"][:] +1j*f["ImB3"][:],
        }  
    return snap


from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

def plot_difference(input_file, reference_file, name,res_xy,res_z):
    B1_in = input_file['B1']
    B2_in = input_file['B2']
    B3_in = input_file['B3']
    B1_ref = reference_file['B1']
    B2_ref = reference_file['B2']
    B3_ref = reference_file['B3']
   
    B1_diff = (B1_in - B1_ref) #/(np.abs(B1_ref)+0.01)
    B2_diff = (B2_in - B2_ref) #/(np.abs(B2_ref)+0.01)
    B3_diff = (B3_in - B3_ref) #/(np.abs(B3_ref)+0.01)

    absB_ref = np.sqrt(np.abs(B1_ref)**2+np.abs(B2_ref)**2+np.abs(B3_ref)**2) 
    absB_in = np.sqrt(np.abs(B1_in)**2+np.abs(B2_in)**2+np.abs(B3_in)**2) 
    absB_diff = np.sqrt(np.abs(B1_diff)**2+np.abs(B2_diff)**2+np.abs(B3_diff)**2) 

    absB_diff_rms = np.sqrt(np.mean(absB_diff**2))

    # Generate figure
    nx = 3
    ny = 1
    fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx, 5 * ny))


    from matplotlib.ticker import AutoMinorLocator
    for a in ax:
        # Ticks on all four sides
        a.tick_params(axis='both',
                      which='major',
                      direction='in',
                      top=True,
                      right=True,
                      length=7,
                      width=1.2,
                      labelsize=16)
    
        a.tick_params(axis='both',
                      which='minor',
                      direction='in',
                      top=True,
                      right=True,
                      length=4,
                      width=1.0)
    
        # Add minor ticks
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())

    ax[2].tick_params(color='w',which='major')
    ax[2].tick_params(color='w',which='minor')


    # Number of color levels
    nlev=100 # 17

    # Colormaps
    # Custom colormaps
    # Select colormaps for r and theta
    cmapr = "viridis" #c1r 
    cmaptheta = "twilight_shifted" #c1r #"RdBu_r"
    a00 = ax[0].pcolormesh(
        periodic_continue(np.abs(absB_in)),
        cmap='plasma',
        vmin=0.0,
        vmax=1.5,
        shading="nearest"  # important!
    )
    a01 = ax[1].pcolormesh(
        periodic_continue(np.abs(absB_ref)),
        cmap='plasma',
        vmin=0.0,
        vmax=1.5,
        shading="nearest"  # important!
    )

    a02 = ax[2].pcolormesh(
        periodic_continue(np.abs(absB_diff)),#/(absB_diff_rms),#/np.abs(absB_ref)),
        cmap=cmapr,
        vmin=0.0,
        vmax=0.1,
        shading="nearest"  # important!
    )
     # Add ticks and labels
    locs = [res_xy/4 , res_xy/2, 3 * res_xy/4, res_xy]
    labels = [r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$',r'$2\pi$']
    #for ii in range(ny):
    #    #ax[ii, 0].set_ylabel(r"$y$",fontsize=20)
    #    ax[ii, 0].set_yticks(locs, labels,fontsize=15)
    ax[0].set_yticks(locs, labels,fontsize=25)
    #for ii in range(ny):
    #    for jj in range(nx):
    #        #ax[ii, jj].set_xlabel(r"$x$",fontsize=20)
    #        ax[ii, jj].set_xticks(locs, labels,fontsize=15)
    #        #ax[ii, jj].set_aspect("equal")
    #        ax[ii, jj].set_aspect("equal", adjustable="box")
    #        ax[ii, jj].set_xlim(0, res_xy)
    #        ax[ii, jj].set_ylim(0, res_xy)
    for jj in range(nx):
        ax[jj].set_xticks(locs, labels,fontsize=25)
        ax[jj].set_aspect("equal", adjustable="box")
        ax[jj].set_xlim(0, res_xy)
        ax[jj].set_ylim(0, res_xy)


    # Add colorbar
    ticks_r = [0.0,0.5,1.0,1.5]
    ticks_err = [0.0,0.05,0.1]
    ticks_phase = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
    ticklabels_phase = [r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"]
    cbar1 = plt.colorbar(a00, ax=ax[0], fraction=0.046, pad=0.04, ticks=ticks_r)
    cbar1.set_ticklabels(ticks_r,fontsize=15)
    #cbar1.set_label(r"$r_1 \;$",fontsize=20,rotation=0)
    ax[0].set_title(r"$|\tilde{B}_{\rm SWIFT}|/B_{\rm rms} \;$",fontsize=30,pad=10)
    cbar2 = plt.colorbar(a01, ax=ax[1], fraction=0.046, pad=0.04, ticks=ticks_r)
    cbar2.set_ticklabels(ticks_r,fontsize=15)
    ax[1].set_title(r"$|\tilde{B}_{\rm Pencil}|/B_{\rm rms} \;$",fontsize=30,pad=10)
    cbar3 = plt.colorbar(a02, ax=ax[2], fraction=0.046, pad=0.04, ticks=ticks_err)
    cbar3.set_ticklabels(ticks_err,fontsize=15)
 
    #cbar2.set_label(r"$r_2 \;$",fontsize=20,rotation=0)
    ax[2].set_title(r"$|\delta \tilde{B}|/B_{\rm rms} \;$",fontsize=30,pad=10)
    #cbar3 = plt.colorbar(a02, ax=ax[2], fraction=0.046, pad=0.04, ticks=ticks_r)
    #cbar3.set_ticklabels(ticks_r,fontsize=15)
    #cbar3.set_label(r"$r_3 \;$",fontsize=20,rotation=0)
    #ax[0,2].set_title(r"$|\delta B_3/B_3| \;$",fontsize=25,pad=10)
    #cbar4 = plt.colorbar(a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar4.set_label(r"$\theta_1 \;$",fontsize=20, rotation=0)
    #cbar4.set_ticklabels(ticklabels_phase,fontsize=15)
    #ax[1,0].set_title(r"$\rm arg(\delta B_1)\;$",fontsize=25,pad=10)
    #cbar5 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar5.set_label(r"$\theta_2 \;$",fontsize=20,rotation=0)
    #cbar5.set_ticklabels(ticklabels_phase,fontsize=15)
    #ax[1,1].set_title(r"$\rm arg(\delta B_2) \;$",fontsize=25,pad=10)
    #cbar6 = plt.colorbar(a12, ax=ax[1, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar6.set_label(r"$\theta_3 \;$",fontsize=20,rotation=0)
    #cbar6.set_ticklabels(ticklabels_phase,fontsize=15)
    #ax[1,2].set_title(r"$\rm arg(\delta B_3) \;$",fontsize=25,pad=10)
    ##cbar7 = plt.colorbar(a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar7.set_label(r"$\tilde{\theta}_1(x,y) \;$")
    ##cbar8 = plt.colorbar(a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar8.set_label(r"$\tilde{\theta}_2(x,y) \;$")
    ##cbar9 = plt.colorbar(a22, ax=ax[2, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
    ##cbar9.set_label(r"$\tilde{\theta}_3(x,y) \;$")

    fig.tight_layout()

    plt.savefig('pixel_absB_'+name+'.png', dpi=220)


if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Generate ICs for ABC flow")

#    parser.add_argument(
#        "-in",
#        "--input_snapshots",
#        help="input folder",
#        default='./input/*.h5',
#        type=str,
#    )
#    parser.add_argument(
#        "-ref",
#        "--reference_snapshots",
#        help="reference folder",
#        default="./reference/*.h5",
#        type=str,
#    )
    parser.add_argument(
        "-out",
        "--output_imagename",
        help="output image name",
        default='abs_delta_B_rms.png',
        type=str,
    )
 
    args = parser.parse_args()

    data        = {'SWIFT':{'16':'./input/slice_SWIFT_016_16_32.h5',
                            '32':'./input/slice_SWIFT_032_32_64.h5',
                            '64':'./input/slice_SWIFT_064_64_128.h5',
                   },
                   'Pencil':{'16':'./reference/slice_Pencil_016_16_32.h5',
                                 '32':'./reference/slice_Pencil_032_32_64.h5',
                                 '64':'./reference/slice_Pencil_064_64_128.h5',
                      }}

    resolution = ['16','32','64']
    resolution_val = [16,32,64]
    keys = ['B1','B2','B3']
    key_names = [r'$B_1$',r'$B_2$',r'$B_3$']
    colors = ['black','black','black']
    markers = ['x','+','o']
    markers_sizes = [140,200,200]
 
    #keys = ['r1','r2','r3','theta1','theta2','theta3']
    #key_names = [r'$r_1$',r'$r_2$',r'$r_3$',r'$\theta_1$',r'$\theta_2$',r'$\theta_3$']
    #colors = ['black','black','black','black','black','black']
    #markers = ['x','+','o','^','>','<']
    #markers_sizes = [70,100,100,100,100,100]
    results = []
    for i in range(len(resolution)):
        res = resolution[i]
        input_file = read_h5(data['SWIFT'][res])
        reference_file = read_h5(data['Pencil'][res])
        plot_difference(input_file, reference_file, res, res_xy=resolution_val[i], res_z=32)
        res_col=[]
        for key in keys:
            res_col.append(compute_L2(input_file[key], reference_file[key]))
        results.append(res_col)
    results = np.array(results).T

    nx = 1
    ny = 1
    lx = 1.61
    ly = 1.0
    fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx * lx, 5 * ny * ly))

 
    from matplotlib.ticker import AutoMinorLocator
    ax.tick_params(axis='both',
                  which='major',
                  direction='in',
                  top=True,
                  right=True,
                  length=7,
                  width=1.2,
                  labelsize=16)
    
    ax.tick_params(axis='both',
                  which='minor',
                  direction='in',
                  top=True,
                  right=True,
                  length=4,
                  width=1.0)
    
    # Add minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


       
    delta_B2 = 0.0
    for i in range(len(keys)):
        delta_B2 += results[i]**2
        data_x = np.array(resolution,dtype=int)
    data_y = np.sqrt(delta_B2) 

    ax.plot(data_x,data_y,marker='x',markersize=7,linestyle='',color='red',markeredgewidth=2)

    #for i in range(len(keys)):
    #    name = key_names[i]
    #    data_y = results[i]
    #    data_x = np.array(resolution,dtype=int)
    #    ax.scatter(data_x,data_y,marker=markers[i],label=name,s=markers_sizes[i],color=colors[i])

    xlin = np.logspace(np.log10(18),np.log10(58),10)
    plt.plot(xlin,0.27*(xlin/16)**(-1),linestyle='dashed',color='black',label=r'$N^{-1}$')
    plt.plot(xlin,0.22*(xlin/16)**(-2),linestyle='dashdot',color='black',label=r'$N^{-2}$')

    ax.text(x=54,y=1e-1, s='$N_\perp^{-1}$',fontsize=20,color='black')
    ax.text(x=46,y=1.4e-2, s='$N_\perp^{-2}$',fontsize=20,color='black')

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_yticks([1e0,1e-1,1e-2],[r'$10^{0}$',r'$10^{-1}$',r'$10^{-2}$'],fontsize=25)
    ax.set_xticks(data_x,data_x,fontsize=25)
    #ax.set_ylabel('L2 norm',fontsize=20)
    ax.set_xlabel(r'$N_\perp$',fontsize=30)
    ax.set_ylabel(r'$\langle |\delta \tilde B| \rangle_{\rm rms} / B_{\rm rms}$',fontsize=30)
    #ax.set_title(r'$\langle |\delta \tilde B| \rangle_{\rm rms} / B_{\rm rms}$',fontsize=30)
    #ax.legend(fontsize=14)
    #ax.grid()

    fig.tight_layout()
    plt.savefig(args.output_imagename)


