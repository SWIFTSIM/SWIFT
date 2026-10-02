import unyt
import matplotlib.pyplot as plt
import argparse
import scipy as sp

from swiftsimio import load
from swiftsimio.visualisation.volume_render import render_gas
from swiftsimio.visualisation.slice import slice_gas
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pencil as pc


import pandas as pd
import os
from glob import glob
#import numpy as np
#import matplotlib.pyplot as plt
import h5py

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.volume_render import render_gas
import matplotlib.ticker as mticker
from swiftsimio import load_statistics
#run_directory = './'

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",    # Computer Modern
})

# Load snapshot

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


def make_2D_slice_z(data,
                    height,
                    res_xy
                    ):

    """
    Making SPH B slice of the box at specific height
    Parameters: data: SPH snapshot data
                height: slice height
    Returns: 
                dict: B - magnetic field array of size (N_x, N_y, 3), z - slice height
    """

    
    Lbox = data.metadata.boxsize

    # compensate half-grid shift
    dx = Lbox[0]/(2*res_xy)
    dy = Lbox[1]/(2*res_xy)

    # select region
    visualise_region_xy = [
    0*unyt.cm-dx,
    Lbox[0]-dx,
    0*unyt.cm-dy,
    Lbox[1]-dy,
    ]

    common_arguments_xy = dict(
        data=data,
        resolution=res_xy,
        parallel=True,
        region=visualise_region_xy,
        z_slice=height,  
        periodic=True,
    )
    mass_map_ij  = slice_gas(**common_arguments_xy, project="masses")
    mass_weighted_Bx_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_Bx")
    mass_weighted_By_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_By")
    mass_weighted_Bz_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_Bz")
    mass_weighted_vx_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_vx")
    mass_weighted_vy_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_vy")
    mass_weighted_vz_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_vz")
 
    Bx_map_ij  = mass_weighted_Bx_map_ij /mass_map_ij 
    By_map_ij  = mass_weighted_By_map_ij /mass_map_ij 
    Bz_map_ij  = mass_weighted_Bz_map_ij /mass_map_ij 

    vx_map_ij  = mass_weighted_vx_map_ij /mass_map_ij 
    vy_map_ij  = mass_weighted_vy_map_ij /mass_map_ij 
    vz_map_ij  = mass_weighted_vz_map_ij /mass_map_ij 

    # Invert y axis to convert ij map to xy map 
    Bx_map_xy = Bx_map_ij #np.flipud(Bx_map_ij)
    By_map_xy = By_map_ij #np.flipud(By_map_ij)
    Bz_map_xy = Bz_map_ij #np.flipud(Bz_map_ij)

    vx_map_xy = vx_map_ij #np.flipud(Bx_map_ij)
    vy_map_xy = vy_map_ij #np.flipud(By_map_ij)
    vz_map_xy = vz_map_ij #np.flipud(Bz_map_ij)

    B = np.stack((Bx_map_xy, By_map_xy, Bz_map_xy), axis=-1)
    v = np.stack((vx_map_xy, vy_map_xy, vz_map_xy), axis=-1)
    return {'B':B, 'v':v, 'z':height}

def make_SPH_grid(data,
                    res_xy,
                    res_z
                    ):
    """
    Rendering B field from SPH particle disrtibution
    Parameters: data: SPH snapshot data
                res_xy: number of samples along x or y
                res_z: number of samples along x or y
    Returns: 
                Bx,By,Bz: 3 3d arrays of shape (N_x, N_y, N_z)
    """
    
    # Retrieve particle attributes of interest
    B = data.gas.magnetic_flux_densities
    v = data.gas.velocities

    Lbox = data.metadata.boxsize
    # Generate mass weighted maps of quantities of interest
    data.gas.mass_weighted_Bx = data.gas.masses * B[:,0]
    data.gas.mass_weighted_By = data.gas.masses * B[:,1]
    data.gas.mass_weighted_Bz = data.gas.masses * B[:,2]
    data.gas.mass_weighted_vx = data.gas.masses * v[:,0]
    data.gas.mass_weighted_vy = data.gas.masses * v[:,1]
    data.gas.mass_weighted_vz = data.gas.masses * v[:,2]


    slices = np.linspace(0,Lbox[2].value, res_z, endpoint=False)*Lbox[2].units

    Bx = []
    By = []
    Bz = []

    vx = []
    vy = []
    vz = []

    for slice in slices:
        Bx.append(make_2D_slice_z(data, slice, res_xy)['B'][:,:,0])
        By.append(make_2D_slice_z(data, slice, res_xy)['B'][:,:,1])
        Bz.append(make_2D_slice_z(data, slice, res_xy)['B'][:,:,2])
        vx.append(make_2D_slice_z(data, slice, res_xy)['v'][:,:,0])
        vy.append(make_2D_slice_z(data, slice, res_xy)['v'][:,:,1])
        vz.append(make_2D_slice_z(data, slice, res_xy)['v'][:,:,2])

    # make with shape N_x, N_y, N_z
    Bx = np.array(Bx).transpose(1, 2, 0)
    By = np.array(By).transpose(1, 2, 0)
    Bz = np.array(Bz).transpose(1, 2, 0)
    vx = np.array(vx).transpose(1, 2, 0)
    vy = np.array(vy).transpose(1, 2, 0)
    vz = np.array(vz).transpose(1, 2, 0)


    return {'B':[Bx,By,Bz],'v':[vx,vy,vz]}


def read_hdf5_B(input_filename, res_xy, res_z):
    # Load snapshot
    data = load(input_filename)

    time = np.round(data.metadata.time,2)
    print("Plotting at ", time)

    Lbox = data.metadata.boxsize

    uu = data.gas.velocities
    vrms = np.sqrt(np.mean(uu[:,0]**2+uu[:,1]**2+uu[:,2]**2))
    print('vrms',vrms)

    grid_data = make_SPH_grid(data, res_xy, res_z)
    grid_data['Lbox'] = Lbox.value 
 
    return grid_data 


def read_VAR_B(input_filename):
    snapshotid = int(input_filename.split('/')[-1].split('.VAR')[0])
    name = str(snapshotid)+'.VAR'
    datafolder = input_filename.split(name)[0]
    snap = pc.read.var(datadir=datafolder,ivar=snapshotid,trimall=True)
    aa = snap.aa
    #uu = snap.uu.transpose(3, 2, 1,0)
    #uu2 = np.sum(uu * uu, axis=-1)
    #vrms = np.sqrt(np.mean(uu2))
    #print('vrms',vrms)
    pad = 6  # depends on stencil (likely 3 for Pencil)
    aa = np.pad(
        aa,
        ((0,0), (pad,pad), (pad,pad), (pad,pad)),
        mode='wrap'
    )
    bb = pc.math.derivatives.div_grad_curl.curl(
        aa,
        dx=snap.dx,
        dy=snap.dy,
        dz=snap.dz
    )
    bb = bb[:, pad:-pad, pad:-pad, pad:-pad]
    bb = bb.transpose(0, 3, 2, 1)
    #uu = uu[:, pad:-pad, pad:-pad, pad:-pad]
    uu = snap.uu
    uu = uu.transpose(0, 3, 2, 1)
    #Bx,By,Bz = bb
    Lbox = [2*np.pi,2*np.pi, 4*np.pi]
   
    grid_data = {'B': bb, 'v': uu, 'Lbox':Lbox}  

    del aa, bb, uu, snap, pad, datafolder, name, snapshotid
    return grid_data

def compute_quantities_3DFFT(grid_data, res_xy, res_z):

    Bx,By,Bz = grid_data['B']
    vx,vy,vz = grid_data['v']
    Lbox = grid_data['Lbox']
    dx,dy,dz = Lbox[0]/res_xy, Lbox[1]/res_xy, Lbox[2]/res_z

    # Get vector potential
    dV = dx * dy * dz
    Bx_k = np.fft.fftn(Bx) * dV
    By_k = np.fft.fftn(By) * dV
    Bz_k = np.fft.fftn(Bz) * dV
    kx1d = 2*np.pi * np.fft.fftfreq(res_xy, d=dx)
    ky1d = 2*np.pi * np.fft.fftfreq(res_xy, d=dy)
    kz1d = 2*np.pi * np.fft.fftfreq(res_z, d=dz)
    kx, ky, kz = np.meshgrid(kx1d, ky1d, kz1d, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    k2[0, 0, 0] = 1.0
    # projection to divergence-free field
    k_dot_B = kx*Bx_k + ky*By_k + kz*Bz_k
    Bx_k -= kx * k_dot_B / k2
    By_k -= ky * k_dot_B / k2
    Bz_k -= kz * k_dot_B / k2
    Ax_k = 1j * (ky * Bz_k - kz * By_k) / k2
    Ay_k = 1j * (kz * Bx_k - kx * Bz_k) / k2
    Az_k = 1j * (kx * By_k - ky * Bx_k) / k2
    Ax_k[0,0,0] = 0
    Ay_k[0,0,0] = 0
    Az_k[0,0,0] = 0
    Ax = np.fft.ifftn(Ax_k).real / dV
    Ay = np.fft.ifftn(Ay_k).real / dV
    Az = np.fft.ifftn(Az_k).real / dV
    del Ax_k,Ay_k,Az_k
    Jx_k = 1j * (ky * Bz_k - kz * By_k)
    Jy_k = 1j * (kz * Bx_k - kx * Bz_k)
    Jz_k = 1j * (kx * By_k - ky * Bx_k)
    Jx = np.fft.ifftn(Jx_k).real / dV
    Jy = np.fft.ifftn(Jy_k).real / dV
    Jz = np.fft.ifftn(Jz_k).real / dV
    del Jx_k,Jy_k,Jz_k

    vx_k = np.fft.fftn(vx) * dV
    vy_k = np.fft.fftn(vy) * dV
    vz_k = np.fft.fftn(vz) * dV
    k_dot_v = kx*vx_k + ky*vy_k + kz*vz_k
    vx_k -= kx * k_dot_v / k2
    vy_k -= ky * k_dot_v / k2
    vz_k -= kz * k_dot_v / k2
    vx = np.fft.ifftn(vx_k).real / dV
    vy = np.fft.ifftn(vy_k).real / dV
    vz = np.fft.ifftn(vz_k).real / dV
    del vx_k,vy_k,vz_k

    del Bx_k,By_k,Bz_k,kx,ky,kz
    B = np.stack((Bx,By,Bz),axis=-1)
    del Bx,By,Bz
    A = np.stack((Ax,Ay,Az),axis=-1)
    del Ax,Ay,Az
    J = np.stack((Jx,Jy,Jz),axis=-1)
    del Jx,Jy,Jz
    v = np.stack((vx,vy,vz),axis=-1)
    del vx,vy,vz

    B2 = np.sum(B * B, axis=-1)
    J2 = np.sum(J * J, axis=-1)
    v2 = np.sum(v * v, axis=-1)
    A2 = np.sum(A * A, axis=-1)

    # plane averaged B
    B_bar = np.mean(B,axis = (0,1))
    B_bar2 = np.sum(B_bar * B_bar, axis=-1)

    # Rms values
    Brms = np.sqrt(np.mean(B2))
    Jrms = np.sqrt(np.mean(J2))
    vrms = np.sqrt(np.mean(v2))
    B_bar_rms = np.sqrt(np.mean(B_bar2))
    Arms = np.sqrt(np.mean(A2))

    # dot products (vectorized, no loop)
    JdotB = np.sum(J * B, axis=-1)
    AdotB = np.sum(A * B, axis=-1)
    #BdotdBdt = np.sum(B * dBdt, axis=-1)

    # JdotB rms
    JdotB_rms = np.sqrt(np.mean(JdotB*JdotB))

    # cross product (computed once, vectorized)
    JcrossB = np.cross(J, B)

    # calculate RMS growth rate
    #growth_rate_Brms = np.mean(BdotdBdt)/Brms**2

    # rms ratio
    B_bar_rms_Norm = B_bar_rms/Brms

    # WL = v · (J × B)
    WL = np.sum(v * JcrossB, axis=-1)

    # final quantities
    JcdotB_avrgNorm = np.mean(JdotB) / (Brms * Jrms)
    JcdotB2_avrgNorm = np.sqrt(np.mean(JdotB**2)) / (Brms * Jrms)

    JcrossB2_avrgNorm = np.sqrt(np.mean(np.sum(JcrossB**2, axis=-1))) / (Brms * Jrms)

    WL_avrg_norm = np.mean(WL) / (vrms * Jrms * Brms)

    J2B2_avrgNorm = np.mean(J2 * B2) / (Brms**2 * Jrms**2)
    AdotB_avrg_Norm = np.mean(AdotB) / (Brms * Brms)

    Jrms_over_Brms = Jrms/Brms
    Jrms_over_Brmsk=Jrms_over_Brms
    Arms_over_Brms = Arms/Brms
    JdotB_avrg_Norm = np.mean(JdotB) / (Brms * Brms)
    AdotB_avrg_Norm = np.mean(AdotB) / (Brms * Brms)
      
    results = {'Bbarrms/Brms':B_bar_rms_Norm,'Jrms/Brms':Jrms_over_Brmsk,'JB/JrmsBrms': JcdotB_avrgNorm, 'JBrms/JrmsBrms':JcdotB2_avrgNorm, 'JxB/JrmsBrms':JcrossB2_avrgNorm,'WL/(urmsBrmsJrms)':-WL_avrg_norm,'J2B2/(JrmsBrms)^2':J2B2_avrgNorm, 'AB/(Brms)^2':-AdotB_avrg_Norm}
    return results


#snapshot_name = "./RF1D_2nd_128_256_14/2.VAR" #"./rf1dn_128_Run2_c1_02/RobertsFlow_0002.hdf5"
snapshot_name = "./RF1D_10th_256_512_14/2.VAR" #"./rf1dn_128_Run2_c1_02/RobertsFlow_0002.hdf5"
#grid_data = read_hdf5_B(snapshot_name, res_xy = 128, res_z = 256)
grid_data = read_VAR_B(snapshot_name)
quantities =  compute_quantities_3DFFT(grid_data, res_xy=256, res_z = 512)
print(quantities)

#snapshot_name = "./RF1D_2nd_128_256_14_old/2.VAR"
#result = read_VAR_B(snapshot_name)
#print(result)

##########################################


def load_test_run_parameters():
    #run_data = pd.read_csv("test_run_parameters.csv", sep="\t")
    run_data = pd.read_csv("test_run_parameters_rf1dn.csv", sep="\t")
    mask = run_data["Status"] == "done"
    run_data = run_data[mask]
    return run_data

def process_info(run_data, snapshot_name = 'RobertsFlow_0002.hdf5',outputfile='snapshot_results.csv'):

   resulting_list = []

   for i in range(len(run_data)):

        run_data_slice = run_data.iloc[[i]]

        # Load snapshot
        filename = run_directory + str(run_data_slice["Run #"].values[0]) + '/' + snapshot_name  
        data = load(filename)
  
        Np = len(data.gas.masses)+1
        boxsize = data.metadata.boxsize
 
        # make volume rendering
        res_xy = int((Np*(boxsize[0].value/boxsize[2].value))**(1/3))
        res_z = 32 #int(0.5*(Np*(boxsize[2].value/boxsize[0].value)**2)**(1/3))
        #print(res_xy,res_z)

        #x1d = np.linspace(0,boxsize[0].value,resolution)*boxsize[0].units
        #y1d = np.linspace(0,boxsize[1].value,resolution)*boxsize[1].units
        #z1d = np.linspace(0,boxsize[2].value,resolution)*boxsize[2].units
        #X, Y, Z = np.meshgrid(x1d, y1d, z1d, indexing="ij")
        #kb0 = 1.0
        #Bx_cube = (np.sin(kb0 * Y) + np.sin(kb0 * Z))*Bx_cube.units
        #By_cube = (np.cos(kb0 * X) - np.cos(kb0 * Z))*By_cube.units
        #Bz_cube = (np.sin(kb0 * X) + np.cos(kb0 * Y))*Bz_cube.units

        #Bnew = data.gas.magnetic_flux_densities
        #pos = data.gas.coordinates
        #kb0 = 1.0
        #Bnew[:,0] = (np.sin(kb0 * pos[:,1]) + np.sin(kb0 * pos[:,2]))*Bnew.units
        #Bnew[:,1] = (np.cos(kb0 * pos[:,0]) - np.cos(kb0 * pos[:,2]))*Bnew.units
        #Bnew[:,2] = (np.sin(kb0 * pos[:,0]) + np.cos(kb0 * pos[:,1]))*Bnew.units

        #Bx_cube,By_cube,Bz_cube = make_grid(data,
        #                                    Bnew,
        #                                    res_xy = res_xy,
        #                                    res_z = res_z)
        Bx_cube,By_cube,Bz_cube = make_grid(data,
                                            data.gas.magnetic_flux_densities,
                                            res_xy = res_xy,
                                            res_z = res_z)
        dBxdt_cube,dBydt_cube,dBzdt_cube = make_grid(data,
                                            data.gas.magnetic_flux_densitiesdt,
                                            res_xy = res_xy,
                                            res_z = res_z)
        vx_cube,vy_cube,vz_cube = make_grid(data,
                                            data.gas.velocities,
                                            res_xy = res_xy,
                                            res_z = res_z)

        dx = boxsize[0] / res_xy
        dy = boxsize[1] / res_xy
        dz = boxsize[2] / res_z

        del data
        
        # test field
        #x1d = np.linspace(0,boxsize[0].value,resolution)*boxsize[0].units
        #y1d = np.linspace(0,boxsize[1].value,resolution)*boxsize[1].units
        #z1d = np.linspace(0,boxsize[2].value,resolution)*boxsize[2].units
        #X, Y, Z = np.meshgrid(x1d, y1d, z1d, indexing="ij")
        #kb0 = 1.0
        #Bx_cube = (np.sin(kb0 * Y) + np.sin(kb0 * Z))*Bx_cube.units
        #By_cube = (np.cos(kb0 * X) - np.cos(kb0 * Z))*By_cube.units
        #Bz_cube = (np.sin(kb0 * X) + np.cos(kb0 * Y))*Bz_cube.units


        # Get vector potential
        dV = dx * dy * dz
        Bx_k = np.fft.fftn(Bx_cube.value)*Bx_cube.units * dV
        By_k = np.fft.fftn(By_cube.value)*By_cube.units * dV
        Bz_k = np.fft.fftn(Bz_cube.value)*Bz_cube.units * dV
        kx1d = 2*np.pi * np.fft.fftfreq(res_xy, d=dx)
        ky1d = 2*np.pi * np.fft.fftfreq(res_xy, d=dy)
        kz1d = 2*np.pi * np.fft.fftfreq(res_z, d=dz)
        kx, ky, kz = np.meshgrid(kx1d, ky1d, kz1d, indexing="ij")
        k2 = kx**2 + ky**2 + kz**2
        k2[0, 0, 0] = 1.0
        # projection to divergence-free field
        k_dot_B = kx*Bx_k + ky*By_k + kz*Bz_k
        Bx_k -= kx * k_dot_B / k2
        By_k -= ky * k_dot_B / k2
        Bz_k -= kz * k_dot_B / k2
        #Bx_cube = np.fft.ifftn(Bx_k.value).real*Bx_k.units / dV
        #By_cube = np.fft.ifftn(By_k.value).real*By_k.units / dV
        #Bz_cube = np.fft.ifftn(Bz_k.value).real*Bz_k.units / dV
        Ax_k = 1j * (ky * Bz_k - kz * By_k) / k2
        Ay_k = 1j * (kz * Bx_k - kx * Bz_k) / k2
        Az_k = 1j * (kx * By_k - ky * Bx_k) / k2
        Ax_k[0,0,0] = 0
        Ay_k[0,0,0] = 0
        Az_k[0,0,0] = 0
        Ax_cube = np.fft.ifftn(Ax_k.value).real*Ax_k.units / dV
        Ay_cube = np.fft.ifftn(Ay_k.value).real*Ay_k.units / dV
        Az_cube = np.fft.ifftn(Az_k.value).real*Az_k.units / dV
        del Ax_k,Ay_k,Az_k
        Jx_k = 1j * (ky * Bz_k - kz * By_k)
        Jy_k = 1j * (kz * Bx_k - kx * Bz_k)
        Jz_k = 1j * (kx * By_k - ky * Bx_k)
        Jx_cube = np.fft.ifftn(Jx_k.value).real*Jx_k.units / dV
        Jy_cube = np.fft.ifftn(Jy_k.value).real*Jy_k.units / dV
        Jz_cube = np.fft.ifftn(Jz_k.value).real*Jz_k.units / dV
        del Jx_k,Jy_k,Jz_k

        vx_k = np.fft.fftn(vx_cube.value) * vx_cube.units * dV
        vy_k = np.fft.fftn(vy_cube.value) * vy_cube.units * dV
        vz_k = np.fft.fftn(vz_cube.value) * vz_cube.units * dV
        k_dot_v = kx*vx_k + ky*vy_k + kz*vz_k
        vx_k -= kx * k_dot_v / k2
        vy_k -= ky * k_dot_v / k2
        vz_k -= kz * k_dot_v / k2
        vx_cube = np.fft.ifftn(vx_k.value).real * vx_k.units / dV
        vy_cube = np.fft.ifftn(vy_k.value).real * vy_k.units / dV
        vz_cube = np.fft.ifftn(vz_k.value).real * vz_k.units / dV
        del vx_k,vy_k,vz_k

        dBxdt_k = np.fft.fftn(dBxdt_cube.value)*dBxdt_cube.units * dV
        dBydt_k = np.fft.fftn(dBydt_cube.value)*dBydt_cube.units * dV
        dBzdt_k = np.fft.fftn(dBzdt_cube.value)*dBzdt_cube.units * dV

        # projection to divergence-free field
        k_dot_dBdt = kx*dBxdt_k + ky*dBydt_k + kz*dBzdt_k
        dBxdt_k -= kx * k_dot_dBdt / k2
        dBydt_k -= ky * k_dot_dBdt / k2
        dBzdt_k -= kz * k_dot_dBdt / k2
        #dBxdt_cube = np.fft.ifftn(dBxdt_k.value).real*dBxdt_k.units / dV
        #dBydt_cube = np.fft.ifftn(dBydt_k.value).real*dBydt_k.units / dV
        #dBzdt_cube = np.fft.ifftn(dBzdt_k.value).real*dBzdt_k.units / dV
        del dBxdt_k,dBydt_k,dBzdt_k
        del Bx_k,By_k,Bz_k,kx,ky,kz
        B = np.stack((Bx_cube,By_cube,Bz_cube),axis=-1)
        del Bx_cube,By_cube,Bz_cube
        A = np.stack((Ax_cube,Ay_cube,Az_cube),axis=-1)
        del Ax_cube,Ay_cube,Az_cube
        J = np.stack((Jx_cube,Jy_cube,Jz_cube),axis=-1)
        del Jx_cube,Jy_cube,Jz_cube
        v = np.stack((vx_cube,vy_cube,vz_cube),axis=-1)
        del vx_cube,vy_cube,vz_cube
        dBdt = np.stack((dBxdt_cube,dBydt_cube,dBzdt_cube),axis=-1)
        del dBxdt_cube,dBydt_cube,dBzdt_cube

        B2 = np.sum(B * B, axis=-1)
        J2 = np.sum(J * J, axis=-1)
        v2 = np.sum(v * v, axis=-1)
        A2 = np.sum(A * A, axis=-1)

        # plane averaged B
        B_bar = np.mean(B,axis = (0,1))
        B_bar2 = np.sum(B_bar * B_bar, axis=-1)

        # Rms values
        Brms = np.sqrt(np.mean(B2))
        Jrms = np.sqrt(np.mean(J2))
        vrms = np.sqrt(np.mean(v2))
        B_bar_rms = np.sqrt(np.mean(B_bar2))
        Arms = np.sqrt(np.mean(A2))

        # dot products (vectorized, no loop)
        JdotB = np.sum(J * B, axis=-1)
        AdotB = np.sum(A * B, axis=-1)
        BdotdBdt = np.sum(B * dBdt, axis=-1)

        # JdotB rms
        JdotB_rms = np.sqrt(np.mean(JdotB*JdotB))

        # cross product (computed once, vectorized)
        JcrossB = np.cross(J, B)

        # calculate RMS growth rate
        growth_rate_Brms = np.mean(BdotdBdt)/Brms**2

        # rms ratio
        B_bar_rms_Norm = B_bar_rms/Brms

        # WL = v · (J × B)
        WL = np.sum(v * JcrossB, axis=-1)

        # final quantities
        JcdotB_avrgNorm = np.mean(JdotB) / (Brms * Jrms)
        JcdotB2_avrgNorm = np.sqrt(np.mean(JdotB**2)) / (Brms * Jrms)

        JcrossB2_avrgNorm = np.sqrt(np.mean(np.sum(JcrossB**2, axis=-1))) / (Brms * Jrms)

        WL_avrg_norm = np.mean(WL) / (vrms * Jrms * Brms)

        J2B2_avrgNorm = np.mean(J2 * B2) / (Brms**2 * Jrms**2)
        AdotB_avrg_Norm = np.mean(AdotB) / (Brms * Brms)

        Jrms_over_Brms = Jrms/Brms
        Jrms_over_Brmsk=Jrms_over_Brms/0.5
        Arms_over_Brms = Arms/Brms
        JdotB_avrg_Norm = np.mean(JdotB) / (Brms * Brms)
        AdotB_avrg_Norm = np.mean(AdotB) / (Brms * Brms)
       
        print('Jrms_over_Brms',Jrms_over_Brms)
        print('Arms_over_Brms',Arms_over_Brms)
        #print('JdotB_avrg_Norm',JdotB_avrg_Norm)
        #print('AdotB_avrg_Norm',AdotB_avrg_Norm)
        #print(WL_avrg_norm)
        #print(J2B2_avrgNorm)
        #print('AdotB_over_Brms2',AdotB_avrg_Norm)
 
        # load statistics
        filename = run_directory + str(run_data_slice["Run #"].values[0]) + '/statistics.txt'
        tau_max=40
        take_last = int(20 / (5e-2))
        the_statistics = np.loadtxt(filename, usecols=range(42)).T
        Time = np.array(the_statistics[1])
        B = np.array(the_statistics[39])
        B = B / B[0]
        divB = np.abs(np.array(the_statistics[35]))

        mask = Time <= tau_max
        Time = Time[mask]
        B = B[mask]
        divB = divB[mask]

        dBdt = np.diff(B) / np.diff(Time)
        cut_B = B[1:].copy()
        local_growth_rate = dBdt[-take_last:] / cut_B[-take_last:]
        cut_timestamps = Time[-take_last:]

        growth_rate = np.mean(local_growth_rate)
        growth_rate_err_2sigma = 2*np.std(local_growth_rate)
        growth_rate_str = str(np.round(growth_rate,5))+"$\pm$"+str(np.round(growth_rate_err_2sigma,5))
       
        results = {'growth_rate_Brms_stat':np.array(growth_rate),'growth_rate_Brms_snap':growth_rate_Brms.value,'B_bar_rms_Norm':B_bar_rms_Norm.value,'Jrms_over_Brmsk':Jrms_over_Brmsk.value,'JcdotB_avrgNorm': JcdotB_avrgNorm.value, 'JcdotB2_avrgNorm':JcdotB2_avrgNorm.value, 'JcrossB2_avrgNorm':JcrossB2_avrgNorm.value,'WL_avrg_norm':WL_avrg_norm.value,'J2B2_avrgNorm':J2B2_avrgNorm.value, 'AdotB_avrg_Norm':AdotB_avrg_Norm.value, 'B_bar_rms':(B_bar_rms/Brms).value, 'WL':np.mean(WL/(vrms * Brms**2)).value, 'Jrms':(Jrms/Brms).value, 'JdotBrms':(JdotB_rms / Brms **2).value, 'JdotB_mean':(np.mean(JdotB)/ Brms**2).value}
        resulting_list.append(results)

   results_df = pd.DataFrame(resulting_list)
   for key in results_df.keys():
       run_data[key] = results_df[key]
   #run_data['JcdotB_avrgNorm'] = results_df['JcdotB_avrgNorm']
   run_data.to_csv(outputfile,sep="\t")
   return results_df

def plot_results(results_df):

   keys_to_plot = ['growth_rate_Brms_stat','B_bar_rms','WL','Jrms','JdotBrms','JdotB_mean']
   key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}B_{rms}^2 \right)$',r'$ \rm err \left(  J_{rms}/B_{rms} \right)$',r'$\rm err \left( (\vec J \cdot \vec B)_{rms}/B_{rms}^2 \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle \right)/B_{rms}^2$'] 

   #keys_to_plot = ['growth_rate_Brms_stat','B_bar_rms_Norm','WL_avrg_norm','JcdotB_avrgNorm']
   #key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}J_{rms}B_{rms} \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle / J_{rms}B_{rms} \right)$'] 
   Nkeys = len(keys_to_plot)
   ncol = 2
   nrow = int(Nkeys/ncol)
   fig, ax = plt.subplots(nrow, ncol, sharey=True,sharex=True, figsize=(5 * ncol, 5 * nrow))
   # load pencil
   resulting_Pencil = pd.read_csv('Pencil_2nd_results.csv',sep='\t')
   for i in range(nrow):
       for j in range(ncol):
            index = i*ncol+j
            if index<Nkeys:
                column_name = keys_to_plot[index]
                col = results_df[column_name]
                col_P = resulting_Pencil[column_name]
                quantity = np.array([x.item() for x in results_df[column_name]])
                quantity -= quantity[-1]
                Narray = 16*2**np.arange(len(quantity)) 
                #print(Narray)
                ax[i][j].scatter(Narray[:-1],np.abs(quantity)[:-1],marker='x',s=70,color='black', label='SWIFT')
                quantity_P = resulting_Pencil[column_name].to_numpy()
                #quantity_P = np.array([x.item() for x in resulting_Pencil[column_name]]) 
                N_P = resulting_Pencil['Resolution'].to_numpy()
                #N_P = np.array([x.item() for x in resulting_Pencil['Resolution']])
                ax[i][j].scatter(N_P,quantity_P,marker='+',s=100,color='black', label='Pencil 2nd')
                ax[i][j].set_yscale('log')
                ax[i][j].set_xscale('log')
                ax[i][j].set_ylabel(key_names[index],fontsize=20)
                ax[i][j].set_ylim([1e-5,1e0])
                ax[i][j].set_xticks(N_P,N_P,fontsize=20)
                ax[i][j].set_yticks([1e0,1e-1,1e-2,1e-3,1e-4,1e-5],['$10^0$','$10^{-1}$','$10^{-2}$','$10^{-3}$','$10^{-4}$','$10^{-5}$'],fontsize=20)
                ax[i][j].xaxis.set_minor_locator(mticker.NullLocator())
       ax[0][0].legend(fontsize=20)

   plt.savefig('test.png')

#run_data=load_test_run_parameters()[:]
#resulting_list = process_info(run_data, snapshot_name = 'RobertsFlow_0002.hdf5', outputfile='snapshot_results_FT.csv')
#plot_results(resulting_list)
