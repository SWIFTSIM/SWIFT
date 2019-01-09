"""
Makes a movie using sphviewer and ffmpeg.

Edit low_frac and up_frac to focus on a certain view of the box.
The colour map can also be changed via colour_map.

Usage: python3 makeMovie.py CoolingHalo_

Written by James Willis (james.s.willis@durham.ac.uk)
"""

import glob
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool
from tqdm import tqdm
from pylab import sys
from sphviewer.tools import QuickView
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm

def gen_data(filename,pixels,case):
    """
    Data generation function. Reads a snapshot, filters the particle data and
    creates an image using QuickView from sphviewer.
    """

    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates = f["/PartType0/Coordinates"][:,:]
        masses = f["/PartType0/Masses"][:]
        smoothing_lengths = f["/PartType0/SmoothingLength"][:]
        coordinates_star= f["/PartType4/Coordinates"][:,:]
        star_flag = f["/PartType4/NewStarFlag"][:]
        SFR = f["/PartType0/SFR"][:]
    
    # Filter data  
    low_frac=0.48
    up_frac=0.52
    absmaxz = 2 #kpc
    absmaxxy = 10 #kpc
    part_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz)

    sfr_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz) & (SFR > 0 )
    
    masses = masses[part_mask]
    SFR = SFR[sfr_mask]
   
    print(np.sum(masses))
    if case ==0:
        coordinates = coordinates[part_mask]
        smoothing_lengths = smoothing_lengths[part_mask]
        masvar = masses
    elif case==1:
        coordinates = coordinates[part_mask]
        smoothing_lengths = smoothing_lengths[part_mask]
        masvar = masses**2
    elif case==2:
        coordinates = coordinates[sfr_mask]
        smoothing_lengths = smoothing_lengths[sfr_mask]
        masvar = SFR
    # Use sphviewer QuickView to generate an image
    qv = QuickView(
        coordinates,
        mass=masvar,
        r="infinity",
        xsize=pixels,
        ysize=pixels,
        p=0,  # Viewing angle theta
        roll=0,  # Viewing angle phi
        plot=False,
        logscale=False,
        hsml=smoothing_lengths
    )

    image_data = qv.get_image()
    extent = qv.get_extent()
   
    return image_data, extent
    
def getnewstars(filename):
    
    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates= f["/PartType4/Coordinates"][:,:]
        star_flag = f["/PartType4/NewStarFlag"][:]

    absmaxz = 2 #kpc
    absmaxxy = 10 #kpc

    part_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz) & (star_flag == 1)

    return coordinates[:,0][part_mask]-box_size/2., coordinates[:,1][part_mask]-box_size/2.

if __name__ == "__main__":

    # Some global variables to define the whole run
    input_file_name = sys.argv[1]
    pixels = 1000 
    num_procs = 4 # No. of processors to use
    colour_map = "inferno"
    
    # Find the list of files and sort them in numerical order
    file_list = sorted(glob.glob(input_file_name + "*.hdf5"), key=lambda x:
            int(x[len(input_file_name):len(input_file_name)+4]))
    
    # get the image data 
    filename = 'output_0066.hdf5' 

    # Get the Surface density from the snapshot
    image_data,extent = gen_data(filename,200,case=0)
    image_data2 = image_data[:,:]*1e10/(1e2**2)
    minvalue = np.min(image_data2[image_data2!=0])
    image_data2[image_data2==0] = minvalue*.9
    
    # Get the SFR surface density from the snapshot
    image_data,extent = gen_data(filename,200,case=2)
    image_data3 = image_data[:,:] * 1e10/1e9/(.1**2)
    SFR_smoothed = image_data3/image_data2
    plot_SFR_smoothed = np.ones(np.shape(image_data3))*np.min(SFR_smoothed[SFR_smoothed!=0])*.9
    plot_SFR_smoothed[image_data3!=0] = SFR_smoothed[image_data3!=0]
    
    # Define the font used in the plots 
    font = {'color':'white', 'size':12}

    # Make the first plot of the Gas density
    fig = plt.figure(1,figsize=(11,11))
    plt.subplot(2,2,1)
    plt.imshow(np.log10(image_data2), extent=extent, cmap='viridis') #,vmin=-3.5, vmax=0.25)
    plt.title('$\log \Sigma_{gas}$ [ $M_\odot \\rm pc^{-2}$]')
    plt.plot([-7.5,-2.5],[-9.0,-9.0],linewidth=3,color='white')
    plt.text(-7.5,-8.5,'5 kpc',fontdict=font)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off
    cbar = plt.colorbar() 
    #cbar.set_label('Surface density ($M_\odot pc^{-2}$)', rotation=270)
    
    # Make the second plot of the SFR density
    plt.subplot(2,2,2)
    plt.imshow(np.log10(plot_SFR_smoothed), extent=extent, cmap='viridis') #,vmin=-3.5, vmax=0.25)
    plt.title('$\log \Sigma_{\star}$ [ $M_\odot \\rm yr^{-1} kpc^{-2}$]')
    plt.plot([-7.5,-2.5],[-9.0,-9.0],linewidth=3,color='white')
    plt.text(-7.5,-8.5,'5 kpc',fontdict=font)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off
    cbar = plt.colorbar() 

    # Make the third plot of the stars that are formed
    plt.subplot(2,2,3)
    xnewstar, ynewstar = getnewstars(filename)
    newstarmatrix = np.histogram2d(xnewstar,ynewstar,bins=200)[0]
    newstarmatrix[newstarmatrix==0] = 1e-1
    plt.imshow(np.log10(newstarmatrix), extent=extent, cmap='viridis') #,vmin=-3.5, vmax=0.25)
    plt.title('New stars')
    plt.plot([-7.5,-2.5],[-9.0,-9.0],linewidth=3,color='white')
    plt.text(-7.5,-8.5,'5 kpc',fontdict=font)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelleft=False)   # labels along the bottom edge are off
    cbar = plt.colorbar()
    
    # Make the fourth plot of the KS law
    plt.subplot(2,2,4)
    sigma_gas = image_data2.flatten()
    sigma_star = SFR_smoothed.flatten()
    sigma_gas=sigma_gas[sigma_star>0]
    sigma_star=sigma_star[sigma_star>0]

    KShistogram = np.histogram2d(np.log10(sigma_gas),np.log10(sigma_star),bins=50,range=[[-1,4],[-5,2]])
    plt.imshow(np.transpose(KShistogram[0]), extent=np.array([-1,4,-5,2]))
    plt.title('KS law')
    cbar = plt.colorbar()
    cbar.set_ticks([])

    def KS_law(Sigmagas,n=1.4,A=1.515e-4):
        return A * Sigmagas**n

    sigma_gas_range = 10**np.linspace(-1,2,100)
    sigma_star_range = KS_law(sigma_gas_range)

    plt.hist2d(np.log10(sigma_gas),np.log10(sigma_star),range=[[-1,2],[-5,-2]],bins=50)
    plt.plot(np.log10(sigma_gas_range), np.log10(sigma_star_range),'k--',label='EAGLE SF')
    plt.xlabel('$\log \Sigma_{gas}$ [ $M_\odot \\rm pc^{-2}$]')
    plt.ylabel('$\log \Sigma_{\star}$ [ $M_\odot \\rm yr^{-1} kpc^{-2}$]')
    plt.legend()

    plt.savefig('./images/'+filename+'_0.1kpc.png')
    plt.close()
    
    image_data,extent = gen_data(filename,20,case=0)
    image_data2 = image_data[:,:]*1e10/(1e3**2)
    fig = plt.figure(1,figsize=(11,11))
    plt.imshow(np.log10(image_data2), extent=extent, cmap='viridis') #,vmin=-3.5, vmax=0.25)
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    cbar = plt.colorbar() 
    cbar.set_label('Surface density ($M_\odot~pc^{-2}$)', rotation=270)
    plt.savefig('./images/'+filename+'_1kpc.png')
