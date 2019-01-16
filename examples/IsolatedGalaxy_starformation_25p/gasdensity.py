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

def getdensfrtemp(filename):
    
    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates= f["/PartType0/Coordinates"][:,:]
        SFR = f["/PartType0/SFR"][:]
        density = f["/PartType0/Density"][:]
        temperature = f["/PartType0/Temperature"][:]

    absmaxz = 2 #kpc
    absmaxxy = 10 #kpc

    part_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz) 
    SFR = SFR[part_mask]
    density = density[part_mask]
    temperature = temperature[part_mask]

    return density, SFR, temperature

def getbirthdensity(filename):
    
    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates= f["/PartType4/Coordinates"][:,:]
        birthdensity = f["/PartType4/BirthDensity"][:]

    absmaxz = 2 #kpc
    absmaxxy = 10 #kpc

    part_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz) & (birthdensity>0)

    birthdensity = birthdensity[part_mask]

    return birthdensity 
    
def KS_law(Sigmagas,n=1.4,A=1.515e-4):
    return A * Sigmagas**n

def makefullplot(filename,binsize=200):
    # Get the Surface density from the snapshot
    image_data,extent = gen_data(filename,binsize,case=0)
    image_data2 = image_data[:,:]*1e10/(1e2**2)
    surface_density = image_data[:,:]*1e10/(1e2**2)
    if binsize==20:
        image_data2 /= 1e2
        surface_density /= 1e2
    minvalue = np.min(image_data2[image_data2!=0])
    image_data2[image_data2==0] = minvalue*.9
    
    # Get the SFR surface density from the snapshot
    image_data,extent = gen_data(filename,binsize,case=2)
    image_data3 = image_data[:,:] * 1e10/1e9/(.1**2) 
    SFR_smoothed = image_data3/image_data2
    SFR_surface_density = image_data[:,:] * 1e10/1e9/(.1**2)/surface_density
    if binsize==20:
        image_data3 /= 1e2
        SFR_surface_density /= 1e2
        SFR_smoothed /= 1e2
    plot_SFR_smoothed = np.ones(np.shape(image_data3))*np.min(SFR_smoothed[SFR_smoothed!=0])*.9
    plot_SFR_smoothed[image_data3!=0] = SFR_smoothed[image_data3!=0]
    
    # Define the font used in the plots 
    font = {'color':'white', 'size':12}

    # Make the first plot of the Gas density
    fig = plt.figure(1,figsize=(16.5,16.5))
    plt.subplot(3,3,1)
    plt.imshow(np.log10(image_data2), extent=np.array([-10,10,-10,10]), cmap='viridis') #,vmin=-3.5, vmax=0.25)
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
    plt.clim(-1.5,3.0)
    #cbar.set_label('Surface density ($M_\odot pc^{-2}$)', rotation=270)
    
    # Make the second plot of the SFR density
    plt.subplot(3,3,2)
    plt.imshow(np.log10(plot_SFR_smoothed), extent=np.array([-10,10,-10,10]), cmap='viridis') #,vmin=-3.5, vmax=0.25)
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
    plt.clim(-5,-1)

    # Make the third plot of the stars that are formed
    plt.subplot(3,3,3)
    xnewstar, ynewstar = getnewstars(filename)
    newstarmatrix = np.histogram2d(xnewstar,ynewstar,bins=binsize)[0]
    newstarmatrix[newstarmatrix==0] = 1e-1
    plt.imshow(np.log10(newstarmatrix), extent=np.array([-10,10,-10,10]), cmap='viridis') #,vmin=-3.5, vmax=0.25)
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
    plt.clim(-1,3)
    
    # Make the fourth plot of the KS law
    plt.subplot(3,3,4)
    sigma_gas = image_data2.flatten()
    sigma_star = SFR_smoothed.flatten()
    sigma_gas=sigma_gas[sigma_star>0]
    sigma_star=sigma_star[sigma_star>0]
    sigma_gas2 = surface_density.flatten()
    sigma_star2 = SFR_surface_density.flatten()

    #KShistogram = np.histogram2d(np.log10(sigma_gas),np.log10(sigma_star),bins=50,range=[[-1,4],[-5,2]])
    #plt.imshow(np.transpose(KShistogram[0]), extent=np.array([-1,4,-5,2]))
    plt.title('KS law')

    sigma_gas_range = 10**np.linspace(-1,2,100)
    sigma_star_range = KS_law(sigma_gas_range)
    sigma_star_one_range = KS_law(sigma_gas_range,n=0.7)

    if binsize==200:
        plt.hist2d(np.log10(sigma_gas),np.log10(sigma_star),range=[[-1,2],[-5,-2]],bins=50)
    else:
        plt.hist2d(np.log10(sigma_gas),np.log10(sigma_star),range=[[-1,2],[-5,-2]],bins=25)
    #plt.hist2d(np.log10(sigma_gas2),np.log10(sigma_star2),range=[[-1,2],[-5,-2]],bins=50)
    plt.plot(np.log10(sigma_gas_range), np.log10(sigma_star_range),'k--',label='EAGLE SF')
    plt.plot(np.log10(sigma_gas_range), np.log10(sigma_star_one_range),'k-.',label='n=.7')
    plt.plot(np.log10(sigma_gas_range), np.log10(KS_law(sigma_gas_range,n=0.6)),'k:',label='n=.6')
    plt.xlabel('$\log \Sigma_{gas}$ [ $M_\odot \\rm pc^{-2}$]')
    plt.ylabel('$\log \Sigma_{\star}$ [ $M_\odot \\rm yr^{-1} kpc^{-2}$]')
    cbar = plt.colorbar()
    if binsize==200:
        plt.clim(0,700)
    elif binsize==20:
        plt.clim(0,40)
    #plt.scatter(np.log10(sigma_gas),np.log10(sigma_star))
    plt.legend()

    # make a plot of the gas density vs SFR
    plt.subplot(3,3,5)
    critden = .1*((0.0129)/.002)**(-.64)
    mu = 1.3
    density, SFR, temperature = getdensfrtemp(filename)
    density *= 1e10/(1e3)**3*40.4759/mu
    SFR *= 1e10/1e9
    SFR_plot = SFR[:]
    plt.hist2d(np.log10(density[SFR_plot>0]),np.log10(SFR_plot[SFR_plot>0]), bins=50,range=[[np.log10(critden),5],[-5.4,-3.5]])
    #plt.loglog(density*1e10/(1e3)**3*40.4759/mu,SFR_plot*1e10/1e9,'.')
    xx = 10**np.linspace(-1,2,100)
    def sfr_model(density,critden=.1,gammaeff=4./3.,mgas=1.4995e4):
        return 6.9e-11 * mgas * (density/critden)**(gammaeff/5.)
    plt.plot(np.log10(xx),np.log10(1e-5*xx**(4./15.)),'k--',label='Slope of EAGLE SF')
    plt.plot(np.log10(xx),np.log10(sfr_model(xx)),'r--',label='Slope of EAGLE SF')
    plt.xlabel('$\log$ density [$\\rm cm^{-3}$]')
    plt.ylabel('$\log$ Star Formation rate (SFR) [$ M_\odot \\rm yr^{-1} $]')
    plt.xlim(np.log10(critden),5)
    cbar = plt.colorbar()
    plt.legend()

    SFR_tracer = np.log10(SFR_plot[SFR_plot>0])
    
    plt.subplot(3,3,6)
    plt.scatter(np.log10(density[SFR_plot>0]),np.log10(temperature[SFR_plot>0]),s=1,c=SFR_tracer,cmap='viridis')
    plt.axvline(x=np.log10(critden),linestyle='--',color='gray',linewidth=3,label='Threshold density')
    cbar = plt.colorbar()
    plt.clim(-5.4,-3.5)
    plt.scatter(np.log10(density[SFR_plot<=0]),np.log10(temperature[SFR_plot<=0]),s=1,c='gray')
    plt.xlabel('$\log$ density [$\\rm cm^{-3}$]')
    plt.ylabel('$\log$ temperature [$\\rm K$]')
    plt.ylim(2,4)
    plt.xlim(-3.5,5)
    plt.legend()

    plt.subplot(3,3,7)
    birthdensity = getbirthdensity(filename)*404.759
    plt.hist(np.log10(birthdensity),bins=50,density=True,range=(np.log10(critden),5))
    plt.axvline(x=np.log10(critden),linestyle='--',color='gray',linewidth=3,label='Threshold density')
    plt.xlabel('$\log$ density [$\\rm cm^{-3}$]')
    plt.xlim(np.log10(critden),5)
    plt.legend()

    plt.subplot(3,3,8)
    plt.hist(np.log10(sigma_gas))
    plt.xlabel('Surface density')
    
    plt.subplot(3,3,9)
    plt.hist(np.log10(sigma_star))
    plt.xlabel('SFR Surface density')

    if binsize==200:
        plt.savefig('./ksplot/'+filename+'_0.1kpc.png')
    elif binsize==20:
        plt.savefig('./ksplot_1kpc/'+filename+'_0.1kpc.png')
    plt.close()


def getmassandsfr(filename):
    
    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates= f["/PartType0/Coordinates"][:,:]
        SFR = f["/PartType0/SFR"][:]
        coordinates_star = f["/PartType4/Coordinates"][:,:]
        masses_star = f["/PartType4/Masses"][:]
        

    absmaxz = 2 #kpc
    absmaxxy = 10 #kpc

    part_mask = ((coordinates[:,0]-box_size/2.) > -absmaxxy) & ((coordinates[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,1]-box_size/2.) > -absmaxxy) & ((coordinates[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates[:,2]-box_size/2.) > -absmaxz) & ((coordinates[:,2]-box_size/2.) < 
        absmaxz) & (SFR>0) 
    mask = ((coordinates_star[:,0]-box_size/2.) > -absmaxxy) & ((coordinates_star[:,0]-box_size/2.) < 
        absmaxxy) & ((coordinates_star[:,1]-box_size/2.) > -absmaxxy) & ((coordinates_star[:,1]-box_size/2.) < 
        absmaxxy) & ((coordinates_star[:,2]-box_size/2.) > -absmaxz) & ((coordinates_star[:,2]-box_size/2.) < 
        absmaxz)  

    SFR_snap = np.sum(SFR[part_mask])
    mass_star = np.sum(masses_star[mask])
    return SFR_snap, mass_star

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
    #for i in [1,5,10,33,66,90]:
    #    print('output_%04d.hdf5'%i)
    #    makefullplot('output_%04d.hdf5'%i)
    #    makefullplot('output_%04d.hdf5'%i,binsize=20)

    '''
    snapshots = 52 
    mass = np.zeros(snapshots)
    sfr = np.zeros(snapshots)
    for i in range(1,snapshots):
        if i!=6:
            print(i)
            sfr[i], mass[i] = getmassandsfr('output_%04d.hdf5'%i)
        
    plt.plot(np.log10(mass*1e10))
    plt.xlabel('Snapshot number (Myr)')
    plt.ylabel('Mass')
    plt.savefig('Mass increase')
    plt.close()

    plt.plot(np.log10(sfr*10))
    plt.xlabel('Snapshot number (Myr)')
    plt.ylabel('SFR ($M_\odot / \\rm yr$)')
    plt.savefig('SFR')
    plt.close()
    '''

    for i in [1,5,10,33,66,90]:
        if i!=6:
            print('output_%04d.hdf5'%i)
            makefullplot('output_%04d.hdf5'%i)
            makefullplot('output_%04d.hdf5'%i,binsize=20)


