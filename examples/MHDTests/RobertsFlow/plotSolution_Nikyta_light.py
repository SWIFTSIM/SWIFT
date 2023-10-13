from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import numpy as np
import sys
import h5py 

filename = sys.argv[1]

with h5py.File(filename, "r") as handle:
    #print(handle["PartType0"].keys())
    #print(handle["PartType0/ViscosityParameters"][0]) 
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"][0]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
    dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    #mhdeta = handle["/HydroScheme"].attrs["Diffusion Eta"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]




filename = sys.argv[1]
data = load(filename)

#getting values from snapshots

img_res=1000
project="xy"

#print(data.metadata.gas_properties.field_names)
x=data.gas.coordinates[:,0]
y=data.gas.coordinates[:,1]
z=data.gas.coordinates[:,2]
rho=data.gas.densities.value
h=data.gas.smoothing_lengths.value
v=data.gas.velocities.value
P=data.gas.pressures.value
B=data.gas.magnetic_flux_densities.value
#divB=data.gas.magnetic_divergence.value
#curlB=data.gas.magnetic_flux_curl.value
#Fmag=data.gas.fmag.value
#Ftot=data.gas.ftot.value
#Ftot=Ftot-Fmag
#SPH_1=data.gas.sph_1.value
#SPH_grad_1=data.gas.sph_grad_1.value
#dt_div_v=data.gas.velocity_divergence_time_differentials.value
#visc_par=data.gas.viscosity_parameters.value




data.gas.mass_weighted_x=(
    data.gas.masses * x[:]
)
data.gas.mass_weighted_y=(
    data.gas.masses * y[:]
)
data.gas.mass_weighted_z=(
    data.gas.masses * z[:]
)
data.gas.mass_weighted_rho=(
    data.gas.masses * rho[:]
)
data.gas.mass_weighted_h=(
    data.gas.masses * h[:]
)
data.gas.mass_weighted_vx=(
    data.gas.masses * v[:,0]
)
data.gas.mass_weighted_vy=(
    data.gas.masses * v[:,1]
)
data.gas.mass_weighted_vz=(
    data.gas.masses * v[:,2]
)
data.gas.mass_weighted_P=(
    data.gas.masses * P[:]
)
data.gas.mass_weighted_Bx=(
    data.gas.masses * B[:,0]
)
data.gas.mass_weighted_By=(
    data.gas.masses * B[:,1]
)
data.gas.mass_weighted_Bz=(
    data.gas.masses * B[:,2]
)
#data.gas.mass_weighted_divB=(
#    data.gas.masses * divB[:]
#)
#data.gas.mass_weighted_curlBx=(
#    data.gas.masses * curlB[:,0]
#)
#data.gas.mass_weighted_curlBy=(
#    data.gas.masses * curlB[:,1]
#)
#data.gas.mass_weighted_curlBz=(
#    data.gas.masses * curlB[:,2]
#)

#data.gas.mass_weighted_Fmagx=(
#    data.gas.masses * Fmag[:,0]
#)
#data.gas.mass_weighted_Fmagy=(
#    data.gas.masses * Fmag[:,1]
#)
#data.gas.mass_weighted_Fmagz=(
#    data.gas.masses * Fmag[:,2]
#)
#data.gas.mass_weighted_Ftotx=(
#    data.gas.masses * Ftot[:,0]
#)
#data.gas.mass_weighted_Ftoty=(
#    data.gas.masses * Ftot[:,1]
#)
#data.gas.mass_weighted_Ftotz=(
#    data.gas.masses * Ftot[:,2]
#)
#data.gas.mass_weighted_SPH_1=(
#    data.gas.masses * SPH_1[:]
#)
#data.gas.mass_weighted_SPH_grad_1x=(
#    data.gas.masses * SPH_grad_1[:,0]
#)
#data.gas.mass_weighted_SPH_grad_1y=(
#    data.gas.masses * SPH_grad_1[:,1]
#)
#data.gas.mass_weighted_SPH_grad_1z=(
#    data.gas.masses * SPH_grad_1[:,2]
#)
#data.gas.mass_weighted_dt_div_v=(
#    data.gas.masses * dt_div_v[:]
#    )
#data.gas.mass_weighted_visc_par=(
#    data.gas.masses * visc_par[:]
#        )

#center = 0.5 * data.metadata.boxsize

#if project=="xy":
#    rotate_vec = [0,0,1]
#elif project=="yz":
#    rotate_vec = [1,0,0]
#elif project=="zx":
#    rotate_vec = [0,1,0]

#matrix = rotation_matrix_from_vector(rotate_vec, axis='z')

def make_slice(key):
    res=slice_gas(data,z_slice=0.01 * data.metadata.boxsize[2],resolution=img_res,project=key,parallel=True)
    return res

divreg=1e-30

masses=make_slice('masses').value
dimy=len(masses)
dimx=len(masses[0])
#print(dimx,dimy)
mass_map = masses.flatten()
mass_map = mass_map + divreg

#print(mass_map)
l=len(mass_map)
#L=np.sqrt(lmm)

v=np.zeros((l,3))
B=np.zeros((l,3))
#curlB=np.zeros((l,3))
#Fmag=np.zeros((l,3))
#Ftot=np.zeros((l,3))
#SPH_grad_1=np.zeros((l,3))

rho=make_slice('mass_weighted_rho').value.flatten()/mass_map
h=make_slice('mass_weighted_h').value.flatten()/mass_map
v[:,0]=make_slice('mass_weighted_vx').value.flatten()/mass_map
v[:,1]=make_slice('mass_weighted_vy').value.flatten()/mass_map
v[:,2]=make_slice('mass_weighted_vz').value.flatten()/mass_map
P=make_slice('mass_weighted_P').value.flatten()/mass_map
B[:,0]=make_slice('mass_weighted_Bx').value.flatten()/mass_map
B[:,1]=make_slice('mass_weighted_By').value.flatten()/mass_map
B[:,2]=make_slice('mass_weighted_Bz').value.flatten()/mass_map
#divB=make_slice('mass_weighted_divB').value.flatten()/mass_map
#curlB[:,0]=make_slice('mass_weighted_curlBx').value.flatten()/mass_map
#curlB[:,1]=make_slice('mass_weighted_curlBy').value.flatten()/mass_map
#curlB[:,2]=make_slice('mass_weighted_curlBz').value.flatten()/mass_map
#Fmag[:,0]=make_slice('mass_weighted_Fmagx').value.flatten()/mass_map
#Fmag[:,1]=make_slice('mass_weighted_Fmagy').value.flatten()/mass_map
#Fmag[:,2]=make_slice('mass_weighted_Fmagz').value.flatten()/mass_map
#Ftot[:,0]=make_slice('mass_weighted_Ftotx').value.flatten()/mass_map
#Ftot[:,1]=make_slice('mass_weighted_Ftoty').value.flatten()/mass_map
#Ftot[:,2]=make_slice('mass_weighted_Ftotz').value.flatten()/mass_map

#SPH_1=make_slice('mass_weighted_SPH_1').value.flatten()/mass_map
#SPH_grad_1[:,0]=make_slice('mass_weighted_SPH_grad_1x').value.flatten()/mass_map
#SPH_grad_1[:,1]=make_slice('mass_weighted_SPH_grad_1y').value.flatten()/mass_map
#SPH_grad_1[:,2]=make_slice('mass_weighted_SPH_grad_1z').value.flatten()/mass_map
#dt_div_v=make_slice('mass_weighted_dt_div_v').value.flatten()/mass_map
#visc_par=make_slice('mass_weighted_visc_par').value.flatten()/mass_map


bb=np.sqrt(B[:,0]**2+B[:,1]**2+B[:,2]**2)
Pmag=(B[:,0]**2+B[:,1]**2+B[:,2]**2)/2
vv = np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2)
#cbcb= np.sqrt(curlB[:,0]**2+curlB[:,1]**2+curlB[:,2]**2)
#ff = np.sqrt(Fmag[:,0]**2+Fmag[:,1]**2+Fmag[:,2]**2)
#fftot = np.sqrt(Ftot[:,0]**2+Ftot[:,1]**2+Ftot[:,2]**2)
#RF=ff/(fftot+ff+divreg)
#cosfftot = (Fmag[:,0]*Ftot[:,0]+Fmag[:,1]*Ftot[:,1]+Fmag[:,2]*Ftot[:,2])/(ff*fftot+divreg)

#fb = Fmag[:,0]*B[:,0]+Fmag[:,1]*B[:,1]+Fmag[:,2]*B[:,2]
#cosfb=fb/(ff*bb+divreg)

#local_flux=bb/(rho**(2/3)+divreg)


#x -= data.metadata.boxsize[0] * 0.5
#y -= data.metadata.boxsize[1] * 0.5
#z -= data.metadata.boxsize[2] * 0.5

reg_err=1.001*1e-2
err_upper_bnd=0.9
above_noise=10

#Estimating errors of divB,curlB,Fmag
#SPH_divB_err = B[:,0]*SPH_grad_1[:,0]+B[:,1]*SPH_grad_1[:,1]+B[:,2]*SPH_grad_1[:,2]
#SPH_curlB_err = np.zeros((len(B),3))
#SPH_curlB_err[:,0] = SPH_grad_1[:,1]*B[:,2]-B[:,2]*SPH_grad_1[:,1]
#SPH_curlB_err[:,1] = SPH_grad_1[:,2]*B[:,0]-B[:,0]*SPH_grad_1[:,2]
#SPH_curlB_err[:,2] = SPH_grad_1[:,0]*B[:,1]-B[:,1]*SPH_grad_1[:,0]
#SPH_abs_curlB_err = np.sqrt(SPH_curlB_err[:,0]**2+SPH_curlB_err[:,1]**2+SPH_curlB_err[:,2]**2)
#SPH_Fmag_err = np.zeros((len(B),3))
#SPH_Fmag_err[:,0] = np.abs(bb*bb/mu0/(rho+divreg)*SPH_grad_1[:,0])+np.abs(SPH_1[:]*Fmag[:,0]) # -corr1-corr2
#SPH_Fmag_err[:,1] = np.abs(bb*bb/mu0/(rho+divreg)*SPH_grad_1[:,1])+np.abs(SPH_1[:]*Fmag[:,1])
#SPH_Fmag_err[:,2] = np.abs(bb*bb/mu0/(rho+divreg)*SPH_grad_1[:,2])+np.abs(SPH_1[:]*Fmag[:,2])
#SPH_abs_Fmag_err = np.sqrt(SPH_Fmag_err[:,0]**2+SPH_Fmag_err[:,1]**2+SPH_Fmag_err[:,2]**2)

#Calculating error metrics
#R0 = np.abs(divB*h)/(bb+np.abs(divB*h)+divreg)
#R1 = np.abs(cosfb+divreg)
#R2 = np.abs(divB)/(cbcb+np.abs(divB)+divreg)

#Condition for quantities to be above noise
#mask_R0 = (np.abs(divB)<above_noise*np.abs(SPH_divB_err))
#mask_R1 = ((ff<above_noise*SPH_abs_Fmag_err) | (RF<0.1))
#mask_R2 = ((mask_R0) | (cbcb<SPH_abs_curlB_err))
#R0[mask_R0]=0.0
#R1[mask_R1]=0.0
#R2[mask_R2]=0.0

#Setting up lower bound to display
#R0[R0<reg_err]=reg_err
#R1[R1<reg_err]=reg_err
#R2[R2<reg_err]=reg_err

#Setting up upper bound to display
#R0[R0>err_upper_bnd]=err_upper_bnd
#R1[R1>err_upper_bnd]=err_upper_bnd
#R2[R2>err_upper_bnd]=err_upper_bnd

#print(R0)
#print(R1)
#print(R2)

#print(len(rho))
from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import ticker


# plot everything
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(14, 10))

#print(np.min(x.value),np.max(x.value))
#print(np.min(y.value),np.max(y.value))
#print(np.min(z.value),np.max(z.value))

if project=="xy":
    new_x=np.linspace(np.min(x.value),np.max(x.value),dimx)
    new_y=np.linspace(np.min(y.value),np.max(y.value),dimy)
elif project=="yz":
    new_x=np.linspace(np.min(y.value),np.max(y.value),dimx)
    new_y=np.linspace(np.min(z.value),np.max(z.value),dimy)
elif project=="zx":
    new_x=np.linspace(np.min(z.value),np.max(z.value),dimx)
    new_y=np.linspace(np.min(x.value),np.max(x.value),dimy)

#print(len(x),len(y))

def make_color_levels(cmin,cmax,c_res=10,log_sc=True):
    cmax+=divreg
    if log_sc:
        levmin=int(np.floor(np.log10(cmin)))
        levmax=int(np.ceil(np.log10(cmax)))
        levels=[]
        levels_short=[]
        for i in range(levmin,levmax):
            levels_short+=[10**i]
            for j in range(int(c_res/10),c_res):
                levels+=[(10/c_res*j)*10**i]

    else:
        levels=[cmin+(cmax-cmin)/c_res*k for k in range(c_res)]
        levels_short=[cmin+(cmax-cmin)/c_res*k for k in range(0,c_res,10)]
    return levels, levels_short

def make_density_plot(Q,cmin,cmax,i,j,Q_name,c_res=10,log_sc=True,cmap='viridis'):
    levels,levels_short=make_color_levels(cmin,cmax,c_res,log_sc)
    if log_sc:
        to_plot=ax[i][j].contourf(new_x,new_y,Q.transpose(),levels=np.array(levels),locator=ticker.LogLocator(),cmap = cmap)
    else:
        to_plot=ax[i][j].contourf(new_x,new_y,Q.transpose(),levels=np.array(levels),cmap = cmap)
    fig.colorbar(to_plot,ticks=levels_short)
    ax[i][j].set_ylabel(Q_name)
    return 0

def make_slice_plot(Q,cmin,cmax,i,j,Q_name):
    slice_Q=Q[:,int(len(Q)/2)]
    ax[i][j].plot(x,slice_Q)
    ax[i][j].plot(x,max(slice_Q)*np.ones(len(slice_Q)),"--")
    ax[i][j].set_ylim([cmin,cmax])
    ax[i][j].set_yscale("log")
    ax[i][j].set_ylabel(Q_name)
    return 0
#print(len(Pmag.reshape((dimx,dimy))),len(Pmag.reshape((dimx,dimy))[0]))

Pmag=bb

make_density_plot(Pmag.reshape((dimx,dimy)),np.min(Pmag),1.1*np.max(Pmag),0,0,'B',c_res=100,log_sc=False)
make_density_plot(P.reshape((dimx,dimy)),np.min(P),1.1*np.max(P),1,0,'Pressure',c_res=100,log_sc=False)
make_density_plot(vv.reshape((dimx,dimy)),np.min(vv),1.1*np.max(vv),0,1,'Velocity',c_res=100,log_sc=False)
#make_density_plot(cosfftot.reshape((img_res,img_res)),0.75,1.01,0,1,'new',c_res=100,log_sc=False)
make_density_plot(rho.reshape((dimx,dimy)),np.min(rho),1.1*np.max(rho),1,1,'Density',c_res=100,log_sc=False)
#make_density_plot(R0.reshape((dimx,dimy)),1e-2,1,0,2,'R0',c_res=100)
#make_density_plot(R1.reshape((dimx,dimy)),1e-2,1,1,2,'R1',c_res=100)
#make_density_plot(R2.reshape((dimx,dimy)),1e-2,1,0,3,'R2',c_res=100)
#make_density_plot(np.abs(dt_div_v).reshape((dimx,dimy)),np.min(np.abs(dt_div_v)),1.1*np.max(np.abs(dt_div_v)),1,3,"d (div v)/dt",c_res=100,log_sc=False)
#make_density_plot(np.abs(visc_par).reshape((dimx,dimy)),np.min(np.abs(visc_par)),1.1*np.max(np.abs(visc_par)),0,4,"alpha_visc*B_switch",c_res=100,log_sc=False)
#make_density_plot(local_flux.reshape((dimx,dimy)),np.min(local_flux),1.1*np.max(local_flux),1,4,"B/rho^2/3",c_res=100,log_sc=False)
#make_slice_plot(Pmag.reshape((img_res,img_res)),np.min(Pmag),np.max(Pmag),0,2,'Pmag')
#make_slice_plot(R1.reshape((img_res,img_res)),1e-2,1,1,2,'R1')
#make_slice_plot(R2.reshape((img_res,img_res)),1e-2,1,2,2,'R2')


ax[0,1].set_title(f"t={t:.2e}")
ax[0,0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
#ax[0,2].set_title("y slice")
#ax[2,2].set_xlabel("y")
#ax[2,0].set_xlabel("x")
#ax[2,1].set_xlabel("x")
fig.tight_layout()

plt.savefig(sys.argv[2],dpi=300)

#fig, ax = plt.subplots(1, 2, figsize=(14, 4))
#ax[0].hist(cosfftot,bins=2000)
#ax[1].hist(RF,bins=2000)
#plt.savefig("hist.png",dpi=300)
