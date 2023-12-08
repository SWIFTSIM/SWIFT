from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import numpy as np
import sys
import h5py

filename = sys.argv[1]

slice_height = sys.argv[3]

projection = sys.argv[4]

prefered_color = 'magma'

cpts = 100


to_plot = 'errors'  # 'B' or 'A' or 'errors'

with h5py.File(filename, "r") as handle:
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"][0]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
    #dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    #dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]


data = load(filename)

# getting values from snapshots

img_res = 512
divreg = 1e-30
above_noise=10
reg_err=0.01
err_upper_bnd=0.999

if projection=='yz':
    rotation_matrix = np.array([[0,1,0],[0,0,1],[1,0,0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(rotation_matrix, data.gas.magnetic_flux_densities[i])
        data.gas.magnetic_vector_potentials[i] =np.matmul(rotation_matrix, data.gas.magnetic_vector_potentials[i])
if projection=='xz':
    rotation_matrix = np.array([[1,0,0],[0,0,1],[0,1,0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(rotation_matrix, data.gas.magnetic_flux_densities[i])
        data.gas.magnetic_vector_potentials[i] =np.matmul(rotation_matrix, data.gas.magnetic_vector_potentials[i])

def abs_vec(vec):
	res = np.sqrt(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2)
	return res
def dot_vec(vec1,vec2):
	res = vec1[:,0] * vec2[:,0] + vec1[:,1] * vec2[:,1] + vec1[:,2] * vec2[:,2]
	return res
def cross_vec(vec1,vec2):
	res_vec = np.zeros((len(vec1),3))
	res_vec[:,0] = vec1[:,1]*vec2[:,2]-vec1[:,2]*vec2[:,1]
	res_vec[:,1] = vec1[:,2]*vec2[:,0]-vec1[:,0]*vec2[:,2]
	res_vec[:,2] = vec1[:,0]*vec2[:,1]-vec1[:,1]*vec2[:,0]
	return res_vec
def rms_vec(vec):
	res = np.sqrt(np.mean(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2))
	return res
# see available fields in snapshot
#print(data.metadata.gas_properties.field_names)

# Get physical quantities
x = data.gas.coordinates[:, 0].value
y = data.gas.coordinates[:, 1].value
z = data.gas.coordinates[:, 2].value
rho = data.gas.densities.value
h = data.gas.smoothing_lengths.value
v = data.gas.velocities.value
P = data.gas.pressures.value
B = data.gas.magnetic_flux_densities.value

R0 = data.gas.r0.value
R1 = data.gas.r1.value
R2 = data.gas.r2.value
R3 = data.gas.r3.value

#Get RMS values
Brms = rms_vec(B)
if to_plot=='A':
    A = data.gas.magnetic_vector_potentials.value
    Arms = rms_vec(A)

#Print Brms
#print(f'Brms = {Brms}')

def make_slice(key, slice_frac_z=float(slice_height)):
        res = slice_gas(
        data,
        z_slice = slice_frac_z *data.metadata.boxsize[2],
        resolution=img_res,
        project=key,
        parallel=True,
        periodic=True,
        )
        return res

masses = make_slice("masses").value
dimy = len(masses)
dimx = len(masses[0])
mass_map = masses.flatten()
mass_map = mass_map + divreg

l = len(mass_map)

def prepare_sliced_quantity(quantity, isvec=False): 
	if isvec:
		data.gas.mass_weighted_temp_qx = data.gas.masses * quantity[:,0]
		data.gas.mass_weighted_temp_qy = data.gas.masses * quantity[:,1]
		data.gas.mass_weighted_temp_qz = data.gas.masses * quantity[:,2]
		sliced_quantity = np.zeros((l, 3))
		sliced_quantity[:,0] = make_slice("mass_weighted_temp_qx").value.flatten() / mass_map
		sliced_quantity[:,1] = make_slice("mass_weighted_temp_qy").value.flatten() / mass_map
		sliced_quantity[:,2] = make_slice("mass_weighted_temp_qz").value.flatten() / mass_map
	else:
		data.gas.mass_weighted_temp_q = data.gas.masses * quantity[:]
		sliced_quantity = make_slice("mass_weighted_temp_q").value.flatten() / mass_map
		
	return sliced_quantity


#slicing scalars
rho = prepare_sliced_quantity(rho)
h = prepare_sliced_quantity(h)
P = prepare_sliced_quantity(P)

#slicing vectors
v = prepare_sliced_quantity(v, isvec=True)
B = prepare_sliced_quantity(B, isvec=True)
if to_plot=='A':
	A = prepare_sliced_quantity(A, isvec=True)

#Get sliced error metrics
R0 = prepare_sliced_quantity(R0)
R1 = prepare_sliced_quantity(R1)
R2 = prepare_sliced_quantity(R2)
R3 = prepare_sliced_quantity(R3)

#mask error metrics
R0[R0<reg_err]=reg_err
R1[R1<reg_err]=reg_err
R2[R2<reg_err]=reg_err
R3[R3<reg_err]=reg_err

# plot everything

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import ticker

new_x = np.linspace(0.0, data.metadata.boxsize.value[0], dimx) #np.linspace(np.min(x), np.max(x), dimx)
new_y = np.linspace(0.0,  data.metadata.boxsize.value[1], dimy) #np.linspace(np.min(y), np.max(y), dimy)

def make_color_levels(cmin, cmax, c_res=10, log_sc=True):
    cmax += divreg
    if log_sc:
        levmin = int(np.floor(np.log10(cmin)))
        levmax = int(np.ceil(np.log10(cmax)))
        levels = []
        levels_short = []
        for i in range(levmin, levmax):
            levels_short += [10 ** i]
            for j in range(int(c_res / 10), c_res):
                levels += [(10 / c_res * j) * 10 ** i]

    else:
        levels = [cmin + (cmax - cmin) / c_res * k for k in range(c_res)]
        levels_short = [cmin + (cmax - cmin) / c_res * k for k in range(0, c_res, 10)]
    return levels, levels_short

def make_density_plot(
    Q, cmin, cmax, i, j, Q_name, c_res=10, log_sc=True, cmap="viridis"
):
    levels, levels_short = make_color_levels(cmin, cmax, c_res, log_sc)
    if log_sc:
        to_plot = ax[j].contourf(
            new_x,
            new_y,
            Q.transpose(),
            levels=np.array(levels),
            locator=ticker.LogLocator(),
            cmap=cmap,
            extend="both",
        )
    else:
        to_plot = ax[j].contourf(
            new_x, new_y, Q.transpose(), levels=np.array(levels), cmap=cmap,
            extend="both",
        )

    fig.colorbar(to_plot, ticks=levels_short)
    ax[j].set_ylabel(Q_name)
    ax[j].set_xlim(min(new_x),max(new_x))
    ax[j].set_ylim(min(new_y),max(new_y))
    return 0

#Plot magnetic fields
if to_plot=='A':
    Ax = A[:,0]/Arms
    Ay = A[:,1]/Arms
    Az = A[:,2]/Arms

    fig, ax = plt.subplots(1,3, sharex=True, figsize=(6*3, 5))
    make_density_plot(
            Ax.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            0,
            '$A_x$/$A_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    make_density_plot(
            Ay.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            1,
            '$A_y$/$A_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    make_density_plot(
            Az.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            2,
            '$A_z$/$A_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    ax[0].streamplot(new_x, new_y, np.transpose(Ax.reshape((dimx, dimy))), np.transpose(Ay.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)
    ax[1].streamplot(new_x, new_y, np.transpose(Ax.reshape((dimx, dimy))), np.transpose(Ay.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)
    ax[2].streamplot(new_x, new_y, np.transpose(Ax.reshape((dimx, dimy))), np.transpose(Ay.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)

if to_plot=='B':
    Bx = B[:,0]/Brms
    By = B[:,1]/Brms
    Bz = B[:,2]/Brms 

    fig, ax = plt.subplots(1,3, sharex=True, figsize=(6*3, 5))
    make_density_plot(
            Bx.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            0,
            '$B_x$/$B_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    make_density_plot(
            By.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            1,
            '$B_y$/$B_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    make_density_plot(
            Bz.reshape((dimx, dimy)),
            -1.0,
            1.0,
            0,
            2,
            '$B_z$/$B_{rms}$',
            c_res=cpts,
            log_sc=False,
            cmap=prefered_color
        )
    ax[0].streamplot(new_x, new_y, np.transpose(Bx.reshape((dimx, dimy))), np.transpose(By.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)
    ax[1].streamplot(new_x, new_y, np.transpose(Bx.reshape((dimx, dimy))), np.transpose(By.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)
    ax[2].streamplot(new_x, new_y, np.transpose(Bx.reshape((dimx, dimy))), np.transpose(By.reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)

if to_plot=='errors':

    fig, ax = plt.subplots(1,4, sharex=True, figsize=(24, 5)) #for 3 plts 18 for 4 plts use 24

    make_density_plot(
	R0.reshape((dimx, dimy)),
	reg_err,
	1.0,
	0,
	0,
	"$R_0$",
	c_res=cpts,
	#log_sc=False,
	cmap=prefered_color
	)
    make_density_plot(
        R1.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        1,
        "$R_1$",
        c_res=cpts,
        #log_sc=False,
        cmap=prefered_color
        )
    make_density_plot(
        R2.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        2,
        "$R_2$",
        c_res=cpts,
        #log_sc=False,
        cmap=prefered_color
        )
    make_density_plot(
        R3.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        3,
        "$R_3$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color
        )

ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
ax[1].set_title(f"t={t:.2e}")
ax[2].set_title("$z_{slice}/L_{box}$="+slice_height)
fig.tight_layout()

plt.savefig(sys.argv[2], dpi=100)
