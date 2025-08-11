from glob import glob
from swiftsimio import load
import unyt as u
import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rc
# Set the global font to be DejaVu Sans, size 10 (or any other sans-seri     f font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':10})
# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})



# setting units here 
units_dict = {'length':1.49597871 * 10**13 *u.cm,'density':u.g/u.cm**3,'velocity': u.km/u.s,'momentum': 1.98e33*u.g*u.km/u.s, 'pressure':u.g/(u.cm**1*u.s**2), 'magneticfield':1e-7*u.g / (u.s**2 * u.A),'magneticfieldpertime':1e-7*u.g/(u.s**2 * u.A * 3.156e7 * u.s),'magneticfieldgradient':1e-7*u.g/(u.s**2 * u.A * u.cm), 'time':u.s, 'mass':u.g,'acceleration':u.cm/u.s**2}
#parsec - 3.086e18*u.cm


# free fall time estimate
# Constants
G = 6.67430e-11 * u.m **3 / u.kg / u.s **2 #6.67430e-8
G = G.to(u.cm **3 / u.g / u.s **2)

def abs_vec(vec):
    res = np.sqrt(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2)
    return res

def dot_vec(vec1, vec2):
    res = vec1[:, 0] * vec2[:, 0] + vec1[:, 1] * vec2[:, 1] + vec1[:, 2]* vec2[:, 2]
    return res

def get_projection(vec,direction_vec):
    res = dot_vec(vec,direction_vec)/abs_vec(direction_vec)
    return res

def cross_vec(vec1, vec2):
    res_vec = np.zeros((len(vec1), 3))
    res_vec[:, 0] = vec1[:, 1] * vec2[:, 2] - vec1[:, 2] * vec2[:, 1]
    res_vec[:, 1] = vec1[:, 2] * vec2[:, 0] - vec1[:, 0] * vec2[:, 2]
    res_vec[:, 2] = vec1[:, 0] * vec2[:, 1] - vec1[:, 1] * vec2[:, 0]
    return res_vec


def extract_variables_of_interest(data):

    # Get snapshot parameters
    time = data.metadata.time
    z = data.metadata.z
    a = data.metadata.a

    pids = data.gas.particle_ids

    # Get physical quantities
    coordinates = data.gas.coordinates.to_physical().astype('float64')
    densities = data.gas.densities.to_physical().astype('float64')
    smoothing_lengths = data.gas.smoothing_lengths.to_physical().astype('float64')
    velocities = data.gas.velocities.to_physical().astype('float64')
    pressures = data.gas.pressures.to_physical().astype('float64')
    masses = data.gas.masses.to_physical().astype('float64')

    momentums = velocities * masses[:, np.newaxis]

    magnetic_flux_densities = data.gas.magnetic_flux_densities.to_physical().astype('float64')
    magnetic_flux_densitiesdt = data.gas.magnetic_flux_densitiesdt.to_physical().astype('float64')
    magnetic_divergences = data.gas.magnetic_divergences.to_physical().astype('float64')

#    r0 = data.gas.r0.to_physical().astype('float64')
#    r1 = data.gas.r1.to_physical().astype('float64')
#    r2 = data.gas.r2.to_physical().astype('float64')
#    r3 = data.gas.r3.to_physical().astype('float64')

    # force tracking variables

    total_force = data.gas.total_force.to_physical().astype('float64')
    
    lorentz_isotropic_force = data.gas.lorentz_isotropic_force.to_physical().astype('float64')
    lorentz_anisotropic_force = data.gas.lorentz_anisotropic_force.to_physical().astype('float64')
    monopole_correction_force = data.gas.monopole_correction_force.to_physical().astype('float64')


    # force magnitude in the velocity direction
    total_force_v = get_projection(total_force.value,velocities.value)*total_force.units
    lorentz_isotropic_force_v = get_projection(lorentz_isotropic_force.value,velocities.value)*lorentz_isotropic_force.units
    lorentz_anisotropic_force_v = get_projection(lorentz_anisotropic_force.value,velocities.value)*lorentz_anisotropic_force.units
    total_lorentz_force_v = lorentz_isotropic_force_v+lorentz_anisotropic_force_v
    monopole_correction_force_v = get_projection(monopole_correction_force.value,velocities.value)*monopole_correction_force.units
    total_hydro_force_v = total_force_v - total_lorentz_force_v - monopole_correction_force_v
    # B field source tracking variables
  
    stretching_bsource = data.gas.stretching_bsource.to_physical().astype('float64')

    dedner_bsource = data.gas.dedner_bsource.to_physical().astype('float64')

    artificial_resistivity_bsource = data.gas.artificial_resistivity_bsource.to_physical().astype('float64')

    physical_resistivity_bsource = data.gas.physical_resistivity_bsource.to_physical().astype('float64')

    # force magnitude in the velocity direction
    stretching_bsource_B = get_projection(stretching_bsource.value,magnetic_flux_densities.value)*stretching_bsource.units
    dedner_bsource_B = get_projection(dedner_bsource.value,magnetic_flux_densities.value)*dedner_bsource.units
    artificial_resistivity_bsource_B = get_projection(artificial_resistivity_bsource.value,magnetic_flux_densities.value)*artificial_resistivity_bsource.units
    physical_resistivity_bsource_B = get_projection(physical_resistivity_bsource.value,magnetic_flux_densities.value)*physical_resistivity_bsource.units
 

    # Center of mass position tracking
    CM_frame_coordinates = coordinates
    for k in range(3):
        CM_k = np.sum(coordinates[:,k].value * masses.value)/np.sum(masses.value)*coordinates.units
        CM_frame_coordinates -= CM_k

    #for k in range(3):
    #    lorentz_isotropic_force[:,k] /= total_force_abs
    #    lorentz_anisotropic_force[:,k] /= total_force_abs
    #    monopole_correction_force[:,k] /= total_force_abs
    # working with units for snapshot parameters:
    time.to(units_dict['time'])

    # working with units:
    coordinates = coordinates.to(units_dict['length'])
    CM_frame_coordinates = CM_frame_coordinates.to(units_dict['length'])
    densities = densities.to(units_dict['density'])
    smoothing_lengths = smoothing_lengths.to(units_dict['length'])
    velocities = velocities.to(units_dict['velocity'])
    momentums = momentums.to(units_dict['momentum'])
    pressures = pressures.to(units_dict['pressure'])
    magnetic_flux_densitiesdt = magnetic_flux_densitiesdt.to(units_dict['magneticfieldpertime'])
    magnetic_divergences = magnetic_divergences.to(units_dict['magneticfieldgradient']) 
    masses = masses.to(units_dict['mass'])
 
    total_force_v = total_force_v.to(units_dict['acceleration'])
    lorentz_isotropic_force_v = lorentz_isotropic_force_v.to(units_dict['acceleration'])
    lorentz_anisotropic_force_v = lorentz_anisotropic_force_v.to(units_dict['acceleration'])
    total_lorentz_force_v = total_lorentz_force_v.to(units_dict['acceleration']) 
    monopole_correction_force_v = monopole_correction_force_v.to(units_dict['acceleration'])
    total_hydro_force_v = total_hydro_force_v.to(units_dict['acceleration'])

    r0_no_cut = magnetic_divergences * smoothing_lengths / (abs_vec(magnetic_flux_densities.value)*magnetic_flux_densities.units)

    mu0 = 1.25663706127e-1 * u.g * u.cm / (u.s**2 * u.A**2)
    plasma_beta = pressures / ((abs_vec(magnetic_flux_densities.value)*magnetic_flux_densities.units)**2/(2*mu0))

    stretching_bsource = stretching_bsource.to(units_dict['magneticfieldpertime'])
    dedner_bsource = data.gas.dedner_bsource.to(units_dict['magneticfieldpertime'])
    artificial_resistivity_bsource = data.gas.artificial_resistivity_bsource.to(units_dict['magneticfieldpertime'])
    physical_resistivity_bsource = data.gas.physical_resistivity_bsource.to(units_dict['magneticfieldpertime'])

    CM_frame_z = np.abs(CM_frame_coordinates[:,2])
    data_dict = {'pids':pids,'coordinates':coordinates,'densities':densities,'smoothing_lengths':smoothing_lengths,'velocities':velocities,'momentums':momentums,'pressures':pressures,'magnetic_flux_densities':magnetic_flux_densities,'CM_frame_coordinates':CM_frame_coordinates,'CM_frame_z':CM_frame_z,'plasma_beta':plasma_beta}
    
    data_dict = {'pids':pids,'coordinates':coordinates,'densities':densities,'smoothing_lengths':smoothing_lengths,'velocities':velocities,'pressures':pressures,'magnetic_flux_densities':magnetic_flux_densities,'magnetic_flux_densitiesdt':magnetic_flux_densitiesdt,'magnetic_divergences':magnetic_divergences,'r0_no_cut':r0_no_cut,'lorentz_isotropic_force_v':lorentz_isotropic_force_v,'lorentz_anisotropic_force_v':lorentz_anisotropic_force_v,'total_lorentz_force_v':total_lorentz_force_v, 'monopole_correction_force_v':monopole_correction_force_v, 'total_force_v':total_force_v,'total_hydro_force_v':total_hydro_force_v,'stretching_bsource_B':stretching_bsource_B,'dedner_bsource_B':dedner_bsource_B,'artificial_resistivity_bsource_B':artificial_resistivity_bsource_B,'physical_resistivity_bsource_B':physical_resistivity_bsource_B,'CM_frame_coordinates':CM_frame_coordinates,'CM_frame_z':CM_frame_z,'plasma_beta':plasma_beta}
    
    values = {}
    units = {}
    for key in data_dict.keys():
        values = values | {key:data_dict[key].value}

    snapshot_parameters = {'time':time.value, 'z':z, 'a':a}

    return {'values':values,'parameters':snapshot_parameters}

def updoad_shapshots_data_from_folder(folderaddr):
    addr_book = glob(folderaddr+'/**/*.hdf5',recursive=True)
    addr_book = sorted(addr_book)
    names = [addr.split('/')[-1].split('.hdf5')[0] for addr in addr_book]
    #print(names)
    snapshots_data = []
    for addr in addr_book:
        snapshot = load(addr)
        snapshots_data+=[extract_variables_of_interest(snapshot)]
    return snapshots_data, names

def upload_particle_history(all_data,pid):
    particle_data = []
    for snapshot in all_data:
        particle_snapshot_dict = snapshot['parameters']
        indx = np.argwhere(snapshot['values']['pids']==pid)[0][0]
        for key in snapshot['values'].keys():
            particle_snapshot_dict = particle_snapshot_dict | {key:snapshot['values'][key][indx]}
        particle_data+=[particle_snapshot_dict]

    # transform list of dicts to single 
    # Dynamically extract keys from the first dictionary
    keys = particle_data[0].keys()

    # Transform the list of dictionaries to a single dictionary
    result = {key: np.array([entry[key] for entry in particle_data]) for key in keys}
    result = result | {'pid':pid}
    return result

def identify_largest_quantity_particles(all_snapshots,quantity,isvec=False,region_cut=True):
    largest_quantity_particles = []
    for snap in all_snapshots:
        quantity_values = snap['values'][quantity]
        zcoords = np.abs(snap['values']['CM_frame_coordinates'][:,2])
        rcoords = np.sqrt(snap['values']['CM_frame_coordinates'][:,0]**2+snap['values']['CM_frame_coordinates'][:,1]**2)
        velocities = snap['values']['velocities']
        vabs = np.sqrt(velocities[:,0]**2+velocities[:,1]**2+velocities[:,2]**2) 
        mask = np.array([True for i in range(len(quantity_values))])
        if region_cut:
            zcut = (zcoords > 0.015*0.015) & (zcoords < 0.15*0.015) # z>0.015 Rcloud
            rcut = (rcoords < 0.15*0.015)
  
            mask = (zcut) & (rcut)

        #quantity_values = quantity_values[zcut]
        if isvec:
            quantity_values = np.sqrt(quantity_values[:,0]**2+quantity_values[:,1]**2+quantity_values[:,2]**2)
        pids = snap['values']['pids']
        quantity_values_cut = quantity_values[mask]
        #print('snapshot at t=',snap['parameters']['time'])
        indx = np.argmax(quantity_values_cut)
        #print('particle #',pids[mask][indx])
        #print('quantity value = ',quantity_values[pids[mask][indx]])
        largest_quantity_particles.append(pids[mask][indx])
    return largest_quantity_particles

def plot_quatities_for_particle_vs_time(particle_history,quantity_list,outputfileaddr='./quantities.png',time_var = 'time',timeinterval=[1.0,1.5],title='quantities vs time',customtimeinterval = False,customquantityinterval = False, quantityinterval=[1e-1,1e1],logscale=True):

    fig, ax = plt.subplots(1, len(quantity_list), figsize=(6.5*len(quantity_list), 6))
    for i in range(len(quantity_list)):
        quantity = particle_history[quantity_list[i]]
        times = particle_history[time_var]
        if len(quantity.shape)>1:
            quantity = np.sqrt(quantity[:,0]**2+quantity[:,1]**2+quantity[:,2]**2)

        mask_nonzero = quantity!=0.0
        #quantity=quantity[mask_nonzero] 
        #times = times[mask_nonzero]
        mask_in_range = (times>timeinterval[0]) & (times<timeinterval[-1])
        #print(len(times[mask_in_range]))
        #if len(quantity[mask_in_range])!=0:
            #ax[i].axvline(x=1.378579,ymin = 0.1*np.min(quantity[mask_in_range]),ymax = 10*np.max(quantity[mask_in_range]),color='red',linestyle='dashed')
            #ax[i].axvline(x=1.378579,ymin = 1e-30,ymax = 1e30,color='red',linestyle='dashed')

        #handle only positive and positive-negative values:
        if logscale:
            mask = quantity<0
            if np.sum(mask)>0:
                # Separate positive and negative values
                positive = np.where(quantity > 0, quantity, np.nan)  # Use NaN for negative values
                negative = np.where(quantity < 0, quantity, np.nan)  # Use NaN for positive values
                ax[i].plot(times, positive,color='red',marker='.')
                ax[i].plot(times,-negative,color='blue',marker='.')
            else:
                ax[i].plot(times,quantity,color='black',marker='.')
        else:
            ax[i].plot(times,quantity,color='black',marker='.')
        ax[i].set_ylabel(quantity_list[i],fontsize=20)
        ax[i].set_xlabel('$t/t_{ff}$',fontsize=20)
        if logscale:
            ax[i].set_yscale('log')
        ax[i].tick_params(axis='both', labelsize=20)
        if customtimeinterval:
            ax[i].set_xlim(timeinterval)
        if customquantityinterval:
            ax[i].set_ylim(quantityinterval)
        ax[i].grid()
    fig.suptitle('Particle #'+str(particle_history['pid']),fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)


    return 0

Lcut = 200 # Lcut in A.U.
def plot_quantity_vs_time(all_snapshots,quantity,outputfileaddr,isvec=False,region_cut=True,label=''):
    quantity_array = []
    rmsq_array = []
    mzq_array = []
    times_array = []
    npart_array = []
    for snap in all_snapshots:
        quantity_values = snap['values'][quantity]
        #total_quantity = np.sum(quantity_values,axis=0)
        #rms_quantity = np.sqrt(np.mean(quantity_values[:,0]**2+quantity_values[:,1]**2+quantity_values[:,2]**2))
        print(len(quantity_values))   
        if region_cut:
            z = snap['values']['CM_frame_coordinates'][:,2]
            y = snap['values']['CM_frame_coordinates'][:,1]
            x = snap['values']['CM_frame_coordinates'][:,0]
            mask = ((np.abs(z)<=Lcut) & (np.abs(y)<=Lcut) & (np.abs(x)<=Lcut))
            quantity_values = quantity_values[mask]

        total_quantity = np.sum(quantity_values,axis=0)
        rms_quantity = np.sqrt(np.mean(quantity_values[:,0]**2+quantity_values[:,1]**2+quantity_values[:,2]**2))
        mz_quantity = np.mean(np.abs(quantity_values[:,2]))
        time = snap['parameters']['time']
        quantity_array.append(total_quantity)
        rmsq_array.append(rms_quantity)
        mzq_array.append(mz_quantity)
        times_array.append(time)    
        npart_array.append(len(quantity_values))
    
    quantity_array = np.array(quantity_array)
    rmsq_array = np.array(rmsq_array)
    mzq_array = np.array(mzq_array)
    times_array = np.array(times_array)
    npart_array = np.array(npart_array)


    fig, ax = plt.subplots(1, 2, figsize=(6.5*2, 6))
    if isvec:
        quantity_abs = np.sqrt(quantity_array[:,0]**2+quantity_array[:,1]**2+quantity_array[:,2]**2)
    else:
        quantity_abs = np.abs(quantity_array)

    #ax[0].plot(times_array,quantity_array[:,2],color='black',marker='.', label='$p_z$')
    #ax[0].plot(times_array,quantity_array[:,0],color='blue' ,marker='.', label='$p_x$')
    #ax[0].plot(times_array,quantity_array[:,1],color='red' ,marker='.', label='$p_y$')
    #ax[0].plot(times_array,np.abs(quantity_array[:,2]),color='black',marker='.', label='$p_z$')
    pq = quantity_array[:,2]#/npart_array
    mask = pq<0
    lpq = '$mean p_z$'
    if np.sum(mask)>0:
        # Separate positive and negative values
        positive = np.where(pq > 0, pq, np.nan)  # Use NaN for negative values
        negative = np.where(pq < 0, pq, np.nan)  # Use NaN for positive values
        ax[0].plot(times_array, positive,color='red', marker='.',label = 'mean $p_z>0$')
        ax[0].plot(times_array,-negative,color='blue',marker='.',label = 'mean $p_z<0$')
    else:
        ax[0].plot(times_array,pq,color='black',marker='.',label = 'mean $p_z>0$')

    #ax[0].plot(times_array, rmsq_array, color='green', marker='.',label='$p_{rms}$')
    #ax[0].plot(times_array, mzq_array, color='black', marker='.',label='mean $|p_z|$')

    #ax[0].plot(times_array,np.abs(quantity_array[:,0]),color='blue' ,marker='.', label='$p_x$')
    #ax[0].plot(times_array,np.abs(quantity_array[:,1]),color='red'  ,marker='.', label='$p_y$')
    
    ax[0].legend()
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].set_ylabel('$p_z^{tot}$, [$M_{sol}\cdot km / s$]',fontsize=20)
    ax[0].set_xlabel('$t/t_{ff}$',fontsize=20)    
    ax[0].set_ylim([1e-5,1e0])
    #ax[0].set_ylim([1e-8,1e-4])

    ax[1].plot(times_array,npart_array,color='black',marker='.')
    ax[1].grid()
    ax[1].set_yscale('log')
    ax[1].set_ylabel(f'$N_p$'+label ,fontsize=20) #, Lcut = {Lcut} A.U.',fontsize=20)
    ax[1].set_xlabel('$t/t_{ff}$',fontsize=20)    
    ax[1].set_ylim([1e2,1e5])

    #fig.suptitle(quantity,fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)
    return 

# 34013   

folder = './MCC_ODI2_control'

snapshots_history, snapshot_names = updoad_shapshots_data_from_folder(folder)

#plot_quantity_vs_time(snapshots_history,'momentums',outputfileaddr=folder+'/momentums_noRegCut.png',isvec=True,region_cut=False, label = ' All')
#plot_quantity_vs_time(snapshots_history,'momentums',outputfileaddr=folder+'/momentums_cut.png',isvec=True,region_cut=True, label = f' $Lbox=${2*Lcut}')

# plot forces for a given particle
pid = 10000
particle_history = upload_particle_history(snapshots_history,pid)
B_vs_time = plot_quatities_for_particle_vs_time(particle_history,['CM_frame_z','densities','velocities','magnetic_flux_densities'],outputfileaddr='./tracked_particles/output_rho_v_B_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])
#errors_vs_time = plot_quatities_for_particle_vs_time(particle_history,['r0','r1','r2'],outputfileaddr='./tracked_particles/output_r0_r1_r2_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])
divergences_and_dedner_vs_time = plot_quatities_for_particle_vs_time(particle_history,['magnetic_flux_densitiesdt','magnetic_divergences'],outputfileaddr='./tracked_particles/output_divB_dBdt_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])
force_vs_time = plot_quatities_for_particle_vs_time(particle_history,['total_force','lorentz_isotropic_force','lorentz_anisotropic_force','monopole_correction_force'],outputfileaddr='./tracked_particles/forces_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e0,1e13])
bsource_vs_time = plot_quatities_for_particle_vs_time(particle_history,['stretching_bsource','dedner_bsource'],outputfileaddr='./tracked_particles/bsources_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e-17,1e7])


#largest_qty_particles = identify_largest_quantity_particles(snapshots_history,quantity = 'magnetic_flux_densities', isvec = True,region_cut=False )

#largest_qty_particles = identify_largest_quantity_particles(snapshots_history,quantity = 'monopole_correction_force', isvec = True,region_cut=False )

##largest_qty_particles = identify_largest_quantity_particles(snapshots_history,quantity = 'velocities', isvec = True,region_cut=False)


#largest_qty_particles = identify_largest_quantity_particles(snapshots_history,quantity = 'r0_no_cut', isvec = False,region_cut=False)

##unique_pids = np.unique(largest_qty_particles[:])

##for pid in unique_pids:
  ##  print(f'making history plots for #{pid}')
  ##  particle_history = upload_particle_history(snapshots_history,pid)

    ##B_vs_time = plot_quatities_for_particle_vs_time(particle_history,['CM_frame_z','densities','velocities','magnetic_flux_densities'],outputfileaddr='./tracked_particles/output_rho_v_B_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

    ##errors_vs_time = plot_quatities_for_particle_vs_time(particle_history,['r0','r1','r2','r0_no_cut'],outputfileaddr='./tracked_particles/output_r0_r1_r2_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

    ##divergences_and_dedner_vs_time = plot_quatities_for_particle_vs_time(particle_history,['magnetic_flux_densitiesdt','magnetic_divergences'],outputfileaddr='./tracked_particles/output_divB_dBdt_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

    ##force_vs_time = plot_quatities_for_particle_vs_time(particle_history,['total_hydro_force_v','total_lorentz_force_v','monopole_correction_force_v'],outputfileaddr='./tracked_particles/forces_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e-7,2e-3], logscale=True)

    ##monopole_force = plot_quatities_for_particle_vs_time(particle_history,['monopole_correction_force_v','plasma_beta'],outputfileaddr='./tracked_particles/monopole_and_beta_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = False, quantityinterval=[1e-7,2e-3], logscale=True)

    ##bsource_vs_time = plot_quatities_for_particle_vs_time(particle_history,['stretching_bsource_B','dedner_bsource_B'],outputfileaddr='./tracked_particles/bsources_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e-7,1e1],logscale=True)
#

#snap_nr = 1380-1200
#pid = largest_v_particles[snap_nr]
#print(snapshot_names[snap_nr])
#particle_history = upload_particle_history(snapshots_history,pid)



#largest_magnetic_force = identify_largest_quantity_particles(snapshots_history,quantity = 'monopole_correction_force', isvec = True )

#largest_F_pid = largest_magnetic_force[-11]
#particle_history = upload_particle_history(snapshots_history,largest_F_pid)

#B_vs_time = plot_quatities_for_particle_vs_time(particle_history,['CM_frame_z','densities','velocities','magnetic_flux_densities'],outputfileaddr='./tracked_particles/output_rho_v_B_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

#errors_vs_time = plot_quatities_for_particle_vs_time(particle_history,['r0','r1','r2'],outputfileaddr='./tracked_particles/output_r0_r1_r2_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

#divergences_and_dedner_vs_time = plot_quatities_for_particle_vs_time(particle_history,['magnetic_flux_densitiesdt','magnetic_divergences'],outputfileaddr='./tracked_particles/output_divB_dBdt_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2])

#force_vs_time = plot_quatities_for_particle_vs_time(particle_history,['total_force','lorentz_isotropic_force','lorentz_anisotropic_force','monopole_correction_force'],outputfileaddr='./tracked_particles/forces_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e0,1e13])

#bsource_vs_time = plot_quatities_for_particle_vs_time(particle_history,['stretching_bsource','dedner_bsource'],outputfileaddr='./tracked_particles/bsources_'+str(pid)+'.png',time_var = 'time',customtimeinterval = False,timeinterval = [0.0,0.2],customquantityinterval = True, quantityinterval=[1e-17,1e7])
#print(B_vs_time)
