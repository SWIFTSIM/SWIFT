from glob import glob
from swiftsimio import load
import unyt as u
import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rc

# Set the global font to be DejaVu Sans, size 10 (or any other sans-seri     f font of your choice!)
rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
# Set the font used for MathJax - more on this later
rc("mathtext", **{"default": "regular"})


# setting units here
units_dict = {
    "length": 3.08e21 * u.cm,
    "density": 1.98e33 * u.g / (3.08e21 * u.cm) ** 3,
    "velocity": u.km / u.s,
    "momentum": 1.98e33 * u.g * u.km / u.s,
    "pressure": 1e10 * (1.98e33 * u.g) / ((3.08e21 * u.cm) * (3.15e16 * u.s) ** 2),
    "magneticfield": 1e-7 * u.g / (u.s ** 2 * u.statA),
    "magneticfieldpertime": 1e-7 * u.g / (u.s ** 2 * u.statA * 3.156e7 * u.s),
    "magneticfieldgradient": 1e-7 * u.g / (u.s ** 2 * u.statA * u.cm),
    "time": 3.15e16 * u.s,
    "mass": 1.98e33 * u.g,
    "acceleration": u.cm / u.s ** 2,
}
# parsec - 3.086e18*u.cm


# free fall time estimate
# Constants
G = 6.67430e-11 * u.m ** 3 / u.kg / u.s ** 2  # 6.67430e-8
G = G.to(u.cm ** 3 / u.g / u.s ** 2)
# Parameters (taken from Hopkins 2016)
Rcloud = 4.628516371e16 * u.cm
M = 1.99e33 * u.g  # total mass of the sphere
rhocloud0 = M / (4 / 3 * np.pi * Rcloud ** 3)
rhouniform = rhocloud0 / 360
t_ff_lit = np.sqrt(3 / (2 * np.pi * G * rhocloud0))
t_ff_wiki = np.sqrt(3 * np.pi / (32 * G * rhocloud0))
t_ff_script = 3e4 * 3.156e7 * u.s


def abs_vec(vec):
    res = np.sqrt(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2)
    return res


def dot_vec(vec1, vec2):
    res = vec1[:, 0] * vec2[:, 0] + vec1[:, 1] * vec2[:, 1] + vec1[:, 2] * vec2[:, 2]
    return res


def get_projection(vec, direction_vec):
    res = dot_vec(vec, direction_vec) / abs_vec(direction_vec)
    return res


def cross_vec(vec1, vec2):
    res_vec = np.zeros((len(vec1), 3))
    res_vec[:, 0] = vec1[:, 1] * vec2[:, 2] - vec1[:, 2] * vec2[:, 1]
    res_vec[:, 1] = vec1[:, 2] * vec2[:, 0] - vec1[:, 0] * vec2[:, 2]
    res_vec[:, 2] = vec1[:, 0] * vec2[:, 1] - vec1[:, 1] * vec2[:, 0]
    return res_vec


def extract_variables_of_interest(data):

    # set constant
    mu0 = 1.25663706127e-1 * u.g * u.cm / (u.s ** 2 * u.statA ** 2)

    # Get snapshot parameters
    time = data.metadata.time
    z = data.metadata.z
    a = data.metadata.a

    pids = data.gas.particle_ids

    # Get physical quantities
    coordinates = data.gas.coordinates.to_physical().astype("float64")
    densities = data.gas.densities.to_physical().astype("float64")
    smoothing_lengths = data.gas.smoothing_lengths.to_physical().astype("float64")
    velocities = data.gas.velocities.to_physical().astype("float64")
    pressures = data.gas.pressures.to_physical().astype("float64")
    masses = data.gas.masses.to_physical().astype("float64")

    momentums = velocities * masses[:, np.newaxis]

    magnetic_flux_densities = data.gas.magnetic_flux_densities.to_physical().astype(
        "float64"
    )
    magnetic_flux_densitiesdt = data.gas.magnetic_flux_densitiesdt.to_physical().astype(
        "float64"
    )
    magnetic_divergences = data.gas.magnetic_divergences.to_physical().astype("float64")

    magnetic_pressures = dot_vec(magnetic_flux_densities, magnetic_flux_densities) / (
        2 * mu0
    )
    dynamic_pressures = densities * dot_vec(velocities, velocities)

    # Center of mass position tracking
    CM_frame_coordinates = coordinates
    for k in range(3):
        CM_k = (
            np.sum(coordinates[:, k].value * masses.value)
            / np.sum(masses.value)
            * coordinates.units
        )
        CM_frame_coordinates -= CM_k

    # working with units:
    coordinates = coordinates.to(units_dict["length"])
    CM_frame_coordinates = CM_frame_coordinates.to(units_dict["length"])
    densities = densities.to(units_dict["density"])
    smoothing_lengths = smoothing_lengths.to(units_dict["length"])
    velocities = velocities.to(units_dict["velocity"])
    momentums = momentums.to(units_dict["momentum"])
    pressures = pressures.to(units_dict["pressure"])
    magnetic_flux_densities = magnetic_flux_densities.to(
        units_dict["magneticfield"]
    )
    magnetic_flux_densitiesdt = magnetic_flux_densitiesdt.to(
        units_dict["magneticfieldpertime"]
    )
    magnetic_divergences = magnetic_divergences.to(units_dict["magneticfieldgradient"])
    masses = masses.to(units_dict["mass"])
    magnetic_pressures = magnetic_pressures.to(units_dict["pressure"])
    dynamic_pressures = dynamic_pressures.to(units_dict["pressure"])

    magnetic_flux_densities_norm = abs_vec(magnetic_flux_densities)

    time = time.to(units_dict["time"])

    r0_no_cut = (
        magnetic_divergences
        * smoothing_lengths
        / (abs_vec(magnetic_flux_densities.value) * magnetic_flux_densities.units)
    )

    CM_frame_z = np.abs(CM_frame_coordinates[:, 2])
    data_dict = {
        "pids": pids,
        "coordinates": coordinates,
        "masses": masses,
        "densities": densities,
        "smoothing_lengths": smoothing_lengths,
        "velocities": velocities,
        "momentums": momentums,
        "pressures": pressures,
        "magnetic_flux_densities": magnetic_flux_densities,
        "CM_frame_coordinates": CM_frame_coordinates,
        "CM_frame_z": CM_frame_z,
        "magnetic_divergences":magnetic_divergences,
        "magnetic_flux_densities_norm": magnetic_flux_densities_norm,
        "magnetic_pressures": magnetic_pressures,
        "dynamic_pressures": dynamic_pressures,
    }

    values = {}
    units = {}
    for key in data_dict.keys():
        values = values | {key: data_dict[key].value}

    snapshot_parameters = {"time": time.value, "z": z, "a": a}

    # Retrieve some information about the simulation run
    artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
    dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
    dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
    dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
    eta = data.metadata.hydro_scheme["Resistive Eta"]
    git = data.metadata.code["Git Revision"]
    gitBranch = data.metadata.code["Git Branch"]
    hydroScheme = data.metadata.hydro_scheme["Scheme"]
    kernel = data.metadata.hydro_scheme["Kernel function"]
    neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]

    metadata = {
        "artDiffusion": artDiffusion,
        "dedHyp": dedHyp,
        "dedHypDivv": dedHypDivv,
        "dedPar": dedPar,
        "eta": eta,
        "git": git,
        "gitBranch": gitBranch,
        "hydroScheme": hydroScheme,
        "kernel": kernel,
        "neighbours": neighbours,
    }

    return {"values": values, "parameters": snapshot_parameters, "metadata": metadata}


def updoad_shapshots_data_from_folder(folderaddr):
    addr_book = glob(folderaddr + "/**/*.hdf5", recursive=True)
    addr_book = sorted(addr_book)
    names = [addr.split("/")[-1].split(".hdf5")[0] for addr in addr_book]
    # print(names)
    snapshots_data = []
    for addr in addr_book:
        snapshot = load(addr)
        snapshots_data += [extract_variables_of_interest(snapshot)]
    return snapshots_data, names


def upload_particle_history(all_data, pid):
    particle_data = []
    for snapshot in all_data:
        particle_snapshot_dict = snapshot["parameters"]
        indx = np.argwhere(snapshot["values"]["pids"] == pid)[0][0]
        for key in snapshot["values"].keys():
            particle_snapshot_dict = particle_snapshot_dict | {
                key: snapshot["values"][key][indx]
            }
        particle_data += [particle_snapshot_dict]

    # transform list of dicts to single
    # Dynamically extract keys from the first dictionary
    keys = particle_data[0].keys()

    # Transform the list of dictionaries to a single dictionary
    result = {key: np.array([entry[key] for entry in particle_data]) for key in keys}
    result = result | {"pid": pid}
    return result


def identify_largest_quantity_particles(
    all_snapshots, quantity, isvec=False, region_cut=True
):
    largest_quantity_particles = []
    for snap in all_snapshots:
        quantity_values = snap["values"][quantity]
        zcoords = np.abs(snap["values"]["CM_frame_coordinates"][:, 2])
        rcoords = np.sqrt(
            snap["values"]["CM_frame_coordinates"][:, 0] ** 2
            + snap["values"]["CM_frame_coordinates"][:, 1] ** 2
        )
        velocities = snap["values"]["velocities"]
        vabs = np.sqrt(
            velocities[:, 0] ** 2 + velocities[:, 1] ** 2 + velocities[:, 2] ** 2
        )
        mask = np.array([True for i in range(len(quantity_values))])
        if region_cut:
            zcut = (zcoords > 0.015 * 0.015) & (
                zcoords < 0.15 * 0.015
            )  # z>0.015 Rcloud
            rcut = rcoords < 0.15 * 0.015

            mask = (zcut) & (rcut)

        # quantity_values = quantity_values[zcut]
        if isvec:
            quantity_values = np.sqrt(
                quantity_values[:, 0] ** 2
                + quantity_values[:, 1] ** 2
                + quantity_values[:, 2] ** 2
            )
        pids = snap["values"]["pids"]
        quantity_values_cut = quantity_values[mask]
        # print('snapshot at t=',snap['parameters']['time'])
        indx = np.argmax(quantity_values_cut)
        # print('particle #',pids[mask][indx])
        # print('quantity value = ',quantity_values[pids[mask][indx]])
        largest_quantity_particles.append(pids[mask][indx])
    return largest_quantity_particles


def plot_quatities_for_particle_vs_time(
    particle_history,
    quantity_list,
    outputfileaddr="./quantities.png",
    time_var="time",
    timeinterval=[0, 3],
    title="quantities vs time",
    customtimeinterval=False,
    customquantityinterval=False,
    quantityinterval=[1e-1, 1e1],
    logscale=True,
):

    fig, ax = plt.subplots(1, len(quantity_list), figsize=(6.5 * len(quantity_list), 6))
    for i in range(len(quantity_list)):
        quantity = particle_history[quantity_list[i]]
        times = particle_history[time_var]
        if len(quantity.shape) > 1:
            quantity = np.sqrt(
                quantity[:, 0] ** 2 + quantity[:, 1] ** 2 + quantity[:, 2] ** 2
            )

        mask_nonzero = quantity != 0.0
        # quantity=quantity[mask_nonzero]
        # times = times[mask_nonzero]
        mask_in_range = (times > timeinterval[0]) & (times < timeinterval[-1])
        # print(len(times[mask_in_range]))
        # if len(quantity[mask_in_range])!=0:
        # ax[i].axvline(x=1.378579,ymin = 0.1*np.min(quantity[mask_in_range]),ymax = 10*np.max(quantity[mask_in_range]),color='red',linestyle='dashed')
        # ax[i].axvline(x=1.378579,ymin = 1e-30,ymax = 1e30,color='red',linestyle='dashed')

        # handle only positive and positive-negative values:
        if logscale:
            mask = quantity < 0
            if np.sum(mask) > 0:
                # Separate positive and negative values
                positive = np.where(
                    quantity > 0, quantity, np.nan
                )  # Use NaN for negative values
                negative = np.where(
                    quantity < 0, quantity, np.nan
                )  # Use NaN for positive values
                ax[i].plot(times, positive, color="red", marker=".")
                ax[i].plot(times, -negative, color="blue", marker=".")
            else:
                ax[i].plot(times, quantity, color="black", marker=".")
        else:
            ax[i].plot(times, quantity, color="black", marker=".")
        ax[i].set_ylabel(quantity_list[i], fontsize=20)
        ax[i].set_xlabel("$time$", fontsize=20)
        if logscale:
            ax[i].set_yscale("log")
        ax[i].tick_params(axis="both", labelsize=20)
        if customtimeinterval:
            ax[i].set_xlim(timeinterval)
        if customquantityinterval:
            ax[i].set_ylim(quantityinterval)
        ax[i].grid()
    fig.suptitle("Particle #" + str(particle_history["pid"]), fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)

    return 0


Lcut = 15  # Lcut in kPc.


def plot_quantity_vs_time(
    all_snapshots, quantity, outputfileaddr, isvec=False, region_cut=True, label=""
):
    quantity_array = []
    rmsq_array = []
    mzq_array = []
    times_array = []
    npart_array = []
    for snap in all_snapshots:
        quantity_values = snap["values"][quantity]
        Np = len(quantity_values)
        # total_quantity = np.sum(quantity_values,axis=0)
        # rms_quantity = np.sqrt(np.mean(quantity_values[:,0]**2+quantity_values[:,1]**2+quantity_values[:,2]**2))
        print(len(quantity_values))
        #if region_cut:
        #    z = snap["values"]["CM_frame_coordinates"][:, 2]
        #    y = snap["values"]["CM_frame_coordinates"][:, 1]
        #    x = snap["values"]["CM_frame_coordinates"][:, 0]
        #    mask = (np.abs(z) <= Lcut) & (np.abs(y) <= Lcut) & (np.abs(x) <= Lcut)
        rcut = 15
        zcut = 0.5
        if region_cut:
            z = snap["values"]["CM_frame_coordinates"][:, 2] * units_dict["length"]
            y = snap["values"]["CM_frame_coordinates"][:, 1] * units_dict["length"]
            x = snap["values"]["CM_frame_coordinates"][:, 0] * units_dict["length"]
            r = np.sqrt(x ** 2 + y ** 2)
            mask = (r <= rcut) & (np.abs(z) <= zcut)
         
            quantity_values = quantity_values[mask]

        total_quantity = np.sum(quantity_values, axis=0)
        rms_quantity = np.sqrt(
            np.mean(
                quantity_values[:, 0] ** 2
                + quantity_values[:, 1] ** 2
                + quantity_values[:, 2] ** 2
            )
        )
        mz_quantity = np.mean(np.abs(quantity_values[:, 2]))
        time = snap["parameters"]["time"]
        quantity_array.append(total_quantity)
        rmsq_array.append(rms_quantity)
        mzq_array.append(mz_quantity)
        times_array.append(time)
        npart_array.append(len(quantity_values))

        # Retrieve some information about the simulation run
        artDiffusion = snap["metadata"]["artDiffusion"]
        dedHyp = snap["metadata"]["dedHyp"]
        dedHypDivv = snap["metadata"]["dedHypDivv"]
        dedPar = snap["metadata"]["dedPar"]
        eta = snap["metadata"]["eta"]
        git = snap["metadata"]["git"]
        gitBranch = snap["metadata"]["gitBranch"]
        hydroScheme = snap["metadata"]["hydroScheme"]
        kernel = snap["metadata"]["kernel"]
        neighbours = snap["metadata"]["neighbours"]

    quantity_array = np.array(quantity_array)
    rmsq_array = np.array(rmsq_array)
    mzq_array = np.array(mzq_array)
    times_array = np.array(times_array)
    npart_array = np.array(npart_array)

    fig, ax = plt.subplots(1, 3, figsize=(5.5 * 3, 5))
    if isvec:
        quantity_abs = np.sqrt(
            quantity_array[:, 0] ** 2
            + quantity_array[:, 1] ** 2
            + quantity_array[:, 2] ** 2
        )
    else:
        quantity_abs = np.abs(quantity_array)

    # ax[0].plot(times_array,quantity_array[:,2],color='black',marker='.', label='$p_z$')
    # ax[0].plot(times_array,quantity_array[:,0],color='blue' ,marker='.', label='$p_x$')
    # ax[0].plot(times_array,quantity_array[:,1],color='red' ,marker='.', label='$p_y$')
    # ax[0].plot(times_array,np.abs(quantity_array[:,2]),color='black',marker='.', label='$p_z$')
    # pq = quantity_array[:,2]#/npart_array
    # mask = pq<0
    # lpq = '$mean p_z$'
    # if np.sum(mask)>0:
    #    # Separate positive and negative values
    #    positive = np.where(pq > 0, pq, np.nan)  # Use NaN for negative values
    #    negative = np.where(pq < 0, pq, np.nan)  # Use NaN for positive values
    #    ax[0].plot(times_array, positive,color='red', marker='.',label = 'mean $p_z>0$')
    #    ax[0].plot(times_array,-negative,color='blue',marker='.',label = 'mean $p_z<0$')
    # else:
    # ax[0].plot(times_array,quantity_abs,color='black',marker='.',label = 'B')

    ax[0].plot(times_array, rmsq_array, color="green", marker=".", label="$B_{rms}$")
    # ax[0].plot(times_array, mzq_array, color='black', marker='.',label='mean $|p_z|$')

    # ax[0].plot(times_array,np.abs(quantity_array[:,0]),color='blue' ,marker='.', label='$p_x$')
    # ax[0].plot(times_array,np.abs(quantity_array[:,1]),color='red'  ,marker='.', label='$p_y$')

    ax[0].legend()
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].set_ylabel("$|B|$, [$\mu G$]", fontsize=20)
    ax[0].set_xlabel("$time [Gyr]$", fontsize=20)
    #ax[0].set_ylim([0, 10.0])
    #ax[0].set_xlim([0, 3.0])
    ax[0].set_ylim([1e-6,1e2])

    ax[1].plot(times_array, npart_array, color="black", marker=".")
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_ylabel(
        f"$N_p$" + label, fontsize=20
    )  # , Lcut = {Lcut} A.U.',fontsize=20)
    ax[1].set_xlabel("$time [Gyr]$", fontsize=20)
    ax[1].set_ylim([1e1, 1e5])
    #ax[1].set_xlim([0, 3.0])

    # add panel with infromation about the run
    text_common_args = dict(
        fontsize=10, ha="center", va="center", transform=ax[2].transAxes
    )

    ax[2].text(0.5, 0.8, "Cooling halo with spin, $Np=%.0f$ " % Np, **text_common_args)
    ax[2].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    ax[2].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    ax[2].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    ax[2].text(
        0.5,
        0.4,
        kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
        **text_common_args,
    )
    ax[2].text(
        0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
    )
    ax[2].text(
        0.5,
        0.2,
        "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
        **text_common_args,
    )
    ax[2].text(
        0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
    )
    ax[2].axis("off")

    # fig.suptitle(quantity,fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)
    return


def plot_pressure_vs_time(
    all_snapshots,
    outputfileaddr,
    isvec=False,
    region_cut=True,
    label="",
    rcut=15,
    zcut=0.5,
):
    pressures_array = []
    magnetic_pressures_array = []
    dynamic_pressures_array = []
    times_array = []
    npart_array = []
    rcut *= units_dict["length"]
    zcut *= units_dict["length"]
    for snap in all_snapshots:
        pressures_values = snap["values"]["pressures"] * units_dict["pressure"]
        magnetic_pressures_values = (
            snap["values"]["magnetic_pressures"] * units_dict["pressure"]
        )
        dynamic_pressures_values = (
            snap["values"]["dynamic_pressures"] * units_dict["pressure"]
        )
        densities_values = snap["values"]["densities"] * units_dict["density"]
        masses_values = snap["values"]["masses"] * units_dict["mass"]
        Np = len(pressures_values)

        if region_cut:
            z = snap["values"]["CM_frame_coordinates"][:, 2] * units_dict["length"]
            y = snap["values"]["CM_frame_coordinates"][:, 1] * units_dict["length"]
            x = snap["values"]["CM_frame_coordinates"][:, 0] * units_dict["length"]
            r = np.sqrt(x ** 2 + y ** 2)
            mask = (r <= rcut) & (np.abs(z) <= zcut)
            pressures_values = pressures_values[mask]
            magnetic_pressures_values = magnetic_pressures_values[mask]
            dynamic_pressures_values = dynamic_pressures_values[mask]
            densities_values = densities_values[mask]
            masses_values = masses_values[mask]

        # calculate volume weighted averages:
        volume = np.sum(masses_values / densities_values)
        pressures_vv = (
            np.sum(pressures_values * masses_values / densities_values) / volume
        )
        magnetic_pressures_vv = (
            np.sum(magnetic_pressures_values * masses_values / densities_values)
            / volume
        )
        dynamic_pressures_vv = (
            np.sum(dynamic_pressures_values * masses_values / densities_values) / volume
        )

        pressures_vv = pressures_vv.to(units_dict["pressure"]).value
        magnetic_pressures_vv = magnetic_pressures_vv.to(units_dict["pressure"]).value
        dynamic_pressures_vv = dynamic_pressures_vv.to(units_dict["pressure"]).value

        pressures_array.append(pressures_vv)
        magnetic_pressures_array.append(magnetic_pressures_vv)
        dynamic_pressures_array.append(dynamic_pressures_vv)

        time = snap["parameters"]["time"]
        times_array.append(time)
        npart_array.append(len(pressures_values))

        # Retrieve some information about the simulation run
        artDiffusion = snap["metadata"]["artDiffusion"]
        dedHyp = snap["metadata"]["dedHyp"]
        dedHypDivv = snap["metadata"]["dedHypDivv"]
        dedPar = snap["metadata"]["dedPar"]
        eta = snap["metadata"]["eta"]
        git = snap["metadata"]["git"]
        gitBranch = snap["metadata"]["gitBranch"]
        hydroScheme = snap["metadata"]["hydroScheme"]
        kernel = snap["metadata"]["kernel"]
        neighbours = snap["metadata"]["neighbours"]

    pressures_array = np.array(pressures_array)
    magnetic_pressures_array = np.array(magnetic_pressures_array)
    dynamic_pressures_array = np.array(dynamic_pressures_array)
    times_array = np.array(times_array)
    npart_array = np.array(npart_array)

    fig, ax = plt.subplots(1, 3, figsize=(5.5 * 3, 5))
    ax[0].plot(
        times_array, pressures_array, color="red", marker=".", label=r"$P_{thermal}$"
    )
    ax[0].plot(
        times_array,
        magnetic_pressures_array,
        color="blue",
        marker=".",
        label=r"$P_{mag}$",
    )
    ax[0].plot(
        times_array,
        dynamic_pressures_array,
        color="green",
        marker=".",
        label=r"$\rho v^2$",
    )

    ax[0].legend()
    ax[0].grid()
    ax[0].set_yscale("log")
    ax[0].set_ylabel("$P$, [$ 10^{10}$ $M_{sol}$ $kpc^{-1}$ $Gyr^{-2}$]", fontsize=20)
    ax[0].set_xlabel("time $[Gyr]$", fontsize=20)
    ax[0].set_ylim([1e-4, 1e2])
    #ax[0].set_xlim([0, 4.0])
    # ax[0].set_ylim([1e-8,1e-4])

    ax[1].plot(times_array, npart_array, color="black", marker=".")
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_ylabel(r"$N_p^{cut}$", fontsize=20)  # , Lcut = {Lcut} A.U.',fontsize=20)
    ax[1].set_xlabel("time $[Gyr]$", fontsize=20)
    ax[1].set_ylim([1e1, 1e5])
    #ax[1].set_xlim([0, 4.0])

    # add panel with infromation about the run
    text_common_args = dict(
        fontsize=10, ha="center", va="center", transform=ax[2].transAxes
    )

    ax[2].text(0.5, 0.8, "Cooling halo with spin, $Np=%.0f$ " % Np, **text_common_args)
    ax[2].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    ax[2].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    ax[2].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    ax[2].text(
        0.5,
        0.4,
        kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
        **text_common_args,
    )
    ax[2].text(
        0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
    )
    ax[2].text(
        0.5,
        0.2,
        "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
        **text_common_args,
    )
    ax[2].text(
        0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
    )
    ax[2].axis("off")

    # fig.suptitle(quantity,fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)
    return


def plot_B_vs_time(
    all_snapshots, outputfileaddr, isvec=False, region_cut=True, label="", rcut=15, zcut=0.5
):
    magnetic_flux_densities_array = []
    times_array = []
    npart_array = []
    total_mass_array = []
    rcut *= units_dict["length"]
    zcut *= units_dict["length"]
    for snap in all_snapshots:
        magnetic_flux_densities_values = (
            snap["values"]["magnetic_flux_densities"] * units_dict["magneticfield"]
        )
        magnetic_flux_densities_values_sq = dot_vec(
            magnetic_flux_densities_values, magnetic_flux_densities_values
        )
        densities_values = snap["values"]["densities"] * units_dict["density"]
        masses_values = snap["values"]["masses"] * units_dict["mass"]
        Np = len(masses_values)
 
        if region_cut:
            z = snap["values"]["CM_frame_coordinates"][:, 2] * units_dict["length"]
            y = snap["values"]["CM_frame_coordinates"][:, 1] * units_dict["length"]
            x = snap["values"]["CM_frame_coordinates"][:, 0] * units_dict["length"]
            r = np.sqrt(x ** 2 + y ** 2)
            mask = (r <= rcut) & (np.abs(z) <= zcut)
            magnetic_flux_densities_values_sq = magnetic_flux_densities_values_sq[mask]
            densities_values = densities_values[mask]
            masses_values = masses_values[mask]

        # calculate volume weighted averages:
        volume = np.sum(masses_values / densities_values)
        magnetic_flux_densities_sq_vv = (
            np.sum(magnetic_flux_densities_values_sq * masses_values / densities_values)
            / volume
        )
        magnetic_flux_densities_vv = np.sqrt(magnetic_flux_densities_sq_vv)

        magnetic_flux_densities_vv = magnetic_flux_densities_vv.to(
            units_dict["magneticfield"]
        ).value

        magnetic_flux_densities_array.append(magnetic_flux_densities_vv)

        time = snap["parameters"]["time"]
        times_array.append(time)
        npart_array.append(len(masses_values))
        total_mass_array.append(np.sum(masses_values.to(units_dict["mass"]).value))

        # Retrieve some information about the simulation run
        artDiffusion = snap["metadata"]["artDiffusion"]
        dedHyp = snap["metadata"]["dedHyp"]
        dedHypDivv = snap["metadata"]["dedHypDivv"]
        dedPar = snap["metadata"]["dedPar"]
        eta = snap["metadata"]["eta"]
        git = snap["metadata"]["git"]
        gitBranch = snap["metadata"]["gitBranch"]
        hydroScheme = snap["metadata"]["hydroScheme"]
        kernel = snap["metadata"]["kernel"]
        neighbours = snap["metadata"]["neighbours"]

    magnetic_flux_densities_array = np.array(magnetic_flux_densities_array)
    times_array = np.array(times_array)
    npart_array = np.array(npart_array)
    total_mass_array = np.array(total_mass_array)

    fig, ax = plt.subplots(1, 4, figsize=(5.5 * 4, 5))
    ax[0].plot(times_array, magnetic_flux_densities_array, color="black", marker=".")
    print(magnetic_flux_densities_array)
    # ax[0].legend()
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].set_ylabel("$B_{rms}$, $[\mu G]$", fontsize=20)
    ax[0].set_xlabel("time $[Gyr]$", fontsize=20)
    #ax[0].set_ylim([0, 3.5])
    #ax[0].set_xlim([0, 3.0])
    ax[0].set_ylim([1e-4,1e1])

    ax[1].plot(times_array, npart_array, color="black", marker=".")
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_ylabel(r"$N_p^{cut}$", fontsize=20)  # , Lcut = {Lcut} A.U.',fontsize=20)
    ax[1].set_xlabel("time $[Gyr]$", fontsize=20)
    ax[1].set_ylim([1e1, 1e6])
    #ax[1].set_xlim([0, 3.0])

    ax[2].plot(times_array, total_mass_array, color="black", marker=".")
    ax[2].grid()
    ax[2].set_yscale("log")
    ax[2].set_ylabel(
        r"$M_{tot} [M_{sol}]$", fontsize=20
    )  # , Lcut = {Lcut} A.U.',fontsize=20)
    ax[2].set_xlabel("time $[Gyr]$", fontsize=20)
    ax[2].set_ylim([1e8, 1e11])
    ax[2].set_xlim([0, 3.0])

    # add panel with infromation about the run
    text_common_args = dict(
        fontsize=10, ha="center", va="center", transform=ax[3].transAxes
    )

    ax[3].text(0.5, 0.8, "Cooling halo with spin, $Np=%.0f$ " % Np, **text_common_args)
    ax[3].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    ax[3].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    ax[3].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    ax[3].text(
        0.5,
        0.4,
        kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
        **text_common_args,
    )
    ax[3].text(
        0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
    )
    ax[3].text(
        0.5,
        0.2,
        "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
        **text_common_args,
    )
    ax[3].text(
        0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
    )
    ax[3].axis("off")

    # fig.suptitle(quantity,fontsize=20)
    fig.tight_layout()
    plt.savefig(outputfileaddr)
    plt.close(fig)
    return


# 34013

folder = "./MHD_eta=1.595_noB_antisym_fd"

snapshots_history, snapshot_names = updoad_shapshots_data_from_folder(folder)

# plot_quantity_vs_time(snapshots_history,'magnetic_flux_densities_norm',outputfileaddr=folder+'/momentums_noRegCut.png',isvec=True,region_cut=False, label = ' All')

plot_pressure_vs_time(
    snapshots_history,
    outputfileaddr=folder + "/P_cut.png",
    isvec=True,
    region_cut=True,
    rcut=15,
    zcut=15.0,
    label=f" $Lbox=${2*Lcut} kPc",
)
plot_B_vs_time(
    snapshots_history,
    outputfileaddr=folder + "/B_cut.png",
    isvec=True,
    region_cut=True,
    rcut=15,
    zcut=1.0,
    label=f" $Lbox=${2*Lcut} kPc",
)

#fast_particles = identify_largest_quantity_particles(
#    snapshots_history, "velocities", isvec=False, region_cut=True
#)
#print('Fast particles: ',fast_particles)

#particle_history = upload_particle_history(snapshots_history,38770)

#plot_quatities_for_particle_vs_time(
#    particle_history,
#    ["CM_frame_coordinates","CM_frame_z","velocities","magnetic_flux_densities","magnetic_divergences"],
#    outputfileaddr=folder+"/quantities.png",
#    time_var="time",
#    timeinterval=[0, 1],
#    title="quantities vs time",
#    customtimeinterval=False,
#    customquantityinterval=False,
#    quantityinterval=[1e-1, 1e1],
#    logscale=True,
#)
