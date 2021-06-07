import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import swiftsimio
import sys
import glob
import unyt
import numpy as np

try:
    plt.style.use("../../../tools/stylesheets/mnras.mplstyle")
except:
    print("Can't find Matplotlib stylesheet.")

## read command line arguments
snapshot_name = sys.argv[1]

################# Read in observational data

data_schaye = np.genfromtxt(
    "./datasets/schaye_et_al_2000_thermal_history.dat", skip_header=2
)
data_walther = np.genfromtxt(
    "./datasets/walther_et_al_2019_thermal_history.dat", skip_header=2
)

data_schaye = data_schaye.T
data_walther = data_walther.T

schaye_z_lower_error = data_schaye[0] - data_schaye[1]
schaye_z_upper_error = data_schaye[2] - data_schaye[0]

schaye_T_lower_error = np.log10(data_schaye[3] * 1.0e4) - np.log10(
    data_schaye[4] * 1.0e4
)
schaye_T_upper_error = np.log10(data_schaye[5] * 1.0e4) - np.log10(
    data_schaye[3] * 1.0e4
)
walther_T_lower_error = np.log10(data_walther[1] * 1.0e4) - np.log10(
    data_walther[2] * 1.0e4
)
walther_T_upper_error = np.log10(data_walther[3] * 1.0e4) - np.log10(
    data_walther[1] * 1.0e4
)

############### Read in simulation data
density_units = "Solar_Mass / Mpc**3"
temperature_units = "K"

## First, get list of all snapshots
reg_exp = "%s*.hdf5" % snapshot_name
snap_list = glob.glob(reg_exp)
snap_data = sorted(
    [swiftsimio.load(snap) for snap in snap_list], key=lambda x: x.metadata.z
)

z = np.empty(len(snap_data), dtype=np.float32)
T_mean = unyt.unyt_array(np.empty_like(z), units=temperature_units)
T_std = unyt.unyt_array(np.empty_like(z), units=temperature_units)
rho_mean = unyt.unyt_array(np.empty_like(z), units=density_units)
rho_std = unyt.unyt_array(np.empty_like(z), units=density_units)

## loop through list
for index, data in enumerate(snap_data):
    # Stick redshift in a list
    z[index] = data.metadata.z

    # Convert units to something sensible so `np.mean` doesn't freak out
    data.gas.temperatures.convert_to_units(temperature_units)
    data.gas.densities.convert_to_units(density_units)

    # Get mean and standard deviation of temperature
    T_mean[index] = np.mean(data.gas.temperatures)
    T_std[index] = np.std(data.gas.temperatures)

    # Get mean and standard deviation of density
    rho_mean[index] = np.mean(data.gas.densities)
    rho_std[index] = np.std(data.gas.densities)


## Put Density into units of mean baryon density today

# first calculate rho_bar_0 from snapshot metadata
data = snap_data[0]
cosmology = data.metadata.cosmology
H0 = unyt.unyt_quantity.from_astropy(cosmology.H0)
Omega_bar = unyt.unyt_quantity(cosmology.Ob0)

### now calculate rho_bar_0 and divide through
rho_bar_0 = 3.0 * H0 ** 2 / (8.0 * np.pi * unyt.G) * Omega_bar
rho_mean /= rho_bar_0
rho_std /= rho_bar_0

### from the first snapshot, get code information
code_info = data.metadata.code
git_branch = code_info["Git Branch"].decode("UTF-8")
git_revision = code_info["Git Revision"].decode("UTF-8")

hydro_metadata = data.metadata.hydro_scheme
scheme = hydro_metadata["Scheme"].decode("UTF-8")
params = data.metadata.parameters

subgrid_metadata = data.metadata.subgrid_scheme
cooling_model = subgrid_metadata["Cooling Model"].decode("UTF-8")
chemistry_model = subgrid_metadata["Chemistry Model"].decode("UTF-8")

if cooling_model == "EAGLE":
    z_r_H = float(params["EAGLECooling:H_reion_z"])
    H_heat_input = float(params["EAGLECooling:H_reion_eV_p_H"])
    z_r_He_centre = float(params["EAGLECooling:He_reion_z_centre"])
    z_r_He_sigma = float(params["EAGLECooling:He_reion_z_sigma"])
    He_heat_input = float(params["EAGLECooling:He_reion_eV_p_H"])

if chemistry_model == "EAGLE":
    metallicity = float(params["EAGLEChemistry:init_abundance_metal"])
else:
    metallicity = "Unknown"

# Make plot of temperature evolution  --------------------------------
fig, ax = plt.subplots(constrained_layout=True)

# Plot sim properties
if cooling_model == "EAGLE":
    ax.axvline(z_r_H, color="k", linestyle="--", alpha=0.5, lw=0.7, zorder=-1)
    ax.text(
        z_r_H + 0.1, 3.55, "H reion.", rotation=90, alpha=0.5, fontsize=7, va="bottom"
    )
    ax.axvline(z_r_He_centre, color="k", linestyle="--", alpha=0.5, lw=0.7, zorder=-1)
    ax.text(
        z_r_He_centre + 0.1,
        3.55,
        "HeII reion.",
        rotation=90,
        alpha=0.5,
        fontsize=7,
        va="bottom",
    )

# Plot observational data
ax.errorbar(
    data_schaye[0],
    np.log10(data_schaye[3] * 1.0e4),
    xerr=[schaye_z_lower_error, schaye_z_upper_error],
    yerr=[schaye_T_lower_error, schaye_T_upper_error],
    fmt="s",
    mec="0.3",
    color="0.3",
    markersize=4,
    markeredgewidth=0.5,
    linewidth=0.5,
    mfc="w",
    label="Schaye et al. (2000)",
)
ax.errorbar(
    data_walther[0],
    np.log10(data_walther[1] * 1.0e4),
    yerr=[walther_T_lower_error, walther_T_upper_error],
    fmt=".",
    mec="0.3",
    color="0.3",
    markersize=7,
    markeredgewidth=0.5,
    linewidth=0.5,
    mfc="w",
    label="Walther et al. (2019)",
)

# Plot simulation
ax.plot(z, np.log10(T_mean))

# Legend
ax.legend(loc="upper right", frameon=True)
ax.text(
    0.025,
    0.975,
    "SWIFT %s \nCooling model: %s \n$\\rho=%.3f$ $\\Omega_{b}\\rho_{crit,0}$\n$Z=%s$"
    % (git_revision, cooling_model, rho_mean[-1], metallicity),
    va="top",
    ha="left",
    bbox={"fc": "white", "ec": "none"},
    transform=plt.gca().transAxes,
    zorder=1,
)

ax.set_xlim(0, 12.2)
ax.set_ylim(3.5, 4.85)
ax.set_xlabel("Redshift", labelpad=-0.5)
ax.set_ylabel(r"$\log_{10}(T/K)$", labelpad=0)
fig.savefig("Temperature_evolution.png")


# Make plot of density evolution  --------------------------------
fig, ax = plt.subplots(constrained_layout=True)

ax.text(0.2, 1.011, "SWIFT %s" % git_revision, va="top", ha="left", fontsize=8)

ax.fill_between(z, rho_mean - rho_std, rho_mean + rho_std, alpha=0.5)
ax.plot(z, rho_mean)

ax.axhline(y=1.0, linestyle="--", color="k", alpha=0.5, lw=1.0)

ax.set_xlim(0.0, 12.2)
ax.set_ylabel(r"$\delta_b = \rho / \Omega_b\rho_{crit,0}$", labelpad=0.0)
ax.set_ylim(0.988, 1.012)
ax.set_yticks([0.99, 1.0, 1.01])
ax.set_xlabel("Redshift", labelpad=-0.5)
fig.savefig("Density_evolution.png")
