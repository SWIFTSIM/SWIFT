import os
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import swiftsimio as sw
from unyt.physical_constants import mh, kboltz_cgs as k_B
import unyt
import re

# --- GLOBAL STATE VARIABLES ---
GLOBAL_DATA_LIST = []
GLOBAL_LINES = {}
GLOBAL_TITLES = {}
GLOBAL_SUPERTITLE = None
GLOBAL_LOG_FILE = None
GLOBAL_ISM_DENSITY = None
GLOBAL_LIFETIME = None


def get_gas_mu(data: sw.SWIFTDataset) -> np.array:
    """
    Return the mean molecular weight of the gas.
    """
    # Use unyt's physical constants
    from unyt.physical_constants import mh

    # Get the cooling model
    cooling = data.metadata.subgrid_scheme["Cooling Model"]

    # Use a variable for better readability
    gas = data.gas

    # Get rho
    rho = gas.densities.to(unyt.g / unyt.cm ** 3)

    # hydrogen mass in gram
    mh.convert_to_cgs()
    mH_in_g = mh

    if cooling == b"Grackle3":  # grackle3
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nH2I = gas.h2_i * rho / (2 * mH_in_g)
        nH2II = gas.h2_ii * rho / (2 * mH_in_g)
        nHDI = gas.hdi * rho / (3 * mH_in_g)

        nel = nHII + nHeII + 2 * nHeIII + nH2II
        mu = (
            (nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHDI * 3
        ) / (nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nHDI + nel)
        return mu

    elif cooling == b"Grackle2":  # grackle2
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nH2I = gas.h2_i * rho / (2 * mH_in_g)
        nH2II = gas.h2_ii * rho / (2 * mH_in_g)

        nel = nHII + nHeII + 2 * nHeIII + nH2II
        mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2) / (
            nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nel
        )
        return mu

    elif cooling == b"Grackle1":  # grackle1
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nel = nHII + nHeII + 2 * nHeIII
        mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4) / (
            nHI + nHII + nHeI + nHeII + nHeIII + nel
        )
        return mu

    else:  # Grackle0
        from unyt.physical_constants import kboltz_cgs as k_B_cgs

        gamma = data.metadata.gas_gamma[0]

        # Get internal energy
        u = data.gas.internal_energies
        u = u.to_physical()
        u = u.to(unyt.erg / unyt.g)

        # Get hydrigen fraction
        if data.metadata.policy["cooling"] == 1:
            H_frac = float(
                data.metadata.parameters["GrackleCooling:HydrogenFractionByMass"]
            )
        else:
            H_frac = 1.0

        # Compute T/mu
        T_over_mu = (gamma - 1.0) * u.value * mH_in_g.value / k_B_cgs.value
        T_trans = 1.1e4
        mu_trans = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))

        # Determine if we are ionized or not
        mu = np.ones(np.size(u))
        mask_ionized = T_over_mu > (T_trans + 1) / mu_trans
        mask_neutral = T_over_mu < (T_trans + 1) / mu_trans

        # Give the right mu
        mu[mask_ionized] = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))
        mu[mask_neutral] = 4.0 / (1.0 + 3.0 * H_frac)

        return mu


def get_gas_temperatures(data: sw.SWIFTDataset) -> np.array:
    """
        Compute the temperature of the gas.
    """

    # Use unyt's physical constants
    from unyt.physical_constants import mh, kboltz_cgs as k_B

    # Convert to cgs
    mh.convert_to_cgs()

    # Get the cooling model
    cooling = data.metadata.subgrid_scheme["Cooling Model"]

    # Get internal energy and convert to physical units in cgs
    u = data.gas.internal_energies
    u = u.to_physical()
    u = u.to(unyt.erg / unyt.g)

    # Get gamma and compute mu
    gamma = data.metadata.gas_gamma[0]

    mu = get_gas_mu(data)

    # Finally compute the temperature
    if cooling == b"Grackle3" or cooling == b"Grackle2" or cooling == b"Grackle1":
        T = mu * (gamma - 1.0) * u * mh / k_B
    else:
        a = (gamma - 1.0) * (mu * mh) / k_B * u
        T = np.where((u.value > 0), a.value, 0) * unyt.kelvin
    return T


def get_data(filename):
    """
    Read the data and convert the properties.
    Returns: (r, v_r, densities, temperatures, time) or None if loading fails.
    """
    from unyt.physical_constants import mh, kboltz_cgs as k_B

    global GLOBAL_ISM_DENSITY, GLOBAL_LIFETIME

    # read files
    try:
        data = sw.load(filename)
    except Exception as e:
        print(f"Error loading file {filename} with swiftsimio: {e}")
        return None

    if data is None:
        print(f"The data file {filename} provided is not suitable or returned None.")
        return None

    # convert units
    data.gas.coordinates.convert_to_units("kpc")
    data.gas.densities.convert_to_units("kg/cm**3")
    data.gas.temperatures = get_gas_temperatures(data)
    data.gas.temperatures.convert_to_cgs()
    data.gas.velocities.convert_to_units("km/s")
    data.metadata.time.convert_to_units("Myr")

    # Get the center
    center = data.stars.coordinates

    # adjust the data
    coordinates = data.gas.coordinates - center

    # Calculate log density (H/cm^3) - using unyt's mass of hydrogen
    # and converting density to the desired unit before logging.
    density_H_cm3 = data.gas.densities.value / mh
    densities = np.log10(density_H_cm3)

    # Calculate log temperature
    temperatures = np.log10(data.gas.temperatures)
    velocities = data.gas.velocities
    time = data.metadata.time
    mu = get_gas_mu(data)
    mu = np.mean(mu)

    if time == 0.0:
        GLOBAL_ISM_DENSITY = np.mean(density_H_cm3.value)
        data.stars.birth_times.convert_to_units("myr")
        GLOBAL_LIFETIME += data.stars.birth_times.value[0]

    # Radial coordinate and radial velocity calculation (copied exactly)
    r = np.sqrt(
        coordinates[:, 0] ** 2 + coordinates[:, 1] ** 2 + coordinates[:, 2] ** 2
    )

    dot_product = np.sum(coordinates * velocities, axis=1)

    r_nonzero = np.where(r != 0, r, 1.0 * unyt.kpc)
    v_r = dot_product / r_nonzero
    v_r = np.where(r != 0, v_r, 0.0 * unyt.km / unyt.s)

    # print("r",r,"v_r",v_r,"densities",densities,"temperature",temperatures,"mu",mu,"time",time)
    return r, v_r, densities, temperatures, mu, time


def extract_stellar_feedback(log_file_path: str):
    """
    Reads a text file and extracts the total cumulative energy and mass ejected 
    from stars by summing up the rates (Energy[erg/yr] and Mass_ejected[Msol/yr])
    found in each matching log line.

    Expected line format: 
    [TIME] ... Energy[erg/yr]=<VALUE_E> Mass_ejected[Msol/yr]=<VALUE_M>
    
    Returns: (total_energy, total_mass) in the cumulative units (erg, Msol), 
             or (None, None) if the file cannot be read.
    """

    global GLOBAL_LIFETIME

    if not os.path.exists(log_file_path):
        print(f"Error: Stellar feedback log file not found at path: {log_file_path}")
        return None, None

    # Regex pattern to capture the number (float, including scientific notation)
    # This pattern captures numbers like 123.4, 1.23e50, -0.1e-10
    float_pattern = r"([+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?)"

    # Define compiled patterns specific to the user's log format (rates)
    # 1. Energy pattern targets: "Energy[erg/yr]=" followed by the float
    energy_pattern = re.compile(r"Energy\[erg/yr\]=" + float_pattern)

    # 2. Mass pattern targets: "Mass_ejected[Msol/yr]=" followed by the float
    mass_pattern = re.compile(r"Mass_ejected\[Msol/yr\]=" + float_pattern)

    # 3. Lifetime pattern targets: "lifetime_myr=" followed by the float
    lifetime_pattern = re.compile(r"lifetime_myr=" + float_pattern)

    # Initialize quantities
    energy_rate = 0.0
    mass_loss = 0.0

    print(
        f"\nSearching for stellar feedback data (cumulative rates) in {log_file_path}..."
    )

    try:
        with open(log_file_path, "r") as f:
            for line in f:

                # Reset temporary variables for the current line check
                current_energy = None
                current_mass = None

                # Search for energy rate
                match_e = energy_pattern.search(line)
                if match_e:
                    current_energy = float(match_e.group(1))

                # Search for mass rate
                match_m = mass_pattern.search(line)
                if match_m:
                    current_mass = float(match_m.group(1))

                # Search for star death
                match_time = lifetime_pattern.search(line)
                if match_time:
                    GLOBAL_LIFETIME = float(match_time.group(1))

                # Check if BOTH quantities were found on this single line
                if current_energy is not None and current_mass is not None:
                    energy_rate = current_energy
                    mass_loss = current_mass

                if GLOBAL_LIFETIME is not None:
                    break

    except Exception as e:
        print(f"Error reading or parsing feedback file: {e}")
        return None, None, None

    # if lifetime not found, it means the star is still alive
    if GLOBAL_LIFETIME is None:
        GLOBAL_LIFETIME = 0.0

    # Check if any significant data was found
    if energy_rate > 0.0 or mass_loss > 0.0:
        print(
            f"Extracted Cumulative Feedback: Energy = {energy_rate:.3e} erg/yr, Mass = {mass_loss:.3e} Msol/yr, Star_Lifetime = {GLOBAL_LIFETIME:.3e} Myr"
        )
    else:
        print(
            "Warning: Log file read, but no positive energy or mass ejection rates were found."
        )

    return energy_rate, mass_loss


def calculate_snowplow_phase(energy_rate, time, mu, wind_velocity, mass_loss):
    """
    Calculate the theoretical radius, velocity, density and tempereture
    of the bubble's shockwave for the constant energy's snowplow phase and constant's momentum snowplow phase.
    """

    global GLOBAL_ISM_DENSITY

    yr_to_s = 1 * unyt.yr
    yr_to_s.convert_to_cgs()
    yr_to_s = yr_to_s.value
    cm_to_kpc = 1 * unyt.cm
    cm_to_kpc.convert_to_units("kpc")
    cm_to_kpc = cm_to_kpc.value
    unyt.mh.convert_to_cgs()
    mH_in_g = unyt.mh.value

    energy_rate.convert_to_cgs()
    time.convert_to_cgs()
    # Some variable used in the book "Introduction to Stellar Winds-Cambridge University Press (1999)"
    L36 = energy_rate / 1e36
    L36.convert_to_cgs()
    t6 = time / 1e6
    t6.convert_to_units("yr")
    v3 = wind_velocity / 1000
    v3.convert_to_units("km/s")

    if time == 0.0:
        print("energy", energy_rate, "L36", L36, "t6", t6, "v3", v3)
        # log(t=0) case...
        return 0, 0, 0, 0, v3.value, v3.value, 0, GLOBAL_ISM_DENSITY

    # Formula in papers (arxiv.org/pdf/2211.15705)
    R_e = (
        0.88
        * (energy_rate.value / mu / GLOBAL_ISM_DENSITY / mH_in_g) ** (1 / 5)
        * (time.value) ** (3 / 5)
    )
    R_e *= cm_to_kpc
    R_time_end = (
        0.88
        * (energy_rate.value / mu / GLOBAL_ISM_DENSITY / mH_in_g) ** (1 / 5)
        * (GLOBAL_LIFETIME * 10 ** 6 * yr_to_s) ** (3 / 5)
        * cm_to_kpc
    )
    R_p = (
        (3 / 2 / np.pi * mass_loss * wind_velocity / mu / GLOBAL_ISM_DENSITY / mH_in_g)
        ** (1 / 4)
        * (time.value) ** (1 / 2)
        * cm_to_kpc
    ).value  # + R_time_end

    # Formula from the book "Introduction to Stellar Winds-Cambridge University Press (1999)"
    R_b_e = (
        28 * (L36 / mu / GLOBAL_ISM_DENSITY) ** (1 / 5) * t6.value ** (3 / 5)
    ).value / 1000
    R_b_time_end = (
        28
        * (L36 / mu / GLOBAL_ISM_DENSITY) ** (1 / 5)
        * (GLOBAL_LIFETIME * 10 ** 6 * yr_to_s) ** (3 / 5)
    ).value * cm_to_kpc
    R_b_p = (
        5.0 * (L36 * v3 / GLOBAL_ISM_DENSITY) ** (1 / 4) * (t6.value) ** (1 / 2)
    ).value / 1000  # + R_b_time_end

    v_e = (17 * (L36 / mu / GLOBAL_ISM_DENSITY) ** (1 / 5) * t6.value ** (-2 / 5)).value
    v_time_end = (
        17
        * (L36 / mu / GLOBAL_ISM_DENSITY) ** (1 / 5)
        * (GLOBAL_LIFETIME * 10 ** 6 * yr_to_s) ** (-2 / 5)
    ).value
    v_p = (
        7.9 * (L36 / v3 / GLOBAL_ISM_DENSITY) ** (1 / 4) * (t6.value) ** (-1 / 2)
    ).value  # + v_time_end

    T = np.log10(
        1.6e6 * L36 ** (8 / 35) * GLOBAL_ISM_DENSITY ** (2 / 35) * t6.value ** (-6 / 35)
    )

    n = np.log10(
        0.01
        * L36 ** (6 / 35)
        * GLOBAL_ISM_DENSITY ** (19 / 35)
        * t6.value ** (-22 / 35)
    )

    return R_e, R_p, R_b_e, R_b_p, v_e, v_p, T, n


# ----------------------------------------------------------------------------------
# --- CORE ANIMATION & PLOTTING LOGIC ---
# ----------------------------------------------------------------------------------


def update_plot_scatter(frame_index):
    """
    Drawing function called by FuncAnimation. Updates the data in the plot lines
    for the movie generation.
    """
    global GLOBAL_SUPERTITLE, GLOBAL_LINES, GLOBAL_DATA_LIST, GLOBAL_LIFETIME, GLOBAL_TEXT
    r, v_raw, log_rho_raw, log_T_raw, _, time = GLOBAL_DATA_LIST[frame_index]

    # Ensure time is in Myr for display, as unyt object are passed by reference in python....
    time.convert_to_units("Myr")

    GLOBAL_LINES["velocity"].set_data(r.value, v_raw.value)
    GLOBAL_LINES["density"].set_data(r.value, log_rho_raw)
    GLOBAL_LINES["temperature"].set_data(r.value, log_T_raw)

    # Show the "THE STAR DIED!" text if the star's lifetime has been reached
    if time > GLOBAL_LIFETIME and GLOBAL_LIFETIME != 0.0:
        GLOBAL_TEXT.set_visible(True)

    GLOBAL_SUPERTITLE.set_text(f"Time: {time.value:.4f} Myr")

    return (
        [line for line in GLOBAL_LINES.values()] + [GLOBAL_SUPERTITLE] + [GLOBAL_TEXT]
    )


def write_quantities(
    coordinates,
    velocities,
    densities,
    temperatures,
    mu,
    frame_index,
    time,
    energy_rate,
    wind_velocity,
    mass_loss,
):
    """
    Analyze the data to find peaks and log them to a file.
    """
    global GLOBAL_LOG_FILE, GLOBAL_ISM_DENSITY, GLOBAL_LIFETIME

    frame_str = f"{frame_index:04d}"

    # --- ANALYSIS PART (Peak Finding and theoretical Logging) ---
    peak_data = f"{time.value:.4f}"
    r_val = coordinates.value  # r [kpc]

    # Density Peak
    # Find the index of the maximum density value
    rho_max = np.max(densities)
    idx_rho_max = np.argmin(
        np.where(densities == rho_max, r_val, np.inf)
    )  # In case of multiple maxima, take the smallest radius, because at the begining the all have the same values
    r_rho_max = r_val[idx_rho_max]
    peak_data += f"\t{r_rho_max:.4f}\t{densities[idx_rho_max]:.4f}"

    # Velocity Peak (Using argmax on the magnitude of the velocity array values)
    idx_v_max = np.argmax(np.abs(velocities.value))
    r_v_max = r_val[idx_v_max]
    peak_data += f"\t{r_v_max:.4f}\t{velocities[idx_v_max].value:.4f}"

    # Velocity at shockwave
    peak_data += f"\t{velocities[idx_rho_max].value:.4f}"

    # Temperature Peak
    idx_temp_max = np.argmax(temperatures)
    r_temp_max = r_val[idx_temp_max]
    peak_data += f"\t{r_temp_max:.4f}\t{temperatures[idx_temp_max]:.4f}"

    # Temperature at shockwave
    temp_at_rho_peak = temperatures[idx_rho_max]
    temp_at_v_peak = temperatures[idx_v_max]
    peak_data += f"\t{temp_at_rho_peak:.4f}\t{temp_at_v_peak:.4f}"

    # theoretical part
    R_e, R_p, R_b_e, R_b_p, v_e, v_p, T_e, n_e, = calculate_snowplow_phase(
        energy_rate, time, np.mean(mu), wind_velocity, mass_loss
    )

    peak_data += f"\t{R_e:.4f}\t{R_p:.4f}\t{R_b_e:.4f}\t{R_b_p:.4f}\t{v_e:.4f}\t{v_p:.4f}\t{T_e:.4f}\t{n_e:.4f}"

    peak_data += f"\t{GLOBAL_LIFETIME:.4f}"

    if GLOBAL_LOG_FILE:
        GLOBAL_LOG_FILE.write(peak_data + "\n")


# ----------------------------------------------------------------------------------
# --- MAIN EXECUTION FUNCTION ---
# ----------------------------------------------------------------------------------


def main(file_list, log_path):
    # Use global variables
    global GLOBAL_DATA_LIST, GLOBAL_LINES, GLOBAL_TITLES, GLOBAL_TEXT, GLOBAL_SUPERTITLE, GLOBAL_LOG_FILE, GLOBAL_ISM_DENSITY, GLOBAL_LIFETIME

    print(
        f"Loading files and getting graph limits from {len(file_list)} potential files..."
    )

    # Initialize Limits
    R_MAX = 0.0
    Y_MIN_MAX = {
        "velocity": [np.inf, -np.inf],
        "density": [np.inf, -np.inf],
        "temperature": [np.inf, -np.inf],
    }

    data_list_temp = []  # Temporary list to store loaded data

    # Get stellar feedback values from the swift log file
    energy_rate, mass_loss = extract_stellar_feedback(log_path)

    # convert to unyt quantities
    energy_rate = energy_rate * unyt.erg / unyt.yr
    mass_loss = mass_loss * unyt.Msun / unyt.yr

    wind_velocity = np.sqrt(2 * energy_rate / mass_loss)
    wind_velocity.convert_to_units("km/s")
    print(
        f"Get Energy rate: {energy_rate}, Mass-Loss: {mass_loss}, and calculated wind velocity: {wind_velocity:.3f}"
    )

    # PASS 1: DATA LOADING AND LIMITS CALCULATION
    for file_path in tqdm(file_list, desc="Limit and Data Loading"):
        data_tuple = get_data(file_path)

        # Only proceed if data loading was successful
        if data_tuple is not None:
            r, v_r, densities, temperatures, _, _ = data_tuple

            # Calculate Limits
            R_MAX = max(R_MAX, r.max().value)
            Y_MIN_MAX["velocity"][0] = min(Y_MIN_MAX["velocity"][0], v_r.min().value)
            Y_MIN_MAX["velocity"][1] = max(Y_MIN_MAX["velocity"][1], v_r.max().value)
            Y_MIN_MAX["density"][0] = min(Y_MIN_MAX["density"][0], densities.min())
            Y_MIN_MAX["density"][1] = max(Y_MIN_MAX["density"][1], densities.max())
            Y_MIN_MAX["temperature"][0] = min(
                Y_MIN_MAX["temperature"][0], temperatures.min()
            )
            Y_MIN_MAX["temperature"][1] = max(
                Y_MIN_MAX["temperature"][1], temperatures.max()
            )

            # Store the successfully loaded data
            data_list_temp.append(data_tuple)

    # Store valid data globally
    GLOBAL_DATA_LIST = data_list_temp

    # spacing limit for readability
    for name in Y_MIN_MAX.keys():
        if Y_MIN_MAX[name][0] > 0:
            Y_MIN_MAX[name][0] *= 0.95
        else:
            Y_MIN_MAX[name][0] *= 1.05

        if Y_MIN_MAX[name][1] > 0:
            Y_MIN_MAX[name][1] *= 1.05
        else:
            Y_MIN_MAX[name][1] *= 0.95

    R_MAX *= 1.05

    print(f"Successfully loaded data for {len(GLOBAL_DATA_LIST)} frames.")

    if not GLOBAL_DATA_LIST:
        print("No suitable data loaded. Exiting.")
        return

    # --- Open Log File ---
    log_filename = "radial_peak_positions.txt"
    try:
        GLOBAL_LOG_FILE = open(log_filename, "w")
        # Write header: Time, R_RhoMax, RhoMax, R_TMax, TMax, R_VMax, VMax
        header = "# Time [Myr]\tR_RhoMax [kpc]\tRhoMax [H/cm^3]\tR_VMax [kpc]\tVMax [km/s]\tV_R_RhoMAX [km/s]\tR_TMax [kpc]\tTMax [K]\tT_R_RhoMAX [K]\tT_R_vMAX [K]\tR_e [kpc]\tR_p [kpc]\tR_b_e [kpc]\tR_b_p [kpc]\tv_e [km/s]\tv_p [kpc]\tT_e [k]\tRho_e [H/cm^3]\tt_end [Myr]\n"
        GLOBAL_LOG_FILE.write(header)
        print(f"Logging peak data to {log_filename}")
    except Exception as e:
        print(f"Warning: Could not open log file {log_filename}: {e}")
        GLOBAL_LOG_FILE = None

    # --- Initialize Figure for FuncAnimation (Movie) ---
    fig, axes = plt.subplots(3, 1, figsize=(6, 12), dpi=300, sharex=True)
    plt.tight_layout(rect=[0.07, 0.03, 1, 0.9])
    plt.subplots_adjust(hspace=0.05)
    GLOBAL_SUPERTITLE = fig.suptitle(f"Time: 0.00 Myr", fontsize=14, y=0.98)
    GLOBAL_TEXT = fig.text(
        0.5,  # x position (center)
        0.92,  # y position (just above the top plot)
        "THE STAR DIED !",  # Text content
        fontsize=18,
        color="red",
        fontweight="bold",
        ha="center",  # Horizontal alignment
        visible=False,  # Hidden initially
    )

    quantities_setup = {
        "velocity": "Radial velocities [km/s]",
        "density": "log(Rho) [H/cm^3]",
        "temperature": "log(T) [K]",
    }

    for i, (name, label) in enumerate(quantities_setup.items()):
        ax = axes[i]

        ax.set_ylabel(label)
        ax.set_xlim(0, R_MAX)
        ax.set_ylim(Y_MIN_MAX[name][0], Y_MIN_MAX[name][1])

        # Only show X-axis labels and ticks on the BOTTOM plot
        if i < 2:
            ax.tick_params(axis="x", labelbottom=False)
            ax.set_xlabel("")  # Remove the label from the top two plots
        else:
            ax.set_xlabel("r [kpc]")  # Only set the label on the bottom plot

        line, = ax.plot([], [], "k.", markersize=1.0, alpha=0.5)

        GLOBAL_LINES[name] = line
        # GLOBAL_TITLES[name] = ax.set_title(f"Time: 0.00 Myr")

    # PASS 2: MANUAL FRAME GENERATION, PNG SAVING, AND LOGGING
    print("\nLogging peak data...")

    for i, data_tuple in tqdm(
        enumerate(GLOBAL_DATA_LIST),
        desc="Frame Gen/Logging",
        total=len(GLOBAL_DATA_LIST),
    ):
        r, v_r, densities, temperatures, mu, time = data_tuple

        write_quantities(
            r,
            v_r,
            densities,
            temperatures,
            mu,
            i,
            time,
            energy_rate,
            wind_velocity,
            mass_loss,
        )

    print("\nPeak data logging complete.")

    # --- Close Log File ---
    if GLOBAL_LOG_FILE:
        GLOBAL_LOG_FILE.close()
        print(f"Log file {log_filename} closed.")

    # --- Generate Movie (Requires PIL/Pillow writer) ---
    print("\nGenerating the GIF Movie...")

    output_filename = "radial_profile_movie.gif"
    writer = "pillow"  # Uses pillow (PIL) writer to avoid ffmpeg dependency

    animate = FuncAnimation(
        fig,
        update_plot_scatter,
        frames=len(GLOBAL_DATA_LIST),
        interval=200,
        blit=False,
        repeat=False,
    )

    try:
        animate.save(output_filename, writer=writer, dpi=300)
        print(f"\nMovie saved successfully as {output_filename}")
    except Exception as e:
        print(f"\nAn error occurred during saving the animation: {e}")
        print("HINT: Ensure the 'Pillow' library is installed (pip install Pillow).")


# --- ARGPARSE AND ENTRY POINT ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate radial profile analysis, movie, and log peak positions from simulation snapshot files."
    )
    # Define the command-line argument for the data directory
    parser.add_argument(
        "data_dir",
        type=str,
        help="Path to the directory containing the simulation. /!\ not the snap/ subdirectory, but the directory before that (because the script also search for the log file there).",
    )
    # Define the command-line argument for the file pattern
    parser.add_argument(
        "--pattern",
        type=str,
        default="snapshot_*.hdf5",
        help="Glob pattern for snapshot files (e.g., 'snapshot_*.hdf5'). Default is 'snapshot_*.hdf5'.",
    )

    args = parser.parse_args()

    # Dynamically find all files matching the pattern in the given directory
    search_path = os.path.join(args.data_dir + "snap", args.pattern)
    # Sort the files numerically to ensure the time evolution is correct
    file_list = sorted(glob.glob(search_path))

    log_path = args.data_dir + "output.log"
    if not file_list:
        print(
            f"Error: No files found matching pattern '{args.pattern}' in directory '{args.data_dir}'."
        )
        print("Please check the directory path and file pattern.")
    else:
        main(file_list, log_path)
