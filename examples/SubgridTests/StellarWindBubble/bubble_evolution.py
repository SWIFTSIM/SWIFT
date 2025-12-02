import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import unyt as u
from scipy.signal import savgol_filter


def plot_peak_quantities_vs_time(log_file_path):
    """
    Reads the radial_peak_positions.txt file and plots the maximum
    values (Density, Temperature, Velocity) as a function of Time.
    """
    if not os.path.exists(log_file_path):
        print(f"Error: Log file not found at path: {log_file_path}")
        return

    print(f"Reading data from {log_file_path}...")

    # 1. Load Data
    # Use genfromtxt to handle the header line and load the data.
    # The first column (Time) is used for the x-axis.
    # The peak values are in columns 3, 5, and 7.
    try:
        data = np.genfromtxt(
            log_file_path,
            skip_header=1,
            dtype=float,
            usecols=(
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
            ),  # Columns: Time, R_RhoMax, RhoMax, R_TMax, TMax, T_R_RhoMAX, T_R_vMAX, R_vmax, VMax, V_R_RhoMax, R_e. R_p, R_b_e, R_b_p, v_e, v_p, T_e, Rho_e, time_end
        )
    except Exception as e:
        print(
            f"Error loading data from file. Ensure the file format is correct (Time, R_RhoMax, RhoMax, R_RhoMin, RhoMin, R_TMax, TMax, T_R_RhoMAX, T_R_vMAX, R_vmax, VMax, V_R_RhoMax, R_e. R_p, R_b_e, R_b_p, v_e, v_p, T_e, Rho_e, time_end): {e}"
        )
        return

    # Check if any data was loaded
    if data.ndim == 1 and data.size == 0:
        print("Error: The log file is empty or only contains the header.")
        return
    elif data.ndim == 1:
        # Handle case where there is only one line of data (single time step)
        data = data.reshape(1, -1)

    # Unpack columns
    time = data[:, 0]
    R_rho_max = data[:, 1]
    rho_max = data[:, 2]
    R_rho_min = data[:, 3]
    rho_min = data[:, 4]
    R_v_max = data[:, 5]
    v_max = data[:, 6]
    v_R_peak = data[:, 7]
    R_T_max = data[:, 8]
    T_max = data[:, 9]
    T_R_peak = data[:, 10]
    T_v_peak = data[:, 11]
    R_e = data[:, 12]
    R_p = data[:, 13]
    R_b_e = data[:, 14]
    R_b_p = data[:, 15]
    v_e = data[:, 16]
    v_p = data[:, 17]
    T_e = data[:, 18]
    rho_e = data[:, 19]
    time_end = data[0, 20]  # Assuming time_end is constant and taken from the first row

    # -------------------------------------------- Bubble's Shell Radius VS Time -----------------------

    fig3, axes3 = plt.subplots(1, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplot

    index_wind_reach_ISM = 0
    # Clean up the values when the wind have not yet reached the ISM
    indices_to_compare_for_cleaning = int(len(R_rho_max) / 10)
    for i in range(len(R_rho_max)):
        if i < len(R_rho_max) / 10:
            if R_rho_max[i] > R_rho_max[indices_to_compare_for_cleaning]:
                R_rho_max[i] = 0.0
            else:
                index_wind_reach_ISM = i
                print(f"Wind reached the ISM at indices {i}")
                break
        else:
            break

    # Plot the quantity vs. time
    axes3.plot(
        time, R_rho_max, color="k", linestyle="", marker=".", label="Simulation data"
    )

    axes3.plot(
        time,
        R_e,
        color="r",
        linestyle="-",
        linewidth=2,
        label="Snowplow energy conservative phase",
    )
    axes3.plot(
        time,
        R_p,
        color="y",
        linestyle="-",
        linewidth=2,
        label="Snowplow momentum conservative phase",
    )
    axes3.plot(
        time,
        R_b_e,
        color="b",
        linestyle="--",
        linewidth=2,
        label="Book's Snowplow energy conservative phase",
    )
    axes3.plot(
        time,
        R_b_p,
        color="g",
        linestyle="--",
        linewidth=2,
        label="Book's Snowplow momentum conservative phase",
    )

    # Formatting
    axes3.set_ylabel(f"Radius [kpc]", fontsize=12)
    axes3.grid(True, linestyle="--", alpha=0.6)
    if time_end != 0.0:
        axes3.vlines(
            time_end,
            ymin=axes3.viewLim.get_points()[0][1],
            ymax=axes3.viewLim.get_points()[1][1],
            linestyle="--",
            label="Star's death",
        )

    # Final axis labels and output
    axes3.set_xlabel("Time [Myr]", fontsize=12)
    axes3.set_title("Evolution of Bubble's Shell Radius Over Time", fontsize=16, y=0.99)
    axes3.legend()

    # 4. Save Figure
    output_filename = "Bubble_radius_vs_time.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")

    # ------------------------------------ Calculating the shockwave velocity --------------------------
    v_shockwave = np.zeros(len(R_rho_max))

    # --- Savitzky-Golay Parameters ---
    # WINDOW_LENGTH: Must be an odd integer. Choose a size large enough to cover the
    # time-scale of the noise, but small enough to resolve physical transitions.
    WINDOW_LENGTH = 35  # (e.g., covers 35 data points/snaps)

    # POLY_ORDER: Degree of the local polynomial fit (e.g., 2 or 3).
    POLY_ORDER = 2

    # 1. Calculate the derivative of R with respect to the array index (i.e., dR/di)
    V_raw_scaled = savgol_filter(
        R_rho_max,
        window_length=WINDOW_LENGTH,
        polyorder=POLY_ORDER,
        deriv=1,  # returns the derivative
    )

    # 2. Calculate the derivative of time with respect to the array index (i.e., dt/di)
    dt_di = np.gradient(time)  # This is effectively dt at each point

    # 3. Scale to get the true derivative V = (dR/di) / (dt/di) = dR/dt
    V_shockwave_smooth_values = V_raw_scaled / dt_di

    # 4. Apply units and convert
    V_shockwave_smooth = V_shockwave_smooth_values * u.kpc / u.Myr
    V_shockwave_smooth.convert_to_units(u.km / u.s)
    V_shockwave_smooth = V_shockwave_smooth.value  # Final smooth velocity array in km/s

    # ------------------------------------ Quantity VS time --------------------------------------
    fig1, axes1 = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplots

    # Define quantities for easy iteration and plotting
    quantities = [
        {
            "data": v_max,
            "ax_idx": 0,
            "title": "Maximum Velocity vs. Time",
            "ylabel": "Max $V_r$ [km/s]",
            "color": "darkblue",
        },
        {
            "data": rho_max,
            "ax_idx": 1,
            "title": "Maximum Density vs. Time",
            "ylabel": "Max log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "data": T_max,
            "ax_idx": 2,
            "title": "Maximum Temperature vs. Time",
            "ylabel": "Max log(T) [K]",
            "color": "darkorange",
        },
    ]

    # 3. Plotting Loop
    for q in quantities:
        ax = axes1[q["ax_idx"]]

        # Plot the quantity vs. time
        ax.plot(time, q["data"], color=q["color"], linestyle="-", linewidth=2)

        # Formatting
        ax.set_ylabel(q["ylabel"], fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)

    # Final axis labels and output
    axes1[-1].set_xlabel("Time [Myr]", fontsize=12)
    fig1.suptitle("Evolution of Peak Radial Quantities Over Time", fontsize=16, y=0.99)

    # 4. Save Figure
    output_filename = "peak_quantities_vs_time.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")

    # ------------------------------------ Quantity (at shockwave) VS time --------------------------------------
    fig1, axes1 = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplots

    # Define quantities for easy iteration and plotting
    quantities = [
        {
            "data": v_max,
            "ax_idx": 0,
            "title": "Velocity at shockwave vs. Time",
            "ylabel": "$V_r$ [km/s]",
            "color": "darkblue",
        },
        {
            "data": rho_max,
            "ax_idx": 1,
            "title": "Density at shockwave vs. Time",
            "ylabel": "log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "data": T_v_peak,
            "ax_idx": 2,
            "title": "Temperature at shockwave vs. Time",
            "ylabel": "log(T) [K]",
            "color": "darkorange",
        },
    ]

    # 3. Plotting Loop
    for q in quantities:
        ax = axes1[q["ax_idx"]]
        # Plot the quantity vs. time
        ax.plot(
            time,
            q["data"],
            color=q["color"],
            linestyle="-",
            linewidth=2,
            label="Simulation data",
        )
        if time_end != 0.0:
            ax.vlines(
                time_end,
                ymin=ax.viewLim.get_points()[0][1],
                ymax=ax.viewLim.get_points()[1][1],
                color="lightblue",
                linestyle="--",
                label="Star's death",
            )

        # Formatting
        ax.set_ylabel(q["ylabel"], fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend()

    # Final axis labels and output
    axes1[2].set_xlabel("Time [Myr]", fontsize=12)
    fig1.suptitle("Evolution of showckwave Over Time", fontsize=16, y=0.99)

    # 4. Save Figure
    output_filename = "shockwave_vs_time.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")

    # -------------------------------------------- Quantity VS radius -----------------------

    fig2, axes2 = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplots

    # Define quantities for easy iteration and plotting
    quantities = [
        {
            "x_axis": R_v_max,
            "data": v_max,
            "ax_idx": 0,
            "title": "Maximum Velocity vs. Radius",
            "ylabel": "Max $V_r$ [km/s]",
            "color": "darkblue",
        },
        {
            "x_axis": R_rho_max,
            "data": rho_max,
            "ax_idx": 1,
            "title": "Maximum Density vs. Radius",
            "ylabel": "Max log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "x_axis": R_T_max,
            "data": T_max,
            "ax_idx": 2,
            "title": "Maximum Temperature vs. Radius",
            "ylabel": "Max log(T) [K]",
            "color": "darkorange",
        },
    ]

    # 3. Plotting Loop
    for q in quantities:
        ax = axes2[q["ax_idx"]]

        # Plot the quantity vs. time
        ax.plot(q["x_axis"], q["data"], color=q["color"], linestyle="", marker=".")

        # Formatting
        ax.set_ylabel(q["ylabel"], fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)

    # Final axis labels and output
    axes2[-1].set_xlabel("Radius [kpc]", fontsize=12)
    fig2.suptitle(
        "Evolution of Peak Radial Quantities Over Radius", fontsize=16, y=0.99
    )

    # 4. Save Figure
    output_filename = "peak_quantities_vs_radius.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")

    # -------------------------------------------- Quantity (at shockwave) VS radius -----------------------

    fig2, axes2 = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplots

    # Define quantities for easy iteration and plotting
    quantities = [
        {
            "x_axis": R_rho_max,
            "data": v_R_peak,
            "ax_idx": 0,
            "title": "Velocity at shockwave vs. Radius",
            "ylabel": "$V_r$ [km/s]",
            "color": "darkblue",
        },
        {
            "x_axis": R_rho_max,
            "data": rho_max,
            "ax_idx": 1,
            "title": "Density at shockwave vs. Radius",
            "ylabel": "log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "x_axis": R_rho_max,
            "data": T_R_peak,
            "ax_idx": 2,
            "title": "Temperature at shockwave vs. Radius",
            "ylabel": "log(T) [K]",
            "color": "darkorange",
        },
    ]

    # 3. Plotting Loop
    for q in quantities:
        ax = axes2[q["ax_idx"]]

        # Plot the quantity vs. time
        ax.plot(q["x_axis"], q["data"], color=q["color"], linestyle="", marker=".")

        # Formatting
        ax.set_ylabel(q["ylabel"], fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)

    # Final axis labels and output
    axes2[-1].set_xlabel("Radius [kpc]", fontsize=12)
    fig2.suptitle("Evolution of Shockwave Over Radius", fontsize=16, y=0.99)

    # 4. Save Figure
    output_filename = "shockwave_vs_radius.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")

    # -------------------------------------------- Theory vs Simulation -------------------------------------------

    fig4, axes4 = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplots

    # Define quantities for easy iteration and plotting
    quantities = [
        {
            "data": V_shockwave_smooth,
            "ax_idx": 0,
            "theory": v_e,
            "title": "Shockwave Velocity vs. Time",
            "ylabel": "$V_{shock}$ [km/s]",
            "color": "purple",
        },
        {
            "data": rho_min,
            "theory": rho_e,
            "ax_idx": 1,
            "title": "Density of hot compressed wind region vs. Time",
            "ylabel": "log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "data": T_max,
            "theory": T_e,
            "ax_idx": 2,
            "title": "Temperature of hot compressed wind region vs. Time",
            "ylabel": "log(T) [K]",
            "color": "darkorange",
        },
    ]

    # 3. Plotting Loop
    for q in quantities:
        ax = axes4[q["ax_idx"]]
        # Plot the quantity vs. time
        ax.plot(
            time,
            q["data"],
            color=q["color"],
            linestyle="-",
            linewidth=2,
            label="Simulation data",
        )
        ax.plot(
            time,
            q["theory"],
            color=q["color"],
            linestyle="--",
            linewidth=2,
            label="Snowplow energy conservative phase",
        )
        if q["ax_idx"] == (0):
            ax.plot(
                time,
                v_p,
                color="purple",
                linestyle="-.",
                linewidth=2,
                label="Snowplow momentum conservative phase",
            )
        if time_end != 0.0:
            ax.vlines(
                time_end,
                ymin=ax.viewLim.get_points()[0][1],
                ymax=ax.viewLim.get_points()[1][1],
                color="lightblue",
                linestyle="--",
                label="Star's death",
            )

        # Formatting
        ax.set_ylabel(q["ylabel"], fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend()

    # Final axis labels and output
    axes4[2].set_xlabel("Time [Myr]", fontsize=12)
    fig4.suptitle(
        "Evolution of theoretical quantities VS the simulation data",
        fontsize=16,
        y=0.99,
    )

    # 4. Save Figure
    output_filename = "theory_vs_simulation.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nSuccessfully generated and saved plot to: {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot the peak quantities (Density, Temperature, Velocity) as a function of time from a log file."
    )
    # Define the command-line argument for the log file path
    parser.add_argument(
        "log_file_path",
        type=str,
        help="Path to the radial_peak_positions.txt file generated by the analysis script.",
    )

    args = parser.parse_args()
    plot_peak_quantities_vs_time(args.log_file_path)
