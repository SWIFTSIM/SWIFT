import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


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
            ),  # Columns: Time, R_RhoMax, RhoMax, R_TMax, TMax, T_R_RhoMAX, T_R_vMAX, R_vmax, VMax, V_R_RhoMax, R_e. R_p, R_b_e, R_b_p, v_e, v_p, T_e, Rho_e, time_end
        )
    except Exception as e:
        print(
            f"Error loading data from file. Ensure the file format is correct (Time, R_RhoMax, RhoMax, R_TMax, TMax, T_R_RhoMAX, T_R_vMAX, R_vmax, VMax, V_R_RhoMax, R_e. R_p, R_b_e, R_b_p, v_e, v_p, T_e, Rho_e, time_end): {e}"
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
    R_v_max = data[:, 3]
    v_max = data[:, 4]
    v_R_peak = data[:, 5]
    R_T_max = data[:, 6]
    T_max = data[:, 7]
    T_R_peak = data[:, 8]
    T_v_peak = data[:, 9]
    R_e = data[:, 10]
    R_p = data[:, 11]
    R_b_e = data[:, 12]
    R_b_p = data[:, 13]
    v_e = data[:, 14]
    v_p = data[:, 15]
    T_e = data[:, 16]
    rho_e = data[:, 17]
    time_end = data[0, 18]  # Assuming time_end is constant and taken from the first row

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
            "theory": v_e,
            "ax_idx": 0,
            "title": "Velocity at shockwave vs. Time",
            "ylabel": "$V_r$ [km/s]",
            "color": "darkblue",
        },
        {
            "data": rho_max,
            "theory": rho_e,
            "ax_idx": 1,
            "title": "Density at shockwave vs. Time",
            "ylabel": "log($\\rho$) [H/cm$^3$]",
            "color": "darkred",
        },
        {
            "data": T_v_peak,
            "theory": T_e,
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
        ax.plot(
            time,
            q["theory"],
            color=q["color"],
            linestyle="--",
            linewidth=2,
            label="Snowplow energy conservative phase",
        )
        if q["ax_idx"] == 0:
            ax.plot(
                time,
                v_p,
                color="blue",
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
    axes1[-1].set_xlabel("Time [Myr]", fontsize=12)
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

    # -------------------------------------------- Bubble's Shell Radius VS Time -----------------------

    fig3, axes3 = plt.subplots(1, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05)  # Reduce space between subplot

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

    # ax.set_title(q["title"], fontsize=14)

    # Final axis labels and output
    axes3.set_xlabel("Time [Myr]", fontsize=12)
    axes3.set_title("Evolution of Bubble's Shell Radius Over Time", fontsize=16, y=0.99)

    axes3.legend()

    # 4. Save Figure
    output_filename = "Bubble_radius_vs_time.png"
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
