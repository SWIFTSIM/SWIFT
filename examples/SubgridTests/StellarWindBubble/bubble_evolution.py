import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import unyt as u


def plot_peak_quantities_vs_time(log_file_path, args):
    """
    Reads radial_peak_positions.txt and plots simulation data vs. 
    two sets of theoretical models (Book vs. Article).
    """
    if not os.path.exists(log_file_path):
        print(f"Error: Log file not found at path: {log_file_path}")
        return

    print(f"Reading data from {log_file_path}...")

    # 1. Load Data
    try:
        data = np.genfromtxt(log_file_path, skip_header=1, dtype=float)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # --- UNPACKING THE COLUMNS ---
    # (Mapping based on common snowplow comparison structures)
    time = data[:, 0]
    R_rho_max = data[:, 1]
    rho_min = data[:, 2]
    t_max = data[:, 3]

    # Theory columns
    R_e_book, R_p_book = data[:, 6], data[:, 7]
    R_e_art, R_p_art = data[:, 8], data[:, 9]

    v_e, v_p = data[:, 10], data[:, 11]

    rho_theory = data[:, 5]
    T_theory = data[:, 4]

    # --- DATA CLEANING ---
    clean_mask = R_rho_max > 0
    t_c = time[clean_mask]
    r_c = R_rho_max[clean_mask]

    # --- RADIUS VS TIME (LINEAR & LOG-LOG) ---
    plt.rcParams.update({"font.size": 16})
    fig_r, ax = plt.subplots(1, 1, figsize=(10, 6))
    plt.subplots_adjust(hspace=0.25)

    if not args.log:

        # Linear Plot
        ax.plot(time, R_rho_max, "k.", label="Simulation (Shock)", alpha=0.4)
        ax.plot(time, R_e_book, "r-", label="Book (E-cons)")
        ax.plot(time, R_e_art, "r--", label="Article (E-cons)")
        ax.plot(time, R_p_book, "y-", label="Book (P-cons)")
        ax.plot(time, R_p_art, "y--", label="Article (P-cons)")

        ax.grid(True, alpha=0.3)
    else:

        # Log-Log Plot
        log_mask = (t_c > 0) & (r_c > 0)
        ax.loglog(t_c[log_mask], r_c[log_mask], "k.", label="Simulation")
        ax.loglog(time[time > 0], R_e_book[time > 0], "r-", label="Book (E-cons)")
        ax.loglog(time[time > 0], R_e_art[time > 0], "r--", label="Article (E-cons)")
        ax.loglog(time[time > 0], R_p_book[time > 0], "y-", label="Book (P-cons)")
        ax.loglog(time[time > 0], R_p_art[time > 0], "y--", label="Article (P-cons)")

        ax.grid(True, which="both", ls="-", alpha=0.2)

    ax.set_ylabel("Radius [pc]")
    ax.set_xlabel("Time [Myr]")
    ax.legend(loc="best", ncol=2)
    ax.set_title("Bubble Shell Radius: Book vs. Article Models vs Simulation")

    plt.savefig("Bubble_radius_comparison.png", dpi=300)

    # --- SHOCK VELOCITY CALCULATION ---
    # Define the Window Size
    WINDOW = 4 if len(time) > 4 else max(len(time) // 2, 1)

    if len(time) > 1:
        # Apply Moving Average to the Radius (R_rho_max)
        weights = np.ones(WINDOW) / WINDOW
        R_smooth = np.convolve(R_rho_max, weights, mode="same")

        # Calculate Velocity (dR/dt)
        dR_dt = np.gradient(R_smooth, time)

        # Physical Constraint & Unit Conversion
        # Velocity cannot be negative in an expanding bubble
        v_calc = np.maximum(dR_dt, 0)

        # Convert from pc/Myr to km/s
        v_shock_smooth = (v_calc * u.pc / u.Myr).to("km/s").value
    else:
        v_shock_smooth = np.zeros_like(time)

    # --- THEORY VS SIMULATION COMPARISON (Velocity, Density, Temp) ---
    fig_comp, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    plt.subplots_adjust(hspace=0.1)

    # Velocity Comparison
    axes[0].plot(time, v_shock_smooth, "purple", label="Simulation", lw=2)
    axes[0].plot(time, v_e, "r-", label="Book $v_e$")
    axes[0].plot(time, v_p, "r--", label="Book $v_p$")
    axes[0].set_ylabel("Velocity [km/s]")

    # Density Comparison
    axes[1].plot(time, rho_min, "darkred", label="Simulation")
    axes[1].plot(time, rho_theory, "k--", label="Theory")
    axes[1].set_ylabel("log($\\rho$) [H/cm$^3$]")

    # Temperature Comparison
    axes[2].plot(time, t_max, "darkorange", label="Simulation")
    axes[2].plot(time, T_theory, "darkorange", linestyle="--", label="Theory")
    axes[2].set_ylabel("log(T) [K]")
    axes[2].set_xlabel("Time [Myr]")

    for ax in axes:
        ax.legend(fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Late Phase Analysis: Book vs Article Comparison", y=0.92)
    plt.savefig("theory_vs_simulation_detailed.png", dpi=300)

    print(
        "Plots generated: 'Bubble_radius_comparison.png' and 'theory_vs_simulation_detailed.png'"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("log_file_path", help="Path to radial_peak_positions.txt")
    parser.add_argument(
        "--log", action="store_true", dest="log", help="make a log-log plot"
    )
    args = parser.parse_args()
    plot_peak_quantities_vs_time(args.log_file_path, args)
