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
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# --- GLOBAL SETTINGS ---
NUM_BINS = 500  
LAST_R_PEAK = 0.0 

def get_gas_mu(data):
    cooling = data.metadata.subgrid_scheme["Cooling Model"]
    gas = data.gas
    rho = gas.densities.to(unyt.g / unyt.cm ** 3)
    mH_in_g = mh.to(unyt.g)
    if cooling in [b"Grackle3", b"Grackle2", b"Grackle1"]:
        nHI, nHII = gas.hi * rho / mH_in_g, gas.hii * rho / mH_in_g
        nHeI, nHeII, nHeIII = gas.he_i * rho / (4*mh), gas.he_ii * rho / (4*mh), gas.he_iii * rho / (4*mh)
        nel = nHII + nHeII + 2 * nHeIII
        num = (nHI + nHII) + (nHeI + nHeII + nHeIII) * 4
        den = (nHI + nHII + nHeI + nHeII + nHeIII + nel)
        return num / den
    else:
        gamma = data.metadata.gas_gamma[0]
        u = data.gas.internal_energies.to_physical().to(unyt.erg / unyt.g)
        H_frac = float(data.metadata.parameters.get("GrackleCooling:HydrogenFractionByMass", 1.0))
        T_mu = (gamma - 1.0) * u.value * mH_in_g.value / k_B.value
        mu = np.where(T_mu > 1.1e4 / (4.0 / (8.0 - 5.0 * (1.0 - H_frac))), 4.0 / (8.0 - 5.0 * (1.0 - H_frac)), 4.0 / (1.0 + 3.0 * H_frac))
        return mu

def get_gas_temperatures(data):
    mH_in_g = mh.to(unyt.g)
    u = data.gas.internal_energies.to_physical().to(unyt.erg / unyt.g)
    gamma = data.metadata.gas_gamma[0]
    mu = get_gas_mu(data)
    T = mu * (gamma - 1.0) * u * mH_in_g / k_B
    return T.to(unyt.K)

def get_data(filename):
    try:
        data = sw.load(filename)
    except: return None
    data.gas.coordinates.convert_to_units("pc")
    data.gas.densities.convert_to_units("g/cm**3")
    data.gas.temperatures = get_gas_temperatures(data)
    data.gas.velocities.convert_to_units("km/s")
    data.metadata.time.convert_to_units("Myr")
    
    center = data.stars.coordinates[0] if len(data.stars.coordinates) > 0 else unyt.unyt_array([0,0,0], "pc")
    coords = data.gas.coordinates - center
    r = np.sqrt(np.sum(coords**2, axis=1))
    v_r = np.sum(coords * data.gas.velocities, axis=1) / np.where(r > 0, r, 1.0 * unyt.pc)
    
    log_rho = np.log10((data.gas.densities / mh).to("1/cm**3").value)
    log_T = np.log10(data.gas.temperatures.value)
    mu_val = np.mean(get_gas_mu(data))
    
    return r, v_r, log_rho, log_T, mu_val, data.metadata.time

def extract_feedback(log_path: str):
    """Parses output.log for Energy and Mass rates."""
    e_rate, m_rate = 0.0, 0.0
    if not os.path.exists(log_path): 
        print(f"Warning: {log_path} not found. Using zero for theory.")
        return 0.0, 0.0
    fp = r"([+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?)"
    with open(log_path, "r") as f:
        for line in f:
            me = re.search(r"Energy\[erg/yr\]=" + fp, line)
            mec = re.search(r"Energy_per_progenitor_mass\[erg/yr/Msol\]=" + fp, line)
            mm = re.search(r"Mass_ejected\[Msol/yr\]=" + fp, line)
            mmc = re.search(r"Mass_ejected_per_progenitor_mass\[Msol/yr/Msol\]=" + fp, line)
            massc = re.search(r"init_mass\[M_odot\]=" + fp,line)
            if me: e_rate = float(me.group(1))
            if mec and massc: e_rate = float(mec.group(1)) * float(massc.group(1))
            if mm: m_rate = float(mm.group(1))
            if mmc and massc: m_rate = float(mmc.group(1)) * float(massc.group(1))
            if e_rate > 0.0 and m_rate > 0.0:
                break
    return e_rate, m_rate

# Find the radial position of the peak at a certain time, knowing the previous peak position
def find_shock_binned(r, rho, time_val):
    global LAST_R_PEAK
    if time_val <= 0:
        LAST_R_PEAK = 0.0
        return 0.0
    
    r_max = np.max(r.value)
    r_bins = np.linspace(0, r_max, NUM_BINS)
    bin_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    indices = np.digitize(r.value, r_bins) - 1
    
    # Use binning to get the spherically symmetric value 
    rho_binned = np.zeros(NUM_BINS-1)
    for i in range(NUM_BINS-1):
        mask = (indices == i)
        if np.any(mask): 
            rho_binned[i] = np.mean(rho[mask])
        else: 
            rho_binned[i] = rho_binned[i-1] if i > 0 else np.min(rho)
    
    # Smooth the value by fitting gaussian
    rho_sm = gaussian_filter1d(rho_binned, sigma=2)
    
    # Limit search to previous peak + 15% box size to prevent jumping on another peak
    search_limit = LAST_R_PEAK + (0.15 * r_max) if LAST_R_PEAK > 0 else 0.2 * r_max
    valid_mask = bin_centers < search_limit
    search_rho = np.where(valid_mask, rho_sm, np.min(rho_sm))

    peaks, _ = find_peaks(search_rho, prominence=0.01)
    
    # If a peak is found, update LAST_R_PEAK and return it. Otherwise, return the last known peak position.
    if len(peaks) > 0:
        found_r = bin_centers[peaks[-1]]
        if LAST_R_PEAK > 0 and found_r > LAST_R_PEAK * 1.5:
             return LAST_R_PEAK 
        LAST_R_PEAK = found_r
        return found_r
    
    return LAST_R_PEAK

# Calculate theoretical values, compare with actual data and write all the result in a txt file
def write_line(f, r, log_rho, log_T, mu, time, e_rate, v_wind, n0):
    r_peak = find_shock_binned(r, log_rho, time.value)
    
    rho_min = np.min(log_rho) if len(log_rho) > 0 else np.log10(n0)

    t_max = np.max(log_T) if len(log_T) > 0 else 0.0
    
    # Theory Scaling (see "Henny J. G. L. M. Lamers, Joseph P. Cassinelli - Introduction to Stellar Winds-Cambridge University Press (1999)")
    L36 = (e_rate.to("erg/s").value) / 1e36
    t6 = time.to("yr").value / 1e6
    cm_to_pc = 1 * unyt.cm
    cm_to_pc.convert_to_units("pc")
    if t6 > 0 and L36 > 0:
        # Theoretical formula from the book
        R_e = 28 * (L36 / mu / n0)**(1/5) * t6**(3/5)                               # Radius vs time for energy-conserving phase
        R_p = 16 * (L36 / (v_wind.value/1000) / n0)**(1/4) * (t6)**(1/2)            # Radius vs time for momentum-conserving phase
        v_e = 17 * (L36 / mu / n0)**(1/5) * t6**(-2/5)                              # Velocity vs time for energy-conserving phase
        v_p = 7.9 * (L36 / (v_wind.value/1000) / n0)**(1/4) * t6**(-1/2)            # Velocity vs time for momentum-conserving phase
        T = 1.6e6 * L36**(8/35) * n0**(2/35) * t6**(-6/35)                          # Temperature vs time for energy-conserving phase
        n = 0.01 * L36**(6/35) * n0**(19/35) * t6**(-22/35) #+ n0                   # Density vs time for energy-conserving phase

        # Theoretical formula from the article (see "lahén et al. (2023) - Formation of star clusters and enrichment by massive stars in simulations of low-metallicity galaxies with a fully sampled initial stellar mass function")
        R_e_art = 0.88 * (e_rate.to("erg/s").value / (mu * n0 * mh.to("g").value))**(1/5) * time.to("s").value**(3/5) * cm_to_pc.value
        R_p_art = (3 / (2 * np.pi)**(1/4) * (2 * e_rate.to("erg/s").value / (v_wind.to("cm/s").value)) / (n0 * mu * mh.to("g").value))**(1/4) * time.to("s").value**(1/2) * cm_to_pc.value
    else: 
        R_e = R_p = v_e = v_p = R_e_art = R_p_art = 0.0
        n = n0
        T = 1.0e4
    if t6 < 0 or L36 < 0:
        print(f"Warning: Non-physical values for theory at time={time.value:.3f} Myr, Energy={e_rate:.2e}. Setting theory to zero.")


    cols = [time.value, r_peak, rho_min, 
            t_max, np.log10(T), np.log10(n), 
            R_e, R_p, R_e_art, 
            R_p_art, v_e, v_p]
    
    f.write("\t".join([f"{x:.6f}" for x in cols]) + "\n")
    return r_peak

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir")
    args = parser.parse_args()
    
    files = sorted(glob.glob(os.path.join(args.data_dir, "snap", "snapshot_*.hdf5")))
    log_path = os.path.join(args.data_dir, "output.log")
    
    # Real extraction
    e_val, m_val = extract_feedback(log_path)
    e_rate = e_val * unyt.erg / unyt.yr
    m_loss = m_val * unyt.Msun / unyt.yr
    v_wind = np.sqrt(2 * e_rate / m_loss).to("km/s") if m_val > 0 else 0 * unyt.km/unyt.s
    
    print(f"Extracted Feedback: Energy={e_rate:.2e}, MassLoss={m_loss:.2e}, WindVel={v_wind:.2f}")

    d0 = sw.load(files[0])
    n0 = np.mean((d0.gas.densities / mh).to("1/cm**3").value)

    # The loop where all required data are loaded
    data_list = []
    print("Pre-loading snapshots...")
    for f in tqdm(files):
        d = get_data(f)
        if d: data_list.append(d)

    # Preparing the gif 
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    lines = [ax.plot([], [], 'k.', markersize=0.5, alpha=0.3)[0] for ax in axes]
    vlines = [ax.axvline(0, color='red', linestyle='--', alpha=0.8) for ax in axes]
    
    mid = len(data_list)//2
    axes[0].set_ylim(data_list[mid][1].min().value * 1.2, data_list[mid][1].max().value * 1.2)
    axes[1].set_ylim(np.min(data_list[mid][2])-0.5, np.max(data_list[mid][2])+1.0)
    axes[2].set_ylim(0,7)
    axes[2].set_xlim(0, np.max([d[0].max().value for d in data_list]))

    axes[2].set_xlabel("Radius [pc]")
    axes[0].set_ylabel("Velocity [km/s]")
    axes[1].set_ylabel("log(Density) [H/cm$^3$]")
    axes[2].set_ylabel("log(Temperature) [K]")

    with open("radial_peak_positions.txt", "w") as f_log:
        f_log.write("# Time R_peak Rho_max R_min Rho_min R_Vmax Vmax V_at_r R_Tmax Tmax ... \n")

        # method to update the gif's frames and write the txt file in the same time
        def update(frame):
            r, v_r, log_rho, log_T, mu, time = data_list[frame]
            r_peak = write_line(f_log, r, log_rho, log_T, mu, time, e_rate, v_wind, n0)
            lines[0].set_data(r.value, v_r.value)
            lines[1].set_data(r.value, log_rho)
            lines[2].set_data(r.value, log_T)
            for vl in vlines: vl.set_xdata([r_peak, r_peak])
            fig.suptitle(f"Time: {time.value:.3f} Myr | R_shock: {r_peak:.3f} pc")
            return lines + vlines

        # Animate the gif using pillow by default
        ani = FuncAnimation(fig, update, frames=len(data_list), blit=True)
        ani.save("bubble_evolution.gif", writer="pillow", fps=10)

    print("Job done")

if __name__ == "__main__":
    main()