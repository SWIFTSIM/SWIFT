###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

import h5py
import numpy as np


def seed_flaws(num_part, volume, m_flaw, k_flaw, hardcoded_max_flaws, seed, chunk_size):
    """Assign Weibull-distributed activation thresholds to particles, adding flaws until each particle has at least one flaw."""
    rng = np.random.default_rng(seed)
    activation_thresholds = np.zeros((num_part, hardcoded_max_flaws))

    # Use chunks to avoid memory issues because of storing (num_draws, num_part) arrays
    chunks = []
    seen = np.zeros(num_part, dtype=bool)
    num_seen = 0

    # Loop over chunks until all particles have at least one flaw
    while num_seen < num_part:
        # Draw chunk worth of particle IDs
        chunk = rng.integers(0, num_part, size=chunk_size)

        # Find all the unique particle IDs that haven't been drawn before
        unique_in_chunk = np.unique(chunk)
        newly_seen = (~seen[unique_in_chunk]).sum()

        # If all particles have been drawn, finish up
        if num_seen + newly_seen >= num_part:
            # Mark which draws in this chunk are the first occurrence of each particle
            _, first_occurrences = np.unique(chunk, return_index=True)
            new_in_chunk = np.zeros(len(chunk), dtype=bool)
            new_in_chunk[first_occurrences] = True
            new_in_chunk &= ~seen[
                chunk
            ]  # exclude particles already seen before this chunk

            # Find the cutoff point where the final missing ID gets drawn for the first time
            cumulative_new = np.cumsum(new_in_chunk)
            cutoff = np.argmax(cumulative_new == (num_part - num_seen))

            # Add this chunk, up to the cutoff, to the list keeping track of IDs
            chunks.append(chunk[: cutoff + 1])
            break

        # Update quantities keeping track of the IDs that we've seen and move onto next chunk
        seen[unique_in_chunk] = True
        num_seen += newly_seen
        chunks.append(chunk)

    # Array of drawn indices
    sample_indices = np.concatenate(chunks)

    # Threshold for the nth flaw
    flaw_numbers = np.arange(1, len(sample_indices) + 1)
    thresholds = (flaw_numbers / (k_flaw * volume)) ** (1.0 / m_flaw)

    # For each draw, how many flaws did that particle already have
    order = np.argsort(sample_indices, kind="stable")
    prior_flaw_counts = np.empty(len(sample_indices), dtype=np.int32)
    prior_flaw_counts[order] = np.concatenate(
        [np.arange(count) for count in np.bincount(sample_indices, minlength=num_part)]
    )

    # Fill in the thresholds and number of flaws for each particle in the right orther
    activation_thresholds[sample_indices, prior_flaw_counts] = thresholds
    N_flaws = np.bincount(sample_indices, minlength=num_part)

    return N_flaws, activation_thresholds


def generate_target(
    L_target,
    boxsize,
    vol,
    numPart_grid,
    target_radius,
    rho_target,
    mat_id_target,
    u_target,
    m_flaw_target,
    k_flaw_target,
    hardcoded_max_flaws,
    seed,
    chunk_size,
):
    # Cubic grid
    coords = np.linspace(0, boxsize, L_target, endpoint=False)
    gx, gy, gz = np.meshgrid(coords, coords, coords, indexing="ij")
    A2_pos_grid = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    # Carve out sphere
    centre = boxsize / 2
    r2 = (
        (A2_pos_grid[:, 0] - centre) ** 2
        + (A2_pos_grid[:, 1] - centre) ** 2
        + (A2_pos_grid[:, 2] - centre) ** 2
    )
    A2_pos_target = A2_pos_grid[r2 <= target_radius**2]
    num_target = len(A2_pos_target)
    print(f"Target particles: {num_target}")

    # Fill qunatities
    A2_vel_target = np.zeros((num_target, 3))
    A1_h_target = np.full(num_target, (boxsize / L_target) * 1.2348)
    A1_m_target = np.full(num_target, vol * rho_target / numPart_grid)
    A1_mat_target = np.full(num_target, mat_id_target)
    A1_rho_target = np.full(num_target, rho_target)
    A1_u_target = np.full(num_target, u_target)

    # Add flaws
    flaw_volume = num_target * vol / numPart_grid
    A1_nflaws_target, A2_thresh_target = seed_flaws(
        num_target,
        flaw_volume,
        m_flaw_target,
        k_flaw_target,
        hardcoded_max_flaws,
        seed,
        chunk_size,
    )
    print(f"Max flaws of a single particle: {A1_nflaws_target.max()}")

    return (
        A2_pos_target,
        A2_vel_target,
        A1_h_target,
        A1_m_target,
        A1_mat_target,
        A1_rho_target,
        A1_u_target,
        A1_nflaws_target,
        A2_thresh_target,
    )


def generate_impactor(
    L_target,
    boxsize,
    vol,
    numPart_grid,
    impactor_radius,
    rho_impactor,
    mat_id_impactor,
    u_impactor,
    vx_impactor,
    impact_angle,
    target_radius,
    hardcoded_max_flaws,
):
    # Reuse the same grid resolution as the target
    coords = np.linspace(0, boxsize, L_target, endpoint=False)
    gx, gy, gz = np.meshgrid(coords, coords, coords, indexing="ij")
    A2_pos_grid = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    # Carve out sphere
    centre = boxsize / 2
    r2 = (
        (A2_pos_grid[:, 0] - centre) ** 2
        + (A2_pos_grid[:, 1] - centre) ** 2
        + (A2_pos_grid[:, 2] - centre) ** 2
    )
    A2_pos_impactor = A2_pos_grid[r2 <= impactor_radius**2].copy()
    num_impactor = len(A2_pos_impactor)
    print(f"Impactor particles: {num_impactor}")

    # Fill qunatities
    A2_vel_impactor = np.zeros((num_impactor, 3))
    A1_h_impactor = np.full(num_impactor, (boxsize / L_target) * 1.2348)
    A1_m_impactor = np.full(num_impactor, vol * rho_impactor / numPart_grid)
    A1_mat_impactor = np.full(num_impactor, mat_id_impactor)
    A1_rho_impactor = np.full(num_impactor, rho_impactor)
    A1_u_impactor = np.full(num_impactor, u_impactor)

    # Impactor velocity
    A2_vel_impactor[:, 0] = vx_impactor

    # Add offset position
    initial_offset = 0.2 * boxsize
    impact_parameter = (target_radius + impactor_radius) * np.sin(
        np.deg2rad(impact_angle)
    )
    A2_pos_impactor[:, 0] -= initial_offset
    A2_pos_impactor[:, 1] += impact_parameter

    # Impactor has no flaws
    A1_nflaws_impactor = np.zeros(num_impactor)
    A2_thresh_impactor = np.zeros((num_impactor, hardcoded_max_flaws))

    return (
        A2_pos_impactor,
        A2_vel_impactor,
        A1_h_impactor,
        A1_m_impactor,
        A1_mat_impactor,
        A1_rho_impactor,
        A1_u_impactor,
        A1_nflaws_impactor,
        A2_thresh_impactor,
    )


def save(fileOutputName, boxsize, target, impactor):

    # Combine target and impactor arrays
    pos = np.vstack([target[0], impactor[0]])
    vel = np.vstack([target[1], impactor[1]])
    h = np.concatenate([target[2], impactor[2]])
    m = np.concatenate([target[3], impactor[3]])
    mat = np.concatenate([target[4], impactor[4]])
    rho = np.concatenate([target[5], impactor[5]])
    u = np.concatenate([target[6], impactor[6]])
    nflaws = np.concatenate([target[7], impactor[7]])
    thresholds = np.vstack([target[8], impactor[8]])
    ids = np.arange(1, len(h) + 1)

    # Save
    numPart = len(h)
    with h5py.File(fileOutputName, "w") as f:
        hdr = f.create_group("Header")
        hdr.attrs["BoxSize"] = [boxsize, boxsize, boxsize]
        hdr.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
        hdr.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        hdr.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
        hdr.attrs["Time"] = 0.0
        hdr.attrs["NumFilesPerSnapshot"] = 1
        hdr.attrs["MassTable"] = [0.0] * 6
        hdr.attrs["Flag_Entropy_ICs"] = [0] * 6
        hdr.attrs["Dimension"] = 3

        units = f.create_group("Units")
        units.attrs["Unit length in cgs (U_L)"] = 100.0
        units.attrs["Unit mass in cgs (U_M)"] = 1000.0
        units.attrs["Unit time in cgs (U_t)"] = 1.0
        units.attrs["Unit current in cgs (U_I)"] = 1.0
        units.attrs["Unit temperature in cgs (U_T)"] = 1.0

        grp = f.create_group("PartType0")
        grp.create_dataset("Coordinates", data=pos, dtype="d")
        grp.create_dataset("Velocities", data=vel, dtype="f")
        grp.create_dataset("Masses", data=m.reshape(-1, 1), dtype="f")
        grp.create_dataset("Density", data=rho.reshape(-1, 1), dtype="f")
        grp.create_dataset("SmoothingLength", data=h.reshape(-1, 1), dtype="f")
        grp.create_dataset("InternalEnergy", data=u.reshape(-1, 1), dtype="f")
        grp.create_dataset("ParticleIDs", data=ids.reshape(-1, 1), dtype="L")
        grp.create_dataset("MaterialIDs", data=mat.reshape(-1, 1), dtype="i")
        grp.create_dataset("NumFlaws", data=nflaws.reshape(-1, 1), dtype="i")
        grp.create_dataset("ActivationThresholds", data=thresholds, dtype="f")


if __name__ == "__main__":
    # Parameters  (Nakamura & Fujiwara 1991, Benz & Asphaug 1994)

    # Box and resolution params
    L_target = 300  # Number of particles on side of grid
    boxsize = 20 / 100
    vol = boxsize**3
    numPart_grid = L_target * L_target * L_target

    # Target
    target_radius = 3 / 100  # cm
    rho_target = 2700
    mat_id_target = 1000
    u_target = 0.074074  # To give pressure = 400 Pa

    # Flaw parameters from B&A 1994
    m_flaw_target = 8
    k_flaw_target = 5e28 * 100**3  # cm^-3 => m^-3
    seed = 12345
    chunk_size = 100000  # to avoid memory issues with large arrays
    hardcoded_max_flaws = 40

    # Impactor
    rho_impactor = 1180
    mass_impactor = 0.2 / 1000  # g => kg
    impactor_radius = (mass_impactor / (4 / 3 * np.pi * rho_impactor)) ** (1 / 3)
    mat_id_impactor = 1001
    u_impactor = 0.074074  # To give pressure = 400 Pa
    vx_impactor = 3.2e5 / 100  # cm => m
    impact_angle = 30

    fileOutputName = "impact.hdf5"
    # -------------------------------------

    target = generate_target(
        L_target,
        boxsize,
        vol,
        numPart_grid,
        target_radius,
        rho_target,
        mat_id_target,
        u_target,
        m_flaw_target,
        k_flaw_target,
        hardcoded_max_flaws,
        seed,
        chunk_size,
    )
    impactor = generate_impactor(
        L_target,
        boxsize,
        vol,
        numPart_grid,
        impactor_radius,
        rho_impactor,
        mat_id_impactor,
        u_impactor,
        vx_impactor,
        impact_angle,
        target_radius,
        hardcoded_max_flaws,
    )

    save(fileOutputName, boxsize, target, impactor)
