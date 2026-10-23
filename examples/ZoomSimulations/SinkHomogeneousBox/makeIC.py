################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#
################################################################################
"""Generate the ICs of the zoom homogeneous box (see README)."""

import argparse

import h5py
import numpy as np

# Constants in cgs
M_PROTON = 1.67262192369e-24
G_NEWTON = 6.67430e-8
MSUN = 1.988409870698051e33
KPC = 3.0856775814913673e21
MPC = 1e3 * KPC

# Internal units: 1e10 Msun, 1 kpc, 1 km/s
UNIT_MASS = 1e10 * MSUN
UNIT_LENGTH = KPC
UNIT_VELOCITY = 1e5
UNIT_TIME = UNIT_LENGTH / UNIT_VELOCITY

# SPH constants of the Wendland C2 kernel and space_stretch (src/space.h)
ETA_RESOLUTION = 1.2348
KERNEL_GAMMA = 1.936492
SPACE_STRETCH = 1.1


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--level", type=int, default=3, help="Number of gas particles 2^(3 level)"
    )
    parser.add_argument(
        "--lJ", type=float, default=0.25, help="Jeans wavelength in gas box units"
    )
    parser.add_argument(
        "--rho", type=float, default=0.1, help="Gas density in atom/cm^3"
    )
    parser.add_argument(
        "--mass", type=float, default=50.0, help="Gas particle mass in Msun"
    )
    parser.add_argument(
        "--embed",
        type=int,
        default=8,
        help="Side of the periodic box in gas box units (power of 2)",
    )
    parser.add_argument(
        "--dm-ratio",
        type=float,
        default=5.36,
        help="Dark matter to gas density ratio in the gas box",
    )
    parser.add_argument(
        "--delta",
        type=float,
        default=0.0,
        help="Density contrast rho_in / rho_out - 1 of the gas box",
    )
    parser.add_argument(
        "--transition-width",
        type=float,
        default=0.25,
        help="Width of the tanh density profile in gas box units",
    )
    parser.add_argument("--h-max", type=float, default=0.08, help="SPH:h_max in kpc")
    parser.add_argument(
        "--pad-factor",
        type=float,
        default=1.3,
        help="ZoomRegion:region_pad_factor",
    )
    parser.add_argument(
        "--max-zoom-cdim",
        type=int,
        default=32,
        help="Maximal number of zoom top-level cells per axis",
    )
    parser.add_argument(
        "--softening-ratio-background",
        type=float,
        default=0.2,
        help="Gravity:softening_ratio_background (printed softenings only)",
    )
    parser.add_argument(
        "--cosmological",
        action="store_true",
        help="Write comoving ICs for a run with --cosmology",
    )
    parser.add_argument(
        "--omega-m", type=float, default=0.3111, help="Omega_m (with --cosmological)"
    )
    parser.add_argument(
        "--hubble-param",
        type=float,
        default=0.6766,
        help="Hubble parameter h (with --cosmological)",
    )
    parser.add_argument(
        "--a-end-factor",
        type=float,
        default=2.0,
        help="Cosmology:a_end in units of a_begin (with --cosmological)",
    )
    parser.add_argument(
        "--softening",
        type=float,
        default=0.005,
        help="Physical softening of gas and PartType1 in kpc (with --cosmological)",
    )
    parser.add_argument(
        "--velocity-sigma",
        type=float,
        default=1.0,
        help="Random velocity dispersion in km/s (with --cosmological only)",
    )
    parser.add_argument("--seed", type=int, default=1, help="Random seed")
    parser.add_argument(
        "-o", dest="output", default="ICs_zoom_homogeneous_box.hdf5", help="IC file"
    )
    parser.add_argument(
        "--zoom-params-file",
        default="zoom_runtime_params.txt",
        help="File for the runtime parameters passed with -P by run.sh",
    )
    opt = parser.parse_args()

    if opt.level < 2:
        parser.error("--level must be at least 2")
    if opt.embed < 2 or opt.embed & (opt.embed - 1):
        parser.error("--embed must be a power of 2, at least 2")
    if opt.dm_ratio < 0 or opt.delta < 0 or opt.transition_width <= 0:
        parser.error("--dm-ratio and --delta must be >= 0, the width > 0")
    return opt


def gas_and_dm_box(opt, rng):
    """Paired gas and dark matter particles in [0, L)^3, in cgs."""
    n_side = 2**opt.level
    N = n_side**3
    rho_gas = opt.rho * M_PROTON
    m_gas = opt.mass * MSUN
    L = (N * m_gas / rho_gas) ** (1.0 / 3.0)
    dx = L / n_side

    # Velocity dispersion of the original HomogeneousBox example
    k_J = 2.0 * np.pi / (opt.lJ * L)
    sigma = np.sqrt(4.0 * np.pi * G_NEWTON * rho_gas) / k_J

    # Split each parent point as in space_generate_gas() (src/space.c)
    m_dm = opt.dm_ratio * m_gas
    gas_fraction = m_gas / (m_gas + m_dm)
    parents = rng.random((N, 3)) * L
    pos_gas = np.mod(parents - 0.5 * dx * (1.0 - gas_fraction), L)
    pos_dm = np.mod(parents + 0.5 * dx * gas_fraction, L)

    return dict(
        N=N,
        L=L,
        dx=dx,
        rho_gas=rho_gas,
        m_gas=m_gas,
        m_dm=m_dm,
        sigma=sigma,
        pos_gas=pos_gas,
        pos_dm=pos_dm if m_dm > 0 else np.zeros((0, 3)),
    )


def background_shells(opt, box, rng):
    """PartType2 on jittered grids in nested cube shells around the gas box."""
    L, dx = box["L"], box["dx"]
    n_side = 2**opt.level
    rho_in = box["rho_gas"] * (1.0 + opt.dm_ratio)
    rho_out = rho_in / (1.0 + opt.delta)
    width = opt.transition_width * L
    centre = 0.5 * opt.embed * L

    idx = (np.arange(n_side) + 0.5) / n_side
    gx, gy, gz = np.meshgrid(idx, idx, idx, indexing="ij")
    unit = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel())) - 0.5
    inner = np.all(np.abs(unit) < 0.25, axis=1)
    unit = unit[~inner]

    pos, mass, shell = [], [], []
    n_shells = int(round(np.log2(opt.embed)))
    for k in range(1, n_shells + 1):
        side = 2**k * L
        spacing = 2**k * dx
        cell_centres = unit * side
        r = np.max(np.abs(cell_centres), axis=1)
        rho = rho_out + (rho_in - rho_out) * (1.0 - np.tanh((r - 0.5 * L) / width))
        jitter = rng.uniform(-0.45, 0.45, size=cell_centres.shape) * spacing
        pos.append(centre + cell_centres + jitter)
        mass.append(rho * spacing**3)
        shell.append(np.full(len(r), k))

    return np.concatenate(pos), np.concatenate(mass), np.concatenate(shell)


def zoom_geometry(opt, L):
    """Background cdim and zoom depth for the given gas box side (cgs)."""
    box_size = opt.embed * L
    bkg_cdim = int(np.floor(opt.embed / opt.pad_factor))
    if bkg_cdim % 2 == 0:
        bkg_cdim -= 1
    if bkg_cdim < 1:
        raise ValueError("--embed too small for --pad-factor")
    bkg_width = box_size / bkg_cdim
    min_zoom_width = SPACE_STRETCH * KERNEL_GAMMA * opt.h_max * KPC
    depth = 0
    while (
        bkg_width / 2 ** (depth + 1) > min_zoom_width
        and 2 ** (depth + 1) <= opt.max_zoom_cdim
    ):
        depth += 1
    if depth < 1:
        raise ValueError(
            "Zoom cells narrower than 1.1 kernel_gamma h_max: lower --h-max "
            "or --pad-factor, or raise --level"
        )
    return bkg_cdim, depth, bkg_width, bkg_width / 2**depth, min_zoom_width


def write_ics(opt, filename, a, box_size, groups):
    """Write the particle groups with h5py. Lengths are divided by a."""
    n_part = [0] * 6
    for ptype, data in groups.items():
        n_part[ptype] = len(data["ParticleIDs"])

    with h5py.File(filename, "w") as f:
        grp = f.create_group("/Header")
        grp.attrs["BoxSize"] = [box_size / a] * 3
        grp.attrs["NumPart_Total"] = n_part
        grp.attrs["NumPart_Total_HighWord"] = [0] * 6
        grp.attrs["NumPart_ThisFile"] = n_part
        grp.attrs["Time"] = a if opt.cosmological else 0.0
        grp.attrs["NumFileOutputsPerSnapshot"] = 1
        grp.attrs["MassTable"] = [0.0] * 6
        grp.attrs["Flag_Entropy_ICs"] = [0] * 6
        grp.attrs["Dimension"] = 3

        grp = f.create_group("/Units")
        grp.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH
        grp.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS
        grp.attrs["Unit time in cgs (U_t)"] = UNIT_TIME
        grp.attrs["Unit current in cgs (U_I)"] = 1.0
        grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

        for ptype, data in groups.items():
            grp = f.create_group(f"/PartType{ptype}")
            for name, (values, dtype) in data_fields(data, a).items():
                grp.create_dataset(name, data=values, dtype=dtype)


def data_fields(data, a):
    fields = {
        "Coordinates": (data["Coordinates"] / a, "d"),
        "Velocities": (data["Velocities"], "f"),
        "Masses": (data["Masses"], "f"),
        "ParticleIDs": (data["ParticleIDs"], "L"),
    }
    if "SmoothingLength" in data:
        fields["SmoothingLength"] = (data["SmoothingLength"] / a, "f")
        fields["InternalEnergy"] = (data["InternalEnergy"], "f")
        fields["Densities"] = (data["Densities"] * a**3, "f")
    return fields


def main():
    opt = parse_options()
    rng = np.random.default_rng(opt.seed)

    box = gas_and_dm_box(opt, rng)
    L, N = box["L"], box["N"]
    box_size = opt.embed * L
    offset = 0.5 * (opt.embed - 1) * L
    centre = 0.5 * box_size

    bkg_cdim, depth, bkg_width, zoom_width, min_zoom_width = zoom_geometry(opt, L)
    mesh = max(32, 2 ** (opt.level + 1))
    pos_bg, m_bg, shell_bg = background_shells(opt, box, rng)

    # Checks
    if np.any(np.max(np.abs(pos_bg - centre), axis=1) < 0.5 * L):
        raise RuntimeError("PartType2 particle inside the gas box")
    if np.any(pos_bg < 0) or np.any(pos_bg >= box_size):
        raise RuntimeError("PartType2 particle outside the periodic box")

    h = min(ETA_RESOLUTION * box["dx"], opt.h_max * KPC)
    N_dm = len(box["pos_dm"])
    N_bg = len(m_bg)

    # Cosmological ICs: comoving mean matter density equal to Omega_m rho_crit_0
    M_total = N * box["m_gas"] + N_dm * box["m_dm"] + m_bg.sum()
    rho_mean = M_total / box_size**3
    a = 1.0
    if opt.cosmological:
        H0 = 100.0 * opt.hubble_param * 1e5 / MPC
        rho_crit_0 = 3.0 * H0**2 / (8.0 * np.pi * G_NEWTON)
        a = (opt.omega_m * rho_crit_0 / rho_mean) ** (1.0 / 3.0)

    to_mass = 1.0 / UNIT_MASS
    to_length = 1.0 / UNIT_LENGTH
    groups = {
        0: {
            "Coordinates": (box["pos_gas"] + offset) * to_length,
            "Masses": np.full(N, box["m_gas"] * to_mass),
            "ParticleIDs": np.arange(1, N + 1),
            "SmoothingLength": np.full(N, h * to_length),
            "InternalEnergy": np.full(N, box["sigma"] ** 2 / UNIT_VELOCITY**2),
            "Densities": np.full(N, box["rho_gas"] * to_mass / to_length**3),
        },
        2: {
            "Coordinates": pos_bg * to_length,
            "Masses": m_bg * to_mass,
            "ParticleIDs": np.arange(N + N_dm + 1, N + N_dm + N_bg + 1),
        },
    }
    if N_dm > 0:
        groups[1] = {
            "Coordinates": (box["pos_dm"] + offset) * to_length,
            "Masses": np.full(N_dm, box["m_dm"] * to_mass),
            "ParticleIDs": np.arange(N + 1, N + N_dm + 1),
        }
    groups = dict(sorted(groups.items()))

    # SWIFT stops on a zero RMS velocity in cosmological runs
    # (engine_recompute_displacement_constraint), so these ICs get a small
    # random velocity
    sigma_v = opt.velocity_sigma if opt.cosmological else 0.0
    for data in groups.values():
        shape = data["Coordinates"].shape
        data["Velocities"] = (
            rng.normal(0.0, sigma_v, shape) if sigma_v > 0 else np.zeros(shape)
        )
    write_ics(opt, opt.output, a, box_size * to_length, groups)

    runtime = [
        f"ZoomRegion:bkg_top_level_cells:{bkg_cdim}",
        f"ZoomRegion:zoom_top_level_depth:{depth}",
        f"ZoomRegion:region_pad_factor:{opt.pad_factor}",
        f"Gravity:mesh_side_length:{mesh}",
    ]
    if opt.cosmological:
        omega_b = opt.omega_m / (1.0 + opt.dm_ratio)
        runtime += [
            f"Cosmology:Omega_cdm:{opt.omega_m - omega_b:.8f}",
            f"Cosmology:Omega_b:{omega_b:.8f}",
            f"Cosmology:Omega_lambda:{1.0 - opt.omega_m:.8f}",
            f"Cosmology:h:{opt.hubble_param}",
            f"Cosmology:a_begin:{a:.8e}",
            f"Cosmology:a_end:{opt.a_end_factor * a:.8e}",
            f"Gravity:comoving_DM_softening:{opt.softening / a:.8e}",
            f"Gravity:comoving_baryon_softening:{opt.softening / a:.8e}",
            f"SPH:h_max:{opt.h_max / a:.8e}",
            f"Snapshots:scale_factor_first:{a:.8e}",
            "Snapshots:delta_time:1.02",
            f"Statistics:scale_factor_first:{a:.8e}",
            "Statistics:delta_time:1.01",
        ]
    with open(opt.zoom_params_file, "w") as f:
        f.write(" ".join(f"-P {r}" for r in runtime) + "\n")

    # Report
    pc = KPC * 1e-3
    print(
        f"{opt.output} written ({'comoving, a = %.6e' % a if opt.cosmological else 'physical'})"
    )
    print(
        f"Gas box side L            : {L / pc:.2f} pc, periodic box {box_size / pc:.2f} pc"
    )
    print(
        f"Gas box                   : [{offset / pc:.2f}, {(offset + L) / pc:.2f}) pc on each axis"
    )
    print(f"Mean gas spacing          : {box['dx'] / pc:.3f} pc, h = {h / pc:.3f} pc")
    print(f"Velocity dispersion       : {box['sigma'] / 1e5:.4f} km/s")
    print(f"PartType0 (gas)           : {N} x {box['m_gas'] / MSUN:.4g} Msun")
    print(f"PartType1 (dark matter)   : {N_dm} x {box['m_dm'] / MSUN:.4g} Msun")
    for k in np.unique(shell_bg):
        sel = shell_bg == k
        m_k = m_bg[sel]
        eps = opt.softening_ratio_background * (m_k / rho_mean) ** (1.0 / 3.0)
        print(
            f"PartType2 shell {k}         : {sel.sum()} particles, "
            f"{m_k.min() / MSUN:.4g} to {m_k.max() / MSUN:.4g} Msun, "
            f"Plummer softening {eps.min() / pc:.3g} to {eps.max() / pc:.3g} pc "
            f"(ratio {opt.softening_ratio_background})"
        )
    M_in = N * box["m_gas"] + N_dm * box["m_dm"]
    M_out = m_bg.sum()
    V_in = L**3
    V_out = box_size**3 - V_in
    print(
        f"Mass inside the gas box   : {M_in / MSUN:.6e} Msun, density {M_in / V_in / M_PROTON:.6f} m_p/cm^3"
    )
    print(
        f"Mass outside the gas box  : {M_out / MSUN:.6e} Msun, density {M_out / V_out / M_PROTON:.6f} m_p/cm^3"
    )
    print(
        f"Total mass                : {M_total / MSUN:.6e} Msun, mean density {rho_mean / M_PROTON:.6f} m_p/cm^3"
    )
    print("PartType2 inside gas box  : 0 (checked)")
    print(
        f"Zoom grid                 : bkg_top_level_cells {bkg_cdim} "
        f"(width {bkg_width / pc:.1f} pc), zoom_top_level_depth {depth} "
        f"(zoom cell {zoom_width / pc:.1f} pc > 1.1 kernel_gamma h_max = "
        f"{min_zoom_width / pc:.1f} pc), mesh_side_length {mesh}, "
        f"written to {opt.zoom_params_file}"
    )
    if opt.cosmological:
        print(f"Cosmological runtime parameters (also in {opt.zoom_params_file}):")
        for r in runtime[3:]:
            print(f"  {r}")


if __name__ == "__main__":
    main()
