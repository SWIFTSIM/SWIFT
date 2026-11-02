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
import sys

import h5py
import numpy as np

# Constants in cgs
M_PROTON = 1.67262192369e-24
G_NEWTON = 6.67430e-8
MSUN = 1.988409870698051e33
KPC = 3.0856775814913673e21
MPC = 1e3 * KPC
PC = 1e-3 * KPC

# Internal units: 1e10 Msun, 1 kpc, 1 km/s
UNIT_MASS = 1e10 * MSUN
UNIT_LENGTH = KPC
UNIT_VELOCITY = 1e5
UNIT_TIME = UNIT_LENGTH / UNIT_VELOCITY

# SPH constants of the Wendland C2 kernel and space_stretch (src/space.h)
ETA_RESOLUTION = 1.2348
KERNEL_GAMMA = 1.936492
SPACE_STRETCH = 1.1
GAS_GAMMA = 5.0 / 3.0

# Gravity defaults of SWIFT (src/gravity_properties.h)
A_SMOOTH = 1.25
R_CUT_MAX = 4.5

# Limits of zoom_region_init() (src/zoom_region/zoom_init.c)
ZOOM_REGIONS_ERROR = 64
ZOOM_REGIONS_WARNING = 16

# Smallest zoom region / (region_pad_factor E) that makeIC.py accepts silently
REGRID_MARGIN = 1.1


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--level",
        "--levelmax",
        dest="level",
        type=int,
        default=5,
        help="levelmax: 2^level gas particles per dimension in the gas box",
    )
    parser.add_argument(
        "--levelmin",
        type=int,
        default=None,
        help="Coarsest background level, spacing L / 2^levelmin (default level - 3)",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=8,
        help="Cells of each level around the next finer one (MUSIC padding)",
    )
    box = parser.add_mutually_exclusive_group()
    box.add_argument(
        "--box-size",
        type=float,
        default=None,
        help="Box side in kpc, or in comoving Mpc with --cosmological",
    )
    box.add_argument(
        "--embed",
        type=float,
        default=None,
        help="Box side in gas box units (default 8)",
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
        help="Width of the tanh profile in gas box units",
    )
    parser.add_argument(
        "--h-max", type=float, default=0.08, help="SPH:h_max in physical kpc"
    )
    parser.add_argument(
        "--pad-factor",
        type=float,
        default=1.3,
        help="ZoomRegion:region_pad_factor",
    )
    parser.add_argument(
        "--max-bkg-cdim",
        type=int,
        default=91,
        help="Cap on ZoomRegion:bkg_top_level_cells",
    )
    parser.add_argument(
        "--max-zoom-cdim",
        type=int,
        default=32,
        help="Cap on the zoom top-level cells per axis",
    )
    parser.add_argument(
        "--max-mesh",
        type=int,
        default=512,
        help="Cap on Gravity:mesh_side_length",
    )
    parser.add_argument(
        "--max-background-particles",
        type=float,
        default=1e8,
        help="Cap on the number of PartType2 particles",
    )
    parser.add_argument(
        "--softening-ratio-background",
        type=float,
        default=0.05,
        help="Gravity:softening_ratio_background",
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
        help="Physical gas and PartType1 softening in kpc (with --cosmological)",
    )
    parser.add_argument(
        "--velocity-fraction",
        type=float,
        default=0.1,
        help="Velocity dispersion in units of the sound speed (with --cosmological)",
    )
    parser.add_argument("--seed", type=int, default=1, help="Random seed")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the layout and the zoom grid only",
    )
    parser.add_argument(
        "-o", dest="output", default="ICs_zoom_homogeneous_box.hdf5", help="IC file"
    )
    parser.add_argument(
        "--zoom-params-file",
        default="zoom_runtime_params.txt",
        help="Runtime parameters for swift -P",
    )
    opt = parser.parse_args()

    if opt.level < 2:
        parser.error("--level must be at least 2")
    if opt.levelmin is None:
        opt.levelmin = opt.level - 3
    if opt.levelmin >= opt.level:
        parser.error("--levelmin must be below --level")
    if opt.padding < 1:
        parser.error("--padding must be at least 1")
    if opt.box_size is None and opt.embed is None:
        opt.embed = 8.0
    if (opt.box_size is not None and opt.box_size <= 0) or (
        opt.embed is not None and opt.embed <= 1
    ):
        parser.error("--box-size must be > 0 and --embed > 1")
    if opt.dm_ratio < 0 or opt.delta < 0 or opt.transition_width <= 0:
        parser.error("--dm-ratio and --delta must be >= 0, the width > 0")
    if opt.velocity_fraction <= 0 and opt.cosmological:
        parser.error("--velocity-fraction must be > 0 with --cosmological")
    return opt


def gas_box(opt):
    """Side, spacing, densities and masses of the gas box, in cgs."""
    n_side = 2**opt.level
    N = n_side**3
    rho_gas = opt.rho * M_PROTON
    m_gas = opt.mass * MSUN
    L = (N * m_gas / rho_gas) ** (1.0 / 3.0)
    m_dm = opt.dm_ratio * m_gas
    rho_in = rho_gas * (1.0 + opt.dm_ratio)

    # Velocity dispersion of the original HomogeneousBox example
    k_J = 2.0 * np.pi / (opt.lJ * L)
    sigma = np.sqrt(4.0 * np.pi * G_NEWTON * rho_gas) / k_J

    return dict(
        N=N,
        L=L,
        dx=L / n_side,
        rho_gas=rho_gas,
        rho_in=rho_in,
        rho_out=rho_in / (1.0 + opt.delta),
        m_gas=m_gas,
        m_dm=m_dm,
        sigma=sigma,
    )


def gas_and_dm_positions(box, rng):
    """Paired gas and dark matter positions in [0, L)^3, in cgs."""
    L, N, dx = box["L"], box["N"], box["dx"]

    # Split each parent point as in space_generate_gas() (src/space.c)
    gas_fraction = box["m_gas"] / (box["m_gas"] + box["m_dm"])
    parents = rng.random((N, 3)) * L
    pos_gas = np.mod(parents - 0.5 * dx * (1.0 - gas_fraction), L)
    pos_dm = np.mod(parents + 0.5 * dx * gas_fraction, L)
    return pos_gas, pos_dm if box["m_dm"] > 0 else np.zeros((0, 3))


def margin(h, padding):
    """Cells added on each side of a hole of h cells: the padding, plus one
    when needed for a multiple of 4 cells, so that the region edges stay on
    the lattice of the next coarser level."""
    return padding if (h // 2 + padding) % 2 == 0 else padding + 1


def level_layout(opt, L, box_target):
    """Background levels from fine to coarse, as MUSIC nests them.

    Level l has the spacing L / 2^l and surrounds the next finer level. The
    coarsest level that fits fills the box. Also returns the coarsest level
    whose padded region still fits, the levelmin that fills the box.
    """
    levels = []
    n_inner = 2**opt.level
    level = opt.level - 1
    while True:
        s = L / 2**level
        h = n_inner // 2
        m = margin(h, opt.padding)
        if (h + 2 * m) * s > box_target * (1.0 + 1e-12):
            break
        if level >= opt.levelmin:
            levels.append(dict(level=level, s=s, h=h, m=m, n=h + 2 * m))
        n_inner = h + 2 * m
        level -= 1
    levelmin_fill = level + 1

    if not levels:
        first = opt.level - 1
        h = 2 ** (opt.level - 1)
        n = h + 2 * margin(h, opt.padding)
        raise SystemExit(
            f"Error: the box ({box_target / PC:.1f} pc) is smaller than the "
            f"padded level {first} ({n * L / 2**first / PC:.1f} pc). "
            "Increase the box or lower --padding."
        )

    outer = levels[-1]
    m_box = max(outer["m"], int(round((box_target / outer["s"] - outer["h"]) / 2)))
    outer["m"] = m_box
    outer["n"] = outer["h"] + 2 * m_box
    outer["fill"] = True
    for lev in levels:
        lev["count"] = lev["n"] ** 3 - lev["h"] ** 3
    return levels, outer["n"] * outer["s"], levelmin_fill


def level_cells(lev):
    """Cell centres of one background level, relative to the box centre."""
    n, h, s = lev["n"], lev["h"], lev["s"]
    idx = (np.arange(n) + 0.5) / n
    gx, gy, gz = np.meshgrid(idx, idx, idx, indexing="ij")
    unit = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel())) - 0.5
    inner = np.all(np.abs(unit) < 0.5 * h / n, axis=1)
    return unit[~inner] * (n * s)


def level_masses(opt, box, lev, cell_centres):
    """Particle masses of one level: the density profile times the cell volume."""
    L = box["L"]
    rho_in, rho_out = box["rho_in"], box["rho_out"]
    r = np.max(np.abs(cell_centres), axis=1)
    x = (r - 0.5 * L) / (opt.transition_width * L)
    return (rho_out + (rho_in - rho_out) * (1.0 - np.tanh(x))) * lev["s"] ** 3


def background_mass(opt, box, levels):
    """Total PartType2 mass, without drawing the particles."""
    L = box["L"]
    total = 0.0
    for lev in levels:
        volume = (lev["n"] ** 3 - lev["h"] ** 3) * lev["s"] ** 3
        inner_edge = 0.5 * lev["h"] * lev["s"] - 0.5 * L
        if (
            box["rho_in"] == box["rho_out"]
            or inner_edge > 20 * opt.transition_width * L
        ):
            total += box["rho_out"] * volume
        else:
            total += level_masses(opt, box, lev, level_cells(lev)).sum()
    return total


def background_particles(opt, box, levels):
    """PartType2 at the cell centres of each level, in cgs.

    No jitter: its random pull on the gas box grows with the spacing.
    """
    centre = 0.5 * levels[-1]["n"] * levels[-1]["s"]
    pos, mass = [], []
    for lev in levels:
        cell_centres = level_cells(lev)
        pos.append(centre + cell_centres)
        mass.append(level_masses(opt, box, lev, cell_centres))
    return np.concatenate(pos), np.concatenate(mass)


def void_geometry(box_size, padded, bkg_cdim):
    """Zoom region cells, side and zoom regions of zoom_get_void_geometry()."""
    bkg_width = box_size / bkg_cdim
    i_low = np.floor((0.5 * box_size - 0.5 * padded) / bkg_width)
    i_up = np.floor((0.5 * box_size + 0.5 * padded) / bkg_width)
    k = int(i_up - i_low + 1)
    region = k * bkg_width
    return k, region, int(np.ceil(region / padded))


def background_cdim(opt, box_size, padded):
    """zoom_optimise_bkg_cells() without its floor of 8, capped (odd)."""
    bkg_cdim = int(np.floor(box_size / padded))
    if bkg_cdim % 2 == 0:
        bkg_cdim -= 1
    if bkg_cdim >= 1:
        width = box_size / bkg_cdim
        if np.floor((0.5 * box_size - 0.5 * padded) / width) != np.floor(
            (0.5 * box_size + 0.5 * padded) / width
        ):
            bkg_cdim -= 2
    cap = opt.max_bkg_cdim if opt.max_bkg_cdim % 2 else opt.max_bkg_cdim - 1
    bkg_cdim = min(bkg_cdim, cap)
    if bkg_cdim < 3:
        raise SystemExit(
            f"Error: the box is only {box_size / padded:.2f} padded zoom regions "
            "wide; the background grid needs at least 3 cells. Increase the box."
        )
    return bkg_cdim


def zoom_depth(opt, bkg_cdim, bkg_width, k, min_zoom_width):
    """Deepest zoom grid with cells wider than min_zoom_width, at most
    max_zoom_cdim cells per axis and not narrower than a mesh cell."""
    depth = 0
    while (
        bkg_width / 2 ** (depth + 1) > min_zoom_width
        and k * 2 ** (depth + 1) <= opt.max_zoom_cdim
        and bkg_cdim * 2 ** (depth + 1) <= opt.max_mesh
    ):
        depth += 1
    return depth


def mesh_size(box_size, bkg_cdim, bkg_width, depth):
    """Smallest power of 2 with mesh cells not wider than a zoom cell
    (space_regrid() stops otherwise) and r_cut_max within a background cell."""
    mesh = 16
    while (
        R_CUT_MAX * A_SMOOTH * box_size / mesh > bkg_width
        or mesh < bkg_cdim * 2**depth
    ):
        mesh *= 2
    return mesh


def zoom_grid(opt, box_size, extent, h_max):
    """Background cdim, zoom depth and mesh, mirroring zoom_region_init()."""
    errors, warnings = [], []
    padded = opt.pad_factor * extent
    bkg_cdim = background_cdim(opt, box_size, padded)
    bkg_width = box_size / bkg_cdim

    k, region, nr_regions = void_geometry(box_size, padded, bkg_cdim)
    if nr_regions >= ZOOM_REGIONS_ERROR:
        needed = bkg_cdim + 2
        while void_geometry(box_size, padded, needed)[2] >= ZOOM_REGIONS_ERROR:
            needed += 2
        errors.append(
            f"{nr_regions} zoom regions in the void region (SWIFT stops at "
            f"{ZOOM_REGIONS_ERROR}): bkg_top_level_cells {bkg_cdim} is capped by "
            f"--max-bkg-cdim (run.sh: max_bkg_cdim); this box needs at least "
            f"{needed} cells"
        )
    elif nr_regions >= ZOOM_REGIONS_WARNING:
        warnings.append(
            f"{nr_regions} zoom regions in the void region (SWIFT warns at "
            f"{ZOOM_REGIONS_WARNING})"
        )

    min_zoom_width = SPACE_STRETCH * KERNEL_GAMMA * h_max
    depth = zoom_depth(opt, bkg_cdim, bkg_width, k, min_zoom_width)
    if depth < 1:
        errors.append(
            f"no zoom depth >= 1 gives zoom cells wider than 1.1 kernel_gamma "
            f"h_max = {min_zoom_width / PC:.1f} pc with at most "
            f"{opt.max_zoom_cdim} cells per axis and a mesh of at most "
            f"{opt.max_mesh} (background cell {bkg_width / PC:.1f} pc)"
        )
        depth = 1
    zoom_width = bkg_width / 2**depth
    if region < extent:
        errors.append("zoom region smaller than the high-resolution particles")
    if region / extent < REGRID_MARGIN * opt.pad_factor:
        warnings.append(
            f"the zoom region regrids when the extent grows by "
            f"{100 * (region / (opt.pad_factor * extent) - 1):.1f} % "
            f"(zoom region / region_pad_factor); keep at least "
            f"{100 * (REGRID_MARGIN - 1):.0f} %: change the box or region_pad_factor"
        )
    if region / extent >= 2.0 * opt.pad_factor:
        warnings.append(
            f"real padding {region / extent:.2f} is at least twice "
            f"region_pad_factor {opt.pad_factor} (SWIFT warns)"
        )

    mesh = mesh_size(box_size, bkg_cdim, bkg_width, depth)
    if mesh > opt.max_mesh:
        mesh = opt.max_mesh
        warnings.append(
            f"mesh capped at {mesh}: r_cut_max is larger than a background cell"
        )
    if mesh < bkg_cdim * 2**depth:
        errors.append(
            f"mesh {mesh} coarser than the zoom cells ({bkg_cdim * 2**depth} "
            "zoom cells across the box)"
        )
    r_s = A_SMOOTH * box_size / mesh

    return dict(
        bkg_cdim=bkg_cdim,
        bkg_width=bkg_width,
        k=k,
        region=region,
        nr_regions=nr_regions,
        depth=depth,
        zoom_cdim=k * 2**depth,
        zoom_width=zoom_width,
        min_zoom_width=min_zoom_width,
        mesh=mesh,
        r_s=r_s,
        r_cut=R_CUT_MAX * r_s,
        errors=errors,
        warnings=warnings,
    )


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


def rho_crit_0(opt):
    H0 = 100.0 * opt.hubble_param * 1e5 / MPC
    return 3.0 * H0**2 / (8.0 * np.pi * G_NEWTON)


def physical_box_target(opt, box):
    """Requested physical box side in cgs."""
    L = box["L"]
    if opt.embed is not None:
        return opt.embed * L
    if not opt.cosmological:
        return opt.box_size * KPC

    # Mean density of the box = Omega_m rho_crit_0 at a_begin, box comoving
    B = opt.box_size * MPC
    excess = (box["rho_in"] - box["rho_out"]) * (L / B) ** 3
    a3 = (opt.omega_m * rho_crit_0(opt) - excess) / box["rho_out"]
    if a3 <= 0:
        raise SystemExit("Error: the comoving box is too small for this overdensity")
    return a3 ** (1.0 / 3.0) * B


def print_levels(opt, box, levels, rho_mean, n_bg_total):
    L, dx = box["L"], box["dx"]
    unit = "pc" if not opt.cosmological else "pc (physical)"
    print(
        f"levelmax {opt.level}  : gas box L = {L / PC:.2f} {unit}, spacing "
        f"d = {dx / PC:.3f} pc, {box['N']} gas and {box['N'] if box['m_dm'] > 0 else 0} "
        "PartType1 particles"
    )
    for lev in levels:
        m = lev["s"] ** 3 * box["rho_out"]
        eps = opt.softening_ratio_background * (m / rho_mean) ** (1.0 / 3.0)
        kind = "fills the box" if lev.get("fill") else f"padding {lev['m']} cells"
        print(
            f"level {lev['level']:3d}   : spacing {lev['s'] / dx:g} d = "
            f"{lev['s'] / PC:.4g} pc, region {lev['n']} cells = "
            f"{lev['n'] * lev['s'] / L:g} L ({kind}), {lev['count']} PartType2 "
            f"of {m / MSUN:.4g} Msun (rho_out), softening {eps / PC:.3g} pc"
        )
    print(
        f"PartType2 total : {n_bg_total} particles on {len(levels)} levels "
        f"({opt.level - 1} to {levels[-1]['level']})"
    )


def main():
    opt = parse_options()
    rng = np.random.default_rng(opt.seed)

    box = gas_box(opt)
    L, N = box["L"], box["N"]
    target = physical_box_target(opt, box)
    levels, box_size, levelmin_fill = level_layout(opt, L, target)
    n_bg_total = sum(lev["count"] for lev in levels)
    offset = 0.5 * box_size - 0.5 * L
    centre = 0.5 * box_size
    h_max = opt.h_max * KPC
    h = min(ETA_RESOLUTION * box["dx"], h_max)

    if n_bg_total > opt.max_background_particles:
        print_levels(opt, box, levels, box["rho_out"], n_bg_total)
        raise SystemExit(
            f"Error: {n_bg_total:.3e} PartType2 particles (limit "
            f"--max-background-particles {opt.max_background_particles:.3g}). "
            f"--levelmin {levelmin_fill} fills this box with padded levels."
        )

    M_hires = N * (box["m_gas"] + box["m_dm"])
    particles = None
    if opt.dry_run:
        extent = L
        M_total = M_hires + background_mass(opt, box, levels)
    else:
        pos_gas, pos_dm = gas_and_dm_positions(box, rng)
        pos_bg, m_bg = background_particles(opt, box, levels)
        pos_gas += offset
        pos_dm += offset
        hires = np.concatenate((pos_gas, pos_dm))
        extent = np.max(hires.max(axis=0) - hires.min(axis=0))

        if np.any(np.max(np.abs(pos_bg - centre), axis=1) < 0.5 * L):
            raise RuntimeError("PartType2 particle inside the gas box")
        if np.any(pos_bg < 0) or np.any(pos_bg >= box_size):
            raise RuntimeError("PartType2 particle outside the periodic box")
        M_total = M_hires + m_bg.sum()
        particles = (pos_gas, pos_dm, pos_bg, m_bg)

    rho_mean = M_total / box_size**3
    if opt.cosmological:
        a = (opt.omega_m * rho_crit_0(opt) / rho_mean) ** (1.0 / 3.0)
    else:
        a = 1.0
    grid = zoom_grid(opt, box_size, extent, h_max)

    # Report
    mode = f"comoving, a_begin = {a:.6e}, z = {1 / a - 1:.2f}"
    print(
        f"{opt.output}{' (dry run, not written)' if opt.dry_run else ''}: "
        f"{mode if opt.cosmological else 'physical'}"
    )
    if opt.box_size:
        requested = f"requested {opt.box_size:.6g}"
    else:
        requested = f"requested {opt.embed:g} gas boxes"
    if opt.cosmological:
        print(
            f"Box             : {box_size / a / MPC:.6g} comoving Mpc ({requested}), "
            f"{box_size / KPC:.6g} physical kpc at a_begin, "
            f"{box_size / L:.6g} gas boxes"
        )
        print(
            f"Gas box         : {L / KPC:.6g} physical kpc = "
            f"{L / a / KPC:.6g} comoving kpc"
        )
    else:
        print(
            f"Box             : {box_size / KPC:.6g} kpc ({requested}), "
            f"{box_size / L:.6g} gas boxes"
        )
    print(
        f"Levels          : levelmax {opt.level}, levelmin {opt.levelmin}, "
        f"padding {opt.padding}, coarsest level used {levels[-1]['level']}, "
        f"levelmin that fills this box with padded levels {levelmin_fill}"
    )
    print_levels(opt, box, levels, rho_mean, n_bg_total)
    print(f"Mean spacing    : {box['dx'] / PC:.3f} pc, h = {h / PC:.3f} pc")
    c_s = np.sqrt(GAS_GAMMA * (GAS_GAMMA - 1.0)) * box["sigma"]
    print(
        f"Gas             : u = sigma^2, sigma = {box['sigma'] / 1e5:.4f} km/s, "
        f"sound speed {c_s / 1e5:.4f} km/s"
    )
    t_ff = np.sqrt(3.0 * np.pi / (32.0 * G_NEWTON * box["rho_in"]))
    print(f"Free-fall time  : {t_ff / UNIT_TIME:.4f} internal units (gas box)")
    print(
        f"Densities       : inside {box['rho_in'] / M_PROTON:.6f}, outside "
        f"{box['rho_out'] / M_PROTON:.6f}, mean {rho_mean / M_PROTON:.6f} m_p/cm^3"
        " (physical)"
    )
    print(f"Total mass      : {M_total / MSUN:.6e} Msun")

    eps_dm = opt.softening if opt.cosmological else None
    print(
        f"Softening       : PartType2 ratio {opt.softening_ratio_background} of the "
        "mean spacing at its mass (values above)"
        + (f"; gas and PartType1 {eps_dm * 1e3:g} pc" if eps_dm else "")
    )
    print(
        f"Zoom grid       : extent E = {extent / PC:.2f} pc, bkg_top_level_cells "
        f"{grid['bkg_cdim']} (cell {grid['bkg_width'] / PC:.6g} pc, cap "
        f"{opt.max_bkg_cdim}), zoom region {grid['k']} background cells = "
        f"{grid['region'] / PC:.6g} pc (real padding {grid['region'] / extent:.3f}, "
        f"{grid['nr_regions']} zoom regions in the void region)"
    )
    print(
        f"                  zoom_top_level_depth {grid['depth']}, zoom cdim "
        f"{grid['zoom_cdim']}, zoom cell {grid['zoom_width'] / PC:.6g} pc "
        f"(> 1.1 kernel_gamma h_max = {grid['min_zoom_width'] / PC:.1f} pc), "
        f"region_pad_factor {opt.pad_factor}"
    )
    print(
        f"                  mesh_side_length {grid['mesh']} (mesh cell "
        f"{box_size / grid['mesh'] / PC:.6g} pc, r_s {grid['r_s'] / PC:.6g} pc, "
        f"r_cut_max {grid['r_cut'] / PC:.6g} pc)"
    )
    for w in grid["warnings"]:
        print(f"Warning: {w}")
    if grid["errors"]:
        for e in grid["errors"]:
            print(f"Error: {e}", file=sys.stderr)
        raise SystemExit("Error: the zoom grid does not fit the rules (see above)")
    if particles is None:
        return
    pos_gas, pos_dm, pos_bg, m_bg = particles

    print(f"PartType2 inside the gas box: 0 (checked); in the box: {len(m_bg)}")

    to_mass = 1.0 / UNIT_MASS
    to_length = 1.0 / UNIT_LENGTH
    N_dm, N_bg = len(pos_dm), len(m_bg)
    groups = {
        0: {
            "Coordinates": pos_gas * to_length,
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
            "Coordinates": pos_dm * to_length,
            "Masses": np.full(N_dm, box["m_dm"] * to_mass),
            "ParticleIDs": np.arange(N + 1, N + N_dm + 1),
        }
    groups = dict(sorted(groups.items()))

    # SWIFT stops on a zero RMS velocity in cosmological runs
    # (engine_recompute_displacement_constraint), so these ICs get a small
    # random velocity
    sigma_v = opt.velocity_fraction * c_s / UNIT_VELOCITY if opt.cosmological else 0.0
    for data in groups.values():
        shape = data["Coordinates"].shape
        data["Velocities"] = (
            rng.normal(0.0, sigma_v, shape) if sigma_v > 0 else np.zeros(shape)
        )
    write_ics(opt, opt.output, a, box_size * to_length, groups)
    if opt.cosmological:
        print(f"Velocities      : Gaussian, {sigma_v:.4g} km/s per component")

    runtime = [
        f"ZoomRegion:bkg_top_level_cells:{grid['bkg_cdim']}",
        f"ZoomRegion:zoom_top_level_depth:{grid['depth']}",
        f"ZoomRegion:region_pad_factor:{opt.pad_factor}",
        f"Gravity:mesh_side_length:{grid['mesh']}",
        f"Gravity:softening_ratio_background:{opt.softening_ratio_background}",
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
            f"Statistics:scale_factor_first:{a:.8e}",
            "Statistics:delta_time:1.01",
        ]
    with open(opt.zoom_params_file, "w") as f:
        f.write(" ".join(f"-P {r}" for r in runtime) + "\n")
    print(f"Runtime parameters written to {opt.zoom_params_file}:")
    for r in runtime:
        print(f"  {r}")


if __name__ == "__main__":
    main()
