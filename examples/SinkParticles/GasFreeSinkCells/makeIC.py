################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
"""Generate ICs with a gas-free, sink-only corner cell.

The box reproduces the HomogeneousBox uniform gas background (same unit
system, gas density, and particle mass) everywhere except inside one cubic
corner region, sized to coincide exactly with one top-level cell (see
``Scheduler:max_top_level_cells`` in ``params.yml``). That corner region
holds zero gas particles and a tight cluster of sink particles, placed
in the ICs.

The sink cluster is deliberately dense enough that:

* its gravitational particle count exceeds ``Scheduler:cell_split_size``
  (raised to 800 in ``params.yml``, above the ~500 gas particles per
  top-level cell in the bulk of the box, so only the sink-only corner
  cell keeps splitting) by a large factor, forcing the tree builder to
  recursively split the gas-free corner cell (``src/space_split.c``
  splits on ``gcount`` alone when self-gravity is on, regardless of
  ``hydro.count``);
* pair separations within the cluster are well below the sink cut-off
  radius (``GEARSink:cut_off_radius``), so sink-sink mergers are actually
  attempted once the simulation runs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import h5py
import numpy as np
from astropy import units
from astropy.constants import G, m_p
from numpy.typing import NDArray

# Standard internal unit system (identical to HomogeneousBox).
UNIT_MASS_IN_CGS = 1.988409870698051e43  # 10^10 Msun in grams
UNIT_LENGTH_IN_CGS = 3.0856775814913673e21  # 1 kpc in centimeters
UNIT_VELOCITY_IN_CGS = 1e5  # km/s in centimeters per second
UNIT_CURRENT_IN_CGS = 1.0  # Amperes
UNIT_TEMP_IN_CGS = 1.0  # Kelvin
UNIT_TIME_IN_CGS = UNIT_LENGTH_IN_CGS / UNIT_VELOCITY_IN_CGS

UNIT_MASS = UNIT_MASS_IN_CGS * units.g
UNIT_LENGTH = UNIT_LENGTH_IN_CGS * units.cm
UNIT_VELOCITY = UNIT_VELOCITY_IN_CGS * units.cm / units.s

# Kernel support radius factor for the Wendland C2 kernel in 3D (must match
# `kernel_gamma` for WENDLAND_C2_KERNEL/HYDRO_DIMENSION_3D in
# src/kernel_hydro.h). Only used to size the placeholder sink smoothing
# length; GEARSink overwrites it from `cut_off_radius` at start-up whenever
# `use_fixed_cut_off_radius` is set (src/sink/GEAR/sink.h::sink_first_init_sink).
WENDLAND_C2_KERNEL_GAMMA_3D = 1.936492

# `space_stretch` in src/space.h: safety factor applied to h_max*kernel_gamma
# when space_regrid.c derives the top-level cell count from the largest
# smoothing length in the box.
SPACE_STRETCH = 1.10


@dataclass
class Geometry:
    """Store the derived geometric quantities of the box, in internal units.

    Parameters
    ----------
    box_size : float
        Side length of the periodic cubic box.
    corner_side : float
        Side length of the gas-free cubic corner region, placed at the
        box origin so that it coincides with the [0, 0, 0] top-level cell.
    cluster_center : NDArray[np.float64]
        Coordinates of the centre of the sink cluster.
    cluster_side : float
        Side length of the cube in which sink particles are drawn.
    """

    box_size: float
    corner_side: float
    cluster_center: NDArray[np.float64]
    cluster_side: float


def parse_args() -> argparse.Namespace:
    """Parse the command-line options for the IC generator.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate a periodic box with a uniform gas background "
        "and one gas-free corner cell packed with sink particles."
    )

    parser.add_argument(
        "--lJ",
        type=float,
        default=0.250,
        help="Jeans wavelength in box size unit, used only to set the gas "
        "internal energy (same convention as HomogeneousBox).",
    )
    parser.add_argument(
        "--rho",
        type=float,
        default=0.1,
        help="Mean gas density in atom/cm3.",
    )
    parser.add_argument(
        "--mass",
        type=float,
        default=50.0,
        help="Gas particle mass in solar masses.",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=5,
        help="Nominal resolution level: N = (2**level)**3 gas particles "
        "before carving out the sink-only corner region.",
    )
    parser.add_argument(
        "--corner-frac",
        type=float,
        default=0.25,
        help="Side length of the gas-free corner cube, as a fraction of "
        "the box size. Must equal 1 / Scheduler:max_top_level_cells in "
        "params.yml so the corner maps exactly onto one top-level cell.",
    )
    parser.add_argument(
        "--max-top-level-cells",
        type=int,
        default=4,
        help="Scheduler:max_top_level_cells in params.yml. Used here only "
        "to check that the requested geometry actually yields this many "
        "top-level cells per dimension (see check_cdim_matches_request).",
    )
    parser.add_argument(
        "--n-sinks",
        type=int,
        default=2000,
        help="Number of sink particles placed in the corner cluster "
        "(should be several times Scheduler:cell_split_size, set to 800 "
        "in params.yml here, to force recursive splitting of the "
        "gas-free cell).",
    )
    parser.add_argument(
        "--cluster-frac",
        type=float,
        default=0.12,
        help="Side length of the cube containing the sink cluster, as a "
        "fraction of the corner cube side length.",
    )
    parser.add_argument(
        "--cluster-offset-frac",
        type=float,
        default=0.3,
        help="Position of the cluster centre along each axis, as a "
        "fraction of the corner cube side length. Kept away from 0.5 so "
        "the cluster sits inside a single octant at every recursive tree "
        "split, exercising several levels of sink-only progeny.",
    )
    parser.add_argument(
        "--cut-off-radius",
        type=float,
        default=5e-3,
        help="GEARSink cut_off_radius in internal length units (kpc). "
        "Must match GEARSink:cut_off_radius in params.yml. Used here only "
        "to size the cluster and to sanity-check the merger geometry.",
    )
    parser.add_argument(
        "--sink-mass-factor",
        type=float,
        default=5.0,
        help="Sink particle mass, as a multiple of the gas particle mass.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed of the random number generator.",
    )
    parser.add_argument(
        "-o",
        dest="outputfilename",
        type=str,
        default="ICs_gas_free_sink_cells.hdf5",
        help="Output filename.",
    )

    return parser.parse_args()


def compute_box_size(rho_atom_cm3: float, mass_msun: float, n_gas: int) -> float:
    """Compute the box side length giving the requested mean gas density.

    Parameters
    ----------
    rho_atom_cm3 : float
        Mean gas number density, in hydrogen atoms per cm3.
    mass_msun : float
        Mass of a single gas particle, in solar masses.
    n_gas : int
        Nominal number of gas particles used to size the box (as if the
        whole box, corner included, were filled with gas).

    Returns
    -------
    float
        Box side length, in internal length units (kpc).
    """
    rho = rho_atom_cm3 * m_p / units.cm**3
    m = mass_msun * units.Msun
    total_mass = n_gas * m
    box_size = (total_mass / rho) ** (1.0 / 3.0)
    return box_size.to(UNIT_LENGTH).value


def build_geometry(args: argparse.Namespace, box_size: float) -> Geometry:
    """Derive the corner and cluster geometry from the command-line options.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    box_size : float
        Box side length, in internal length units.

    Returns
    -------
    Geometry
        The derived geometric quantities describing the gas-free corner
        cell and the sink cluster it contains.
    """
    corner_side = args.corner_frac * box_size
    cluster_side = args.cluster_frac * corner_side
    offset = args.cluster_offset_frac * corner_side
    cluster_center = np.array([offset, offset, offset])

    return Geometry(
        box_size=box_size,
        corner_side=corner_side,
        cluster_center=cluster_center,
        cluster_side=cluster_side,
    )


def predict_top_level_cdim(
    box_size: float, h_max: float, max_top_level_cells: int
) -> int:
    """Predict the number of top-level cells per dimension SWIFT will use.

    Reproduces the ``cdim`` derivation in ``src/space_regrid.c``: the
    requested ``Scheduler:max_top_level_cells`` is only honoured if the
    largest smoothing length in the box is small enough; otherwise SWIFT
    silently falls back to fewer, larger top-level cells.

    Parameters
    ----------
    box_size : float
        Side length of the (cubic) box, in internal length units.
    h_max : float
        Largest smoothing length of any particle in the ICs, in internal
        length units.
    max_top_level_cells : int
        Requested ``Scheduler:max_top_level_cells``.

    Returns
    -------
    int
        The number of top-level cells per dimension SWIFT would actually
        build for this box.
    """
    tol = max(1.0 - 1.0 / (max_top_level_cells * max_top_level_cells), 0.99)
    cell_min = tol * box_size / max_top_level_cells
    denom = max(h_max * WENDLAND_C2_KERNEL_GAMMA_3D * SPACE_STRETCH, cell_min)
    return min(max_top_level_cells, int(box_size / denom))


def generate_background_gas(
    args: argparse.Namespace, geometry: Geometry, n_gas_nominal: int
) -> NDArray[np.float64]:
    """Generate the uniform gas background, excluding the corner region.

    Positions are drawn uniformly in the full box and particles falling
    inside the gas-free corner cube are rejected, so the surviving
    particles are exactly uniform (in a Poisson-sampling sense) over the
    complement of the corner cube.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    geometry : Geometry
        The box and corner geometry.
    n_gas_nominal : int
        Number of candidate positions to draw before rejection.

    Returns
    -------
    NDArray[np.float64]
        Gas positions with shape (n_gas, 3), where n_gas is the number of
        particles surviving rejection of the corner cube region.
    """
    rng = np.random.default_rng(args.seed)

    candidate_pos = rng.random((n_gas_nominal, 3)) * geometry.box_size
    inside_corner = np.all(candidate_pos < geometry.corner_side, axis=1)
    pos = candidate_pos[~inside_corner]

    return pos


def generate_sink_cluster(
    args: argparse.Namespace, geometry: Geometry, rng: np.random.Generator
) -> NDArray[np.float64]:
    """Draw sink positions in a tight cube inside the gas-free corner.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    geometry : Geometry
        The box and corner geometry.
    rng : numpy.random.Generator
        Random number generator, shared with the gas sampling so that
        `--seed` controls the whole IC reproducibly.

    Returns
    -------
    NDArray[np.float64]
        Sink positions, shape (n_sinks, 3).
    """
    half_side = 0.5 * geometry.cluster_side
    low = geometry.cluster_center - half_side
    offsets = rng.random((args.n_sinks, 3)) * geometry.cluster_side
    return low + offsets


def write_ic_file(
    filename: str,
    geometry: Geometry,
    gas_pos: NDArray[np.float64],
    gas_mass: float,
    gas_h: float,
    gas_u: float,
    gas_rho: float,
    sink_pos: NDArray[np.float64],
    sink_mass: float,
    sink_h: float,
) -> None:
    """Write the gas and sink particle groups to an HDF5 IC file.

    Parameters
    ----------
    filename : str
        Output HDF5 file path.
    geometry : Geometry
        The box geometry (used for the header BoxSize).
    gas_pos : NDArray[np.float64]
        Gas particle positions, shape (n_gas, 3).
    gas_mass : float
        Gas particle mass, in internal mass units.
    gas_h : float
        Gas particle smoothing length, in internal length units.
    gas_u : float
        Gas particle internal energy, in internal energy units.
    gas_rho : float
        Gas particle density, in internal density units.
    sink_pos : NDArray[np.float64]
        Sink particle positions, shape (n_sinks, 3).
    sink_mass : float
        Sink particle mass, in internal mass units.
    sink_h : float
        Sink particle smoothing length, in internal length units.
    """
    n_gas = gas_pos.shape[0]
    n_sinks = sink_pos.shape[0]

    gas_ids = np.arange(n_gas, dtype=np.uint64)
    sink_ids = np.arange(n_gas, n_gas + n_sinks, dtype=np.uint64)

    with h5py.File(filename, "w") as f:
        # Header. Arrays use the current swift_type_count = 7 ordering:
        # [gas, dark_matter, dark_matter_background, sink, stars,
        #  black_hole, neutrino] (src/part_type.h).
        grp = f.create_group("/Header")
        grp.attrs["BoxSize"] = [geometry.box_size] * 3
        grp.attrs["NumPart_Total"] = [n_gas, 0, 0, n_sinks, 0, 0, 0]
        grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0, 0]
        grp.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, n_sinks, 0, 0, 0]
        grp.attrs["Time"] = 0.0
        grp.attrs["NumFileOutputsPerSnapshot"] = 1
        grp.attrs["MassTable"] = [0.0] * 7
        grp.attrs["Flag_Entropy_ICs"] = [0] * 7
        grp.attrs["Dimension"] = 3

        # Units.
        grp = f.create_group("/Units")
        grp.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH_IN_CGS
        grp.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS_IN_CGS
        grp.attrs["Unit time in cgs (U_t)"] = UNIT_TIME_IN_CGS
        grp.attrs["Unit current in cgs (U_I)"] = UNIT_CURRENT_IN_CGS
        grp.attrs["Unit temperature in cgs (U_T)"] = UNIT_TEMP_IN_CGS

        # Gas particles.
        grp = f.create_group("/PartType0")
        grp.create_dataset("Coordinates", data=gas_pos, dtype="d")
        grp.create_dataset("Velocities", data=np.zeros((n_gas, 3)), dtype="f")
        grp.create_dataset("Masses", data=np.full(n_gas, gas_mass), dtype="f")
        grp.create_dataset("SmoothingLength", data=np.full(n_gas, gas_h), dtype="f")
        grp.create_dataset("InternalEnergy", data=np.full(n_gas, gas_u), dtype="f")
        grp.create_dataset("ParticleIDs", data=gas_ids, dtype="L")
        grp.create_dataset("Densities", data=np.full(n_gas, gas_rho), dtype="f")

        # Sink particles: read directly by
        # src/sink/GEAR/sink_io.h::sink_read_particles (Coordinates,
        # Velocities, Masses, ParticleIDs are compulsory; SmoothingLength
        # and BirthTime are optional but provided for robustness).
        grp = f.create_group("/PartType3")
        grp.create_dataset("Coordinates", data=sink_pos, dtype="d")
        grp.create_dataset("Velocities", data=np.zeros((n_sinks, 3)), dtype="f")
        grp.create_dataset("Masses", data=np.full(n_sinks, sink_mass), dtype="f")
        grp.create_dataset("SmoothingLength", data=np.full(n_sinks, sink_h), dtype="f")
        grp.create_dataset("ParticleIDs", data=sink_ids, dtype="L")
        grp.create_dataset("BirthTime", data=np.zeros(n_sinks), dtype="f")


def main() -> None:
    """Generate and write the gas-free sink-cell test ICs."""
    args = parse_args()

    n_gas_nominal = (2**args.level) ** 3
    box_size = compute_box_size(args.rho, args.mass, n_gas_nominal)
    geometry = build_geometry(args, box_size)

    # Gas particle properties (mirrors HomogeneousBox/makeIC.py).
    rho = (args.rho * m_p / units.cm**3).to(UNIT_MASS / UNIT_LENGTH**3).value
    gas_mass = (args.mass * units.Msun).to(UNIT_MASS).value
    gas_h = 3.0 * box_size / n_gas_nominal ** (1.0 / 3.0)

    lJ = args.lJ * box_size * UNIT_LENGTH
    kJ = 2.0 * np.pi / lJ
    sigma = (
        np.sqrt(4.0 * np.pi * G * (args.rho * m_p / units.cm**3))
        / kJ
    )
    sigma = sigma.to(UNIT_VELOCITY).value
    gas_u = sigma**2

    rng = np.random.default_rng(args.seed)
    gas_pos = generate_background_gas(args, geometry, n_gas_nominal)
    sink_pos = generate_sink_cluster(args, geometry, rng)

    sink_mass = args.sink_mass_factor * gas_mass
    sink_h = args.cut_off_radius / WENDLAND_C2_KERNEL_GAMMA_3D

    predicted_cdim = predict_top_level_cdim(
        box_size, max(gas_h, sink_h), args.max_top_level_cells
    )
    if predicted_cdim != args.max_top_level_cells:
        raise ValueError(
            f"Requested {args.max_top_level_cells} top-level cells per "
            f"dimension, but the gas smoothing length ({gas_h:.6g}) is too "
            f"large relative to the box size ({box_size:.6g}): SWIFT would "
            f"actually build {predicted_cdim}. The gas-free corner cube "
            "would then no longer coincide with a single top-level cell. "
            "Lower --level (finer resolution -> smaller h) or "
            "--max-top-level-cells to match."
        )

    write_ic_file(
        args.outputfilename,
        geometry,
        gas_pos,
        gas_mass,
        gas_h,
        gas_u,
        rho,
        sink_pos,
        sink_mass,
        sink_h,
    )

    n_gas = gas_pos.shape[0]
    mean_gas_sep = box_size / n_gas_nominal ** (1.0 / 3.0)
    cluster_mean_sep = geometry.cluster_side / args.n_sinks ** (1.0 / 3.0)

    print(f"{args.outputfilename} saved.")
    print(f"Predicted top-level cdim per dimension : {predicted_cdim}")
    print(f"Box size (internal units)              : {box_size:.6g}")
    print(f"Gas particles                          : {n_gas} (nominal {n_gas_nominal})")
    print(f"Mean gas inter-particle separation      : {mean_gas_sep:.6g}")
    print(f"Corner cube side (gas-free region)      : {geometry.corner_side:.6g}")
    print(f"Sink particles                          : {args.n_sinks}")
    print(f"Sink cluster side                       : {geometry.cluster_side:.6g}")
    print(f"Sink cluster centre                     : {geometry.cluster_center}")
    print(f"Sink cutoff radius (GEARSink)            : {args.cut_off_radius:.6g}")
    print(f"Sink cluster mean nearest-neighbour scale: {cluster_mean_sep:.6g}")
    print(f"Sink mass                               : {sink_mass:.6g}")
    print(f"Gas particle mass                       : {gas_mass:.6g}")


if __name__ == "__main__":
    main()
