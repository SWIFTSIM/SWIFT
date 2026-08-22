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
"""
Shared reader for the Data/Radiation group of a GEAR yields_table HDF5 file
(e.g. POPIIsw.h5), used by the SubgridRadiation examples' *_analytic_check.py
scripts to look up the ionizing photon rate Q_H(mass), the bolometric
luminosity L(mass) and the mean excess photon energy above 13.6 eV, from the
actual pychem-generated table -- instead of the standalone piecewise-fit
Python ports these scripts used before the table migration.

Replicates the read/interpolate chain src/feedback/GEAR/radiation.c actually
runs, not pychem's own generation-time interpolation convention:

  * The table is read exactly like radiation_read_grid_metadata() +
    radiation_build_tables(): the group-level "dimensionality"/"m0"/"dm"/"nm"
    attributes and the Luminosity/Q_H/DotEExcess datasets, all still plain
    CGS doubles (see radiation.c's own doxygen for why no SWIFT-internal
    unit conversion or RADIATION_DOT_N_ION_TABLE_SCALING is needed here: this
    module stays in CGS throughout, matching the pre-conversion checker
    scripts' own convention, and those two factors would cancel out again on
    the way back to CGS in any case).
  * Interpolation is SWIFT's own semi-log scheme (interpolation.h's
    interpolate_1d_init()/interpolate_1d()): linear in log10(mass), linear
    in the raw physical value -- NOT pychem's own log10(mass)-vs-log10(value)
    "np.interp" convention. The two schemes were measured against each
    other for this migration's in-scope examples and found to differ by
    <0.2%, well under the 1% trigger the plan set for requiring a separate
    log-log implementation here (see
    .claude/dev/logs/2026-08-22_0253_phase2-parity-harness.md, item 2) --
    so replicating SWIFT's own convention, not pychem's, is deliberate.
  * Same two-stage structure as radiation_build_tables(): the native table
    (evenly spaced in log10(mass), "m0"/"dm"/"nm") is first resampled onto
    "radiation_interpolation_size_mass" points evenly spaced in log10(mass)
    between the table's own Data/IMF "Mmin"/"Mmax" (NOT the native table's
    own mass range -- this is a deliberate SWIFT design choice, see
    radiation_build_tables()'s doxygen), then a query at an arbitrary mass
    linearly interpolates that resampled grid. Both stages clamp to the
    edge value outside their domain (boundary_condition_const) -- never
    extrapolate, never raise.

Only 1D ("M", fits-mode) tables are supported: every SubgridRadiation
example this module is used by loads a fits-mode table. A 2D ("M,Z")
table raises, mirroring radiation_check_is_1d()'s abort in the C code.

The table path is never hardcoded: yields_table_path_from_snapshot() reads
GEARFeedback:yields_table from the run's own Parameters group (the resolved
value SWIFT actually used, exactly as parsed from params.yml, including
any command-line override) -- the same, deliberately uncanonicalized,
cwd-relative path SWIFT's own H5Fopen call resolves it against (see
PLAN_radiation_tables.md finding #12); this module does not attempt to fix
that, only to read the file the run itself read.
"""

import os

import h5py
import numpy as np

#: Default for GEARFeedback:radiation_interpolation_size_mass
#: (src/feedback/GEAR/radiation.c, examples/parameter_example.yml).
DEFAULT_INTERPOLATION_SIZE_MASS = 200


def _decode(raw):
    """Decode an HDF5 attribute value that may come back as bytes."""
    if isinstance(raw, (bytes, np.bytes_)):
        return raw.decode()
    if isinstance(raw, np.ndarray):
        raw = raw.item()
        return raw.decode() if isinstance(raw, (bytes, np.bytes_)) else raw
    return raw


def yields_table_path_from_snapshot(snapshot_file):
    """Return the raw GEARFeedback:yields_table string this run actually
    used, read from the snapshot's own Parameters group -- exactly as
    parsed from params.yml (relative or absolute, whatever the run was
    given), not yet resolved to an existing file; see
    resolve_yields_table_path().

    @param snapshot_file One of the run's snapshot files (any one holds the
    same Parameters group).
    """
    with h5py.File(snapshot_file, "r") as h:
        return _decode(h["Parameters"].attrs["GEARFeedback:yields_table"])


def resolve_yields_table_path(raw_path, snapshot_file, script_dir=None):
    """Resolve GEARFeedback:yields_table's raw, uncanonicalized path
    (PLAN_radiation_tables.md finding #12: SWIFT stores it exactly as
    given, relative to whatever cwd the run itself used) to an existing
    file, without hardcoding a location.

    Tried in order: as given, relative to the current working directory
    (matches SWIFT's own H5Fopen resolution, the common case when this
    check script runs right after run.sh in the same directory); relative
    to the snapshot file's own directory and that directory's parent
    (covers a snapshot glob pointed at an archived/moved run directory,
    e.g. run.sh's `mv snap $run_name`, where the table file was left
    behind one level up rather than moved with the snapshots); relative to
    `script_dir`, the calling check script's own directory (covers the 5
    examples with a table checked in next to the script itself).

    @param raw_path The string from yields_table_path_from_snapshot().
    @param snapshot_file One of the run's snapshot files (for the
    snapshot-directory candidates).
    @param script_dir Directory of the calling check script, or None to
    skip that candidate.
    @return The first candidate that exists.
    @raise FileNotFoundError naming every path tried, if none exists.
    """
    if os.path.isabs(raw_path):
        candidates = [raw_path]
    else:
        snap_dir = os.path.dirname(os.path.abspath(snapshot_file))
        candidates = [
            raw_path,
            os.path.join(snap_dir, raw_path),
            os.path.join(os.path.dirname(snap_dir), raw_path),
        ]
        if script_dir is not None:
            candidates.append(os.path.join(script_dir, raw_path))

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate

    raise FileNotFoundError(
        f"Could not find the radiation table 'GEARFeedback:yields_table="
        f"{raw_path}'. Tried: {candidates}"
    )


def open_radiation_table(
    files,
    script_dir=None,
    default_interpolation_size_mass=DEFAULT_INTERPOLATION_SIZE_MASS,
):
    """Open the RadiationTable for a run, from its own recorded
    GEARFeedback:yields_table path -- never a hardcoded filename.

    @param files The run's sorted snapshot filenames (any one holds the
    Parameters group this reads).
    @param script_dir Directory of the calling check script (see
    resolve_yields_table_path()), or None to skip that fallback candidate.
    @param default_interpolation_size_mass Used when the snapshot predates
    GEARFeedback:radiation_interpolation_size_mass (see
    interpolation_size_mass_from_snapshot()).
    """
    raw_path = yields_table_path_from_snapshot(files[0])
    resolved_path = resolve_yields_table_path(raw_path, files[0], script_dir)
    interpolation_size_mass = interpolation_size_mass_from_snapshot(
        files[0], default_interpolation_size_mass
    )
    return RadiationTable(resolved_path, interpolation_size_mass)


def interpolation_size_mass_from_snapshot(
    snapshot_file, default=DEFAULT_INTERPOLATION_SIZE_MASS
):
    """Return GEARFeedback:radiation_interpolation_size_mass this run used,
    or `default` if the snapshot predates that key (e.g. an old archived
    run made before the table migration)."""
    with h5py.File(snapshot_file, "r") as h:
        raw = h["Parameters"].attrs.get(
            "GEARFeedback:radiation_interpolation_size_mass"
        )
    return int(_decode(raw)) if raw is not None else default


def _resample_1d(
    log_mass_min_out, log_mass_max_out, n_out, log_data_xmin, step_size, data_native
):
    """Direct port of interpolate_1d_init() (interpolation.h), with
    boundary_condition_const: resample `data_native` (evenly spaced in
    log10(mass), starting at log_data_xmin with step step_size) onto n_out
    points evenly spaced in log10(mass) in [log_mass_min_out,
    log_mass_max_out].

    @return (resampled data array of length n_out, dx between output points).
    """
    n_data = len(data_native)
    dx = (log_mass_max_out - log_mass_min_out) / (n_out - 1)
    out = np.empty(n_out, dtype=np.float64)

    for i in range(n_out):
        log_x = log_mass_min_out + i * dx
        x_j = (log_x - log_data_xmin) / step_size

        if x_j < 0:
            # boundary_condition_const, left edge.
            out[i] = data_native[0]
            continue
        if x_j >= n_data:
            # boundary_condition_const, right edge: interpolate_1d_init()
            # literally reuses the previous OUTPUT sample here, not
            # data_native[-1] -- replicated verbatim (the two agree once the
            # output grid has crossed the native table's edge, since x_j is
            # monotonically increasing in i).
            out[i] = out[i - 1]
            continue

        j = int(x_j)
        if j >= n_data - 1:
            out[i] = data_native[n_data - 1]
        else:
            f = x_j - j
            out[i] = (1.0 - f) * data_native[j] + f * data_native[j + 1]

    return out, dx


def _interp_1d(resampled, log_mass_min_out, dx, log_mass):
    """Direct port of interpolate_1d() (interpolation.h), with
    boundary_condition_const: query the resampled grid at an arbitrary
    log10(mass), clamping to the edge value outside its domain."""
    n = len(resampled)
    i = (log_mass - log_mass_min_out) / dx
    if i < 0:
        return float(resampled[0])
    if i >= n - 1:
        return float(resampled[n - 1])
    idx = int(i)
    frac = i - idx
    return float(resampled[idx] * (1.0 - frac) + resampled[idx + 1] * frac)


class RadiationTable:
    """Reads a GEAR yields_table's Data/Radiation group and reproduces
    radiation_get_luminosities_from_raw() / radiation_get_ionization_rate_
    from_raw() / radiation_get_mean_excess_photon_energy_HI_from_raw()
    (src/feedback/GEAR/radiation.c), in plain CGS.

    @param h5_path Path to the yields_table HDF5 file (e.g. POPIIsw.h5),
    resolved relative to the current working directory -- see
    yields_table_path_from_snapshot().
    @param interpolation_size_mass Number of points in the resampled mass
    grid (GEARFeedback:radiation_interpolation_size_mass; default 200,
    matching radiation.c's own default).
    """

    def __init__(
        self, h5_path, interpolation_size_mass=DEFAULT_INTERPOLATION_SIZE_MASS
    ):
        with h5py.File(h5_path, "r") as h:
            if "Data/Radiation" not in h:
                raise RuntimeError(
                    f"'{h5_path}' has no Data/Radiation group -- it needs "
                    f"regenerating with the current pychem radiation "
                    f"generator (see PLAN_radiation_tables.md Phase 1c)."
                )
            rad = h["Data/Radiation"]
            imf = h["Data/IMF"]

            dimensionality = _decode(rad.attrs["dimensionality"])
            if dimensionality != "M":
                raise NotImplementedError(
                    f"Data/Radiation has dimensionality '{dimensionality}': "
                    f"only a 1D ('M', mass-only, fits-mode) table is "
                    f"supported here, mirroring radiation_check_is_1d()'s "
                    f"abort in radiation.c -- no SubgridRadiation example "
                    f"exercises a 2D ('M,Z') table yet."
                )

            log_data_xmin = float(rad.attrs["m0"])
            step_size = float(rad.attrs["dm"])

            luminosity_native = rad["Luminosity"][:].astype(np.float64)
            q_h_native = rad["Q_H"][:].astype(np.float64)
            dot_e_excess_native = rad["DotEExcess"][:].astype(np.float64)

            # Output grid bounds: the model's own IMF mass range
            # (Data/IMF's Mmin/Mmax), not the table's own native range --
            # matches radiation_build_tables()'s deliberate choice
            # (sm->imf.mass_min/max, not the file's grid extent).
            mass_min = float(imf.attrs["Mmin"])
            mass_max = float(imf.attrs["Mmax"])

        log_mass_min_out = np.log10(mass_min)
        log_mass_max_out = np.log10(mass_max)

        self._log_mass_min_out = log_mass_min_out
        self._luminosity, self._dx = _resample_1d(
            log_mass_min_out,
            log_mass_max_out,
            interpolation_size_mass,
            log_data_xmin,
            step_size,
            luminosity_native,
        )
        self._q_h, _ = _resample_1d(
            log_mass_min_out,
            log_mass_max_out,
            interpolation_size_mass,
            log_data_xmin,
            step_size,
            q_h_native,
        )
        self._dot_e_excess, _ = _resample_1d(
            log_mass_min_out,
            log_mass_max_out,
            interpolation_size_mass,
            log_data_xmin,
            step_size,
            dot_e_excess_native,
        )

    def luminosity_erg_s(self, mass_msun):
        """Bolometric luminosity of a single star of the given mass, cgs
        erg/s (radiation_get_luminosities_from_raw())."""
        log_m = np.log10(mass_msun)
        return _interp_1d(self._luminosity, self._log_mass_min_out, self._dx, log_m)

    def ionizing_photon_rate_s(self, mass_msun):
        """Q_H, the ionizing photon emission rate of a single star of the
        given mass, cgs photons/s
        (radiation_get_ionization_rate_from_raw())."""
        log_m = np.log10(mass_msun)
        return _interp_1d(self._q_h, self._log_mass_min_out, self._dx, log_m)

    def mean_excess_photon_energy_erg(self, mass_msun):
        """Mean excess photon energy above the 13.6 eV HI ionization
        threshold, Q_H-weighted, cgs erg
        (radiation_get_mean_excess_photon_energy_HI_from_raw()): the ratio
        of the interpolated DotEExcess and Q_H tables, not a direct read of
        the file's own MeanExcessPhotonEnergyHI dataset -- radiation.c
        never reads that dataset (see the Phase 2 parity harness log entry
        for the correction that established this).

        @return 0 if this mass produces no ionizing photons (Q_H(mass) <= 0),
        matching the C getter's degenerate-ratio guard.
        """
        log_m = np.log10(mass_msun)
        q_h = _interp_1d(self._q_h, self._log_mass_min_out, self._dx, log_m)
        if q_h <= 0.0:
            return 0.0
        dot_e = _interp_1d(self._dot_e_excess, self._log_mass_min_out, self._dx, log_m)
        return dot_e / q_h
