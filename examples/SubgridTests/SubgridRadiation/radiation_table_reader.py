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
"""Read the Data/Radiation group of a GEAR yields_table HDF5 file.

Shared by the SubgridRadiation examples' ``*_analytic_check.py`` scripts to
look up the ionizing photon rate Q_H(mass), the bolometric luminosity
L(mass), and the mean excess photon energy above 13.6 eV, from the actual
pychem-generated table (e.g. ``POPIIsw.h5``) -- instead of the standalone
piecewise-fit Python ports these scripts used before the table migration.

Replicates the read/interpolate chain ``src/feedback/GEAR/radiation.c``
actually runs, not pychem's own generation-time interpolation convention:

* The table is read exactly like ``radiation_read_grid_metadata()`` +
  ``radiation_build_tables()``: the group-level ``dimensionality``/``m0``/
  ``dm``/``nm`` attributes and the Luminosity/Q_H/DotEExcess datasets, all
  still plain CGS doubles (see ``radiation.c``'s own doxygen for why no
  SWIFT-internal unit conversion or ``RADIATION_DOT_N_ION_TABLE_SCALING`` is
  needed here: this module stays in CGS throughout, matching the
  pre-conversion checker scripts' own convention, and those two factors
  would cancel out again on the way back to CGS in any case).
* Interpolation is SWIFT's own log-log scheme, matching pychem's own
  log10(mass)-vs-log10(value) ``np.interp`` convention: ``interpolation.h``'s
  ``interpolate_1d_init()``/``interpolate_1d()`` are linear interpolators,
  unchanged, fed log10(value) (floored at :data:`LOG_FLOOR`, matching
  ``radiation.c``'s ``RADIATION_LOG_FLOOR_CGS`` bit-for-bit) instead of the
  raw physical value, and ``radiation_get_luminosities_from_raw()`` /
  ``_ionization_rate_from_raw()`` / ``_mean_excess_photon_energy_HI_from_raw()``
  each exponentiate the interpolated result back before returning. Ported
  here as :func:`_safe_log10` before :func:`_resample_1d`, and ``10 **
  value`` after :func:`_interp_1d`.
* The C getters get an exact-0.0 return "for free" in the threshold-
  crossing cell, via IEEE-754 float32 underflow of the floored value
  (SWIFT's internal-unit conversion divides every value down into
  float32's representable range first); that trick does NOT transfer to
  this module (it stays in plain CGS -- Q_H reaches ~1e49 -- so narrowing
  a value that large through float32 overflows to ``inf`` rather than
  underflowing to ``0``). Note pychem's own table-generation-time
  interpolator has no equivalent zero-passthrough at all (a floored entry
  there comes back as its floor value, not an exact 0); this module's
  zero-indicator table (a companion array resampled/interpolated
  alongside each log-valued one) is this module's own mechanism for
  matching the C getters' behaviour, not a port of anything in pychem.
  See :func:`_exp10_with_zero_passthrough`'s own docstring.
* Same two-stage structure as ``radiation_build_tables()``: the native
  table (evenly spaced in log10(mass), ``m0``/``dm``/``nm``) is first
  resampled onto ``radiation_interpolation_size_mass`` points evenly spaced
  in log10(mass) between the table's own Data/IMF ``Mmin``/``Mmax`` (NOT
  the native table's own mass range -- this is a deliberate SWIFT design
  choice, see ``radiation_build_tables()``'s doxygen), then a query at an
  arbitrary mass linearly interpolates that resampled grid. Both stages
  clamp to the edge value outside their domain (``boundary_condition_const``)
  -- never extrapolate, never raise.

Only 1D (``"M"``, fits-mode) tables are supported: every SubgridRadiation
example this module is used by loads a fits-mode table. A 2D (``"M,Z"``)
table raises, mirroring ``radiation_check_is_1d()``'s abort in the C code.

The table path is never hardcoded: :func:`yields_table_path_from_snapshot`
reads ``GEARFeedback:yields_table`` from the run's own Parameters group (the
resolved value SWIFT actually used, exactly as parsed from ``params.yml``,
including any command-line override) -- the same, deliberately
uncanonicalized, cwd-relative path SWIFT's own ``H5Fopen`` call resolves it
against; this module does not attempt to fix that, only to read the file the
run itself read.
"""

from __future__ import annotations

import os
from typing import Sequence

import h5py
import numpy as np

#: Default for GEARFeedback:radiation_interpolation_size_mass
#: (src/feedback/GEAR/radiation.c, examples/parameter_example.yml).
DEFAULT_INTERPOLATION_SIZE_MASS = 200

#: Floor applied before taking log10() of a native Data/Radiation value, so
#: a genuinely-zero table entry (Q_H/DotEExcess below the source table's
#: own ionization threshold) does not produce log10(0) = -inf. Matches
#: pychem's own floor bit-for-bit (PyChemInitTable/libparsec_radiation.py's
#: `_LOG_FLOOR`) and src/feedback/GEAR/radiation.c's
#: `RADIATION_LOG_FLOOR_CGS`.
LOG_FLOOR = 1e-300


def _safe_log10(data_native: np.ndarray) -> np.ndarray:
    """Compute log10 of a native Data/Radiation array with pychem's floor.

    Matches ``radiation_read_cgs_array()``'s ``log_data_internal``
    computation (radiation.c) -- this module stays in CGS throughout (see
    the module docstring), so there is no unit-conversion term to subtract
    here, unlike the C function it mirrors.

    Parameters
    ----------
    data_native : np.ndarray
        Native CGS values (Luminosity, Q_H, or DotEExcess), non-negative.

    Returns
    -------
    np.ndarray
        ``log10(clip(data_native, LOG_FLOOR, None))``.
    """
    return np.log10(np.clip(data_native, LOG_FLOOR, None))


def _native_zero_mask(data_native_cgs: np.ndarray) -> np.ndarray:
    """Build a 0/1 indicator of exactly-zero native CGS table entries.

    1.0 where a native Data/Radiation CGS entry is exactly 0 (Q_H/
    DotEExcess below the source table's own ionization threshold), else
    0.0. Resampled/interpolated through the same
    :func:`_resample_1d`/:func:`_interp_1d` machinery as the log-valued
    tables themselves -- see :func:`_exp10_with_zero_passthrough`.

    Parameters
    ----------
    data_native_cgs : np.ndarray
        Native CGS values (Q_H or DotEExcess).

    Returns
    -------
    np.ndarray
        1.0/0.0 indicator array, float64, same length as `data_native_cgs`.
    """
    return (data_native_cgs <= 0.0).astype(np.float64)


def _exp10_with_zero_passthrough(log_value: float, zero_indicator: float) -> float:
    """Exponentiate a log10(value), forcing exact 0.0 near any native zero.

    Forces an exact 0.0 return whenever the query's interpolation bracket
    touches at least one native-zero table entry, not only when every
    contributing point was zero: the source data jumps by hundreds of
    decades across the ionization threshold, so a bracket straddling it
    interpolates two physically unrelated log10 values and produces
    garbage on exponentiation, not merely an inaccurate number. This
    module stays in plain CGS throughout (Q_H up to ~1e49), so it cannot
    rely on radiation.c's float32-underflow shortcut the way that C code
    does for a query entirely inside an all-zero run; the C shortcut does
    not help in the threshold-crossing cell either, which is why this
    function's own explicit check exists.

    `zero_indicator` is :func:`_interp_1d` applied to the SAME resampled
    grid as `log_value`, but built from :func:`_native_zero_mask` instead
    of a log-valued table: linearly interpolating a strict {0.0, 1.0}
    array (through both :func:`_resample_1d`'s build-time resample and
    this function's own query) is exactly 0.0 only when every native
    point the interpolation drew from was itself nonzero (IEEE-754
    double: a weighted sum of zeros stays exactly 0.0), so ``> 0.0`` is an
    exact check for "at least one contributing point was zero", not a
    fuzzy threshold.

    Parameters
    ----------
    log_value : float
        log10(value), as returned by :func:`_interp_1d` on a log-valued
        table.
    zero_indicator : float
        :func:`_interp_1d` on the matching zero-mask table (see
        :func:`_native_zero_mask`); any value greater than 0.0 forces an
        exact 0.0 return.

    Returns
    -------
    float
        0.0 if `zero_indicator` > 0.0, else ``10 ** log_value``.
    """
    if zero_indicator > 0.0:
        return 0.0
    return 10.0**log_value


def _decode(raw: object) -> object:
    """Decode an HDF5 attribute value that may come back as bytes.

    Parameters
    ----------
    raw : object
        The raw value returned by h5py for an attribute.

    Returns
    -------
    object
        `raw` decoded to `str` if it was `bytes`/`np.bytes_` (directly, or
        wrapped in a 0-d `np.ndarray`); otherwise `raw` unchanged.
    """
    if isinstance(raw, (bytes, np.bytes_)):
        return raw.decode()
    if isinstance(raw, np.ndarray):
        raw = raw.item()
        return raw.decode() if isinstance(raw, (bytes, np.bytes_)) else raw
    return raw


def yields_table_path_from_snapshot(snapshot_file: str) -> str:
    """Return the raw yields_table path a run actually used.

    Reads ``GEARFeedback:yields_table`` from the snapshot's own Parameters
    group -- exactly as parsed from ``params.yml`` (relative or absolute,
    whatever the run was given), not yet resolved to an existing file; see
    :func:`resolve_yields_table_path`.

    Parameters
    ----------
    snapshot_file : str
        One of the run's snapshot files (any one holds the same Parameters
        group).

    Returns
    -------
    str
        The raw, unresolved ``GEARFeedback:yields_table`` string.
    """
    with h5py.File(snapshot_file, "r") as h:
        return _decode(h["Parameters"].attrs["GEARFeedback:yields_table"])


def resolve_yields_table_path(
    raw_path: str, snapshot_file: str, script_dir: str | None = None
) -> str:
    """Resolve a raw, uncanonicalized yields_table path to an existing file.

    SWIFT stores ``GEARFeedback:yields_table`` exactly as given in
    ``params.yml``, relative to whatever working directory the run itself
    used; this resolves that same string against a handful of plausible
    locations without hardcoding any one of them.

    Tried in order: as given, relative to the current working directory
    (matches SWIFT's own H5Fopen resolution, the common case when this
    check script runs right after run.sh in the same directory); relative
    to the snapshot file's own directory and that directory's parent
    (covers a snapshot glob pointed at an archived/moved run directory,
    e.g. run.sh's `mv snap $run_name`, where the table file was left
    behind one level up rather than moved with the snapshots); relative to
    `script_dir`, the calling check script's own directory (covers the
    examples with a table checked in next to the script itself).

    Parameters
    ----------
    raw_path : str
        The string from :func:`yields_table_path_from_snapshot`.
    snapshot_file : str
        One of the run's snapshot files (for the snapshot-directory
        candidates).
    script_dir : str or None, optional
        Directory of the calling check script, or None to skip that
        candidate.

    Returns
    -------
    str
        The first candidate that exists.

    Raises
    ------
    FileNotFoundError
        If none of the candidate paths exists; names every path tried.
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
    files: Sequence[str],
    script_dir: str | None = None,
    default_interpolation_size_mass: int = DEFAULT_INTERPOLATION_SIZE_MASS,
) -> "RadiationTable":
    """Open the RadiationTable for a run, from its own recorded table path.

    Reads ``GEARFeedback:yields_table`` from the run's own Parameters group
    -- never a hardcoded filename.

    Parameters
    ----------
    files : Sequence[str]
        The run's sorted snapshot filenames (any one holds the Parameters
        group this reads).
    script_dir : str or None, optional
        Directory of the calling check script (see
        :func:`resolve_yields_table_path`), or None to skip that fallback
        candidate.
    default_interpolation_size_mass : int, optional
        Used when the snapshot predates
        ``GEARFeedback:radiation_interpolation_size_mass`` (see
        :func:`interpolation_size_mass_from_snapshot`).

    Returns
    -------
    RadiationTable
        The opened table, resampled onto this run's own interpolation grid.
    """
    raw_path = yields_table_path_from_snapshot(files[0])
    resolved_path = resolve_yields_table_path(raw_path, files[0], script_dir)
    interpolation_size_mass = interpolation_size_mass_from_snapshot(
        files[0], default_interpolation_size_mass
    )
    return RadiationTable(resolved_path, interpolation_size_mass)


def interpolation_size_mass_from_snapshot(
    snapshot_file: str, default: int = DEFAULT_INTERPOLATION_SIZE_MASS
) -> int:
    """Return the run's radiation_interpolation_size_mass, or a default.

    Parameters
    ----------
    snapshot_file : str
        One of the run's snapshot files.
    default : int, optional
        Value to return if the snapshot predates
        ``GEARFeedback:radiation_interpolation_size_mass`` (e.g. an old
        archived run made before the table migration).

    Returns
    -------
    int
        The run's ``GEARFeedback:radiation_interpolation_size_mass``, or
        `default`.
    """
    with h5py.File(snapshot_file, "r") as h:
        raw = h["Parameters"].attrs.get(
            "GEARFeedback:radiation_interpolation_size_mass"
        )
    return int(_decode(raw)) if raw is not None else default


def _resample_1d(
    log_mass_min_out: float,
    log_mass_max_out: float,
    n_out: int,
    log_data_xmin: float,
    step_size: float,
    data_native: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Resample a native log-mass-uniform array onto an output grid.

    Direct port of ``interpolate_1d_init()`` (interpolation.h), with
    ``boundary_condition_const``: resample `data_native` (evenly spaced in
    log10(mass), starting at `log_data_xmin` with step `step_size`) onto
    `n_out` points evenly spaced in log10(mass) in
    ``[log_mass_min_out, log_mass_max_out]``.

    Parameters
    ----------
    log_mass_min_out : float
        Lower bound of the output grid, in log10(mass).
    log_mass_max_out : float
        Upper bound of the output grid, in log10(mass).
    n_out : int
        Number of points in the output grid.
    log_data_xmin : float
        log10(mass) of `data_native`'s first sample.
    step_size : float
        Spacing, in log10(mass), between consecutive `data_native` samples.
    data_native : np.ndarray
        The native table values, evenly spaced in log10(mass).

    Returns
    -------
    out : np.ndarray
        The resampled array, of length `n_out`.
    dx : float
        Spacing, in log10(mass), between consecutive `out` samples.
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


def _interp_1d(
    resampled: np.ndarray, log_mass_min_out: float, dx: float, log_mass: float
) -> float:
    """Query a resampled grid at an arbitrary log10(mass).

    Direct port of ``interpolate_1d()`` (interpolation.h), with
    ``boundary_condition_const``: clamps to the edge value outside the
    grid's domain.

    Parameters
    ----------
    resampled : np.ndarray
        The output of :func:`_resample_1d`.
    log_mass_min_out : float
        Lower bound of `resampled`'s domain, in log10(mass).
    dx : float
        Spacing, in log10(mass), between consecutive `resampled` samples.
    log_mass : float
        log10(mass) to query.

    Returns
    -------
    float
        The (possibly edge-clamped) linearly interpolated value.
    """
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
    """Read a GEAR yields_table's Data/Radiation group and interpolate it.

    Reproduces ``radiation_get_luminosities_from_raw()`` /
    ``radiation_get_ionization_rate_from_raw()`` /
    ``radiation_get_mean_excess_photon_energy_HI_from_raw()``
    (src/feedback/GEAR/radiation.c), in plain CGS.

    Parameters
    ----------
    h5_path : str
        Path to the yields_table HDF5 file (e.g. ``POPIIsw.h5``), resolved
        relative to the current working directory -- see
        :func:`yields_table_path_from_snapshot`.
    interpolation_size_mass : int, optional
        Number of points in the resampled mass grid
        (``GEARFeedback:radiation_interpolation_size_mass``; default 200,
        matching radiation.c's own default).
    """

    def __init__(
        self,
        h5_path: str,
        interpolation_size_mass: int = DEFAULT_INTERPOLATION_SIZE_MASS,
    ) -> None:
        with h5py.File(h5_path, "r") as h:
            if "Data/Radiation" not in h:
                raise RuntimeError(
                    f"'{h5_path}' has no Data/Radiation group -- it needs "
                    f"regenerating with the current pychem radiation table "
                    f"generator (run pychem's yields-table generator on the "
                    f"source chimie parameters, then point "
                    f"GEARFeedback:yields_table at the regenerated file)."
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

            luminosity_raw = rad["Luminosity"][:].astype(np.float64)
            q_h_raw = rad["Q_H"][:].astype(np.float64)
            dot_e_excess_raw = rad["DotEExcess"][:].astype(np.float64)

            # log10(value), pychem-floored (_safe_log10()): the tables below
            # are resampled/interpolated in log-value space, matching
            # radiation_build_tables()'s own log-log storage. Luminosity is
            # never 0 in this table (unlike Q_H/DotEExcess), so it gets no
            # zero-mask companion -- see _exp10_with_zero_passthrough().
            luminosity_native = _safe_log10(luminosity_raw)
            q_h_native = _safe_log10(q_h_raw)
            dot_e_excess_native = _safe_log10(dot_e_excess_raw)
            q_h_zero_native = _native_zero_mask(q_h_raw)
            dot_e_excess_zero_native = _native_zero_mask(dot_e_excess_raw)

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
        self._q_h_zero, _ = _resample_1d(
            log_mass_min_out,
            log_mass_max_out,
            interpolation_size_mass,
            log_data_xmin,
            step_size,
            q_h_zero_native,
        )
        self._dot_e_excess_zero, _ = _resample_1d(
            log_mass_min_out,
            log_mass_max_out,
            interpolation_size_mass,
            log_data_xmin,
            step_size,
            dot_e_excess_zero_native,
        )

    def luminosity_erg_s(self, mass_msun: float) -> float:
        """Return the bolometric luminosity of a single star, cgs erg/s.

        Mirrors ``radiation_get_luminosities_from_raw()``:
        ``exp10(interpolate_1d(log10-valued table))``. No zero-passthrough
        is applied: Luminosity is never 0 in this table.

        Parameters
        ----------
        mass_msun : float
            Stellar mass, in Msun.

        Returns
        -------
        float
            Bolometric luminosity, cgs erg/s.
        """
        log_m = np.log10(mass_msun)
        log_lum = _interp_1d(self._luminosity, self._log_mass_min_out, self._dx, log_m)
        return 10.0**log_lum

    def ionizing_photon_rate_s(self, mass_msun: float) -> float:
        """Return Q_H, the ionizing photon emission rate, cgs photons/s.

        Mirrors ``radiation_get_ionization_rate_from_raw()``:
        ``exp10(interpolate_1d(log10-valued table))``, with an exact-0.0
        passthrough near the ionization threshold -- see
        :func:`_exp10_with_zero_passthrough`.

        Parameters
        ----------
        mass_msun : float
            Stellar mass, in Msun.

        Returns
        -------
        float
            Ionizing photon emission rate, cgs photons/s.
        """
        log_m = np.log10(mass_msun)
        log_q_h = _interp_1d(self._q_h, self._log_mass_min_out, self._dx, log_m)
        q_h_zero = _interp_1d(self._q_h_zero, self._log_mass_min_out, self._dx, log_m)
        return _exp10_with_zero_passthrough(log_q_h, q_h_zero)

    def mean_excess_photon_energy_erg(self, mass_msun: float) -> float:
        """Return the Q_H-weighted mean excess photon energy, cgs erg.

        Mean excess photon energy above the 13.6 eV HI ionization
        threshold, for a single star of the given mass: the ratio of the
        interpolated, exponentiated DotEExcess and Q_H tables, not a direct
        read of the file's own MeanExcessPhotonEnergyHI dataset --
        radiation.c never reads that dataset, computing the ratio directly
        instead (mirrors
        ``radiation_get_mean_excess_photon_energy_HI_from_raw()``).

        Parameters
        ----------
        mass_msun : float
            Stellar mass, in Msun.

        Returns
        -------
        float
            Mean excess photon energy, cgs erg, or 0 if this mass produces
            no ionizing photons (Q_H(mass) <= 0, via
            :func:`ionizing_photon_rate_s`'s own zero-passthrough),
            matching the C getter's degenerate-ratio guard.
        """
        q_h = self.ionizing_photon_rate_s(mass_msun)
        if q_h <= 0.0:
            return 0.0
        log_m = np.log10(mass_msun)
        log_dot_e = _interp_1d(
            self._dot_e_excess, self._log_mass_min_out, self._dx, log_m
        )
        dot_e_zero = _interp_1d(
            self._dot_e_excess_zero, self._log_mass_min_out, self._dx, log_m
        )
        dot_e = _exp10_with_zero_passthrough(log_dot_e, dot_e_zero)
        return dot_e / q_h
