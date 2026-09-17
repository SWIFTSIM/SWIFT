#!/usr/bin/env python3
"""
SPH-kernel-smoothed maps of any PartType0 dataset in a HomogeneousBox run.

Generalises plot_lw_fuv_maps_smoothed2.py's four hardcoded panels (FUV/LW
specific energy, dust temperature, H2I mass fraction) to any gas dataset,
via swiftsimio.visualisation.projection: a thick z-slab, projected with
each particle's own SPH smoothing length, either as a mass-weighted
average (default) or as a projected sum (--projection sum, e.g. a column
quantity).

A field is given as NAME or NAME:COMPONENT. COMPONENT selects a vector's
norm (":norm") or component index (":0", ":1", ...), or a named-columns
species/element array's column, by name or index (e.g.
"MetalMassFractions:fe" or "MetalMassFractions:0").

Two modes:

* Single snapshot (default): the newest complete snapshot matching
  --snapshots is used; each field's colour scale is data-driven from that
  frame.
* --time-series: every complete snapshot matched by --snapshots is used.
  One frozen colour scale per field (from the global min/max across the
  series) and one frozen slab (z0, half-thickness) are used for every
  frame, so a given pixel maps to the same physical location and colour
  meaning throughout. A combined summary panel (one small panel per
  snapshot, shared colour scale) is written per field to <outdir>; add
  --frames to also write each snapshot's individual panel there.

Stars (PartType4) and sinks (PartType3) inside the projected slab are
overlaid on every panel by default (--no-stars/--no-sinks to turn either
off); --star-age-max restricts the stars shown to those younger than a
given age in Myr, when the snapshot has a BirthTimes field.

With no --fields, the four original panels are reproduced.
"""

import argparse
import glob
import os
import re
import sys
from dataclasses import dataclass
from typing import Optional

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import swiftsimio as sw
from matplotlib import patheffects
from matplotlib.colors import LogNorm, Normalize
from swiftsimio.objects import cosmo_array
from swiftsimio.visualisation.projection import project_gas

RESOLUTION_DEFAULT = 512
#: Auto scale heuristic: log only if every finite value is positive and the
#: span exceeds this many decades.
LOG_DEX_THRESHOLD = 2

#: --fields default: the four panels of the original hardcoded script.
DEFAULT_FIELDS = ("FUVSpecificEnergies", "LWSpecificEnergies", "DustTemperature", "H2I")

#: Per-field style, matching the original hardcoded script exactly, keyed
#: by the bare HDF5 dataset name (a field's component, if any, does not
#: change its style). Any field not listed here gets fully auto-derived
#: style (see resolve_style).
FIELD_PRESETS = {
    "FUVSpecificEnergies": dict(
        scale="log",
        cmap="inferno",
        floor=1e-2,
        clip_negative=True,
        label="mass-weighted FUVSpecificEnergies [internal units]",
        title="FUV-band (6-11.2 eV) specific energy (SPH-smoothed)",
    ),
    "LWSpecificEnergies": dict(
        scale="log",
        cmap="viridis",
        floor=1e-2,
        clip_negative=True,
        label="mass-weighted LWSpecificEnergies [internal units]",
        title="Lyman-Werner-band (11.2-13.6 eV) specific energy (SPH-smoothed)",
    ),
    "DustTemperature": dict(
        scale="linear",
        cmap="magma",
        floor=None,
        clip_negative=False,
        label="mass-weighted DustTemperature [K]",
        title="Dust temperature (Grackle diagnostic, SPH-smoothed)",
    ),
    "H2I": dict(
        scale="log",
        cmap="cividis",
        floor=1e-12,
        clip_negative=False,
        label="mass-weighted H2I mass fraction",
        title="H2I mass fraction (Grackle, SPH-smoothed)",
    ),
}


@dataclass
class FieldSpec:
    """One parsed --fields token."""

    raw: str
    name: str
    component: Optional[str]

    @property
    def key(self) -> str:
        """Filesystem/dict-safe identifier derived from `raw`."""
        return re.sub(r"[^0-9A-Za-z_]", "_", self.raw)


def parse_field_spec(token: str) -> FieldSpec:
    """Parse one --fields token into a FieldSpec.

    Parameters
    ----------
    token : str
        ``NAME`` or ``NAME:COMPONENT`` (see the module docstring).

    Returns
    -------
    FieldSpec
    """
    if ":" in token:
        name, component = token.split(":", 1)
    else:
        name, component = token, None
    return FieldSpec(raw=token, name=name, component=component)


def parse_kv_list(items: list, value_type=str) -> dict:
    """Parse a list of "KEY=VALUE" strings into a dict.

    Parameters
    ----------
    items : list of str
        Each element must contain "=" once, splitting a field spec's
        `raw` token from an override value.
    value_type : callable, optional
        Applied to the value half of each pair (default: keep as str).

    Returns
    -------
    dict
    """
    result = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"Expected KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        result[key] = value_type(value)
    return result


def available_field_names(gas) -> list:
    """List every PartType0 dataset this snapshot exposes, for an error message.

    A named-columns field (e.g. MetalMassFractions) is listed once with
    its column names, since a bare column index would otherwise be
    meaningless without them.
    """
    gm = gas.group_metadata
    names = []
    for path in gm.field_paths:
        if not path.startswith("PartType0/"):
            continue
        bare = path[len("PartType0/") :]
        cols = gm.named_columns.get(path)
        names.append(f"{bare} (columns: {', '.join(cols)})" if cols else bare)
    return sorted(names)


def resolve_field_array(data, spec: FieldSpec) -> np.ndarray:
    """Read one field spec's per-particle gas values.

    Reads a scalar or vector field through swiftsimio's own attribute
    (whatever HDF5 dataset it is, including a GEAR-custom one swiftsimio
    has no built-in schema entry for). A named-columns species/element
    array (e.g. MetalMassFractions) is read directly from the underlying
    HDF5 file instead, since swiftsimio wraps that case in a per-species
    accessor object rather than a plain array.

    Parameters
    ----------
    data : swiftsimio SWIFTDataset
        The loaded snapshot.
    spec : FieldSpec
        The field (and optional component) to read.

    Returns
    -------
    np.ndarray or unyt_array
        One value per gas particle.
    """
    gas = data.gas
    gm = gas.group_metadata
    full_path = f"PartType0/{spec.name}"
    if full_path not in gm.field_paths:
        raise SystemExit(
            f"Unknown field {spec.name!r}. Available PartType0 fields:\n  "
            + "\n  ".join(available_field_names(gas))
        )

    cols = gm.named_columns.get(full_path)
    if cols:
        with h5py.File(gas.filename, "r") as f:
            raw2d = f["PartType0"][spec.name][:]
        if spec.component is None:
            raise SystemExit(
                f"{spec.name!r} has {len(cols)} columns; select one, e.g. "
                f"'{spec.name}:{cols[0]}' or '{spec.name}:0'. Columns: {cols}"
            )
        if spec.component in cols:
            idx = cols.index(spec.component)
        else:
            try:
                idx = int(spec.component)
            except ValueError:
                raise SystemExit(
                    f"{spec.name!r} has no column {spec.component!r}. Columns: {cols}"
                )
        return raw2d[:, idx].astype(np.float64)

    attr = dict(zip(gm.field_paths, gm.field_names))[full_path]
    arr = getattr(gas, attr)
    if arr.ndim == 1:
        if spec.component is not None:
            raise SystemExit(f"{spec.name!r} is a scalar field; it takes no component.")
        return arr
    if arr.ndim == 2:
        if spec.component == "norm":
            return np.sqrt((arr**2).sum(axis=1))
        if spec.component is not None:
            try:
                idx = int(spec.component)
            except ValueError:
                raise SystemExit(
                    f"{spec.name!r} is a vector field; component must be 'norm' or "
                    f"a column index 0..{arr.shape[1] - 1}, not {spec.component!r}."
                )
            return arr[:, idx]
        raise SystemExit(
            f"{spec.name!r} is a vector field ({arr.shape[1]} components); select "
            f"'{spec.name}:norm' or '{spec.name}:0'."
        )
    raise SystemExit(f"{spec.name!r} has unsupported shape {arr.shape}.")


def project_field(
    data,
    field_array,
    mass_map,
    region,
    resolution: int,
    projection: str,
    clip_negative: bool,
    key: str,
) -> tuple:
    """Project one field's slab values, mass-weighted-average or summed.

    Parameters
    ----------
    data : swiftsimio SWIFTDataset
    field_array : np.ndarray or unyt_array
        Per-particle values from `resolve_field_array`.
    mass_map : unyt_array
        The slab's projected mass map (from project_gas(project="masses")).
    region : cosmo_array
        The projection region (see `make_region`).
    resolution : int
    projection : str
        "average" (mass-weighted mean) or "sum" (raw projected sum).
    clip_negative : bool
        If True, negative values are set to 0 before projecting, so they
        cannot pull an average negative (and off a log scale); the
        clipped particle count/mass fraction is returned for reporting.
    key : str
        Unique per-field identifier for the temporary gas attribute.

    Returns
    -------
    tuple
        (avg, n_clipped, mass_fraction_clipped): the projected 2D array as
        plain floats, and the negative-clip diagnostics (0, 0.0 if
        `clip_negative` is False).
    """
    weight = field_array.copy()
    n_clipped, mfrac_clipped = 0, 0.0
    if clip_negative:
        neg = weight < 0
        n_clipped = int(np.sum(neg))
        if n_clipped:
            mass = data.gas.masses
            mfrac_clipped = float(np.sum(mass[neg]) / np.sum(mass))
        weight[neg] = 0.0

    tmp_name = f"_tmp_{key}_w"
    if projection == "average":
        setattr(data.gas, tmp_name, weight * data.gas.masses)
        proj = project_gas(data, resolution=resolution, project=tmp_name, region=region)
        avg = proj / mass_map
    else:
        setattr(data.gas, tmp_name, weight)
        avg = project_gas(data, resolution=resolution, project=tmp_name, region=region)

    return np.asarray(getattr(avg, "value", avg)), n_clipped, mfrac_clipped


def _snapshot_time_if_complete(path: str, required_names: list):
    """Return the snapshot's Header Time if it holds every required field, else None."""
    try:
        with h5py.File(path, "r") as f:
            gas = f["PartType0"]
            if all(k in gas for k in required_names) and "Header" in f:
                return float(np.atleast_1d(f["Header"].attrs["Time"])[0])
    except (OSError, KeyError):
        pass
    return None


def find_all_snapshots(patterns: list, required_names: list) -> list:
    """Expand one or more globs/paths to every complete snapshot, sorted by time.

    A single-pattern, single-target caller just takes the last element of
    the returned list (the newest one).

    Parameters
    ----------
    patterns : list of str
        Snapshot globs/paths.
    required_names : list of str
        Bare HDF5 dataset names that must be present for a snapshot to
        count as complete (the requested fields' names, without their
        component suffix).
    """
    seen = set()
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if not matches:
            print(f"warning: {pattern!r} matched no files", file=sys.stderr)
        for m in matches:
            if m not in seen:
                seen.add(m)
                paths.append(m)
    if not paths:
        raise SystemExit(f"No snapshot matches {patterns!r}")

    timed = []
    for path in paths:
        t = _snapshot_time_if_complete(path, required_names)
        if t is not None:
            timed.append((t, path))
    if not timed:
        with h5py.File(paths[-1], "r") as f:
            available = sorted(f["PartType0"].keys()) if "PartType0" in f else []
        missing = [n for n in required_names if n not in available]
        raise SystemExit(
            f"No snapshot has every requested field {required_names!r}. "
            f"Field(s) not found in {paths[-1]}: {missing!r}. Available "
            f"PartType0 datasets there:\n  " + "\n  ".join(available)
        )
    timed.sort(key=lambda pair: (pair[0], pair[1]))
    return [path for _, path in timed]


def pick_slab(data, boxsize: float, min_count: int = 3000):
    """Centre a z-slab on the brightest star and widen it until it is populated.

    Falls back to a box-centred slab if this snapshot has no stars, or no
    FUV-luminosity field on its stars (a config without radiation
    feedback), rather than assuming both are present.
    """
    z0 = 0.5 * boxsize
    try:
        stars = data.stars
        if stars.fuvluminosities.size:
            z0 = float(stars.coordinates[np.argmax(stars.fuvluminosities), 2].value)
    except AttributeError:
        pass
    z = data.gas.coordinates[:, 2].value
    for frac in (0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 1.0):
        half = frac * boxsize
        sel = np.abs(z - z0) < half
        if sel.sum() >= min_count:
            return z0, half, int(sel.sum())
    return z0, boxsize, int(z.size)


def freeze_slab(paths: list, boxsize: float) -> tuple:
    """Pick one (z0, half) slab shared by the whole time series.

    Without this, a pixel at index (i, j) would map to a different
    physical (x, y, z) region in each frame, since pick_slab's z0/half
    otherwise track whichever star is currently brightest. Scans
    backwards from the latest snapshot for the first with a star; falls
    back to the latest snapshot's own box-centred pick otherwise.
    """
    chosen = None
    for path in reversed(paths):
        data = sw.load(path)
        try:
            if data.stars.fuvluminosities.size:
                chosen = data
                break
        except AttributeError:
            continue
    if chosen is None:
        chosen = sw.load(paths[-1])
    z0, half, _ = pick_slab(chosen, boxsize)
    return z0, half


def make_region(boxsize_arr, L: float, z0: float, half: float) -> cosmo_array:
    """Build the (x, y, z) projection region covering the full box in x/y and
    the thick slab in z, tagged with the same comoving/cosmo_factor as the
    snapshot's own BoxSize so swiftsimio's unit machinery accepts it."""
    return cosmo_array(
        [0.0, L, 0.0, L, z0 - half, z0 + half],
        units=boxsize_arr.units,
        comoving=boxsize_arr.comoving,
        cosmo_factor=boxsize_arr.cosmo_factor,
    )


def safe_log_norm(vmin: float, vmax: float) -> LogNorm:
    """LogNorm(vmin, vmax), widening a degenerate range so LogNorm never
    receives vmax <= vmin (e.g. an early snapshot with no signal yet)."""
    if vmin is None or not np.isfinite(vmin) or vmin <= 0:
        vmin = 1e-30
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = vmin * 10.0
    return LogNorm(vmin=vmin, vmax=vmax)


def safe_linear_norm(vmin: float, vmax: float) -> Normalize:
    """Normalize(vmin, vmax), widening a degenerate range the same way as safe_log_norm."""
    if vmin is None or not np.isfinite(vmin):
        vmin = 0.0
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = vmin + 1.0
    return Normalize(vmin=vmin, vmax=vmax)


#: A black stroke around every marker keeps it visible regardless of which
#: colormap (or where on it) the marker happens to sit, since no single
#: marker colour stays legible across every colormap's bright end.
_MARKER_STROKE = [patheffects.withStroke(linewidth=1.8, foreground="black")]


def overlay_sources(
    ax,
    stars_xyz: np.ndarray,
    sinks_xyz: np.ndarray,
    z0: float,
    half: float,
    legend: bool = True,
) -> None:
    """Mark the stars and sinks that lie inside the slab.

    `stars_xyz`/`sinks_xyz` are expected pre-filtered by the caller (e.g.
    an empty array when --no-stars/--no-sinks was requested, or a
    --star-age-max cut already applied), so this function only handles
    the slab selection and drawing.
    """
    any_plotted = False
    for arr, marker, colour, label in (
        (stars_xyz, "*", "white", "stars"),
        (sinks_xyz, "o", "cyan", "sinks"),
    ):
        if arr.shape[0] == 0:
            continue
        sel = np.abs(arr[:, 2] - z0) < half
        if sel.sum() == 0:
            continue
        ax.scatter(
            arr[sel, 0],
            arr[sel, 1],
            marker=marker,
            s=45 if marker == "*" else 18,
            facecolors="none" if marker == "o" else colour,
            edgecolors="black" if marker == "*" else colour,
            linewidths=0.8,
            path_effects=_MARKER_STROKE,
            label=f"{label} ({sel.sum()})",
        )
        any_plotted = True
    if legend and any_plotted:
        ax.legend(loc="upper right", fontsize=7, framealpha=0.6)


#: Set on the first snapshot missing BirthTimes, so --time-series with
#: --star-age-max warns once per run rather than once per frame.
_warned_no_birth_times = False


def select_sources(
    data,
    show_stars: bool,
    show_sinks: bool,
    star_age_max_myr: Optional[float],
) -> tuple:
    """Read stars'/sinks' coordinates for overlay, applying CLI display filters.

    Missing PartType3/PartType4 data is handled the same way as absent
    fields elsewhere in this script: an empty array, not an error. A
    `star_age_max_myr` cut is applied only when the snapshot exposes a
    plain `BirthTimes` field; a birth-scale-factor-only snapshot is not a
    simple age computation, so the cut is skipped with a warning (printed
    once per run, not once per frame) rather than guessed at.

    Parameters
    ----------
    data : swiftsimio SWIFTDataset
    show_stars, show_sinks : bool
        --stars/--sinks (or their --no- counterparts).
    star_age_max_myr : float or None
        --star-age-max; None disables the age cut.

    Returns
    -------
    tuple
        (stars_xyz, sinks_xyz), each an (N, 3) array (possibly empty).
    """
    stars_xyz = np.zeros((0, 3))
    if show_stars and hasattr(data, "stars"):
        stars = data.stars
        stars_xyz = np.asarray(stars.coordinates.value)
        if star_age_max_myr is not None and stars_xyz.shape[0]:
            if hasattr(stars, "birth_times"):
                myr_per_internal = float(
                    (1.0 * stars.birth_times.units).to("Myr").value
                )
                age_myr = (
                    float(data.metadata.time.value)
                    - np.asarray(stars.birth_times.value)
                ) * myr_per_internal
                stars_xyz = stars_xyz[age_myr <= star_age_max_myr]
            else:
                global _warned_no_birth_times
                if not _warned_no_birth_times:
                    print(
                        "warning: --star-age-max requested but at least one "
                        "snapshot has no BirthTimes field; showing all stars "
                        "for those snapshots",
                        file=sys.stderr,
                    )
                    _warned_no_birth_times = True

    sinks_xyz = np.zeros((0, 3))
    if show_sinks and hasattr(data, "sinks"):
        sinks_xyz = np.asarray(data.sinks.coordinates.value)

    return stars_xyz, sinks_xyz


def combined_range(entries: list, field_raw: str) -> tuple:
    """Global (min, max, all_positive) of one field's projected maps across frames.

    Parameters
    ----------
    entries : list of dict
        Frames from `compute_frame`.
    field_raw : str
        A FieldSpec's `raw` token, used as the maps dict key.

    Returns
    -------
    tuple
        (vmin, vmax, all_positive); (None, None, False) if every frame's
        map is empty/non-finite.
    """
    mins, maxs, all_positive = [], [], True
    for entry in entries:
        arr = entry["maps"][field_raw]
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            continue
        if np.any(finite <= 0):
            all_positive = False
        mins.append(float(finite.min()))
        maxs.append(float(finite.max()))
    if not mins:
        return None, None, False
    return min(mins), max(maxs), all_positive


def resolve_style(
    spec: FieldSpec,
    vmin,
    vmax,
    all_positive: bool,
    projection: str,
    overrides: dict,
) -> dict:
    """Resolve one field's scale/colormap/floor/label/title.

    A field in `FIELD_PRESETS` (keyed by its bare name) uses that preset
    unless a CLI override matches its `raw` token; any other field gets a
    fully automatic style: log scale only if every value across the
    frame(s) considered is positive and spans more than
    `LOG_DEX_THRESHOLD` decades, else linear.

    Parameters
    ----------
    spec : FieldSpec
    vmin, vmax : float or None
        The value range to style against (one frame in single-snapshot
        mode, the global range in --time-series mode).
    all_positive : bool
        Whether every finite value in that range is > 0.
    projection : str
        "average" or "sum", used only for the auto-generated label.
    overrides : dict
        CLI overrides, keyed "scale"/"cmap"/"floor"/"label" -> {raw: value}.

    Returns
    -------
    dict
        scale, cmap, floor, label, title.
    """
    preset = FIELD_PRESETS.get(spec.name, {})
    auto_scale = (
        "log"
        if (
            all_positive
            and vmin is not None
            and vmin > 0
            and vmax / vmin > 10**LOG_DEX_THRESHOLD
        )
        else "linear"
    )
    scale = overrides["scale"].get(spec.raw, preset.get("scale", auto_scale))

    cmap = overrides["cmap"].get(spec.raw, preset.get("cmap", "viridis"))

    default_floor = preset.get("floor")
    if default_floor is None and scale == "log" and vmin is not None and vmin > 0:
        default_floor = vmin
    floor = overrides["floor"].get(spec.raw, default_floor)

    weighting_desc = (
        "mass-weighted average of" if projection == "average" else "projected sum of"
    )
    default_label = preset.get("label", f"{weighting_desc} {spec.raw}")
    label = overrides["label"].get(spec.raw, default_label)

    title = preset.get("title", f"{spec.raw} (SPH-smoothed, {projection})")

    return dict(scale=scale, cmap=cmap, floor=floor, label=label, title=title)


def resolve_clip_negative(spec: FieldSpec, overrides: dict) -> bool:
    """Whether to clip this field's negative values to 0 before weighting."""
    if spec.raw in overrides["clip_negative_off"]:
        return False
    if spec.raw in overrides["clip_negative_on"]:
        return True
    return FIELD_PRESETS.get(spec.name, {}).get("clip_negative", False)


def compute_frame(
    path: str,
    specs: list,
    clip_flags: dict,
    projection: str,
    resolution: int,
    z0: float = None,
    half: float = None,
    show_stars: bool = True,
    show_sinks: bool = True,
    star_age_max_myr: Optional[float] = None,
) -> dict:
    """Load one snapshot and compute every requested field's SPH-projected map.

    When `z0`/`half` are given (time-series mode), that frozen slab is
    reused instead of re-picking one from this snapshot's own brightest
    star, so every frame shares the same physical slab. Only plain numpy
    arrays are kept on the returned entry (not the swiftsimio dataset
    itself), so caching one entry per snapshot for a long series stays
    bounded in memory and does not hold the HDF5 file open.

    `show_stars`/`show_sinks`/`star_age_max_myr` are display filters only
    (see `select_sources`); they never affect slab selection, which is
    fixed before this function is even called in time-series mode.
    """
    data = sw.load(path)
    L = float(data.metadata.boxsize[0].value)
    if z0 is None or half is None:
        z0, half, _ = pick_slab(data, L)

    region = make_region(data.metadata.boxsize, L, z0, half)
    mass_map = project_gas(data, resolution=resolution, project="masses", region=region)

    z = data.gas.coordinates[:, 2].value
    slab_sel = np.abs(z - z0) < half
    n_slab = int(slab_sel.sum())

    maps = {}
    clip_notes = {}
    for spec in specs:
        field_array = resolve_field_array(data, spec)
        avg, n_clipped, mfrac_clipped = project_field(
            data,
            field_array,
            mass_map,
            region,
            resolution,
            projection,
            clip_flags[spec.raw],
            spec.key,
        )
        maps[spec.raw] = avg
        clip_notes[spec.raw] = (
            (n_clipped, mfrac_clipped) if clip_flags[spec.raw] else None
        )

    stars_xyz, sinks_xyz = select_sources(
        data, show_stars, show_sinks, star_age_max_myr
    )

    return {
        "path": path,
        "tag": os.path.splitext(os.path.basename(path))[0],
        "time_myr": float(data.metadata.time.to("Myr").value),
        "time_internal": float(data.metadata.time.value),
        "L": L,
        "z0": z0,
        "half": half,
        "n_slab": n_slab,
        "stars_xyz": stars_xyz,
        "sinks_xyz": sinks_xyz,
        "maps": maps,
        "clip_notes": clip_notes,
    }


def _draw_image_panel(ax, avg: np.ndarray, L: float, norm, cmap: str, floor):
    """Draw one projected map onto `ax`; returns the mappable for a colorbar."""
    masked = np.ma.masked_invalid(avg)
    if floor is not None:
        masked = np.ma.masked_less_equal(masked, floor)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("0.72")
    im = ax.imshow(
        masked.T, origin="lower", extent=[0.0, L, 0.0, L], norm=norm, cmap=cmap_obj
    )
    ax.set_facecolor("0.72")
    ax.set_aspect("equal")
    return im


def draw_map(
    avg: np.ndarray,
    L: float,
    stars_xyz: np.ndarray,
    sinks_xyz: np.ndarray,
    z0: float,
    half: float,
    n_slab: int,
    t_internal: float,
    t_myr: float,
    norm,
    cmap: str,
    title: str,
    cbar_label: str,
    outpath: str,
    floor,
    note: str,
) -> None:
    """Draw one smoothed-projection panel and write it to disk."""
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    im = _draw_image_panel(ax, avg, L, norm, cmap, floor)
    ax.set_xlabel("x [internal length units, kpc]")
    ax.set_ylabel("y [internal length units, kpc]")
    ax.set_title(
        f"{title}\n"
        f"t = {t_internal:.4f} internal = {t_myr:.1f} Myr; "
        f"SPH-projected slab |z - {z0:.3f}| < {half:.3f} "
        f"({n_slab} gas particles selected by centre; kernel overlap extends beyond)",
        fontsize=9,
    )
    overlay_sources(ax, stars_xyz, sinks_xyz, z0, half)
    if note:
        ax.text(
            0.02,
            0.02,
            note,
            transform=ax.transAxes,
            fontsize=7,
            color="black",
            va="bottom",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.5),
        )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(cbar_label, fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"wrote {outpath}")


def field_stats(avg: np.ndarray) -> tuple:
    """(min, max, mean) of an array's finite values, or three NaNs if none are finite."""
    finite = avg[np.isfinite(avg)]
    if finite.size:
        return float(finite.min()), float(finite.max()), float(finite.mean())
    return float("nan"), float("nan"), float("nan")


def render_frame(
    entry: dict,
    specs: list,
    styles: dict,
    outdir: str,
    suffix: str,
    projection: str,
) -> dict:
    """Draw every requested field's panel for one computed entry.

    Parameters
    ----------
    entry : dict
        One frame from `compute_frame`.
    specs : list of FieldSpec
    styles : dict
        {spec.raw: {scale, cmap, floor, label, title, norm}}, with `norm`
        added by the caller (data-driven per frame, or frozen across a
        time series -- see `main`).
    outdir : str
    suffix : str
        Appended to the output filename, before ".png".
    projection : str
        Only used for the printed min/max/mean summary line.

    Returns
    -------
    dict
        {"outputs": {spec.raw: path}, "stats": {spec.raw: (min, max, mean)}}.
    """
    tag = entry["tag"]
    L, z0, half, n_slab = entry["L"], entry["z0"], entry["half"], entry["n_slab"]
    stars_xyz, sinks_xyz = entry["stars_xyz"], entry["sinks_xyz"]
    t_internal, t_myr = entry["time_internal"], entry["time_myr"]
    outputs, stats = {}, {}

    for spec in specs:
        avg = entry["maps"][spec.raw]
        style = styles[spec.raw]
        clip_note = entry["clip_notes"][spec.raw]

        note_lines = []
        if style["floor"] is not None:
            n_below = int(np.sum((avg <= style["floor"]) | ~np.isfinite(avg)))
            note_lines.append(
                f"grey: {n_below}/{avg.size} pixels at/below floor {style['floor']:g}"
            )
        if clip_note is not None:
            n_clipped, mfrac_clipped = clip_note
            note_lines.append(
                f"negative values clipped to 0 before weighting: {n_clipped} "
                f"particles ({mfrac_clipped:.2%} of gas mass)"
            )

        outpath = os.path.join(outdir, f"{tag}_{spec.key}{suffix}.png")
        draw_map(
            avg,
            L,
            stars_xyz,
            sinks_xyz,
            z0,
            half,
            n_slab,
            t_internal,
            t_myr,
            style["norm"],
            style["cmap"],
            style["title"],
            style["label"],
            outpath,
            style["floor"],
            "\n".join(note_lines),
        )
        outputs[spec.raw] = outpath
        stats[spec.raw] = field_stats(avg)

    return {"outputs": outputs, "stats": stats}


def make_summary_panel(
    entries: list,
    spec: FieldSpec,
    z0: float,
    half: float,
    style: dict,
    outpath: str,
) -> None:
    """One combined figure: one small panel per snapshot, sharing the
    colour scale and colorbar, so the field's evolution is visible at a
    glance before diving into the individual per-frame figures."""
    n = len(entries)
    ncols = min(n, 6)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(2.6 * ncols, 2.6 * nrows), squeeze=False
    )
    mappable = None
    for i, entry in enumerate(entries):
        ax = axes[i // ncols][i % ncols]
        im = _draw_image_panel(
            ax,
            entry["maps"][spec.raw],
            entry["L"],
            style["norm"],
            style["cmap"],
            style["floor"],
        )
        mappable = im
        overlay_sources(
            ax, entry["stars_xyz"], entry["sinks_xyz"], z0, half, legend=False
        )
        ax.set_title(f"t={entry['time_myr']:.1f} Myr\n{entry['tag']}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(
        f"{style['title']} -- time series ({n} snapshots, shared colour scale)",
        fontsize=11,
    )
    if mappable is not None:
        cb = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.8)
        cb.set_label(style["label"], fontsize=9)
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"wrote {outpath}")


def main(argv: list) -> int:
    """Build the requested fields' SPH-smoothed maps for the given snapshot(s)."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        default=list(DEFAULT_FIELDS),
        help="Fields to plot, as NAME or NAME:COMPONENT (see module docstring). "
        "Default: the original four panels.",
    )
    parser.add_argument(
        "--projection",
        choices=("average", "sum"),
        default="average",
        help="'average' (default): mass-weighted mean per pixel. 'sum': raw "
        "projected sum per pixel (e.g. a column quantity).",
    )
    parser.add_argument("--scale", nargs="+", default=[], metavar="FIELD=log|linear")
    parser.add_argument("--cmap", nargs="+", default=[], metavar="FIELD=CMAP")
    parser.add_argument("--floor", nargs="+", default=[], metavar="FIELD=VALUE")
    parser.add_argument("--label", nargs="+", default=[], metavar="FIELD=TEXT")
    parser.add_argument(
        "--clip-negative",
        nargs="+",
        default=[],
        dest="clip_negative_on",
        metavar="FIELD",
        help="Force-clip these fields' negative values to 0 before weighting.",
    )
    parser.add_argument(
        "--no-clip-negative",
        nargs="+",
        default=[],
        dest="clip_negative_off",
        metavar="FIELD",
        help="Force-disable negative-value clipping for these fields.",
    )
    parser.add_argument("--resolution", type=int, default=RESOLUTION_DEFAULT)
    parser.add_argument(
        "--stars",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Overlay PartType4 stars inside the slab (default: on). --no-stars turns it off.",
    )
    parser.add_argument(
        "--sinks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Overlay PartType3 sinks inside the slab (default: on). --no-sinks turns it off.",
    )
    parser.add_argument(
        "--star-age-max",
        type=float,
        default=None,
        metavar="MYR",
        help="Show only stars younger than this age in Myr. Requires the snapshot's "
        "stars to have a BirthTimes field; ignored (with a warning) otherwise.",
    )
    parser.add_argument(
        "-s",
        "--snapshots",
        nargs="+",
        default=["snap/snapshot_0252.hdf5"],
        help="One or more snapshot globs/paths. Default mode (no --time-series): "
        "the newest complete snapshot among all of them is used. With "
        "--time-series: every complete snapshot matched by any of them is "
        "used, sorted by simulation time.",
    )
    parser.add_argument("-o", "--outdir", default=".", help="Output directory.")
    parser.add_argument(
        "--time-series",
        action="store_true",
        help="Build every requested field's map for every complete snapshot "
        "matched by --snapshots instead of just the newest one. One frozen "
        "colour scale per field and one frozen slab are used for every "
        "frame -- see freeze_slab(). Writes one combined "
        "<outdir>/<field>_smoothed_timeseries_summary.png panel per field "
        "(one small panel per snapshot, shared colour scale). Add --frames "
        "to also write each snapshot's individual "
        "<outdir>/<tag>_<field>_smoothed_timeseries.png panel.",
    )
    parser.add_argument(
        "--frames",
        action="store_true",
        help="In --time-series mode, also write each snapshot's individual "
        "per-field panel, not just the combined summary panel. Off by "
        "default, since a long series' individual frames are usually not needed.",
    )
    args = parser.parse_args(argv)

    specs = [parse_field_spec(token) for token in args.fields]
    required_names = sorted({spec.name for spec in specs})

    overrides = dict(
        scale=parse_kv_list(args.scale),
        cmap=parse_kv_list(args.cmap),
        floor=parse_kv_list(args.floor, value_type=float),
        label=parse_kv_list(args.label),
        clip_negative_on=set(args.clip_negative_on),
        clip_negative_off=set(args.clip_negative_off),
    )
    clip_flags = {spec.raw: resolve_clip_negative(spec, overrides) for spec in specs}
    os.makedirs(args.outdir, exist_ok=True)

    if not args.time_series:
        path = find_all_snapshots(args.snapshots, required_names)[-1]
        print(f"using {path}")
        entry = compute_frame(
            path,
            specs,
            clip_flags,
            args.projection,
            args.resolution,
            show_stars=args.stars,
            show_sinks=args.sinks,
            star_age_max_myr=args.star_age_max,
        )
        print(
            f"t = {entry['time_internal']:.5f} internal = {entry['time_myr']:.1f} Myr; "
            f"slab centre z0 = {entry['z0']:.4f}, half-thickness = {entry['half']:.4f}, "
            f"{entry['n_slab']} gas particles (by centre)"
        )
        styles = {}
        for spec in specs:
            arr = entry["maps"][spec.raw]
            finite = arr[np.isfinite(arr)]
            vmin, vmax = (
                (float(finite.min()), float(finite.max()))
                if finite.size
                else (None, None)
            )
            all_positive = bool(finite.size) and bool(np.all(finite > 0))
            style = resolve_style(
                spec, vmin, vmax, all_positive, args.projection, overrides
            )
            style["norm"] = (
                safe_log_norm(style["floor"], vmax)
                if style["scale"] == "log"
                else safe_linear_norm(vmin, vmax)
            )
            styles[spec.raw] = style
        render_frame(entry, specs, styles, args.outdir, "_smoothed", args.projection)
        return 0

    paths = find_all_snapshots(args.snapshots, required_names)
    print(f"time series: {len(paths)} complete snapshots")

    with h5py.File(paths[-1], "r") as f:
        boxsize = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
    z0, half = freeze_slab(paths, boxsize)
    print(
        f"frozen slab centre z0 = {z0:.4f}, half-thickness = {half:.4f} (reused for every frame)"
    )

    entries = [
        compute_frame(
            path,
            specs,
            clip_flags,
            args.projection,
            args.resolution,
            z0=z0,
            half=half,
            show_stars=args.stars,
            show_sinks=args.sinks,
            star_age_max_myr=args.star_age_max,
        )
        for path in paths
    ]

    styles = {}
    print("global colour ranges (from the projected maps):")
    for spec in specs:
        vmin, vmax, all_positive = combined_range(entries, spec.raw)
        style = resolve_style(
            spec, vmin, vmax, all_positive, args.projection, overrides
        )
        style["norm"] = (
            safe_log_norm(style["floor"], vmax)
            if style["scale"] == "log"
            else safe_linear_norm(vmin, vmax)
        )
        styles[spec.raw] = style
        print(
            f"  {spec.raw}: [{vmin:.6g}, {vmax:.6g}] scale={style['scale']}"
            if vmin is not None
            else f"  {spec.raw}: empty"
        )

    print(
        f"{'snapshot':>18} {'t_Myr':>8} {'field':>28} {'min':>12} {'max':>12} {'mean':>12}"
    )
    for entry in entries:
        if args.frames:
            stats = render_frame(
                entry,
                specs,
                styles,
                args.outdir,
                "_smoothed_timeseries",
                args.projection,
            )["stats"]
        else:
            stats = {spec.raw: field_stats(entry["maps"][spec.raw]) for spec in specs}
        for field_raw, (vmin, vmax, vmean) in stats.items():
            print(
                f"{entry['tag']:>18} {entry['time_myr']:8.2f} {field_raw:>28} {vmin:12.4g} {vmax:12.4g} {vmean:12.4g}"
            )

    for spec in specs:
        outpath = os.path.join(
            args.outdir, f"{spec.key}_smoothed_timeseries_summary.png"
        )
        make_summary_panel(entries, spec, z0, half, styles[spec.raw], outpath)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
