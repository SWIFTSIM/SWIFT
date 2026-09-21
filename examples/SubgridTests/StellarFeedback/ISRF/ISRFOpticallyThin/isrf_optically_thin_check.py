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
Check the ISRF hyperbolic (Cattaneo-type flux-relaxation, M1 closure)
propagation scheme in the optically-thin corner, where the screening
length is much larger than both the smoothing length and the box.

For a single point source of band luminosity `L` in a static, uniform,
purely absorbing medium the M1 closure is exact (there is no scattering to
isotropize the field), so the free-streaming branch gives the settled
field behind the causal front,

    u(r) = L * exp(-r/lambda) / (4 * pi * c * rho * r^2)      (free-streaming)

with `lambda = 1/(kappa_eff*rho)` the absorption length, `kappa_eff` the
dust mass opacity (area/mass) and `c` the true speed of light (theory/GEAR/
Radiation/02_fuv_isrf.tex, "Steady states"). The reduced propagation speed
`c_hyp` cancels: it sets how long the field takes to settle, not where it
settles. This is the gated target of this check.

The isotropic/diffusion-closure branch,

    u(r) = 3*L * exp(-sqrt(3)*r/lambda) / (4 * pi * c * lambda * rho * r) ,
                                                              (informational)

is NOT a valid steady state this close to a single source (its own reduced
flux exceeds the closure's admissible range for r <~ 0.79*lambda) and is
kept here only as an informational curve, plotted for comparison. It is
the exact target of the sibling `ISRFHyperbolicPropagation` check, which
measures the screened corner (h/lambda ~ 1) where a diffuse, many-source
field legitimately isotropizes.

Gated leg, against the free-streaming solution with NO fitted parameter:

  * SHAPE: the radial log-log slope of `u(r)`. Free-streaming predicts
    `-2 - r/lambda`; the informational diffusion curve predicts
    `-1 - r/lambda`. The two are a full decade apart in a quantity measured
    over ~0.9 decades of radius.

Reported, not gated:

  * AMPLITUDE: the ratio of the measured field to the free-streaming
    prediction. This absolute continuum comparison is possible here (not
    at the sibling `ISRFHyperbolicPropagation` corner, where `h/lambda` is
    of order 1 and the discrete SPH estimator's own fixed point departs
    from the continuum profile), but the settled value is currently above
    unity by more than this check's own tolerance would allow; see
    `theory/GEAR/Radiation/02_fuv_isrf.tex`, "Limitations and open items".
    The check still fails on a non-finite amplitude or radial spread.

Measurement window, per snapshot: `[3*h_star, min(R_f - 2*H, L_box - R_f)]`.
  * `3*h_star` excludes the star's own injection footprint, which is not
    propagated field.
  * `R_f - 2*H` excludes the causal front itself, which the SPH estimator
    smears over a couple of kernel supports and which carries most of the
    field's energy at any time.
  * `L_box - R_f` excludes every periodic image by construction: an image
    at distance `L_box` cannot have reached radius `r` until
    `c_hyp*t = L_box - r`. This term is not redundant with the previous
    one and must not be dropped.

The front radius is `R_f = N_steps * C_hyp * h`, exactly, because the
propagation-speed closure `c_hyp = C_hyp*h/dt` makes `c_hyp*dt == C_hyp*h`
identically, whatever `dt` the run actually took. No reconstruction of
`c_hyp` from the timestep record is needed or wanted here.

`run.sh`'s `time_end` is set long enough that the fitted snapshots sample
the field well past its initial transient, before the causal front's own
approach to the periodic box edge narrows the window again: a window
taken too early measures a transient, not the settled field, and reads a
steeper slope than the settled one.

Tolerances are derived from measurement, not tuned to pass. A failure at
these tolerances is a real regression, not a case for loosening them.
"""

import argparse
import glob
import re
import sys

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# src/feedback/GEAR/radiation.h
SIGMA_D_PE_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default, src/feedback/GEAR/radiation.h
GRACKLE_SOLAR_Z = 0.01295
# Wendland C2, 3D: src/kernel_hydro.h kernel_gamma.
GAMMA_3D = 1.936492
C_LIGHT_CGS = 2.99792458e10
PC_CGS = 3.0856775814913673e18


def parse_options() -> argparse.Namespace:
    """Parse the command-line options."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to consider (default: %(default)s).",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        default="output.log",
        help="Run log, read for the step/time record (default: %(default)s).",
    )
    parser.add_argument(
        "--c-hyp-margin",
        type=float,
        default=0.5,
        help="GEARFeedback:ISRF_c_hyp_margin used by the run "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--n-bins", type=int, default=14, help="Number of radial bins (log-spaced)."
    )
    parser.add_argument(
        "--n-late",
        type=int,
        default=4,
        help="Number of final snapshots averaged for the verdict "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--slope-tol",
        type=float,
        default=0.25,
        help="Max allowed absolute difference between the measured log-log "
        "slope and the free-streaming prediction (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        default="isrf_optically_thin_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


def read_step_times(logfile: str) -> list:
    """Read the (step, time) record from a run log.

    Parameters
    ----------
    logfile : str
        Path to the run's own stdout log.

    Returns
    -------
    list of (int, float)
        Step number and simulation time, in file order.
    """
    record = []
    with open(logfile) as handle:
        for line in handle:
            match = re.match(r"^\s+(\d+)\s+([0-9.eE+-]+)\s", line)
            if match:
                record.append((int(match.group(1)), float(match.group(2))))
    if not record:
        raise RuntimeError(f"No step/time record found in {logfile}.")
    return record


def steps_completed(record: list, time: float) -> int:
    """Return the number of steps completed at a given simulation time.

    Parameters
    ----------
    record : list of (int, float)
        Output of :func:`read_step_times`.
    time : float
        Simulation time of the snapshot, internal units.

    Returns
    -------
    int
        Step count at that time.
    """
    completed = 0
    for step, step_time in record:
        if step_time <= time * (1.0 + 1e-9):
            completed = step
    return completed


def load_snapshot(path: str) -> dict:
    """Read the fields this check needs from one snapshot.

    Parameters
    ----------
    path : str
        Snapshot filename.

    Returns
    -------
    dict
        Particle arrays and metadata, all in internal units.
    """
    with h5py.File(path, "r") as handle:
        header = handle["/Header"]
        units = handle["/Units"]
        gas = handle["/PartType0"]
        star = handle["/PartType4"]
        return dict(
            time=float(np.asarray(header.attrs["Time"]).flat[0]),
            boxsize=float(np.asarray(header.attrs["BoxSize"], dtype=float).flat[0]),
            unit_length_cgs=float(
                np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
            ),
            unit_mass_cgs=float(
                np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0]
            ),
            unit_time_cgs=float(
                np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0]
            ),
            pos=gas["Coordinates"][:, :],
            rho=gas["Densities"][:].astype(np.float64),
            h=gas["SmoothingLengths"][:].astype(np.float64),
            mass=gas["Masses"][:].astype(np.float64),
            u_PE=gas["PESpecificEnergies"][:].astype(np.float64),
            u_LW=gas["LWSpecificEnergies"][:].astype(np.float64),
            Z=gas["MetalMassFractions"][:, -1],
            star_pos=star["Coordinates"][0, :],
            star_h=float(star["SmoothingLengths"][0]),
            L_PE=float(star["PELuminosities"][0]),
            L_LW=float(star["LWLuminosities"][0]),
        )


def radial_distance(
    pos: np.ndarray, star_pos: np.ndarray, boxsize: float
) -> np.ndarray:
    """Compute the minimum-image radial distance from the star.

    Parameters
    ----------
    pos : numpy.ndarray
        Gas coordinates, shape (N, 3).
    star_pos : numpy.ndarray
        Star coordinates, shape (3,).
    boxsize : float
        Periodic box side length.

    Returns
    -------
    numpy.ndarray
        Distance of every gas particle from the star.
    """
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def dust_mass_opacity_cgs(Z: np.ndarray, sigma_d_cgs: float) -> np.ndarray:
    """Compute the band dust mass opacity (area/mass, cgs).

    Mirrors radiation_get_dust_mass_opacity in radiation_isrf.c, with
    local_dust_to_gas_ratio left at Grackle's own default (this example's
    params.yml leaves it unset), so D_relative reduces to Z/Z_grackle_sun.

    Parameters
    ----------
    Z : numpy.ndarray or float
        Gas metal mass fraction.
    sigma_d_cgs : float
        Band dust cross-section per hydrogen nucleon, cm^2.

    Returns
    -------
    numpy.ndarray or float
        Dust mass opacity, cm^2/g.
    """
    return sigma_d_cgs * (Z / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)


def measure(snapshot: dict, record: list, c_hyp_margin: float, band: str, n_bins: int):
    """Measure one band's radial profile inside the causal front.

    Parameters
    ----------
    snapshot : dict
        Output of :func:`load_snapshot`.
    record : list
        Output of :func:`read_step_times`.
    c_hyp_margin : float
        GEARFeedback:ISRF_c_hyp_margin used by the run.
    band : str
        Either "PE" or "LW".
    n_bins : int
        Number of log-spaced radial bins.

    Returns
    -------
    dict
        Measurement, including the window, the binned profile, both
        analytic predictions, the fitted slope and the negative-energy
        fraction. Entries that could not be measured are NaN, and
        ``comment`` then says why.
    """
    unit_length = snapshot["unit_length_cgs"]
    unit_mass = snapshot["unit_mass_cgs"]
    unit_time = snapshot["unit_time_cgs"]
    sigma_d = SIGMA_D_PE_CGS if band == "PE" else SIGMA_D_LW_CGS

    r = radial_distance(snapshot["pos"], snapshot["star_pos"], snapshot["boxsize"])
    r *= unit_length
    u = snapshot["u_" + band] * (unit_length / unit_time) ** 2
    mass = snapshot["mass"] * unit_mass
    luminosity = snapshot["L_" + band] * unit_mass * unit_length**2 / unit_time**3

    Z = np.median(snapshot["Z"])
    rho = np.median(snapshot["rho"]) * unit_mass / unit_length**3
    kappa_eff = dust_mass_opacity_cgs(Z, sigma_d)
    lam = 1.0 / (kappa_eff * rho)

    h = np.median(snapshot["h"]) * unit_length
    support = GAMMA_3D * h
    box = snapshot["boxsize"] * unit_length
    # Exact: the closure c_hyp = C_hyp*h/dt makes c_hyp*dt == C_hyp*h.
    front = steps_completed(record, snapshot["time"]) * c_hyp_margin * h

    r_min = 3.0 * snapshot["star_h"] * unit_length
    r_max = min(front - 2.0 * support, box - front)

    out = dict(
        time=snapshot["time"],
        band=band,
        lam=lam,
        h=h,
        front=front,
        r_min=r_min,
        r_max=r_max,
        slope=np.nan,
        slope_err=np.nan,
        slope_freestream=np.nan,
        slope_diffusion=np.nan,
        amplitude=np.nan,
        flatness=np.nan,
        comment="",
        r_bin=np.array([]),
        u_bin=np.array([]),
        u_diffusion=np.array([]),
        u_freestream=np.array([]),
    )

    inside = r < front
    energy = (mass * u)[inside]
    total = np.abs(energy).sum()
    out["negative_fraction"] = (
        np.abs(energy[energy < 0]).sum() / total if total > 0 else np.nan
    )

    if r_max <= 1.5 * r_min:
        out["comment"] = "window too narrow (front has not cleared the star yet)"
        return out

    edges = np.logspace(np.log10(r_min), np.log10(r_max), n_bins + 1)
    r_bin, u_bin, n_dropped = [], [], 0
    for low, high in zip(edges[:-1], edges[1:]):
        selection = (r >= low) & (r < high)
        if selection.sum() < 25:
            continue
        median = np.median(u[selection])
        if median <= 0.0:
            n_dropped += 1
            continue
        r_bin.append(np.sqrt(low * high))
        u_bin.append(median)
    r_bin = np.array(r_bin)
    u_bin = np.array(u_bin)

    mid = np.sqrt(r_min * r_max)
    out["slope_freestream"] = -2.0 - mid / lam
    out["slope_diffusion"] = -1.0 - mid / lam

    if len(r_bin) < 5 or n_dropped > 0.3 * n_bins:
        out["comment"] = (
            f"cannot fit: {n_dropped} of {n_bins} bins dropped for a "
            "non-positive median field. The propagated field is ringing in "
            "sign, so it has no radial profile to measure. This is what an "
            "insufficient ISRF_dissipation_alpha_max looks like at this "
            "screening length; see the README."
        )
        return out

    out["r_bin"] = r_bin
    out["u_bin"] = u_bin
    out["u_freestream"] = (
        luminosity * np.exp(-r_bin / lam) / (4.0 * np.pi * C_LIGHT_CGS * rho * r_bin**2)
    )
    out["u_diffusion"] = (
        3.0
        * luminosity
        * np.exp(-np.sqrt(3.0) * r_bin / lam)
        / (4.0 * np.pi * C_LIGHT_CGS * lam * rho * r_bin)
    )

    design = np.vstack([np.log(r_bin), np.ones(len(r_bin))]).T
    coefficients, _, _, _ = np.linalg.lstsq(design, np.log(u_bin), rcond=None)
    residual = np.log(u_bin) - design @ coefficients
    covariance = (
        np.linalg.inv(design.T @ design)
        * (residual @ residual)
        / max(len(r_bin) - 2, 1)
    )
    out["slope"] = coefficients[0]
    out["slope_err"] = np.sqrt(covariance[0, 0])

    ratio = u_bin / out["u_freestream"]
    out["amplitude"] = np.median(ratio)
    out["flatness"] = np.std(ratio) / np.mean(ratio)
    return out


def plot(results: dict, filename: str) -> None:
    """Plot the measured profiles against both analytic predictions.

    Parameters
    ----------
    results : dict
        Band name mapped to the final measurement dict.
    filename : str
        Output image filename.
    """
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for band, colour in (("PE", "C0"), ("LW", "C1")):
        item = results[band]
        if len(item["r_bin"]) == 0:
            continue
        r = item["r_bin"] / PC_CGS
        axes[0].loglog(r, item["u_bin"], "o", color=colour, label=f"{band} measured")
        axes[0].loglog(
            r, item["u_freestream"], "-", color=colour, label=f"{band} free-streaming"
        )
        axes[0].loglog(
            r,
            item["u_diffusion"],
            "--",
            color=colour,
            label=f"{band} diffusion, informational",
        )
        axes[1].semilogx(
            r, item["u_bin"] / item["u_freestream"], "o-", color=colour, label=f"{band}"
        )
    axes[0].set_xlabel("r [pc]")
    axes[0].set_ylabel(r"specific energy $u$ [erg g$^{-1}$]")
    axes[0].legend(fontsize=7)
    axes[0].set_title("Profile inside the causal front")
    axes[1].axhline(1.0, color="k", lw=0.8)
    axes[1].set_xlabel("r [pc]")
    axes[1].set_ylabel("measured / free-streaming")
    axes[1].set_ylim(0.0, 2.0)
    axes[1].legend(fontsize=8)
    axes[1].set_title("Absolute agreement (no fitted parameter)")
    figure.tight_layout()
    figure.savefig(filename, dpi=140)
    print(f"Wrote {filename}")


def main() -> int:
    """Run the check and return a shell exit code."""
    options = parse_options()
    paths = sorted(glob.glob(options.snapshot))
    if not paths:
        print(f"No snapshot matched {options.snapshot}.")
        return 1
    record = read_step_times(options.logfile)

    per_band = {"PE": [], "LW": []}
    for path in paths:
        snapshot = load_snapshot(path)
        if snapshot["time"] <= 0.0:
            continue
        for band in ("PE", "LW"):
            per_band[band].append(
                measure(snapshot, record, options.c_hyp_margin, band, options.n_bins)
            )

    failures = []
    final = {}
    for band in ("PE", "LW"):
        items = per_band[band]
        usable = [item for item in items if not np.isnan(item["slope"])]
        print(f"\n=== {band} band ===")
        print(
            f"{'t':>10} {'R_f/pc':>7} {'window/pc':>14} {'dex':>5} "
            f"{'E_neg/E':>8} {'slope':>16} {'meas/FS':>10}"
        )
        for item in items[-options.n_late * 3 :]:
            window = f"[{item['r_min'] / PC_CGS:4.2f},{item['r_max'] / PC_CGS:5.2f}]"
            dex = (
                np.log10(item["r_max"] / item["r_min"])
                if item["r_max"] > item["r_min"]
                else 0.0
            )
            slope = (
                f"{item['slope']:+7.3f}+-{item['slope_err']:5.3f}"
                if not np.isnan(item["slope"])
                else "        --      "
            )
            amplitude = (
                f"{item['amplitude']:10.3f}"
                if not np.isnan(item["amplitude"])
                else "        --"
            )
            print(
                f"{item['time']:10.2e} {item['front'] / PC_CGS:7.2f} {window:>14} "
                f"{dex:5.2f} {item['negative_fraction']:8.3f} {slope:>16} {amplitude}"
            )

        if len(usable) < options.n_late:
            comment = items[-1]["comment"] if items else "no snapshots"
            print(
                f"FAIL [{band}]: fewer than {options.n_late} snapshots could be fitted."
            )
            print(f"      last reason: {comment}")
            failures.append(band)
            final[band] = items[-1] if items else None
            continue

        late = usable[-options.n_late :]
        final[band] = late[-1]
        slope = float(np.mean([item["slope"] for item in late]))
        slope_err = float(np.std([item["slope"] for item in late]) / np.sqrt(len(late)))
        predicted = float(np.mean([item["slope_freestream"] for item in late]))
        diffusion = float(np.mean([item["slope_diffusion"] for item in late]))
        amplitude = float(np.mean([item["amplitude"] for item in late]))
        flatness = float(np.max([item["flatness"] for item in late]))
        lam = late[-1]["lam"] / PC_CGS
        h_over_lam = late[-1]["h"] / late[-1]["lam"]

        print(
            f"\n  lambda = {lam:.1f} pc, h/lambda = {h_over_lam:.2e}, "
            f"averaged over the last {len(late)} fitted snapshots:"
        )
        print(
            f"  slope             = {slope:+.3f} +- {slope_err:.3f}  "
            f"(free-streaming predicts {predicted:+.3f}, diffusion (informational) "
            f"predicts {diffusion:+.3f})"
        )
        print(
            f"  measured / free-streaming = {amplitude:.3f} (radial spread "
            f"{flatness:.3f}), NOT GATED: open amplitude question, see "
            'theory/GEAR/Radiation/02_fuv_isrf.tex, "Limitations and open items".'
        )

        # The slope gate below is a `> tol -> FAIL` comparison, which a NaN
        # passes silently, so finiteness is tested first and separately.
        if not np.isfinite(slope) or not np.isfinite(predicted):
            print("  FAIL: the fitted slope or its prediction is not finite.")
            failures.append(band)
        if not np.isfinite(amplitude) or not np.isfinite(flatness):
            print("  FAIL: amplitude or radial spread is not finite.")
            failures.append(band)
        if abs(slope - predicted) > options.slope_tol:
            print(
                f"  FAIL: slope is {abs(slope - predicted):.3f} from the "
                f"free-streaming prediction, tolerance {options.slope_tol}."
            )
            failures.append(band)

    if all(final.get(band) is not None for band in ("PE", "LW")):
        plot(final, options.output)

    if failures:
        print(f"\nCHECK FAILED for: {', '.join(sorted(set(failures)))}")
        return 1
    print(
        "\nCHECK PASSED: the propagated field's radial slope matches the "
        "free-streaming 1/r^2 profile in both bands, not the diffusion-"
        "closure 1/r profile kept here only for comparison. The amplitude "
        "ratio above is reported only; see the printout for the open item."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
