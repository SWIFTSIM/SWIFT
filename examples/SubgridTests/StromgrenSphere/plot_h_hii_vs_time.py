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
Plot r_hii(t) per star from a run's snapshots.

Two star fields are involved and neither alone gives a clean curve:
- `HIIRegionRadii` (stars_io.h) is the LIVE value (h_hii * kernel_gamma),
  read fresh from the star each snapshot dump. It is reset to 0 the moment
  the star dies or ages out of HII eligibility
  (feedback_will_do_feedback, feedback_common.c), so a raw plot of this
  field alone shows every star's curve crash to zero at death.
- `FinalHIIRegionRadii` (tracers_io.h) is a ONE-TIME terminal snapshot of
  the same quantity, written at that same instant, and held constant
  afterwards.

This script stitches the two per star, per snapshot: use the live value
while it is nonzero, otherwise fall back to the final value -- giving a
curve that rises during the star's active HII lifetime and then holds flat
after death/ineligibility, instead of dropping back to zero.

Usage:
    python3 plot_h_hii_vs_time.py [-s snap/snapshot] [-o out.png]
"""
import argparse
import glob

import numpy as np

try:
    import swiftsimio as sw
except ImportError:
    sw = None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_per_star_r_hii(snapshot_glob):
    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    # star_id -> (list of times [Myr], list of r_hii [pc], birth mass [Msun])
    per_star = {}

    for f in files:
        data = sw.load(f)
        ids = np.asarray(data.stars.particle_ids.value)
        if len(ids) == 0:
            continue

        t_myr = float(data.metadata.time.to("Myr").value)
        r_live = data.stars.hiiregion_radii.to("pc").value
        r_final = data.stars.final_hiiregion_radii.to("pc").value
        masses = data.stars.masses.to("Msun").value

        # Stitch: live value while the star is actively growing an HII
        # region; fall back to the terminal snapshot once it has died or
        # aged out (live field reset to 0 at that point).
        r_eff = np.where(r_live > 0.0, r_live, r_final)

        for sid, r, m in zip(ids, r_eff, masses):
            sid = int(sid)
            entry = per_star.setdefault(sid, {"t": [], "r": [], "mass": m})
            entry["t"].append(t_myr)
            entry["r"].append(float(r))

    return per_star


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot-glob",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to read (default: snap/snapshot_*.hdf5)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="h_hii_vs_time.png",
        help="Output plot filename.",
    )
    args = parser.parse_args()

    per_star = read_per_star_r_hii(args.snapshot_glob)

    fig, ax = plt.subplots()
    for sid, entry in sorted(per_star.items()):
        t = np.asarray(entry["t"])
        r = np.asarray(entry["r"])
        if np.all(r == 0.0):
            # Star never ionized anything (e.g. too low-mass, or eligible
            # window closed before the first HII rebuild) -- skip clutter.
            continue
        ax.plot(t, r, "o-", ms=3, label=f"id={sid} ({entry['mass']:.1f} Msun)")

    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"$r_{hii}$ [pc]")
    ax.legend(fontsize="small")
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
