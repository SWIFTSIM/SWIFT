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
Compare the simulated growth of the HII region radius against the STARBENCH
semi-empirical D-type expansion formula (Bisbas et al. 2015, MNRAS 453,
1324, arXiv:1507.05621), Eqns 8, 11, 28 and 29 -- the same target Hu et al.
(2017, MNRAS 471, 2151, Appendix A) and Smith et al. (2021, MNRAS 506,
3882) validate their own D-type test against, rather than the classical
Spitzer (1978) solution alone (which Spitzer over-predicts against at late
times -- the whole reason STARBENCH exists).

  Eqn 8  (Raga-I):  the exact pressure-driven thin-shell ODE (Spitzer's
                     closed-form solution is what you get by dropping the
                     small mu_i*T_o/(mu_o*T_i) term from this).
  Eqn 11 (Raga-II): Hosokawa & Inutsuka (2006)'s ODE, including the inertia
                     of the swept-up shell (2nd order in R).
  Eqn 28: R_SB = R_II + f_SB * (R_I - R_II)
  Eqn 29: f_SB = 1 - 0.733 * exp(-t / 1 Myr)

Default config matches Hu et al. (2017) Appendix A's D-type convergence
test: 4 co-located 19.2 Msun sources (Q_H~2.5e48/s each, ~1e49/s combined),
n_H=100/cm3, t_end=8 Myr, at this code's own Z=0 ionized/neutral
temperature floors (T_i~47500 K, T_o~1e3 K via SPH:minimal_temperature),
not the papers' flat T_i=1e4 K -- the comparison is evaluated at whatever
T_i/T_o this run actually used (--T-ionized-K/--T-neutral-K), not
hardcoded to the papers' own numbers. See this example's README for the
box-size derivation and why Z=0 is used here instead of the
Z/Zsun~0.231 workaround used elsewhere in this project.

r_hii is read the same way as stromgren_analytic_check.py's corrected
measure: the current position of every gas particle *any* of the n_stars
co-located sources has ever tagged (HIIStarIDs), not any star's frozen
HIIRegionRadii (h_hii) -- generalizing the single-star version to Hu et
al.'s own r_HII definition ("the maximum radius where a gas particle with
an ionization fraction x_H+ > 0.95 can be found", i.e. a single front
radius for the whole source cluster, not per-star).

Usage:
    python3 hu_smith_analytic_check.py [-s snap/snapshot] [-o out.png]
"""

import argparse
import glob

import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy.integrate import solve_ivp

try:
    import swiftsimio as sw
except ImportError:
    sw = None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c). Copied verbatim from
    stromgren_analytic_check.py."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


# -----------------------------------------------------------------------------
# STARBENCH semi-empirical formula (Bisbas et al. 2015)
# -----------------------------------------------------------------------------
def starbench_curve(
    t_myr, R_St_pc, c_i_km_s, c_o_km_s, mu_i=0.5, mu_o=1.0, T_i_K=1e4, T_o_K=1e3
):
    """Solve Eqns 8 and 11 numerically and blend via Eqns 28-29.

    mu_i=0.5 (fully ionized pure H, proton+electron) and mu_o=1.0 (neutral
    atomic H) are the paper's implicit values -- confirmed by cross-checking
    against its own stated "ratio ~ 1/200" for the early-phase test
    (T_o=100 K) and its quoted c_i/c_o for the late-phase test; do not
    change these without re-deriving them the same way.

    @param t_myr Array of times (Myr) to evaluate at, must start > 0.
    @param R_St_pc Initial Stromgren radius (pc).
    @param c_i_km_s, c_o_km_s Ionized/neutral isothermal sound speeds (km/s).
    @return (R_I, R_II, R_SB) arrays in pc, same shape as t_myr.
    """
    pc_to_km = 3.0857e13
    Myr_to_s = 3.15576e13
    c_i = c_i_km_s * Myr_to_s / pc_to_km  # pc / Myr
    c_o = c_o_km_s * Myr_to_s / pc_to_km

    def rhs_raga1(t, y):
        R = y[0]
        ratio = (R_St_pc / R) ** 0.75
        term2 = (mu_i * T_o_K) / (mu_o * T_i_K) * (R_St_pc / R) ** (-0.75)
        return [c_i * (ratio - term2)]

    def rhs_raga2_conservative(t, y):
        # w = R^3 * dR/dt (proportional to the shell's physical momentum);
        # avoids the (3/R)(dR/dt)^2 term's numerical stiffness in the naive
        # R'' + (3/R)R'^2 = ... form.
        R, w = y
        dRdt = w / R**3
        dwdt = 3.0 * R_St_pc**1.5 * c_i**2 * np.sqrt(R) - 3.0 * c_o**2 * R**2
        return [dRdt, dwdt]

    t0 = min(t_myr.min(), 1e-4)
    t_span = [t0, t_myr.max()]

    sol1 = solve_ivp(rhs_raga1, t_span, [R_St_pc], t_eval=t_myr, rtol=1e-10, atol=1e-12)
    if sol1.status != 0:
        raise RuntimeError(f"Raga-I (Eqn 8) integration failed: {sol1.message}")

    # v(0) = c_i matches Raga-I's own initial expansion rate at R=R_St; a
    # larger initial velocity (e.g. sqrt(4/3)*c_i, Eqn 13's implied slope)
    # sends this ODE into an unstable oscillation that crashes through R=0
    # well before 3 Myr -- verified empirically, not just in principle.
    w0 = R_St_pc**3 * c_i
    sol2 = solve_ivp(
        rhs_raga2_conservative,
        t_span,
        [R_St_pc, w0],
        t_eval=t_myr,
        method="Radau",
        rtol=1e-10,
        atol=1e-12,
    )
    if sol2.status != 0:
        raise RuntimeError(f"Raga-II (Eqn 11) integration failed: {sol2.message}")

    R_I = sol1.y[0]
    R_II = sol2.y[0]
    f_SB = 1.0 - 0.733 * np.exp(-t_myr / 1.0)  # Eqn 29
    R_SB = R_II + f_SB * (R_I - R_II)  # Eqn 28
    return R_I, R_II, R_SB


# -----------------------------------------------------------------------------
# Q_H(mass): direct port of radiation_get_individual_star_ionizing_photon_
# emission_rate_fit() and the radius/luminosity fits it calls
# (src/feedback/GEAR/radiation.c). Copied verbatim from
# stromgren_analytic_check.py (not imported, to keep this example
# self-contained) -- keep both in sync if the fit ever changes.
# -----------------------------------------------------------------------------
def star_radius(mass_msun):
    """Empirical mass-radius relation (radiation_get_individual_star_radius)."""
    if mass_msun < 1.0:
        return (mass_msun**0.8) * const.R_sun
    elif mass_msun < 8.0:
        return (mass_msun**0.57) * const.R_sun
    else:
        return (mass_msun**0.5) * const.R_sun


def star_luminosity(mass_msun):
    """Empirical mass-luminosity relation (radiation_get_individual_star_luminosity)."""
    if mass_msun < 0.43:
        lum_sol = 0.185 * mass_msun**2
    elif mass_msun < 2.0:
        lum_sol = mass_msun**4
    elif mass_msun < 54.0:
        lum_sol = 1.5 * mass_msun**3.5
    else:
        lum_sol = 32000.0 * mass_msun
    return lum_sol * const.L_sun


def ionizing_photon_rate(mass_msun):
    """
    Q_H, the ionizing photon emission rate for a single star of the given
    mass, from a blackbody fit to its (R, L) -- direct port of
    radiation_get_individual_star_ionizing_photon_emission_rate_fit().
    """
    R = star_radius(mass_msun)
    L = star_luminosity(mass_msun)
    if R <= 0 or L <= 0:
        return 0.0 / u.s

    R_in_Rsun = (R / const.R_sun).decompose().value
    L_in_Lsun = (L / const.L_sun).decompose().value

    T = 5780.0 * (L_in_Lsun / R_in_Rsun**2) ** 0.25 * u.K

    E_threshold = 13.605 * u.eV
    x_0 = (E_threshold / (const.k_B * T)).decompose().value

    if x_0 > 45.0:
        return 0.0 / u.s

    photon_integral_sum = 0.0
    for n in range(1, 6):
        exp_term = np.exp(-n * x_0)
        if exp_term < 1e-10:
            break
        photon_integral_sum += (exp_term / n) * (x_0**2 + (2.0 * x_0) / n + 2.0 / n**2)

    prefactor = 15.0 / np.pi**4
    N_dot_ion = (L / (const.k_B * T)) * prefactor * photon_integral_sum
    return N_dot_ion.to(1 / u.s)


# -----------------------------------------------------------------------------
# Read the simulated r_hii(t) from the snapshots. Generalizes the
# single-star ever-tagged measure (stromgren_analytic_check.py,
# starbench_analytic_check.py) to n_stars co-located sources: a gas
# particle counts if it was ever tagged by *any* of them (np.isin, not a
# single ID match), matching Hu et al.'s own r_HII definition -- one front
# radius for the whole source cluster, not per-star.
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    times, r_hii = [], []
    star_mass_msun = None
    n_stars_found = None
    n_H_atom_cc = None
    boxsize_kpc = None

    for f in files:
        data = sw.load(f)
        if boxsize_kpc is None:
            boxsize_kpc = float(np.min(data.metadata.boxsize.to("kpc").value))
        if len(data.stars.hiiregion_radii) == 0:
            continue

        times.append(data.metadata.time.to("Myr").value)

        star_ids = data.stars.particle_ids.value.astype(np.int64)
        ever_tagged = np.isin(data.gas.hiistar_ids, star_ids)
        if np.any(ever_tagged):
            # Co-located sources: use the mean position as the common center.
            cluster_pos = data.stars.coordinates.to("kpc").value.mean(axis=0)
            gas_pos = data.gas.coordinates[ever_tagged].to("kpc").value
            box_kpc = data.metadata.boxsize.to("kpc").value
            dx = gas_pos - cluster_pos
            dx -= box_kpc * np.round(dx / box_kpc)
            r_hii.append(float(np.max(np.linalg.norm(dx, axis=1))))
        else:
            r_hii.append(float(np.max(data.stars.hiiregion_radii).to("kpc").value))

        if star_mass_msun is None:
            n_stars_found = len(data.stars.masses)
            masses = data.stars.masses.to("Msun").value
            if not np.allclose(masses, masses[0], rtol=1e-6):
                raise RuntimeError(
                    f"Expected {n_stars_found} equal-mass co-located sources, "
                    f"found masses {masses} Msun -- this script sums Q_H(mass) "
                    f"per identical star, not for a mixed population."
                )
            star_mass_msun = float(masses[0])
        if n_H_atom_cc is None:
            rho_g_cm3 = float(np.mean(data.gas.densities).to("g/cm**3").value)
            X_H_raw = data.metadata.parameters.get(
                "GrackleCooling:HydrogenFractionByMass"
            )
            X_H = float(X_H_raw) if X_H_raw is not None else 0.716
            n_H_atom_cc = (rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p).to(1 / u.cm**3)

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        star_mass_msun,
        n_stars_found,
        n_H_atom_cc,
        boxsize_kpc * u.kpc,
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot-glob", default="snap/snapshot_*.hdf5")
    parser.add_argument("-o", "--output", default="hu_smith_check.png")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.1,
        help="Relative-error tolerance for a pass (default: 0.1, "
        "tighter than the Spitzer check's 0.3 since STARBENCH "
        "agrees with high-resolution sims to <~2%%).",
    )
    parser.add_argument(
        "--T-ionized-K",
        type=float,
        default=47500.0,
        dest="T_ionized_K",
        help="Applied ionized-gas temperature floor (default: 47500, this code's own "
        "Z=0 floor). Pass 1e4 if this run instead used the Z/Zsun~0.231 workaround "
        "to hit Hu/Smith's own convention.",
    )
    parser.add_argument("--T-neutral-K", type=float, default=1e3, dest="T_neutral_K")
    args = parser.parse_args()

    t_sim, r_sim, star_mass_msun, n_stars, n_H, boxsize = read_simulated_r_hii(
        args.snapshot_glob
    )
    box_half_width = 0.5 * boxsize

    # n_stars identical co-located sources: sum Q_H per star, not
    # Q_H(n_stars * mass) -- the fit is nonlinear in mass.
    Q_H = n_stars * ionizing_photon_rate(star_mass_msun)
    alpha_B = alpha_b_hui_gnedin(args.T_ionized_K * u.K)
    R_St = ((3 * Q_H / (4 * np.pi * alpha_B * n_H**2)) ** (1 / 3.0)).to(u.pc)
    c_i = (np.sqrt(const.k_B * args.T_ionized_K * u.K / (0.5 * const.m_p))).to(
        u.km / u.s
    )
    c_o = (np.sqrt(const.k_B * args.T_neutral_K * u.K / const.m_p)).to(u.km / u.s)

    print(f"Star mass (each)   : {star_mass_msun:.3f} Msun x {n_stars} sources")
    print(f"n_H                : {n_H:.4g}")
    print(f"Q_H (combined)     : {Q_H.to(1/u.s).value:.3e} 1/s")
    print(f"R_St               : {R_St:.4g}")
    print(f"c_i = {c_i:.4g}, c_o = {c_o:.4g}")
    print(f"Box half-width     : {box_half_width.to(u.pc):.4g}")

    t_grid = np.linspace(max(t_sim.min().value, 1e-4), t_sim.max().value, 400) * u.Myr
    R_I, R_II, R_SB = starbench_curve(
        t_grid.to(u.Myr).value,
        R_St.to(u.pc).value,
        c_i.value,
        c_o.value,
        T_i_K=args.T_ionized_K,
        T_o_K=args.T_neutral_K,
    )

    # Interpolate the STARBENCH curve onto the simulation's own snapshot
    # times for the comparison-point error (avoid re-solving the ODE once
    # per snapshot).
    R_SB_at_t_sim = np.interp(t_sim.to(u.Myr).value, t_grid.to(u.Myr).value, R_SB)

    r_sim_pc = r_sim.to(u.pc).value
    box_valid = R_SB_at_t_sim <= box_half_width.to(u.pc).value
    alive = r_sim_pc > 0
    ok = box_valid & alive
    if not np.any(ok):
        raise RuntimeError("No box-valid, alive snapshot to compare at.")
    last = np.where(ok)[0][-1]
    t_c = t_sim[last]
    r_c = r_sim_pc[last]
    R_SB_c = R_SB_at_t_sim[last]
    rel_error = abs(r_c - R_SB_c) / R_SB_c
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"\nComparison time: t={t_c:.4g}  r_sim={r_c:.4g} pc  "
        f"R_SB={R_SB_c:.4g} pc  rel_error={rel_error:.2%}  [{verdict}]"
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(
        t_sim.to(u.Myr),
        r_sim_pc,
        "o-",
        color="#d62728",
        markersize=4,
        label="Simulation (ever-tagged $r_{\\rm HII}$)",
    )
    ax.plot(t_grid, R_I, ":", color="grey", label="Raga-I (Eqn 8)")
    ax.plot(t_grid, R_II, "-.", color="grey", label="Raga-II (Eqn 11)")
    ax.plot(t_grid, R_SB, "--", color="black", label="STARBENCH (Eqns 28-29)")
    ax.axhline(
        box_half_width.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        alpha=0.6,
        label="Box half-width",
    )
    ax.axvline(
        t_c.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.4,
        label="Comparison point",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"HII region radius $r_{\rm HII}$ [pc]")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
