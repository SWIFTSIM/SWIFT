#!/usr/bin/env python3
"""Check the luminosity and the H2 dissociation rate of one star in thin gas.

One star sits in an isolated, dust-free box of H2-poor gas. The check reads the
luminosity L of the star from the radiation table, not from the transported field.
Three tests gate the run. Two other numbers are reported without a bar.

Luminosity (gate). The star's LWLuminosities against the table L_LW as SWIFT reads
it. For Z = 0, below the first table row, SWIFT takes the first row and interpolates
log10 L linearly in log M between the two nodes of its own mass grid around the star's
mass (interpolate_2d, clamped on the Z axis). SWIFT's mass grid runs from the IMF minimum to the IMF
maximum mass of the table in GEARRadiation:interpolation_size_mass steps. The check
follows that rule where this grid equals the mass grid of the table file. Where the two
differ, the run is void (exit 4). The bar is 8 float32 ulp of log10 L (a chosen margin
on the one ulp that the float32 readout cannot resolve), as a relative error. The check
tests the readout of the table, not whether the boundary rule is the one the model
wants. A star with Z inside the table range, or a mass outside it, is a void run.

Rate identity (gate, test 4). Per particle ln(k_meas / (kappa_D c rho u_diss)), with
k_meas = -d ln x2 / dt from two snapshots, kappa_D = 1.4425e-7 cm^2 erg^-1 (the
Draine-field rate per unit LW flux), c rho u_diss the field the cooling receives
(LWDissociatingSpecificEnergies, averaged over the interval; it already holds the
line transmission). Its median per shell must lie within the timing term plus
2 se. The timing term is the shell median of |ln(u_diss,b / u_diss,a)|: it covers
whether the cooling reads the field of the step before or after. A second shielding of
the transported rate by the library multiplies k_meas by a factor far below 1 and
fails it.

Stored-energy fraction Q (reported, no bar). With no dust the energy that crosses a
sphere around the star per step is L dt. By the field formula of the transport
(theory, Eq. packet-field) the volume mean over a shell of c rho u, with u the
particle's own LW energy before line absorption, over L / (4 pi r^2) is

    Q(r) = sum_i V_i c rho_i u_i / sum_i V_i L / (4 pi r_i^2),   V_i = m_i / rho_i.

Q depends on the cone power of the transport. The reported Q is printed with its
bootstrap error and not gated. The ratio of the measured dissociation to the
kernel-mean thin limit (2b) is reported the same way: it carries the same Q.

Geometry bar on Q (gate, test Q). The bar is two-sided and comes from the geometry of
the run, not from a clean run. A model of the transport on the gas of the run gives the
expected Q of three radial shells, and the scatter of that expectation.

Model. Symbols: r [pc] distance to the star; V [pc^3] volume m / rho of a particle;
h [pc] smoothing length of a particle; h_s [pc] smoothing length of the star;
gamma = 1.936492 (kernel support over h); p (dimensionless) the cone power
GEARFeedback:ISRF_packet_cone_power; mu (dimensionless) the cosine between the
direction of the light of the sender (radial from the star) and the vector to the
receiver; l [pc] the weighted mean advance r mu of the hops of one particle; E
(dimensionless) the energy of a particle in units of the energy emitted per step.
  1. Hop weight from sender i to receiver j: V_j W(r_ij / (gamma H)) mu^p / H^3, with W
     the Wendland C2 kernel, H = max(h_i, h_j), for r_ij < gamma h_i or gamma h_j and
     mu > 1e-3 (radiation_packet.c, hop weights). The weights of one sender are
     normalised to 1 and form the matrix Phi.
  2. Emission per step: the gas with r < 0.95 gamma h_s gets the tent weight
     V (1 - r / (0.95 gamma h_s)) (star_gas_weights.h). The weights are summed in the 12
     HEALPix cells (nside 1) of a frame rotated at random each step. A cell whose sum is
     below the cone floor 0.5 (pi / 3) (0.95 gamma h_s)^3 / 12 sends the missing fraction
     of its 1 / 12 out of the box. The delivered light, 1 minus that escaped fraction, is
     shared over all receivers in proportion to their weight (star_gas_global_share in
     star_gas_weights.h). This gives the emission vector e.
  3. One step: E <- Phi E + e (unlimited light speed: one hop per step). The stored
     energy of a particle gives Q = 4 pi r^2 E l / V (the field formula of the transport).
  4. A snapshot sees E after the steps that ended at or before its time. Q of a shell is
     the volume mean, then the mean over the snapshots of the check window.
The expectation is the mean of that shell Q over N_SEQUENCE independent random frame
sequences (fixed seeds). The bar is the expectation +- 3 sigma, with sigma the standard
deviation over the sequences of the window mean. The run fails if one shell lies
outside. No datum of the Q of the run enters. The shells are fixed multiples of the
median gas h: BAR_SHELLS_H.

Q does not depend on the number of reception slots. With one star every share that
reaches a particle has the star as origin, so the arrivals of a particle fall in one
slot. The bar does not test the slot count.

What a PASS proves: the star's table luminosity matches the SWIFT value, the H2
dissociation follows the field that the cooling is given with no further factor (no
second shielding), and the stored light of the thin box equals the geometry model to 3
sigma per shell (a scale error of Q of about 1 per cent is rejected). What it does not
prove: that the model itself is right where it was not tested. The Q bar is a limit of
the model: it is gated only at resolution level 5 and a star of 29.7 Msun, and other
levels, star masses and slot counts are void or untested. A loss of light that grows with the H2
column or the dust is not tested (the box is thin and dust-free), nor are errors that
appear only with a limited light speed or several stars. The rate identity alone does
not see a scale error of the LW field, because the rate and the field scale together;
the Q bar does.

Exit status: 0 pass, 1 a failed bar, 2 a missing file or field, 3 a value that is not
finite, 4 a void Q bar with no other failure (a failed gate always gives 1; the
luminosity and the rate identity run at every level and star mass).

The --plant-* options change the data as a known error would, to show that the check
fails when it should. They act on the arrays read from the snapshots, not on SWIFT.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]

C_LIGHT = 2.99792458e10
MSUN = 1.988409870698051e33
PARSEC = 3.0856775814913673e18
KAPPA_D = 1.4425206587e-7
KB = 1.380649e-16
MH = 1.67262192369e-24
GRAV = 6.6743e-8
GAMMA_KERNEL = 1.936492
N_LIT = 4
N_MIN = 100
N_BOOT = 300
N_SEQUENCE = 24
BAR_LEVEL = 5
BAR_STAR_MASS_MSUN = 29.7
BAR_SIGMA = 3.0
BAR_SHELLS_H = (5.6, 7.7, 9.2, 11.2)
MU_MIN = 1e-3
STAR_EDGE = 0.95
CONE_FLOOR_FRACTION = 0.5
N_CONE = 12
SPECIES = ["HI", "HII", "HM", "H2I", "H2II", "HeI", "HeII", "HeIII", "DI", "DII", "HDI"]
SPECIES_A = [1, 1, 1, 2, 2, 4, 4, 4, 2, 2, 3]
FIELDS = [
    "LWOwnSpecificEnergies",
    "LWOwnDissociatingSpecificEnergies",
    "LWSpecificEnergies",
    "LWDissociatingSpecificEnergies",
    "LWLineTransmissions",
]

EXIT_PASS, EXIT_FAIL, EXIT_MISSING, EXIT_NONFINITE, EXIT_VOID = 0, 1, 2, 3, 4


def wg_alpha(temp: Arr, n_tot: Arr) -> Arr:
    """Return the exponent of the shielding fit, Eq. 7 of Wolcott-Green and Haiman (2019)."""
    tc = np.clip(temp, 1e2, 8e3)
    nc = np.clip(n_tot, 1e-8, 1e7)
    return (0.8711 * np.log10(tc) - 1.928) * np.exp(-0.2856 * np.log10(nc)) + (
        -0.9639 * np.log10(tc) + 3.892
    )


def wg_factor(col: Arr, temp: Arr, n_tot: Arr) -> Arr:
    """Return the shielding factor of the library at the H2 column col [cm^-2]."""
    x = 2e-15 * col
    b5 = 1e-5 * np.sqrt(KB * temp / MH)
    s = 0.965 / (1.0 + x / b5) ** wg_alpha(temp, n_tot) + 0.035 * np.exp(
        -8.5e-4 * np.sqrt(1.0 + x)
    ) / np.sqrt(1.0 + x)
    return np.minimum(s, 1.0)


def kernel_mean_inverse_square(r: Arr, big_h: Arr) -> Arr:
    """Return the mean of 1/|x|^2 over the Wendland C2 kernel of support big_h around a point at r.

    The mean of 1/|x|^2 over a sphere of radius s around a point at distance r from the
    origin is ln((r + s) / |r - s|) / (2 r s). The kernel weight is (1 - q)^4 (1 + 4 q),
    q = s / big_h, on the volume. It needs r > big_h.
    """
    nodes, weights = np.polynomial.legendre.leggauss(64)
    q = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    kernel = (1.0 - q) ** 4 * (1.0 + 4.0 * q) * q**2
    s = q[None, :] * big_h[:, None]
    rr = r[:, None]
    sphere = np.log((rr + s) / np.abs(rr - s)) / (2.0 * rr * s)
    return (sphere * (kernel * w)[None, :]).sum(axis=1) / (kernel * w).sum()


def random_rotation(rng: np.random.Generator) -> Arr:
    """Return a uniform random rotation matrix from a random unit quaternion."""
    quat = rng.normal(size=4)
    a, b, c, d = quat / np.linalg.norm(quat)
    return np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],
        ]
    )


def healpix_nside1(v: Arr) -> NDArray[np.int64]:
    """Return the HEALPix ring pixel (nside 1) of the direction vectors v, shape (N, 3)."""
    v = v / np.linalg.norm(v, axis=1)[:, None]
    z = v[:, 2]
    tt = (np.arctan2(v[:, 1], v[:, 0]) % (2.0 * np.pi)) / (0.5 * np.pi)
    out = np.zeros(len(v), dtype=np.int64)
    eq = np.abs(z) <= 2.0 / 3.0
    t1 = 0.5 + tt[eq]
    t2 = 0.75 * z[eq]
    jp = np.floor(t1 - t2).astype(np.int64)
    jm = np.floor(t1 + t2).astype(np.int64)
    ir = 2 + jp - jm
    ip = ((jp + jm + (1 - (ir & 1))) // 2) % 4
    out[eq] = (ir - 1) * 4 + ip
    pol = ~eq
    tp = tt[pol] - np.floor(tt[pol])
    tmp = np.sqrt(3.0 * (1.0 - np.abs(z[pol])))
    jp = np.floor(tp * tmp).astype(np.int64)
    jm = np.floor((1.0 - tp) * tmp).astype(np.int64)
    ir = jp + jm + 1
    ip = np.floor(tt[pol] * ir).astype(np.int64) % (4 * ir)
    out[pol] = np.where(z[pol] > 0, 2 * ir * (ir - 1) + ip, 12 - 2 * ir * (ir + 1) + ip)
    return out


def hop_matrix(pos: Arr, h: Arr, vol: Arr, star: Arr, power: float) -> tuple:
    """Return the hop matrix Phi (sparse, receiver by sender) and the advance l of each sender.

    pos, h and vol are in pc, pc and pc^3. The light of every sender is radial from the star.
    """
    from scipy.sparse import csr_matrix
    from scipy.spatial import cKDTree

    pairs = cKDTree(pos).query_pairs(GAMMA_KERNEL * h.max(), output_type="ndarray")
    i = np.concatenate([pairs[:, 0], pairs[:, 1]])
    j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    dx = pos[j] - pos[i]
    dist = np.linalg.norm(dx, axis=1)
    reach = (dist < GAMMA_KERNEL * h[i]) | (dist < GAMMA_KERNEL * h[j])
    radial = pos - star
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    mu = np.einsum("ij,ij->i", dx, radial[i]) / dist
    keep = reach & (mu > MU_MIN) & (dist > 0.0)
    i, j, dist, mu = i[keep], j[keep], dist[keep], mu[keep]
    big_h = np.maximum(h[i], h[j])
    q = np.clip(dist / (GAMMA_KERNEL * big_h), 0.0, 1.0)
    weight = vol[j] * (1.0 - q) ** 4 * (1.0 + 4.0 * q) / big_h**3 * mu**power
    total = np.bincount(i, weight, len(pos))
    advance = np.bincount(i, weight * dist * mu, len(pos))
    length = np.where(total > 0.0, advance / np.where(total > 0.0, total, 1.0), 0.0)
    phi = csr_matrix((weight / total[i], (j, i)), shape=(len(pos), len(pos)))
    return phi, length


def shell_mean(vol: Arr, r: Arr, q: Arr, lo: float, hi: float) -> float:
    """Return the volume mean of q over the particles with lo <= r < hi."""
    sel = (r >= lo) & (r < hi)
    return float((vol[sel] * q[sel]).sum() / vol[sel].sum())


def model_window_q(
    pos: Arr,
    vol: Arr,
    star: Arr,
    star_h: float,
    hop: tuple,
    step_of_snap: list[int],
    shells: list[tuple[float, float]],
    seed: int,
) -> Arr:
    """Return the model Q of each shell, averaged over the snapshots of the window.

    One random frame sequence. step_of_snap lists the number of steps done at each snapshot.
    """
    rng = np.random.default_rng(seed)
    phi, length = hop
    offset = pos - star
    r = np.linalg.norm(offset, axis=1)
    support = STAR_EDGE * GAMMA_KERNEL * star_h
    weight = np.where(r < support, vol * (1.0 - r / support), 0.0)
    sel = np.where(weight > 0.0)[0]
    floor = CONE_FLOOR_FRACTION * (np.pi / 3.0) * support**3 / N_CONE
    factor = 4.0 * np.pi * r**2 / vol * length
    energy = np.zeros(len(pos))
    wanted = set(step_of_snap)
    rows = []
    for step in range(1, max(step_of_snap) + 1):
        cell = healpix_nside1(offset[sel] @ random_rotation(rng).T)
        cone = np.bincount(cell, weight[sel], N_CONE)
        escape = np.where(
            floor > cone,
            np.where(cone > 0.0, (1.0 - cone / floor) / N_CONE, 1.0 / N_CONE),
            0.0,
        ).sum()
        emit = np.zeros(len(pos))
        emit[sel] = (1.0 - escape) * weight[sel] / cone.sum()
        energy = phi @ energy + emit
        if step in wanted:
            q = factor * energy
            rows.append([shell_mean(vol, r, q, lo, hi) for lo, hi in shells])
    return np.mean(rows, axis=0)


def bar_shells(s: dict) -> list[tuple[float, float]]:
    """Return the radial shells of the Q bar [pc], fixed multiples of the median gas h."""
    h_med = float(np.median(s["h"])) / PARSEC
    edges = [x * h_med for x in BAR_SHELLS_H]
    return list(zip(edges[:-1], edges[1:]))


def geometry_bar(s: dict, step_of_snap: list[int]) -> dict:
    """Return the model centre and standard deviation of Q in each bar shell.

    s is the last snapshot: its gas is the geometry of the model.
    """
    pos, h = s["pos"] / PARSEC, s["h"] / PARSEC
    vol = s["mass"] / s["rho"] / PARSEC**3
    star = s["star_pos"] / PARSEC
    shells = bar_shells(s)
    hop = hop_matrix(pos, h, vol, star, s["cone_power"])
    seqs = np.array(
        [
            model_window_q(
                pos, vol, star, s["star_h"] / PARSEC, hop, step_of_snap, shells, 100 + k
            )
            for k in range(N_SEQUENCE)
        ]
    )
    return {"centre": seqs.mean(axis=0), "sigma": seqs.std(axis=0, ddof=1)}


def measured_window_q(snaps: list[dict], window: list[int], l_table: float) -> Arr:
    """Return the measured Q of each bar shell, averaged over the snapshots of the window."""
    s = snaps[-1]
    r = np.linalg.norm(s["pos"] - s["star_pos"], axis=1)
    vol = s["mass"] / s["rho"]
    shells = [(lo * PARSEC, hi * PARSEC) for lo, hi in bar_shells(s)]
    rows = []
    for k in window:
        q = (
            C_LIGHT
            * snaps[k]["rho"]
            * snaps[k]["LWOwnSpecificEnergies"]
            * 4.0
            * np.pi
            * r**2
            / l_table
        )
        rows.append([shell_mean(vol, r, q, lo, hi) for lo, hi in shells])
    return np.mean(rows, axis=0)


def bar_void_reasons(s: dict) -> list[str]:
    """Return why the Q geometry bar is not gated for this run (empty if it is).

    The bar is gated only at resolution level BAR_LEVEL and for a star of
    BAR_STAR_MASS_MSUN.
    """
    out = []
    level = np.log2(len(s["mass"])) / 3.0
    if abs(level - BAR_LEVEL) > 1e-9:
        out.append(
            f"{len(s['mass'])} gas particles: level {level:.3g}, not {BAR_LEVEL}"
        )
    mass = s["star_mass"] / MSUN
    if abs(mass - BAR_STAR_MASS_MSUN) > 0.005 * BAR_STAR_MASS_MSUN:
        out.append(f"star mass {mass:.4g} Msun, not {BAR_STAR_MASS_MSUN}")
    return out


def window_steps(snaps: list[dict], steps: Arr) -> tuple[list[int], list[int]]:
    """Return the window snapshot indices and the number of steps done at each of them.

    The first row of timesteps.txt is step 0, the initial state: it is not a step.
    """
    n = len(snaps)
    window = list(range(max(12, n // 3), n - 1, 2))
    tol = 1e-6 * float(steps[-1])
    done = [int((steps[1:] <= snaps[k]["time_code"] + tol).sum()) for k in window]
    return window, done


def read_snapshot(path: Path) -> dict:
    """Read the gas and star fields of one snapshot in cgs."""
    with h5py.File(path, "r") as f:
        um = float(f["Units"].attrs["Unit mass in cgs (U_M)"][0])
        ul = float(f["Units"].attrs["Unit length in cgs (U_L)"][0])
        ut = float(f["Units"].attrs["Unit time in cgs (U_t)"][0])
        g = f["PartType0"]
        order = np.argsort(g["ParticleIDs"][:])
        out: dict = {"id": g["ParticleIDs"][:][order].astype(float)}
        out["pos"] = g["Coordinates"][:][order].astype(float) * ul
        out["mass"] = g["Masses"][:][order].astype(float) * um
        out["h"] = g["SmoothingLengths"][:][order].astype(float) * ul
        out["rho"] = g["Densities"][:][order].astype(float) * um / ul**3
        out["u"] = g["InternalEnergies"][:][order].astype(float) * (ul / ut) ** 2
        for k in SPECIES:
            out[k] = g[k][:][order].astype(float)
        for k in FIELDS:
            out[k] = (
                g[k][:][order].astype(float) * (ul / ut) ** 2
                if "Specific" in k
                else g[k][:][order].astype(float)
            )
        out["z"] = g["SmoothedMetalMassFractions"][:][order][:, -1].astype(float)
        out["time_code"] = float(f["Header"].attrs["Time"][0])
        out["t_unit"] = ut
        out["box"] = float(f["Header"].attrs["BoxSize"][0]) * ul
        s = f["PartType4"] if "PartType4" in f else None
        out["n_star"] = 0 if s is None else int(s["Coordinates"].shape[0])
        if out["n_star"] > 0:
            out["star_pos"] = s["Coordinates"][0].astype(float) * ul
            out["star_mass"] = float(s["Masses"][0]) * um
            out["star_z"] = float(s["MetalMassFractions"][0][-1])
            out["star_L"] = float(s["LWLuminosities"][0]) * um * ul**2 / ut**3
            out["star_h"] = float(s["SmoothingLengths"][0]) * ul
        else:
            out["star_pos"] = np.zeros(3)
            out["star_mass"] = out["star_z"] = out["star_L"] = 1.0
            out["star_h"] = 1.0
        par = f["Parameters"].attrs
        out["interp_size_mass"] = int(
            float(par.get("GEARRadiation:interpolation_size_mass", 0))
        )
        out["solver"] = float(par.get("GeajayCooling:solver", np.nan))
        out["periodic"] = float(par.get("InitialConditions:periodic", np.nan))
        raw = par.get("GEARFeedback:ISRF_packet_light_speed", b"")
        out["light_speed"] = raw.decode() if isinstance(raw, bytes) else str(raw)
        out["cone_power"] = float(
            par.get("GEARFeedback:ISRF_packet_cone_power", np.nan)
        )
    return out


def step_table(path: Path) -> tuple[Arr, Arr]:
    """Return the end time of every step and its star-update count."""
    rows = [
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return np.array([float(r[1]) for r in rows]), np.array([float(r[9]) for r in rows])


def last_update_time(steps: Arr, t_snap: float) -> float:
    """Return the end of the last step at or before a snapshot time (code time).

    Species are not drifted to the snapshot time: they hold the value of the last step end.
    """
    tol = 1e-6 * float(steps[-1])
    return float(steps[steps <= t_snap + tol].max())


def hydrogen_state(s: dict) -> dict:
    """Return x2 and the number density of H nuclei of a snapshot."""
    h_tot = s["HI"] + s["HII"] + s["HM"] + s["H2I"] + s["H2II"]
    n_tot = sum(s[k] / a for k, a in zip(SPECIES, SPECIES_A))
    n_elec = s["HII"] - s["HM"] + s["H2II"] / 2.0
    return {
        "x2": 0.5 * s["H2I"] / h_tot,
        "n_h": h_tot * s["rho"] / MH,
        "n_tot": (n_tot + n_elec) * s["rho"] / MH,
        "hii_ratio": s["HII"] / h_tot,
    }


def jeans_length(temp: Arr, rho: Arr, mu: Arr) -> Arr:
    """Return the Jeans length [cm]."""
    return np.sqrt((5.0 / 3.0) * np.pi * KB * temp / (GRAV * mu * MH * rho))


def table_luminosity(
    path: Path, mass_msun: float, z: float, interp_size: int
) -> tuple[float, float]:
    """Return the LW luminosity [erg/s] that SWIFT reads from the table, and its bar.

    SWIFT resamples the table onto a uniform mass grid from the IMF minimum to the
    IMF maximum mass (Data/IMF Mmin and Mmax) in interpolation_size_mass nodes
    (src/feedback/GEAR/radiation_table_io.c), and interpolates log L bilinearly in
    log Z and log M. A metallicity below the first row of the table (here Z = 0)
    is clamped to the first row, and log10 L is interpolated linearly in log M
    between the two nodes of the resampled grid around the star's mass
    (interpolate_axis_bracket_clamped in src/feedback/GEAR/interpolation.h, with
    edge_interpolates_inside_axis = 1 for the radiation tables). This function
    follows that rule, with the index and its fraction formed in float32 as SWIFT
    does, when the resampled grid equals the grid of the file.

    The bar is 8 float32 ulp of log10 L as a relative error: 8 is a chosen margin
    on the one ulp, eps * |log10 L| * ln 10, which the float32 readout cannot resolve.

    Parameters
    ----------
    path
        The radiation table.
    mass_msun
        The star mass [Msun].
    z
        The star metallicity.
    interp_size
        GEARRadiation:interpolation_size_mass of the run.

    Raises
    ------
    ValueError
        If the mass is outside the table, the metallicity is inside its range or
        the resampled mass grid differs from the grid of the file.
    """
    with h5py.File(path, "r") as f:
        g = f["Data/Radiation"]
        a = g.attrs
        m0, dm, n_mass = float(a["m0"]), float(a["dm"]), int(a["nm"])
        imf = f["Data/IMF"].attrs
        mass_min, mass_max = np.float32(imf["Mmin"]), np.float32(imf["Mmax"])
        z_first = float(g["Metallicity"][0])
        lum = g["L_LW"][0]
    if z >= z_first:
        raise ValueError(
            f"star metallicity {z:g} is not below the first table row {z_first:g}"
        )
    grid_min = np.log10(mass_min)
    grid_step = (np.log10(mass_max) - grid_min) / np.float32(interp_size - 1)
    if (
        interp_size != n_mass
        or abs(grid_min - m0) > 1e-5
        or abs(grid_step - dm) > 1e-5 * dm
    ):
        raise ValueError(
            f"the resampled mass grid ({interp_size} nodes from {float(grid_min):.5f}, "
            f"step {float(grid_step):.6f}) differs from the grid of the file "
            f"({n_mass} nodes from {m0:.5f}, step {dm:.6f})"
        )
    log_m = np.float32(np.log10(mass_msun))
    index = int((log_m - grid_min) / grid_step)
    if index < 0 or index >= n_mass - 1:
        raise ValueError(f"log10 mass {float(log_m):.4f} is outside the table")
    frac = float(np.float32((log_m - grid_min) / grid_step) - np.float32(index))
    value = float(
        10.0
        ** (
            (1.0 - frac) * np.log10(float(lum[index]))
            + frac * np.log10(float(lum[index + 1]))
        )
    )
    bar = 8.0 * float(np.finfo(np.float32).eps) * abs(np.log10(value)) * np.log(10.0)
    return value, bar


def boot_mean(v: Arr, w: Arr, rng: np.random.Generator) -> float:
    """Return the bootstrap standard error of the weighted mean."""
    n = v.size
    idx = rng.integers(0, n, size=(N_BOOT, n))
    means = (v[idx] * w[idx]).sum(axis=1) / w[idx].sum(axis=1)
    return float(np.std(means))


def boot_median(v: Arr, rng: np.random.Generator) -> float:
    """Return the bootstrap standard error of the median."""
    n = v.size
    idx = rng.integers(0, n, size=(N_BOOT, n))
    return float(np.std(np.median(v[idx], axis=1)))


def lit_mask(snaps: list[dict], first: int, last: int) -> Arr:
    """Return the particles whose own LW energy is positive in every snapshot first..last."""
    ok = np.ones(snaps[0]["id"].size, dtype=bool)
    for k in range(first, last + 1):
        ok &= snaps[k]["LWOwnSpecificEnergies"] > 0.0
    return ok


def set_x2(s: dict, x2_new: Arr) -> None:
    """Set the H2 abundance of a snapshot to x2_new and move the mass difference to HI."""
    old = s["H2I"].copy()
    s["H2I"] = old * x2_new / hydrogen_state(s)["x2"]
    s["HI"] = s["HI"] + (old - s["H2I"])


def plant_lw_scale(snaps: list[dict], scale: float) -> None:
    """Scale every LW field by scale and the H2 history as a lowered field would."""
    x0 = hydrogen_state(snaps[0])["x2"]
    for k, s in enumerate(snaps):
        for name in FIELDS[:-1]:
            s[name] = s[name] * scale
        if k > 0:
            set_x2(s, x0 * (hydrogen_state(s)["x2"] / x0) ** scale)


def plant_double_shielding(snaps: list[dict]) -> None:
    """Slow the H2 decay of each interval by the library factor S at the Jeans column."""
    x_plant = hydrogen_state(snaps[0])["x2"].copy()
    for j in range(len(snaps) - 1):
        a = hydrogen_state(snaps[j])
        b = hydrogen_state(snaps[j + 1])
        mu = snaps[j]["rho"] / (a["n_tot"] * MH)
        temp = (2.0 / 3.0) * snaps[j]["u"] * mu * MH / KB
        column = x_plant * a["n_h"] * jeans_length(temp, snaps[j]["rho"], mu)
        factor = wg_factor(column, temp, a["n_tot"])
        x_plant = x_plant * (b["x2"] / a["x2"]) ** factor
        set_x2(snaps[j + 1], x_plant)


def luminosity_stats(s_last: dict, table: Path) -> dict:
    """Return the snapshot luminosity, the table luminosity and the bar."""
    l_table, bar = table_luminosity(
        table,
        s_last["star_mass"] / MSUN,
        s_last["star_z"],
        s_last["interp_size_mass"],
    )
    diff = s_last["star_L"] / l_table - 1.0
    return {
        "snapshot_erg_s": s_last["star_L"],
        "table_erg_s": l_table,
        "relative_difference": diff,
        "bar": bar,
        "ok": abs(diff) <= bar,
    }


def interval_rates(a: dict, b: dict, steps: Arr) -> dict:
    """Return the per-particle rates of the interval between two snapshots."""
    ha, hb = hydrogen_state(a), hydrogen_state(b)
    t_a = last_update_time(steps, a["time_code"])
    t_b = last_update_time(steps, b["time_code"])
    dt = (t_b - t_a) * a["t_unit"]
    u_a, u_b = a["LWDissociatingSpecificEnergies"], b["LWDissociatingSpecificEnergies"]
    with np.errstate(divide="ignore", invalid="ignore"):
        k_meas = -np.log(hb["x2"] / ha["x2"]) / dt
        t_line = 0.5 * (a["LWLineTransmissions"] + b["LWLineTransmissions"])
        rho_mid = 0.5 * (a["rho"] + b["rho"])
        return {
            "dt": dt,
            "k_meas": k_meas,
            "t_line": t_line,
            "cu_eff": k_meas / (KAPPA_D * t_line),
            "ratio4": k_meas / (KAPPA_D * C_LIGHT * rho_mid * 0.5 * (u_a + u_b)),
            "timing": np.abs(np.log(u_b / u_a)),
        }


def field_stats(v: Arr, q: Arr, rng: np.random.Generator) -> dict:
    """Return the volume mean of Q and its bootstrap error in one shell (no bar)."""
    return {
        "q_mean": float((v * q).sum() / v.sum()),
        "q_median": float(np.median(q)),
        "q_se": boot_mean(q, v, rng),
    }


def dissociation_stats(v: Arr, eff: Arr, ref: Arr, rng: np.random.Generator) -> dict:
    """Return the dissociation over its thin limit in one shell (no bar).

    The ratio of sums is the mean of eff / ref weighted by v * ref.
    """
    ratio = eff / ref
    weight = v * ref
    return {
        "d_mean": float((weight * ratio).sum() / weight.sum()),
        "d_se": boot_mean(ratio, weight, rng),
    }


def rate_identity_stats(ratio4: Arr, timing: Arr, rng: np.random.Generator) -> dict:
    """Return the median of ln(ratio4) in one shell and its bar (test 4)."""
    ln4 = np.log(ratio4)
    median = float(np.median(ln4))
    se = boot_median(ln4, rng)
    tim = float(np.median(timing))
    bar = tim + 2.0 * se
    return {
        "r4_median": median,
        "r4_se": se,
        "r4_timing": tim,
        "r4_bar": bar,
        "r4_ok": bool(abs(median) <= bar),
    }


def evaluate(snaps: list[dict], steps: Arr, table: Path, seed: int = 1) -> dict:
    """Return the luminosity, field, dissociation and rate statistics of the run."""
    rng = np.random.default_rng(seed)
    s_last = snaps[-1]
    n_snap = len(snaps)
    lum = luminosity_stats(s_last, table)
    r = np.linalg.norm(s_last["pos"] - s_last["star_pos"], axis=1)
    h_med = float(np.median(s_last["h"]))
    r_hi = 0.5 * s_last["box"] - 2.0 * h_med
    edges = np.arange(3.0 * GAMMA_KERNEL * h_med, r_hi, h_med)
    j_last = n_snap - 2
    qualified_field = lit_mask(snaps, n_snap - N_LIT, n_snap - 1)
    qualified_rate = lit_mask(snaps, j_last - N_LIT + 1, j_last + 1)
    in_domain = (r > 3.0 * GAMMA_KERNEL * s_last["h"]) & (r < r_hi)
    l_table = lum["table_erg_s"]
    vol = s_last["mass"] / s_last["rho"]
    cu_own = C_LIGHT * s_last["rho"] * s_last["LWOwnSpecificEnergies"]
    q_part = cu_own * 4.0 * np.pi * r**2 / l_table
    rates = interval_rates(snaps[j_last], snaps[j_last + 1], steps)
    ref_kernel = (
        l_table
        / (4.0 * np.pi)
        * kernel_mean_inverse_square(r, GAMMA_KERNEL * s_last["h"])
    )
    usable = (
        np.isfinite(rates["k_meas"]) & (rates["k_meas"] > 0.0) & (rates["t_line"] > 0.0)
    )
    shells = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        base = in_domain & (r >= lo) & (r < hi)
        sel_f = base & qualified_field
        sel_r = base & qualified_rate & usable
        row: dict = {
            "r_pc": float(0.5 * (lo + hi) / PARSEC),
            "n_field": int(sel_f.sum()),
            "n_rate": int(sel_r.sum()),
        }
        if sel_f.sum() >= N_MIN:
            row.update(field_stats(vol[sel_f], q_part[sel_f], rng))
        if sel_r.sum() >= N_MIN:
            row.update(
                dissociation_stats(
                    vol[sel_r], rates["cu_eff"][sel_r], ref_kernel[sel_r], rng
                )
            )
            row.update(
                rate_identity_stats(rates["ratio4"][sel_r], rates["timing"][sel_r], rng)
            )
        shells.append(row)
    return {
        "luminosity": lum,
        "shells": shells,
        "dt_s": rates["dt"],
        "interval": [j_last, j_last + 1],
        "n_lit_field": int(qualified_field.sum()),
    }


def finite_report(snaps: list[dict]) -> list[str]:
    """Return one message per non-finite array of any snapshot."""
    out = []
    for i, s in enumerate(snaps):
        for k, v in s.items():
            if (
                isinstance(v, np.ndarray)
                and v.dtype.kind == "f"
                and not np.all(np.isfinite(v))
            ):
                out.append(
                    f"snapshot {i}: {int(np.sum(~np.isfinite(v)))} non-finite values in {k}"
                )
    return out


def void_reasons(snaps: list[dict], steps: Arr, star_updates: Arr) -> list[str]:
    """Return the reasons why the run is void for this check."""
    out = []
    s = snaps[-1]
    if s["n_star"] != 1:
        out.append(f"{s['n_star']} stars, need 1")
    if s["solver"] != 1.0:
        out.append(f"GeajayCooling:solver is {s['solver']}, not 1")
    if s["periodic"] != 0.0:
        out.append("the box is periodic: images of the star add light")
    if s["light_speed"] != "unlimited":
        out.append(f"light speed is {s['light_speed']!r}, not 'unlimited'")
    if float(np.max(s["z"])) > 0.0:
        out.append("the gas has metals: dust absorbs light")
    h = hydrogen_state(s)
    if float(np.max(h["hii_ratio"])) > 1e-6:
        out.append("the gas is ionised above 1e-6: H- formation is not negligible")
    if np.any(star_updates[1:] < 1):
        out.append("the star does not update at every gas step")
    return out


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("snapdir", type=Path)
    p.add_argument(
        "--radiation-table",
        type=Path,
        required=True,
        help="the radiation table of the run",
    )
    p.add_argument("--timesteps", type=Path, default=None)
    p.add_argument(
        "--plant-lw-scale",
        type=float,
        default=None,
        help="plant: scale the LW fields and the H2 rate",
    )
    p.add_argument(
        "--plant-double-shielding",
        action="store_true",
        help="plant: shield the transported rate once more",
    )
    p.add_argument(
        "--plant-luminosity",
        type=float,
        default=None,
        help="plant: scale the star luminosity of the snapshot",
    )
    p.add_argument(
        "--plant-interp-size",
        type=int,
        default=None,
        help="plant: set GEARRadiation:interpolation_size_mass of the snapshot",
    )
    p.add_argument(
        "--plant-star-skip",
        action="store_true",
        help="plant: the star makes no update in one step",
    )
    p.add_argument(
        "--plant-q-scale",
        type=float,
        default=None,
        help="plant: scale the stored LW energy that gives Q (nan: make it not finite)",
    )
    p.add_argument(
        "--plant-nan",
        action="store_true",
        help="plant: one NaN in the density of the last snapshot",
    )
    p.add_argument("--json", type=Path, default=None)
    o = p.parse_args()

    files = sorted(o.snapdir.glob("snapshot_*.hdf5"))
    timesteps = o.timesteps or (o.snapdir.parent / "timesteps.txt")
    if len(files) < 12 or not timesteps.exists() or not o.radiation_table.exists():
        print(
            f"MISSING: {len(files)} snapshots, timesteps {timesteps.exists()}, table {o.radiation_table.exists()}"
        )
        return EXIT_MISSING
    try:
        snaps = [read_snapshot(f) for f in files]
    except (KeyError, OSError) as err:
        print(f"MISSING: {err}")
        return EXIT_MISSING
    steps, star_updates = step_table(timesteps)
    keep = [0] + [
        k
        for k in range(1, len(snaps))
        if snaps[k]["time_code"] > snaps[k - 1]["time_code"] * (1.0 + 1e-9)
    ]
    # The snapshot written at the end of the run holds the species of one step earlier than its time says.
    snaps = [snaps[k] for k in keep][:-1]
    if o.plant_nan:
        snaps[-1]["rho"][0] = np.nan
    bad = finite_report(snaps)
    if bad:
        print("NON-FINITE:\n  " + "\n  ".join(bad))
        return EXIT_NONFINITE
    if o.plant_interp_size is not None:
        snaps[-1]["interp_size_mass"] = o.plant_interp_size
    if o.plant_star_skip:
        star_updates[len(star_updates) // 2] = 0.0
    void = void_reasons(snaps, steps, star_updates)
    if void:
        print("VOID:\n  " + "\n  ".join(void))
        return EXIT_VOID
    if o.plant_lw_scale is not None:
        plant_lw_scale(snaps, o.plant_lw_scale)
        print(f"PLANT: LW fields and H2 rate x {o.plant_lw_scale}")
    if o.plant_double_shielding:
        plant_double_shielding(snaps)
        print("PLANT: transported rate shielded a second time at the Jeans column")

    if o.plant_q_scale is not None:
        for s in snaps:
            s["LWOwnSpecificEnergies"] = s["LWOwnSpecificEnergies"] * o.plant_q_scale
        print(f"PLANT: stored LW energy (Q) x {o.plant_q_scale}")

    if o.plant_luminosity is not None:
        snaps[-1]["star_L"] *= o.plant_luminosity
        print(f"PLANT: star luminosity x {o.plant_luminosity}")

    try:
        res = evaluate(snaps, steps, o.radiation_table)
    except ValueError as err:
        print(f"VOID: {err}")
        return EXIT_VOID
    lum = res["luminosity"]
    for key in ("snapshot_erg_s", "table_erg_s", "relative_difference", "bar"):
        if not np.isfinite(lum[key]):
            print(f"NON-FINITE luminosity {key}")
            return EXIT_NONFINITE
    print(
        f"luminosity: snapshot {lum['snapshot_erg_s']:.4e} erg/s, table {lum['table_erg_s']:.4e} erg/s, "
        f"difference {lum['relative_difference']:+.2e}, bar {lum['bar']:.2e}  {'ok' if lum['ok'] else 'FAIL'}"
    )
    print(
        f"interval {res['interval']}, dt {res['dt_s'] / 3.15576e7:.4f} yr, lit particles {res['n_lit_field']}"
    )
    print(
        "Q and 2b are reported without a bar here; the luminosity, the rate identity (4) and the Q geometry bar gate."
    )
    print(
        f"{'r[pc]':>6} {'N':>6} {'Q mean':>8} {'Q med':>8} {'se':>7} | {'2b':>8} {'se':>7} | {'ln r4':>10} {'bar':>9}"
    )
    status = EXIT_PASS if lum["ok"] else EXIT_FAIL
    void_q = bar_void_reasons(snaps[-1])
    if not void_q:
        window, done = window_steps(snaps, steps)
        bar = geometry_bar(snaps[-1], done)
        q_win = measured_window_q(snaps, window, lum["table_erg_s"])
        if not (np.all(np.isfinite(q_win)) and np.all(np.isfinite(bar["centre"]))):
            print("NON-FINITE geometry bar or measured Q")
            return EXIT_NONFINITE
    n_tested = 0
    for sh in res["shells"]:
        cells = [f"{sh['r_pc']:6.2f}", f"{sh['n_field']:6d}"]
        for key in (
            "q_mean",
            "q_median",
            "q_se",
            "d_mean",
            "d_se",
            "r4_median",
            "r4_bar",
        ):
            if key in sh and not np.isfinite(sh[key]):
                print(f"NON-FINITE {key} in shell {sh['r_pc']:.2f}")
                return EXIT_NONFINITE
        fail = []
        if "r4_ok" in sh:
            n_tested += 1
            if not sh["r4_ok"]:
                fail.append("4")
        if fail:
            status = EXIT_FAIL
        q = (
            f"{sh['q_mean']:8.4f} {sh['q_median']:8.4f} {sh['q_se']:7.4f}"
            if "q_mean" in sh
            else " " * 25
        )
        d = f"{sh['d_mean']:8.4f} {sh['d_se']:7.4f}" if "d_mean" in sh else " " * 16
        r4 = (
            f"{sh['r4_median']:10.4e} {sh['r4_bar']:9.3e}"
            if "r4_ok" in sh
            else " " * 20
        )
        print(
            f"{' '.join(cells)} {q} | {d} | {r4}  {'FAIL ' + ','.join(fail) if fail else 'ok'}"
        )
    if n_tested == 0:
        print("VOID: no shell holds enough lit particles")
        void_q = void_q + ["no shell holds enough lit particles"]
    if void_q:
        print("VOID: the Q geometry bar is not gated here:\n  " + "\n  ".join(void_q))
    else:
        print(
            f"Q geometry bar: window snapshots {window[0]}..{window[-1]} step 2, cone power "
            f"{snaps[-1]['cone_power']:g}, model {N_SEQUENCE} frame sequences, bar centre +- {BAR_SIGMA:g} sigma"
        )
        print(
            f"{'r[pc]':>10} {'Q window':>9} {'centre':>8} {'sigma':>7} {'(Q-c)/sig':>10}"
        )
        if not (np.all(np.isfinite(bar["sigma"])) and np.all(bar["sigma"] > 0.0)):
            print("NON-FINITE or not positive sigma of the geometry bar")
            return EXIT_NONFINITE
        for (lo, hi), q, c, sg in zip(
            bar_shells(snaps[-1]), q_win, bar["centre"], bar["sigma"]
        ):
            z = (q - c) / sg
            bad = abs(z) > BAR_SIGMA
            if bad:
                status = EXIT_FAIL
            print(
                f"{lo:4.2f}-{hi:4.2f} {q:9.4f} {c:8.4f} {sg:7.4f} {z:10.2f}  {'FAIL Q' if bad else 'ok'}"
            )
    if status == EXIT_PASS and void_q:
        status = EXIT_VOID
    print("RESULT:", {EXIT_PASS: "PASS", EXIT_VOID: "VOID"}.get(status, "FAIL"))
    if o.json:
        o.json.write_text(json.dumps({"exit": status, **res}, indent=1, default=float))
    return status


if __name__ == "__main__":
    sys.exit(main())
