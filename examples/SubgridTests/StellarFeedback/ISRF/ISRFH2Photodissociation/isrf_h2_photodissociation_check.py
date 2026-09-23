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
"""Check the H2 photodissociation rate the LW band drives through Grackle.

Compares the measured H2 abundance history of every gas particle against the
shielded photodissociation rate predicted from that particle's own LW energy
density, temperature, density and H2 abundance. The prediction chain is
reproduced here term by term, in CGS throughout.

Step 1: the rate the LW band hands to Grackle
---------------------------------------------
``src/feedback/GEAR/radiation.h`` converts the propagated LW specific energy
into a photodissociation rate as a band-averaged cross-section times the LW
photon flux,

    k_diss,0 = sigma_H2 * c * rho * u_LW / E_LW ,                        (1)

with ``sigma_H2 = RADIATION_SIGMA_H2_LW_CGS``,
``E_LW = RADIATION_LW_PHOTON_ENERGY_EV`` (both read from that header at run
time, since only their quotient is calibrated) and ``u_LW`` the snapshot's
``LWSpecificEnergies`` in CGS. This script takes ``u_LW`` from the snapshots
rather than predicting it from a transport solution: the field the chemistry
actually saw is an input to this test, not one of its claims, so a transport
amplitude error cannot make this check pass or fail.

Equation (1) is the convention this code uses. The literature convention is a
free-space rate per Habing field, ``k_diss,0 = 3.3e-11 * G_0 s^-1``
(Draine and Bertoldi 1996, their Table 2 unshielded rate, with
``G_0 = c * rho * (u_PE + u_LW) / 1.6e-3 erg s^-1 cm^-2``, Habing 1968).
The two differ because (1) counts only the LW band while ``G_0`` is the
band-summed 6 to 13.6 eV field. Both are reported; only (1) is gated.

Step 2: the shielding factor Grackle applies
--------------------------------------------
Grackle multiplies the supplied rate by a self-shielding factor
(``solve_rate_cool_g.F``, the ``iH2shield > 0`` block). The functional form is
that of Draine and Bertoldi (1996), their Eq. 37,

    f_shield = 0.965 / (1 + x/b5)**a
             + 0.035 * exp(-8.5e-4 * sqrt(1 + x)) / sqrt(1 + x) ,        (2)
    x  = N_H2 / 5e14 cm^-2 ,
    b5 = sqrt(2 k_B T / m_H) / (1e5 cm s^-1) ,

capped at 1. Draine and Bertoldi (1996) use ``a = 2``; Wolcott-Green,
Haiman and Bryan (2011) revised it to ``a = 1.1``; the version Grackle runs
uses the temperature- and density-dependent exponent of Wolcott-Green and
Haiman (2019),

    a = (0.8711 log10(T') - 1.928) * exp(-0.2856 log10(n'))
        + (-0.9639 log10(T') + 3.892) ,                                  (3)

with ``T'`` clamped to [100, 8000] K and ``n'`` the total number density
clamped to 1e7 cm^-3. Note that ``b5`` uses the unclamped temperature while
(3) uses the clamped one, as in the Fortran.

Step 3: the column the shielding factor is evaluated on
-------------------------------------------------------
Grackle does not integrate a column along a line of sight. It forms a local
estimate from the H2 density and a shielding length selected by
``GrackleCooling:H2_self_shielding``,

    N_H2 = 2 * n_H2 * l_shield .                                         (4)

The factor 2 is Grackle's species-density convention: its ``H2I`` field is a
mass density, so ``H2I / m_H`` is ``2 n_H2``, not ``n_H2``. SWIFT stores the
same quantity as the mass fraction ``H2I`` in the snapshots, so
``n_H2 = H2I * rho / (2 m_H)``.

Mode 3 uses the local Jeans length,

    l_shield = sqrt(gamma * pi * k_B * T / (G * mu * m_H * rho)) ,        (5)

and mode 2 uses the length SWIFT supplies per particle, the kernel support
radius

    l_shield = gamma_K * h ,                                             (5b)

with ``h`` the snapshot's ``SmoothingLengths`` and ``gamma_K`` the kernel's
support-to-smoothing ratio (``--kernel-gamma``, 1.936492 for the Wendland C2
kernel in 3D). With Eq. (4), mode 2 is the H2 column through a path of
``2 gamma_K h``. With ``GrackleCooling:H2_self_shielding_path:
kernel_radius`` (``--h2-self-shielding-path``) SWIFT supplies half that
length, a path of ``gamma_K h``. SWIFT rejects mode 1 (Sobolev-like) at
start-up: it reads six neighbouring grid points that do not exist when
Grackle is called on one particle. Mode 0 disables shielding altogether,
``f_shield = 1``, and is the unshielded reference.

The mean molecular weight follows Grackle's own definition
(``cool1d_multi_g.F``),

    1/mu = HI + HII + e + (HeI + HeII + HeIII)/4
           + HM + (H2I + H2II)/2 + Z/16 ,                                (6)

on the snapshot mass fractions, and the temperature from
``T = (gamma - 1) u mu m_H / k_B``. Grackle applies a small correction to
``T`` for the H2 rotational degrees of freedom; at the H2 abundances used
here it is below the pass bars and is not reproduced.

Step 4: what is compared
------------------------
With formation suppressed (zero metallicity removes dust, and the initial
ionization is set low enough that the H- route is negligible; both are
verified below from the run's own numbers), each particle obeys

    d ln x_H2 / dt = -k_diss,0 * f_shield ,                              (7)

so between the first and last snapshot of the measurement window

    ln[x_H2(t_1) / x_H2(t_0)] = -integral k_diss,0 f_shield dt .          (8)

The right-hand side is integrated with the trapezoid rule over the snapshot
times, using each snapshot's own ``k_diss,0`` and ``f_shield``. Each particle
gets its own window. It opens at the first snapshot at which the LW field has
reached that particle (the field starts at zero and propagates out from the
star, so arrival times differ by several snapshots across the box). In the
optically thin configuration it closes at the last snapshot at which ``x_H2``
is still ``--floor-margin`` (default 100) times above the lowest value the
particle reaches in the run: at that point the residual formation rate that
eventually balances destruction is below 1 percent of it, so Eq. (7) holds
to better than the pass bar. In the shielded configuration ``x_H2`` drops by
about a percent over the run and the window closes at the last snapshot.

Pass bars
---------
Both bars come from the expected numerical error of the comparison, not from
what the run happens to produce.

* Optically thin configuration: for each particle, take the first snapshot
  at which the predicted integral of Eq. (8) reaches ``--thin-efolds``
  (default 2). The median, over particles, of the relative error on
  ``x_H2`` there must be at or below ``--thin-tol`` (default 0.08), and
  ``f_shield`` must stay at or above ``--thin-f-min`` (default 0.90) so the
  leg really is the unshielded limit. Comparing at a fixed number of
  predicted e-folds keeps the bar independent of the window length: a
  fractional rate error ``delta`` becomes an ``x_H2`` error of about
  ``exp(N delta) - 1`` after ``N`` e-folds. The rate error budget, summed
  linearly: the trapezoid rule on a field that changes by tens of percent
  per snapshot, about 1e-2; a lag of up to one gas step between the field
  the chemistry used and the field the snapshot records, about 3e-2 at
  this run's field growth rate; float32 storage, 1e-5. That is 4e-2, and
  ``exp(2 * 0.04) - 1 = 0.0833`` at two e-folds, hence 0.08. The median
  ratio of measured to predicted e-folds over each full window is printed
  as well.

* Self-shielded configuration: the measured shielding ratio

      r = -ln[x_H2(t_1)/x_H2(t_0)] / integral k_diss,0 dt                (9)

  must agree with the predicted, ``f_shield``-weighted equivalent
  ``integral k_diss,0 f_shield dt / integral k_diss,0 dt``, particle by
  particle, with the median of the per-particle quotient within
  ``--thick-tol`` (default 0.25) of 1, and ``f_shield`` must stay at or
  below ``--thick-f-max`` (default 1e-2) so the leg really is strongly
  shielded. The bar is looser than the thin one because ``r`` is a ratio of
  a small ln-drop to a large integral: at ``f_shield ~ 1e-3`` the measured
  ln-drop is of order 1e-2, so the same absolute few-times-1e-2
  trapezoid/lag error (see the thin budget above) is a larger relative
  error on the numerator.

Rate normalisation (reported, never gated)
------------------------------------------

Eq. (1) over ``3.3e-11 G_0`` is ``6.229 * u_LW / (u_PE + u_LW)`` for a mean
LW photon energy of ``E_LW``: the constant is
``sigma_H2 * 1.6e-3 erg s^-1 cm^-2 / (E_LW * 3.3e-11 s^-1)``. It is 1 only
for an LW fraction of 0.161. The Draine (1978) field, which the
``3.3e-11 G_0`` rate assumes, has an LW fraction of 0.148 (ratio 0.92). The
script prints the measured ratio and ``6.229`` times the measured LW
fraction; ``draine_spectrum`` uses a star with the Draine LW fraction.

Mode-choice diagnostic (reported, never gated)
----------------------------------------------
The quantity that decides which ``H2_self_shielding`` mode to use is not a
pass bar but the ratio of the mode's own column, Eq. (4), to the H2 column
actually present between the particle and the star,

    N_H2,geom(r) = integral_0^r n_H2 dl ,                               (10)

evaluated on the run's own radial H2 profile. In a uniform medium the ratio
is ``2 l_shield / r``: the factor 2 is the convention of Eq. (4), the rest is
how far the local length is from the true path. Printed per radial shell,
alongside the shielding length itself, for a shielding mode only; it is
meaningful while the H2 profile is still close to uniform, as in the
self-shielded configuration, and is not printed for mode 0.
"""

import argparse
import glob
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RADIATION_H = (
    Path(__file__).resolve().parents[5] / "src" / "feedback" / "GEAR" / "radiation.h"
)


def read_radiation_h_constant(name: str) -> float:
    """Read a #define'd float constant's value out of radiation.h.

    Parameters
    ----------
    name : str
        The macro name, e.g. ``"RADIATION_SIGMA_H2_LW_CGS"``.

    Returns
    -------
    float
        The macro's value.
    """
    text = RADIATION_H.read_text()
    match = re.search(rf"^#define\s+{re.escape(name)}\s+([0-9.eE+-]+)", text, re.M)
    if match is None:
        raise ValueError(f"Could not find #define {name} in {RADIATION_H}")
    return float(match.group(1))


# The cross-section and the photon energy are calibrated together: only their
# quotient is constrained, so a copy of one that drifts from the other
# rescales every dissociation rate. Read both from the header instead.
SIGMA_H2_LW_CGS: float = read_radiation_h_constant("RADIATION_SIGMA_H2_LW_CGS")
LW_PHOTON_ENERGY_EV: float = read_radiation_h_constant("RADIATION_LW_PHOTON_ENERGY_EV")
HABING_FLUX_CGS: float = 1.6e-3
# Draine and Bertoldi (1996), unshielded free-space rate per Habing field
DB96_UNSHIELDED_RATE_CGS: float = 3.3e-11
# src/physical_constants_cgs.h
C_LIGHT_CGS: float = 2.99792458e10
K_BOLTZMANN_CGS: float = 1.380649e-16
M_H_CGS: float = 1.67262192369e-24
NEWTON_G_CGS: float = 6.67430e-8
ELECTRON_VOLT_CGS: float = 1.602176634e-12
PARSEC_CGS: float = 3.0856775814913673e18
# Grackle's H2 column normalisation, Eq. (2)
N_H2_NORM_CGS: float = 5.0e14
# Wendland C2 support-to-smoothing ratio in 3D, src/kernel_hydro.h
KERNEL_GAMMA_WENDLAND_C2: float = 1.936492
# Eq. (1) over 3.3e-11 G_0, per unit LW fraction
RATE_RATIO_PER_LW_FRACTION: float = (
    SIGMA_H2_LW_CGS
    * HABING_FLUX_CGS
    / (LW_PHOTON_ENERGY_EV * ELECTRON_VOLT_CGS * DB96_UNSHIELDED_RATE_CGS)
)
# Mean atomic weight Grackle assigns to the metal field, cool1d_multi_g.F
MU_METAL: float = 16.0
GAMMA: float = 5.0 / 3.0


def parse_options() -> argparse.Namespace:
    """Parse command-line options.

    Returns
    -------
    argparse.Namespace
        The parsed options.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to read (default: %(default)s)",
    )
    parser.add_argument(
        "--config",
        choices=["thin", "thick", "draine_spectrum"],
        required=True,
        help="Configuration the run was launched with; selects the pass bar "
        "(draine_spectrum uses the thin bar)",
    )
    parser.add_argument(
        "--h2-self-shielding",
        type=int,
        default=3,
        help="GrackleCooling:H2_self_shielding the run used: 0, 2 or 3 "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--kernel-gamma",
        type=float,
        default=KERNEL_GAMMA_WENDLAND_C2,
        help="Kernel support-to-smoothing ratio, mode 2 only "
        "(default: %(default)s, Wendland C2 in 3D)",
    )
    parser.add_argument(
        "--h2-self-shielding-path",
        choices=["kernel_diameter", "kernel_radius"],
        default="kernel_diameter",
        help="GrackleCooling:H2_self_shielding_path the run used, mode 2 only "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--thin-tol",
        type=float,
        default=0.08,
        help="Max median relative error on x_H2 at the end of the window, "
        "thin configuration (default: %(default)s)",
    )
    parser.add_argument(
        "--thin-efolds",
        type=float,
        default=2.0,
        help="Thin configuration: predicted e-folds at which x_H2 is compared "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--thin-f-min",
        type=float,
        default=0.90,
        help="Min f_shield for the thin configuration to count as unshielded "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--thick-tol",
        type=float,
        default=0.25,
        help="Max relative error on the measured shielding ratio, thick "
        "configuration (default: %(default)s)",
    )
    parser.add_argument(
        "--thick-f-max",
        type=float,
        default=1e-2,
        help="Max f_shield for the thick configuration to count as strongly "
        "shielded (default: %(default)s)",
    )
    parser.add_argument(
        "--field-arrival-fraction",
        type=float,
        default=0.05,
        help="A particle's LW field counts as arrived once k_diss,0 exceeds "
        "this fraction of its own maximum (default: %(default)s)",
    )
    parser.add_argument(
        "--floor-margin",
        type=float,
        default=100.0,
        help="Thin configuration: close each particle's window while x_H2 is "
        "still this many times above the lowest value it reaches "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--n-shells",
        type=int,
        default=6,
        help="Number of radial shells for the printed profiles "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default="isrf_h2_photodissociation_check.png",
        help="Output figure name (default: %(default)s)",
    )
    return parser.parse_args()


def read_snapshot(filename: str) -> Dict[str, np.ndarray]:
    """Read the gas fields this check needs from one snapshot.

    Parameters
    ----------
    filename : str
        Path to the snapshot.

    Returns
    -------
    dict of str to numpy.ndarray
        Physical CGS gas quantities, plus the scalar snapshot time in seconds
        under key ``time`` and the star position in cm under ``star_position``.
    """
    with h5py.File(filename, "r") as handle:
        units = handle["/Units"].attrs
        length_cgs = float(np.atleast_1d(units["Unit length in cgs (U_L)"]).ravel()[0])
        mass_cgs = float(np.atleast_1d(units["Unit mass in cgs (U_M)"]).ravel()[0])
        time_cgs = float(np.atleast_1d(units["Unit time in cgs (U_t)"]).ravel()[0])
        gas = handle["/PartType0"]
        names = [
            "Coordinates",
            "Densities",
            "InternalEnergies",
            "LWSpecificEnergies",
            "PESpecificEnergies",
            "SmoothingLengths",
            "ParticleIDs",
            "HI",
            "HII",
            "HeI",
            "HeII",
            "HeIII",
            "e",
            "HM",
            "H2I",
            "H2II",
        ]
        raw = {name: gas[name][:].astype(np.float64) for name in names}
        metals = gas["MetalMassFractions"][:].astype(np.float64)
        named_columns = "/SubgridScheme/NamedColumns/MetalMassFractions"
        metal_names = (
            [name.decode() for name in handle[named_columns][:]]
            if named_columns in handle
            else []
        )
        star_position = handle["/PartType4/Coordinates"][:].astype(np.float64)
        snapshot_time = float(np.atleast_1d(handle["/Header"].attrs["Time"]).ravel()[0])

    energy_cgs = (length_cgs / time_cgs) ** 2
    # The per-element columns are part of the total in the "Metals" column.
    if metals.ndim == 2:
        metallicity = metals[:, metal_names.index("Metals")]
    else:
        metallicity = metals

    density = raw["Densities"] * mass_cgs / length_cgs**3
    # Grackle's own mean molecular weight, Eq. (6).
    inverse_mu = (
        raw["HI"]
        + raw["HII"]
        + raw["e"]
        + 0.25 * (raw["HeI"] + raw["HeII"] + raw["HeIII"])
        + raw["HM"]
        + 0.5 * (raw["H2I"] + raw["H2II"])
        + metallicity / MU_METAL
    )
    mu = 1.0 / inverse_mu
    temperature = (
        (GAMMA - 1.0)
        * raw["InternalEnergies"]
        * energy_cgs
        * mu
        * M_H_CGS
        / K_BOLTZMANN_CGS
    )

    return {
        "time": snapshot_time * time_cgs,
        "star_position": star_position[0] * length_cgs,
        "position": raw["Coordinates"] * length_cgs,
        "smoothing_length": raw["SmoothingLengths"] * length_cgs,
        "ids": raw["ParticleIDs"],
        "density": density,
        "mu": mu,
        "temperature": temperature,
        "u_LW": raw["LWSpecificEnergies"] * energy_cgs,
        "u_PE": raw["PESpecificEnergies"] * energy_cgs,
        "H2I_fraction": raw["H2I"],
        "n_H2": raw["H2I"] * density / (2.0 * M_H_CGS),
    }


def unshielded_rate(density: np.ndarray, u_LW: np.ndarray) -> np.ndarray:
    """Compute the unshielded H2 photodissociation rate, Eq. (1).

    Parameters
    ----------
    density : numpy.ndarray
        Gas mass density, g cm^-3.
    u_LW : numpy.ndarray
        LW-band specific energy, erg g^-1.

    Returns
    -------
    numpy.ndarray
        Rate in s^-1.
    """
    photon_energy = LW_PHOTON_ENERGY_EV * ELECTRON_VOLT_CGS
    return SIGMA_H2_LW_CGS * C_LIGHT_CGS * density * u_LW / photon_energy


def habing_field(density: np.ndarray, u_PE: np.ndarray, u_LW: np.ndarray) -> np.ndarray:
    """Compute the band-summed field strength in Habing units.

    Parameters
    ----------
    density : numpy.ndarray
        Gas mass density, g cm^-3.
    u_PE : numpy.ndarray
        PE-band specific energy, erg g^-1.
    u_LW : numpy.ndarray
        LW-band specific energy, erg g^-1.

    Returns
    -------
    numpy.ndarray
        G_0, dimensionless.
    """
    return C_LIGHT_CGS * density * (u_PE + u_LW) / HABING_FLUX_CGS


def jeans_shielding_length(
    temperature: np.ndarray, density: np.ndarray, mu: np.ndarray
) -> np.ndarray:
    """Compute Grackle's mode-3 shielding length, Eq. (5).

    Parameters
    ----------
    temperature : numpy.ndarray
        Gas temperature, K.
    density : numpy.ndarray
        Gas mass density, g cm^-3.
    mu : numpy.ndarray
        Mean molecular weight.

    Returns
    -------
    numpy.ndarray
        Local Jeans length in cm.
    """
    return np.sqrt(
        GAMMA
        * np.pi
        * K_BOLTZMANN_CGS
        * temperature
        / (NEWTON_G_CGS * mu * M_H_CGS * density)
    )


def shielding_factor(
    column: np.ndarray,
    temperature: np.ndarray,
    number_density: np.ndarray,
) -> np.ndarray:
    """Compute Grackle's H2 self-shielding factor, Eqs. (2) and (3).

    Parameters
    ----------
    column : numpy.ndarray
        H2 column as Grackle forms it, Eq. (4), cm^-2.
    temperature : numpy.ndarray
        Gas temperature, K.
    number_density : numpy.ndarray
        Total particle number density, cm^-3.

    Returns
    -------
    numpy.ndarray
        f_shield, capped at 1.
    """
    clamped_temperature = np.clip(temperature, 1e2, 8e3)
    clamped_density = np.minimum(number_density, 1e7)
    exponent = (0.8711 * np.log10(clamped_temperature) - 1.928) * np.exp(
        -0.2856 * np.log10(clamped_density)
    ) + (-0.9639 * np.log10(clamped_temperature) + 3.892)
    x = column / N_H2_NORM_CGS
    b5 = np.sqrt(2.0 * K_BOLTZMANN_CGS * temperature / M_H_CGS) / 1e5
    factor = 0.965 / (1.0 + x / b5) ** exponent + 0.035 * np.exp(
        -8.5e-4 * np.sqrt(1.0 + x)
    ) / np.sqrt(1.0 + x)
    return np.minimum(factor, 1.0)


def build_history(
    filenames: List[str],
    self_shielding_mode: int,
    kernel_gamma: float,
    path_in_kernel_radii: float,
) -> Dict[str, np.ndarray]:
    """Assemble the per-particle time series this check compares.

    Particles are matched across snapshots by ``ParticleIDs``, so a changed
    write order cannot silently scramble the histories.

    Parameters
    ----------
    filenames : list of str
        Snapshot paths in time order.
    self_shielding_mode : int
        ``GrackleCooling:H2_self_shielding`` the run used.
    kernel_gamma : float
        Kernel support-to-smoothing ratio, for the mode-2 length.
    path_in_kernel_radii : float
        Mode-2 H2 column path in kernel support radii: 2 for
        ``kernel_diameter``, 1 for ``kernel_radius``.

    Returns
    -------
    dict of str to numpy.ndarray
        Arrays of shape ``(n_snapshots, n_particles)`` for ``x_H2``,
        ``rate``, ``f_shield``, ``column``, ``habing``, ``length``,
        ``temperature`` and ``lw_fraction``, plus ``time`` of shape ``(n_snapshots,)``,
        ``radius`` and ``n_H2`` of shape ``(n_particles,)`` taken at the
        first snapshot.
    """
    first = read_snapshot(filenames[0])
    order = np.argsort(first["ids"])
    reference_ids = first["ids"][order]

    series: Dict[str, List[np.ndarray]] = {
        key: []
        for key in (
            "x_H2",
            "rate",
            "f_shield",
            "column",
            "habing",
            "length",
            "temperature",
            "lw_fraction",
        )
    }
    times: List[float] = []

    for filename in filenames:
        snapshot = read_snapshot(filename)
        index = np.argsort(snapshot["ids"])
        if not np.array_equal(snapshot["ids"][index], reference_ids):
            raise RuntimeError(f"{filename} has a different particle set")
        take = lambda name: snapshot[name][index]

        density = take("density")
        temperature = take("temperature")
        mu = take("mu")
        n_H2 = take("n_H2")
        if self_shielding_mode == 2:
            length = (
                0.5 * path_in_kernel_radii * kernel_gamma * take("smoothing_length")
            )
        elif self_shielding_mode in (0, 3):
            length = jeans_shielding_length(temperature, density, mu)
        else:
            raise RuntimeError(
                f"H2_self_shielding={self_shielding_mode} is not a mode SWIFT "
                "runs. Use 0, 2 or 3."
            )
        column = 2.0 * n_H2 * length
        number_density = density / (mu * M_H_CGS)
        if self_shielding_mode == 0:
            factor = np.ones_like(column)
        else:
            factor = shielding_factor(column, temperature, number_density)
        u_LW = take("u_LW")
        u_total = take("u_PE") + u_LW

        series["x_H2"].append(take("H2I_fraction"))
        series["rate"].append(unshielded_rate(density, take("u_LW")))
        series["f_shield"].append(factor)
        series["column"].append(column)
        series["habing"].append(habing_field(density, take("u_PE"), take("u_LW")))
        series["length"].append(length)
        series["temperature"].append(temperature)
        series["lw_fraction"].append(
            np.divide(u_LW, u_total, out=np.zeros_like(u_LW), where=u_total > 0.0)
        )
        times.append(snapshot["time"])

    radius = np.linalg.norm(first["position"][order] - first["star_position"], axis=1)
    history = {key: np.array(value) for key, value in series.items()}
    history["time"] = np.array(times)
    history["radius"] = radius
    history["n_H2_initial"] = first["n_H2"][order]
    return history


def particle_windows(
    history: Dict[str, np.ndarray],
    arrival_fraction: float,
    floor_margin: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Choose, per particle, the snapshot range over which Eq. (8) holds.

    Parameters
    ----------
    history : dict of str to numpy.ndarray
        Output of :func:`build_history`.
    arrival_fraction : float
        Fraction of the particle's own peak rate that counts as field arrival.
    floor_margin : float
        Factor above the particle's lowest x_H2 at which the window closes;
        0 keeps the window open to the last snapshot.

    Returns
    -------
    tuple of numpy.ndarray
        First and last snapshot index of each particle's window, inclusive,
        and a boolean mask of the particles the field reached at all.
    """
    rate = history["rate"]
    x_H2 = history["x_H2"]
    n_snapshots = x_H2.shape[0]

    peak = rate.max(axis=0)
    arrived = rate >= arrival_fraction * peak[np.newaxis, :]
    reached = arrived.any(axis=0) & (peak > 0.0)
    start = np.argmax(arrived, axis=0)

    lowest = x_H2.min(axis=0)
    above_floor = x_H2 > floor_margin * lowest[np.newaxis, :]
    after_start = np.arange(n_snapshots)[:, np.newaxis] >= start[np.newaxis, :]
    usable = above_floor & after_start
    last_usable = n_snapshots - 1 - np.argmax(usable[::-1], axis=0)
    end = np.where(usable.any(axis=0), last_usable, start)
    return start, end, reached


def check_gate(name: str, value: float, bad: bool, message: str) -> Optional[str]:
    """Build a failure string for one pass/fail gate, failing closed on NaN/inf.

    A gate written as a bare ``value > bar`` (or ``< bar``) comparison lets a
    non-finite ``value`` slip through silently: NaN compares ``False``
    against every bound, so ``NaN > bar`` and ``NaN < bar`` are both
    ``False`` and a bare comparison never fires. This helper tests
    finiteness first, so a non-finite value always fails regardless of
    ``bad``.

    Parameters
    ----------
    name : str
        Short name of the gated quantity, used only in the non-finite
        message.
    value : float
        The measured value being gated.
    bad : bool
        Whether the value fails its bound; only consulted when ``value`` is
        finite.
    message : str
        Failure message to use when ``value`` is finite and ``bad``.

    Returns
    -------
    str or None
        A failure message, or ``None`` if the gate passes.
    """
    if not np.isfinite(value):
        return f"{name} is non-finite ({value!r})"
    if bad:
        return message
    return None


def cumulative_integral(times: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Trapezoid-integrate a per-particle time series from the first snapshot.

    Parameters
    ----------
    times : numpy.ndarray
        Snapshot times, s.
    values : numpy.ndarray
        Values of shape ``(n_snapshots, n_particles)``.

    Returns
    -------
    numpy.ndarray
        Running integral of the same shape, zero at the first snapshot, so the
        integral over snapshots ``[i, j]`` is ``result[j] - result[i]``.
    """
    segments = 0.5 * (values[1:] + values[:-1]) * np.diff(times)[:, np.newaxis]
    return np.vstack([np.zeros((1, values.shape[1])), np.cumsum(segments, axis=0)])


def report_shells(
    history: Dict[str, np.ndarray], index: int, n_shells: int, shielded: bool
) -> None:
    """Print the radial profiles at one snapshot, including the mode diagnostic.

    Parameters
    ----------
    history : dict of str to numpy.ndarray
        Output of :func:`build_history`.
    index : int
        Snapshot index to report.
    n_shells : int
        Number of equal-width radial shells.
    shielded : bool
        Whether a shielding mode was active. Without one the column plays no
        role, so the column ratio is not printed.
    """
    radius = history["radius"]
    edges = np.linspace(0.0, radius.max(), n_shells + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    # Eq. (4) divided by the shielding length gives back 2 n_H2.
    n_H2 = 0.5 * history["column"][index] / history["length"][index]

    shell_n_H2 = np.full(n_shells, np.nan)
    for i in range(n_shells):
        mask = (radius >= edges[i]) & (radius < edges[i + 1])
        if np.any(mask):
            shell_n_H2[i] = np.median(n_H2[mask])
    # Eq. (10), from the star to each shell centre.
    geometric = np.cumsum(shell_n_H2 * np.diff(edges)) - 0.5 * shell_n_H2 * np.diff(
        edges
    )

    print(
        f"{'r [pc]':>9} {'G_0':>11} {'k_diss,0 [1/s]':>15} "
        f"{'l_shield [pc]':>14} {'N_H2 [1/cm2]':>13} {'f_shield':>10} "
        f"{'N_H2/N_geom':>12}"
    )
    for i in range(n_shells):
        mask = (radius >= edges[i]) & (radius < edges[i + 1])
        if not np.any(mask):
            continue
        column = np.median(history["column"][index][mask])
        ratio = column / geometric[i] if shielded and geometric[i] > 0.0 else np.nan
        print(
            f"{centres[i] / PARSEC_CGS:9.3f} "
            f"{np.median(history['habing'][index][mask]):11.4g} "
            f"{np.median(history['rate'][index][mask]):15.4g} "
            f"{np.median(history['length'][index][mask]) / PARSEC_CGS:14.4g} "
            f"{column:13.4g} "
            f"{np.median(history['f_shield'][index][mask]):10.4g} "
            f"{ratio:12.4g}"
        )


def make_figure(
    history: Dict[str, np.ndarray],
    predicted: np.ndarray,
    measured: np.ndarray,
    config: str,
    output: str,
) -> None:
    """Plot the per-particle comparison and the median rate history.

    Parameters
    ----------
    history : dict of str to numpy.ndarray
        Output of :func:`build_history`.
    predicted : numpy.ndarray
        Predicted quantity per qualifying particle: e-folds (thin) or
        shielding ratio (thick).
    measured : numpy.ndarray
        Measured counterpart of ``predicted``.
    config : str
        ``thin`` or ``thick``.
    output : str
        Figure filename.
    """
    times = history["time"]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    label = "e-folds" if config == "thin" else "shielding ratio"
    axes[0].loglog(predicted, measured, ".", markersize=2, alpha=0.4)
    bounds = [
        min(predicted.min(), measured.min()),
        max(predicted.max(), measured.max()),
    ]
    axes[0].loglog(bounds, bounds, "k--", linewidth=1)
    axes[0].set_xlabel(f"predicted {label}")
    axes[0].set_ylabel(f"measured {label}")

    axes[1].semilogy(
        times, np.median(history["rate"], axis=1), "o-", label=r"$k_{\rm diss,0}$"
    )
    axes[1].semilogy(
        times,
        np.median(history["rate"] * history["f_shield"], axis=1),
        "s-",
        label=r"$k_{\rm diss,0}\,f_{\rm shield}$",
    )
    axes[1].set_xlabel("time [s]")
    axes[1].set_ylabel(r"median rate [s$^{-1}$]")
    axes[1].legend()

    figure.tight_layout()
    figure.savefig(output, dpi=140)
    print(f"Wrote {output}")


def main() -> int:
    """Run the check.

    Returns
    -------
    int
        0 on pass, 1 on failure.
    """
    options = parse_options()
    filenames = sorted(glob.glob(options.snapshot))
    if len(filenames) < 4:
        print(f"Need at least 4 snapshots, found {len(filenames)}")
        return 1

    path_in_kernel_radii = (
        1.0 if options.h2_self_shielding_path == "kernel_radius" else 2.0
    )
    history = build_history(
        filenames,
        options.h2_self_shielding,
        options.kernel_gamma,
        path_in_kernel_radii,
    )
    gated_config = "thin" if options.config == "draine_spectrum" else options.config
    times = history["time"]
    n_particles = history["x_H2"].shape[1]
    columns = np.arange(n_particles)

    floor_margin = options.floor_margin if gated_config == "thin" else 0.0
    start, end, reached = particle_windows(
        history, options.field_arrival_fraction, floor_margin
    )

    # Grackle receives the rate clamped at zero.
    rate = np.maximum(history["rate"], 0.0)
    shielded = cumulative_integral(times, rate * history["f_shield"])
    unshielded = cumulative_integral(times, rate)
    shielded_integral = shielded[end, columns] - shielded[start, columns]
    unshielded_integral = unshielded[end, columns] - unshielded[start, columns]
    measured_drop = -np.log(
        history["x_H2"][end, columns] / history["x_H2"][start, columns]
    )

    qualifies = reached & (end - start >= 2) & (unshielded_integral >= 0.5)
    n_qualifying = int(qualifies.sum())
    window_f_shield = np.array(
        [
            np.median(history["f_shield"][start[i] : end[i] + 1, i])
            for i in np.flatnonzero(qualifies)
        ]
    )
    print(
        f"Snapshots: {len(filenames)} ({times[0]:.4g} to {times[-1]:.4g} s); "
        f"qualifying particles: {n_qualifying} of {n_particles}"
    )
    if n_qualifying < 0.1 * n_particles:
        print("FAIL: fewer than 10 percent of the particles have a usable window")
        return 1

    median_f_shield = float(np.median(window_f_shield))
    median_window = float(np.median(end[qualifies] - start[qualifies]))
    print(f"Median window length: {median_window:.0f} snapshots")
    print(f"Median f_shield over the windows: {median_f_shield:.4g}")
    print(
        f"Median integral k_diss,0 dt: "
        f"{np.median(unshielded_integral[qualifies]):.4g} (unshielded e-folds)"
    )
    print(
        f"Median integral k_diss,0 f_shield dt: "
        f"{np.median(shielded_integral[qualifies]):.4g} (predicted e-folds)"
    )
    print(f"Median measured e-folds: {np.median(measured_drop[qualifies]):.4g}")
    lit = history["habing"][-1] > 0.0
    measured_rate_ratio = np.median(
        rate[-1][lit] / (DB96_UNSHIELDED_RATE_CGS * history["habing"][-1][lit])
    )
    lw_fraction = np.median(history["lw_fraction"][-1][lit])
    print(
        f"Cross-check, median G_0 at the last snapshot: "
        f"{np.median(history['habing'][-1]):.4g}; Draine and Bertoldi (1996) "
        f"3.3e-11 G_0 = {DB96_UNSHIELDED_RATE_CGS * np.median(history['habing'][-1]):.4g} "
        f"1/s vs Eq. (1) {np.median(rate[-1]):.4g} 1/s"
    )
    print(
        f"Rate normalisation: Eq. (1) / (3.3e-11 G_0) median "
        f"{measured_rate_ratio:.4g}; LW fraction u_LW/(u_PE+u_LW) "
        f"{lw_fraction:.4g}, times {RATE_RATIO_PER_LW_FRACTION:.4g} = "
        f"{RATE_RATIO_PER_LW_FRACTION * lw_fraction:.4g} (Draine 1978 field: "
        f"LW fraction 0.148, ratio 0.92)"
    )
    print()
    print(f"Radial profiles at the last snapshot:")
    report_shells(
        history,
        len(filenames) - 1,
        options.n_shells,
        options.h2_self_shielding != 0,
    )
    print()

    failures: List[str] = []
    if gated_config == "thin":
        # First snapshot, inside each window, at which the predicted
        # integral reaches the reference number of e-folds.
        snapshot_index = np.arange(len(times))[:, np.newaxis]
        inside = (snapshot_index >= start[np.newaxis, :]) & (
            snapshot_index <= end[np.newaxis, :]
        )
        reached_target = inside & (
            shielded - shielded[start, columns][np.newaxis, :] >= options.thin_efolds
        )
        compared = qualifies & reached_target.any(axis=0)
        at = np.argmax(reached_target, axis=0)
        predicted_at = shielded[at, columns] - shielded[start, columns]
        measured_at = -np.log(
            history["x_H2"][at, columns] / history["x_H2"][start, columns]
        )
        relative_error = np.abs(np.exp(measured_at - predicted_at) - 1.0)[compared]
        median_error = float(np.median(relative_error))
        efold_ratio = measured_drop[qualifies] / shielded_integral[qualifies]
        print(
            f"THIN: {int(compared.sum())} particles reach "
            f"{options.thin_efolds:g} predicted e-folds; relative error on "
            f"x_H2 there, median {median_error:.4g}, 16th to 84th percentile "
            f"{np.percentile(relative_error, 16):.4g} to "
            f"{np.percentile(relative_error, 84):.4g} (bar {options.thin_tol})"
        )
        print(
            f"THIN: measured over predicted e-folds, full windows, median "
            f"{np.median(efold_ratio):.4g}"
        )
        if int(compared.sum()) < 0.1 * n_particles:
            failures.append(
                f"fewer than 10 percent of the particles reach "
                f"{options.thin_efolds:g} predicted e-folds"
            )
        for message in (
            check_gate(
                "THIN relative error",
                median_error,
                median_error > options.thin_tol,
                f"relative error {median_error:.4g} above {options.thin_tol}",
            ),
            check_gate(
                "THIN f_shield",
                median_f_shield,
                median_f_shield < options.thin_f_min,
                f"f_shield {median_f_shield:.4g} below {options.thin_f_min}: "
                "this run is not in the unshielded limit",
            ),
        ):
            if message is not None:
                failures.append(message)
        predicted_plot = shielded_integral[qualifies]
        measured_plot = measured_drop[qualifies]
    else:
        measured_ratio = measured_drop[qualifies] / unshielded_integral[qualifies]
        predicted_ratio = shielded_integral[qualifies] / unshielded_integral[qualifies]
        quotient = measured_ratio / predicted_ratio
        median_quotient = float(np.median(quotient))
        relative_error = abs(median_quotient - 1.0)
        print(
            f"THICK: shielding ratio measured {np.median(measured_ratio):.4g} vs "
            f"predicted {np.median(predicted_ratio):.4g}; per-particle "
            f"quotient median {median_quotient:.4g}, 16th to 84th percentile "
            f"{np.percentile(quotient, 16):.4g} to "
            f"{np.percentile(quotient, 84):.4g} (bar {options.thick_tol})"
        )
        for message in (
            check_gate(
                "THICK shielding ratio",
                relative_error,
                relative_error > options.thick_tol,
                f"shielding ratio off by {relative_error:.4g}, above "
                f"{options.thick_tol}",
            ),
            check_gate(
                "THICK f_shield",
                median_f_shield,
                median_f_shield > options.thick_f_max,
                f"f_shield {median_f_shield:.4g} above {options.thick_f_max}: "
                "this run is not strongly shielded",
            ),
        ):
            if message is not None:
                failures.append(message)
        predicted_plot = predicted_ratio
        measured_plot = np.maximum(measured_ratio, 1e-12)

    make_figure(history, predicted_plot, measured_plot, gated_config, options.output)

    print()
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
