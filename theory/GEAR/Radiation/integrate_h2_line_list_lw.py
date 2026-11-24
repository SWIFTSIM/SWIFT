"""Integrate the H2 line list over the Draine field to get sigma_H2/E_LW.

Standalone, runs in well under a minute and prints the first value of the
three-computation derivation in 04_isrf_packet_transport.tex. The coefficient is

    kappa = sum_l sigma_int,l N_lambda(lambda_l) / F_LW,

with sigma_int,l the integrated dissociation cross section of band l
(cm^2 A), N_lambda the photon flux per unit wavelength of the field
(photons s^-1 cm^-2 A^-1) and F_LW the energy flux of the field over the LW
band (11.2 to 13.6 eV, the band edges of radiation.h). Only bands inside the
LW band enter, since the code carries only the LW energy.

The field is the Draine (1978) fit of Sternberg et al. (2014), their Eq. 2.
The band list is the Leiden photodissociation database file h2.pd
(Heays et al. 2017), which holds the H2 data of Abgrall et al. (2000): 32
Lyman and Werner vibronic bands summed over the rotational levels of a gas at
300 K.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

H_PLANCK_CGS = 6.62607015e-27  # erg s
C_LIGHT_CGS = 2.99792458e10  # cm s^-1
ELECTRON_VOLT_CGS = 1.602176634e-12  # erg
LW_LOWER_EDGE_EV = 11.2
LW_UPPER_EDGE_EV = 13.6

# (wavelength in Angstrom, integrated cross section in cm^2 Angstrom), h2.pd.
H2_BANDS = (
    (912.076, 1.73028e-17),
    (915.991, 5.19425e-19),
    (916.526, 7.47491e-22),
    (918.501, 2.71922e-17),
    (925.283, 2.79348e-17),
    (930.344, 7.55216e-19),
    (931.326, 2.39592e-22),
    (932.265, 4.1567e-17),
    (939.792, 4.12305e-17),
    (946.57, 2.17386e-17),
    (947.815, 5.97175e-17),
    (948.538, 4.13045e-24),
    (955.761, 5.79599e-17),
    (964.677, 5.28225e-17),
    (965.781, 5.81254e-18),
    (966.824, 1.49778e-25),
    (973.381, 6.93884e-17),
    (982.845, 7.08941e-17),
    (986.649, 6.56527e-19),
    (987.586, 7.36158e-27),
    (992.768, 9.35668e-17),
    (1003.29, 7.79871e-17),
    (1010.46, 3.94255e-27),
    (1010.48, 1.19154e-19),
    (1014.21, 5.38714e-17),
    (1026.09, 1.06602e-17),
    (1038.07, 5.01481e-18),
    (1050.98, 2.92367e-19),
    (1064.67, 3.48684e-20),
    (1078.68, 1.59918e-21),
    (1094.39, 7.21152e-24),
    (1109.84, 8.9413e-26),
)


def edge_wavelength_angstrom(energy_ev: float) -> float:
    """Convert a photon energy to a wavelength.

    Parameters
    ----------
    energy_ev : float
        Photon energy, eV.

    Returns
    -------
    float
        Wavelength, Angstrom.
    """
    return H_PLANCK_CGS * C_LIGHT_CGS / (energy_ev * ELECTRON_VOLT_CGS) * 1.0e8


def draine_photon_flux_per_angstrom(wavelength_angstrom: float) -> float:
    """Evaluate the Draine field of Sternberg et al. (2014), Eq. 2.

    Parameters
    ----------
    wavelength_angstrom : float
        Wavelength, Angstrom.

    Returns
    -------
    float
        Photon flux per unit wavelength, photons s^-1 cm^-2 A^-1.
    """
    lam = wavelength_angstrom
    # Photons s^-1 cm^-2 Hz^-1 sr^-1, then 4 pi sr and |dnu/dlambda|.
    intensity = (1.068e-3 / lam - 1.719 / lam**2 + 6.853e2 / lam**3) / (4.0 * np.pi)
    dnu_dlambda = C_LIGHT_CGS / (lam * 1.0e-8) / lam
    return 4.0 * np.pi * intensity * dnu_dlambda


def lw_energy_flux() -> float:
    """Integrate the energy flux of the Draine field over the LW band.

    Returns
    -------
    float
        Energy flux, erg s^-1 cm^-2.
    """
    lam_low = edge_wavelength_angstrom(LW_UPPER_EDGE_EV)
    lam_high = edge_wavelength_angstrom(LW_LOWER_EDGE_EV)

    def energy_per_angstrom(lam: float) -> float:
        photon_energy = H_PLANCK_CGS * C_LIGHT_CGS / (lam * 1.0e-8)
        return draine_photon_flux_per_angstrom(lam) * photon_energy

    return quad(energy_per_angstrom, lam_low, lam_high, limit=200)[0]


def line_list_dissociation_rate() -> float:
    """Sum the in-band bands of the line list over the Draine field.

    Returns
    -------
    float
        Free-space dissociation rate, s^-1.
    """
    lam_low = edge_wavelength_angstrom(LW_UPPER_EDGE_EV)
    lam_high = edge_wavelength_angstrom(LW_LOWER_EDGE_EV)
    return sum(
        sigma * draine_photon_flux_per_angstrom(lam)
        for lam, sigma in H2_BANDS
        if lam_low <= lam <= lam_high
    )


def main() -> None:
    """Print the line-list rate, the LW energy flux and their ratio."""
    rate = line_list_dissociation_rate()
    flux = lw_energy_flux()
    print(f"k_diss (line list, Draine field) = {rate:.4e} s^-1")
    print(f"F_LW (Draine field, 11.2 to 13.6 eV) = {flux:.4e} erg s^-1 cm^-2")
    print(f"kappa = k_diss / F_LW = {rate / flux:.4e} cm^2 erg^-1")


if __name__ == "__main__":
    main()
