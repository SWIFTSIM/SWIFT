/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_DEFAULT_NEUTRINO_H
#define SWIFT_DEFAULT_NEUTRINO_H

/* Riemann function zeta(3) */
#define M_ZETA_3 1.2020569031595942853997

/* Compute the ratio of macro particle mass in internal mass units to
 * the mass of one microscopic neutrino in eV.
 *
 * @param cosmo The #cosmology used for this run.
 * @param internal_units The system of units used internally.
 * @param physical_constants The #phys_const used for this run.
 * @param volume The volume occupied by neutrino particles
 * @param nr_nuparts The number of macro neutrino particles
 */
INLINE static double neutrino_mass_factor(
    const struct cosmology *cosmo, const struct unit_system *internal_units,
    const struct phys_const *physical_constants, double volume,
    double nr_nuparts) {
  /* Some constants */
  const double k_b = physical_constants->const_boltzmann_k;
  const double hbar = physical_constants->const_planck_hbar;
  const double c = physical_constants->const_speed_light_c;
  const double eV = physical_constants->const_electron_volt;
  const double eV_mass = eV / (c * c);  // 1 eV/c^2 in internal mass units
  const double prefactor = (1.5 * M_ZETA_3) / (M_PI * M_PI);
  const double T_nu = cosmo->T_nu_0;

  /* Count the number of neutrino flavours according to multiplicity */
  double flavours = 0.;
  for (int i = 0; i < cosmo->N_nu; i++) {
    flavours += cosmo->deg_nu[i];
  }

  /* Compute the comoving number density per flavour */
  const double kThc = k_b * T_nu / (hbar * c);
  const double n = prefactor * kThc * kThc * kThc;

  /* Compute the conversion factor */
  const double mass_factor = nr_nuparts / (flavours * n * volume);

  /* Convert to eV */
  const double mass_factor_eV = mass_factor / eV_mass;

  return mass_factor_eV;
}

#endif /* SWIFT_DEFAULT_NEUTRINO_H */
