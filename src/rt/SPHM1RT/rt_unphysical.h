/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c)  2022 Tsang Keung Chan (chantsangkeung@gmail.com)
 *                2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_UNPHYSICAL_SPHM1RT_H
#define SWIFT_RT_UNPHYSICAL_SPHM1RT_H

/**
 * @file src/rt/SPHM1RT/rt_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

#include "rt_properties.h"

/**
 * @brief check for and correct if needed unphysical
 * values for a radiation state.
 *
 * @param energy_density pointer to the radiation energy density
 * @param flux pointer to radiation flux (3 dimensional)
 * @param e_old energy density before change to check. Set = 0 if not available
 * @param cred reduced speed of light in (comoving) code unit
 * @param loc indicating the location
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_state(
    float* energy_density, float flux[3], const float e_old, const float cred, 
    const char* loc) {

  /* Check for negative energies */
  /* Note to self for printouts: Maximal allowable F = E * c.
   * In some cases, e.g. while cooling, we don't modify the fluxes,
   * so you can get an estimate of what the photon energy used to be
   * by dividing the printed out fluxes by the speed of light in
   * code units */
  /* this check will never trigger unless you manually remove the
   * -ffinite-math-only flag for gcc */
  if (isinf(*energy_density) || isnan(*energy_density))
    error("Got inf/nan radiation energy case | %.6e | %.6e %.6e %.6e",
          *energy_density, flux[0], flux[1], flux[2]);

  if (*energy_density <= 0.f) {
    *energy_density = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
    return;
  }

  /* Check for too high fluxes */
  double fluxdouble[3];
  fluxdouble[0] = (double)(flux[0]);
  fluxdouble[1] = (double)(flux[1]);
  fluxdouble[2] = (double)(flux[2]);
  const double flux2double = fluxdouble[0] * fluxdouble[0] +
                             fluxdouble[1] * fluxdouble[1] +
                             fluxdouble[2] * fluxdouble[2];
  if (isinf(flux2double) || isnan(flux2double))
    error("Got inf/nan in flux2 | %.6e| %.6e %.6e %.6e %s", flux2double, flux[0],
          flux[1], flux[2], loc);

  const double flux_normdouble = (flux2double == 0.0) ? 0.0 : sqrt(flux2double);
  const float flux_norm = (float)(flux_normdouble);

  const float flux_norm_inv = (flux_norm == 0.f) ? 0.f : 1.f / flux_norm;
  const float flux_max = cred * *energy_density;
  float flux_diff = flux_norm - flux_max;

  if (flux_norm != 0.f) {
    if (flux_diff > 0.f) {
      const float correct = flux_max * flux_norm_inv;
      flux[0] *= correct;
      flux[1] *= correct;
      flux[2] *= correct;
    }
  }
}

/**
 * @brief check whether gas species abundances have physical
 * values and correct small errors if necessary.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_elem_spec(
    struct part* restrict p, const struct rt_props* rt_props) {

  /*************************/
  /* check mass fraction   */
  /*************************/

  /* check individual mass fraction */
  char name[10];
  double test_mass_fraction;
  for (int elem = 0; elem < rt_chemistry_element_count; elem++) {
    test_mass_fraction = p->rt_data.tchem.metal_mass_fraction[elem];
    if (isinf(test_mass_fraction) || isnan(test_mass_fraction))
      error("Got inf/nan test_mass_fraction %d | with value %.6e ", elem,
            test_mass_fraction);
    if (test_mass_fraction < 0.f) {
      sprintf(name, "%s",
              rt_chemistry_get_element_name((enum rt_chemistry_element)elem));
      error("Error: Got negative mass fraction in %s", name);
    }
  }

  /* check total mass fraction */
  float mass_fraction_tot = 0.f;
  for (int elem = 0; elem < rt_chemistry_element_count; elem++) {
    mass_fraction_tot += p->rt_data.tchem.metal_mass_fraction[elem];
  }
  /* Make sure we sum up to 1. */
  if (fabsf(mass_fraction_tot - 1.f) > 1e-3)
    error("Got total mass fraction = %.6g", mass_fraction_tot);

  /*********************/
  /* check abundaces   */
  /*********************/

  double test_abundance;
  for (int spec = 0; spec < rt_species_count; spec++) {
    test_abundance = p->rt_data.tchem.abundances[spec];
    if (isinf(test_abundance) || isnan(test_abundance))
      error("Got inf/nan test_abundance %d | with value %.6e ", spec,
            test_abundance);
    if (test_abundance < 0.f) {
      if (test_abundance < -1e4) {
        sprintf(name, "%s", rt_get_species_name((enum rt_cooling_species)spec));
        message("WARNING: Got negative abundance in %s", name);
      }
      p->rt_data.tchem.abundances[spec] = 0.f;
    }
  }

  /* normalized total abundances summed over many species */
  float abundance_tot_nor = 0.f;
  /* first check Hydrogen */
  abundance_tot_nor = p->rt_data.tchem.abundances[rt_sp_HI] +
                      p->rt_data.tchem.abundances[rt_sp_HII];
  /* Make sure we sum up to 1. */
  if (fabsf(abundance_tot_nor - 1.f) > 1e-3)
    error("Got total abundances of hydrogen gas = %.6g", abundance_tot_nor);

  /* second check Helium */
  abundance_tot_nor = p->rt_data.tchem.abundances[rt_sp_HeI] +
                      p->rt_data.tchem.abundances[rt_sp_HeII] +
                      p->rt_data.tchem.abundances[rt_sp_HeIII];
  /* Helium is more tricky. We need to check AHe=nHe/nH. */
  /* expected total abundance = MHe/MH / am_He * am_H */
  float abundance_tot_exp;
  abundance_tot_exp =
      p->rt_data.tchem.metal_mass_fraction[rt_chemistry_element_He] /
      p->rt_data.tchem.metal_mass_fraction[rt_chemistry_element_H] *
      rt_props->atomicmass[rt_chemistry_element_H] *
      rt_props->atomicmass_inv[rt_chemistry_element_He];
  /* Make sure we sum up to expected. */
  if (fabsf(abundance_tot_nor - abundance_tot_exp) > 1e-4)
    error("Got total abundances of helium gas = %.6g, expect = %.6g",
          abundance_tot_nor, abundance_tot_exp);

  /* third check electron density */
  abundance_tot_exp = p->rt_data.tchem.abundances[rt_sp_HII] +
                      p->rt_data.tchem.abundances[rt_sp_HeII] +
                      2.0f * p->rt_data.tchem.abundances[rt_sp_HeIII];
  /* Make sure we sum up to the expected. */
  if (fabsf(p->rt_data.tchem.abundances[rt_sp_elec] - abundance_tot_exp) >
      1e-3 * abundance_tot_exp)
    error("Got total abundances of electron = %.6g; expected = %.6g",
          p->rt_data.tchem.abundances[rt_sp_elec], abundance_tot_exp);
}

#endif /* SWIFT_RT_UNPHYSICAL_SPHM1RT_H */
