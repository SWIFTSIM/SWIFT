/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_IACT_H
#define SWIFT_EAGLE_FEEDBACK_IACT_H

/* Local includes */
#include "random.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xp Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float *dx,
                                    const float hi, const float hj,
                                    struct spart *restrict si,
                                    const struct part *restrict pj,
                                    const struct xpart *restrict xp,
                                    const struct cosmology *restrict cosmo,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.to_collect.ngb_mass += mj;

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  const float rho = hydro_get_comoving_density(pj);
  if (rho != 0.f)
    si->feedback_data.to_collect.enrichment_weight_inv += wi / rho;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xp Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  const struct spart *restrict si,
                                  struct part *restrict pj,
                                  struct xpart *restrict xp,
                                  const struct cosmology *restrict cosmo,
                                  const integertime_t ti_current) {

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);

  /* Compute weighting for distributing feedback quantities */
  float Omega_frac;
  if (rho_j != 0.f) {
    Omega_frac = si->feedback_data.to_distribute.enrichment_weight * wi / rho_j;
  } else {
    Omega_frac = 0.f;
  }

  /* Update particle mass */
  const double current_mass = hydro_get_mass(pj);
  const double delta_mass = si->feedback_data.to_distribute.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Update total metallicity */
  const double current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const double delta_metal_mass_total =
      si->feedback_data.to_distribute.total_metal_mass * Omega_frac;
  const double new_metal_mass_total =
      current_metal_mass_total + delta_metal_mass_total;

  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total / new_mass;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.to_distribute.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] = new_metal_mass / new_mass;
  }

  /* Update iron mass fraction from SNIa  */
  const double current_iron_from_SNIa_mass =
      pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  const double delta_iron_from_SNIa_mass =
      si->feedback_data.to_distribute.Fe_mass_from_SNIa * Omega_frac;
  const double new_iron_from_SNIa_mass =
      current_iron_from_SNIa_mass + delta_iron_from_SNIa_mass;

  pj->chemistry_data.iron_mass_fraction_from_SNIa =
      new_iron_from_SNIa_mass / new_mass;

  /* Update mass fraction from SNIa  */
  const double current_mass_from_SNIa =
      pj->chemistry_data.mass_from_SNIa * current_mass;
  const double delta_mass_from_SNIa =
      si->feedback_data.to_distribute.mass_from_SNIa * Omega_frac;
  const double new_mass_from_SNIa =
      current_mass_from_SNIa + delta_mass_from_SNIa;

  pj->chemistry_data.mass_from_SNIa = new_mass_from_SNIa / new_mass;

  /* Update metal mass fraction from SNIa */
  const double current_metal_mass_from_SNIa =
      pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const double delta_metal_mass_from_SNIa =
      si->feedback_data.to_distribute.metal_mass_from_SNIa * Omega_frac;
  const double new_metal_mass_from_SNIa =
      current_metal_mass_from_SNIa + delta_metal_mass_from_SNIa;

  pj->chemistry_data.metal_mass_fraction_from_SNIa =
      new_metal_mass_from_SNIa / new_mass;

  /* Update mass fraction from SNII  */
  const double current_mass_from_SNII =
      pj->chemistry_data.mass_from_SNII * current_mass;
  const double delta_mass_from_SNII =
      si->feedback_data.to_distribute.mass_from_SNII * Omega_frac;
  const double new_mass_from_SNII =
      current_mass_from_SNII + delta_mass_from_SNII;

  pj->chemistry_data.mass_from_SNII = new_mass_from_SNII / new_mass;

  /* Update metal mass fraction from SNII */
  const double current_metal_mass_from_SNII =
      pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const double delta_metal_mass_from_SNII =
      si->feedback_data.to_distribute.metal_mass_from_SNII * Omega_frac;
  const double new_metal_mass_from_SNII =
      current_metal_mass_from_SNII + delta_metal_mass_from_SNII;

  pj->chemistry_data.metal_mass_fraction_from_SNII =
      new_metal_mass_from_SNII / new_mass;

  /* Update mass fraction from AGB  */
  const double current_mass_from_AGB =
      pj->chemistry_data.mass_from_AGB * current_mass;
  const double delta_mass_from_AGB =
      si->feedback_data.to_distribute.mass_from_AGB * Omega_frac;
  const double new_mass_from_AGB = current_mass_from_AGB + delta_mass_from_AGB;

  pj->chemistry_data.mass_from_AGB = new_mass_from_AGB / new_mass;

  /* Update metal mass fraction from AGB */
  const double current_metal_mass_from_AGB =
      pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  const double delta_metal_mass_from_AGB =
      si->feedback_data.to_distribute.metal_mass_from_AGB * Omega_frac;
  const double new_metal_mass_from_AGB =
      current_metal_mass_from_AGB + delta_metal_mass_from_AGB;

  pj->chemistry_data.metal_mass_fraction_from_AGB =
      new_metal_mass_from_AGB / new_mass;

  /* /\* Update momentum *\/ */
  /* for (int i = 0; i < 3; i++) { */
  /*   pj->v[i] += si->feedback_data.to_distribute.mass * Omega_frac * */
  /*               (si->v[i] - pj->v[i]); */
  /* } */

  /* Energy feedback */
  // d_energy += si->feedback_data.to_distribute.d_energy *
  // Omega_frac;

  /* if (feedback_props->continuous_heating) { */
  /*   // We're doing ONLY continuous heating */
  /*   d_energy += si->feedback_data.to_distribute.num_SNIa * */
  /*               feedback_props->total_energy_SNe * Omega_frac * */
  /*               si->mass_init; */
  /* } else { */

  /* /\* Add contribution from thermal and kinetic energy of ejected material */
  /*    (and continuous SNIa feedback) *\/ */
  /* u_init = hydro_get_physical_internal_energy(pj, xp, cosmo); */
  /* du = d_energy / hydro_get_mass(pj); */
  /* hydro_set_physical_internal_energy(pj, xp, cosmo, u_init + du); */
  /* hydro_set_drifted_physical_internal_energy(pj, cosmo, u_init + du); */

  /* Get the SNII feedback properties */
  const float prob = si->feedback_data.to_distribute.SNII_heating_probability;
  const float delta_u = si->feedback_data.to_distribute.SNII_delta_u;

  /* Are we doing some SNII feedback? */
  if (prob > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(si->id + pj->id, ti_current,
                                            random_number_stellar_feedback);
    /* Are we lucky? */
    if (rand < prob) {

      /* Compute new energy of this particle */
      const float u_init = hydro_get_physical_internal_energy(pj, xp, cosmo);
      const float u_new = u_init + delta_u;

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xp, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      message(
          "We did some heating! id %llu star id %llu probability %.5e "
          "random_num %.5e du %.5e du/ini %.5e",
          pj->id, si->id, prob, rand, delta_u, delta_u / u_init);
    }
  }
}

#endif /* SWIFT_EAGLE_FEEDBACK_IACT_H */
