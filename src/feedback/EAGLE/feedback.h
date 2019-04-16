/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_EAGLE_H
#define SWIFT_FEEDBACK_EAGLE_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const float age,
                               const float dt);

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.density_weighted_frac_normalisation_inv = 0.f;
  sp->feedback_data.ngb_mass = 0.f;
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {

  feedback_init_spart(sp);
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.to_distribute.mass = 0.f;

  /* Zero the number of SN */
  sp->feedback_data.to_distribute.num_SNIa = 0.f;
  sp->feedback_data.to_distribute.num_SNII = 0.f;
  sp->feedback_data.to_distribute.num_SNe = 0.f;

  /* Zero the enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++)
    sp->feedback_data.to_distribute.metal_mass[i] = 0.f;
  sp->feedback_data.to_distribute.total_metal_mass = 0.f;
  sp->feedback_data.to_distribute.mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.to_distribute.d_energy = 0.f;

  /* Zero the SNII feedback probability */
  sp->feedback_data.to_distribute.SNII_heating_probability = 0.f;

  /* Zero the SNII feedback energy */
  sp->feedback_data.to_distribute.SNII_delta_u = 0.f;
}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_propss The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
__attribute__((always_inline)) INLINE static void feedback_evolve_spart(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const double star_age_beg_step, const double dt) {

  /* Compute amount of enrichment and feedback that needs to be done in this
   * step */
  compute_stellar_evolution(feedback_props, cosmo, sp, us, star_age_beg_step,
                            dt);

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= sp->feedback_data.to_distribute.mass;
}

#endif /* SWIFT_FEEDBACK_EAGLE_H */
