/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_SINK_SETTERS_H
#define SWIFT_GEAR_SINK_SETTERS_H

#include "sink_part.h"

/**
 * @file src/sink/GEAR/sink_setters.h
 * @brief Setters functions for GEAR sink scheme to avoid exposing
 * implementation details to the outer world. Keep the code clean and lean.
 */

/**
 * @brief Set the birth time/scale-factor of a #sink particle.
 *
 * @param sink The #sink.
 * @param birth_time Birth time of the #sink.
 * @param birth_scale_factor Birth scale-factor of the star.
 * @param with_cosmology If we run with cosmology.
 */

__attribute__((always_inline)) INLINE void
sink_set_sink_birth_time_or_scale_factor(struct sink *restrict sink,
                                         const float birth_time,
                                         const float birth_scale_factor,
                                         const int with_cosmology) {
  if (with_cosmology) {
    sink->birth_scale_factor = birth_scale_factor;
  } else {
    sink->birth_time = birth_time;
  }
}

/**
 * @brief Update the target mass of the sink particle.
 *
 * @param sink the #sink particle.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param e The #engine
 * @param star_counter The star loop counter.
 */
INLINE static void sink_update_target_mass(struct sink *sink,
                                           const struct sink_props *sink_props,
                                           const struct engine *e,
                                           int star_counter) {

  float random_number = random_unit_interval_part_ID_and_index(
      sink->id, star_counter, e->ti_current, random_number_sink_formation);

  const struct feedback_props *feedback_props = e->feedback_props;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_sink_total_iron_mass_fraction_for_feedback(sink);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then the sink generates first star particles. */
  const int is_first_star = metal < threshold;
  const struct stellar_model *model;
  double minimal_discrete_mass_Msun;

  /* Take the correct values if your are a first star or not. */
  if (!is_first_star) /* (metal >= threshold)*/ {
    model = &feedback_props->stellar_model;
    minimal_discrete_mass_Msun = sink_props->minimal_discrete_mass_Msun;
  } else {
    model = &feedback_props->stellar_model_first_stars;
    minimal_discrete_mass_Msun =
        sink_props->minimal_discrete_mass_first_stars_Msun;
  }

  const struct initial_mass_function *imf = &model->imf;

  if (random_number < imf->sink_Pc) {
    /* We are dealing with the continous part of the IMF. */
    sink->target_mass_Msun = imf->stellar_particle_mass_Msun;
    sink->target_type = star_population_continuous_IMF;
  } else {
    /* We are dealing with the discrete part of the IMF. */
    random_number = random_unit_interval_part_ID_and_index(
        sink->id, star_counter + 1, e->ti_current,
        random_number_sink_formation);
    double m = initial_mass_function_sample_power_law(
        minimal_discrete_mass_Msun, imf->mass_max, imf->exp[imf->n_parts - 1],
        random_number);
    sink->target_mass_Msun = m;
    sink->target_type = single_star;
  }

  /* Also store the mass of the IMF into the sink. */
  double dummy;
  ;
  initial_mass_function_compute_Mc_Md_Mtot(imf, &dummy, &dummy,
                                           &sink->mass_IMF);

  /* Convert from M_sun to internal units. */
  sink->mass_IMF *= e->physical_constants->const_solar_mass;
}
#endif /* SWIFT_GEAR_SINK_SETTERS_H */
