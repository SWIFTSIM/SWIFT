/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.reavz@epfl.ch)
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
#ifndef SWIFT_AGORA_FEEDBACK_PROPERTIES_H
#define SWIFT_AGORA_FEEDBACK_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"
#include "inline.h"

/**
 * @brief Properties of the EAGLE feedback model.
 */
struct feedback_props {

  /*! Supernovae energy effectively deposited */
  float supernovae_efficiency;

  /*! Energy released per CCSN in code unit */
  double energy_per_CCSN;

  /*! Supernovae explosion time */
  double supernovae_explosion_time;

  /*! CCSNe per solar mass  */
  double ccsne_per_solar_mass;

  /*! Ejected mass per CCSN */
  double ejected_mass_per_CCSN;

  /*! Ejected Fe mass per CCSN */
  double ejected_Fe_mass_per_CCSN;

  /*! Ejected metal mass per CCSN */
  double ejected_metal_mass_per_CCSN;
};

/**
 * @brief Print the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void feedback_props_print(
    const struct feedback_props *feedback_props) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* Print the feedback properties */
  message("Supernovae efficiency = %.2g",
          feedback_props->supernovae_efficiency);

  message("Energy relased per CCSN (code units) = %.2g",
          feedback_props->energy_per_CCSN);

  message("Supernovae explosition time (code units) = %.2g",
          feedback_props->supernovae_explosion_time);

  message("CCSNe per solar mass = %.2g", feedback_props->ccsne_per_solar_mass);

  message("Ejected mass per CCSN (code units) = %.2g",
          feedback_props->ejected_mass_per_CCSN);

  message("Ejected Fe mass per CCSN (code units) = %.2g",
          feedback_props->ejected_Fe_mass_per_CCSN);

  message("Ejected metal mass per CCSN (code units) = %.2g",
          feedback_props->ejected_metal_mass_per_CCSN);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * Nothing to do here for the no feedback model.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void feedback_props_init(struct feedback_props *fp,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct hydro_props *hydro_props,
                                       const struct cosmology *cosmo) {

  /* Supernovae energy efficiency */
  double e_efficiency =
      parser_get_param_double(params, "AGORAFeedback:supernovae_efficiency");
  fp->supernovae_efficiency = e_efficiency;

  /* Supernovae energy */
  double e_feedback =
      parser_get_param_double(params, "AGORAFeedback:energy_in_erg_per_CCSN");
  e_feedback /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  fp->energy_per_CCSN = e_feedback;

  /* Supernovae explosion time */
  double e_time = parser_get_param_double(
      params, "AGORAFeedback:supernovae_explosion_time_myr");
  fp->supernovae_explosion_time = e_time * phys_const->const_year * 1e6;

  /* SCCSNe per solar mass */
  fp->ccsne_per_solar_mass =
      parser_get_param_double(params, "AGORAFeedback:ccsne_per_solar_mass");

  /* Ejected gas mass per CCSN */
  fp->ejected_mass_per_CCSN = parser_get_param_double(
      params, "AGORAFeedback:ejected_mass_in_solar_mass_per_CCSN");
  fp->ejected_mass_per_CCSN *= phys_const->const_solar_mass;

  /* Ejected Fe mass per CCSN */
  fp->ejected_Fe_mass_per_CCSN = parser_get_param_double(
      params, "AGORAFeedback:ejected_Fe_mass_in_solar_mass_per_CCSN");
  fp->ejected_Fe_mass_per_CCSN *= phys_const->const_solar_mass;

  /* Ejected metal mass per CCSN */
  fp->ejected_metal_mass_per_CCSN = parser_get_param_double(
      params, "AGORAFeedback:ejected_metal_mass_in_solar_mass_per_CCSN");
  fp->ejected_metal_mass_per_CCSN *= phys_const->const_solar_mass;

  /* Print the stellar properties */
  feedback_props_print(fp);
}

#endif /* SWIFT_AGORA_FEEDBACK_PROPERTIES_H */
