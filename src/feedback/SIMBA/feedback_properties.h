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
#ifndef SWIFT_SIMBA_FEEDBACK_PROPERTIES_H
#define SWIFT_SIMBA_FEEDBACK_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"
#include "inline.h"

/**
 * @brief Properties of the EAGLE feedback model.
 */
struct feedback_props {
  /* Parameters for calculating rotational velocity of galaxy */
  float simba_host_galaxy_mass_norm;
  float simba_v_circ_exp;

  /* Parameters for calculating ejection velocity */
  float galsf_firevel;
  float galsf_firevel_slope;
  float scale_factor_norm;
  float vwvf_scatter; // ALEXEI: rename this variable to something intelligible.

  /* Delay time */
  float simba_delay_time; //ALEXEI use this variable to read in value from yml file, to be copied later to individual particles so that feedback iact function can access it. Perhaps move somewhere else
};

/**
 * @brief Initialize the global properties of the feedback scheme.
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

  /* Initialize parameters for calculating rotational velocity of galaxy */
  fp->simba_host_galaxy_mass_norm = parser_get_param_float(params, "SIMBAFeedback:galaxy_mass_norm"); // 102.329 ALEXEI: guide values added in until figured out what are appropriate values.
  fp->simba_v_circ_exp = parser_get_param_float(params, "SIMBAFeedback:v_circ_exp"); // 0.26178;

  /* Initialize parameters for calculating ejection velocity */
  fp->galsf_firevel = parser_get_param_float(params,"SIMBAFeedback:galsf_vel"); // Look in config file 1.6
  fp->galsf_firevel_slope = parser_get_param_float(params, "SIMBAFeedback:galsf_vel_slope"); // 0.12;
  fp->scale_factor_norm = parser_get_param_float(params,"SIMBAFeedback:scale_factor_norm"); // 200.;
  fp->vwvf_scatter = parser_get_param_float(params,"SIMBAFeedback:wind_scatter"); // 0.1;

  /* read in delay time */
  fp->simba_delay_time = parser_get_param_float(params,"SIMBAFeedback:delay_time");

}

#endif /* SWIFT_SIMBA_FEEDBACK_PROPERTIES_H */
