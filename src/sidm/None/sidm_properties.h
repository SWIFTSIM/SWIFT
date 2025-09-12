/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_NONE_SIDM_PROPERTIES_H
#define SWIFT_NONE_SIDM_PROPERTIES_H

/* Local header */
#include "parser.h"

/**
 * @brief Properties of SIDM in the Basic model.
 */
struct sidm_props {

  /* ----- Basic neighbour search properties ------ */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /*! Are we using a fixed cutoff radius? (all smoothing length calculations are
   * disabled if so) */
  char use_fixed_r_cut;
};

/**
 * @brief Initialise the SIDM properties from the parameter file.
 *
 * @param sip The #sidm_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void sidm_props_init(struct sidm_props *sip,
                                   const struct phys_const *phys_const,
                                   const struct unit_system *us,
                                   struct swift_params *params,
                                   const struct hydro_props *hydro_props,
                                   const struct cosmology *cosmo) {

  /* We don't use a fixed cutoff radius in this model. This property must always
   * be present, as we use it to skip smoothing length iterations in
   * runner_ghost if a fixed cutoff is being used. */
  sip->use_fixed_r_cut = 0;

  /* Read in the basic neighbour search properties or default to the hydro
     ones if the user did not provide any different values */

  /* Kernel properties */
  sip->eta_neighbours = parser_get_opt_param_float(
      params, "BasicSIDM:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  sip->h_tolerance = parser_get_opt_param_float(params, "BasicSIDM:h_tolerance",
                                                hydro_props->h_tolerance);

  /* Get derived properties */
  sip->target_neighbours = pow_dimension(sip->eta_neighbours) * kernel_norm;
  const float delta_eta = sip->eta_neighbours * (1.f + sip->h_tolerance);
  sip->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(sip->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  sip->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BasicSIDM:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BasicSIDM:max_volume_change", -1);
  if (max_volume_change == -1)
    sip->log_max_h_change = hydro_props->log_max_h_change;
  else
    sip->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));
}

/**
 * @brief Write an sidm_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the sidm properties struct
 * @param stream the file stream
 */
INLINE static void sidm_struct_dump(const struct sidm_props *props,
                                    FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct sidm_props), 1, stream,
                       "sidm props", "SIDM props");
}

/**
 * @brief Restore a sidm_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the sidm properties struct
 * @param stream the file stream
 */
INLINE static void sidm_struct_restore(const struct sidm_props *props,
                                       FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct sidm_props), 1, stream, NULL,
                      "SIDM props");
}

#endif /* SWIFT_NONE_SIDM_PROPERTIES_H */
