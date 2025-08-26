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
#ifndef SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H
#define SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H

/* Standard headers */
#include <float.h>

/* Local headers */
#include "gravity_properties.h"
#include "hydro_properties.h"
#include "inline.h"
#include "kernel_hydro.h"

#define df_props_default_eta_neighbours 1.2348        // THIS NEEDS TO CHANGE "48 Ngb" with the cubic spline kernel
#define df_props_default_max_iterations 30
#define df_props_default_volume_change 1.4f
#define df_props_default_h_max FLT_MAX
#define df_props_default_h_min_ratio 0.f
#define df_props_default_h_tolerance 1e-4

/**
 * @brief Properties of the DF model.
 */
struct df_props {

  /* ------ Smoothing lengths parameters ---------- */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal smoothing length (internal units) */
  float h_max;

  /*! Minimal smoothing length expressed as ratio to softening length */
  float h_min_ratio;

  /*! Minimal smoothing length (internal units) */
  float h_min;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

};

/**
 * @brief Initialize the global properties of the df scheme.
 *
 * Nothing to do here for the no feedback model.
 *
 * @param dfp The #df_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void df_props_init(struct df_props *dfp,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct hydro_props *hydro_props,
                                       const struct cosmology *cosmo) {

  /* ------ Smoothing lengths parameters ---------- */

  /* Kernel properties */
  dfp->eta_neighbours = parser_get_opt_param_float(params, "DynamicalFriction:resolution_eta", df_props_default_eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  dfp->h_tolerance = parser_get_opt_param_float(params, "DynamicalFriction:h_tolerance",
                                              df_props_default_h_tolerance);

  /* Get derived properties */
  dfp->target_neighbours = pow_dimension(dfp->eta_neighbours) * kernel_norm;
  const float delta_eta = dfp->eta_neighbours * (1.f + dfp->h_tolerance);
  dfp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(dfp->eta_neighbours)) *
      kernel_norm;

  /* Maximal smoothing length */
  dfp->h_max = parser_get_opt_param_float(params, "DynamicalFriction:h_max",
                                        df_props_default_h_max);

  /* Minimal smoothing length ratio to softening */
  dfp->h_min_ratio = parser_get_opt_param_float(params, "DynamicalFriction:h_min_ratio",
                                              df_props_default_h_min_ratio);

  /* Temporarily set the minimal softening to 0. */
  dfp->h_min = 0.f;

  /* Number of iterations to converge h */
  dfp->max_smoothing_iterations = parser_get_opt_param_int(
      params, "DynamicalFriction:max_ghost_iterations", df_props_default_max_iterations);

  if (dfp->max_smoothing_iterations <= 10)
    error("The number of smoothing length iterations should be > 10");

}

/**
 * @brief Update the global properties of the hydro scheme for that time-step.
 *
 * @param p The properties to update.
 * @param gp The properties of the gravity scheme.
 * @param cosmo The cosmological model.
 */
void df_props_update(struct df_props *dfp, const struct gravity_props *gp,
                        const struct cosmology *cosmo) {

  /* Update the minimal allowed smoothing length
   *
   * We follow Gadget here and demand that the kernel support (h * gamma)
   * is a fixed fraction of the radius at which the softened forces
   * recover a Newtonian behaviour (i.e. 2.8 * Plummer equivalent softening
   * in the case of a cubic spline kernel). */
  dfp->h_min = dfp->h_min_ratio * gp->epsilon_baryon_cur / kernel_gamma;
}

/**
 * @brief Write a #stars_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
INLINE static void dynamical_friction_struct_dump(const struct df_props *p,
                                           FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct df_props), 1, stream,
                       "dfprops", "dynamical friction props");
}

/**
 * @brief Restore a df_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void dynamical_friction_struct_restore(
    const struct df_props *props, FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct df_props), 1,
                      stream, NULL, "dynamical friction props");
}

#endif /* SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H */
