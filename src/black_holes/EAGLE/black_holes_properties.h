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
#ifndef SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H
#define SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"

/**
 * @brief Properties of black holes and AGN feedback in the EAGEL model.
 */
struct black_holes_props {

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

  /* ----- Initialisation properties  ------ */

  /* ----- Properties of the accretion model ------ */

  /*! Maximal fraction of the Eddington rate allowed. */
  float f_Edd;

  /*! Radiative efficiency of the black holes. */
  float epsilon_r;

  /*! Feedback coupling efficiency of the black holes. */
  float epsilon_f;

  /* ---- Properties of the feedback model ------- */

  /*! Temperature increase induced by AGN feedback (Kelvin) */
  float AGN_delta_T_desired;

  /*! Number of gas neighbours to heat in a feedback event */
  float num_ngbs_to_heat;

  /* ---- Common conversion factors --------------- */

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;
};

/**
 * @brief Initialise the black hole properties from the parameter file.
 *
 * For the basic black holes neighbour finding properties we use the
 * defaults from the hydro scheme if the users did not provide specific
 * values.
 *
 * @param bp The #black_holes_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void black_holes_props_init(struct black_holes_props *bp,
                                          const struct phys_const *phys_const,
                                          const struct unit_system *us,
                                          struct swift_params *params,
                                          const struct hydro_props *hydro_props,
                                          const struct cosmology *cosmo) {

  /* Read in the basic neighbour search properties or default to the hydro
     ones if the user did not provide any different values */

  /* Kernel properties */
  bp->eta_neighbours = parser_get_opt_param_float(
      params, "BlackHoles:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  bp->h_tolerance = parser_get_opt_param_float(params, "BlackHoles:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  bp->target_neighbours = pow_dimension(bp->eta_neighbours) * kernel_norm;
  const float delta_eta = bp->eta_neighbours * (1.f + bp->h_tolerance);
  bp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(bp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  bp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BlackHoles:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Accretion parameters ---------------------------------- */

  bp->f_Edd = parser_get_param_float(params, "EAGLEAGN:max_eddington_fraction");
  bp->epsilon_r =
      parser_get_param_float(params, "EAGLEAGN:radiative_efficiency");
  bp->epsilon_f =
      parser_get_param_float(params, "EAGLEAGN:coupling_efficiency");

  /* Feedback parameters ---------------------------------- */

  bp->AGN_delta_T_desired =
      parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_K");

  bp->num_ngbs_to_heat =
      parser_get_param_float(params, "EAGLEAGN:AGN_num_ngb_to_heat");

  /* Common conversion factors ----------------------------- */

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);
}

#endif /* SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H */
