/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* This object's header. */
#include "hydro_properties.h"

/* Standard headers */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "dimension.h"
#include "error.h"
#include "hydro.h"
#include "kernel_hydro.h"

#define hydro_props_default_max_iterations 30
#define hydro_props_default_volume_change 1.4f
#define hydro_props_default_h_max FLT_MAX
#define hydro_props_default_h_tolerance 1e-4

void hydro_props_init(struct hydro_props *p,
                      const struct swift_params *params) {

  /* Kernel properties */
  p->eta_neighbours = parser_get_param_float(params, "SPH:resolution_eta");

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  p->h_tolerance = parser_get_opt_param_float(params, "SPH:h_tolerance",
                                              hydro_props_default_h_tolerance);

  /* Get derived properties */
  p->target_neighbours = pow_dimension(p->eta_neighbours) * kernel_norm;
  const float delta_eta = p->eta_neighbours * (1.f + p->h_tolerance);
  p->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(p->eta_neighbours)) *
      kernel_norm;

#ifdef SHADOWFAX_SPH
  /* change the meaning of target_neighbours and delta_neighbours */
  p->target_neighbours = 1.0f;
  p->delta_neighbours = 0.0f;
#endif

  /* Maximal smoothing length */
  p->h_max = parser_get_opt_param_float(params, "SPH:h_max",
                                        hydro_props_default_h_max);

  /* Number of iterations to converge h */
  p->max_smoothing_iterations = parser_get_opt_param_int(
      params, "SPH:max_ghost_iterations", hydro_props_default_max_iterations);

  /* Time integration properties */
  p->CFL_condition = parser_get_param_float(params, "SPH:CFL_condition");
  const float max_volume_change = parser_get_opt_param_float(
      params, "SPH:max_volume_change", hydro_props_default_volume_change);
  p->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));
}

void hydro_props_print(const struct hydro_props *p) {

#if defined(EOS_IDEAL_GAS)
  message("Equation of state: Ideal gas.");
#elif defined(EOS_ISOTHERMAL_GAS)
  message(
      "Equation of state: Isothermal with internal energy "
      "per unit mass set to %f.",
      const_isothermal_internal_energy);
#endif

  message("Adiabatic index gamma: %f.", hydro_gamma);

  message("Hydrodynamic scheme: %s in %dD.", SPH_IMPLEMENTATION,
          (int)hydro_dimension);

  message("Hydrodynamic kernel: %s with eta=%f (%.2f neighbours).", kernel_name,
          p->eta_neighbours, p->target_neighbours);

  message("Hydrodynamic relative tolerance in h: %.5f (+/- %.4f neighbours).",
          p->h_tolerance, p->delta_neighbours);

  message("Hydrodynamic integration: CFL parameter: %.4f.", p->CFL_condition);

  message(
      "Hydrodynamic integration: Max change of volume: %.2f "
      "(max|dlog(h)/dt|=%f).",
      pow_dimension(expf(p->log_max_h_change)), p->log_max_h_change);

  if (p->h_max != hydro_props_default_h_max)
    message("Maximal smoothing length allowed: %.4f", p->h_max);

  if (p->max_smoothing_iterations != hydro_props_default_max_iterations)
    message("Maximal iterations in ghost task set to %d (default is %d)",
            p->max_smoothing_iterations, hydro_props_default_max_iterations);
}

#if defined(HAVE_HDF5)
void hydro_props_print_snapshot(hid_t h_grpsph, const struct hydro_props *p) {

  io_write_attribute_f(h_grpsph, "Adiabatic index", hydro_gamma);
  io_write_attribute_i(h_grpsph, "Dimension", (int)hydro_dimension);
  io_write_attribute_s(h_grpsph, "Scheme", SPH_IMPLEMENTATION);
  io_write_attribute_s(h_grpsph, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpsph, "Kernel target N_ngb", p->target_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel delta N_ngb", p->delta_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel eta", p->eta_neighbours);
  io_write_attribute_f(h_grpsph, "Smoothing length tolerance", p->h_tolerance);
  io_write_attribute_f(h_grpsph, "Maximal smoothing length", p->h_max);
  io_write_attribute_f(h_grpsph, "CFL parameter", p->CFL_condition);
  io_write_attribute_f(h_grpsph, "Volume log(max(delta h))",
                       p->log_max_h_change);
  io_write_attribute_f(h_grpsph, "Volume max change time-step",
                       pow_dimension(expf(p->log_max_h_change)));
  io_write_attribute_i(h_grpsph, "Max ghost iterations",
                       p->max_smoothing_iterations);
}
#endif
