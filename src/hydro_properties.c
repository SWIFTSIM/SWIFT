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
#include "equation_of_state.h"
#include "error.h"
#include "hydro.h"
#include "kernel_hydro.h"

#define hydro_props_default_max_iterations 30
#define hydro_props_default_volume_change 1.4f
#define hydro_props_default_h_max FLT_MAX
#define hydro_props_default_h_tolerance 1e-4
#define hydro_props_default_init_temp 0.f
#define hydro_props_default_min_temp 0.f
#define hydro_props_default_H_fraction 0.76

/**
 * @brief Initialize the global properties of the hydro scheme.
 *
 * @param p The #hydro_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 */
void hydro_props_init(struct hydro_props *p,
                      const struct phys_const *phys_const,
                      const struct unit_system *us,
                      struct swift_params *params) {

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
  p->eta_neighbours = 1.0f;
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

  /* Initial temperature */
  p->initial_temperature = parser_get_opt_param_float(
      params, "SPH:initial_temperature", hydro_props_default_init_temp);

  /* Initial temperature */
  p->minimal_temperature = parser_get_opt_param_float(
      params, "SPH:minimal_temperature", hydro_props_default_min_temp);

  if ((p->initial_temperature != 0.) &&
      (p->initial_temperature < p->minimal_temperature))
    error("Initial temperature lower than minimal allowed temperature!");

  /* Hydrogen mass fraction */
  p->hydrogen_mass_fraction = parser_get_opt_param_double(
      params, "SPH:H_mass_fraction", hydro_props_default_H_fraction);

  /* Compute the initial energy (Note the temp. read is in internal units) */
  double u_init = phys_const->const_boltzmann_k / phys_const->const_proton_mass;
  u_init *= p->initial_temperature;
  u_init *= hydro_one_over_gamma_minus_one;

  /* Correct for hydrogen mass fraction */
  double mean_molecular_weight;
  if (p->initial_temperature *
          units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE) >
      1e4)
    mean_molecular_weight = 4. / (8. - 5. * (1. - p->hydrogen_mass_fraction));
  else
    mean_molecular_weight = 4. / (1. + 3. * p->hydrogen_mass_fraction);

  p->initial_internal_energy = u_init / mean_molecular_weight;

  /* Compute the minimal energy (Note the temp. read is in internal units) */
  double u_min = phys_const->const_boltzmann_k / phys_const->const_proton_mass;
  u_min *= p->minimal_temperature;
  u_min *= hydro_one_over_gamma_minus_one;

  /* Correct for hydrogen mass fraction */
  if (p->minimal_temperature *
          units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE) >
      1e4)
    mean_molecular_weight = 4. / (8. - 5. * (1. - p->hydrogen_mass_fraction));
  else
    mean_molecular_weight = 4. / (1. + 3. * p->hydrogen_mass_fraction);

  p->minimal_internal_energy = u_min / mean_molecular_weight;
}

/**
 * @brief Print the global properties of the hydro scheme.
 *
 * @param p The #hydro_props.
 */
void hydro_props_print(const struct hydro_props *p) {

  /* Print equation of state first */
  eos_print(&eos);

  /* Now SPH */
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

  if (p->initial_temperature != hydro_props_default_init_temp)
    message("Initial gas temperature set to %f", p->initial_temperature);

  if (p->minimal_temperature != hydro_props_default_min_temp)
    message("Minimal gas temperature set to %f", p->minimal_temperature);
}

#if defined(HAVE_HDF5)
void hydro_props_print_snapshot(hid_t h_grpsph, const struct hydro_props *p) {

  eos_print_snapshot(h_grpsph, &eos);

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
  io_write_attribute_f(h_grpsph, "Minimal temperature", p->minimal_temperature);
  io_write_attribute_f(h_grpsph, "Initial temperature", p->initial_temperature);
  io_write_attribute_f(h_grpsph, "Initial energy per unit mass",
                       p->initial_internal_energy);
  io_write_attribute_f(h_grpsph, "Hydrogen mass fraction",
                       p->hydrogen_mass_fraction);
}
#endif

/**
 * @brief Write a hydro_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void hydro_props_struct_dump(const struct hydro_props *p, FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct hydro_props), 1, stream,
                       "hydroprops", "hydro props");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void hydro_props_struct_restore(const struct hydro_props *p, FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct hydro_props), 1, stream, NULL,
                      "hydro props");
}
