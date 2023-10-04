/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "gravity_properties.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "mhd.h"
#include "parser.h"
#include "pressure_floor.h"
#include "units.h"

#define hydro_props_default_max_iterations 30
#define hydro_props_default_volume_change 1.4f
#define hydro_props_default_h_max FLT_MAX
#define hydro_props_default_h_min_ratio 0.f
#define hydro_props_default_h_tolerance 1e-4
#define hydro_props_default_init_temp 0.f
#define hydro_props_default_min_temp 0.f
#define hydro_props_default_H_ionization_temperature 1e4

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

  /* ------ Smoothing lengths parameters ---------- */

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

  /* Minimal smoothing length ratio to softening */
  p->h_min_ratio = parser_get_opt_param_float(params, "SPH:h_min_ratio",
                                              hydro_props_default_h_min_ratio);

  /* Temporarily set the minimal softening to 0. */
  p->h_min = 0.f;

  /* Number of iterations to converge h */
  p->max_smoothing_iterations = parser_get_opt_param_int(
      params, "SPH:max_ghost_iterations", hydro_props_default_max_iterations);

  if (p->max_smoothing_iterations <= 10)
    error("The number of smoothing length iterations should be > 10");

  /* ------ Neighbour number definition ------------ */

  /* Non-conventional neighbour number definition */
  p->use_mass_weighted_num_ngb =
      parser_get_opt_param_int(params, "SPH:use_mass_weighted_num_ngb", 0);

  if (p->use_mass_weighted_num_ngb) {
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH) || defined(SHADOWFAX_SPH)
    error("Can't use alternative neighbour definition with this scheme!");
#endif
  }

  /* ------ Time integration parameters ------------ */

  /* Time integration properties */
  p->CFL_condition = parser_get_param_float(params, "SPH:CFL_condition");
  const float max_volume_change = parser_get_opt_param_float(
      params, "SPH:max_volume_change", hydro_props_default_volume_change);
  p->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Initial temperature */
  p->initial_temperature = parser_get_opt_param_float(
      params, "SPH:initial_temperature", hydro_props_default_init_temp);

  if (p->initial_temperature < 0.f)
    error("ERROR: Initial temperature set to a negative value!!!");

  /* ------ Temperature parameters ----------------- */

  /* Minimal temperature */
  p->minimal_temperature = parser_get_opt_param_float(
      params, "SPH:minimal_temperature", hydro_props_default_min_temp);

  if (p->minimal_temperature < 0.f)
    error("ERROR: Minimal temperature set to a negative value!!!");

  if ((p->initial_temperature != 0.) &&
      (p->initial_temperature < p->minimal_temperature))
    error("Initial temperature lower than minimal allowed temperature!");

  /* Neutral to ionized Hydrogen transition temperature */
  p->hydrogen_ionization_temperature =
      parser_get_opt_param_double(params, "SPH:H_ionization_temperature",
                                  hydro_props_default_H_ionization_temperature);

  /* Hydrogen mass fraction */
  const float default_H_fraction =
      1. - phys_const->const_primordial_He_fraction;
  p->hydrogen_mass_fraction = parser_get_opt_param_double(
      params, "SPH:H_mass_fraction", default_H_fraction);

  /* Mean molecular mass for neutral gas */
  p->mu_neutral = 4. / (1. + 3. * p->hydrogen_mass_fraction);

  /* Mean molecular mass for fully ionised gas */
  p->mu_ionised = 4. / (8. - 5. * (1. - p->hydrogen_mass_fraction));

  /* Initialise the implementation-dependent viscosity parameters
   * (see hydro/SCHEME/hydro_parameters.h for this implementation) */
  viscosity_init(params, us, phys_const, &(p->viscosity));

  /* Same for the thermal diffusion parameters */
  diffusion_init(params, us, phys_const, &(p->diffusion));

  /* Compute the initial energy (Note the temp. read is in internal units)
   */
  /* u_init = k_B T_init / (mu m_p (gamma - 1)) */
  double u_init = phys_const->const_boltzmann_k / phys_const->const_proton_mass;
  u_init *= p->initial_temperature;
  u_init *= hydro_one_over_gamma_minus_one;

  /* Correct for hydrogen mass fraction (mu) */
  double mean_molecular_weight;
  if (p->initial_temperature > p->hydrogen_ionization_temperature)
    mean_molecular_weight = 4. / (8. - 5. * (1. - p->hydrogen_mass_fraction));
  else
    mean_molecular_weight = 4. / (1. + 3. * p->hydrogen_mass_fraction);

  p->initial_internal_energy = u_init / mean_molecular_weight;

  /* Compute the minimal energy (Note the temp. read is in internal units) */
  /* u_min = k_B T_min / (mu m_p (gamma - 1)) */
  double u_min = phys_const->const_boltzmann_k / phys_const->const_proton_mass;
  u_min *= p->minimal_temperature;
  u_min *= hydro_one_over_gamma_minus_one;

  /* Correct for hydrogen mass fraction (mu) */
  if (p->minimal_temperature > p->hydrogen_ionization_temperature)
    mean_molecular_weight = 4. / (8. - 5. * (1. - p->hydrogen_mass_fraction));
  else
    mean_molecular_weight = 4. / (1. + 3. * p->hydrogen_mass_fraction);

  p->minimal_internal_energy = u_min / mean_molecular_weight;

  /* ------ Particle splitting parameters ---------- */

  /* Are we doing particle splitting? */
  p->particle_splitting =
      parser_get_opt_param_int(params, "SPH:particle_splitting", 0);

  if (p->particle_splitting) {
    p->particle_splitting_mass_threshold =
        parser_get_param_float(params, "SPH:particle_splitting_mass_threshold");

    p->generate_random_ids = parser_get_opt_param_int(
        params, "SPH:particle_splitting_generate_random_ids", 0);
  }
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

  if (p->use_mass_weighted_num_ngb)
    message("Neighbour number definition: Mass-weighted.");
  else
    message("Neighbour number definition: Unweighted.");

  if (p->h_max != hydro_props_default_h_max)
    message("Maximal smoothing length allowed: %.4f", p->h_max);

  message("Maximal time-bin difference between neighbours: %d",
          time_bin_neighbour_max_delta_bin);

  if (p->max_smoothing_iterations != hydro_props_default_max_iterations)
    message("Maximal iterations in ghost task set to %d (default is %d)",
            p->max_smoothing_iterations, hydro_props_default_max_iterations);

  if (p->initial_temperature != hydro_props_default_init_temp)
    message("Initial gas temperature set to %f", p->initial_temperature);

  if (p->minimal_temperature != hydro_props_default_min_temp)
    message("Minimal gas temperature set to %f", p->minimal_temperature);

  if (p->particle_splitting)
    message("Splitting particles with mass > %e U_M",
            p->particle_splitting_mass_threshold);
  else
    message("No particle splitting");

  /* Print out the implementation-dependent viscosity parameters
   * (see hydro/SCHEME/hydro_parameters.h for this implementation) */
  viscosity_print(&(p->viscosity));

  /* Same for the diffusion */
  diffusion_print(&(p->diffusion));

  /* Same for MHD */
  message("MHD scheme: %s.", MHD_IMPLEMENTATION);

  // MATTHIEU: Temporary location for this planetary SPH i/o business.

#ifdef PLANETARY_SPH
#ifdef PLANETARY_SPH_NO_BALSARA
  message("Planetary SPH: Balsara switch DISABLED");
#else
  message("Planetary SPH: Balsara switch ENABLED");
#endif
#endif
}

#if defined(HAVE_HDF5)
void hydro_props_print_snapshot(hid_t h_grpsph, const struct hydro_props *p) {

  eos_print_snapshot(h_grpsph, &eos);
  pressure_floor_print_snapshot(h_grpsph);

  io_write_attribute_i(h_grpsph, "Dimension", (int)hydro_dimension);
  if (p->use_mass_weighted_num_ngb) {
    io_write_attribute_s(h_grpsph, "Neighbour number definition",
                         "mass-weighted");
  } else {
    io_write_attribute_s(h_grpsph, "Neighbour number definition", "unweighted");
  }
  io_write_attribute_s(h_grpsph, "Scheme", SPH_IMPLEMENTATION);
  io_write_attribute_s(h_grpsph, "MHD Scheme", MHD_IMPLEMENTATION);
  io_write_attribute_s(h_grpsph, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpsph, "Kernel target N_ngb", p->target_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel delta N_ngb", p->delta_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel eta", p->eta_neighbours);
  io_write_attribute_f(h_grpsph, "Kernel gamma", kernel_gamma);
  io_write_attribute_f(h_grpsph, "Smoothing length tolerance", p->h_tolerance);
  io_write_attribute_f(h_grpsph, "Maximal smoothing length [internal units]",
                       p->h_max);
  io_write_attribute_f(h_grpsph, "CFL parameter", p->CFL_condition);
  io_write_attribute_f(h_grpsph, "Volume log(max(delta h))",
                       p->log_max_h_change);
  io_write_attribute_f(h_grpsph, "Volume max change time-step",
                       pow_dimension(expf(p->log_max_h_change)));
  io_write_attribute_i(h_grpsph, "Max ghost iterations",
                       p->max_smoothing_iterations);
  io_write_attribute_f(h_grpsph, "Minimal temperature", p->minimal_temperature);
  io_write_attribute_f(h_grpsph,
                       "Minimal energy per unit mass [internal units]",
                       p->minimal_internal_energy);
  io_write_attribute_f(h_grpsph, "Initial temperature", p->initial_temperature);
  io_write_attribute_f(h_grpsph,
                       "Initial energy per unit mass [internal units]",
                       p->initial_internal_energy);
  io_write_attribute_f(h_grpsph, "Hydrogen mass fraction",
                       p->hydrogen_mass_fraction);
  io_write_attribute_f(h_grpsph, "Hydrogen ionization transition temperature",
                       p->hydrogen_ionization_temperature);
  io_write_attribute_i(h_grpsph,
                       "Maximal time-bin difference between neighbours",
                       time_bin_neighbour_max_delta_bin);

  /* Write out the implementation-dependent viscosity parameters
   * (see hydro/SCHEME/hydro_parameters.h for this implementation) */
  viscosity_print_snapshot(h_grpsph, &(p->viscosity));

  /* Same for the diffusion */
  diffusion_print_snapshot(h_grpsph, &(p->diffusion));
}
#endif

/**
 * @brief Initialises a hydro_props struct with somewhat useful values for
 *        the automated test suite. This is not intended for production use,
 *        but rather to fill for the purposes of mocking.
 *
 * @param p the struct
 */
void hydro_props_init_no_hydro(struct hydro_props *p) {

  p->eta_neighbours = 1.2348;
  p->h_tolerance = hydro_props_default_h_tolerance;
  p->target_neighbours = pow_dimension(p->eta_neighbours) * kernel_norm;
  const float delta_eta = p->eta_neighbours * (1.f + p->h_tolerance);
  p->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(p->eta_neighbours)) *
      kernel_norm;
  p->h_max = hydro_props_default_h_max;
  p->h_min = 0.f;
  p->h_min_ratio = hydro_props_default_h_min_ratio;
  p->max_smoothing_iterations = hydro_props_default_max_iterations;
  p->CFL_condition = 0.1;
  p->log_max_h_change = logf(powf(1.4, hydro_dimension_inv));

  /* These values are inconsistent and in a production run would probably lead
     to a crash. Again, this function is intended for mocking use in unit tests
     and is _not_ to be used otherwise! */
  p->minimal_temperature = hydro_props_default_min_temp;
  p->minimal_internal_energy = 0.f;
  p->initial_temperature = hydro_props_default_init_temp;
  p->initial_internal_energy = 0.f;

  p->hydrogen_mass_fraction = 0.755;
  p->hydrogen_ionization_temperature =
      hydro_props_default_H_ionization_temperature;

  /* Initialise the implementation-dependent viscosity parameters
   * (see hydro/SCHEME/hydro_parameters.h for this implementation) */
  viscosity_init_no_hydro(&(p->viscosity));

  /* Same for the diffusion */
  diffusion_init_no_hydro(&(p->diffusion));
}

/**
 * @brief Update the global properties of the hydro scheme for that time-step.
 *
 * @param p The properties to update.
 * @param gp The properties of the gravity scheme.
 * @param cosmo The cosmological model.
 */
void hydro_props_update(struct hydro_props *p, const struct gravity_props *gp,
                        const struct cosmology *cosmo) {

  /* Update the minimal allowed smoothing length
   *
   * We follow Gadget here and demand that the kernel support (h * gamma)
   * is a fixed fraction of the radius at which the softened forces
   * recover a Newtonian behaviour (i.e. 2.8 * Plummer equivalent softening
   * in the case of a cubic spline kernel). */
  p->h_min = p->h_min_ratio * gp->epsilon_baryon_cur / kernel_gamma;
}

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
