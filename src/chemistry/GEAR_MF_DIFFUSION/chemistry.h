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
 ******************************************************************************/
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_H

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry.h
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chemistry_flux.h"
#include "chemistry_getters.h"
#include "chemistry_gradients.h"
#include "chemistry_setters.h"
#include "chemistry_struct.h"
#include "chemistry_timesteps.h"
#include "error.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* Some constants */
#define FILTERING_SMOOTHING_FACTOR 0.8
#define DEFAULT_DIFFUSION_NORMALISATION 1
#define DEFAULT_PSI_RIEMANN_SOLVER 0.1
#define DEFAULT_EPSILON_RIEMANN_SOLVER 0.5
#define DEFAULT_USE_SUPERTIMESTEPPING 0
#define DEFAULT_N_SUBSTEPS 5
#define DEFAULT_NU_SUPERTIMESTEPPPING 0.04
#define DEFAULT_C_CFL_CHEMISTRY_SUPERTIMESTEPPPING 0.4

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  message("Chemistry function is 'Gear MF hyperbolic diffusion'.");
#else
  message("Chemistry function is 'Gear MF diffusion'.");
#endif
}

/**
 * @brief Read the solar abundances and scale with them the initial
 * metallicities.
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_scale_initial_metallicities(
    struct swift_params* parameter_file, struct chemistry_global_data* data) {

#ifndef HAVE_HDF5
  error("Cannot scale the solar abundances without HDF5");
#endif
  /* Scale the initial metallicities */
  char txt[DESCRIPTION_BUFFER_SIZE] = "Scaling initial metallicities by:";
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    data->initial_metallicities[i] *= data->solar_abundances[i];
    char tmp[10];
    sprintf(tmp, " %.2g", data->solar_abundances[i]);
    strcat(txt, tmp);
  }

  if (engine_rank == 0) {
    message("%s", txt);
  }
}

/**
 * @brief Read the solar abundances and scale with them the initial
 * metallicities.
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_read_solar_abundances(
    struct swift_params* parameter_file, struct chemistry_global_data* data) {
#if defined(HAVE_HDF5)

  /* Get the yields table */
  char filename[DESCRIPTION_BUFFER_SIZE];
  parser_get_param_string(parameter_file, "GEARFeedback:yields_table",
                          filename);

  /* Open file. */
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s.\n", filename);

  /* Open group. */
  hid_t group_id = H5Gopen(file_id, "Data", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group Data.\n");

  /* Read the data */
  io_read_array_attribute(group_id, "SolarMassAbundances", FLOAT,
                          data->solar_abundances, GEAR_CHEMISTRY_ELEMENT_COUNT);

  /* Close group */
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file.");

#else
  message("Cannot read the solar abundances without HDF5");
#endif
}

/**
 * @brief Get the name of the element i.
 *
 * @param sm The #stellar_model.
 * @param i The element indice.
 */
static INLINE const char* chemistry_get_element_name(
    const struct chemistry_global_data* data, int i) {

  return data->elements_name + i * GEAR_LABELS_SIZE;
}

/**
 * @brief Get the index of the element .
 *
 * @param sm The #stellar_model.
 * @param element_name The element name.
 */
static INLINE int chemistry_get_element_index(
    const struct chemistry_global_data* data, const char* element_name) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (strcmp(chemistry_get_element_name(data, i), element_name) == 0)
      return i;
  }
  error("Chemical element %s not found !", element_name);

  return -1;
}

/**
 * @brief Read the name of all the elements present in the tables.
 * It is nearly a copy/paste of stellar_evolution_read_elements
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_read_elements(struct swift_params* params,
                                           struct chemistry_global_data* data) {

  /* Read the elements from the parameter file. */
  int nval = -1;
  char** elements;
  parser_get_param_string_array(params, "GEARFeedback:elements", &nval,
                                &elements);

  /* Check that we have the correct number of elements. */
  if (nval != GEAR_CHEMISTRY_ELEMENT_COUNT - 1) {
    error(
        "You need to provide %i elements but found %i. "
        "If you wish to provide a different number of elements, "
        "you need to compile with --with-chemistry=GEAR_N where N "
        "is the number of elements + 1.",
        GEAR_CHEMISTRY_ELEMENT_COUNT, nval);
  }

  /* Copy the elements into the stellar model. */
  for (int i = 0; i < nval; i++) {
    if (strlen(elements[i]) >= GEAR_LABELS_SIZE) {
      error("Element name '%s' too long", elements[i]);
    }
    strcpy(data->elements_name + i * GEAR_LABELS_SIZE, elements[i]);
  }

  /* Cleanup. */
  parser_free_param_string_array(nval, elements);

  /* Add the metals to the end. */
  strcpy(data->elements_name +
             (GEAR_CHEMISTRY_ELEMENT_COUNT - 1) * GEAR_LABELS_SIZE,
         "Metals");

  /* Check the elements */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    for (int j = i + 1; j < GEAR_CHEMISTRY_ELEMENT_COUNT; j++) {
      const char* el_i = chemistry_get_element_name(data, i);
      const char* el_j = chemistry_get_element_name(data, j);
      if (strcmp(el_i, el_j) == 0) {
        error("You need to provide each element only once (%s).", el_i);
      }
    }
  }

  /* Check that iron is at index 0 */
  if (chemistry_get_element_index(data, "Fe") != 0)
    error("Element Fe must be at index 0 !");
}

/**
 * @brief Initialises the chemistry properties.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* read parameters */
  const float initial_metallicity = parser_get_param_float(
      parameter_file, "GEARChemistry:initial_metallicity");

  if (initial_metallicity < 0) {
    message("Setting the initial metallicity from the snapshot.");
  } else {
    message("Setting the initial metallicity from the parameter file.");
  }

  /* Set the initial metallicities */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    data->initial_metallicities[i] = initial_metallicity;
  }

  /* Check if need to scale the initial metallicity */
  const int scale_metallicity = parser_get_opt_param_int(
      parameter_file, "GEARChemistry:scale_initial_metallicity", 0);

  /* Scale the metallicities if required */
  if (scale_metallicity) {
    /* Read the solar abundances */
    chemistry_read_solar_abundances(parameter_file, data);
    chemistry_read_elements(parameter_file, data);

    /* Scale the solar abundances */
    chemistry_scale_initial_metallicities(parameter_file, data);
  }
  /* We do not care about the solar abundances without feedback */
#ifdef FEEDBACK_GEAR
  else {
    chemistry_read_solar_abundances(parameter_file, data);
    chemistry_read_elements(parameter_file, data);
  }
#endif

  /***************************************************************************/
  data->diffusion_coefficient = parser_get_opt_param_double(
      parameter_file, "GEARChemistry:diffusion_coefficient",
      DEFAULT_DIFFUSION_NORMALISATION);

  char temp[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(parameter_file, "GEARChemistry:diffusion_mode", temp);
  if (strcmp(temp, "Isotropic_constant") == 0)
    data->diffusion_mode = isotropic_constant;
  else if (strcmp(temp, "Smagorinsky") == 0)
    data->diffusion_mode = isotropic_smagorinsky;
  else if (strcmp(temp, "Gradient") == 0)
    data->diffusion_mode = anisotropic_gradient;
  else
    error(
        "The chemistry diffusion mode must be Isotropic_constant, Smagorinsky "
        "or Gradient "
        " not %s",
        temp);

  /***************************************************************************/
  /* Read parameters for the Riemann solver */
  data->hll_riemann_solver_psi = parser_get_opt_param_float(
      parameter_file, "GEARChemistry:hll_riemann_solver_psi",
      DEFAULT_PSI_RIEMANN_SOLVER);

  if ((data->hll_riemann_solver_psi < 0)) {
    error("hll_riemann_solver_psi must be positive!");
  }

  data->hll_riemann_solver_epsilon = parser_get_opt_param_float(
      parameter_file, "GEARChemistry:hll_riemann_solver_epsilon",
      DEFAULT_EPSILON_RIEMANN_SOLVER);

  if ((data->hll_riemann_solver_epsilon > 1) ||
      (data->hll_riemann_solver_epsilon < 0)) {
    error("hll_riemann_solver_epsilon must respect 0 <= epsilon <= 1!");
  }

  /***************************************************************************/
  /* Supertimestepping */
  data->use_supertimestepping = parser_get_opt_param_int(
      parameter_file, "GEARChemistry:use_supertimestepping",
      DEFAULT_USE_SUPERTIMESTEPPING);

  data->N_substeps = parser_get_opt_param_int(
      parameter_file, "GEARChemistry:N_substeps", DEFAULT_N_SUBSTEPS);

  data->nu = parser_get_opt_param_float(parameter_file,
                                        "GEARChemistry:nu_supertimestepping",
                                        DEFAULT_NU_SUPERTIMESTEPPPING);

  data->C_CFL_chemistry = parser_get_opt_param_float(
      parameter_file, "GEARChemistry:C_CFL_chemistry",
      DEFAULT_C_CFL_CHEMISTRY_SUPERTIMESTEPPPING);

  /***************************************************************************/
  /* Hyperbolic diffusion */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  if (data->diffusion_mode == isotropic_constant) {
    data->tau =
        parser_get_param_double(parameter_file, "GEARChemistry:tau");
  }

  /* This is used for testing only */
  data->riemann_solver = parser_get_opt_param_int(
      parameter_file, "GEARChemistry:riemann_solver", 1);
#endif

  /***************************************************************************/
  /* Print the parameters we use */
  if (engine_rank == 0) {
    message("Diffusion mode:             %u", data->diffusion_mode);
    message("Diffusion coefficient:      %e", data->diffusion_coefficient);
    message("HLL Riemann solver psi:     %e", data->hll_riemann_solver_psi);
    message("HLL Riemann solver epsilon: %e", data->hll_riemann_solver_epsilon);
    message("Use supertimestepping:      %d", data->use_supertimestepping);
    message("N_substeps:                 %d", data->N_substeps);
    message("nu:                         %e", data->nu);
  }
}

/**
 * @brief Computes the chemistry-related timestep constraints.
 *
 * Parabolic constraint proportional to h^2 or supertimestepping.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param chem_data The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* chem_data,
    const struct part* restrict p) {

  if (chem_data->use_supertimestepping) {
    /* For supertimestepping, use the advective timestep eq (D1) */
    return chemistry_compute_advective_supertimestep(p, chem_data, cosmo);
  } else {
    /* Without supertimestepping, use the parabolic timestep eq (15) */
    return chemistry_compute_parabolic_timestep(p, chem_data, cosmo);
  }
}

/**
 * @brief Do supertimestepping following Alexiades, Amiez and Gremaud (1996)
 * and Hopkins (2017) modifications for individual timestep scheme.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 * @param dt_part Minimal timestep for the other physical processes (hydro,
 * MHD, rt, gravity, ...).
 * @param time_base The system's minimal time-step.
 * @param ti_current The current time on the integer time-line.
 */
__attribute__((always_inline)) INLINE static float chemistry_supertimestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, struct part* restrict p,
    const float dt_part, const double time_base, const double ti_current) {

  struct chemistry_part_data* chd = &p->chemistry_data;

  /* Do not use supertimestepping in the fake timestep */
  if (cd->use_supertimestepping && p->time_bin != 0) {

    /* If we need to start a new cycle, compute the new explicit timestep */
    if (chd->timesteps.current_substep == 0) {
      /* Get the explicit timestep from eq (15) */
      chd->timesteps.explicit_timestep =
          chemistry_compute_parabolic_timestep(p, cd, cosmo);
    }

    /* Compute the current substep */
    double substep = chemistry_compute_subtimestep(
        p, cd, chd->timesteps.current_substep + 1);
    double timestep_to_use = substep;

    if (dt_part <= substep || chd->timesteps.current_substep > 0) {
      /* If the part timestep is smaller, other constraints beat the
         substep. We can safely iterate the substep.
         If we are in the middle of a supertimestepping cycle, iterate. */
      ++chd->timesteps.current_substep;

      /* Increment substep and loop if it cycles fully */
      if (chd->timesteps.current_substep >= cd->N_substeps) {
        chd->timesteps.current_substep = 0;
      }
    } else {
      /* If we are at the beginning of a new cycle (current_substep = 0) and
      the next substep matters for starting a new cycle (dt_part > substep),
      think whether to start a cycle */

      /* Apply cosmology correction (This is 1 if non-cosmological) */
      const double dt_pred = substep * cosmo->time_step_factor;

      /* Convert to integer time.
         Note: This function is similar to Gizmo code checking for
         synchronization. */
      const integertime_t dti_allowed = chemistry_make_integer_timestep(
          dt_pred, p->time_bin, p->limiter_data.min_ngb_time_bin, ti_current,
          1.0 / time_base);

      /* Now convert this back to a physical timestep */
      const timebin_t bin_allowed = get_time_bin(dti_allowed);
      const double dt_allowed =
          get_timestep(bin_allowed, time_base) / cosmo->time_step_factor;

      if (substep > 1.5 * dt_allowed) {
        /* The next allowed timestep is too small to fit the big step part of
        the supertimestepping cycle. Do not waste our timestep (which will put
        us into a lower bin and defeat the supertimestep purpose of using bigger
        timesteps) and use the safe parabolic timestep until the next desired
        time-bin synchs up */
        timestep_to_use = chemistry_compute_parabolic_timestep(p, cd, cosmo);
      } else {
        /* We can jump up in bins to use our super-step => begin the cycle */
        ++chd->timesteps.current_substep;

        /* Increment substep and loop if it cycles fully */
        if (chd->timesteps.current_substep >= cd->N_substeps) {
          chd->timesteps.current_substep = 0;
        }
      } /* substep > 1.5*dt_allowed */
    }
    return timestep_to_use;
  } else { /* No supertimestepping */
    return FLT_MAX;
  }
}

/**
 * @brief Finishes geometry and density calculations.
 *
 * This function requires the #hydro_end_density to have been called.
 *
 * @param p The particle to act upon.
 * @param cd The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {

  /*****************************************/
  /* Finish computations on the filtered quantities */
  /* Multiply by the smoothing factor */
  p->chemistry_data.filtered.rho_v[0] *= FILTERING_SMOOTHING_FACTOR;
  p->chemistry_data.filtered.rho_v[1] *= FILTERING_SMOOTHING_FACTOR;
  p->chemistry_data.filtered.rho_v[2] *= FILTERING_SMOOTHING_FACTOR;

  /* Add self term (is it needed for rho? the formula does not include it) */
  /* p->chemistry_data.filtered.rho += hydro_get_mass(p) * kernel_root; */
  p->chemistry_data.filtered.rho_v[0] +=
      hydro_get_comoving_density(p) * p->v[0];
  p->chemistry_data.filtered.rho_v[1] +=
      hydro_get_comoving_density(p) * p->v[1];
  p->chemistry_data.filtered.rho_v[2] +=
      hydro_get_comoving_density(p) * p->v[2];

  /* Insert missing h-factor to rho. For rho*v, notice that the h-factors were
     already given in the density loop since they depend on bar{h_ij}. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  p->chemistry_data.filtered.rho *= h_inv_dim;

  /*****************************************/
  /* Finish computations on the MF geometry quantities */

  /* Store a copy of the matrix to avoid modifying it when computing its
     inverse */
  float matrix_E[3][3];
  matrix_E[0][0] = p->geometry.matrix_E[0][0];
  matrix_E[0][1] = p->geometry.matrix_E[0][1];
  matrix_E[0][2] = p->geometry.matrix_E[0][2];

  matrix_E[1][0] = p->geometry.matrix_E[1][0];
  matrix_E[1][1] = p->geometry.matrix_E[1][1];
  matrix_E[1][2] = p->geometry.matrix_E[1][2];

  matrix_E[2][0] = p->geometry.matrix_E[2][0];
  matrix_E[2][1] = p->geometry.matrix_E[2][1];
  matrix_E[2][2] = p->geometry.matrix_E[2][2];

  const float condition_number_E_inv =
      matrix_E[0][0] * matrix_E[0][0] + matrix_E[0][1] * matrix_E[0][1] +
      matrix_E[0][2] * matrix_E[0][2] + matrix_E[1][0] * matrix_E[1][0] +
      matrix_E[1][1] * matrix_E[1][1] + matrix_E[1][2] * matrix_E[1][2] +
      matrix_E[2][0] * matrix_E[2][0] + matrix_E[2][1] * matrix_E[2][1] +
      matrix_E[2][2] * matrix_E[2][2];

  if (invert_dimension_by_dimension_matrix(matrix_E) != 0) {
    /* something went wrong in the inversion; force bad condition number */
    p->chemistry_data.geometry_condition_number =
        const_gizmo_max_condition_number + 1.0f;
  } else {
    const float condition_number_E =
        matrix_E[0][0] * matrix_E[0][0] + matrix_E[0][1] * matrix_E[0][1] +
        matrix_E[0][2] * matrix_E[0][2] + matrix_E[1][0] * matrix_E[1][0] +
        matrix_E[1][1] * matrix_E[1][1] + matrix_E[1][2] * matrix_E[1][2] +
        matrix_E[2][0] * matrix_E[2][0] + matrix_E[2][1] * matrix_E[2][1] +
        matrix_E[2][2] * matrix_E[2][2];

    p->chemistry_data.geometry_condition_number =
        hydro_dimension_inv *
        sqrtf(condition_number_E * condition_number_E_inv);
  }

  /* Check that the metal masses are physical */
  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {
    const double m_metal_old = p->chemistry_data.metal_mass[g];

#ifdef SWIFT_DEBUG_CHECKS
    if (p->geometry.volume == 0.) {
      error("Volume is 0!");
    }
#endif
    chemistry_check_unphysical_state(&p->chemistry_data.metal_mass[g],
                                     m_metal_old, hydro_get_mass(p),
                                     /*callloc=*/3, /*element*/ g);
  }

  /* Sanity check on the total metal mass */
  chemistry_check_unphysical_total_metal_mass(p, 0);
}

/**
 * @brief Finishes the gradient calculations.
 *
 * Just a wrapper around chemistry_gradients_finalize, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param p The particle to act upon.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_gradient(
    struct part* p, const struct chemistry_global_data* cd) {
  chemistry_gradients_finalise(p, cd);
}

/**
 * @brief Prepare a particle for the force calculations.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_prepare_force(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo, const float dt_alpha, const float dt_therm,
    const struct chemistry_global_data* cd) {
  p->chemistry_data.flux_dt = dt_therm;

  /* Update the diffusion coefficient for the new loop */
  p->chemistry_data.kappa =
      chemistry_compute_diffusion_coefficient(p, cd, cosmo);

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  p->chemistry_data.tau = chemistry_compute_physical_tau(p, cd, cosmo);
#endif
}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 * @param with_cosmology Are we running with the cosmology?
 * @param time Current time of the simulation.
 * @param dt Time step (in physical units).
 * @param chem_data The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* restrict p, const struct cosmology* cosmo,
    const int with_cosmology, const double time, const double dt,
    const struct chemistry_global_data* cd) {

  /* Update active particles */
  struct chemistry_part_data* chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    double flux;
    chemistry_get_fluxes(p, i, &flux);

    /* Update the conserved variable */
    chd->metal_mass[i] += flux;

    /* Update the diffused metal mass */
    chd->diffused_metal_mass[i] += flux;
  }

  /* Reset the fluxes now that they have been applied */
  chemistry_reset_chemistry_fluxes(p);

  /* Invalidate the particle time-step. It is considered to be inactive until
     dt is set again in hydro_prepare_force() */
  chd->flux_dt = -1.0f;

  /* Element-wise sanity checks */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    const double m_metal_old = chd->metal_mass[i];
    chemistry_check_unphysical_state(&chd->metal_mass[i], m_metal_old,
                                     hydro_get_mass(p), /*callloc=*/2,
                                     /*element*/ i);
  }

  /* Sanity check on the total metal mass */
  chemistry_check_unphysical_total_metal_mass(p, 1);

  /* Store the density of the current timestep for the next timestep */
  p->chemistry_data.rho_prev = chemistry_get_comoving_density(p);
  p->chemistry_data.filtered.rho_prev = p->chemistry_data.filtered.rho;

  /* Take care of the case where \bar{rho_prev} = 0 */
  if (p->chemistry_data.filtered.rho_prev == 0) {
    p->chemistry_data.filtered.rho_prev = p->chemistry_data.rho_prev;
  }
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0
 * ngbs.
 *
 * Note: When the #part has 0 ngbs, the chemistry_end_density() is not called.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_has_no_neighbours(struct part* restrict p,
                                 struct xpart* restrict xp,
                                 const struct chemistry_global_data* cd,
                                 const struct cosmology* cosmo) {

  /* Re-set geometry problematic values */
  p->geometry.volume = 1.0f;
  p->geometry.matrix_E[0][0] = 1.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 1.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 1.0f;
  p->chemistry_data.geometry_condition_number = hydro_dimension_inv * 1.f;

  /* Reset the centroid to disable MFV velocity corrections for this particle */
  fvpm_reset_centroids(p);

  /* Set density loop variables to meaningful values */
  p->chemistry_data.filtered.rho = hydro_get_comoving_density(p);
  p->chemistry_data.filtered.rho_v[0] = hydro_get_comoving_density(p) * p->v[0];
  p->chemistry_data.filtered.rho_v[1] = hydro_get_comoving_density(p) * p->v[1];
  p->chemistry_data.filtered.rho_v[2] = hydro_get_comoving_density(p) * p->v[2];
}

/**
 * @brief Prepares a particle for the diffusionx calculations.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various diffusion tasks.
 *
 * Warning: DON'T call p->rho here. It is basically 0.
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  cpd->filtered.rho = 0.0;
  cpd->filtered.rho_v[0] = 0.0;
  cpd->filtered.rho_v[1] = 0.0;
  cpd->filtered.rho_v[2] = 0.0;

  /* Init the gradient for the next loops */
  chemistry_gradients_init(p);

  /* Initialize time step criterion variables */
  cpd->timestepvars.vmax = 0.;
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param cosmo The #cosmology.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* cd, struct part* restrict p,
    struct xpart* restrict xp) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (cd->initial_metallicities[i] < 0) {
      /* Use the value from the IC. We are reading the metal mass fraction. */
      p->chemistry_data.metal_mass[i] *= hydro_get_mass(p);
    } else {
      /* Use the value from the parameter file */
      p->chemistry_data.metal_mass[i] =
          cd->initial_metallicities[i] * hydro_get_mass(p);
    }
  }

  /* Init the part chemistry data */
  chemistry_init_part(p, cd);

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  p->geometry.wcorr = 1.0f;

  /* Supertimestepping ----*/
  /* Set the substep to 0 so that we can compute everything in timestep for the
     first time. */
  p->chemistry_data.timesteps.current_substep = 0;

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.hyperbolic_flux[i].F_diff_pred[0] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].F_diff_pred[1] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].F_diff_pred[2] = 0.0;

    p->chemistry_data.hyperbolic_flux[i].F_diff[0] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].F_diff[1] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].F_diff[2] = 0.0;

    p->chemistry_data.hyperbolic_flux[i].dF_dt[0] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].dF_dt[1] = 0.0;
    p->chemistry_data.hyperbolic_flux[i].dF_dt[2] = 0.0;

    p->chemistry_data.timestepvars.vmax = 0.0;
  }
#endif
}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data* data, struct spart* restrict sp) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Bug fix (26.07.2024): Check that the initial me metallicities are non
       negative. */
    if (data->initial_metallicities[i] >= 0) {
      /* Use the value from the parameter file */
      sp->chemistry_data.metal_mass_fraction[i] =
          data->initial_metallicities[i];
    }
    /* else : Use the value from the IC. We are reading the metal mass
     fraction. So do not overwrite the metallicities */
  }
}

/* Add chemistry first init sink ? */

/**
 * @brief Sets the chemistry properties of the sink particles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sink Pointer to the sink particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_sink(
    const struct chemistry_global_data* data, struct sink* restrict sink) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Use the value from the parameter file */
    if (data->initial_metallicities[i] >= 0) {
      sink->chemistry_data.metal_mass_fraction[i] =
          data->initial_metallicities[i];
    }
    /* else : read the metallicities from the ICs. */
  }
}

/**
 * @brief Copies the chemistry properties of the gas particle over to the
 * star particle.
 *
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 */
INLINE static void chemistry_copy_star_formation_properties(
    struct part* p, const struct xpart* xp, struct spart* sp) {

  /* gas mass after update */
  float mass = hydro_get_mass(p);

  /* Store the chemistry struct in the star particle */
  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    sp->chemistry_data.metal_mass_fraction[k] =
        p->chemistry_data.metal_mass[k] / mass;

    /* Remove the metals taken by the star. */
    p->chemistry_data.metal_mass[k] *= mass / (mass + sp->mass);
  }
}

/**
 * @brief Copies the chemistry properties of the sink particle over to the
 * stellar particle.
 *
 * @param sink the sink particle with its properties.
 * @param sp the new star particles.
 */
INLINE static void chemistry_copy_sink_properties_to_star(struct sink* sink,
                                                          struct spart* sp) {

  /* Store the chemistry struct in the star particle */
  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    sp->chemistry_data.metal_mass_fraction[k] =
        sink->chemistry_data.metal_mass_fraction[k];
  }
}

/**
 * @brief Copies the chemistry properties of the gas particle over to the
 * sink particle.
 *
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sink the new created star particle with its properties.
 */
INLINE static void chemistry_copy_sink_properties(const struct part* p,
                                                  const struct xpart* xp,
                                                  struct sink* sink) {

  /* Store the chemistry struct in the star particle */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    sink->chemistry_data.metal_mass_fraction[i] =
        p->chemistry_data.metal_mass[i] / hydro_get_mass(p);
  }
}

/**
 * @brief Add the chemistry data of a sink particle to a sink.
 *
 * @param si_data The sink data to add to.
 * @param sj_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_sink_to_sink(
    struct sink* si, const struct sink* sj, const double mi_old) {

  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    double mk = si->chemistry_data.metal_mass_fraction[k] * mi_old +
                sj->chemistry_data.metal_mass_fraction[k] * sj->mass;

    si->chemistry_data.metal_mass_fraction[k] = mk / si->mass;
  }
}

/**
 * @brief Add the chemistry data of a gas particle to a sink.
 *
 * @param s The #sink to add to.
 * @param p The gas #part to use.
 * @param ms_old The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_sink(
    struct sink* s, const struct part* p, const double ms_old) {
  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    double mk = s->chemistry_data.metal_mass_fraction[k] * ms_old +
                p->chemistry_data.metal_mass[k];

    s->chemistry_data.metal_mass_fraction[k] = mk / s->mass;
  }
}

/**
 * @brief Transfer chemistry data of a gas particle to a black hole.
 *
 * In contrast to `chemistry_add_part_to_bpart`, only a fraction of the
 * masses stored in the gas particle are transferred here. Absolute masses
 * of the gas particle are adjusted as well.
 * Black holes don't store fractions so we need to add element masses.
 *
 * Nothing to do here.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param nibble_mass The mass to be removed from the gas particle.
 * @param nibble_fraction The fraction of the (original) mass of the gas
 *        particle that is removed.
 */
__attribute__((always_inline)) INLINE static void
chemistry_transfer_part_to_bpart(struct chemistry_bpart_data* bp_data,
                                 struct chemistry_part_data* p_data,
                                 const double nibble_mass,
                                 const double nibble_fraction) {
  error("No BH yet in GEAR");
}

/**
 * @brief Add the chemistry data of a black hole to another one.
 *
 * Nothing to do here.
 *
 * @param bp_data The black hole data to add to.
 * @param swallowed_data The black hole data to use.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_bpart_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_bpart_data* swallowed_data) {
  error("No BH yet in GEAR");
}

/**
 * @brief Split the metal content of a particle into n pieces
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void chemistry_split_part(
    struct part* p, const double n) {
  error("Loic: to be implemented");
}

/**
 * @brief Extra chemistry operations to be done during the drift.
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The current cosmological model.
 * @param chem_data The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_predict_extra(
    struct part* p, struct xpart* xp, float dt_drift, float dt_therm,
    const struct cosmology* cosmo,
    const struct chemistry_global_data* chem_data) {

  struct chemistry_part_data* chd = &p->chemistry_data;

  /* Predict rho_prev using the density predicted by the hydro. */
  chd->rho_prev = p->rho;

  /* Predict filtered density. Notice that h was predicted before by
     hydro_predict_extra() */
  float h_bar_inv = 1 / (p->h * kernel_gamma);
  const float w1 = p->force.h_dt * h_bar_inv * dt_drift;

  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    const float expf_approx =
        approx_expf(w2); /* 4th order expansion of exp(w) */
    chd->filtered.rho *= expf_approx;
  } else {
    const float expf_exact = expf(w2);
    chd->filtered.rho *= expf_exact;
  }
  chd->filtered.rho_prev = chd->filtered.rho;

  /* Update diffusion coefficient */
  chd->kappa = chemistry_compute_diffusion_coefficient(p, chem_data, cosmo);

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  /* Compute the predicted flux */
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    chd->hyperbolic_flux[m].F_diff_pred[0] =
        chd->hyperbolic_flux[m].F_diff[0] +
        0.5*dt_therm * chd->hyperbolic_flux[m].dF_dt[0];
    chd->hyperbolic_flux[m].F_diff_pred[1] =
        chd->hyperbolic_flux[m].F_diff[1] +
        0.5*dt_therm * chd->hyperbolic_flux[m].dF_dt[1];
    chd->hyperbolic_flux[m].F_diff_pred[2] =
        chd->hyperbolic_flux[m].F_diff[2] +
        0.5*dt_therm * chd->hyperbolic_flux[m].dF_dt[2];
  }
#endif

  /* Update inactive particles that are drifted */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    double flux;
    chemistry_get_fluxes(p, i, &flux);

    /* Update the conserved variable */
    chd->metal_mass[i] += flux;

    /* Update the diffused metal mass */
    chd->diffused_metal_mass[i] += flux;
  }

  /* Reset the fluxes now that they have been applied */
  chemistry_reset_chemistry_fluxes(p);

  /* Invalidate the particle time-step. It is considered to be inactive until
     dt is set again in hydro_prepare_force() */
  chd->flux_dt = -1.0f;

  /* Element-wise sanity checks */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    const double m_metal_old = chd->metal_mass[i];
    chemistry_check_unphysical_state(&chd->metal_mass[i], m_metal_old,
                                     hydro_get_mass(p), /*callloc=*/10,
                                     /*element*/ i);
  }

  /* Sanity check on the total metal mass */
  chemistry_check_unphysical_total_metal_mass(p, 10);


}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_H */
