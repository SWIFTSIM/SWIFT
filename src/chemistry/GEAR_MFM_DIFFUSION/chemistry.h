/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Roduit Darwin (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_H

/**
 * @file src/chemistry/GEAR_MFM_DIFFUSION/chemistry.h
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "chemistry_gradients.h"
#include "chemistry_flux.h"
#include "error.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief TODO
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static double
chemistry_part_compute_diffusion_coefficient(struct part* restrict p) {
  return kernel_gamma2*p->h*p->h;
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
    sp->chemistry_data.metal_mass_fraction[k] = p->chemistry_data.metal_mass[k] / mass;

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
        p->chemistry_data.metal_mass[i] / p->mass;
  }
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {

  message("Chemistry function is 'Gear MFM diffusion'.");
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
}

/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Reset the smoothed metallicity */
    cpd->smoothed_metal_mass_fraction[i] = 0.f;
  }

  /* Reset the geometry matrix */
  cpd->geometry.volume = 0.0f;
  cpd->geometry.matrix_E[0][0] = 0.0f;
  cpd->geometry.matrix_E[0][1] = 0.0f;
  cpd->geometry.matrix_E[0][2] = 0.0f;
  cpd->geometry.matrix_E[1][0] = 0.0f;
  cpd->geometry.matrix_E[1][1] = 0.0f;
  cpd->geometry.matrix_E[1][2] = 0.0f;
  cpd->geometry.matrix_E[2][0] = 0.0f;
  cpd->geometry.matrix_E[2][1] = 0.0f;
  cpd->geometry.matrix_E[2][2] = 0.0f;

  /* Update the diffusion coefficient for the new loops */
  /* TODO: This might be moved to the end of density or gradient, if we want to
     use a shear tensor using grad velocity */
  p->chemistry_data.kappa = chemistry_part_compute_diffusion_coefficient(p);

  /* Init the gradient */
  chemistry_gradients_init(p);
}

/**
 * @brief Finishes the smooth metal calculation.
 *
 * Multiplies the smoothed metallicity and number of neighbours by the
 * appropiate constants and add the self-contribution term.
 *
 * This function requires the #hydro_end_density to have been called.
 *
 * @param p The particle to act upon.
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float ihdim = pow_dimension(h_inv);
  const float factor = pow_dimension(h_inv) / p->rho; /* 1 / h^d * rho */

  struct chemistry_part_data* cpd = &p->chemistry_data;

  /* First finish the smoothed metallicites */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] += cpd->metal_mass[i] * kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    cpd->smoothed_metal_mass_fraction[i] *= factor;
  }

  /* Final operation on the geometry. */
  /* We multiply with the smoothing kernel normalization ih3 and calculate the
   * volume */
  const float volume_inv = ihdim * (p->chemistry_data.geometry.volume + kernel_root);
  const float volume = 1.0f / volume_inv;
  p->chemistry_data.geometry.volume = volume;

  /* we multiply with the smoothing kernel normalization */
  p->chemistry_data.geometry.matrix_E[0][0] *= ihdim;
  p->chemistry_data.geometry.matrix_E[0][1] *= ihdim;
  p->chemistry_data.geometry.matrix_E[0][2] *= ihdim;
  p->chemistry_data.geometry.matrix_E[1][0] *= ihdim;
  p->chemistry_data.geometry.matrix_E[1][1] *= ihdim;
  p->chemistry_data.geometry.matrix_E[1][2] *= ihdim;
  p->chemistry_data.geometry.matrix_E[2][0] *= ihdim;
  p->chemistry_data.geometry.matrix_E[2][1] *= ihdim;
  p->chemistry_data.geometry.matrix_E[2][2] *= ihdim;

  /* Check the condition number to see if we have a stable geometry. */
  const float condition_number_E =
      p->chemistry_data.geometry.matrix_E[0][0] * p->chemistry_data.geometry.matrix_E[0][0] +
      p->chemistry_data.geometry.matrix_E[0][1] * p->chemistry_data.geometry.matrix_E[0][1] +
      p->chemistry_data.geometry.matrix_E[0][2] * p->chemistry_data.geometry.matrix_E[0][2] +
      p->chemistry_data.geometry.matrix_E[1][0] * p->chemistry_data.geometry.matrix_E[1][0] +
      p->chemistry_data.geometry.matrix_E[1][1] * p->chemistry_data.geometry.matrix_E[1][1] +
      p->chemistry_data.geometry.matrix_E[1][2] * p->chemistry_data.geometry.matrix_E[1][2] +
      p->chemistry_data.geometry.matrix_E[2][0] * p->chemistry_data.geometry.matrix_E[2][0] +
      p->chemistry_data.geometry.matrix_E[2][1] * p->chemistry_data.geometry.matrix_E[2][1] +
      p->chemistry_data.geometry.matrix_E[2][2] * p->chemistry_data.geometry.matrix_E[2][2];

  float condition_number = 0.0f;
  if (invert_dimension_by_dimension_matrix(p->chemistry_data.geometry.matrix_E) != 0) {
    /* something went wrong in the inversion; force bad condition number */
    condition_number = const_gizmo_max_condition_number + 1.0f;
  } else {
    const float condition_number_Einv =
        p->chemistry_data.geometry.matrix_E[0][0] * p->chemistry_data.geometry.matrix_E[0][0] +
        p->chemistry_data.geometry.matrix_E[0][1] * p->chemistry_data.geometry.matrix_E[0][1] +
        p->chemistry_data.geometry.matrix_E[0][2] * p->chemistry_data.geometry.matrix_E[0][2] +
        p->chemistry_data.geometry.matrix_E[1][0] * p->chemistry_data.geometry.matrix_E[1][0] +
        p->chemistry_data.geometry.matrix_E[1][1] * p->chemistry_data.geometry.matrix_E[1][1] +
        p->chemistry_data.geometry.matrix_E[1][2] * p->chemistry_data.geometry.matrix_E[1][2] +
        p->chemistry_data.geometry.matrix_E[2][0] * p->chemistry_data.geometry.matrix_E[2][0] +
        p->chemistry_data.geometry.matrix_E[2][1] * p->chemistry_data.geometry.matrix_E[2][1] +
        p->chemistry_data.geometry.matrix_E[2][2] * p->chemistry_data.geometry.matrix_E[2][2];

    condition_number =
        hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);
  }

  if (condition_number > const_gizmo_max_condition_number &&
      p->chemistry_data.geometry.wcorr > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
          const_gizmo_max_condition_number, condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            condition_number, const_gizmo_max_condition_number, p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    p->chemistry_data.geometry.wcorr = const_gizmo_w_correction_factor * p->chemistry_data.geometry.wcorr;
  }
 
  /* Check that the metal masses are physical */
  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {
    double n_metal_old = p->chemistry_data.metal_mass[g];

#ifdef SWIFT_DEBUG_CHECKS
    if (volume == 0.) {
      error("Volume is 0!");
    }
#endif

    /* Sanity checks */
    chemistry_check_unphysical_state(&p->chemistry_data.metal_mass[g], n_metal_old, /*callloc=*/0);
  }
}

/**
 * @brief Finishes the gradient calculation.
 *
 * TODO: Add it at the right place in the code.
 *
 * Just a wrapper around chemistry_gradients_finalize, which can be an empty method,
 * in which case no gradients are used.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_gradient(
    struct part* p) {
  chemistry_gradients_finalise(p);

  /* Do not reset the gradients here. We need nabla_otimes_q to compute the
     chemistry timestep */
}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* restrict p, const struct cosmology* cosmo,
    const int with_cosmology, const double time, const double dt) {}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void chemistry_prepare_force(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo, const float dt_alpha,
    const float dt_therm) {
  p->chemistry_data.flux_dt = dt_therm;
}


/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
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

  /* Set the smoothed fractions with the non smoothed fractions */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.smoothed_metal_mass_fraction[i] =
        p->chemistry_data.metal_mass[i] / hydro_get_mass(p);
  }

  /* Update geometry data */
  p->chemistry_data.geometry.volume = 1.0f;
  p->chemistry_data.geometry.matrix_E[0][0] = 1.0f;
  p->chemistry_data.geometry.matrix_E[0][1] = 0.0f;
  p->chemistry_data.geometry.matrix_E[0][2] = 0.0f;
  p->chemistry_data.geometry.matrix_E[1][0] = 0.0f;
  p->chemistry_data.geometry.matrix_E[1][1] = 1.0f;
  p->chemistry_data.geometry.matrix_E[1][2] = 0.0f;
  p->chemistry_data.geometry.matrix_E[2][0] = 0.0f;
  p->chemistry_data.geometry.matrix_E[2][1] = 0.0f;
  p->chemistry_data.geometry.matrix_E[2][2] = 1.0f;
}

/**
 * @brief Computes the chemistry-related time-step constraint.
 *
 * No constraints in the GEAR model (no diffusion) --> FLT_MAX
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, const struct part* restrict p) {

  /* K = kappa * I_3. The norm is the Froebenius norm */
  /* const float norm_K = sqrtf(3*p->chemistry_data.kappa*p->chemistry_data.kappa); */
  /* const float norm_q = p->chemistry_data.kappa; */
  return kernel_gamma2 * p->h*p->h / (sqrtf(3) * fabs(p->chemistry_data.kappa));
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param cosmo The #cosmology.
 * @param data The global chemistry information.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (data->initial_metallicities[i] < 0) {
      /* Use the value from the IC. We are reading the metal mass fraction. */
      p->chemistry_data.metal_mass[i] *= hydro_get_mass(p);
    } else {
      /* Use the value from the parameter file */
      p->chemistry_data.metal_mass[i] =
          data->initial_metallicities[i] * hydro_get_mass(p);
    }
  }
  p->chemistry_data.kappa = chemistry_part_compute_diffusion_coefficient(p);
  chemistry_init_part(p, data);

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  p->chemistry_data.geometry.wcorr = 1.0f;
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
 * @brief Add the chemistry data of a sink particle to a sink.
 *
 * @param si_data The black hole data to add to.
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
 * @param sp_data The sink data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
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
  error("To be implemented.");
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
  error("Loic: to be implemented");
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
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * This is unused in GEAR. --> return NULL
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(const struct part* restrict p) {
  error("Not implemented");
  return NULL;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * This is unused in GEAR. --> return 0
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct part* restrict p) {
  error("Not implemented");
  return 0.f;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_star_total_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data
      .metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the total iron mass fraction of the
 * star particle to be used in feedback/enrichment related routines.
 * We assume iron to be stored at index 0.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_star_total_iron_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction[0];
}

/**
 * @brief Returns the total iron mass fraction of the
 * sink particle to be used in feedback/enrichment related routines.
 * We assume iron to be stored at index 0.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_sink_total_iron_mass_fraction_for_feedback(
    const struct sink* restrict sink) {

  return sink->chemistry_data.metal_mass_fraction[0];
}

/**
 * @brief Returns the abundances (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part* restrict p) {

  return p->chemistry_data
      .metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT - 1]/p->mass;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
 *
 * @TODO: This must be changed
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_metal_mass_fraction_for_cooling(const struct part* restrict p) {

  return p->chemistry_data.smoothed_metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data
      .metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT - 1]/p->mass;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
 *
 * TODO: Take care of this...
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.smoothed_metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in the stats related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part* restrict p) {

  return p->chemistry_data.metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart* restrict sp) {

  return sp->chemistry_data
             .metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1] *
         sp->mass;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * black hole particle to be used in the stats related routines.
 *
 * @param bp Pointer to the BH particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart* restrict bp) {
  error("Not implemented");
  return 0.f;
}


/**
 * @brief Extra operations done during the kick. This needs to be
 * done before the particle mass is updated in the hydro_kick_extra.
 *
 * @TODO
 *
 * @param p Particle to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 */
__attribute__((always_inline)) INLINE static void chemistry_kick_extra(
    struct part* p, float dt_therm, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {

  if (p->chemistry_data.flux_dt > 0.0f) {

    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
      double flux;
      chemistry_part_get_fluxes(p, i, &flux);

      /* Update the conserved variable */
      p->chemistry_data.metal_mass[i] += flux;
    }

    /* Reset the fluxes, so that they do not get used again in kick1 */
    chemistry_part_reset_chemistry_fluxes(p);

    /* Invalidate the particle time-step. It is considered to be inactive until
       dt is set again in hydro_prepare_force() */
    p->chemistry_data.flux_dt = -1.0f;
  } else if (p->chemistry_data.flux_dt == 0.0f) {
    /* something tricky happens at the beginning of the simulation: the flux
       exchange is done for all particles, but using a time step of 0. This
       in itself is not a problem. However, it causes some issues with the
       initialisation of flux.dt for inactive particles, since this value will
       remain 0 until the particle is active again, and its flux.dt is set to
       the actual time step in hydro_prepare_force(). We have to make sure it
       is properly set to -1 here, so that inactive particles are indeed found
       to be inactive during the flux loop. */
    p->chemistry_data.flux_dt = -1.0f;
  }

  /* Sanity checks. We don't want negative metal masses. */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    const double n_metal_old = p->chemistry_data.metal_mass[i];
    chemistry_check_unphysical_state(&p->chemistry_data.metal_mass[i], n_metal_old, /*callloc=*/2);
  }

  /* Reset wcorr */
  p->chemistry_data.geometry.wcorr = 1.0f;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_H */