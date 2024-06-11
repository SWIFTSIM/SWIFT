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
#ifndef SWIFT_CHEMISTRY_GEAR_H
#define SWIFT_CHEMISTRY_GEAR_H

/**
 * @file src/chemistry/GEAR/chemistry.h
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

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
        p->chemistry_data.smoothed_metal_mass_fraction[k];

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
        p->chemistry_data.smoothed_metal_mass_fraction[i];
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

  message("Chemistry function is 'Gear'.");
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
  const float factor = pow_dimension(h_inv) / p->rho; /* 1 / h^d * rho */

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] += cpd->metal_mass[i] * kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    cpd->smoothed_metal_mass_fraction[i] *= factor;
  }
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
  return FLT_MAX;
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

  chemistry_init_part(p, data);
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
    sp->chemistry_data.metal_mass_fraction[i] = data->initial_metallicities[i];
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

  /* gas mass */
  const float mass = hydro_get_mass(p);

  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    double mk = s->chemistry_data.metal_mass_fraction[k] * ms_old +
                p->chemistry_data.smoothed_metal_mass_fraction[k] * mass;

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
      .smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
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
      .smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
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

#endif /* SWIFT_CHEMISTRY_GEAR_H */
