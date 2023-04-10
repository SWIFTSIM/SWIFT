/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_AGORA_H
#define SWIFT_CHEMISTRY_AGORA_H

/**
 * @file src/chemistry/none/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cosmology.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Get the name of the element i.
 *
 * @param sm The #stellar_model.
 * @param i The element indice.
 */
static INLINE const char* chemistry_get_element_name(
    const struct chemistry_global_data* data, int i) {

  return data->elements_name + i * AGORA_LABELS_SIZE;
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

  float mass = hydro_get_mass(p);

  /* Store the chemistry struct in the star particle */
  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    sp->chemistry_data.metal_mass_fraction[i] =
        p->chemistry_data.smoothed_metal_mass_fraction[i];

    /* Remove the metals taken by the star. */
    p->chemistry_data.metal_mass[i] *= mass / (mass + sp->mass);
  }
}

/**
 * @brief Initialises the chemistry properties.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The global chemistry information (to be filled).
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* read parameters */
  const float initial_metallicity = parser_get_param_float(
      parameter_file, "AGORAChemistry:initial_metallicity");

  if (initial_metallicity < 0) {
    message("Setting the initial metallicity from the snapshot.");
  } else {
    message("Setting the initial metallicity from the parameter file.");
  }

  /* Set the initial metallicities */
  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    data->initial_metallicities[i] = initial_metallicity;
  }

  /* Check if need to scale the initial metallicity */
  const int scale_metallicity = parser_get_opt_param_int(
      parameter_file, "AGORAChemistry:scale_initial_metallicity", 0);

  /* Scale the metallicities if required */
  if (scale_metallicity) {

    /* Set element name (Metals goes to the end). */
    strcpy(data->elements_name + (0) * AGORA_LABELS_SIZE, "Fe");
    strcpy(data->elements_name +
               (AGORA_CHEMISTRY_ELEMENT_COUNT - 1) * AGORA_LABELS_SIZE,
           "Metals");

    /* Read the solar abundances. */
    data->solar_abundances[0] = parser_get_opt_param_double(
        parameter_file, "AGORAChemistry:solar_abundance_Fe", 0.001771);
    data->solar_abundances[AGORA_CHEMISTRY_ELEMENT_COUNT - 1] =
        parser_get_opt_param_double(
            parameter_file, "AGORAChemistry:solar_abundance_Metals", 0.02);

    /* Scale the metallicities */
    for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
      data->initial_metallicities[i] *= data->solar_abundances[i];
    }
  }
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

  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Reset the smoothed metallicity */
    cpd->smoothed_metal_mass_fraction[i] = 0.f;
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

  message("Chemistry function is 'AGORA'.");
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param cd The global chemistry information.
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

  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] += cpd->metal_mass[i] * kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    cpd->smoothed_metal_mass_fraction[i] *= factor;
  }
}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 * @param with_cosmology Are we running with the cosmology?
 * @param time Current time of the simulation.
 * @param dt Time step (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* restrict p, const struct cosmology* cosmo,
    const int with_cosmology, const double time, const double dt) {}

/**
 * @brief Computes the chemistry-related time-step constraint.
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
  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.smoothed_metal_mass_fraction[i] =
        p->chemistry_data.metal_mass[i] / hydro_get_mass(p);
  }
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information used for this run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
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
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data* data, struct spart* restrict sp) {

  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
    sp->chemistry_data.metal_mass_fraction[i] = data->initial_metallicities[i];
  }
}

/**
 * @brief Initialise the chemistry properties of a black hole with
 * the chemistry properties of the gas it is born from.
 *
 * Nothing to do here.
 *
 * @param bp_data The black hole data to initialise.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_bpart_from_part(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {
  error("To be implemented.");
}

/**
 * @brief Add the chemistry data of a gas particle to a black hole.
 *
 * Nothing to do here.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {
  error("To be implemented.");
}

/**
 * @brief Transfer chemistry data of a gas particle to a black hole.
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
  error("To be implemented.");
}

/**
 * @brief Add the chemistry data of a sink particle to a sink.
 *
 * Nothing to do here.
 *
 * @param si_data The black hole data to add to.
 * @param sj_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_sink_to_sink(
    struct chemistry_sink_data* si_data,
    const struct chemistry_sink_data* sj_data) {
  error("To be implemented.");
}

/**
 * @brief Add the chemistry data of a gas particle to a sink.
 *
 * Nothing to do here.
 *
 * @param sp_data The sink data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_sink(
    struct chemistry_sink_data* sp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {
  error("To be implemented.");
}

/**
 * @brief Split the metal content of a particle into n pieces
 *
 * Nothing to do here.
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void chemistry_split_part(
    struct part* p, const double n) {
  error("To be implemented.");
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct part* restrict p) {
  error("To be implemented.");
  return 0.f;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * No metallicity treatment here -> return NULL array.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(const struct part* restrict p) {
  error("To be implemented.");
  return NULL;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data
      .metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * No metallicity treatment here -> return NULL array.
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
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part* restrict p) {

  return p->chemistry_data
      .smoothed_metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
 *
 * No metallicity treatment here -> return NULL array.
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
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data
      .smoothed_metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
 *
 * No metallicity treatment here -> return NULL array.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.smoothed_metal_mass_fraction;
}

/**
 * @brief Returns the total metal mass of the
 * gas particle to be used in the stats related routines.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part* restrict p) {

  return p->chemistry_data.metal_mass[AGORA_CHEMISTRY_ELEMENT_COUNT - 1];
}

/**
 * @brief Returns the total metal mass of the
 * star particle to be used in the stats related routines.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart* restrict sp) {

  return sp->chemistry_data
             .metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT - 1] *
         sp->mass;
}

/**
 * @brief Returns the total metal mass of the
 * black hole particle to be used in the stats related routines.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart* restrict bp) {

  error("To be implemented.");
  return 0.f;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in the luminosity calculations.
 *
 * No metallicity treatment here -> return 0.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_luminosity(
    const struct spart* restrict sp) {

  error("To be implemented.");
  return 0.f;
}

#endif /* SWIFT_CHEMISTRY_AGORA_H */
