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
#ifndef SWIFT_CHEMISTRY_GEAR_H
#define SWIFT_CHEMISTRY_GEAR_H

/**
 * @file src/chemistry/GEAR/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Compute the metal mass fraction
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static float
chemistry_metal_mass_fraction(const struct part* restrict p,
                              const struct xpart* restrict xp) {
  return p->chemistry_data.Z;
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
  data->initial_metallicity = parser_get_opt_param_float(
      parameter_file, "GearChemistry:InitialMetallicity", -1);
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

  for (int i = 0; i < chemistry_element_count; i++) {
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
  const float m = p->mass;

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < chemistry_element_count; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] +=
        m * cpd->metal_mass_fraction[i] * kernel_root;

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
    struct part* restrict p, const struct cosmology* cosmo) {}

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
  error("Needs implementing!");
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
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  p->chemistry_data.Z = data->initial_metallicity;
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

  error(
      "MATTHIEU: Loic this is a new function. I don't know whether you"
      " want something here.");
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
  error("Loic: to be implemented");
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
  error("Loic: to be implemented");
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

#endif /* SWIFT_CHEMISTRY_GEAR_H */
