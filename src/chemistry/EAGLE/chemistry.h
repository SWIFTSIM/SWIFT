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
#ifndef SWIFT_CHEMISTRY_EAGLE_H
#define SWIFT_CHEMISTRY_EAGLE_H

/**
 * @file src/chemistry/EAGLE/chemistry.h
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
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen",
      "Neon",     "Magnesium", "Silicon", "Iron"};

  return chemistry_element_names[elem];
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

  cpd->smoothed_metal_mass_fraction_total = 0.f;
  cpd->smoothed_iron_mass_fraction_from_SNIa = 0.f;
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
  const float m = hydro_get_mass(p);

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int i = 0; i < chemistry_element_count; i++) {
    /* Final operation on the density (add self-contribution). */
    cpd->smoothed_metal_mass_fraction[i] +=
        m * cpd->metal_mass_fraction[i] * kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    cpd->smoothed_metal_mass_fraction[i] *= factor;
  }

  /* Smooth mass fraction of all metals */
  cpd->smoothed_metal_mass_fraction_total +=
      m * cpd->metal_mass_fraction_total * kernel_root;
  cpd->smoothed_metal_mass_fraction_total *= factor;

  /* Smooth iron mass fraction from SNIa */
  cpd->smoothed_iron_mass_fraction_from_SNIa +=
      m * cpd->iron_mass_fraction_from_SNIa * kernel_root;
  cpd->smoothed_iron_mass_fraction_from_SNIa *= factor;
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

  /* Just make all the smoothed fields default to the un-smoothed values */
  struct chemistry_part_data* cpd = &p->chemistry_data;

  /* Total metal mass fraction */
  cpd->smoothed_metal_mass_fraction_total = cpd->metal_mass_fraction_total;

  /* Iron frac from SNIa */
  cpd->smoothed_iron_mass_fraction_from_SNIa =
      cpd->iron_mass_fraction_from_SNIa;

  /* Individual metal mass fractions */
  for (int i = 0; i < chemistry_element_count; i++) {
    cpd->smoothed_metal_mass_fraction[i] = cpd->metal_mass_fraction[i];
  }
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
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

  // Add initialization of all other fields in chemistry_part_data struct.
  if (data->initial_metal_mass_fraction_total != -1) {
    p->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem)
      p->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
  }
  chemistry_init_part(p, data);
}

/**
 * @brief Initialises the chemistry properties.
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

  /* Read the total metallicity */
  data->initial_metal_mass_fraction_total = parser_get_opt_param_float(
      parameter_file, "EAGLEChemistry:init_abundance_metal", -1);

  if (data->initial_metal_mass_fraction_total != -1) {
    /* Read the individual mass fractions */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      char buffer[50];
      sprintf(buffer, "EAGLEChemistry:init_abundance_%s",
              chemistry_get_element_name((enum chemistry_element)elem));

      data->initial_metal_mass_fraction[elem] =
          parser_get_param_float(parameter_file, buffer);
    }
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

  message("Chemistry model is 'EAGLE' tracking %d elements.",
          chemistry_element_count);
}

#endif /* SWIFT_CHEMISTRY_EAGLE_H */
