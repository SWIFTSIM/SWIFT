/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "stellar_evolution.h"

/* Include local headers */
#include "hdf5_functions.h"
#include "initial_mass_function.h"
#include "lifetime.h"
#include "random.h"
#include "stellar_evolution_struct.h"
#include "supernovae_ia.h"
#include "supernovae_ii.h"

#include <math.h>
#include <stddef.h>

/**
 * @brief Print the stellar model.
 *
 * @param sm The #stellar_model.
 */
void stellar_model_print(const struct stellar_model* sm) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* Print the type of yields */
  message("Discrete yields? %i", sm->discrete_yields);

  /* Print the sub properties */
  initial_mass_function_print(&sm->imf);
  lifetime_print(&sm->lifetime);
  supernovae_ia_print(&sm->snia);
  supernovae_ii_print(&sm->snii);
}

/**
 * @brief Compute the integer number of supernovae from the floating number.
 *
 * @param sp The particle to act upon
 * @param number_supernovae_f Floating number of supernovae during this step.
 * @param ti_begin The #integertime_t at the begining of the step.
 * @param random_type The categorie of random.
 *
 * @return The integer number of supernovae.
 */
int stellar_evolution_compute_integer_number_supernovae(
    struct spart* restrict sp, float number_supernovae_f,
    const integertime_t ti_begin, enum random_number_type random_type) {

  const int number_supernovae_i = floor(number_supernovae_f);

  /* Get the random number for the decimal part */
  const float rand_sn = random_unit_interval(sp->id, ti_begin, random_type);

  /* Get the fraction part */
  const float frac_sn = number_supernovae_f - number_supernovae_i;

  /* Get the integer number of SN */
  return number_supernovae_i + ((rand_sn < frac_sn) ? 1 : 0);
}

/**
 * @brief Compute the feedback properties.
 *
 * @param sp The particle to act upon
 * @param sm The #stellar_model structure.
 * @param phys_const The physical constants in the internal unit system.
 * @param log_m_beg_step Mass of a star ending its life at the begining of the
 * step (log10(solMass))
 * @param log_m_end_step Mass of a star ending its life at the end of the step
 * (log10(solMass))
 * @param m_beg_step Mass of a star ending its life at the begining of the step
 * (solMass)
 * @param m_end_step Mass of a star ending its life at the end of the step
 * (solMass)
 * @param number_snia_f (Floating) Number of SNIa produced by the stellar particle.
 * @param number_snii_f (Floating) Number of SNII produced by the stellar particle.
 *
 */
void stellar_evolution_compute_continuous_feedback_properties(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const, const float log_m_beg_step,
    const float log_m_end_step, const float m_beg_step, const float m_end_step,
    const float m_init, const float number_snia_f, const float number_snii_f) {

  /* Compute the mass ejected */
  /* SNIa */
  const float mass_snia =
    supernovae_ia_get_ejected_mass_processed(&sm->snia) * number_snia_f;

  /* SNII */
  const float mass_frac_snii = supernovae_ii_get_ejected_mass_fraction_processed(
          &sm->snii, log_m_end_step, log_m_beg_step); 

  /* Sum the contributions from SNIa and SNII */
  sp->feedback_data.mass_ejected = mass_frac_snii * sp->birth.mass
    + mass_snia * phys_const->const_solar_mass;

  if (sp->mass <= sp->feedback_data.mass_ejected) {
    error("Stars cannot have negative mass. (%g <= %g). Initial mass = %g",
          sp->mass, sp->feedback_data.mass_ejected, sp->sf_data.birth_mass);
  }

  /* Update the mass */
  sp->mass -= sp->feedback_data.mass_ejected;

  /* Now deal with the metals */

  /* Get the SNIa yields */
  const float* snia_yields = supernovae_ia_get_yields(&sm->snia);

  /* Compute the SNII yields */
  float snii_yields[GEAR_CHEMISTRY_ELEMENT_COUNT];
  supernovae_ii_get_yields(&sm->snii, log_m_end_step, log_m_beg_step,
                           snii_yields);

  /* Compute the mass fraction of non processed elements */
  const float non_processed = supernovae_ii_get_ejected_mass_fraction(
      &sm->snii, log_m_end_step, log_m_beg_step);

  /* Set the yields */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Compute the mass fraction of metals */
    sp->feedback_data.metal_mass_ejected[i] =
      /* Supernovae II yields */
      snii_yields[i] +
      /* Gas contained in stars initial metallicity */
      chemistry_get_metal_mass_fraction_for_feedback(sp)[i] * non_processed;

    /* Convert it to total mass */
    sp->feedback_data.metal_mass_ejected[i] *= sp->sf_data.birth_mass;

    sp->feedback_data.metal_mass_ejected[i] *= sp->sf_data.birth_mass;

    /* Add the Supernovae Ia */
    sp->feedback_data.metal_mass_ejected[i] +=
      snia_yields[i] * number_snia_f * phys_const->const_solar_mass;

  }
}

/**
 * @brief Compute the feedback properties.
 *
 * @param sp The particle to act upon
 * @param sm The #stellar_model structure.
 * @param phys_const The physical constants in the internal unit system.
 * @param log_m_beg_step Mass of a star ending its life at the begining of the
 * step (log10(solMass))
 * @param log_m_end_step Mass of a star ending its life at the end of the step
 * (log10(solMass))
 * @param m_beg_step Mass of a star ending its life at the begining of the step
 * (solMass)
 * @param m_end_step Mass of a star ending its life at the end of the step
 * (solMass)
 * @param number_snia Number of SNIa produced by the stellar particle.
 * @param number_snii Number of SNII produced by the stellar particle.
 *
 */
void stellar_evolution_compute_discrete_feedback_properties(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const, const float log_m_beg_step,
    const float log_m_end_step, const float m_beg_step, const float m_end_step,
    const float m_init, const int number_snia, const int number_snii) {

  /* Get the normalization to the average */
  const float normalization =
    number_snii == 0 ? 0. :
    number_snii / (supernovae_ii_get_number_per_unit_mass(&sm->snii, m_end_step, m_beg_step) *
		   m_init);

  /* Compute the mass ejected */
  /* SNIa */
  const float mass_snia =
    (number_snia == 0) ? 0 :
    (supernovae_ia_get_ejected_mass_processed(&sm->snia) * number_snia);

  /* SNII */
  const float mass_snii =
      normalization * supernovae_ii_get_ejected_mass_fraction_processed(
                          &sm->snii, log_m_end_step, log_m_beg_step);

  sp->feedback_data.mass_ejected = mass_snia + mass_snii;

  /* Transform into internal units */
  sp->feedback_data.mass_ejected *= phys_const->const_solar_mass;

  if (sp->mass <= sp->feedback_data.mass_ejected) {
    error("Stars cannot have negative mass. (%g <= %g). Initial mass = %g",
          sp->mass, sp->feedback_data.mass_ejected, sp->sf_data.birth_mass);
  }

  /* Update the mass */
  sp->mass -= sp->feedback_data.mass_ejected;

  /* Get the SNIa yields */
  const float* snia_yields = supernovae_ia_get_yields(&sm->snia);

  /* Compute the SNII yields (without the normalization) */
  float snii_yields[GEAR_CHEMISTRY_ELEMENT_COUNT];
  supernovae_ii_get_yields(&sm->snii, log_m_end_step, log_m_beg_step,
                           snii_yields);

  /* Compute the mass fraction of non processed elements */
  const float non_processed =
      normalization * supernovae_ii_get_ejected_mass_fraction(
                          &sm->snii, log_m_end_step, log_m_beg_step);

  /* Set the yields */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Compute the mass fraction of metals */
    sp->feedback_data.metal_mass_ejected[i] =
        /* Supernovae Ia yields */
        snia_yields[i] * number_snia +
        /* Supernovae II yields */
        normalization * snii_yields[i] +
        /* Gas contained in stars initial metallicity */
        chemistry_get_metal_mass_fraction_for_feedback(sp)[i] * non_processed;

    /* Convert it to total mass */
    sp->feedback_data.metal_mass_ejected[i] *= sp->sf_data.birth_mass;
  }
}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function compute the SN rate and yields before sending
 * this information to a different MPI rank.
 *
 * Here I am using Myr-solar mass units internally in order to
 * avoid numerical errors.
 *
 * @param sp The particle to act upon
 * @param sm The #stellar_model structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The physical constants in the internal unit system.
 * @param ti_begin The #integertime_t at the begining of the step.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
void stellar_evolution_evolve_spart(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct phys_const* phys_const, const integertime_t ti_begin,
    const double star_age_beg_step, const double dt) {

  /* Convert the inputs */
  const double conversion_to_myr = phys_const->const_year * 1e6;
  const double star_age_beg_step_myr = star_age_beg_step / conversion_to_myr;
  const double dt_myr = dt / conversion_to_myr;

  /* Get the metallicity */
  const float metallicity =
      chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Compute masses range */
  const float log_m_beg_step =
      star_age_beg_step == 0.
          ? FLT_MAX
          : lifetime_get_log_mass_from_lifetime(
                &sm->lifetime, log10(star_age_beg_step_myr), metallicity);
  const float log_m_end_step = lifetime_get_log_mass_from_lifetime(
      &sm->lifetime, log10(star_age_beg_step_myr + dt_myr), metallicity);

  const float m_beg_step =
      star_age_beg_step == 0. ? FLT_MAX : pow(10, log_m_beg_step);
  const float m_end_step = pow(10, log_m_end_step);

  /* Check if the star can produce a supernovae */
  const int can_produce_snia =
      supernovae_ia_can_explode(&sm->snia, m_end_step, m_beg_step);
  const int can_produce_snii =
      supernovae_ii_can_explode(&sm->snii, m_end_step, m_beg_step);

  /* Is it possible to generate a supernovae? */
  if (!can_produce_snia && !can_produce_snii) return;

  /* Compute the initial mass */
  const float m_init = sp->sf_data.birth_mass / phys_const->const_solar_mass;

  /* Compute number of SNIa */
  const float number_snia_f = can_produce_snia ?
    supernovae_ia_get_number_per_unit_mass(&sm->snia, m_end_step, m_beg_step) * m_init : 0;

  /* Compute number of SNII */
  const float number_snii_f = can_produce_snii ?
    supernovae_ii_get_number_per_unit_mass(&sm->snii, m_end_step, m_beg_step) * m_init : 0;

  /* Does this star produce a supernovae? */
  if (number_snia_f == 0 && number_snii_f == 0) return;

  sp->feedback_data.number_sn = number_snia_f + number_snii_f;

  /* Compute the properties of the feedback (e.g. yields) */
  if (sm->discrete_yields) {
    /* Get the integer number of supernovae */
    const int number_snia = stellar_evolution_compute_integer_number_supernovae(
      sp, number_snia_f, ti_begin, random_number_stellar_feedback_1);

    /* Get the integer number of supernovae */
    const int number_snii = stellar_evolution_compute_integer_number_supernovae(
      sp, number_snii_f, ti_begin, random_number_stellar_feedback_2);

    stellar_evolution_compute_discrete_feedback_properties(
        sp, sm, phys_const, log_m_beg_step, log_m_end_step, m_beg_step,
        m_end_step, m_init, number_snia, number_snii);

  } else {
    stellar_evolution_compute_continuous_feedback_properties(
        sp, sm, phys_const, log_m_beg_step, log_m_end_step, m_beg_step,
        m_end_step, m_init, number_snia_f, number_snii_f);
  }
}

/**
 * @brief Get the name of the element i.
 *
 * @param sm The #stellar_model.
 * @param i The element indice.
 */
const char* stellar_evolution_get_element_name(const struct stellar_model* sm,
                                               int i) {

  return sm->elements_name + i * GEAR_LABELS_SIZE;
}

/**
 * @brief Read the name of all the elements present in the tables.
 *
 * @param sm The #stellar_model.
 */
void stellar_evolution_read_elements(struct stellar_model* sm,
                                     struct swift_params* params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data", &file_id, &group_id);

  /* Read the elements */
  io_read_string_array_attribute(group_id, "elts", sm->elements_name,
                                 GEAR_CHEMISTRY_ELEMENT_COUNT,
                                 GEAR_LABELS_SIZE);

  /* Check that we received correctly the metals */
  if (strcmp(stellar_evolution_get_element_name(
                 sm, GEAR_CHEMISTRY_ELEMENT_COUNT - 1),
             "Metals") != 0) {
    error(
        "The chemistry table should contain the metals in the last column "
        "(found %s)",
        stellar_evolution_get_element_name(sm,
                                           GEAR_CHEMISTRY_ELEMENT_COUNT - 1));
  }

  /* Print the name of the elements */
  char txt[GEAR_CHEMISTRY_ELEMENT_COUNT * (GEAR_LABELS_SIZE + 2)] = "";
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (i != 0) {
      strcat(txt, ", ");
    }
    strcat(txt, stellar_evolution_get_element_name(sm, i));
  }

  if (engine_rank == 0) {
    message("Chemistry elements: %s", txt);
  }

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Initialize the global properties of the stellar evolution scheme.
 *
 * @param sm The #stellar_model.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void stellar_evolution_props_init(struct stellar_model* sm,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* us,
                                  struct swift_params* params,
                                  const struct cosmology* cosmo) {

  /* Read the list of elements */
  stellar_evolution_read_elements(sm, params);

  /* Use the discrete yields approach? */
  sm->discrete_yields =
      parser_get_param_int(params, "GEARFeedback:discrete_yields");

  /* Initialize the initial mass function */
  initial_mass_function_init(&sm->imf, phys_const, us, params);

  /* Initialize the lifetime model */
  lifetime_init(&sm->lifetime, phys_const, us, params);

  /* Initialize the supernovae Ia model */
  supernovae_ia_init(&sm->snia, phys_const, us, params, sm);

  /* Initialize the supernovae II model */
  supernovae_ii_init(&sm->snii, phys_const, us, params, sm);
}

/**
 * @brief Write a stellar_evolution struct to the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything has been copied in the
 * feedback.
 *
 * @param sm the struct
 * @param stream the file stream
 */
void stellar_evolution_dump(const struct stellar_model* sm, FILE* stream) {

  /* Dump the initial mass function */
  initial_mass_function_dump(&sm->imf, stream, sm);

  /* Dump the lifetime model */
  lifetime_dump(&sm->lifetime, stream, sm);

  /* Dump the supernovae Ia model */
  supernovae_ia_dump(&sm->snia, stream, sm);

  /* Dump the supernovae II model */
  supernovae_ii_dump(&sm->snii, stream, sm);
}

/**
 * @brief Restore a stellar_evolution struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything has been copied in the
 * feedback.
 *
 * @param sm the struct
 * @param stream the file stream
 */
void stellar_evolution_restore(struct stellar_model* sm, FILE* stream) {

  /* Restore the initial mass function */
  initial_mass_function_restore(&sm->imf, stream, sm);

  /* Restore the lifetime model */
  lifetime_restore(&sm->lifetime, stream, sm);

  /* Restore the supernovae Ia model */
  supernovae_ia_restore(&sm->snia, stream, sm);

  /* Restore the supernovae II model */
  supernovae_ii_restore(&sm->snii, stream, sm);
}

/**
 * @brief Clean the allocated memory.
 *
 * @param sm the #stellar_model.
 */
void stellar_evolution_clean(struct stellar_model* sm) {

  initial_mass_function_clean(&sm->imf);
  lifetime_clean(&sm->lifetime);
  supernovae_ia_clean(&sm->snia);
  supernovae_ii_clean(&sm->snii);
}
