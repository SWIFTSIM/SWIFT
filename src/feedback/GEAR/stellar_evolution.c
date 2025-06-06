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
#include "exp10.h"
#include "hdf5_functions.h"
#include "initial_mass_function.h"
#include "lifetime.h"
#include "random.h"
#include "stellar_evolution_struct.h"
#include "supernovae_ia.h"
#include "supernovae_ii.h"

#include <math.h>
#include <stddef.h>

#define DEFAULT_STAR_MINIMAL_GRAVITY_MASS_MSUN 1e-1

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
  return number_supernovae_i + (rand_sn < frac_sn);
}

/**
 * @brief Update the #spart mass from the supernovae ejected mass.
 *
 * This function deals with each star_type.
 *
 * Note: This function is called by
 * stellar_evolution_compute_discrete_feedback_properties() and
 * stellar_evolution_compute_continuous_feedback_properties().
 *
 * @param sp The particle to act upon
 * @param sm The #stellar_model structure.
 */
void stellar_evolution_sn_apply_ejected_mass(struct spart* restrict sp,
                                             const struct stellar_model* sm) {
  /* If a star is a discrete star */
  if (sp->star_type == single_star) {
    const int null_mass = (sp->mass == sp->feedback_data.mass_ejected);
    const int negative_mass = (sp->mass < sp->feedback_data.mass_ejected);

    if (null_mass) {
      message("Star %lld (m_star = %e, m_ej = %e) completely exploded!", sp->id,
              sp->mass, sp->feedback_data.mass_ejected);
      /* If the star ejects all its mass (for very massive stars), give it a
         zero mass so that we know it has exploded.
         We do not remove the star from the simulation to keep track of its
         properties, e.g. to check the IMF sampling (with sinks).

         Bug fix (28.04.2024): The mass of the star should not be set to
         0 because of gravity. So, we give some minimal value. */
      sp->mass = sm->discrete_star_minimal_gravity_mass;

      /* If somehow the star has a negative mass, we have a problem. */
    } else if (negative_mass) {
      error(
          "(Discrete star) Negative mass (m_star = %e, m_ej = %e), skipping "
          "current star: %lli",
          sp->mass, sp->feedback_data.mass_ejected, sp->id);
      /* Reset everything */
      sp->feedback_data.number_snia = 0;
      sp->feedback_data.number_snii = 0;
      sp->feedback_data.mass_ejected = 0;

      /* Reset energy to avoid injecting anything in the
         runner_iact_nonsym_feedback_apply() */
      sp->feedback_data.energy_ejected = 0;
      return;
    } else {
      /* Update the mass */
      sp->mass -= sp->feedback_data.mass_ejected;
    }

    /* If the star is the continuous part of the IMF or the enteire IMF */
  } else {
    /* Check if we can ejected the required amount of elements. */
    const int negative_mass = (sp->mass <= sp->feedback_data.mass_ejected);
    if (negative_mass) {
      warning(
          "(Continuous star) Negative mass (m_star = %e, m_ej = %e), skipping "
          "current star: %lli",
          sp->mass, sp->feedback_data.mass_ejected, sp->id);
      /* Reset everything */
      sp->feedback_data.number_snia = 0;
      sp->feedback_data.number_snii = 0;
      sp->feedback_data.mass_ejected = 0;

      /* Reset energy to avoid injecting anything in the
         runner_iact_nonsym_feedback_apply() */
      sp->feedback_data.energy_ejected = 0;
      return;
    }
    /* Update the mass */
    sp->mass -= sp->feedback_data.mass_ejected;
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
 * @param m_init Birth mass of the stellar particle (solMass).
 * @param number_snia_f (Floating) Number of SNIa produced by the stellar
 * particle.
 * @param number_snii_f (Floating) Number of SNII produced by the stellar
 * particle.
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
  const float mass_frac_snii =
      supernovae_ii_get_ejected_mass_fraction_processed_from_integral(
          &sm->snii, log_m_end_step, log_m_beg_step);

  /* Sum the contributions from SNIa and SNII */
  sp->feedback_data.mass_ejected = mass_frac_snii * sp->sf_data.birth_mass +
                                   mass_snia * phys_const->const_solar_mass;

  /* Removes the ejected mass from the star */
  stellar_evolution_sn_apply_ejected_mass(sp, sm);

  /* Now deal with the metals */

  /* Get the SNIa yields */
  const float* snia_yields = supernovae_ia_get_yields(&sm->snia);

  /* Compute the SNII yields */
  float snii_yields[GEAR_CHEMISTRY_ELEMENT_COUNT];
  supernovae_ii_get_yields_from_integral(&sm->snii, log_m_end_step,
                                         log_m_beg_step, snii_yields);

  /* Compute the mass fraction of non processed elements */
  const float non_processed =
      supernovae_ii_get_ejected_mass_fraction_non_processed_from_integral(
          &sm->snii, log_m_end_step, log_m_beg_step);

  /* Set the yields */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    /* Compute the mass fraction of metals */
    sp->feedback_data.metal_mass_ejected[i] =
        /* Supernovae II yields */
        snii_yields[i] +
        /* Gas contained in stars initial metallicity */
        chemistry_get_star_metal_mass_fraction_for_feedback(sp)[i] *
            non_processed;

    /* Convert it to total mass */
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
 * @param m_beg_step Mass of a star ending its life at the begining of the step
 * (solMass)
 * @param m_end_step Mass of a star ending its life at the end of the step
 * (solMass)
 * @param m_init Birth mass in solMass.
 * @param number_snia Number of SNIa produced by the stellar particle.
 * @param number_snii Number of SNII produced by the stellar particle.
 *
 */
void stellar_evolution_compute_discrete_feedback_properties(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const, const float m_beg_step,
    const float m_end_step, const float m_init, const int number_snia,
    const int number_snii) {

  /* Limit the mass within the imf limits */
  const float m_beg_lim = min(m_beg_step, sm->imf.mass_max);
  const float m_end_lim = max(m_end_step, sm->imf.mass_min);

  /* Compute the average mass */
  const float m_avg = 0.5 * (m_beg_lim + m_end_lim);
  const float log_m_avg = log10(m_avg);

  /* Compute the mass ejected */
  /* SNIa */
  const float mass_snia =
      supernovae_ia_get_ejected_mass_processed(&sm->snia) * number_snia;

  /* SNII */
  const float mass_snii =
      supernovae_ii_get_ejected_mass_fraction_processed_from_raw(&sm->snii,
                                                                 log_m_avg) *
      m_avg * number_snii;

  sp->feedback_data.mass_ejected = mass_snia + mass_snii;

  /* Transform into internal units */
  sp->feedback_data.mass_ejected *= phys_const->const_solar_mass;

  /* Removes the ejected mass from the star */
  stellar_evolution_sn_apply_ejected_mass(sp, sm);

  /* Get the SNIa yields */
  const float* snia_yields = supernovae_ia_get_yields(&sm->snia);

  /* Compute the SNII yields */
  float snii_yields[GEAR_CHEMISTRY_ELEMENT_COUNT];
  supernovae_ii_get_yields_from_raw(&sm->snii, log_m_avg, snii_yields);

  /* Compute the mass fraction of non processed elements */
  const float non_processed =
      supernovae_ii_get_ejected_mass_fraction_non_processed_from_raw(&sm->snii,
                                                                     log_m_avg);

  /* Set the yields */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {

    /* Compute the mass fraction of metals */
    sp->feedback_data.metal_mass_ejected[i] =
        /* Supernovae II yields */
        snii_yields[i] +
        /* Gas contained in stars initial metallicity */
        chemistry_get_star_metal_mass_fraction_for_feedback(sp)[i] *
            non_processed;

    /* Convert it to total mass */
    sp->feedback_data.metal_mass_ejected[i] *= m_avg * number_snii;

    /* Supernovae Ia yields */
    sp->feedback_data.metal_mass_ejected[i] += snia_yields[i] * number_snia;

    /* Convert everything in code units */
    sp->feedback_data.metal_mass_ejected[i] *= phys_const->const_solar_mass;
  }
}

/**
 * @brief Evolve an individual star represented by a #spart.
 *
 * This function compute the SN rate and yields before sending
 * this information to a different MPI rank.
 * It also compute the supernovae energy to be released by the
 * star.
 *
 * Here I am using Myr-solar mass units internally in order to
 * avoid numerical errors.
 *
 * Note: This function treats the case of single/individual stars.
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
void stellar_evolution_evolve_individual_star(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct phys_const* phys_const, const integertime_t ti_begin,
    const double star_age_beg_step, const double dt) {

  /* Check that this function is called for single_star only. */
  if (sp->star_type != single_star) {
    error("This function can only be called for single/individual star!");
  }

  /* Convert the inputs */
  const double conversion_to_myr = phys_const->const_year * 1e6;
  const double star_age_end_step_myr =
      (star_age_beg_step + dt) / conversion_to_myr;
  const double star_age_beg_step_myr = star_age_beg_step / conversion_to_myr;

  /* Get the metallicity */
  const float metallicity =
      chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);

  const float log_mass = log10(sp->mass / phys_const->const_solar_mass);
  const float lifetime_myr = pow(10, lifetime_get_log_lifetime_from_mass(
                                         &sm->lifetime, log_mass, metallicity));

  /* if the lifetime is outside the interval */
  if ((lifetime_myr < star_age_beg_step_myr) ||
      (lifetime_myr > star_age_end_step_myr))
    return;

  message(
      "(%lld) lifetime_myr=%g %g star_age_beg_step=%g star_age_end_step=%g "
      "(%g)",
      sp->id, lifetime_myr, lifetime_myr * conversion_to_myr,
      star_age_beg_step_myr, star_age_end_step_myr,
      sp->mass / phys_const->const_solar_mass);

  /* This is needed by stellar_evolution_compute_discrete_feedback_properties(),
     but this is not used inside the function. */
  const float m_init = 0;

  /* Get the integer number of supernovae */
  const int number_snia = 0;
  const int number_snii = 1;

  /* Save the number of supernovae */
  sp->feedback_data.number_snia = 0;
  sp->feedback_data.number_snii = number_snii;

  /* this is needed for  stellar_evolution_compute_discrete_feedback_properties
   */
  const float m_beg_step = sp->mass / phys_const->const_solar_mass;
  const float m_end_step = sp->mass / phys_const->const_solar_mass;
  const float m_avg = 0.5 * (m_beg_step + m_end_step);

  /* Compute the yields */
  stellar_evolution_compute_discrete_feedback_properties(
      sp, sm, phys_const, m_beg_step, m_end_step, m_init, number_snia,
      number_snii);

  /* Compute the supernovae energy associated to the stellar particle */

  const float energy_conversion =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) / 1e51;

  /* initialize */
  sp->feedback_data.energy_ejected = 0;

  /* snia contribution */
  const float snia_energy = sm->snia.energy_per_supernovae;
  sp->feedback_data.energy_ejected +=
      sp->feedback_data.number_snia * snia_energy;

  /* snii contribution */
  const float snii_energy =
      supernovae_ii_get_energy_from_progenitor_mass(&sm->snii, m_avg) /
      energy_conversion;
  sp->feedback_data.energy_ejected +=
      sp->feedback_data.number_snii * snii_energy;
}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function compute the SN rate and yields before sending
 * this information to a different MPI rank.
 * It also compute the supernovae energy to be released by the
 * star.
 *
 * Here I am using Myr-solar mass units internally in order to
 * avoid numerical errors.
 *
 * Note: This function treats the case of particles representing the whole IMF
 * (star_type = star_population) and the particles representing only the
 * continuous part of the IMF (star_type = star_population_continuous_IMF).
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

  /* Check that this function is called for populations of stars and not
     individual stars. */
  if (sp->star_type == single_star) {
    error(
        "This function can only be called for sparts representing stars "
        "populations!");
  }

  /* Convert the inputs */
  const double conversion_to_myr = phys_const->const_year * 1e6;
  const double star_age_beg_step_myr = star_age_beg_step / conversion_to_myr;
  const double dt_myr = dt / conversion_to_myr;

  /* Get the metallicity */
  const float metallicity =
      chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);

  /* Compute masses range */
  const float log_m_beg_step =
      star_age_beg_step == 0.
          ? FLT_MAX
          : lifetime_get_log_mass_from_lifetime(
                &sm->lifetime, log10(star_age_beg_step_myr), metallicity);
  const float log_m_end_step = lifetime_get_log_mass_from_lifetime(
      &sm->lifetime, log10(star_age_beg_step_myr + dt_myr), metallicity);

  float m_beg_step = star_age_beg_step == 0. ? FLT_MAX : exp10(log_m_beg_step);
  float m_end_step = exp10(log_m_end_step);

  /* Limit the mass interval to the IMF boundaries */
  m_end_step = max(m_end_step, sm->imf.mass_min);
  m_beg_step = min(m_beg_step, sm->imf.mass_max);

  /* Here we are outside the IMF, i.e., both masses are too large or too small
   */
  if (m_end_step >= m_beg_step) return;

  /* Star particles representing only the continuous part of the IMF need a
  special treatment. They do not contain stars above the mass that separate the
  IMF into two parts (variable called minimal_discrete_mass_Msun in the sink
  module). So, if m_end_step > minimal_discrete_mass_Msun, you don't do
  feedback. Note that the sm structure contains different information for the
  'first stars' and the 'late stars'. The right sm data is passed to this
  function so we do not need any special treatment here. */
  if (sp->star_type == star_population_continuous_IMF) {
    /* If it's not time yet for feedback, exit. Notice that both masses are in
      solar mass. */
    if (m_end_step > sm->imf.minimal_discrete_mass_Msun) {
      return;
    }

    /* If we are in a case where
                m_beg_step > minimal_discrete_mass_Msun > m_end_step,
       then we need to be careful. We don't want feedback from the discrete
       part, only the continuous part. Hence, we need to update m_beg_step.
    */
    if (m_beg_step > sm->imf.minimal_discrete_mass_Msun) {
      m_beg_step = sm->imf.minimal_discrete_mass_Msun;
    }
  }

  /* Check if the star can produce a supernovae */
  const int can_produce_snia =
      supernovae_ia_can_explode(&sm->snia, m_end_step, m_beg_step);
  const int can_produce_snii =
      supernovae_ii_can_explode(&sm->snii, m_end_step, m_beg_step);

  /* Is it possible to generate a supernovae? */
  if (!can_produce_snia && !can_produce_snii) return;

  /* Compute the initial mass. The initial mass is different if the star
     particle is of type 'star_population' or
     'star_population_continuous_IMF'. The function call treats both cases. */
  const float m_init =
      stellar_evolution_compute_initial_mass(sp, sm, phys_const);

  /* Then, for 'star_population_continuous_IMF', everything remain the same as
     with the "old" 'star_population'! */

  /* Compute number of SNIa */
  float number_snia_f = 0;
  if (can_produce_snia) {
    number_snia_f = supernovae_ia_get_number_per_unit_mass(
                        &sm->snia, m_end_step, m_beg_step) *
                    m_init;
  }

  /* Compute number of SNII */
  float number_snii_f = 0;
  if (can_produce_snii) {
    number_snii_f = supernovae_ii_get_number_per_unit_mass(
                        &sm->snii, m_end_step, m_beg_step) *
                    m_init;
  }

  /* Does this star produce a supernovae? */
  if (number_snia_f == 0 && number_snii_f == 0) return;

  /* Compute the properties of the feedback (e.g. yields) */
  if (sm->discrete_yields) {
    /* Get the integer number of supernovae */
    const int number_snia = stellar_evolution_compute_integer_number_supernovae(
        sp, number_snia_f, ti_begin, random_number_stellar_feedback_1);

    /* Get the integer number of supernovae */
    const int number_snii = stellar_evolution_compute_integer_number_supernovae(
        sp, number_snii_f, ti_begin, random_number_stellar_feedback_2);

    /* Do we have a supernovae? */
    if (number_snia == 0 && number_snii == 0) return;

    /* Save the number of supernovae */
    sp->feedback_data.number_snia = number_snia;
    sp->feedback_data.number_snii = number_snii;

    /* Compute the yields */
    stellar_evolution_compute_discrete_feedback_properties(
        sp, sm, phys_const, m_beg_step, m_end_step, m_init, number_snia,
        number_snii);

  } else {
    /* Save the number of supernovae */
    sp->feedback_data.number_snia = number_snia_f;
    sp->feedback_data.number_snii = number_snii_f;

    /* Compute the yields */
    stellar_evolution_compute_continuous_feedback_properties(
        sp, sm, phys_const, log_m_beg_step, log_m_end_step, m_beg_step,
        m_end_step, m_init, number_snia_f, number_snii_f);
  }

  /* Compute the supernovae energy associated to the stellar particle */

  /* Compute the average mass (in solar mass) */
  const float m_avg = 0.5 * (m_beg_step + m_end_step);
  const float energy_conversion =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) / 1e51;

  /* initialize */
  sp->feedback_data.energy_ejected = 0;

  /* snia contribution */
  const float snia_energy = sm->snia.energy_per_supernovae;
  sp->feedback_data.energy_ejected +=
      sp->feedback_data.number_snia * snia_energy;

  /* snii contribution */
  const float snii_energy =
      supernovae_ii_get_energy_from_progenitor_mass(&sm->snii, m_avg) /
      energy_conversion;
  sp->feedback_data.energy_ejected +=
      sp->feedback_data.number_snii * snii_energy;
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
 * @brief Get the index of the element .
 *
 * @param sm The #stellar_model.
 * @param element_name The element name.
 */
int stellar_evolution_get_element_index(const struct stellar_model* sm,
                                        const char* element_name) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (strcmp(stellar_evolution_get_element_name(sm, i), element_name) == 0)
      return i;
  }
  error("Chemical element %s not found !", element_name);

  return -1;
}

/**
 * @brief Get the solar abundance of the element .
 *
 * @param sm The #stellar_model.
 * @param element_name The element name.
 */
float stellar_evolution_get_solar_abundance(const struct stellar_model* sm,
                                            const char* element_name) {

  int element_index = stellar_evolution_get_element_index(sm, element_name);
  float solar_abundance = sm->solar_abundances[element_index];

  return solar_abundance;
}

/**
 * @brief Read the name of all the elements present in the tables.
 *
 * @param sm The #stellar_model.
 * @param params The #swift_params.
 */
void stellar_evolution_read_elements(struct stellar_model* sm,
                                     struct swift_params* params) {

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
    strcpy(sm->elements_name + i * GEAR_LABELS_SIZE, elements[i]);
  }

  /* Cleanup. */
  parser_free_param_string_array(nval, elements);

  /* Add the metals to the end. */
  strcpy(
      sm->elements_name + (GEAR_CHEMISTRY_ELEMENT_COUNT - 1) * GEAR_LABELS_SIZE,
      "Metals");

  /* Check the elements */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    for (int j = i + 1; j < GEAR_CHEMISTRY_ELEMENT_COUNT; j++) {
      const char* el_i = stellar_evolution_get_element_name(sm, i);
      const char* el_j = stellar_evolution_get_element_name(sm, j);
      if (strcmp(el_i, el_j) == 0) {
        error("You need to provide each element only once (%s).", el_i);
      }
    }
  }
}

/**
 * @brief Read the solar abundances.
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
void stellar_evolution_read_solar_abundances(struct stellar_model* sm,
                                             struct swift_params* params) {

#if defined(HAVE_HDF5)

  /* Get the yields table */
  char filename[DESCRIPTION_BUFFER_SIZE];
  parser_get_param_string(params, "GEARFeedback:yields_table", filename);

  /* Open file. */
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s.\n", filename);

  /* Open group. */
  hid_t group_id = H5Gopen(file_id, "Data", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group Data.\n");

  /* Read the data */
  io_read_array_attribute(group_id, "SolarMassAbundances", FLOAT,
                          sm->solar_abundances, GEAR_CHEMISTRY_ELEMENT_COUNT);

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
 * @brief Initialize the global properties of the stellar evolution scheme.
 *
 * @param sm The #stellar_model.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
void stellar_evolution_props_init(struct stellar_model* sm,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* us,
                                  struct swift_params* params,
                                  const struct cosmology* cosmo) {

  /* Read the list of elements */
  stellar_evolution_read_elements(sm, params);

  /* Read the solar abundances */
  stellar_evolution_read_solar_abundances(sm, params);

  /* Use the discrete yields approach? */
  sm->discrete_yields =
      parser_get_param_int(params, "GEARFeedback:discrete_yields");

  /* Initialize the initial mass function */
  initial_mass_function_init(&sm->imf, phys_const, us, params,
                             sm->yields_table);

  /* Initialize the lifetime model */
  lifetime_init(&sm->lifetime, phys_const, us, params, sm->yields_table);

  /* Initialize the supernovae Ia model */
  supernovae_ia_init(&sm->snia, phys_const, us, params, sm);

  /* Initialize the supernovae II model */
  supernovae_ii_init(&sm->snii, params, sm, us);

  /* Initialize the minimal gravity mass for the stars */
  /* const float default_star_minimal_gravity_mass_Msun = 1e-1; */
  sm->discrete_star_minimal_gravity_mass = parser_get_opt_param_float(
      params, "GEARFeedback:discrete_star_minimal_gravity_mass_Msun",
      DEFAULT_STAR_MINIMAL_GRAVITY_MASS_MSUN);

  /* Convert from M_sun to internal units */
  sm->discrete_star_minimal_gravity_mass *= phys_const->const_solar_mass;

  if (engine_rank == 0) {
    message("discrete_star_minimal_gravity_mass: (internal units)          %e",
            sm->discrete_star_minimal_gravity_mass);
  }
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

/**
 * @brief Computes the initial mass of a #spart. This function distinguishes
 * between the stellar particle representing a whole IMF and the stellar
 * particles representing only the continuous part.
 *
 * @param sp The particle for which we compute the initial mass.
 * @param sm The #stellar_model structure.
 * @param phys_const the physical constants in internal units.
 * @param (return) m_init Initial mass of the star particle (in M_sun).
 */
float stellar_evolution_compute_initial_mass(
    const struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const) {

  const struct initial_mass_function* imf = &sm->imf;
  switch (sp->star_type) {
    case star_population:
      return sp->sf_data.birth_mass / phys_const->const_solar_mass;
    case star_population_continuous_IMF: {
      double M_IMF_tot, M_d_dummy, M_c_dummy;
      initial_mass_function_compute_Mc_Md_Mtot(imf, &M_c_dummy, &M_d_dummy,
                                               &M_IMF_tot);
      /* No need to convert from internal units to M_sun because the masses are
         already in solar masses (to avoid numerical errors) */
      return M_IMF_tot;
    }
    case single_star:
      return sp->sf_data.birth_mass / phys_const->const_solar_mass;
    default: {
      error("This star_type (%d) is not implemented!", sp->star_type);
      return -1.0;
    }
  }
}
