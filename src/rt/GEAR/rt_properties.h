/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_PROPERTIES_GEAR_H
#define SWIFT_RT_PROPERTIES_GEAR_H

#include "rt_parameters.h"

/**
 * @file src/rt/GEAR/rt_properties.h
 * @brief Main header file for the 'GEAR' radiative transfer scheme
 * properties.
 */

/**
 * @brief Properties of the 'GEAR' radiative transfer model
 */
struct rt_props {

  /* Are we running with hydro or star controlled injection?
   * This is added to avoid #ifdef macros as far as possible */
  int hydro_controlled_injection;

  /* Do we need to run a conversion after the zeroth
   * step, but before the first step? */
  int convert_stars_after_zeroth_step;
  int convert_parts_after_zeroth_step;

  /* Are we using constant stellar emission rates? */
  int use_const_emission_rates;

  /* Frequency bin edges for photon groups
   * Includes 0 as leftmost edge, doesn't include infinity as
   * rightmost bin edge*/
  float photon_groups[RT_NGROUPS];

  /* Global constant stellar emission rates */
  double stellar_const_emission_rates[RT_NGROUPS];

  /* CFL condition */
  float CFL_condition;

  /* do we set initial ionization mass fractions manually? */
  int set_initial_ionization_mass_fractions;
  int set_equilibrium_initial_ionization_mass_fractions;

  /* initial mass fractions for ionization species */
  /* the following are required for manually setting exact values */
  float mass_fraction_HI_init;
  float mass_fraction_HII_init;
  float mass_fraction_HeI_init;
  float mass_fraction_HeII_init;
  float mass_fraction_HeIII_init;
  float number_density_electrons_init; /* todo: do we need this? */

  /* Hydrogen and Helium mass fractions of the non-metal portion of the gas */
  float hydrogen_mass_fraction;
  float helium_mass_fraction;

  /* Skip thermochemistry? For testing/debugging only! */
  int skip_thermochemistry;

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Do extended tests where we assume that all parts
   * have spart neighbours? */
  /* skip this for GEAR */
  /* int debug_do_all_parts_have_stars_checks; */

  /* radiation emitted by stars this step. This is not really a property,
   * but a placeholder to sum up a global variable. It's being reset
   * every timestep. */
  int debug_radiation_emitted_this_step;

  /* total radiation emitted by stars. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_emitted_tot;

  /* radiation absorbed by gas this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  int debug_radiation_absorbed_this_step;

  /* total radiation absorbed by gas. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_absorbed_tot;

  /* Interactions of a star with gas during injection prep this step. This is
   * not really a property, but a placeholder to sum up a global variable */
  int debug_star_injection_prep_iacts_with_parts_this_step;

  /* Interactions of a star with gas during injection prep. This is not
   * really a property, but a placeholder to sum up a global variable */
  unsigned long long debug_star_injection_prep_iacts_with_parts_tot;

  /* Interactions of a star with gas during injection prep this step. This is
   * not really a property, but a placeholder to sum up a global variable */
  int debug_part_injection_prep_iacts_with_stars_this_step;

  /* Interactions of a star with gas during injection prep. This is not
   * really a property, but a placeholder to sum up a global variable */
  unsigned long long debug_part_injection_prep_iacts_with_stars_tot;

  /* Total radiation energy in the gas. It's being reset every step. */
  float debug_total_radiation_conserved_energy[RT_NGROUPS];
  float debug_total_radiation_energy_density[RT_NGROUPS];
  float debug_total_star_emitted_energy[RT_NGROUPS];

  /* Files to write energy budget to after every step */
  FILE* conserved_energy_filep;
  FILE* energy_density_filep;
  FILE* star_emitted_energy_filep;
#endif
};

/**
 * @brief Print the RT model.
 *
 * @param rtp The #rt_props
 */
__attribute__((always_inline)) INLINE static void rt_props_print(
    const struct rt_props* rtp) {

  /* Only the master print */
  if (engine_rank != 0) return;

  message("Radiative transfer scheme: '%s'", RT_IMPLEMENTATION);
  char messagestring[200] = "Using photon frequency bins: [0.";
  char freqstring[20];
  for (int g = 1; g < RT_NGROUPS; g++) {
    sprintf(freqstring, ", %.3g", rtp->photon_groups[g]);
    strcat(messagestring, freqstring);
  }
  strcat(messagestring, "]");
  message("%s", messagestring);

  if (rtp->use_const_emission_rates) {
    strcpy(messagestring, "Using constant stellar emission rates: [ ");
    for (int g = 0; g < RT_NGROUPS; g++) {
      sprintf(freqstring, "%.3g ", rtp->stellar_const_emission_rates[g]);
      strcat(messagestring, freqstring);
    }
    strcat(messagestring, "]");
    message("%s", messagestring);
  }

  if (rtp->set_equilibrium_initial_ionization_mass_fractions)
    message(
        "Setting initial mass fractions for "
        "ionizing species assuming ionization equilibrium");
  if (rtp->set_initial_ionization_mass_fractions)
    message(
        "Using initial ionization mass fractions specified in parameter file");
  if (rtp->skip_thermochemistry)
    message("WARNING: Thermochemistry will be skipped.");
}

/**
 * @brief Initialize the global properties of the RT scheme.
 *
 * @param rtp The #rt_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void rt_props_init(
    struct rt_props* rtp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    struct cosmology* cosmo) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  rtp->hydro_controlled_injection = 1;
#else
  rtp->hydro_controlled_injection = 0;
#endif

  /* Make sure we reset debugging counters correctly after
   * zeroth step. */
#ifdef SWIFT_RT_DEBUG_CHECKS
  rtp->convert_parts_after_zeroth_step = 1;
  rtp->convert_stars_after_zeroth_step = 1;
#else
  rtp->convert_parts_after_zeroth_step = 0;
  rtp->convert_stars_after_zeroth_step = rtp->hydro_controlled_injection;
#endif

  /* Read in photon frequency group properties */
  /* ----------------------------------------- */
  if (RT_NGROUPS <= 0) {
    error(
        "You need to run GEAR-RT with at least 1 photon group, "
        "you have %d",
        RT_NGROUPS);
  } else if (RT_NGROUPS == 1) {
    rtp->photon_groups[0] = 0.f;
  } else {
    float frequencies[RT_NGROUPS - 1];
    parser_get_param_float_array(params, "GEARRT:photon_groups_Hz",
                                 RT_NGROUPS - 1, frequencies);
    float Hz_internal = units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);
    float Hz_internal_inv = 1.f / Hz_internal;
    for (int g = 0; g < RT_NGROUPS - 1; g++) {
      rtp->photon_groups[g + 1] = frequencies[g] * Hz_internal_inv;
    }
    rtp->photon_groups[0] = 0.f;
  }

  /* Are we using constant emission rates? */
  /* ------------------------------------- */
  rtp->use_const_emission_rates = parser_get_opt_param_int(
      params, "GEARRT:use_const_emission_rates", /* default = */ 0);

  if (rtp->use_const_emission_rates) {
    double emission_rates[RT_NGROUPS];
    parser_get_param_double_array(params, "GEARRT:star_emission_rates_LSol",
                                  RT_NGROUPS, emission_rates);
    const double unit_power = units_cgs_conversion_factor(us, UNIT_CONV_POWER);
    const double unit_power_inv = 1. / unit_power;
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->stellar_const_emission_rates[g] = emission_rates[g] * unit_power_inv;
    }
  } else {
    /* kill the run for now */
    error("GEAR-RT can't run without constant stellar emission rates for now.");
  }

  /* get reduced speed of light factor */
  /* --------------------------------- */
  const float f_r = parser_get_param_float(params, "GEARRT:f_reduce_c");
  rt_params.reduced_speed_of_light = phys_const->const_speed_light_c * f_r;
  rt_params.reduced_speed_of_light_inverse =
      1.f / rt_params.reduced_speed_of_light;

  /* get CFL condition */
  /* ----------------- */
  const float CFL = parser_get_param_float(params, "GEARRT:CFL_condition");
  rtp->CFL_condition = CFL;

  /* Get thermochemistry set-up */
  /* -------------------------- */
  rtp->hydrogen_mass_fraction =
      parser_get_param_float(params, "GEARRT:hydrogen_mass_fraction");
  rtp->helium_mass_fraction = 1.f - rtp->hydrogen_mass_fraction;

  /* Are we manually overwriting initial mass fractions of H and He? */
  rtp->set_initial_ionization_mass_fractions = parser_get_opt_param_int(
      params, "GEARRT:set_initial_ionization_mass_fractions",
      /* default = */ 0);
  if (rtp->set_initial_ionization_mass_fractions) {
    /* Read in mass fractions */
    rtp->mass_fraction_HI_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HI");
    rtp->mass_fraction_HII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HII");
    rtp->mass_fraction_HeI_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeI");
    rtp->mass_fraction_HeII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeII");
    rtp->mass_fraction_HeIII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeIII");

    /* Temporary check neglecting metals. Make sure we sum up to 1. */
    const float h_sum =
        rtp->mass_fraction_HI_init + rtp->mass_fraction_HII_init;
    if (fabsf(h_sum - rtp->hydrogen_mass_fraction) > 1e-4)
      error(
          "Inconsistent Hydrogen mass fractions: XH_tot %.6g != XH %.6g + XH+ "
          "%.6g",
          rtp->hydrogen_mass_fraction, rtp->mass_fraction_HI_init,
          rtp->mass_fraction_HII_init);

    const float he_sum = rtp->mass_fraction_HeI_init +
                         rtp->mass_fraction_HeII_init +
                         rtp->mass_fraction_HeIII_init;
    if (fabsf(he_sum - rtp->helium_mass_fraction) > 1e-4)
      error(
          "Inconsistent Helium mass fractions: XHe_tot %.6g != XHe %.6g + XH+ "
          "%.6g + XHe++ %.6g",
          rtp->helium_mass_fraction, rtp->mass_fraction_HeI_init,
          rtp->mass_fraction_HeII_init, rtp->mass_fraction_HeIII_init);

    const float mass_fraction_sum = h_sum + he_sum;
    if (fabsf(mass_fraction_sum - 1.f) > 1e-5)
      error("Ionizing species mass fraction sums up to %.6f, I expect 1.0",
            mass_fraction_sum);
  } else {
    /* Initialize properties to deliberately bogus values */
    rtp->mass_fraction_HI_init = -1.f;
    rtp->mass_fraction_HII_init = -1.f;
    rtp->mass_fraction_HeI_init = -1.f;
    rtp->mass_fraction_HeII_init = -1.f;
    rtp->mass_fraction_HeIII_init = -1.f;
  }

  /* Are we setting up initial mass fractions in equilibrium? */
  rtp->set_equilibrium_initial_ionization_mass_fractions =
      parser_get_opt_param_int(
          params, "GEARRT:set_equilibrium_initial_ionization_mass_fractions",
          /* default = */ 0);

  /* Mark that we need some conversion after the first step now */
  if (rtp->set_equilibrium_initial_ionization_mass_fractions ||
      rtp->set_initial_ionization_mass_fractions)
    rtp->convert_parts_after_zeroth_step = 1;

  if (rtp->set_equilibrium_initial_ionization_mass_fractions &&
      rtp->set_initial_ionization_mass_fractions)
    error(
        "Can't use equilibrium initial ionization mass fractions "
        "simultaneously with manually set mass fractions. Pick one.");

  /* Are we skipping thermochemistry? */
  rtp->skip_thermochemistry = parser_get_opt_param_int(
      params, "GEARRT:skip_thermochemistry", /* default = */ 0);

#ifdef SWIFT_RT_DEBUG_CHECKS
  rtp->debug_radiation_emitted_tot = 0ULL;
  rtp->debug_radiation_absorbed_tot = 0ULL;
  rtp->debug_star_injection_prep_iacts_with_parts_tot = 0LL;
  rtp->debug_part_injection_prep_iacts_with_stars_tot = 0LL;
  for (int g = 0; g < RT_NGROUPS; g++)
    rtp->debug_total_star_emitted_energy[g] = 0.f;

  /* Open up files for energy budgets */
  rtp->conserved_energy_filep = fopen("RT_conserved_energy_budget.txt", "w");
  if (rtp->conserved_energy_filep == NULL)
    error("Couldn't open RT conserved energy budget file to write in");
  rtp->energy_density_filep = fopen("RT_energy_density_budget.txt", "w");
  if (rtp->energy_density_filep == NULL)
    error("Couldn't open RT energy density budget file to write in");
  rtp->star_emitted_energy_filep = fopen("RT_star_injected_energy.txt", "w");
  if (rtp->star_emitted_energy_filep == NULL)
    error("Couldn't open RT star energy budget file to write in");

  if (rtp->use_const_emission_rates) {
    FILE* files[3] = {rtp->conserved_energy_filep, rtp->energy_density_filep,
                      rtp->star_emitted_energy_filep};
    for (int f = 0; f < 3; f++) {
      fprintf(files[f], "# Emission rates: ");
      const double solar_luminosity = 3.826e33; /* erg/s */
      for (int g = 0; g < RT_NGROUPS; g++)
        fprintf(files[f], "%12.6e ",
                rtp->stellar_const_emission_rates[g] * solar_luminosity);
      fprintf(files[f], "\n");
    }
  }
#endif

  /* Finishers */
  /* --------- */

  /* After initialisation, print params to screen */
  rt_props_print(rtp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Radiative transfer initialized");
  }
}

/**
 * @brief Write an RT properties struct to the given FILE as a
 * stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_dump(
    const struct rt_props* props, FILE* stream) {

  restart_write_blocks((void*)props, sizeof(struct rt_props), 1, stream,
                       "RT props", "RT properties struct");
}

/**
 * @brief Restore an RT properties struct from the given FILE as
 * a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
}

#endif /* SWIFT_RT_PROPERTIES_GEAR_H */
