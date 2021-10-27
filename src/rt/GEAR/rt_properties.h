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

  if (RT_NGROUPS <= 0) {
    error(
        "You need to run GEAR-RT with at least 1 photon group, "
        "you have %d",
        RT_NGROUPS);
  } else if (RT_NGROUPS == 1) {
    rtp->photon_groups[0] = 0.f;
  } else {
    /* Read in parameters */
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
  rtp->use_const_emission_rates = parser_get_opt_param_float(
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
  const float f_r = parser_get_param_float(params, "GEARRT:f_reduce_c");
  rt_params.reduced_speed_of_light = phys_const->const_speed_light_c * f_r;
  rt_params.reduced_speed_of_light_inverse =
      1.f / rt_params.reduced_speed_of_light;

  /* get CFL condition */
  const float CFL = parser_get_param_float(params, "GEARRT:CFL_condition");
  rtp->CFL_condition = CFL;

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
