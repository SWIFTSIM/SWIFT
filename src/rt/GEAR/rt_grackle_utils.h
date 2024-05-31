/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GRACKLE_UTILS_H
#define SWIFT_RT_GRACKLE_UTILS_H

/* skip deprecation warnings. I cleaned old API calls. */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC

/* need hydro gamma */
#include "hydro.h"

#include <grackle.h>

/* need to rework (and check) code if changed */
#define FIELD_SIZE 1

/**
 * @file src/rt/GEAR/rt_grackle_utils.h
 * @brief Utility and helper functions related to using grackle.
 */

/**
 * @brief Update grackle units during run
 *
 * @param grackle_units grackle units struct
 * @param cosmo cosmology struct
 *
 * NOTE: In the current implementation, this function does nothing.
 * However, there might be use-cases in the future (e.g. switching
 * UV background on or off depending on redshift) that might be
 * needed in the future, which can be implemented into this function.
 */
__attribute__((always_inline)) INLINE void update_grackle_units_cosmo(
    code_units *grackle_units, const struct unit_system *us,
    const struct cosmology *restrict cosmo) {}

/**
 * @brief initialize grackle during rt_props_init
 *
 * @param grackle_units grackle units struct to fill up correctly.
 * @param grackle_chemistry_dat grackle chemistry data struct to fill up
 *correctly.
 * @param hydrogen_mass_fraction global hydrogen mass fraction.
 * @param grackle_verb run grackle in verbose mode?
 * @param case_B_recombination use grackle with case B recombination?
 * @param us #unit_system struct
 **/
__attribute__((always_inline)) INLINE static void rt_init_grackle(
    code_units *grackle_units, chemistry_data *grackle_chemistry_data,
    chemistry_data_storage *grackle_chemistry_rates,
    float hydrogen_mass_fraction, const int grackle_verb,
    const int case_B_recombination, const struct unit_system *us,
    const struct cosmology *restrict cosmo) {

  grackle_verbose = grackle_verb;

  /* Initialize units */
  /* ---------------- */
  /* we assume all quantities to be physical, not comoving */
  grackle_units->a_units = 1.0;
  grackle_units->a_value = 1.0;
  grackle_units->comoving_coordinates = 0;
  grackle_units->density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  grackle_units->length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  grackle_units->time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  /* Set velocity units */
  set_velocity_units(grackle_units);

  /* Chemistry Parameters */
  /* -------------------- */
  /* More details on
   * https://grackle.readthedocs.io/en/grackle-3.2.0/Integration.html#chemistry-data
   */

  if (local_initialize_chemistry_parameters(grackle_chemistry_data) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  /* chemistry on */
  grackle_chemistry_data->use_grackle = 1;
  /* cooling on */
  /* NOTE: without cooling on, it also won't heat... */
  grackle_chemistry_data->with_radiative_cooling = 1;
  /* 6 species atomic H and He */
  grackle_chemistry_data->primordial_chemistry = 1;
  /* No dust processes */
  grackle_chemistry_data->dust_chemistry = 0;
  /* No H2 formation on dust */
  grackle_chemistry_data->h2_on_dust = 0;
  /* metal cooling (uses Cloudy) off (for now) */
  grackle_chemistry_data->metal_cooling = 0;
  /* no cooling below CMB temperature */
  grackle_chemistry_data->cmb_temperature_floor = 1;
  /* UV background off */
  grackle_chemistry_data->UVbackground = 0;
  /* data file - currently not used */
  grackle_chemistry_data->grackle_data_file = "";
  /* adiabatic index */
  grackle_chemistry_data->Gamma = hydro_gamma;
  /* we'll provide grackle with ionization and heating rates from RT */
  grackle_chemistry_data->use_radiative_transfer = 1;

  /* fraction by mass of Hydrogen in the metal-free portion of the gas */
  grackle_chemistry_data->HydrogenFractionByMass = hydrogen_mass_fraction;
  /* Use case B recombination? (On-the-spot approximation) */
  grackle_chemistry_data->CaseBRecombination = case_B_recombination;

  if (local_initialize_chemistry_data(grackle_chemistry_data,
                                      grackle_chemistry_rates,
                                      grackle_units) == 0) {
    error("Error in initialize_chemistry_data");
  }
}

/**
 * @brief fill out a grackle field struct with the relevant (gas) data from a
 *particle
 *
 * @param grackle_fields (return) grackle field to copy into
 * @param density array of particle density
 * @param internal_energy array of particle internal_energy
 * @param species_densities array of species densities of particle (HI, HII,
 *HeI, HeII, HeIII, e-)
 * @param iact_rates array of interaction rates (heating, 3 ioniziation, H2
 *dissociation)
 *
 **/
__attribute__((always_inline)) INLINE static void
rt_get_grackle_particle_fields(grackle_field_data *grackle_fields,
                               gr_float density, gr_float internal_energy,
                               gr_float species_densities[6],
                               gr_float iact_rates[5]) {

  int *dimension = malloc(3 * sizeof(int));
  int *start = malloc(3 * sizeof(int));
  int *end = malloc(3 * sizeof(int));

  dimension[0] = FIELD_SIZE;
  dimension[1] = 0;
  dimension[2] = 0;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  end[0] = FIELD_SIZE - 1;
  end[1] = 0;
  end[2] = 0;

  grackle_fields->grid_dx = 0.;
  grackle_fields->grid_rank = 3;
  grackle_fields->grid_dimension = dimension;
  grackle_fields->grid_start = start;
  grackle_fields->grid_end = end;

  /* Set initial quantities */
  grackle_fields->density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->internal_energy = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->x_velocity = NULL;
  grackle_fields->y_velocity = NULL;
  grackle_fields->z_velocity = NULL;
  /* for primordial_chemistry >= 1 */
  grackle_fields->HI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeIII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->e_density = malloc(FIELD_SIZE * sizeof(gr_float));
  /* for primordial_chemistry >= 2 */
  grackle_fields->HM_density = NULL;
  grackle_fields->H2I_density = NULL;
  grackle_fields->H2II_density = NULL;
  /* for primordial_chemistry >= 3 */
  grackle_fields->DI_density = NULL;
  grackle_fields->DII_density = NULL;
  grackle_fields->HDI_density = NULL;
  /* for metal_cooling = 1 */
  grackle_fields->metal_density = NULL;
  /* for use_dust_density_field = 1 */
  grackle_fields->dust_density = NULL;

  /* volumetric heating rate (provide in units [erg s^-1 cm^-3]) */
  grackle_fields->volumetric_heating_rate = NULL;
  /* specific heating rate (provide in units [egs s^-1 g^-1] */
  grackle_fields->specific_heating_rate = NULL;

  /* radiative transfer ionization / dissociation rate fields (provide in units
   * [1/s]) */
  grackle_fields->RT_HI_ionization_rate = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_HeI_ionization_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_HeII_ionization_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_H2_dissociation_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  /* radiative transfer heating rate field
   * (provide in units [erg s^-1 cm^-3] / nHI in cgs) */
  grackle_fields->RT_heating_rate = malloc(FIELD_SIZE * sizeof(gr_float));

  grackle_fields->H2_self_shielding_length = NULL;
  grackle_fields->H2_custom_shielding_factor = NULL;
  grackle_fields->isrf_habing = NULL;

  for (int i = 0; i < FIELD_SIZE; i++) {

    grackle_fields->density[i] = density;
    grackle_fields->internal_energy[i] = internal_energy;

    grackle_fields->HI_density[i] = species_densities[0];
    grackle_fields->HII_density[i] = species_densities[1];
    grackle_fields->HeI_density[i] = species_densities[2];
    grackle_fields->HeII_density[i] = species_densities[3];
    grackle_fields->HeIII_density[i] = species_densities[4];
    /* e_density = electron density*mh/me = n_e * m_h */
    grackle_fields->e_density[i] = species_densities[5];

    /* grackle_fields->HM_density[i] = species_densities[6]; */
    /* grackle_fields->H2I_density[i] = species_densities[7]; */
    /* grackle_fields->H2II_density[i] = species_densities[8]; */
    /* grackle_fields->DI_density[i] = species_densities[9]; */
    /* grackle_fields->DII_density[i] = species_densities[10]; */
    /* grackle_fields->HDI_density[i] = species_densities[11]; */

    /* grackle_fields->metal_density[i] = 0.0; */
    /* solar metallicity */
    /* grackle_chemistry_data.SolarMetalFractionByMass *
     * grackle_fields->density[i]; */

    /* grackle_fields->x_velocity[i] = 0.0; */
    /* grackle_fields->y_velocity[i] = 0.0; */
    /* grackle_fields->z_velocity[i] = 0.0; */

    /* grackle_fields->volumetric_heating_rate[i] = 0.0; */
    /* grackle_fields->specific_heating_rate[i] = 0.0; */

    grackle_fields->RT_heating_rate[i] = iact_rates[0];
    grackle_fields->RT_HI_ionization_rate[i] = iact_rates[1];
    grackle_fields->RT_HeI_ionization_rate[i] = iact_rates[2];
    grackle_fields->RT_HeII_ionization_rate[i] = iact_rates[3];
    grackle_fields->RT_H2_dissociation_rate[i] = iact_rates[4];
  }
}

/**
 * @brief free arrays allocated in grackle_fields.
 *
 * @param grackle_fields grackle fields to clean up
 *
 **/
__attribute__((always_inline)) INLINE static void rt_clean_grackle_fields(
    grackle_field_data *grackle_fields) {

  free(grackle_fields->grid_dimension);
  free(grackle_fields->grid_start);
  free(grackle_fields->grid_end);

  /* initial quantities */
  free(grackle_fields->density);
  free(grackle_fields->internal_energy);
  /* free(grackle_fields->x_velocity); */
  /* free(grackle_fields->y_velocity); */
  /* free(grackle_fields->z_velocity); */

  /* for primordial_chemistry >= 1 */
  free(grackle_fields->HI_density);
  free(grackle_fields->HII_density);
  free(grackle_fields->HeI_density);
  free(grackle_fields->HeII_density);
  free(grackle_fields->HeIII_density);
  free(grackle_fields->e_density);

  /* for primordial_chemistry >= 2 */
  /* free(grackle_fields->HM_density); */
  /* free(grackle_fields->H2I_density); */
  /* free(grackle_fields->H2II_density); */

  /* for primordial_chemistry >= 3 */
  /* free(grackle_fields->DI_density); */
  /* free(grackle_fields->DII_density); */
  /* free(grackle_fields->HDI_density); */

  /* for metal_cooling = 1 */
  /* free(grackle_fields->metal_density); */

  /* for use_dust_density_field = 1 */
  /* free(grackle_fields->dust_density); */

  /* free(grackle_fields->volumetric_heating_rate); */
  /* free(grackle_fields->specific_heating_rate); */

  free(grackle_fields->RT_HI_ionization_rate);
  free(grackle_fields->RT_HeI_ionization_rate);
  free(grackle_fields->RT_HeII_ionization_rate);
  free(grackle_fields->RT_H2_dissociation_rate);
  free(grackle_fields->RT_heating_rate);

  /* free(grackle_fields->H2_self_shielding_length); */
  /* free(grackle_fields->H2_custom_shielding_factor); */
  /* free(grackle_fields->isrf_habing); */
}

/**
 * @brief Write out all available grackle field data for a given index
 * and setup to a file.
 * This function is intended for debugging.
 *
 * @param fp FILE pointer to write into.
 * @param grackle_fields grackle field data
 * @param grackle_chemistry_data grackle chemistry data.
 * @param grackle_units units used by grackle
 * @param field_index grackle field index to print out.
 **/
__attribute__((always_inline)) INLINE static void
rt_write_grackle_setup_and_field(FILE *fp, grackle_field_data grackle_fields,
                                 chemistry_data *grackle_chemistry_data,
                                 code_units *grackle_units, int field_index) {

  fprintf(fp, "Grackle chemistry parameters:\n");

  fprintf(fp, "use_grackle                       = %d\n",
          grackle_chemistry_data->use_grackle);
  fprintf(fp, "with_radiative_cooling            = %d\n",
          grackle_chemistry_data->with_radiative_cooling);
  fprintf(fp, "primordial_chemistry              = %d\n",
          grackle_chemistry_data->primordial_chemistry);
  fprintf(fp, "dust_chemistry                    = %d\n",
          grackle_chemistry_data->dust_chemistry);
  fprintf(fp, "metal_cooling                     = %d\n",
          grackle_chemistry_data->metal_cooling);
  fprintf(fp, "UVbackground                      = %d\n",
          grackle_chemistry_data->UVbackground);
  fprintf(fp, "grackle_data_file                 = %s\n",
          grackle_chemistry_data->grackle_data_file);
  fprintf(fp, "cmb_temperature_floor             = %d\n",
          grackle_chemistry_data->cmb_temperature_floor);
  fprintf(fp, "Gamma                             = %g\n",
          grackle_chemistry_data->Gamma);
  fprintf(fp, "h2_on_dust                        = %d\n",
          grackle_chemistry_data->h2_on_dust);
  fprintf(fp, "use_dust_density_field            = %d\n",
          grackle_chemistry_data->use_dust_density_field);
  fprintf(fp, "dust_recombination_cooling        = %d\n",
          grackle_chemistry_data->dust_recombination_cooling);
  fprintf(fp, "photoelectric_heating             = %d\n",
          grackle_chemistry_data->photoelectric_heating);
  fprintf(fp, "photoelectric_heating_rate        = %g\n",
          grackle_chemistry_data->photoelectric_heating_rate);
  fprintf(fp, "use_isrf_field                    = %d\n",
          grackle_chemistry_data->use_isrf_field);
  fprintf(fp, "interstellar_radiation_field      = %g\n",
          grackle_chemistry_data->interstellar_radiation_field);
  fprintf(fp, "use_volumetric_heating_rate       = %d\n",
          grackle_chemistry_data->use_volumetric_heating_rate);
  fprintf(fp, "use_specific_heating_rate         = %d\n",
          grackle_chemistry_data->use_specific_heating_rate);
  fprintf(fp, "three_body_rate                   = %d\n",
          grackle_chemistry_data->three_body_rate);
  fprintf(fp, "cie_cooling                       = %d\n",
          grackle_chemistry_data->cie_cooling);
  fprintf(fp, "h2_optical_depth_approximation    = %d\n",
          grackle_chemistry_data->h2_optical_depth_approximation);
  fprintf(fp, "ih2co                             = %d\n",
          grackle_chemistry_data->ih2co);
  fprintf(fp, "ipiht                             = %d\n",
          grackle_chemistry_data->ipiht);
  fprintf(fp, "HydrogenFractionByMass            = %g\n",
          grackle_chemistry_data->HydrogenFractionByMass);
  fprintf(fp, "DeuteriumToHydrogenRatio          = %g\n",
          grackle_chemistry_data->DeuteriumToHydrogenRatio);
  fprintf(fp, "SolarMetalFractionByMass          = %g\n",
          grackle_chemistry_data->SolarMetalFractionByMass);
  fprintf(fp, "local_dust_to_gas_ratio           = %g\n",
          grackle_chemistry_data->local_dust_to_gas_ratio);
  fprintf(fp, "NumberOfTemperatureBins           = %d\n",
          grackle_chemistry_data->NumberOfTemperatureBins);
  fprintf(fp, "CaseBRecombination                = %d\n",
          grackle_chemistry_data->CaseBRecombination);
  fprintf(fp, "TemperatureStart                  = %g\n",
          grackle_chemistry_data->TemperatureStart);
  fprintf(fp, "TemperatureEnd                    = %g\n",
          grackle_chemistry_data->TemperatureEnd);
  fprintf(fp, "NumberOfDustTemperatureBins       = %d\n",
          grackle_chemistry_data->NumberOfDustTemperatureBins);
  fprintf(fp, "DustTemperatureStart              = %g\n",
          grackle_chemistry_data->DustTemperatureStart);
  fprintf(fp, "DustTemperatureEnd                = %g\n",
          grackle_chemistry_data->DustTemperatureEnd);
  fprintf(fp, "Compton_xray_heating              = %d\n",
          grackle_chemistry_data->Compton_xray_heating);
  fprintf(fp, "LWbackground_sawtooth_suppression = %d\n",
          grackle_chemistry_data->LWbackground_sawtooth_suppression);
  fprintf(fp, "LWbackground_intensity            = %g\n",
          grackle_chemistry_data->LWbackground_intensity);
  fprintf(fp, "UVbackground_redshift_on          = %g\n",
          grackle_chemistry_data->UVbackground_redshift_on);
  fprintf(fp, "UVbackground_redshift_off         = %g\n",
          grackle_chemistry_data->UVbackground_redshift_off);
  fprintf(fp, "UVbackground_redshift_fullon      = %g\n",
          grackle_chemistry_data->UVbackground_redshift_fullon);
  fprintf(fp, "UVbackground_redshift_drop        = %g\n",
          grackle_chemistry_data->UVbackground_redshift_drop);
  fprintf(fp, "cloudy_electron_fraction_factor   = %g\n",
          grackle_chemistry_data->cloudy_electron_fraction_factor);
  fprintf(fp, "use_radiative_transfer            = %d\n",
          grackle_chemistry_data->use_radiative_transfer);
  fprintf(fp, "radiative_transfer_coupled_rate_solver = %d\n",
          grackle_chemistry_data->radiative_transfer_coupled_rate_solver);
  fprintf(fp, "radiative_transfer_intermediate_step = %d\n",
          grackle_chemistry_data->radiative_transfer_intermediate_step);
  fprintf(fp, "radiative_transfer_hydrogen_only  = %d\n",
          grackle_chemistry_data->radiative_transfer_hydrogen_only);
  fprintf(fp, "self_shielding_method             = %d\n",
          grackle_chemistry_data->self_shielding_method);
  fprintf(fp, "H2_custom_shielding               = %d\n",
          grackle_chemistry_data->H2_custom_shielding);
  fprintf(fp, "H2_self_shielding                 = %d\n",
          grackle_chemistry_data->H2_self_shielding);

  fprintf(fp, "\nUnits:\n");
  fprintf(fp, "a_units               = %g\n", grackle_units->a_units);
  fprintf(fp, "a_value               = %g\n", grackle_units->a_value);
  fprintf(fp, "comoving_coordinates  = %d\n",
          grackle_units->comoving_coordinates);
  fprintf(fp, "density_units         = %g\n", grackle_units->density_units);
  fprintf(fp, "length_units          = %g\n", grackle_units->length_units);
  fprintf(fp, "time_units            = %g\n", grackle_units->time_units);
  fprintf(fp, "velocity_units        = %g\n", grackle_units->velocity_units);

#define rt_print_grackle_field(v) \
  if (grackle_fields.v != NULL)   \
  fprintf(fp, "grackle_fields." #v " = %g\n", grackle_fields.v[field_index])

  fprintf(fp, "\nGrackle field data:\n");
  rt_print_grackle_field(density);
  rt_print_grackle_field(internal_energy);
  rt_print_grackle_field(HI_density);
  rt_print_grackle_field(HII_density);
  rt_print_grackle_field(HeI_density);
  rt_print_grackle_field(HeII_density);
  rt_print_grackle_field(HeIII_density);
  rt_print_grackle_field(e_density);
  rt_print_grackle_field(HM_density);
  rt_print_grackle_field(H2I_density);
  rt_print_grackle_field(H2II_density);
  rt_print_grackle_field(DI_density);
  rt_print_grackle_field(DII_density);
  rt_print_grackle_field(HDI_density);
  rt_print_grackle_field(metal_density);
  rt_print_grackle_field(x_velocity);
  rt_print_grackle_field(y_velocity);
  rt_print_grackle_field(z_velocity);
  rt_print_grackle_field(volumetric_heating_rate);
  rt_print_grackle_field(specific_heating_rate);
  rt_print_grackle_field(RT_HI_ionization_rate);
  rt_print_grackle_field(RT_HeI_ionization_rate);
  rt_print_grackle_field(RT_HeII_ionization_rate);
  rt_print_grackle_field(RT_H2_dissociation_rate);
  rt_print_grackle_field(RT_heating_rate);

#undef rt_print_grackle_field
}

#endif /* SWIFT_RT_GRACKLE_UTILS_H */
