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
#ifndef SWIFT_CHEMISTRY_STRUCT_KIARA_H
#define SWIFT_CHEMISTRY_STRUCT_KIARA_H

#define FIREHOSE_COOLLIM 0.1f            /* -u_new / u */
#define FIREHOSE_HEATLIM 10.f            /* +u_new / u */
#define FIREHOSE_EPSILON_TOLERANCE 1.e-6 /* Minimum rel. difference to add */

/**
 * @brief The individual elements traced in the KIARA model.
 */
enum chemistry_element {
  chemistry_element_H = 0,
  chemistry_element_He,
  chemistry_element_C,
  chemistry_element_N,
  chemistry_element_O,
  chemistry_element_Ne,
  chemistry_element_Mg,
  chemistry_element_Si,
  chemistry_element_S,
  chemistry_element_Ca,
  chemistry_element_Fe,
  chemistry_element_count
};

#if COOLING_GRACKLE_MODE >= 2
/**
 * @brief The individual elements traced in the Grackle dust model.
 */
enum dust_element {
  dust_element_C,
  dust_element_O,
  dust_element_Mg,
  dust_element_Si,
  dust_element_S,
  dust_element_Ca,
  dust_element_Fe,
  dust_element_count
};
#endif

/**
 * @brief Global chemical abundance information in the KIARA model.
 */
struct chemistry_global_data {

  /*! Fraction of the particle mass in given elements at the start of the run */
  float initial_metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals at the start of the run */
  float initial_metal_mass_fraction_total;

  /*! Is metal diffusion turned on? */
  int diffusion_flag;

  /*! The timestep beta value from Parshikov & Medin 2002 equation 41 */
  float diffusion_beta;

  /*! The minimum time step size in internal units for diffusion */
  float time_step_min;

  /*! A limiter for how much Z/Z_init can be transferred (~0.25) */
  float max_fractional_Z_transfer;

  /*! The metal diffusion coefficient (Smag ~0.23) */
  float C_Smagorinsky;

  /*! Use Firehose wind model (1) or standard decoupled winds (0) */
  int use_firehose_wind_model;

  /*! Firehose wind model maximum density */
  float firehose_ambient_rho_max;

  /*! Firehose wind model minimum thermal energy */
  float firehose_u_floor;

  /*! Firehose wind model Mach number threshold for recoupling */
  float firehose_recoupling_mach;

  /*! Firehose wind model thermal energy ratio threshold for recoupling */
  float firehose_recoupling_u_factor;

  /*! Firehose wind model mixing mass threshold for recoupling */
  float firehose_recoupling_fmix;

  /*! Firehose threshold relative velocity (km/s) above which model is turned
   * off */
  float firehose_max_velocity;

  /*! Firehose maximum fraction of particles' mass that can be mixed in a single
   * step */
  float firehose_max_fmix_per_step;

  /*! Dust sputtering constant */
  float dust_sputtering_const;

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Converts temperature to internal energy */
  float temp_to_u_factor;

  /*! Converst temperature to internal units */
  float T_to_internal;

  /*! Factor to convert km/s to internal units */
  float kms_to_internal;

  /*! Convert internal units to kpc */
  float length_to_kpc;

  /*! Convert internal time to Myr */
  float time_to_Myr;
};

/**
 * @brief Chemical abundances traced by the #part in the KIARA model.
 */
struct chemistry_part_data {

  /*! Fraction of the particle mass in a given element */
  float metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals */
  float metal_mass_fraction_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Fraction of total gas mass in metals coming from SNIa */
  float metal_mass_fraction_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Fraction of total gas mass in metals coming from AGB */
  float metal_mass_fraction_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Fraction of total gas mass in metals coming from SNII */
  float metal_mass_fraction_from_SNII;

  /*! Fraction of total gas mass in Iron coming from SNIa */
  float iron_mass_fraction_from_SNIa;

  /*! Diffusion coefficient */
  float diffusion_coefficient;

  /*! Variation of the total metal mass */
  float dZ_dt_total;

  /*! Variation of the metal mass by element */
  float dZ_dt[chemistry_element_count];

  /*! Velocity shear tensor in internal and physical units. */
  float shear_tensor[3][3];

#if COOLING_GRACKLE_MODE >= 2
  /*! SFR density (physical) within smoothing kernel needed for G0 calculation
   */
  float local_sfr_density;
#endif

  /*! Firehose ambient gas thermal energy */
  float u_ambient;

  /*! Firehose ambient gas density */
  float rho_ambient;

  /*! Weighting factor for ambient thermal energy sum */
  float w_ambient;

  /*! Firehose radius of outflowing stream */
  float radius_stream;

  /*! Firehose initial mass of the stream */
  float exchanged_mass;

  /*! Firehose exchanged mass this step */
  float dm;

  /*! Firehose exchanged metal fractions */
  float dm_Z[chemistry_element_count];

  /*! Firehose exchanged dust mass */
  float dm_dust;

  /*! Firehose exchanged dust mass metals */
  float dm_dust_Z[chemistry_element_count];

  /*! Firehose exchanged internal energy, internal units */
  double du;

  /*! Firehose exchanged velocities, internal units */
  float dv[3];
};

#define chemistry_spart_data chemistry_part_data

/**
 * @brief Chemical abundances traced by the #bpart in the KIARA model.
 */
struct chemistry_bpart_data {

  /*! Mass in a given element */
  float metal_mass[chemistry_element_count];

  /*! Mass in *all* metals */
  float metal_mass_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Metal mass coming from SNIa */
  float metal_mass_from_SNIa;

  /*! Metal mass coming from AGB */
  float metal_mass_from_AGB;

  /*! Metal mass coming from SNII */
  float metal_mass_from_SNII;

  /*! Iron mass coming from SNIa */
  float iron_mass_from_SNIa;

  /*! Metallicity of converted part. */
  float formation_metallicity;
};

/**
 * @brief Chemical abundances traced by the #sink in the KIARA model.
 *
 * Nothing here.
 */
struct chemistry_sink_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_KIARA_H */
