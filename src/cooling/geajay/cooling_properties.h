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
#ifndef SWIFT_COOLING_PROPERTIES_GEAJAY_H
#define SWIFT_COOLING_PROPERTIES_GEAJAY_H

/* Hide the deprecated separate-array functions of geajay. */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC

/* geajay: the chemistry API and the particle API */
#include <geajay.h>
#include <geajay_cosmic_ray.h>
#include <geajay_particle.h>

/**
 * @file src/cooling/geajay/cooling_properties.h
 * @brief Properties of the cooling with the geajay library.
 */

/**
 * @brief Self-shielding of the UV background (GeajayCooling:uv_self_shielding).
 *
 * The values 0 to 3 are those of geajay_uv_self_shielding.
 */
enum cooling_uv_self_shielding {
  /*! Density threshold, handled by SWIFT (option "density_threshold") */
  cooling_uv_self_shielding_density_threshold = -1,
  /*! None */
  cooling_uv_self_shielding_none = 0,
  /*! HI only */
  cooling_uv_self_shielding_HI = 1,
  /*! HI and HeI, HeII thin */
  cooling_uv_self_shielding_HI_HeI = 2,
  /*! HI and HeI, no HeII ionisation and heating */
  cooling_uv_self_shielding_HI_HeI_no_HeII = 3
};

/**
 * @brief Length of the H2 self-shielding column
 * (GeajayCooling:h2_shielding_length_method).
 *
 * The values are those of geajay_h2_shielding_length_method.
 */
enum cooling_h2_shielding_length {
  /*! No H2 self-shielding */
  cooling_h2_shielding_length_none = 0,
  /*! Sobolev-like length (rejected at start-up) */
  cooling_h2_shielding_length_Sobolev = 1,
  /*! Length supplied per particle */
  cooling_h2_shielding_length_given = 2,
  /*! Local Jeans length, computed by geajay */
  cooling_h2_shielding_length_Jeans = 3
};

/* The structure is dumped raw in the restart files */
_Static_assert(sizeof(enum cooling_uv_self_shielding) == sizeof(int) &&
                   sizeof(enum cooling_h2_shielding_length) == sizeof(int),
               "The shielding enums must keep the size of an int.");

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Filename of the Cloudy Table */
  char cloudy_table[200];

  /*! Enable/Disable UV backgroud */
  int with_uv_background;

  /*! Chemistry network */
  int primordial_chemistry;

  /*! Set the three-body reaction rate fit (see geajay_options) */
  int H2_three_body_rate;

  /*! Enable/disable H2 collision-induced emission cooling from Ripamonti & Abel
   * (2004) */
  int H2_cie_cooling;

  /*! Flag to enable H2 formation on dust grains */
  int H2_on_dust;

  /*! Flag to enable heating from H2 photodissociation and UV pumping */
  int H2_photodissociation_heating;

  /*! Compton heating by the X-ray background
   * (GeajayCooling:compton_xray_heating, chemistry_data.Compton_xray_heating).
   */
  int compton_xray_heating;

  /*! The ratio of total dust mass to gas mass in the local Universe. */
  double local_dust_to_gas_ratio;

  /*! Enable/disable CMB temperature floor */
  int cmb_temperature_floor;

  /*! Flag to enable the AGORA redshift dependent CMB floor:
                    T_CMB = T_CMB,0 * (1+z)
      This floor is different from the `cmb_temperature_floor`, which is
      internal to geajay.
  */
  char agora_cmb_temperature_floor;

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  geajay_code_units units;

  /*! geajay chemistry parameters */
  geajay_chemistry_data chemistry_data;

  /*! geajay rate storage (kept alive for as long as the context) */
  geajay_chemistry_data_storage chemistry_rates;

  /*! Chemistry parameters and unit factors, shared by all threads between two
   * calls of cooling_update(). Created in cooling_init_geajay(), not valid
   * after a restart read. */
  geajay_context *context;

  /*! Where the species live in the #xpart, checked against the context. */
  geajay_species_layout layout;

  /*! Cosmic-ray settings (GeajayCooling:cosmic_rays and the options after
   * it). */
  geajay_cosmic_ray_settings cosmic_ray_settings;

  /*! Context of the reduced cosmic-ray network over #context; NULL without
   * cosmic rays. Created in cooling_init_geajay(), not valid after a restart
   * read. */
  geajay_cosmic_ray_reduced_context *cosmic_ray_context;

  /*! Background rules (GeajayCooling:background_mode), a
   * #geajay_background_mode. With "grackle" the context is built from the
   * chemistry as before; with "per_band" from geajay_options. */
  int background_mode;

  /*! Per-band choices of "per_band" (geajay_choices.h); -1 for "default". */
  int ionising_background_model, lw_background_model, isrf_background_model;

  /*! Background G_0 of isrf_background_model "constant" [Habing]; negative:
   * unset. */
  double isrf_background_habing;

  /*! Dust attenuation of the background, a #geajay_background_attenuation. */
  int background_attenuation;

  /*! Constant LW background J21 (Grackle LWbackground_intensity); negative:
   * the library default. */
  double lw_background_j21;

  /*! G_0 is read from the particle field (Grackle use_isrf_field). */
  int isrf_from_field;

  /*! Dust cross-sections per H nucleon of the LW (11.2-13.6 eV) and PE
   * (6-11.2 eV) bands [cm^2]. */
  double dust_cross_section_lw_cgs, dust_cross_section_fuv_cgs;

  /*! Converts rho * path (internal units) to the H nucleon column [cm^-2]. */
  double column_factor;

  /*! Initial cosmic-ray value of every particle, in the unit of the input
   * (GeajayCooling:cosmic_ray_initial_value); 0 with the input "constant". */
  float cosmic_ray_initial_value;

  /*! Enable/Disable metal cooling */
  int with_metal_cooling;

  /*! Arrays of ionization and heating rates are provided */
  int use_radiative_transfer;

  /*! RT_heating_rate (in IU) */
  float RT_heating_rate;

  /*! RT_HI_ionization_rate (in IU) */
  float RT_HI_ionization_rate;

  /*! RT_HeI_ionization_rate (in IU) */
  float RT_HeI_ionization_rate;

  /*! RT_HeII_ionization_rate (in IU) */
  float RT_HeII_ionization_rate;

  /*! RT_H2_dissociation_rate (in IU) */
  float RT_H2_dissociation_rate;

  /*! Volumetric heating rates */
  float volumetric_heating_rates;

  /*! Specific heating rates */
  float specific_heating_rates;

  /*! Hydrogen fraction by mass (negative in the parameter file: the library
   * value, resolved at initialisation) */
  double HydrogenFractionByMass;

  /*! initial nHII to nH ratio (number density ratio) */
  float initial_nHII_to_nH_ratio;

  /*! initial nHeI to nH ratio (number density ratio) */
  float initial_nHeI_to_nH_ratio;

  /*! initial nHeII to nH ratio (number density ratio) */
  float initial_nHeII_to_nH_ratio;

  /*! initial nHeIII to nH ratio (number density ratio) */
  float initial_nHeIII_to_nH_ratio;

  /*! initial nDI to nH ratio (number density ratio) */
  float initial_nDI_to_nH_ratio;

  /*! initial nDII to nH ratio (number density ratio) */
  float initial_nDII_to_nH_ratio;

  /*! initial nHM to nH ratio (number density ratio) */
  float initial_nHM_to_nH_ratio;

  /*! initial nH2I to nH ratio (number density ratio) */
  float initial_nH2I_to_nH_ratio;

  /*! initial nH2II to nH ratio (number density ratio) */
  float initial_nH2II_to_nH_ratio;

  /*! initial nHDI to nH ratio (number density ratio) */
  float initial_nHDI_to_nH_ratio;

  /*! Self shielding method (1 -> 3 for the library's ones, 0 for none and -1
   * for a density threshold) */
  enum cooling_uv_self_shielding self_shielding_method;

  /*! Self shielding threshold */
  float self_shielding_threshold;

  /*! H2 self-shielding mode (chemistry_data.H2_self_shielding):
   * 0 = off, 2 = shielding length supplied per particle, half the column
   * path set by H2_self_shielding_path (one kernel support radius
   * kernel_gamma * h, physical, by default), 3 = local Jeans length
   * (computed internally by geajay).
   * Mode 1 (Sobolev-like, differences neighbouring grid cells) is rejected
   * at start-up, since SWIFT calls geajay on a single particle. Unrelated
   * to self_shielding_method above, which is a separate field
   * gating UV-background self-shielding. */
  enum cooling_h2_shielding_length H2_self_shielding;

  /*! Mode-2 H2 shielding path in kernel support radii kernel_gamma * h: 1
   * for "kernel_radius" (default), 2 for "kernel_diameter"
   * (GeajayCooling:H2_self_shielding_path). The column is N_H2 = 2 n_H2 l, so
   * the length handed to geajay is half this path. */
  float H2_self_shielding_path_in_kernel_radii;

  /*! Maximal number of sub-cycle iterations per solve
   * (GeajayCooling:max_substeps, chemistry_data.max_iterations). */
  int max_step;

  /*! Number of retries of a failed solve
   * (GeajayCooling:subcycle_on_failure). Retry i splits the time-step into
   * 2^i sub-steps. 0 disables the retries. */
  int subcycle_on_failure;

  /*! Largest fraction of the internal energy that one sub-step of the solver
   * may change, in (0, 1] (GeajayCooling:energy_change_limit,
   * chemistry_data.energy_change_limit). Grackle fixes 0.1. */
  double energy_change_limit;

  /*! Duration for switching off cooling after an event (e.g. supernovae) */
  double thermal_time;

  /*! Lookup beyond the density range of the Cloudy table
   * (GeajayCooling:table_density_edge). */
  geajay_table_edge table_density_edge;

  /*! Proton mass [g] */
  double proton_mass_cgs;

  /*! Hydrogen number density range of the Cloudy table [cm^-3], set by
   * cooling_init_geajay(). [0, inf] without a table. */
  double table_density_min_Hpcm3, table_density_max_Hpcm3;

  /*! Physical density (internal units) whose X_H rho / m_p is
   * table_density_max_Hpcm3. */
  double table_edge_density;

  /*! Largest physical density passed to geajay, in internal units
   * (GeajayCooling:maximal_density_Hpcm3). Zero or negative: no cap. */
  double cooling_density_max;

  /*! Debug/test-only: skip cooling_cool_part()'s per-particle energy
      update entirely (leaving hydro's own internal-energy derivative
      untouched), while cooling_init() and everything it resolves (table
      loading, chemistry_data, including local_dust_to_gas_ratio) still
      run completely normally. Never set in a production run. */
  int disable_cooling_for_debugging;
};

#endif /* SWIFT_COOLING_PROPERTIES_GEAJAY_H */
