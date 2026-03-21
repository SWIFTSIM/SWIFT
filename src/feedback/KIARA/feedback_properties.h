/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_KIARA_FEEDBACK_PROPERTIES_H
#define SWIFT_KIARA_FEEDBACK_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "hydro_properties.h"

#define NM 5000
#define NZSN 7
#define NMLF 41
#define NZLF 9
#define NZSN1R 5 /* secondary mass ranges */
#define NZSN1Y 7 /* yields */
#define NMSN 32
#define NXSNall 84

#define SN1E_idx(A, B) ((A) * NZSN1Y + B)
#define LFLT_idx(A, B) ((A) * NMLF + B)
#define SN1R_idx(A, B) ((A) * NM + B)
#define SN2R_idx(A, B) ((A) * NM + B)
#define SWR_idx(A, B) ((A) * NM + B)
#define SN2E_idx(A, B, C) ((A) * NZSN * NM + (B) * NM + C)

#define LINEAR_INTERPOLATION(x1, y1, x2, y2, x) \
  (((y2 - y1) / (x2 - x1)) * (x - x1) + y1)
#define LOG_INTERPOLATION(x, x2, x1) \
  ((log10(x2) - log10(x)) / (log10(x2) - log10(x1)))

enum kiara_metal_boosting {
  kiara_metal_boosting_off,
  kiara_metal_boosting_vwind,
  kiara_metal_boosting_eta,
  kiara_metal_boosting_both
};

/* Chem5 tracks A LOT of elements but we will just map the standard 11 back */
enum chem5_element {
  chem5_element_Z = 0,
  chem5_element_H,
  chem5_element_He,
  chem5_element_Li,
  chem5_element_Be,
  chem5_element_B,
  chem5_element_C,
  chem5_element_N,
  chem5_element_O,
  chem5_element_F,
  chem5_element_Ne,
  chem5_element_Na,
  chem5_element_Mg,
  chem5_element_Al,
  chem5_element_Si,
  chem5_element_P,
  chem5_element_S,
  chem5_element_Cl,
  chem5_element_Ar,
  chem5_element_K,
  chem5_element_Ca,
  chem5_element_Sc,
  chem5_element_Ti,
  chem5_element_V,
  chem5_element_Cr,
  chem5_element_Mn,
  chem5_element_Fe,
  chem5_element_Co,
  chem5_element_Ni,
  chem5_element_Cu,
  chem5_element_Zn,
  chem5_element_Ga,
  chem5_element_Ge,
  chem5_element_unknown,
  chem5_element_count,
  chem5_dummy1,
  chem5_dummy2,
  chem5_NXSN
};

/**
 * @brief Stores the yield tables
 */
struct feedback_tables {
  double *LFLT;
  double *LFLM;
  double *LFLZ;
  double *LFLT2;
  double *SWR;
  double *SN2E;
  double *SN2R;
  double *SN1R;
  double *SNLM;
  double *SNLZ;
  double *SNLZ1R;
  double *SN1E;
  double *SNLZ1Y;
};

/**
 * @brief Properties of the KIARA feedback model.
 */
struct feedback_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we depositing energy from HN directly from Chem5? */
  int with_HN_energy_from_chem5;

  /*! Are we depositing energy from SNII directly from Chem5? */
  int with_SNII_energy_from_chem5;

  /*! Are we depositing energy from SNIa directly from Chem5? */
  int with_SNIa_energy_from_chem5;

  /*! If time since last chemical enrichment is above this value times the
   * current stellar age, recompute */
  float stellar_enrichment_frequency;

  /* ------------ Yield tables    ----------------- */

  struct feedback_tables tables;

  /* Conversion indices from Chem5 to Swift */
  /* 0-Z, 2-He, 6-C , 7-N, 8-O, 10-Ne, 12-Mg, 14-Si, 26-Fe */
  int element_index_conversions[chemistry_element_count];

  /* Location of feedback tables */
  char tables_path[200];

  /* ------------- Conversion factors --------------- */

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! The mass of the sun in g */
  double solar_mass_in_g;

  /*! Conversion factor from internal mass unit to solar mass */
  double solar_mass_to_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Conversion factor from temperature to internal energy */
  double temp_to_u_factor;

  /*! Conversion factor from km/s to cm/s */
  double kms_to_cms;

  /*! Factor to convert km/s to internal units */
  double kms_to_internal;

  /*! Convert internal units to kpc */
  double length_to_kpc;

  /*! Convert internal time to Myr */
  double time_to_Myr;

  /*! Convert internal time to yr */
  double time_to_yr;

  /*! Convert code energy units to cgs */
  double energy_to_cgs;

  /*! Convert temperature in K to internal */
  double T_to_internal;

  /* ------------ Enrichment sampling properties ------------ */

  /*! Star age above which the enrichment will be downsampled (in internal
   * units) */
  double stellar_evolution_age_cut;

  /*! Number of time-steps in-between two enrichment events */
  int stellar_evolution_sampling_rate;

  /* ------------ Kinetic feedback properties --------------- */

  /*! Velocity normalization */
  float FIRE_velocity_normalization;

  /*! FIRE velocity slope */
  float FIRE_velocity_slope;

  /*! Normalization for the mass loading curve */
  float FIRE_eta_normalization;

  /*! The location (in internal mass units) where the break in the
   * mass loading curve occurs */
  float FIRE_eta_break;

  /*! The power-law slope of eta below FIRE_eta_break */
  float FIRE_eta_lower_slope;

  /*! The power-law slope of eta above FIRE_eta_break */
  float FIRE_eta_upper_slope;

  /*! The wind speed of stellar feedback suppressed above this z */
  float wind_velocity_suppression_redshift;

  /*! The mass loading factor of stellar feedback suppressed above this z */
  float wind_eta_suppression_redshift;

  /*! Maxiumum multiple of SNII energy that is available to launch winds */
  float SNII_energy_multiplier;

  /*! For KIARA, the radius from within which to launch wind */
  float kick_radius_over_h;

  /*! For KIARA, the maximum fraction of gas mass within kernel to launch */
  float max_frac_of_kernel_to_launch;

  /*! For KIARA, weight the launch probability by the particles's SFR (0/1) */
  int use_sfr_weighted_launch;

  /*! Flag to set feedback boost at low Z:
   * 0=Off, 1=vwind boost, 2=eta boost, 3=both boost */
  int metal_dependent_vwind;

  /*! The minimum galaxy stellar mass in internal units */
  float minimum_galaxy_stellar_mass;

  /*! The number of star particles when a galaxy is considered resolved */
  int galaxy_particle_resolution_count;

  /*! Floor for the eta suppression factor */
  float eta_suppression_factor_floor;

  /*! Direction to launch wind: 0=random, 1=v x a, 2=outwards */
  float kick_direction_flag;

  /*! Added scatter to the wind velocities */
  float kick_velocity_scatter;

  /*! max decoupling time is (this factor) * current Hubble time */
  float wind_decouple_time_factor;

  /*! Density (cgs) above which recoupling considers it within ISM */
  float recouple_ism_density_nH_cgs;

  /*! Factor (<1) below ISM density below which to recouple */
  float recouple_density_factor;

  /*! The internal energy corresponding to the unheated wind temperature */
  float cold_wind_internal_energy;

  /*! The internal energy corresponding to the heated wind temperature */
  float hot_wind_internal_energy;

  /*! Whether the firehose wind model is on */
  char use_firehose_model;

  /*! Early stellar feedback alpha value from Keller et al 2022 */
  float early_stellar_feedback_alpha;

  /*! Early stellar feedback term for SF efficiency Keller et al 2022 */
  float early_stellar_feedback_epsterm;

  /*! Early stellar feedback t_fb from Keller et al 2022 */
  float early_stellar_feedback_tfb;

  /*! Early stellar feedback t_fb inverse */
  float early_stellar_feedback_tfb_inv;

  /*! Maximum mass growth factor from stellar feedback for a particle */
  float max_mass_increase_factor;

  /*! Maximum internal energy growth factor from stellar feedback to particle */
  float max_energy_increase_factor;

  /*! Minimum internal energy loss factor from momentum exchange */
  // float min_energy_decrease_factor;

  /*! Flag to add heat to destroy cold gas in the ISM from SNIa gas */
  int SNIa_add_heat_to_ISM;

  /*! Avoid floating point errors when comparing internal energies in the ISM */
  float SNIa_add_heat_to_ISM_tolerance;

  /* ------------ Chem5 Default Parameters --------------- */

  /*! Which IMF? Kroupa=0, Chabrier=1, Else=2 */
  int imf;

  /*! Solar H */
  double H_mf;

  /*! Solar He */
  double He_mf;

  /*! Solar Z */
  double Z_mf;

  /*! Solar O */
  double O_mf;

  /*! Solar Fe */
  double Fe_mf;

  /*! IMF parameter */
  double ximf;

  /*! Upper limit for IMF integration */
  double M_u;

  /*! Lower limit for IMF integration */
  double M_l;

  /*! IMF parameter */
  double ximf3;

  /*! >= M_u */
  double M_u3;

  /*! >= M_l */
  double M_l3;

  /*! If set greater than zero, activates Pop3 stars */
  double zmax3;

  /*! Upper limit on IMF integration */
  double M_u2;

  /*! Lower limit on IMF integration */
  double M_l2;

  /*! binary parameter for SNIa */
  double b_rg;

  /*! binary parameter for SNIa */
  double b_ms;

  /*! Energy in supernova (Ia) */
  double E_sn1;

  /*! Energy in supernova (II) */
  double E_sw;

#if COOLING_GRACKLE_MODE >= 2
  /* ------------ Dust Efficiency Tables --------------- */

  /*! dust condensation efficiency for C/O>1 */
  double delta_AGBCOG1[chemistry_element_count];

  /*! dust condensation efficiency for C/O<1 */
  double delta_AGBCOL1[chemistry_element_count];

  /*! dust condensation efficiency from SNII */
  double delta_SNII[chemistry_element_count];

  /*! max fraction of metals locked into dust */
  float max_dust_fraction;

  /*! Rolling value for number of SNe is smoothed over this timescale
   * in Myr (0 for instantaneous) */
  float SNe_smoothing_time_in_Myr;
#endif

  /*! chem5 metal yield multiplier */
  float metal_yield_multiplier;
};

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);

#endif /* SWIFT_KIARA_FEEDBACK_PROPERTIES_H */
