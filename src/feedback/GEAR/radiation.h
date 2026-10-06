/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_GEAR_H
#define SWIFT_RADIATION_GEAR_H

/**
 * @file src/feedback/GEAR/radiation.h
 * @brief Subgrid radiation feedback for GEAR. This file contains functions to
 * compute quantities for the radiation feedback.
 */

#include "../../feedback_properties.h"
#include "cooling_properties.h"
#include "hdf5_functions.h"
#include "hydro.h"
#include "part.h"
#include "physical_constants.h"
#include "radiation_isrf.h"
#include "stellar_evolution_struct.h"
#include "units.h"

/*! Scale of the dot_N_ion tables: they store (raw value / this factor) to fit
    float, and every reader multiplies back. */
#define RADIATION_DOT_N_ION_TABLE_SCALING 1e50

/*! Lifetime of an ionization tag, in rebuild intervals. Must exceed 1 so the
    tag outlives the gap to its star's next rebuild pass. */
#define RADIATION_TAG_LIFETIME_INTERVALS 2.0

/*! Lifetime of an LW/PE illumination episode, in integer timesteps of the
    illuminating star. Must exceed 1 so it outlives the gap to the star's next
    visit. */
#define RADIATION_ISRF_TAG_LIFETIME_INTERVALS 2

/*! Decay length, in particle smoothing lengths, of the negativity-triggered
    artificial-dissipation coefficient. Sets the transient, not the steady
    state (ISRF_dissipation_alpha_max does that). */
#define RADIATION_ISRF_DISSIPATION_DECAY_LENGTH 5.0f

/*! Ceiling, in rebuild cadences, on the interval the per-pass photon budget
    is integrated over. Photons emitted while the star's cell held no gas
    escaped and must not be handed to the next pass. */
#define HII_DT_BACK_MAX_INTERVALS 2.0

/*! Floor applied to a Data/Radiation CGS value before log10(). It matches
    pychem's `_LOG_FLOOR` bit for bit, so a native-zero entry means the same in
    both codes. Applied to the CGS value, before the unit conversion. */
#define RADIATION_LOG_FLOOR_CGS 1e-300

/*! Mean mass per hydrogen nucleon (He included), in proton masses. It pairs
    with the per-hydrogen dust cross-sections below. It is not the run's
    hydrogen fraction: the dust-to-gas mass ratio cancels the composition. */
#define RADIATION_MU_H 1.4

/*! Dust cross-section per hydrogen nucleon, cm^2 (Kim et al. 2023): PE
    (6-11.2 eV) and Lyman-Werner (11.2-13.6 eV) bands. */
#define RADIATION_SIGMA_D_PE_CGS 9e-22
#define RADIATION_SIGMA_D_LW_CGS 1.5e-21

/*! Grackle's solar metal mass fraction (SolarMetalFractionByMass), used to
    scale the dust-to-gas ratio as Grackle does. */
#define RADIATION_GRACKLE_SOLAR_METAL_FRACTION 0.01295

/*! Grackle's default local_dust_to_gas_ratio (Pollack et al. 1994), used when
    GrackleCooling:local_dust_to_gas_ratio is unset. */
#define RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO 0.009387

/*! Habing flux, erg/s/cm^2: G0=1 is this flux over the PE+LW bands. */
#define RADIATION_HABING_FLUX_CGS 1.6e-3

/*! H2 Lyman-Werner photodissociation coefficient sigma_H2/E_LW, cm^2 erg^-1
    (Sternberg et al. 2014, ApJ 790:10). The particle carries only the LW
    energy flux, so the rate depends on this quotient alone. Agreement with the
    anchor is at the 10% level. See
    theory/GEAR/Radiation/verify_sigma_h2_lw_sternberg2014.py. */
#define RADIATION_SIGMA_H2_OVER_E_LW_CGS 1.2847106348798106e-07

/*! Representative Lyman-Werner photon energy, eV (Kim et al. 2023, Table 3).
    Not used by the dissociation rate. */
#define RADIATION_LW_PHOTON_ENERGY_EV 12.2

/*! Effective H2 LW cross section, cm^2: the quotient above times
    #RADIATION_LW_PHOTON_ENERGY_EV in erg. The rate does not read it. It is a
    literal so the check scripts can parse it; a unit test asserts it equals
    the product. */
#define RADIATION_SIGMA_H2_LW_CGS 2.5111667e-18

/*! Metallicity mass fraction at which the diagnostic mean LW photon energy is
    read from a 2D table (#radiation_lw_photon_energy_cgs). The gas particle
    is source-anonymous, so one representative rung is used. */
#define RADIATION_LW_PHOTON_ENERGY_REFERENCE_METALLICITY 0.014

/*! Relative epsilon by which a 2D IMF-integrated query mass is nudged below
    the table's top mass edge, so an exact mass_max query takes the blended
    branch of interpolate_2d(). */
#define RADIATION_2D_EDGE_EPS 1e-5f

/*! PE/LW band lower edges, eV. Photons redshift downward through these
    energies (fixed in physical energy). See
    theory/GEAR/Radiation/02_fuv_isrf.tex. */
#define RADIATION_PE_BAND_LOWER_EDGE_EV 6.0
#define RADIATION_LW_BAND_LOWER_EDGE_EV 11.2

/*! The band lower edges in erg (1 eV = 1.602176634e-12 erg), hardcoded because
    eV is not an internal unit. They enter only a dimensionless ratio. */
#define RADIATION_PE_BAND_LOWER_EDGE_CGS 9.6130598e-12
#define RADIATION_LW_BAND_LOWER_EDGE_CGS 1.7944378e-11

/*! Fallback band-edge weights lambda_E(PE), lambda_E(LW), lambda_N(LW). They
    initialise the #feedback_props band_edge weights and are used when the
    table's Integrated_L_PE or Integrated_L_LW denominator vanishes.
    lambda_E(b) = 1 + Lambda_b * E_lo(b)/<E>_b (theory/GEAR/Radiation/
    02_fuv_isrf.tex). The values are the young-population end of the spectral
    family, which dominates the LW luminosity. */
#define RADIATION_BAND_EDGE_WEIGHT_PE_DEFAULT 2.154
#define RADIATION_BAND_EDGE_WEIGHT_LW_DEFAULT 6.508
#define RADIATION_BAND_EDGE_PHOTON_WEIGHT_LW_DEFAULT 6.0

/**
 * @brief Read-time grid metadata shared by the datasets of a Data/Radiation
 * HDF5 group. Rebuilt on every read, not stored in #radiation.
 */
struct radiation_grid_metadata {
  /*! "M" (mass-only) or "M,Z" (mass x metallicity). */
  char dimensionality[8];

  /*! Is this a 2D ("M,Z") table? */
  int is_2d;

  /*! log10(mass grid minimum). */
  float log_mass_min;

  /*! log10 mass grid step. */
  float mass_step;

  /*! Number of mass grid points. */
  int n_mass;

  /*! Number of metallicity grid points (0 for a 1D table). */
  int n_metallicity;

  /*! Metallicity grid (mass fraction, NULL for 1D). Increasing and positive,
      not log-uniform, so 2D tables keep it as their own axis. */
  float *metallicity;

  /*! Mass-axis boundary condition for "Luminosity" (2D only). The metallicity
      axis always clamps. */
  enum interpolate_boundary_condition edge_policy_luminosity;

  /*! Mass-axis boundary condition for "Q_H" (2D only), read from the generic
      edge_policy_q_h_below/above attributes written by pychem. */
  enum interpolate_boundary_condition edge_policy_q_h;

  /*! Mass-axis boundary condition for "DotEExcess" (2D only). It may differ
      from #edge_policy_q_h. */
  enum interpolate_boundary_condition edge_policy_dot_e_excess;

  /*! Mass-axis boundary condition for "Teff" (2D only, if present). */
  enum interpolate_boundary_condition edge_policy_teff;

  /*! Mass-axis boundary condition for "L_PE" (2D only, if present). It has its
      own attribute pair, separate from #edge_policy_luminosity. */
  enum interpolate_boundary_condition edge_policy_l_pe;

  /*! Mass-axis boundary condition for "L_LW" (2D only, if present). */
  enum interpolate_boundary_condition edge_policy_l_lw;
};

double radiation_get_part_number_hydrogen_atoms(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double radiation_get_part_number_neutral_hydrogen_atoms(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double radiation_get_part_rate_to_fully_ionize(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double radiation_get_part_ionized_internal_energy(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double radiation_get_case_b_recombination_coefficient_cgs(const double T);

double radiation_get_T_collisional_K(const double Z);

void radiation_tag_part_as_ionized(struct part *p, struct xpart *xpj,
                                   long long star_id, double end_time,
                                   float excess_photon_energy_HI,
                                   float photoionization_rate_HI);
void radiation_reset_part_ionized_tag(struct part *p, struct xpart *xpj);
char radiation_is_part_tagged_as_ionized(const struct part *p,
                                         const struct xpart *xpj);
double radiation_get_part_ionized_end_time(const struct part *p,
                                           const struct xpart *xpj);
void radiation_reset_part_ISRF_illumination_tag(struct part *p,
                                                const struct engine *e);
long long radiation_get_part_ionized_star_id(const struct part *p,
                                             const struct xpart *xpj);
float radiation_get_part_excess_photon_energy_HI(const struct part *p,
                                                 const struct xpart *xpj);
float radiation_get_part_photoionization_rate_coefficient(
    const struct part *p, const struct xpart *xpj);
double radiation_get_photoionization_rate_coefficient_from_flux_HI(
    const struct unit_system *us, const double ionizing_flux_HI);
double radiation_get_part_isrf_habing(const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct cosmology *cosmo,
                                      const struct part *p);
double radiation_get_part_LW_dissociation_rate_internal(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct part *p);
void radiation_set_ionizing_photon_rate(struct spart *sp,
                                        double dot_N_ion_total,
                                        int n_HII_pixels);
void radiation_zero_spart_output(struct spart *sp);
void radiation_open_ionizing_photon_budget(struct spart *sp, double dt_back);
void radiation_resync_ionizing_photon_rate_cache(struct spart *sp);
void radiation_consume_ionizing_photons(struct spart *sp, int pixel,
                                        double Delta_N_ion);
float radiation_get_comoving_gas_column_density_at_star(const struct spart *sp);

float radiation_get_star_physical_radiation_pressure(
    const struct spart *sp, const float Delta_t,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo);

/******************************************************************************/
/* Functions to deal with integrated data over an IMF. These functions read,
   interpolate and integrate. */
/******************************************************************************/
/*! Names of the radiation table's optional provenance attributes, in the
    order #radiation.table_source stores them. Defined in
    radiation_table_io.c, which reads them. */
extern const char
    *const radiation_table_source_keys[RADIATION_TABLE_SOURCE_COUNT];

void radiation_print(const struct radiation *rad);
void radiation_init(struct radiation *rad, struct swift_params *params,
                    const struct stellar_model *sm,
                    const struct unit_system *us,
                    const struct phys_const *phys_const);
void radiation_dump(const struct radiation *rad, FILE *stream,
                    const struct stellar_model *sm);
void radiation_restore(struct radiation *rad, FILE *stream,
                       const struct stellar_model *sm,
                       const struct unit_system *us,
                       const struct phys_const *phys_const,
                       const char with_radiation);
void radiation_clean(struct radiation *rad);
void radiation_zero_pointers(struct radiation *rad);

float radiation_get_luminosities_from_integral(const struct radiation *rad,
                                               float log_m1, float log_m2);
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                          float log_m);
double radiation_get_ionization_rate_from_integral(const struct radiation *rad,
                                                   float log_m1, float log_m2);
double radiation_get_ionization_rate_from_raw(const struct radiation *rad,
                                              float log_m);
double radiation_get_mean_excess_photon_energy_HI_from_integral(
    const struct radiation *rad, float log_m1, float log_m2);
double radiation_get_mean_excess_photon_energy_HI_from_raw(
    const struct radiation *rad, float log_m);

float radiation_get_log_metallicity(float Z);
float radiation_get_luminosities_from_raw_2d(const struct radiation *rad,
                                             float log_z, float log_m);
double radiation_get_ionization_rate_from_raw_2d(const struct radiation *rad,
                                                 float log_z, float log_m,
                                                 float star_age_myr);
double radiation_get_mean_excess_photon_energy_HI_from_raw_2d(
    const struct radiation *rad, float log_z, float log_m, float star_age_myr);

float radiation_get_luminosities_from_integral_2d(const struct radiation *rad,
                                                  float log_z, float log_m1,
                                                  float log_m2);
double radiation_get_ionization_rate_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2);
double radiation_get_mean_excess_photon_energy_HI_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2);
float radiation_get_main_sequence_lifetime_inverse_mass_2d(
    const struct radiation *rad, float log_z, float star_age_myr, float m_min);

float radiation_get_star_luminosity(const struct radiation *rad, float log_m,
                                    float log_z);
double radiation_get_star_ionization_rate(const struct radiation *rad,
                                          float log_m, float log_z,
                                          float star_age_myr);
double radiation_get_star_mean_excess_photon_energy_HI(
    const struct radiation *rad, float log_m, float log_z, float star_age_myr);

float radiation_get_teff_from_raw(const struct radiation *rad, float log_m);
float radiation_get_teff_from_raw_2d(const struct radiation *rad, float log_z,
                                     float log_m);
float radiation_get_star_teff(const struct radiation *rad, float log_m,
                              float log_z);

float radiation_get_luminosity_pe_from_raw(const struct radiation *rad,
                                           float log_m);
float radiation_get_luminosity_pe_from_raw_2d(const struct radiation *rad,
                                              float log_z, float log_m);
float radiation_get_star_luminosity_pe(const struct radiation *rad, float log_m,
                                       float log_z);
float radiation_get_luminosity_lw_from_raw(const struct radiation *rad,
                                           float log_m);
float radiation_get_luminosity_lw_from_raw_2d(const struct radiation *rad,
                                              float log_z, float log_m);
double radiation_get_mean_photon_energy_lw_from_raw(const struct radiation *rad,
                                                    float log_m);
double radiation_get_mean_photon_energy_lw_from_raw_2d(
    const struct radiation *rad, float log_z, float log_m);
double radiation_get_star_mean_photon_energy_lw(const struct radiation *rad,
                                                float log_m, float log_z);
double radiation_get_mean_photon_energy_lw_from_integral(
    const struct radiation *rad, float log_m);
double radiation_get_mean_photon_energy_lw_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m);
float radiation_get_star_luminosity_lw(const struct radiation *rad, float log_m,
                                       float log_z);

float radiation_get_luminosity_pe_from_integral(const struct radiation *rad,
                                                float log_m1, float log_m2);
float radiation_get_luminosity_pe_from_integral_2d(const struct radiation *rad,
                                                   float log_z, float log_m1,
                                                   float log_m2);
float radiation_get_luminosity_lw_from_integral(const struct radiation *rad,
                                                float log_m1, float log_m2);
float radiation_get_luminosity_lw_from_integral_2d(const struct radiation *rad,
                                                   float log_z, float log_m1,
                                                   float log_m2);

float radiation_get_luminosity_edge_pe_from_integral(
    const struct radiation *rad, float log_m1, float log_m2);
float radiation_get_luminosity_edge_pe_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2);
float radiation_get_luminosity_edge_lw_from_integral(
    const struct radiation *rad, float log_m1, float log_m2);
float radiation_get_luminosity_edge_lw_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2);

void radiation_read_data(struct radiation *rad, struct swift_params *params,
                         const struct stellar_model *sm,
                         const struct unit_system *us,
                         const struct phys_const *phys_const,
                         const int restart);
void radiation_read_luminosities_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_ionization_rate_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_mean_excess_photon_energy_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_teff_array(struct radiation *rad, hid_t group_id,
                               const struct radiation_grid_metadata *grid,
                               const struct stellar_model *sm,
                               const struct unit_system *us);
void radiation_read_luminosity_pe_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_mean_photon_energy_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm);
void radiation_read_luminosity_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_luminosity_edge_pe_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_luminosity_edge_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_main_sequence_lifetime_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_main_sequence_lifetime_inverse_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);

#endif /* SWIFT_RADIATION_GEAR_H */
