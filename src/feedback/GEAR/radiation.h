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
 * @brief Subgrid radiation feedback for GEAR. This files contains functions to
 * compute quantities for the radiation feedback.
 */

#include "../../feedback_properties.h"
#include "cooling_properties.h"
#include "hdf5_functions.h"
#include "hydro.h"
#include "part.h"
#include "physical_constants.h"
#include "stellar_evolution_struct.h"
#include "units.h"

/*! Rescaling applied when building/reading the ionizing-photon-rate
    (dot_N_ion) interpolation tables: raw values are ~1e48-1e50 photons/s in
    code units, too large for the tables' float storage, so the table holds
    (raw value / this factor) and every reader multiplies back by it. */
#define RADIATION_DOT_N_ION_TABLE_SCALING 1e50

/*! Lifetime granted to an ionization tag, in units of the rebuild interval
    that produced it. Must exceed 1: a tag has to outlive the gap to its own
    star's next rebuild pass, or cooling (which expires tags on the gas
    particle's own, independently-binned timestep; see
    src/cooling/grackle/cooling_gear_subgrid.h) can clear it before the star
    ever gets the chance to renew it. */
#define RADIATION_TAG_LIFETIME_INTERVALS 2.0

/*! Ceiling on the elapsed interval the per-pass photon budget is integrated
    over, in units of the rebuild cadence actually in force. A scheduled pass
    is skipped whenever the star's working-level cell holds no gas, which
    leaves the last-rebuild stamp untouched and lets the elapsed interval grow
    without bound; the photons emitted meanwhile escaped an empty cell rather
    than being stored, so they must not be handed to the next landing pass. */
#define HII_DT_BACK_MAX_INTERVALS 2.0

/*! Floor applied before taking log10() of a Data/Radiation CGS value, so a
    genuinely-zero table entry (Q_H/DotEExcess below the source table's own
    ionization-threshold mass, where pychem defines them to be exactly 0;
    see PyChemInitTable/libradiation.py) does not produce log10(0) = -inf.
    Matches pychem's own floor bit-for-bit (PyChemInitTable/
    libparsec_radiation.py's `_LOG_FLOOR`): SWIFT's radiation interpolation
    must match pychem's own log10(mass)-vs-log10(value) interpolation
    scheme. Applied to the CGS value itself (before the internal-unit
    conversion), not the converted value: this is the number pychem's own
    floor actually clamps. So, a query that is "native-zero" in SWIFT means
    exactly what it means in pychem, independently of SWIFT's unit
    system. */
#define RADIATION_LOG_FLOOR_CGS 1e-300

/*! Relative epsilon a 2D IMF-integrated getter's query mass is nudged below
    the integrated table's own top mass edge before calling interpolate_2d(),
    so an exact-mass_max query deterministically takes the blended (not
    boundary-clamped) branch. See radiation_get_luminosities_from_
    integral_2d()'s own doxygen for why this matters. */
#define RADIATION_2D_EDGE_EPS 1e-5f

/**
 * @brief Transient, read-time-only grid metadata shared by every dataset in
 * a Data/Radiation HDF5 group. Not part of the persistent #radiation
 * struct: rebuilt fresh by radiation_read_data() on every read, including
 * on restart.
 */
struct radiation_grid_metadata {
  /*! "M" (mass-only) or "M,Z" (mass x metallicity), from the group's own
      "dimensionality" attribute. */
  char dimensionality[8];

  /*! Is this a 2D ("M,Z") table? Derived from #dimensionality. */
  int is_2d;

  /*! log10(mass grid minimum), from the group's "m0" attribute. */
  float log_mass_min;

  /*! log10 mass grid step, from the group's "dm" attribute. */
  float mass_step;

  /*! Number of mass grid points, from the group's "nm" attribute. */
  int n_mass;

  /*! Number of metallicity grid points (0 for a 1D table), from the
      group's "nz" attribute. */
  int n_metallicity;

  /*! Metallicity grid values (mass fraction Z, native units; NULL for a
      1D table). Not guaranteed log-uniformly spaced by the file. See
      radiation_read_data()'s own comment on the approximation this forces
      for interpolate_2d_init(). */
  float *metallicity;

  /*! Mass-axis boundary condition for the "Luminosity" dataset (2D tables
      only; boundary_condition_error otherwise), from the group's
      edge_policy_luminosity_below/above attributes. The metallicity axis
      always clamps (boundary_condition_const), matching pychem's own
      "clamp to nearest grid Z, never extrapolated" convention; see
      radiation_build_tables(). */
  enum interpolate_boundary_condition edge_policy_luminosity;

  /*! Mass-axis boundary condition for the "Q_H" dataset (2D tables only).
      Dispatched on the group's "source" attribute between
      edge_policy_q_h_blackbody_* and edge_policy_q_h_parsec_*, matching
      which Q_H variant pychem wrote as the table's primary Q_H. */
  enum interpolate_boundary_condition edge_policy_q_h;

  /*! Mass-axis boundary condition for the "DotEExcess" dataset (2D tables
      only). Always edge_policy_q_h_blackbody_*, regardless of "source":
      pychem computes DotEExcess as Q_H_Blackbody * MeanExcessPhotonEnergyHI
      unconditionally. */
  enum interpolate_boundary_condition edge_policy_dot_e_excess;
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

double radiation_get_part_mean_molecular_weight(
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
long long radiation_get_part_ionized_star_id(const struct part *p,
                                             const struct xpart *xpj);
float radiation_get_part_excess_photon_energy_HI(const struct part *p,
                                                 const struct xpart *xpj);
float radiation_get_part_photoionization_rate_coefficient(
    const struct part *p, const struct xpart *xpj);
double radiation_get_photoionization_rate_coefficient_from_flux_HI(
    const struct unit_system *us, const double ionizing_flux_HI);
void radiation_set_ionizing_photon_rate(struct spart *sp,
                                        double dot_N_ion_total,
                                        int n_HII_pixels);
void radiation_zero_spart_output(struct spart *sp);
void radiation_open_ionizing_photon_budget(struct spart *sp, double dt_back);
void radiation_consume_ionizing_photons(struct spart *sp, int pixel,
                                        double Delta_N_ion);
float radiation_get_comoving_gas_column_density_at_star(const struct spart *sp);

float radiation_get_physical_IR_opacity(const struct spart *sp,
                                        const struct unit_system *us);

float radiation_get_physical_IR_optical_depth(const struct spart *sp,
                                              const struct unit_system *us,
                                              const struct cosmology *cosmo);

float radiation_get_star_physical_radiation_pressure(
    const struct spart *sp, const float Delta_t,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo);

/******************************************************************************/
/* Functions to deal with integrated data over an IMF. These functions read,
   interpolate and integrate. */
/******************************************************************************/
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
float radiation_get_ms_lifetime_inverse_mass_2d(const struct radiation *rad,
                                                float log_z, float star_age_myr,
                                                float m_min);

float radiation_get_star_luminosity(const struct radiation *rad, float log_m,
                                    float log_z);
double radiation_get_star_ionization_rate(const struct radiation *rad,
                                          float log_m, float log_z,
                                          float star_age_myr);
double radiation_get_star_mean_excess_photon_energy_HI(
    const struct radiation *rad, float log_m, float log_z, float star_age_myr);

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
void radiation_read_main_sequence_lifetime_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);
void radiation_read_main_sequence_lifetime_inverse_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us);

#endif /* SWIFT_RADIATION_GEAR_H */
