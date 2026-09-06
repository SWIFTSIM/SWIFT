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
#include "radiation_isrf.h"
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

/*! Lifetime granted to an LW/FUV illumination episode, in units of the
    illuminating star's own integer timestep. Integer, not a double factor
    on a physical time like RADIATION_TAG_LIFETIME_INTERVALS: the
    injection/expiry comparison this feeds (radiation_iact.h,
    radiation_gas.c:radiation_reset_part_LW_FUV_illumination_tag) works
    entirely in integertime_t, so there is no float boundary to round
    against. Must exceed 1 for the same reason as RADIATION_TAG_LIFETIME_
    INTERVALS: a touch has to outlive the gap to the star's next visit, or
    the gas particle's own once-per-step expiry check (feedback_reset_part)
    can clear it first. */
#define RADIATION_LW_FUV_TAG_LIFETIME_INTERVALS 2

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

/*! Non-ionizing FUV band, eV (Habing band, 912-2000 Angstrom). */
#define RADIATION_FUV_BAND_LOW_EV 6.0
#define RADIATION_FUV_BAND_HIGH_EV 11.2

/*! Lyman-Werner band, eV (H2 photodissociation, 912-1108 Angstrom): see
    RADIATION_FUV_BAND_LOW_EV's doxygen. */
#define RADIATION_LW_BAND_LOW_EV 11.2
#define RADIATION_LW_BAND_HIGH_EV 13.6

/*! Boltzmann constant in eV/K, for #radiation_planck_band_fraction's
    dimensionless photon energy x = h*nu / (k_B*T). */
#define RADIATION_BOLTZMANN_K_EV_PER_K 8.617333262e-5

/*! Number of Simpson's-rule sub-intervals #radiation_planck_band_fraction
    integrates the Planck function over; see that function's own doxygen. */
#define RADIATION_PLANCK_QUADRATURE_N 200

/*! Mean mass per hydrogen nucleon (He folded in), for converting a
    per-hydrogen-nucleon dust cross-section (#RADIATION_SIGMA_D_FUV_CGS/
    #RADIATION_SIGMA_D_LW_CGS) into a mass opacity; see
    radiation_get_dust_extinction_factor()'s own doxygen. */
#define RADIATION_MU_H 1.4

/*! Hydrogen atomic mass, g (cgs); see #RADIATION_MU_H. */
#define RADIATION_HYDROGEN_MASS_CGS 1.6726219e-24

/*! Band-specific dust cross-section per hydrogen nucleon, cm^2 (Kim et
    al. 2023, Weingartner & Draine 2001 grain population): 6-11.2 eV
    (FUV) and 11.2-13.6 eV (Lyman-Werner) bands respectively. */
#define RADIATION_SIGMA_D_FUV_CGS 9e-22
#define RADIATION_SIGMA_D_LW_CGS 1.5e-21

/*! Grackle's own solar metal mass fraction, SolarMetalFractionByMass
    (grackle_chemistry_data_fields.def, default 0.01295), not
    radiation_pressure.c's unrelated Z_sun=0.02 (used only for the
    kappa_IR/kappa_NUV fit). By default (chemistry_data.use_dust_density_
    field=0), Grackle computes its own dust-to-gas ratio as
    local_dust_to_gas_ratio * (Z/#RADIATION_GRACKLE_SOLAR_METAL_FRACTION);
    matching that convention keeps our own assumed dust abundance
    consistent with Grackle's dust_chemistry=1-coupled channels for the
    same gas. local_dust_to_gas_ratio only cancels out of our relative
    D(Z)/D(Zsun) scaling when a run leaves it at Grackle's own compiled
    default (#RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO); a run that
    overrides GrackleCooling:local_dust_to_gas_ratio must scale D(Z) by
    the resolved value relative to that default, so it is read
    explicitly at both extinction call sites instead of assumed away. */
#define RADIATION_GRACKLE_SOLAR_METAL_FRACTION 0.01295

/*! Grackle's own compiled default for chemistry_data.local_dust_to_gas_
    ratio (Pollack et al. 1994), resolved from the `-1` sentinel in
    cooling.c when GrackleCooling:local_dust_to_gas_ratio is left unset;
    see #RADIATION_GRACKLE_SOLAR_METAL_FRACTION. */
#define RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO 0.009387

/*! Standard Habing-unit flux normalization, erg/s/cm^2: G0=1 corresponds
    to this flux integrated over the FUV+LW bands. */
#define RADIATION_HABING_FLUX_CGS 1.6e-3

/*! Representative Lyman-Werner photon energy, eV (~12 eV, the band's own
    6-11.2/11.2-13.6 eV midpoint region), for converting an energy flux to
    a photon flux; see radiation_get_part_LW_dissociation_rate_internal()'s
    own doxygen. */
#define RADIATION_LW_PHOTON_ENERGY_EV 12.0

/*! Effective H2 Lyman-Werner-band photodissociation cross section, cm^2:
    an approximation with an implicit assumed spectral shape (the
    band-integrated H2 cross section depends on the spectrum within
    11.2-13.6 eV, not a single atomic-physics constant, so this is a
    single-line-averaged effective value, not re-derivable to arbitrary
    precision). Verified 2026-09-06 against Sternberg, Le Petit, Roueff
    & Le Bourlot (2014, ApJ 790:10): this codebase's own k_diss formula
    (radiation_get_part_LW_dissociation_rate_internal()) compared against
    their free-space photodissociation rate D0=5.8e-11*I_UV (their Eq 5),
    for an equivalent LW-band radiation intensity, gives k_diss/D0=0.898
    (within ~10%, consistent with the ~5-10% uncertainty in converting
    between this codebase's Habing-band convention and Sternberg's
    Draine-normalized I_UV). See
    theory/GEAR/Radiation/verify_sigma_h2_lw_sternberg2014.py for the full
    calculation. */
#define RADIATION_SIGMA_H2_LW_CGS 2.47e-18

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

  /*! Mass-axis boundary condition for the "Q_H" dataset (2D tables only),
      from the group's generic edge_policy_q_h_below/above attributes.
      pychem sets these, at write time, to whichever per-variant edge
      policy (e.g. edge_policy_q_h_blackbody_*, edge_policy_q_h_parsec_*)
      matches the table's own primary Q_H; SWIFT reads them directly and
      never dispatches on the group's "source" attribute itself, so a new
      pychem source mode needs no companion SWIFT change. */
  enum interpolate_boundary_condition edge_policy_q_h;

  /*! Mass-axis boundary condition for the "DotEExcess" dataset (2D tables
      only), from the group's generic
      edge_policy_mean_excess_energy_below/above attributes (mirroring
      #edge_policy_q_h above). DotEExcess's VALUE is (the primary Q_H
      variant "source" selects) times MeanExcessPhotonEnergyHI; this field
      only records the edge policy pychem assigns to
      MeanExcessPhotonEnergyHI/DotEExcess's own primary content, which
      need not match #edge_policy_q_h's variant. */
  enum interpolate_boundary_condition edge_policy_dot_e_excess;

  /*! Mass-axis boundary condition for the "Teff" dataset (2D tables only),
      from the group's generic edge_policy_teff_below/above attributes. */
  enum interpolate_boundary_condition edge_policy_teff;

  /*! Mass-axis boundary condition for the "L_FUV" dataset (2D tables with
      an "L_FUV" dataset only; boundary_condition_error otherwise, matching
      every other edge_policy_* field's convention for a table where the
      corresponding dataset does not apply), from the group's own
      edge_policy_l_fuv_below/above attributes -- a dedicated pair, NOT
      shared with #edge_policy_luminosity, since pychem's L_FUV/L_LW default
      policy ("zero" below the native mass floor, "constant" above) differs
      from Luminosity's own. */
  enum interpolate_boundary_condition edge_policy_l_fuv;

  /*! Mass-axis boundary condition for the "L_LW" dataset (2D tables with an
      "L_LW" dataset only), from the group's own edge_policy_l_lw_below/above
      attributes. See #edge_policy_l_fuv's own doxygen. */
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
void radiation_reset_part_LW_FUV_illumination_tag(struct part *p,
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

float radiation_get_teff_from_raw(const struct radiation *rad, float log_m);
float radiation_get_teff_from_raw_2d(const struct radiation *rad, float log_z,
                                     float log_m);
float radiation_get_star_teff(const struct radiation *rad, float log_m,
                              float log_z);
double radiation_planck_band_fraction(double T_kelvin, double E_low_eV,
                                      double E_high_eV);

float radiation_get_l_fuv_from_raw(const struct radiation *rad, float log_m);
float radiation_get_l_fuv_from_raw_2d(const struct radiation *rad, float log_z,
                                      float log_m);
float radiation_get_star_l_fuv(const struct radiation *rad, float log_m,
                               float log_z);
float radiation_get_l_lw_from_raw(const struct radiation *rad, float log_m);
float radiation_get_l_lw_from_raw_2d(const struct radiation *rad, float log_z,
                                     float log_m);
float radiation_get_star_l_lw(const struct radiation *rad, float log_m,
                              float log_z);

float radiation_get_l_fuv_from_integral(const struct radiation *rad,
                                        float log_m1, float log_m2);
float radiation_get_l_fuv_from_integral_2d(const struct radiation *rad,
                                           float log_z, float log_m1,
                                           float log_m2);
float radiation_get_l_lw_from_integral(const struct radiation *rad,
                                       float log_m1, float log_m2);
float radiation_get_l_lw_from_integral_2d(const struct radiation *rad,
                                          float log_z, float log_m1,
                                          float log_m2);

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
void radiation_read_l_fuv_array(struct radiation *rad, hid_t group_id,
                                const struct radiation_grid_metadata *grid,
                                const struct stellar_model *sm,
                                const struct unit_system *us);
void radiation_read_l_lw_array(struct radiation *rad, hid_t group_id,
                               const struct radiation_grid_metadata *grid,
                               const struct stellar_model *sm,
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
