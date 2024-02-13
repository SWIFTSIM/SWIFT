/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COSMOLOGY_H
#define SWIFT_COSMOLOGY_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "parser.h"
#include "physical_constants.h"
#include "timeline.h"
#include "units.h"

/**
 * @brief Cosmological parameters
 */
struct cosmology {

  /*! Current expansion factor of the Universe */
  double a;

  /*! Inverse of the current expansion factor of the Universe */
  double a_inv;

  /*! Inverse square of the current expansion factor of the Universe */
  double a2_inv;

  /*! Inverse cube of the current expansion factor of the Universe */
  double a3_inv;

  /*! Power of the scale-factor used for internal energy conversion to physical
   */
  double a_factor_internal_energy;

  /*! Power of the scale-factor used for pressure conversion to physical */
  double a_factor_pressure;

  /*! Power of the scale-factor used for sound-speed conversion to physical */
  double a_factor_sound_speed;

  /*! Power of the scale-factor used for relative velocities in visc. terms */
  double a_factor_mu;

  /*! Power of the scale-factor used for epsilon term in the Balsara switch */
  double a_factor_Balsara_eps;

  /*! Power of the scale-factor used for gravity accelerations */
  double a_factor_grav_accel;

  /*! Power of the scale-factor used for hydro accelerations */
  double a_factor_hydro_accel;

  /*! Current redshift */
  double z;

  /*! Hubble constant at the current redshift (in internal units) */
  double H;

  /*! The critical density at the current redshift (in internal physical units)
   */
  double critical_density;

  /*! The critical density at redshift 0 (in internal physical units) */
  double critical_density_0;

  /*! The mean density at the current redshift (in internal physical units) */
  double mean_density;

  /*! The mean matter density at the current redshift (in internal physical
   * units) */
  double mean_density_Omega_m;

  /*! The mean baryonic density at the current redshift (in internal physical
   * units) */
  double mean_density_Omega_b;

  /*! Over-density for virialised haloes at the current redshift
   * from the Bryan & Norman 1998 fit */
  double overdensity_BN98;

  /*! Conversion factor from internal time-step size to cosmological step */
  double time_step_factor;

  /*! Expansion rate at the current redshift (in internal units) */
  double a_dot;

  /*! Time (in internal units) since the Big Bang */
  double time;

  /*! Lookback time (in internal units) */
  double lookback_time;

  /*! Dark-energy equation of state at the current time */
  double w;

  /*! Scale-factor at the previous time-step */
  double a_old;

  /*! Redshit at the previous time-step */
  double z_old;

  /*------------------------------------------------------------------ */

  /*! Starting expansion factor */
  double a_begin;

  /*! Final expansion factor */
  double a_end;

  /*! Time (in internal units) since the Big Bang at the start */
  double time_begin;

  /*! Time (in internal units) since the Big Bang at the end */
  double time_end;

  /*! Conversion factor from integer time-line to \f$ d\log{a} \f$ */
  double time_base;

  /*! Inverse of conversion factor from integer time-line to \f$ d\log{a} \f$ */
  double time_base_inv;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc)) */
  double h;

  /*! Hubble constant at z = 0 (in internal units) */
  double H0;

  /*! Hubble time 1/H0 */
  double Hubble_time;

  /*! Cold Dark Matter density parameter */
  double Omega_cdm;

  /*! Baryon density parameter */
  double Omega_b;

  /*! Cosmological constant density parameter */
  double Omega_lambda;

  /*! Total radiation density parameter (photons and other relics) */
  double Omega_r;

  /*! CMB radiation density parameter (Omega_gamma) */
  double Omega_g;

  /*! Massive neutrino density parameter */
  double Omega_nu;

  /*! Massive neutrino density parameter at z=0 */
  double Omega_nu_0;

  /*! Ultra-relativistic species (e.g. massless neutrinos) density parameter */
  double Omega_ur;

  /*! Curvature density parameter */
  double Omega_k;

  /*! Dark-energy equation of state at z=0 */
  double w_0;

  /*! Dark-energy evolution parameter */
  double w_a;

  /*! CMB temperature at z = 0 implied by cosmology (internal units) */
  double T_CMB_0;

  /*! CMB temperature at z = 0 implied by cosmology (Kelvins) */
  double T_CMB_0_K;

  /*! Neutrino temperature at z = 0 (internal units) */
  double T_nu_0;

  /* Neutrino temperature at z = 0 (electron-volts) */
  double T_nu_0_eV;

  /*! Number of massive neutrino species */
  int N_nu;

  /*! Number of ultra-relativistic species (excluding massive neutrinos) */
  double N_ur;

  /*! Effective number of relativistic species (including massive neutrinos) */
  double N_eff;

  /*! Mass of each massive neutrino species in electron-volts */
  double *M_nu_eV;

  /*! Degeneracy of each massive neutrino species */
  double *deg_nu;

  /*! Sum of massive neutrino degeneracies */
  double deg_nu_tot;

  /*! Log of starting expansion factor for neutrino interpolation tables */
  double log_a_long_begin;

  /*! Log of midpoint expansion factor for neutrino interpolation tables */
  double log_a_long_mid;

  /*! Log of ending expansion factor for neutrino interpolation tables */
  double log_a_long_end;

  /*! Log of starting expansion factor */
  double log_a_begin;

  /*! Log of final expansion factor */
  double log_a_end;

  /*! Speed of light (internal units) */
  double const_speed_light_c;

  /*! Drift factor interpolation table */
  double *drift_fac_interp_table;

  /*! Kick factor (gravity) interpolation table */
  double *grav_kick_fac_interp_table;

  /*! Kick factor (hydro) interpolation table */
  double *hydro_kick_fac_interp_table;

  /*! Kick factor (hydro correction) interpolation table (GIZMO-MFV only) */
  double *hydro_kick_corr_interp_table;

  /*! Time interpolation table */
  double *time_interp_table;

  /*! Scale factor interpolation table */
  double *scale_factor_interp_table;

  /*! Comoving distance interpolation table */
  double *comoving_distance_interp_table;

  /*! Comoving distance from present day (a=1) to a_end */
  double comoving_distance_interp_table_offset;

  /*! Comoving distance from a_start to a_end */
  double comoving_distance_start_to_end;

  /*! Comoving distance inverse interpolation table */
  double *comoving_distance_inverse_interp_table;

  /*! Massive neutrino density interpolation table at early times */
  double *neutrino_density_early_table;

  /*! Massive neutrino density interpolation table at late times */
  double *neutrino_density_late_table;

  /*! Time between Big Bang and first entry in the table */
  double time_interp_table_offset;

  /*! Time between Big Bang and last entry in the table */
  double time_interp_table_max;

  /*! Time at the present-day (a=1) */
  double universe_age_at_present_day;
};

void cosmology_update(struct cosmology *c, const struct phys_const *phys_const,
                      integertime_t ti_current);

double cosmology_get_drift_factor(const struct cosmology *cosmo,
                                  const integertime_t ti_start,
                                  const integertime_t ti_end);
double cosmology_get_grav_kick_factor(const struct cosmology *cosmo,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end);
double cosmology_get_hydro_kick_factor(const struct cosmology *cosmo,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end);
double cosmology_get_therm_kick_factor(const struct cosmology *cosmo,
                                       const integertime_t ti_start,
                                       const integertime_t ti_end);
double cosmology_get_corr_kick_factor(const struct cosmology *cosmo,
                                      const integertime_t ti_start,
                                      const integertime_t ti_end);
double cosmology_get_delta_time(const struct cosmology *c,
                                const integertime_t ti_start,
                                const integertime_t ti_end);
double cosmology_get_neutrino_density(const struct cosmology *c, double a);

double cosmology_get_delta_time_from_scale_factors(const struct cosmology *c,
                                                   const double a_start,
                                                   const double a_end);

double cosmology_get_timebase(struct cosmology *c,
                              const integertime_t ti_current);

double cosmology_get_scale_factor(const struct cosmology *cosmo, double t);

double cosmology_get_comoving_distance(const struct cosmology *c,
                                       const double a);

double cosmology_scale_factor_at_comoving_distance(const struct cosmology *c,
                                                   double r);

double cosmology_get_time_since_big_bang(const struct cosmology *c, double a);
void cosmology_init(struct swift_params *params, const struct unit_system *us,
                    const struct phys_const *phys_const, struct cosmology *c);

void cosmology_init_no_cosmo(struct cosmology *c);

void cosmology_print(const struct cosmology *c);
void cosmology_clean(struct cosmology *c);

#ifdef HAVE_HDF5
void cosmology_write_model(hid_t h_grp, const struct cosmology *c);
#endif

/* Dump/restore. */
void cosmology_struct_dump(const struct cosmology *cosmology, FILE *stream);
void cosmology_struct_restore(int enabled, struct cosmology *cosmology,
                              FILE *stream);

#endif /* SWIFT_COSMOLOGY_H */
