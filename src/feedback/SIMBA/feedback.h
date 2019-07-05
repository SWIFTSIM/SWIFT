/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_SIMBA_H
#define SWIFT_FEEDBACK_SIMBA_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"


/**
 * @brief Calculates speed particles will be kicked based on
 * host galaxy properties 
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_kick_speed(struct spart *sp, const struct feedback_props *feedback_props, const struct cosmology *cosmo) {

  /* Calculate circular velocity based on Baryonic Tully-Fisher relation*/
  const float v_circ = pow(sp->feedback_data.host_galaxy_mass/feedback_props->simba_host_galaxy_mass_norm, feedback_props->simba_v_circ_exp);

  /* checkout what this random number does and how to generate it */
  const float random_num = 1.;

  /* Calculate wind speed */
  // ALEXEI: checkout what the numbers in this equation mean. Maybe possible future simplifications
  sp->feedback_data.to_distribute.v_kick = feedback_props->galsf_firevel 
      * pow(v_circ * cosmo->a /feedback_props->scale_factor_norm,feedback_props->galsf_firevel_slope) 
      * pow(feedback_props->scale_factor_norm,0.12 - feedback_props->galsf_firevel_slope)
      * (1. - feedback_props->vwvf_scatter - 2.*feedback_props->vwvf_scatter*random_num)
      * v_circ;

  // ALEXEI: temporarily set to arbitrary number for testing.
  sp->feedback_data.to_distribute.v_kick = 500.;
  
}

/**
 * @brief Calculates speed particles will be kicked based on
 * host galaxy properties 
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_mass_loading(struct spart *sp, const struct feedback_props *feedback_props) {
  // ALEXEI: check units of host_galaxy_mass with Romeel: is this a mass or mass/msun?
  if (sp->mass < feedback_props->simba_mass_spectrum_break) {
    sp->feedback_data.to_distribute.wind_mass = feedback_props->simba_wind_mass_eta 
        * sp->feedback_data.host_galaxy_mass * pow(sp->mass/feedback_props->simba_mass_spectrum_break,feedback_props->simba_low_mass_power);
  } else {
    sp->feedback_data.to_distribute.wind_mass = feedback_props->simba_wind_mass_eta 
        * sp->feedback_data.host_galaxy_mass * pow(sp->mass/feedback_props->simba_mass_spectrum_break,feedback_props->simba_high_mass_power);
  }
};

/**
 * @brief Calculates amount of extra thermal energy injection required to
 * make up difference between energy injected as wind and total energy
 * injected due to SN.  
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_heating(struct spart *sp, const struct feedback_props *feedback_props) {
  /* Calculate the amount of energy injected in the wind. 
   * This is in terms of internal energy because we don't 
   * know the mass of the particle being kicked yet. */
  float u_wind = 0.5*sp->feedback_data.to_distribute.v_kick*sp->feedback_data.to_distribute.v_kick;

  /* Calculate internal energy contribution from SN */
  float u_SN = feedback_props->SN_energy * sp->feedback_data.host_galaxy_mass 
               / sp->feedback_data.to_distribute.wind_mass;

  // ALEXEI: should this be smoothed metal mass fraction?
  if (sp->chemistry_data.metal_mass_fraction[0] < 1.e-9) {
    u_SN *= exp10(-0.0029*pow(log10(sp->chemistry_data.metal_mass_fraction[0])+9,2.5)+0.417694); // ALEXEI: what are all these numbers?
  } else {
    u_SN *= 2.61634;
  }

  if (u_wind > feedback_props->simba_wind_energy_limit) 
    sp->feedback_data.to_distribute.v_kick *= sqrt(feedback_props->simba_wind_energy_limit*u_SN/u_wind);
  if (feedback_props->simba_wind_energy_limit < 1.f) 
    u_SN *= feedback_props->simba_wind_energy_limit;

  /* Now we can decide if there's any energy left over to distribute */
  sp->feedback_data.to_distribute.u_extra = max(u_SN - u_wind, 0.);
};

inline void compute_feedback_probability(struct spart *sp) {
  const float mass_frac = sp->feedback_data.to_distribute.wind_mass/sp->mass;
  sp->feedback_data.to_distribute.feedback_probability = 1. - exp(-mass_frac);
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {
  // Temporarily set the particle's galaxy host mass artificially.
  sp->feedback_data.host_galaxy_mass = 1.;

}

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 */
__attribute__((always_inline)) INLINE static int feedback_do_feedback(
    const struct spart* sp) {

  return (sp->birth_time != -1.);
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param time The current simulation time (Non-cosmological runs).
 * @param cosmo The cosmological model (cosmological runs).
 * @param with_cosmology Are we doing a cosmological run?
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const float time, const struct cosmology* cosmo,
    const int with_cosmology) {

  return 1;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {
  sp->feedback_data.to_distribute.simba_delay_time = feedback_props->simba_delay_time;
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
__attribute__((always_inline)) INLINE static void feedback_evolve_spart(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const double star_age_beg_step, const double dt) {
  
  /* Calculate the velocity to kick neighbouring particles with */
  compute_kick_speed(sp, feedback_props, cosmo);

  /* Compute wind mass loading */
  compute_mass_loading(sp, feedback_props);

  /* Compute residual heating */
  compute_heating(sp, feedback_props);

  /* Compute probability of doing feedback */
  compute_feedback_probability(sp);

}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
static INLINE void feedback_struct_dump(const struct feedback_props* feedback,
                                        FILE* stream) {}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
static INLINE void feedback_struct_restore(struct feedback_props* feedback,
                                           FILE* stream) {}

#endif /* SWIFT_FEEDBACK_SIMBA_H */
