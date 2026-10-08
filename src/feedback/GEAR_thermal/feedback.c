/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "feedback.h"

/* Local includes */
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "feedback_tracers.h"
#include "hydro_properties.h"
#include "minmax.h"
#include "part.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part *p, struct xpart *xp,
                          const struct engine *e) {

  /* Did the particle receive a supernovae */
  if (!xp->feedback_data.hit_by_SN && !xp->feedback_data.hit_by_winds) return;

  const struct cosmology *cosmo = e->cosmology;
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;

  /* Turn off the cooling only when SN are involved */
  if (xp->feedback_data.hit_by_SN) {
    cooling_set_part_time_cooling_off(p, xp, e->time);
  }

  /* Update mass */
  const float old_mass = hydro_get_mass(p);
  const float new_mass = old_mass + xp->feedback_data.delta_mass;

  if (xp->feedback_data.delta_mass < 0.) {
    error("Delta mass smaller than 0");
  }

  /* Correction when several events reach the particle: the directed momentum
     is scaled down and the thermal energy takes the difference between the
     kinetic energy the events gave alone and the one the particle gets. It
     needs the velocity before the update and changes delta_p and delta_E_th
     in place. */
  float f_corr = 1.0f;
  float u_residual = 0.0f;
  if (e->feedback_props->enable_multiple_SN_momentum_correction_factor &&
      (xp->feedback_data.number_SN + xp->feedback_data.number_winds > 1)) {
    f_corr = feedback_compute_momentum_correction_factor_for_multiple_sn_events(
        p, xp, cosmo);
    u_residual =
        feedback_compute_residual_internal_energy_for_multiple_sn_events(
            p, xp, cosmo, old_mass, new_mass, f_corr);
    for (int i = 0; i < 3; i++)
      xp->feedback_data.delta_p[i] -=
          (1.0f - f_corr) * xp->feedback_data.delta_p_directed[i];
    xp->feedback_data.delta_E_th += u_residual * new_mass;
  }

  hydro_set_mass(p, new_mass);

  xp->feedback_data.delta_mass = 0;

  /* Update the density */
  p->rho *= new_mass / old_mass;

  /* Update internal energy (the final mass shares the energy of all events) */
  const float new_mass_inv = 1.0f / new_mass;
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo) * old_mass *
                  new_mass_inv;
  const float u_new = u + xp->feedback_data.delta_E_th * new_mass_inv;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

  xp->feedback_data.delta_E_th = 0.0f;

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = xp->feedback_data.delta_p[i] * new_mass_inv;

    xp->v_full[i] += dv;
    p->v[i] += dv;

    xp->feedback_data.delta_p[i] = 0;
  }

  /* Update chemistry metal mass and reset feedback variable */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.metal_mass[i] += xp->feedback_data.delta_metal_mass[i];
    xp->feedback_data.delta_metal_mass[i] = 0.0;
  }

  /* The tracers get what the particle received, with the final mass */
  feedback_tracers_pending_update(xp, xp->feedback_data.hit_by_SN,
                                  xp->feedback_data.hit_by_winds, f_corr,
                                  u_residual, new_mass_inv);
  feedback_tracers_pending_reset(xp);

  for (int i = 0; i < 3; i++) xp->feedback_data.delta_p_directed[i] = 0.0f;
  xp->feedback_data.delta_p_norm_2_sum = 0.0f;
  xp->feedback_data.delta_E_kin_events = 0.0f;
  xp->feedback_data.delta_p_hubble_work = 0.0f;
  xp->feedback_data.number_SN = 0;
  xp->feedback_data.number_winds = 0;

  xp->feedback_data.hit_by_SN = 0;
  xp->feedback_data.hit_by_winds = 0;
}

/**
 * @brief Compute the factor that scales the directed momentum down when
 * several feedback events reach one #part in a timestep.
 *
 * It is the factor that makes the kinetic energy of the summed directed
 * momenta equal to the sum of the kinetic energies of each, never above 1
 * (Okamoto et al., https://arxiv.org/pdf/2603.17421). Only the momentum that
 * the winds direct away from the stars enters; the ejecta keep their momentum.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param cosmo The #cosmology.
 */
float feedback_compute_momentum_correction_factor_for_multiple_sn_events(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {

  const float dp[3] = {xp->feedback_data.delta_p_directed[0] * cosmo->a_inv,
                       xp->feedback_data.delta_p_directed[1] * cosmo->a_inv,
                       xp->feedback_data.delta_p_directed[2] * cosmo->a_inv};
  const float dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

  /* No directed momentum, or the events cancelled each other */
  if (dp_norm_2 <= 0.f) return 1.f;

  const float f_corr = sqrtf(xp->feedback_data.delta_p_norm_2_sum / dp_norm_2);
  return (f_corr >= 1.f) ? 1.f : f_corr;
}

/**
 * @brief Add the kinetic energy that one feedback event would give to the
 * #part alone, for the energy balance of overlapping events.
 *
 * The energy is the one in the peculiar frame of the gas: the event changes
 * the momentum of the gas of mass mj and velocity v by dp + dp_ejecta and its
 * mass to new_mass. The work of the directed momentum against the Hubble flow
 * relative to the star is stored apart, since the correction rescales it.
 *
 * @param xp The #xpart.
 * @param mj The mass of the gas particle before the events.
 * @param new_mass The mass of the gas particle after this event.
 * @param v_pec The physical peculiar velocity of the gas particle.
 * @param v_hubble The physical Hubble flow velocity of the gas particle
 * relative to the star.
 * @param dp The physical directed momentum given by the event.
 * @param dp_ejecta The physical momentum of the ejecta in the lab frame.
 */
void feedback_accumulate_kinetic_energy_for_multiple_sn_events(
    struct xpart *xp, const float mj, const float new_mass,
    const float v_pec[3], const float v_hubble[3], const double dp[3],
    const double dp_ejecta[3]) {

  const double dp_tot[3] = {dp[0] + dp_ejecta[0], dp[1] + dp_ejecta[1],
                            dp[2] + dp_ejecta[2]};
  const double v_norm_2 =
      v_pec[0] * v_pec[0] + v_pec[1] * v_pec[1] + v_pec[2] * v_pec[2];
  const double v_dot_dp =
      v_pec[0] * dp_tot[0] + v_pec[1] * dp_tot[1] + v_pec[2] * dp_tot[2];
  const double dp_tot_norm_2 =
      dp_tot[0] * dp_tot[0] + dp_tot[1] * dp_tot[1] + dp_tot[2] * dp_tot[2];

  /* A difference of the energies, written without two large terms */
  const double dE_kin =
      0.5 *
      (-mj * v_norm_2 * (new_mass - mj) + 2.0 * mj * v_dot_dp + dp_tot_norm_2) /
      new_mass;

  xp->feedback_data.delta_E_kin_events += dE_kin;
  xp->feedback_data.delta_p_hubble_work +=
      v_hubble[0] * dp[0] + v_hubble[1] * dp[1] + v_hubble[2] * dp[2];
}

/**
 * @brief Store what one stellar wind event gives to a #part, for the
 * multiple-event correction.
 *
 * Only the momentum that the wind directs away from the star is rescaled by
 * the correction; the ejecta keep their momentum. The values are recomputed
 * from the inputs of runner_iact_nonsym_feedback_apply(), and the function is
 * not inlined, so the injection itself compiles as without the correction.
 *
 * @param xp The #xpart.
 * @param si The star.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param r2 Comoving square distance between the two particles.
 * @param weight The share of the star ejecta given to the #part.
 * @param mj The mass of the gas particle before the events.
 * @param dm_SW The wind mass given to the #part.
 * @param new_mass The mass of the gas particle after this event.
 * @param cosmo The #cosmology.
 */
__attribute__((noinline)) void feedback_accumulate_wind_for_multiple_sn_events(
    struct xpart *xp, const struct spart *si, const float dx[3], const float r2,
    const double weight, const float mj, const double dm_SW,
    const double new_mass, const struct cosmology *cosmo) {

  const float a = cosmo->a;
  const float a_inv = cosmo->a_inv;
  const float a_dot = a * cosmo->H;
  const float r_p = sqrtf(r2) * a;
  const float p_ej = sqrt(2.0 * si->feedback_data.winds.mass_ejected *
                          si->feedback_data.winds.energy_ejected);

  float v_pec[3], v_hubble[3];
  double dp[3], dp_ejecta[3];
  double dp_norm_2 = 0.0;
  for (int i = 0; i < 3; i++) {
    v_pec[i] = xp->v_full[i] * a_inv;
    v_hubble[i] = -a_dot * dx[i];
    dp[i] = -weight * p_ej * dx[i] * a / r_p;
    dp_ejecta[i] = dm_SW * si->v[i] * a_inv;
    dp_norm_2 += dp[i] * dp[i];
    xp->feedback_data.delta_p_directed[i] += dp[i] * a;
  }

  feedback_accumulate_kinetic_energy_for_multiple_sn_events(
      xp, mj, new_mass, v_pec, v_hubble, dp, dp_ejecta);
  xp->feedback_data.delta_p_norm_2_sum += dp_norm_2;
  xp->feedback_data.number_winds += 1;
}

/**
 * @brief Store what one supernova event gives to a #part, for the
 * multiple-event correction.
 *
 * The ejecta carry no directed momentum. Ejecta without energy are not an
 * event, but their kinetic energy enters the balance. Not inlined, like
 * feedback_accumulate_wind_for_multiple_sn_events().
 *
 * @param xp The #xpart.
 * @param si The star.
 * @param mj The mass of the gas particle before the events.
 * @param dm_SN The supernova mass given to the #part.
 * @param new_mass The mass of the gas particle after this event.
 * @param cosmo The #cosmology.
 * @param is_event Does the supernova bring energy?
 */
__attribute__((noinline)) void feedback_accumulate_SN_for_multiple_sn_events(
    struct xpart *xp, const struct spart *si, const float mj,
    const double dm_SN, const double new_mass, const struct cosmology *cosmo,
    const int is_event) {

  const float a_inv = cosmo->a_inv;
  const float v_pec[3] = {xp->v_full[0] * a_inv, xp->v_full[1] * a_inv,
                          xp->v_full[2] * a_inv};
  const float v_zero[3] = {0.f, 0.f, 0.f};
  const double dp_zero[3] = {0.0, 0.0, 0.0};
  const double dp_ejecta[3] = {dm_SN * si->v[0] * a_inv,
                               dm_SN * si->v[1] * a_inv,
                               dm_SN * si->v[2] * a_inv};
  feedback_accumulate_kinetic_energy_for_multiple_sn_events(
      xp, mj, new_mass, v_pec, v_zero, dp_zero, dp_ejecta);
  if (is_event) xp->feedback_data.number_SN += 1;
}

/**
 * @brief Compute the specific internal energy that conserves the energy when
 * several feedback events reach one #part in a timestep.
 *
 * Each event gave its thermal energy as if it were alone. The thermal energy
 * takes the difference between the kinetic energy the events gave alone and
 * the kinetic energy of the summed, rescaled momentum, in the peculiar frame.
 * The correction never removes more thermal energy than the events gave.
 * Called before the velocity update.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param cosmo The #cosmology.
 * @param old_mass The mass of the #part before the events.
 * @param new_mass The mass of the #part after all the events.
 * @param f_corr The momentum correction factor.
 * @return The physical specific internal energy to add.
 */
float feedback_compute_residual_internal_energy_for_multiple_sn_events(
    const struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const float old_mass, const float new_mass, const float f_corr) {

  /* Physical peculiar velocity and momentum in the frame of the gas */
  const double a_inv = cosmo->a_inv;
  double v[3], dp[3];
  for (int i = 0; i < 3; i++) {
    v[i] = xp->v_full[i] * a_inv;
    dp[i] = (xp->feedback_data.delta_p[i] -
             (1.0 - f_corr) * xp->feedback_data.delta_p_directed[i]) *
            a_inv;
  }
  const double v_norm_2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  const double v_dot_dp = v[0] * dp[0] + v[1] * dp[1] + v[2] * dp[2];
  const double dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

  /* Kinetic energy actually given: the gained mass is accelerated to v, then
     the whole particle by dp / new_mass */
  const double dE_kin_actual = 0.5 * (new_mass - old_mass) * v_norm_2 +
                               v_dot_dp + 0.5 * dp_norm_2 / new_mass;

  const double dE_residual =
      xp->feedback_data.delta_E_kin_events +
      (1.0 - f_corr) * xp->feedback_data.delta_p_hubble_work - dE_kin_actual;
  const double u_residual = dE_residual / new_mass;

  /* Do not remove more thermal energy than what the events gave */
  const double u_events = max(xp->feedback_data.delta_E_th, 0.0) / new_mass;
  return (u_residual < -u_events) ? -u_events : u_residual;
}

/**
 * @brief Finishes the #part density calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 */
__attribute__((always_inline)) INLINE void feedback_end_density(
    struct part *p, struct xpart *xp) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here in the GEAR model.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
void feedback_reset_part(struct part *p, struct xpart *xp) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart *sp, const struct engine *e) {

  /* the particle is inactive if its birth_scale_factor or birth_time is
   * negative */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_feedback;
}

/**
 * @brief Should this particle inject anything as supernovae feedback?
 *
 * The supernova branch of runner_iact_nonsym_feedback_apply() runs when the
 * star has supernova energy to give.
 *
 * @param sp The #spart.
 */
int feedback_should_inject_SN_feedback(const struct spart *sp) {
  return sp->feedback_data.supernovae.energy_ejected != 0.f;
}

/**
 * @brief Should this particle inject anything as stellar wind feedback?
 *
 * The wind branch of runner_iact_nonsym_feedback_apply() runs when the star
 * has wind energy to give.
 *
 * @param sp The #spart.
 */
int feedback_should_inject_wind_feedback(const struct spart *sp) {
  return sp->feedback_data.winds.energy_ejected != 0.f;
}

/**
 * @brief Should this particle inject anything as stellar feedback?
 *
 * @param sp The #spart.
 */
int feedback_should_inject_feedback(const struct spart *sp) {
  return feedback_should_inject_SN_feedback(sp) ||
         feedback_should_inject_wind_feedback(sp);
}

/**
 * @brief Get the comoving SPH gas density at the star position.
 *
 * Only valid after feedback_prepare_feedback() has been called for this step.
 *
 * @param sp The #spart.
 */
float feedback_get_comoving_gas_density_at_star(const struct spart *sp) {
  return sp->feedback_data.enrichment_weight;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart *sp) {

  sp->feedback_data.enrichment_weight = 0.f;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 *
 * @param sp The #spart.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props) {}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props) {}

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 *
 * In here, we only need to add the missing coefficients.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  sp->feedback_data.enrichment_weight *= hi_inv_dim;
}
