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
#ifndef SWIFT_RADIATION_IACT_GEAR_H
#define SWIFT_RADIATION_IACT_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_iact.h
 * @brief Subgrid radiation feedback for GEAR: functions called from
 * feedback_iact.h and feedback_prepare_feedback().
 */

#include "chemistry.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "feedback_properties.h"
#include "minmax.h"
#include "radiation.h"
#include "timestep_sync_part.h"
#include "tracers.h"

/**
 * @brief Radiation density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param hydro_props The #hydro_props.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cooling The #cooling_function_data used in the run.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
radiation_iact_nonsym_feedback_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, const struct part *pj, const struct xpart *xpj,
    const struct cosmology *cosmo, const struct feedback_props *fb_props,
    const struct hydro_props *hydro_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current) {

  /* Exit if no radiation policies are enabled */
  if (fb_props->radiation_policy == radiation_policy_none) {
    return;
  }

  const float mj = hydro_get_mass(pj);
  /* Floor avoids a NaN in dx_unit for a coincident pair (r2 == 0). */
  const float r2_min = 1e-6f * hi * hi;
  const float r = sqrtf(max(r2, r2_min));

  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(ui, &wi, &wi_dx);

  /* Unit vector pointing to pj */
  float dx_unit[3];
  for (int k = 0; k < 3; ++k) {
    dx_unit[k] = dx[k] / r;
  }

  float gradW[3];
  for (int k = 0; k < 3; ++k) {
    gradW[k] = wi_dx * dx_unit[k];
  }

  for (int k = 0; k < 3; ++k) {
    si->feedback_data.grad_rho_star[k] += mj * gradW[k];
  }

  /* Same weighting as the gas density at the star, which normalizes it. */
  si->feedback_data.Z_star +=
      chemistry_get_total_metal_mass_fraction_for_feedback(pj) * mj * wi;
}

/**
 * @brief Store the star's feedback time-step for this step.
 *
 * Cached so the per-pair loop avoids the cosmological lookup. @p dt must be
 * the star's own step, as GEAR's feedback_get_enrichment_timestep() returns,
 * in physical time also under cosmology.
 *
 * @param sp The #spart to update.
 * @param dt Length of the star's feedback step, in internal units.
 * @param ti_begin Integer time at the start of that step.
 */
__attribute__((always_inline)) INLINE static void feedback_star_store_timestep(
    struct spart *restrict sp, const double dt, const integertime_t ti_begin) {

  sp->feedback_data.radiation.Delta_t = (float)dt;
#ifdef SWIFT_DEBUG_CHECKS
  sp->feedback_data.radiation.Delta_t_cached_ti_begin = ti_begin;
#else
  (void)ti_begin;
#endif
}

/**
 * @brief Finalize a #spart's radiation-feedback inputs (density gradient,
 * metallicity) and cache its feedback timestep, once per star per step.
 *
 * This is called in the feedback_prepare_feedback(), which is called in the
 * stars ghost task.
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
__attribute__((always_inline)) INLINE static void
feedback_prepare_radiation_feedback(
    struct spart *restrict sp, const struct feedback_props *feedback_props,
    const struct cosmology *cosmo, const struct unit_system *us,
    const struct phys_const *phys_const, const double star_age_beg_step,
    const double dt, const double time, const integertime_t ti_begin,
    const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  sp->feedback_data.grad_rho_star[0] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[1] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[2] *= hi_inv_dim_plus_one;

  /* The gas density is 0 only when the Z_star sum is too, so skipping is
     exact. Z_star needs 1/h^d, not the gradient's 1/h^(d+1). The density is
     already normalized by the caller. */
  const float rho_gas = feedback_get_comoving_gas_density_at_star(sp);
  if (rho_gas > 0.0f) {
    sp->feedback_data.Z_star *= hi_inv_dim / rho_gas;
  }

  feedback_star_store_timestep(sp, dt, ti_begin);
}

/**
 * @brief Radiation feedback of a star on one gas neighbour, for a given
 * share of the star's emission (non-symmetric).
 *
 * Applies radiation pressure and injects the local Lyman-Werner/PE field.
 * The neighbour receives the fraction @p weight of the star's radiation
 * momentum, along -@p dir, and the fraction @p weight of its LW/PE energy.
 * The caller chooses the weights; they must sum to 1 over the neighbours.
 *
 * @param r Comoving distance between the two particles, floored above 0.
 * @param weight Share of the star's emission given to pj.
 * @param dir Vector pointing from pj towards the star (comoving separation
 * or dimensionless weight).
 * @param dir_norm Norm of @p dir, > 0.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param fb_props Properties of the feedback scheme.
 * @param phys_const The physical constants (in internal units).
 * @param us The internal system of units.
 * @param cooling The properties of the cooling scheme.
 * @param ti_current Current integer time
 */
__attribute__((always_inline)) INLINE static void
radiation_iact_nonsym_feedback_apply_weighted(
    const float r, const double weight, const float dir[3],
    const float dir_norm, struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current) {

  const float mj = hydro_get_mass(pj);

  /* Also used to renew the LW/PE illumination window. */
  const integertime_t ti_step = get_integer_timestep(si->time_bin);

  /* Cached once per star by feedback_prepare_radiation_feedback(). */
  const float Delta_t = si->feedback_data.radiation.Delta_t;
#ifdef SWIFT_DEBUG_CHECKS
  if (get_integer_time_begin(ti_current, si->time_bin) !=
      si->feedback_data.radiation.Delta_t_cached_ti_begin)
    error(
        "Stale cached Delta_t: star %lld's step boundary moved since it "
        "was cached.",
        si->id);
#endif

  /* Test the policy bit first: L_bol can be positive from two negative
     factors even with the switch off. */
  if ((fb_props->radiation_policy & radiation_policy_radiation_pressure) &&
      si->feedback_data.radiation.L_bol > 0.0) {
    const float p_rad = radiation_get_star_physical_radiation_pressure(
        si, Delta_t, phys_const, us, cosmo);
    const float delta_p_rad = weight * p_rad;

    /* Along -dir, away from the star; cosmo->a converts to comoving. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.radiation.delta_p[i] -=
          delta_p_rad * dir[i] / dir_norm * cosmo->a;
    }

    /* Tracer. The kick velocity is a magnitude because
       MaxKickVelocityFromRadiationPressure is a running maximum from 0;
       the cumulative momentum carries the sign. */
    tracers_after_radiation_pressure_feedback_part(xpj, delta_p_rad,
                                                   fabsf(delta_p_rad / mj));

    /* Required for feedback_update_part_radiation() to apply the momentum. */
    xpj->feedback_data.hit_by_radiation = 1;
  }

  /* L_band is zero unless with_interstellar_radiation_field is on. */
  if (si->feedback_data.radiation.L_band[ISRF_MOMENT_PE] != 0.0 ||
      si->feedback_data.radiation.L_band[ISRF_MOMENT_LW] != 0.0) {

    const float Z_j = chemistry_get_total_metal_mass_fraction_for_cooling(pj);
    float extinction[ISRF_OPERATOR_COUNT];
    /* Receiver-side: uses pj's own column density, not the source's. */
    const float extinction_path = radiation_get_comoving_extinction_path(
        fb_props, pj, xpj, r, cosmo, phys_const, hydro_props, us, cooling);
    radiation_get_part_ISRF_extinction_factors(
        us, phys_const, cosmo, pj, Z_j, cooling, extinction_path, extinction);

    /* Energy; divided by mj below to get the stored specific energy. */
    double u_inject[ISRF_MOMENT_COUNT];
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
      u_inject[m] = (double)Delta_t * weight *
                    si->feedback_data.radiation.L_band[m] *
                    (double)extinction[radiation_isrf_moment_to_operator[m]];
    }

    if (fb_props->ISRF_propagation) {
      /* Pure accumulation, no reset, so stars on any time bin just add. The
         fold-in happens in radiation_end_force_propagation(). */
      for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
        pj->feedback_data.isrf_moment[m].u_dose_reservoir +=
            (float)(u_inject[m] / (double)mj);
      }
      pj->feedback_data.ISRF_reservoir_end_ti =
          max(pj->feedback_data.ISRF_reservoir_end_ti, ti_current + ti_step);
      pj->feedback_data.ISRF_last_touch_ti = ti_current;
    } else {
      /* Instantaneous field: reset on the first touch this step, then sum
         over the stars. */
      if (pj->feedback_data.ISRF_last_touch_ti != ti_current) {
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
          pj->feedback_data.isrf_moment[m].u = 0.f;
        pj->feedback_data.ISRF_last_touch_ti = ti_current;
      }

      /* #u is double: no narrowing cast. */
      for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
        pj->feedback_data.isrf_moment[m].u += u_inject[m] / (double)mj;
      }
    }

    /* Renew the window on every touch. The expiry check runs once per step
       in feedback_reset_part(). */
    pj->feedback_data.ISRF_illumination_end_ti =
        ti_current + RADIATION_ISRF_TAG_LIFETIME_INTERVALS * ti_step;

    /* Sync on the first touch only: re-syncing every pass drags the whole
       region to the shortest time bin. */
    if (!pj->feedback_data.is_illuminated_ISRF) {
      pj->feedback_data.is_illuminated_ISRF = 1;
      timestep_sync_part(pj);
    }
  }
}

/**
 * @brief Record that a star pass reached a gas neighbour with a zero share of
 * the star's LW/PE emission.
 *
 * Without propagation the field is instantaneous: the neighbour must hold this
 * pass's share (zero), not the value of an earlier pass.
 *
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time
 */
__attribute__((always_inline)) INLINE static void
radiation_iact_nonsym_feedback_apply_zero_share(
    const struct spart *si, struct part *pj,
    const struct feedback_props *fb_props, const integertime_t ti_current) {

  if (fb_props->ISRF_propagation) return;
  if (si->feedback_data.radiation.L_band[ISRF_MOMENT_PE] == 0.0 &&
      si->feedback_data.radiation.L_band[ISRF_MOMENT_LW] == 0.0)
    return;

  if (pj->feedback_data.ISRF_last_touch_ti != ti_current) {
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
      pj->feedback_data.isrf_moment[m].u = 0.f;
    pj->feedback_data.ISRF_last_touch_ti = ti_current;
  }
}

/**
 * @brief Radiation feedback interaction between two particles (non-symmetric),
 * updating the gas particles neighbouring a star particle.
 *
 * Applies radiation pressure and injects the local Lyman-Werner/PE field,
 * shared by the SPH kernel mass weights and directed radially away from the
 * star.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param fb_props Properties of the feedback scheme.
 * @param phys_const The physical constants (in internal units).
 * @param us The internal system of units.
 * @param cooling The properties of the cooling scheme.
 * @param ti_current Current integer time
 */
__attribute__((always_inline)) INLINE static void
radiation_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current) {

  const float mj = hydro_get_mass(pj);
  /* Floor avoids a NaN in the radial kick for a coincident pair (r2 == 0). */
  const float r2_min = 1e-6f * hi * hi;
  const float r = sqrtf(max(r2, r2_min));

  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  const float rho_gas = feedback_get_comoving_gas_density_at_star(si);
  const double si_inv_weight = rho_gas == 0 ? 0. : 1. / rho_gas;
  const double weight = mj * wi * si_inv_weight;

  radiation_iact_nonsym_feedback_apply_weighted(
      r, weight, dx, r, si, pj, xpj, cosmo, hydro_props, fb_props, phys_const,
      us, cooling, ti_current);
}

/**
 * @brief Update the properties of the particle due to radiation feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param mass Mass of the gas the momentum is divided by (GEAR thermal: the
 * mass before any winds or SN).
 */
__attribute__((always_inline)) INLINE static void
feedback_update_part_radiation(struct part *p, struct xpart *xp,
                               const struct engine *e, const float mass) {

  /* Momentum only; internal energy is handled elsewhere, before cooling. */
  if (xp->feedback_data.hit_by_radiation) {
    for (int i = 0; i < 3; i++) {
      const float dv = xp->feedback_data.radiation.delta_p[i] / mass;
      xp->v_full[i] += dv;
      p->v[i] += dv;

      xp->feedback_data.radiation.delta_p[i] = 0;
    }
    xp->feedback_data.hit_by_radiation = 0;
  }
}

#endif /* SWIFT_RADIATION_IACT_GEAR_H */
