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
 * feedback_iact.h and feedback_prepare_feedback(), split out so the
 * mechanical feedback module does not duplicate them.
 */

#include "chemistry.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "feedback_properties.h"
#include "minmax.h"
#include "radiation.h"
#include "timestep_sync_part.h"

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
  const float r = sqrtf(r2);

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

  /* Weighted like enrichment_weight (feedback_iact.h) so both share the
     same kernel normalization in feedback_prepare_radiation_feedback(). */
  si->feedback_data.Z_star +=
      chemistry_get_total_metal_mass_fraction_for_feedback(pj) * mj * wi;
}

/**
 * @brief Store the star's feedback time-step for this step.
 *
 * The pairwise injection loop needs it once per neighbour; caching it here
 * keeps the cosmological lookup off the per-pair path.
 *
 * @p dt already carries the d(ln a) versus proper-time distinction under
 * cosmology, mirroring compute_time() in feedback_common.c. It is the star's
 * own step only because GEAR's feedback_get_enrichment_timestep() returns
 * dt_star unchanged; were that ever to differ, this would silently cache the
 * wrong quantity and the staleness check below would not catch it, because
 * the value would be fresh rather than stale.
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

  /* enrichment_weight is 0 only when Z_star's own sum is too (same
     weighting), so skipping is exact, not an approximation. Uses
     hi_inv_dim, not hi_inv_dim_plus_one: Z_star needs the density
     estimate's 1/h^d, not the gradient's extra 1/h. */
  if (sp->feedback_data.enrichment_weight > 0.0f) {
    sp->feedback_data.Z_star *=
        hi_inv_dim / sp->feedback_data.enrichment_weight;
  }

  feedback_star_store_timestep(sp, dt, ti_begin);
}

/**
 * @brief Radiation feedback interaction between two particles
 * (non-symmetric). Used for updating properties of gas particles neighbouring
 * a star particle.
 *
 * Applies radiation pressure and injects the local Lyman-Werner/PE field.
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
  const float r = sqrtf(r2);

  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  const double si_inv_weight = si->feedback_data.enrichment_weight == 0
                                   ? 0.
                                   : 1. / si->feedback_data.enrichment_weight;
  const double weight = mj * wi * si_inv_weight;

  /* Also reused below to renew the LW/PE illumination window
   * (radiation_reset_part_ISRF_illumination_tag). */
  const integertime_t ti_step = get_integer_timestep(si->time_bin);

  /* Cached once per star by feedback_prepare_radiation_feedback, not
     recomputed per neighbour: shared by radiation pressure and LW/PE
     injection below. */
  const float Delta_t = si->feedback_data.radiation.Delta_t;
#ifdef SWIFT_DEBUG_CHECKS
  if (get_integer_time_begin(ti_current, si->time_bin) !=
      si->feedback_data.radiation.Delta_t_cached_ti_begin)
    error(
        "Stale cached Delta_t: star %lld's step boundary moved since it "
        "was cached.",
        si->id);
#endif

  if (si->feedback_data.radiation.L_bol != 0.0) {
    const float p_rad = radiation_get_star_physical_radiation_pressure(
        si, Delta_t, phys_const, us, cosmo);
    const float delta_p_rad = weight * p_rad;

    /* Radially outwards from the star; * cosmo->a converts to comoving
       units. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.radiation.delta_p[i] -=
          delta_p_rad * dx[i] / r * cosmo->a;
    }

    /* Lifetime-cumulative diagnostic. delta_p_rad is the physical momentum
       magnitude for this pair, before it is projected onto the radial
       direction above; no separate energy channel (radiation pressure only
       deposits momentum). */
    xpj->feedback_data.radiation.cumulative_momentum += delta_p_rad;
    const float kick_velocity_rad = delta_p_rad / mj;
    if (kick_velocity_rad > xpj->feedback_data.radiation.max_kick_velocity) {
      xpj->feedback_data.radiation.max_kick_velocity = kick_velocity_rad;
    }

    /* Matches hit_by_SN/hit_by_winds: without it,
       feedback_update_part_radiation() never applies this momentum. */
    xpj->feedback_data.hit_by_radiation = 1;
  }

  /* Zero unless GEARFeedback:with_interstellar_radiation_field is on
     (L_band is then computed by stellar_evolution.c). */
  if (si->feedback_data.radiation.L_band[ISRF_MOMENT_PE] != 0.0 ||
      si->feedback_data.radiation.L_band[ISRF_MOMENT_LW] != 0.0) {

    const float Z_j = chemistry_get_total_metal_mass_fraction_for_cooling(pj);
    float extinction[ISRF_OPERATOR_COUNT];
    /* Receiver-side, using pj's own column density, not the source; see
       radiation_get_part_ISRF_extinction_factors for the formula. */
    const float extinction_path = radiation_get_comoving_extinction_path(
        fb_props, pj, xpj, r, cosmo, phys_const, hydro_props, us, cooling);
    radiation_get_part_ISRF_extinction_factors(us, cosmo, pj, Z_j, cooling,
                                               extinction_path, extinction);

    /* u_inject is an energy; dividing by mj below converts it to the
       specific energy each moment (or the dose reservoir) stores. Each
       moment takes its operator's extinction through the map. */
    double u_inject[ISRF_MOMENT_COUNT];
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
      u_inject[m] = (double)Delta_t * weight *
                    si->feedback_data.radiation.L_band[m] *
                    (double)extinction[radiation_isrf_moment_to_operator[m]];
    }

    if (fb_props->ISRF_propagation) {
      /* Dose-reservoir accumulator: pure accumulation of the elapsed star
         step's own (unrescaled) deposit, no reset, no first-touch logic, so
         any number of stars on any time bins just add without losing or
         double-counting emission. The rescale/phi fold-in happens once, at
         the receiving particle's own cadence, in
         radiation_end_force_propagation. */
      for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
        pj->feedback_data.isrf_moment[m].u_dose_reservoir +=
            (float)(u_inject[m] / (double)mj);
      }
      pj->feedback_data.ISRF_reservoir_end_ti =
          max(pj->feedback_data.ISRF_reservoir_end_ti, ti_current + ti_step);
      pj->feedback_data.ISRF_last_touch_ti = ti_current;
    } else {
      /* An instantaneous field strength, not an accumulated dose: reset to
         0 on the first touch this step (by any star), so a later read sees
         this step's illumination rather than a total across every step
         since the last cooling call. A later touch this same step (a
         second illuminating star) sums into what the first just wrote. */
      if (pj->feedback_data.ISRF_last_touch_ti != ti_current) {
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
          pj->feedback_data.isrf_moment[m].u = 0.f;
        pj->feedback_data.ISRF_last_touch_ti = ti_current;
      }

      for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
        pj->feedback_data.isrf_moment[m].u += (float)(u_inject[m] / (double)mj);
      }
    }

    /* Renew the illumination window on every touch, first or not: mirrors
       feedback_iact_HII_maintain_ionized_part's per-pass renewal of the HII
       tag's own end_time, so a continuously-illuminated particle's window
       never lapses between touches. The expiry check itself
       (radiation_reset_part_ISRF_illumination_tag) runs once per step in
       feedback_reset_part, not here. */
    pj->feedback_data.ISRF_illumination_end_ti =
        ti_current + RADIATION_ISRF_TAG_LIFETIME_INTERVALS * ti_step;

    /* First-touch-only sync, mirroring feedback_hii_claim_part vs.
       feedback_iact_HII_maintain_ionized_part's claim-vs-maintain split
       (feedback_common.c): do not re-sync an already-illuminated particle
       every pass, or every held particle drags the whole region down to
       the shortest time bin. */
    if (!pj->feedback_data.is_illuminated_ISRF) {
      pj->feedback_data.is_illuminated_ISRF = 1;
      timestep_sync_part(pj);
    }
  }
}

/**
 * @brief Update the properties of the particle due to radiation feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_update_part_radiation(struct part *p, struct xpart *xp,
                               const struct engine *e,
                               const float initial_mass) {

  /* Momentum only; internal energy is handled elsewhere, before cooling. */
  if (xp->feedback_data.hit_by_radiation) {
    for (int i = 0; i < 3; i++) {
      /* Initial mass of the gas, i.e. before any winds or SN. */
      const float dv = xp->feedback_data.radiation.delta_p[i] / initial_mass;
      xp->v_full[i] += dv;
      p->v[i] += dv;

      xp->feedback_data.radiation.delta_p[i] = 0;
    }
    xp->feedback_data.hit_by_radiation = 0;
  }
}

#endif /* SWIFT_RADIATION_IACT_GEAR_H */
