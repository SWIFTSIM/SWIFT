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
#ifndef SWIFT_FEEDBACK_TRACERS_GEAR_THERMAL_H
#define SWIFT_FEEDBACK_TRACERS_GEAR_THERMAL_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "minmax.h"
#include "part.h"
#include "tracers.h"

/* The GEAR tracers get the momentum and the thermal energy that a gas particle
 * received over a step, once the final mass of the particle is known. The
 * events add their values to the pending fields of the #xpart, and
 * feedback_update_part() gives the sums to the tracers, once per channel and
 * per step. The specific energy and the kick speed are then divided by the
 * final mass, not by the mass after each event.
 * With the other tracers (none, EAGLE, FLAMINGO) nothing is stored and their
 * hooks are called by each event. */

#if !defined(TRACERS_GEAR)

/**
 * @brief Give the values of a supernova event to the tracers (not the GEAR
 * ones: there is nothing to store).
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SN(struct xpart *xp, const float dp_mag,
                                const double dE_th, const double new_mass) {
  const float new_mass_inv = new_mass > 0.0 ? (float)(1.0 / new_mass) : 0.0f;
  tracers_after_supernovae_feedback_part(xp, dp_mag, (float)dE_th * new_mass_inv,
                                         dp_mag * new_mass_inv);
}

/**
 * @brief Give the values of a stellar wind event to the tracers (not the GEAR
 * ones: there is nothing to store).
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SW(struct xpart *xp, const float dp_mag,
                                const double dE_th, const double new_mass) {
  const float new_mass_inv = new_mass > 0.0 ? (float)(1.0 / new_mass) : 0.0f;
  tracers_after_stellar_winds_feedback_part(
      xp, dp_mag, (float)dE_th * new_mass_inv, dp_mag * new_mass_inv);
}

/**
 * @brief Give the pending values to the tracers. Does nothing without the GEAR
 * tracers.
 *
 * @param xp The #xpart.
 * @param new_mass_inv The inverse of the mass of the #part after all events.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_update(struct xpart *xp, const float new_mass_inv) {}

/**
 * @brief Reset the pending values of the tracers. Does nothing without the
 * GEAR tracers.
 *
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_reset(struct xpart *xp) {}

#else

/**
 * @brief Keep the values of a supernova event for the tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas: the mass after
 * the event times the velocity change of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event (not used: the final
 * mass is used at the update).
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SN(struct xpart *xp, const float dp_mag,
                                const double dE_th, const double new_mass) {
  xp->feedback_data.tracers_pending.p_sum_SN += dp_mag;
  xp->feedback_data.tracers_pending.p_max_SN =
      max(xp->feedback_data.tracers_pending.p_max_SN, dp_mag);
  xp->feedback_data.tracers_pending.E_th_SN += dE_th;
}

/**
 * @brief Keep the values of a stellar wind event for the tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event (not used).
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SW(struct xpart *xp, const float dp_mag,
                                const double dE_th, const double new_mass) {
  xp->feedback_data.tracers_pending.p_sum_SW += dp_mag;
  xp->feedback_data.tracers_pending.p_max_SW =
      max(xp->feedback_data.tracers_pending.p_max_SW, dp_mag);
  xp->feedback_data.tracers_pending.E_th_SW += dE_th;
}

/**
 * @brief Give to the tracers the momentum and the thermal energy that a gas
 * particle received over the step, with the final mass.
 *
 * The momentum is the sum over the events of the momentum applied in the frame
 * of the gas. The specific thermal energy and the kick speed use the mass of
 * the particle after all the events.
 *
 * Note: This function is called in feedback_update_part(), before the reset of
 * the feedback fields.
 *
 * @param xp The #xpart.
 * @param new_mass_inv The inverse of the mass of the #part after all events.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_update(struct xpart *xp, const float new_mass_inv) {

  const struct feedback_tracers_pending *pending =
      &xp->feedback_data.tracers_pending;

  if (xp->feedback_data.hit_by_SN) {
    tracers_after_supernovae_feedback_part(
        xp, pending->p_sum_SN, (float)pending->E_th_SN * new_mass_inv,
        pending->p_max_SN * new_mass_inv);
  }

  if (xp->feedback_data.hit_by_winds) {
    tracers_after_stellar_winds_feedback_part(
        xp, pending->p_sum_SW, (float)pending->E_th_SW * new_mass_inv,
        pending->p_max_SW * new_mass_inv);
  }
}

/**
 * @brief Reset the pending values of the tracers.
 *
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_reset(struct xpart *xp) {
  xp->feedback_data.tracers_pending.p_sum_SN = 0.0f;
  xp->feedback_data.tracers_pending.p_max_SN = 0.0f;
  xp->feedback_data.tracers_pending.E_th_SN = 0.0f;
  xp->feedback_data.tracers_pending.p_sum_SW = 0.0f;
  xp->feedback_data.tracers_pending.p_max_SW = 0.0f;
  xp->feedback_data.tracers_pending.E_th_SW = 0.0f;
}

#endif /* !defined(TRACERS_GEAR) */

#endif /* SWIFT_FEEDBACK_TRACERS_GEAR_THERMAL_H */
