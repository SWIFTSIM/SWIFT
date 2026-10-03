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
#ifndef SWIFT_FEEDBACK_TRACERS_GEAR_MECHANICAL_H
#define SWIFT_FEEDBACK_TRACERS_GEAR_MECHANICAL_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "minmax.h"
#include "part.h"
#include "tracers.h"

/* The GEAR tracers get the momentum and the thermal energy that a gas particle
 * received, after the multiple-event correction. The events add their values
 * to the pending fields of the #xpart, and feedback_update_part() gives the
 * sum to the tracers, once per channel and per step (not once per event).
 * With the other tracers (none, EAGLE, FLAMINGO) nothing is stored and the
 * tracer hooks are not called by GEAR_mechanical. */

#if !defined(TRACERS_GEAR)

/**
 * @brief Keep the values of a supernova event for the tracers. Does nothing
 * without the GEAR tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dU The physical thermal energy given by the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SN(struct xpart *xp, const float dp_mag,
                                const double dU) {}

/**
 * @brief Keep the values of a stellar wind event for the tracers. Does nothing
 * without the GEAR tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dU The physical thermal energy given by the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SW(struct xpart *xp, const float dp_mag,
                                const double dU) {}

/**
 * @brief Give the pending values to the tracers. Does nothing without the GEAR
 * tracers.
 *
 * @param xp The #xpart.
 * @param f_corr The momentum correction factor.
 * @param u_residual The physical specific residual thermal energy.
 * @param new_mass_inv The inverse of the mass of the #part after all events.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_update(struct xpart *xp, const float f_corr,
                                const float u_residual,
                                const float new_mass_inv) {}

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
 * @param dU The physical thermal energy given by the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SN(struct xpart *xp, const float dp_mag,
                                const double dU) {
  xp->feedback_data.tracers_pending.p_sum_SN += dp_mag;
  xp->feedback_data.tracers_pending.p_max_SN =
      max(xp->feedback_data.tracers_pending.p_max_SN, dp_mag);
  xp->feedback_data.tracers_pending.E_th_SN += dU;
}

/**
 * @brief Keep the values of a stellar wind event for the tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dU The physical thermal energy given by the event.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_add_SW(struct xpart *xp, const float dp_mag,
                                const double dU) {
  xp->feedback_data.tracers_pending.p_sum_SW += dp_mag;
  xp->feedback_data.tracers_pending.p_max_SW =
      max(xp->feedback_data.tracers_pending.p_max_SW, dp_mag);
  xp->feedback_data.tracers_pending.E_th_SW += dU;
}

/**
 * @brief Give to the tracers the momentum and the thermal energy that a gas
 * particle received, once the multiple-event correction is known.
 *
 * The momentum is the one in the frame of the particle (the mass after the
 * event times the velocity change), for supernovae and winds alike. The
 * momentum and the kick speed are rescaled by the momentum correction factor.
 * The thermal energy of each channel is divided by the final mass, and the
 * residual thermal energy is shared between the channels in proportion to the
 * thermal energy that they gave.
 *
 * Note: This function is called in feedback_update_part(), before the reset of
 * the feedback fields.
 *
 * @param xp The #xpart.
 * @param f_corr The momentum correction factor.
 * @param u_residual The physical specific residual thermal energy.
 * @param new_mass_inv The inverse of the mass of the #part after all events.
 */
__attribute__((always_inline)) INLINE static void
feedback_tracers_pending_update(struct xpart *xp, const float f_corr,
                                const float u_residual,
                                const float new_mass_inv) {

  const struct feedback_tracers_pending *pending =
      &xp->feedback_data.tracers_pending;

  const float E_th_SN_pos = max(pending->E_th_SN, 0.0f);
  const float E_th_SW_pos = max(pending->E_th_SW, 0.0f);
  const float E_th_pos = E_th_SN_pos + E_th_SW_pos;

  /* Share of the residual thermal energy given to each channel */
  float share_SN = 0.0f;
  float share_SW = 0.0f;
  if (E_th_pos > 0.0f) {
    share_SN = E_th_SN_pos / E_th_pos;
    share_SW = 1.0f - share_SN;
  } else if (xp->feedback_data.number_SN > 0) {
    share_SN = 1.0f;
  } else {
    share_SW = 1.0f;
  }

  if (xp->feedback_data.number_SN > 0) {
    tracers_after_supernovae_feedback_part(
        xp, f_corr * pending->p_sum_SN,
        pending->E_th_SN * new_mass_inv + share_SN * u_residual,
        f_corr * pending->p_max_SN * new_mass_inv);
  }

  if (xp->feedback_data.number_winds > 0) {
    tracers_after_stellar_winds_feedback_part(
        xp, f_corr * pending->p_sum_SW,
        pending->E_th_SW * new_mass_inv + share_SW * u_residual,
        f_corr * pending->p_max_SW * new_mass_inv);
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

#endif /* SWIFT_FEEDBACK_TRACERS_GEAR_MECHANICAL_H */
