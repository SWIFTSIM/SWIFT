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
#include "../GEAR/feedback_tracers_common.h"
#include "part.h"
#include "tracers.h"

/* GEAR_thermal gives the events to the tracers as GEAR_mechanical does, with
 * the code of ../GEAR/feedback_tracers_common.h. The only difference is that
 * the other tracers (none, EAGLE, FLAMINGO) keep being called by each event,
 * as before. */

/**
 * @brief Give the values of a supernova event to the tracers: they are kept for
 * the update with the GEAR tracers, and given at once to the other tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event.
 */
__attribute__((always_inline)) INLINE static void feedback_tracers_event_SN(
    struct xpart *xp, const float dp_mag, const double dE_th,
    const double new_mass) {
  feedback_tracers_pending_add_SN(xp, dp_mag, dE_th);
#if !defined(TRACERS_GEAR)
  const float new_mass_inv = new_mass > 0.0 ? (float)(1.0 / new_mass) : 0.0f;
  tracers_after_supernovae_feedback_part(
      xp, dp_mag, (float)dE_th * new_mass_inv, dp_mag * new_mass_inv);
#endif
}

/**
 * @brief Give the values of a stellar wind event to the tracers: they are kept
 * for the update with the GEAR tracers, and given at once to the other tracers.
 *
 * @param xp The #xpart.
 * @param dp_mag The physical momentum in the frame of the gas.
 * @param dE_th The physical thermal energy given by the event.
 * @param new_mass The mass of the #part after the event.
 */
__attribute__((always_inline)) INLINE static void feedback_tracers_event_SW(
    struct xpart *xp, const float dp_mag, const double dE_th,
    const double new_mass) {
  feedback_tracers_pending_add_SW(xp, dp_mag, dE_th);
#if !defined(TRACERS_GEAR)
  const float new_mass_inv = new_mass > 0.0 ? (float)(1.0 / new_mass) : 0.0f;
  tracers_after_stellar_winds_feedback_part(
      xp, dp_mag, (float)dE_th * new_mass_inv, dp_mag * new_mass_inv);
#endif
}

#endif /* SWIFT_FEEDBACK_TRACERS_GEAR_THERMAL_H */
