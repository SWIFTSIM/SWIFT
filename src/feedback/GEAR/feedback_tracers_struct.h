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
#ifndef SWIFT_FEEDBACK_TRACERS_STRUCT_GEAR_H
#define SWIFT_FEEDBACK_TRACERS_STRUCT_GEAR_H

/* Config parameters. */
#include <config.h>

/**
 * @brief Momentum and thermal energy of the feedback events of one step, kept
 * for the tracers of a gas particle (GEAR_thermal and GEAR_mechanical). The
 * tracers get them at the update, when the final mass and, for the mechanical
 * model, the multiple-event correction are known. The struct is empty without
 * the GEAR tracers, so it takes no memory.
 */
struct feedback_tracers_pending {
#if defined(TRACERS_GEAR)
  /*! Sum over the supernova events of the physical momentum received, in the
      frame of the particle (mass after the event times the velocity change) */
  float p_sum_SN;

  /*! Largest momentum of one supernova event, same frame */
  float p_max_SN;

  /*! Sum over the supernova events of the physical thermal energy given */
  float E_th_SN;

  /*! Same as p_sum_SN for the stellar winds */
  float p_sum_SW;

  /*! Same as p_max_SN for the stellar winds */
  float p_max_SW;

  /*! Same as E_th_SN for the stellar winds */
  float E_th_SW;
#endif /* defined(TRACERS_GEAR) */
};

#endif /* SWIFT_FEEDBACK_TRACERS_STRUCT_GEAR_H */
