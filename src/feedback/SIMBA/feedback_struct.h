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
#ifndef SWIFT_FEEDBACK_STRUCT_SIMBA_H
#define SWIFT_FEEDBACK_STRUCT_SIMBA_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  union {

    /**
     * @brief Values collected from the gas neighbours.
     */
    struct {

    } to_collect;

    /**
     * @brief Values to be distributed to the gas neighbours.
     */
    struct {

      /* Velocity to update particles with */
      float v_kick;

      /* Remaining energy to distribute as heat */
      float delta_u; // ALEXEI: surely this should be an energy and not an internal energy because don't know particle's mass we're distributing to.

      /* Delay time */
      double simba_delay_time; // ALEXEI: think of a better place to put this

      /* wind mass loading */
      float wind_mass;

      /* Excess energy after kicking particle converted to thermal */
      float u_extra;

    } to_distribute;

    float host_galaxy_mass;
  };
};

#endif /* SWIFT_FEEDBACK_STRUCT_SIMBA_H */
