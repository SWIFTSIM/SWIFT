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
      float v_kick[3];

    } to_distribute;
  };
};

#endif /* SWIFT_FEEDBACK_STRUCT_SIMBA_H */
