/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PART_LOCK_H
#define SWIFT_PART_LOCK_H

/* Config parameters. */
#include "config.h"

/* Local headers */
#include "error.h"
#include "lock.h"

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS

struct swift_particle_lock {};
#define swift_particle_lock_init(p) ;
#define swift_particle_lock_lock(p) ;
#define swift_particle_lock_unlock(p) ;
#else
#define swift_particle_lock_t swift_lock_type
#define swift_particle_lock_init(p) lock_init(&((p)->lock))
#define swift_particle_lock_lock(p) lock_lock(&((p)->lock))
#ifdef SWIFT_DEBUG_CHECKS
#define swift_particle_lock_unlock(p)                                        \
  ({                                                                         \
    if (lock_unlock(&((p)->lock)) != 0) error("Failed to unlock particle!"); \
  })
#else
#define swift_particle_lock_unlock(p)   \
  ({                                    \
    if (lock_unlock(&((p)->lock)) != 0) \
      ;                                 \
  })
#endif /* SWIFT_DEBUG_CHECKS */
#endif /* SWIFT_TASKS_WITHOUT_ATOMICS */

#endif /* SWIFT_PART_LOCK_H */
