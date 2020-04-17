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

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS

/* Place-holders for the case where tasks never access particlesc concurrently
 */

struct swift_particle_lock {};
#define swift_particle_lock_init(p) ;
#define swift_particle_lock_lock(p) ;
#define swift_particle_lock_unlock(p) ;
#define swift_particle_lock_lock_both(p, q) ;
#define swift_particle_lock_unlock_both(p, q) ;
#define swift_particle_lock_is_locked(p) 1
#else

/*! A particle-carried lock for neighbour interactions */
typedef volatile int swift_particle_lock_t;

/**
 * @brief Initialise the lock of a particle.
 *
 * @param p Pointer to the particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_init(p) ((p)->lock = 0)

/**
 * @brief Lock a particle for access within this thread.
 *
 * The thread will spin until the lock can be aquired.
 *
 * @param p Pointer to the particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_lock(p)               \
  ({                                              \
    while (atomic_cas(&((p)->lock), 0, 1) != 0) { \
      ;                                           \
    }                                             \
  })

/**
 * @brief Lock two particles for access within this thread.
 *
 * The thread will spin until the locks can be aquired.
 *
 * To avoid the dining philosopher's problem, we create a
 * hierarchy by using the address of the particles.
 *
 * @param p Pointer to the first particle (part, spart, gpart, bpart,...)
 * @param q Pointer to the second particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_lock_both(p, q) \
  ({                                        \
    if ((void*)p < (void*)q) {              \
      swift_particle_lock_lock(p);          \
      swift_particle_lock_lock(q);          \
    } else {                                \
      swift_particle_lock_lock(q);          \
      swift_particle_lock_lock(p);          \
    }                                       \
  })

#ifdef SWIFT_DEBUG_CHECKS

/**
 * @brief Release the lock of a particle.
 *
 * Raises and error if the lock can't be released.
 *
 * @param p Pointer to the particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_unlock(p)        \
  ({                                         \
    if (atomic_cas(&((p)->lock), 1, 0) != 1) \
      error("Failed to unlock particle!");   \
  })
#else

/**
 * @brief Release the lock of a particle.
 *
 * @param p Pointer to the particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_unlock(p) ({ atomic_cas(&((p)->lock), 1, 0); })
#endif /* SWIFT_DEBUG_CHECKS */

/**
 * @brief Release the lock of two particles.
 *
 * @param p Pointer to the first particle (part, spart, gpart, bpart,...)
 * @param q Pointer to the second particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_unlock_both(p, q) \
  ({                                          \
    swift_particle_lock_unlock(p);            \
    swift_particle_lock_unlock(q);            \
  })

/**
 * @brief Verifies whether a given particle lock has been aquired.
 *
 * Note that there is no way to know whether the lock has been
 * aquired by the calling thread.
 *
 * @param p Pointer to the particle (part, spart, gpart, bpart,...)
 */
#define swift_particle_lock_is_locked(p) ((p)->lock == 1)

#endif /* SWIFT_TASKS_WITHOUT_ATOMICS */

#endif /* SWIFT_PART_LOCK_H */
