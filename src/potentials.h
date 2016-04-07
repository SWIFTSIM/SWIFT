/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#ifndef SWIFT_POTENTIALS_H
#define SWIFT_POTENTIALS_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* External Potential Properties */
struct external_potential {

#ifdef EXTERNAL_POTENTIAL_POINTMASS
  struct {
    double x, y, z;
    double mass;
  } point_mass;
#endif
};

/* Include exteral pointmass potential */
#ifdef EXTERNAL_POTENTIAL_POINTMASS

/**
 * @brief Computes the time-step due to the acceleration from a point mass
 *
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline))
    INLINE static float external_gravity_pointmass_timestep(
        const struct external_potential* potential,
        const struct phys_const* const phys_const,
        const struct gpart* const g) {

  const float G_newton = phys_const->const_newton_G;
  const float dx = g->x[0] - potential->point_mass.x;
  const float dy = g->x[1] - potential->point_mass.y;
  const float dz = g->x[2] - potential->point_mass.z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float drdv = (g->x[0] - potential->point_mass.x) * (g->v_full[0]) +
                     (g->x[1] - potential->point_mass.y) * (g->v_full[1]) +
                     (g->x[2] - potential->point_mass.z) * (g->v_full[2]);
  const float dota_x = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[0] + 3.f * rinv * rinv * drdv * dx);
  const float dota_y = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[1] + 3.f * rinv * rinv * drdv * dy);
  const float dota_z = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[2] + 3.f * rinv * rinv * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return 0.03f * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration of a particle due to a point
 * mass
 *
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_pointmass(
    const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float G_newton = phys_const->const_newton_G;
  const float dx = g->x[0] - potential->point_mass.x;
  const float dy = g->x[1] - potential->point_mass.y;
  const float dz = g->x[2] - potential->point_mass.z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);

  g->a_grav[0] +=
      -G_newton * potential->point_mass.mass * dx * rinv * rinv * rinv;
  g->a_grav[1] +=
      -G_newton * potential->point_mass.mass * dy * rinv * rinv * rinv;
  g->a_grav[2] +=
      -G_newton * potential->point_mass.mass * dz * rinv * rinv * rinv;
}
#endif /* EXTERNAL_POTENTIAL_POINTMASS */

/* Now, some generic functions, defined in the source file */

void potential_init(const struct swift_params* parameter_file,
		    struct UnitSystem* us,
		    struct external_potential* potential);

void potential_print(const struct external_potential* potential);

#endif /* SWIFT_POTENTIALS_H */
