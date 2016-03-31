/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include <float.h>
#include "potentials.h"
/**
 * @brief Computes the gravity time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */

__attribute__((always_inline)) INLINE static float gravity_compute_timestep(
    struct part* p, struct xpart* xp) {

  /* Currently no limit is imposed */
  return FLT_MAX;
}

__attribute__((always_inline)) INLINE static void gravity_init_gpart(struct gpart* g) {
  
  /* zero gravity acceleration */
  g->a_grav_external[0] = 0.f;
  g->a_grav_external[1] = 0.f;
  g->a_grav_external[2] = 0.f;

}



__attribute__((always_inline)) INLINE static float external_gravity_compute_timestep(struct gpart* g) {
  float dtmin = FLT_MAX;
#ifdef EXTERNAL_POTENTIAL_POINTMASS
  dtmin = fminf(dtmin, external_gravity_pointmass_timestep(g));
#endif
  return dtmin;
}

__attribute__((always_inline)) INLINE static void external_gravity(struct gpart *g)
{
#ifdef EXTERNAL_POTENTIAL_POINTMASS
  external_gravity_pointmass(g);
#endif

}
