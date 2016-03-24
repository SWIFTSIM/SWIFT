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



__attribute__((always_inline)) INLINE static float external_gravity_compute_timestep(
    struct gpart* g) {

  /* Currently no limit is imposed */
  const float dx   = g->x[0]-External_Potential_X;
  const float dy   = g->x[1]-External_Potential_Y;
  const float dz   = g->x[2]-External_Potential_Z;
  const float rinv = 1.f / sqrtf(dx*dx + dy*dy + dz*dz);
  const float drdv = (g->x[0]-External_Potential_X) * (g->v_full[0]) + (g->x[1]-External_Potential_Y) * (g->v_full[1]) + (g->x[2]-External_Potential_Z) * (g->v_full[2]); 
  const float dota_x = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[0] + 3.f * rinv * rinv * drdv * dx);
  const float dota_y = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[1] + 3.f * rinv * rinv * drdv * dy);
  const float dota_z = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[2] + 3.f * rinv * rinv * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2    = g->a_grav_external[0] * g->a_grav_external[0] + g->a_grav_external[1] * g->a_grav_external[1] + g->a_grav_external[2] * g->a_grav_external[2];
  
  return 0.03f * sqrtf(a_2/dota_2);
}
