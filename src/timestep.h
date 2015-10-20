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

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */
__attribute__((always_inline)) INLINE static float compute_timestep_hydro(
    struct part* p, struct xpart* xp) {

  /* CFL condition */
  float dt_cfl = const_cfl * p->h / p->force.v_sig;

  /* Limit change in h */
  float dt_h_change = (p->force.h_dt != 0.0f)
                          ? fabsf(const_ln_max_h_change * p->h / p->force.h_dt)
                          : FLT_MAX;

  /* Limit change in u */
  float dt_u_change = (p->force.u_dt != 0.0f)
                          ? fabsf(const_max_u_change * p->u / p->force.u_dt)
                          : FLT_MAX;

  return fminf(dt_cfl, fminf(dt_h_change, dt_u_change));
}

/**
 * @brief Computes the gravity time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */

__attribute__((always_inline)) INLINE static float compute_timestep_grav(
    struct part* p, struct xpart* xp) {

  /* Currently no limit is imposed */
  return FLT_MAX;
}
