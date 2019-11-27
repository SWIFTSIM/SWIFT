/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_MFM_HYDRO_GETTERS_H
#define SWIFT_GIZMO_MFM_HYDRO_GETTERS_H

/**
 * @brief Get a 5-element state vector W containing the primitive hydrodynamic
 * variables.
 *
 * @param p Particle.
 * @param W Pointer to the array in which the result needs to be stored (of size
 * 5 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_part_get_primitive_variables(const struct part *restrict p, float *W) {

  W[0] = p->rho;
  W[1] = p->v[0];
  W[2] = p->v[1];
  W[3] = p->v[2];
  W[4] = p->P;
}

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 */
#define hydro_part_geometry_well_behaved(p) \
  (p->geometry.wcorr > const_gizmo_min_wcorr)

#endif /* SWIFT_GIZMO_MFM_HYDRO_GETTERS_H */
