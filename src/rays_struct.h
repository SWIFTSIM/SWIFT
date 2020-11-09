/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Evgenii Chaikin (chaikin@strw.leidenuniv.nl)
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
#ifndef SWIFT_RAYS_STRUCT_H
#define SWIFT_RAYS_STRUCT_H

/**
 * @brief Fields used in isotropic thermal SN/AGN feedback
 */
struct ray_data {

  /*! The mininum length (arc length or distance) between this
  ray and the gas neighbours in the stellar/BH kernel */
  float min_length;

  /*! The gas particle ID that has the minimal length
   * (arc length or distance) with respect to this ray !*/
  long long id_min_length;

  /*! Gas-particle mass in code units */
  float mass;
};

/**
 * @brief Additional fields used in isotropic kinetic SN feedback
 */
struct ray_data_extra {

  /*! Gas particle's comoving position in
   code units with respect to the star/BH */
  float x[3];

  /*! Gas particle's internal velocity in code units */
  float v[3];
};

#endif /* SWIFT_RAYS_STRUCT_H */
