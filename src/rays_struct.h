/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Evgenii Chaikin (chaikin@strw.leidenuniv.nl)
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
 * @brief The ray type: "thermal" or "kinetic". The kinetic type is then
 * subdivided into "kinetic true" and "kinetic mirr" types because in SNII
 * kinetic feedback, in one event we kick two particles in the opposite
 * directions; hence, we need two rays. The "true" ray is pointing to the 1st
 * particle and the second ray, which mirrors the direction of the first one, is
 * pointing to the 2nd particle.
 */
typedef enum ray_feedback_modes {
  ray_feedback_thermal,      /*< Used in AGN and SNII thermal feedback */
  ray_feedback_kinetic_true, /*< Used in SNII kinetic feedback (1st particle) */
  ray_feedback_kinetic_mirr  /*< Used in SNII kinetic feedback (2nd particle) */
} ray_feedback_type;

/**
 * @brief In SNII kinetic feedback, each ray needs to carry such an enum
 * to decide whether its kick will be allowed or not
 */
enum ray_feedback_kick_permission {
  ray_feedback_kick_non_allowed, /*< Stellar particle is not allowed to kick */
  ray_feedback_kick_allowed      /*< Stellar particle is allowed to kick */
};

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

  /* Does the gas particle this ray is pointing to allow this stellar particle
   * to kick it or not? Note, however, that the kick will only occur if both
   * gas particles in the pair agree to be kicked by this stellar particle. */
  enum ray_feedback_kick_permission status;
};

#endif /* SWIFT_RAYS_STRUCT_H */
