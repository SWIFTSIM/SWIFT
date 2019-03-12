/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_VELOCIRAPTOR_PART_H
#define SWIFT_VELOCIRAPTOR_PART_H

#include "part_type.h"

/**
 * @brief SWIFT/VELOCIraptor particle.
 *
 * This should match the structure Swift::swift_vel_part
 * defined in the file NBodylib/src/NBody/SwiftParticle.h
 * of the VELOCIraptor code.
 */
struct swift_vel_part {

  /*! Particle ID. */
  long long id;

  /*! Particle position. */
  double x[3];

  /*! Particle velocity. */
  float v[3];

  /*! Particle mass. */
  float mass;

  /*! Gravitational potential */
  float potential;

  /*! Internal energy of gas particle */
  float u;

  /*! Temperature of a gas particle */
  float T;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;

  /*! MPI rank on which this #gpart lives on the SWIFT side. */
  int task;

  /*! Index of this #gpart in the global array of this rank on the SWIFT
    side. */
  int index;
};

#endif /* SWIFT_VELOCIRAPTOR_PART_H */
