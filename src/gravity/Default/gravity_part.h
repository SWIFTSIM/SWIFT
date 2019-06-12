/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_GRAVITY_PART_H
#define SWIFT_DEFAULT_GRAVITY_PART_H

/* Gravity particle. */
struct gpart {

  /*! Particle ID. If negative, it is the negative offset of the #part with
     which this gpart is linked. */
  long long id_or_neg_offset;

  /*! Particle position. */
  double x[3];

  /*! Particle velocity. */
  float v_full[3];

  /*! Particle acceleration. */
  float a_grav[3];

  /*! Particle mass. */
  float mass;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;

  /* Particle group ID and size in the FOF. */
  size_t group_id, group_size;

#ifdef SWIFT_DEBUG_CHECKS

  /* Numer of gparts this gpart interacted with */
  long long num_interacted;

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

  /* Has this particle been initialised? */
  int initialised;

#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  /*! Acceleration taken from the mesh only */
  float a_grav_PM[3];

  /*! Potential taken from the mesh only */
  float potential_PM;

  /* Brute-force particle acceleration. */
  double a_grav_exact[3];

  /* Brute-force particle potential. */
  double potential_exact;
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DEFAULT_GRAVITY_PART_H */
