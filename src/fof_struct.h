/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FOF_STRUCT_H
#define SWIFT_FOF_STRUCT_H

/* Config parameters. */
#include <config.h>

#ifdef WITH_FOF

#ifdef WITH_FOF_GALAXIES

/**
 * @brief Particle-carried fields for the FoF galaxies scheme.
 */
struct group_data {

  /*! The gas+stellar mass of the host galaxy */
  float mass;

  /*! The stellar mass of the host galaxy */
  float stellar_mass;

  /*! The specific star formation rate of the host galaxy */
  float ssfr;

};

#endif

/**
 * @brief Particle-carried fields for the FoF scheme.
 */
struct fof_gpart_data {

  /*! Particle group ID */
  size_t group_id;

  /*! Size of the FOF group of this particle */
  size_t group_size;

#ifdef WITH_FOF_GALAXIES
  /*! The gas+stellar mass of the host galaxy */
  float group_mass;

  /*! The stellar mass of the host galaxy */
  float group_stellar_mass;

  /*! The star formation rate of the host galaxy */
  float group_sfr;

  /*! The star formation rate of the particle (duplicate of sf_data.SFR) */
  float part_sfr;

  /*! Is this particle able to form a group? */
  int is_grouppable;
#endif
};

#else

/**
 * @brief Particle-carried fields for the FoF scheme.
 */
struct fof_gpart_data {};


#endif

#endif /* SWIFT_FOF_STRUCT_H */
