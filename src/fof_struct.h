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

/**
 * @brief Particle-carried fields for the FoF scheme.
 */
struct fof_gpart_data {

  /*! Particle group ID */
  size_t group_id;

  /*! Size of the FOF group of this particle */
  size_t group_size;

#ifdef WITH_HALO_FINDER
  
    /*! Particle host halo ID */
  size_t host_id;

  /*! Particle subhalo ID */
  size_t subhalo_id;

  /*! Size of the host of this particle */
  size_t host_size;
  
  /*! Size of the subhalo of this particle */
  size_t subhalo_size;

  /*! Mass of the group of this particle */
  double group_mass;

  /*! Mass of the host of this particle */
  double host_mass;

  /*! Particle previous group ID */
  size_t prev_group_id;

  /*! Particle previous host halo ID */
  size_t prev_host_id;

  /*! Particle previous subhalo ID */
  size_t prev_subhalo_id;
  
#endif
};

#else

/**
 * @brief Particle-carried fields for the FoF scheme.
 */
struct fof_gpart_data {};

#endif

#endif /* SWIFT_FOF_STRUCT_H */
