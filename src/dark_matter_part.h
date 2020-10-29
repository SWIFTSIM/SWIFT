/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_DARK_MATTER_PART_H
#define SWIFT_DARK_MATTER_PART_H

#include "dark_matter_properties.h"


/**
 * @brief Particle fields for the dark matter particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct dmpart {

  /*! Particle ID. */
  long long id_or_neg_offset;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v_full[3];

  /*! Dark matter mass */
  float mass;

  /*! Dark matter density */
  float rho;

  /* Particle cutoff radius. */
  float h;

  /*! Particle time bin */
  timebin_t time_bin;
    
#ifdef SWIFT_DEBUG_CHECKS
    
  /* Time of the last drift */
  integertime_t ti_drift;
    
  /* Time of the last kick */
  integertime_t ti_kick;
    
#endif
    
    /*! Average probability of scattering with another DM part */
    float sidm_probability;
    
    /*! Average relative velocity wrt to neighbours */
    float avg_pair_v;
    
    float time_step_size;
    
    float num_neighbours;


  struct {

      /*! Neighbour number count. */
      float wcount;
      
      /*! Derivative of the neighbour number with respect to h. */
      float wcount_dh;
      
      /*! Derivative of density with respect to h */
      float rho_dh;
      
  } density;

  /*! Add self-interacting DM specific stuff. */
  struct sidm_dmpart_data sidm_data;


} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DARK_MATTER_PART_H */
