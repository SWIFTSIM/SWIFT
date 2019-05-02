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
#ifndef SWIFT_EAGLE_BLACK_HOLE_PART_H
#define SWIFT_EAGLE_BLACK_HOLE_PART_H

/**
 * @brief Particle fields for the black hole particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct bpart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v[3];

  /*! Black hole mass */
  float mass;

  /* Particle cutoff radius. */
  float h;

  /*! Particle time bin */
  timebin_t time_bin;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Subgrid mass of the black hole */
  float subgrid_mass;

  /*! Union for the formation time and formation scale factor */
  union {

    /*! Formation time */
    float formation_time;

    /*! Formation scale factor */
    float formation_scale_factor;
  };

  /*! Energy reservoir for feedback */
  float energy_reservoir;

  /*! Instantaneous accretion rate */
  float accretion_rate;

  /*! Density of the gas surrounding the black hole. */
  float rho_gas;

  /*! Smoothed sound speed of the gas surrounding the black hole. */
  float sound_speed_gas;

  /*! Smoothed velocity (peculiar) of the gas surrounding the black hole */
  float velocity_gas[3];

  /*! Total mass of the gas neighbours. */
  float ngb_mass;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_EAGLE_BLACK_HOLE_PART_H */
