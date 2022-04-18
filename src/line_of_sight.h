/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifndef SWIFT_LOS_H
#define SWIFT_LOS_H

/* Config parameters. */
#include "../config.h"
#include "engine.h"
#include "io_properties.h"

/**
 * @brief Maps the LOS axis geometry to the simulation axis geometry.
 *
 * Sigtlines will always shoot down the los_direction_z,
 * which can map to x,y or z of the simulation geometry.
 *
 * The remainng two axes, los_direction_x/y, then create
 * the plane orthogonal to the LOS direction. The random
 * sightline positions are created on this plane.
 */
enum los_direction { simulation_x_axis, simulation_y_axis, simulation_z_axis };

/**
 * @brief Properties of a single line-of-sight
 */
struct line_of_sight {

  /*! Simulation axis the LOS shoots down. */
  enum los_direction zaxis;

  /*! The two remaining axes defining the plane orthogonal to the sightline. */
  enum los_direction xaxis, yaxis;

  /*! Sightline position along los_direction_x. */
  double Xpos;

  /*! Sightline position along los_direction_y. */
  double Ypos;

  /*! Number of parts in LOS. */
  int particles_in_los_total;

  /*! Number of parts in LOS on this node. */
  int particles_in_los_local;

  /*! Is the simulation periodic? */
  int periodic;

  /*! Dimensions of the space. */
  double dim[3];

  /*! How many top level cells does ths LOS intersect? */
  int num_intersecting_top_level_cells;

  /*! The min--max range to consider for parts in LOS. */
  double range_when_shooting_down_axis[2];
};

/**
 * @brief Properties of the line-of-sight computation
 */
struct los_props {

  /*! Number of sightlines shooting down simulation z axis. */
  int num_along_z;

  /*! Number of sightlines shooting down simulation x axis. */
  int num_along_x;

  /*! Number of sightlines shooting down simulation y axis. */
  int num_along_y;

  /*! Total number of sightlines. */
  int num_tot;

  /*! The min--max range along the simulation x axis random sightlines are
   * allowed. */
  double allowed_losrange_x[2];

  /*! The min--max range along the simulation y axis random sightlines are
   * allowed. */
  double allowed_losrange_y[2];

  /*! The min--max range along the simulation z axis random sightlines are
   * allowed. */
  double allowed_losrange_z[2];

  /*! The min--max range to consider when LOS is shooting down each simulation
   * axis. */
  double range_when_shooting_down_axis[3][2];

  /*! Base name for line of sight HDF5 files. */
  char basename[200];
};

void print_los_info(const struct line_of_sight *Los, const int i);
void do_line_of_sight(struct engine *e);
void los_init(const double dim[3], struct los_props *los_params,
              struct swift_params *params);

void los_struct_dump(const struct los_props *internal_los, FILE *stream);
void los_struct_restore(const struct los_props *internal_los, FILE *stream);

#endif /* SWIFT_LOS_H */
