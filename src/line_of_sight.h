/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* 
 * Maps the LOS axis geometry to the simulation axis geometry.
 *
 * Sigtlines will always shoot down the los_direction_z, 
 * which can map to x,y or z of the simulation geometry.
 *
 * The remainng two axes, los_direction_x/y, then create
 * the plane orthogonal to the LOS direction. The random
 * sightline positions are created on this plane.
 */
enum los_direction { 
   simulation_x_axis = 0,
   simulation_y_axis = 1,
   simulation_z_axis = 2
};

struct line_of_sight {
  /* Simulation axis the LOS shoots down. */
  enum los_direction zaxis;

  /* The two remaining axes defining the plane orthogonal to the sightline. */
  enum los_direction xaxis, yaxis;

  /* Sightline position along los_direction_x. */
  double Xpos;

  /* Sightline position along los_direction_y. */
  double Ypos;

  /* Number of parts in LOS. */
  size_t particles_in_los_total;

  /* Number of parts in LOS on this node. */
  size_t particles_in_los_local;

  /* Is the simulation periodic? */
  int periodic;

  /* Dimensions of the space. */
  double dim[3];

  /* Flag what top level cells this LOS intersects with */
  int *cells_top;

  /* How many top level cells does ths LOS intersect? */
  int num_intersecting_top_level_cells;
};

struct los_props {
  /* Number of sightlines shooting down simulation z axis. */
  int num_along_xy;

  /* Number of sightlines shooting down simulation x axis. */
  int num_along_yz;

  /* Number of sightlines shooting down simulation y axis. */
  int num_along_xz;

  /* Total number of sightlines. */
  int num_tot;

  /* The min--max range along the simulation x axis random sightlines are allowed. */
  double xmin, xmax;

  /* The min--max range along the simulation y axis random sightlines are allowed. */
  double ymin, ymax;

  /* The min--max range along the simulation z axis random sightlines are allowed. */
  double zmin, zmax;

  /* Basename for line of sight HDF5 files. */
  char basename[200];
};

double los_periodic(double x, double dim);
void generate_line_of_sights(struct line_of_sight *Los,
                             const struct los_props *params,
                             const int periodic, const double dim[3]);
void print_los_info(const struct line_of_sight *Los, const int i);
void do_line_of_sight(struct engine *e);
void los_init(double dim[3], struct los_props *los_params,
        struct swift_params *params);
void write_los_hdf5_datasets(hid_t grp, int j, size_t N, const struct part* parts,
                struct engine* e, const struct xpart* xparts);
void write_los_hdf5_dataset(const struct io_props p, size_t N, int j, struct engine* e, hid_t grp);
void write_hdf5_header(hid_t h_file, const struct engine *e, const struct los_props* LOS_params,
        const size_t total_num_parts_in_los);
void create_line_of_sight(const double Xpos, const double Ypos,
        enum los_direction xaxis, enum los_direction yaxis, enum los_direction zaxis,
        const int periodic, const double dim[3], struct line_of_sight *los);
void los_struct_dump(const struct los_props *internal_los,
                            FILE *stream);
void los_struct_restore(const struct los_props *internal_los,
                                       FILE *stream);
void find_intersecting_top_level_cells(const struct engine *e, struct line_of_sight *los);
int does_los_intersect(const struct cell *c, const struct line_of_sight *los);

#endif /* SWIFT_LOS_H */ 
