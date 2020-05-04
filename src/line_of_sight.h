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

/* Config parameters. */
#include "../config.h"
#include "engine.h"
#include "io_properties.h"

struct line_of_sight {
  /* Axis the sightline is shooting down.
   * 0 = simulation x axis.
   * 1 = simulation y axis.
   * 2 = simulation z axis. */
  int zaxis;

  /* The two remaining axes defining the plane orthogonal to the sightline. */
  int xaxis, yaxis;

  /* Random sightline position on the defined xaxis--yaxis plane. */
  double Xpos, Ypos;

  /* Number of parts in LOS. */
  size_t particles_in_los_total;

  /* Number of parts in LOS on this node. */
  size_t particles_in_los_local;

  /* Is the simulation periodic? */
  int periodic;

  /* Dimensions of the space. */
  double dim[3];
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
        const int xaxis, const int yaxis, const int zaxis,
        const int periodic, const double dim[3], struct line_of_sight *los);
void los_struct_dump(const struct los_props *internal_los,
                            FILE *stream);
void los_struct_restore(const struct los_props *internal_los,
                                       FILE *stream);
