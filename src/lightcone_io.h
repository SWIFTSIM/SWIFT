/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_LIGHTCONE_IO_H
#define SWIFT_LIGHTCONE_IO_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <hdf5.h>

/* This object's header. */
#include "lightcone_io.h"

/* Forward declarations */
struct gpart;
struct lightcone_props;


/**
 * @brief Gas particle data for lightcone output
 */
struct lightcone_gas_data {
  long long id;
  double x[3];
};

/**
 * @brief Dark matter particle data for lightcone output
 */
struct lightcone_dark_matter_data {
  long long id;
  double x[3];
};

void lightcone_store_dark_matter(const struct gpart *gp, const double x_cross[3],
                                 struct lightcone_dark_matter_data *data);

void lightcone_write_dark_matter(struct lightcone_props *props, hid_t file_id, int ptype);

/**
 * @brief Dark matter background particle data for lightcone output
 */
struct lightcone_dark_matter_background_data {
  long long id;
  double x[3];
};

/**
 * @brief Star particle data for lightcone output
 */
struct lightcone_stars_data {
  long long id;
  double x[3];
};

/**
 * @brief Black hole particle data for lightcone output
 */
struct lightcone_black_hole_data {
  long long id;
  double x[3];
};

/**
 * @brief Neutrino particle data for lightcone output
 */
struct lightcone_neutrino_data {
  long long id;
  double x[3];
};

#endif
