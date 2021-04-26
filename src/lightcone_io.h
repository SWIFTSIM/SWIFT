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
#include <string.h>

/* This object's header. */
#include "lightcone_io.h"

/* Local headers. */
#include "common_io.h"
#include "error.h"
#include "part_type.h"

/* Forward declarations */
struct gpart;
struct part;
struct xpart;
struct spart;
struct bpart;
struct lightcone_props;

/*
 * Struct to describe an output field in the lightcone
 */
struct lightcone_io_props {
  
  /* Name */
  char name[FIELD_BUFFER_SIZE];

  /* Type of the field */
  enum IO_DATA_TYPE type;

  /* Dimension (1D, 3D, ...) */
  int dimension;
  
  /* Offset to this field in the data struct  */
  size_t offset;

};


inline static struct lightcone_io_props lightcone_io_make_output_field(
    char *name, enum IO_DATA_TYPE type, int dimension, ptrdiff_t offset) {
  
  struct lightcone_io_props r;
  bzero(&r, sizeof(struct lightcone_io_props));
  strcpy(r.name, name);
  r.type = type;
  r.dimension = dimension;
  r.offset = offset;
  return r;
}


/**
 * @brief Gas particle data for lightcone output
 */
struct lightcone_gas_data {
  long long id;
  double x[3];
  double mass;
};

void lightcone_store_gas(const struct gpart *gp, const struct part *p,
                         const struct xpart *xp, const double a_cross,
                         const double x_cross[3], struct lightcone_gas_data *data);


/**
 * @brief Dark matter particle data for lightcone output
 */
struct lightcone_dark_matter_data {
  long long id;
  double x[3];
  double mass;
};

void lightcone_store_dark_matter(const struct gpart *gp, const double a_cross,
                                 const double x_cross[3], struct lightcone_dark_matter_data *data);


/**
 * @brief Star particle data for lightcone output
 */
struct lightcone_stars_data {
  long long id;
  double x[3];
  double mass;
};

void lightcone_store_stars(const struct gpart *gp, const struct spart *sp,
                           const double a_cross, const double x_cross[3],
                           struct lightcone_stars_data *data);


/**
 * @brief Black hole particle data for lightcone output
 */
struct lightcone_black_hole_data {
  long long id;
  double x[3];
  double mass;
};

void lightcone_store_black_hole(const struct gpart *gp, const struct bpart *bp,
                                const double a_cross, const double x_cross[3],
                                struct lightcone_black_hole_data *data);


/**
 * @brief Neutrino particle data for lightcone output
 */
struct lightcone_neutrino_data {
  long long id;
  double x[3];
  double mass;
};

void lightcone_store_neutrino(const struct gpart *gp, const double a_cross,
                              const double x_cross[3],
                              struct lightcone_neutrino_data *data);


void lightcone_write_particles(struct lightcone_props *props, int ptype, hid_t file_id);


inline static size_t lightcone_io_struct_size(int ptype) {
  switch(ptype) {
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
    return sizeof(struct lightcone_dark_matter_data);
  case swift_type_gas:
    return sizeof(struct lightcone_gas_data);
  case swift_type_stars:
    return sizeof(struct lightcone_stars_data);
  case swift_type_black_hole:
    return sizeof(struct lightcone_black_hole_data);
  case swift_type_neutrino:
    return sizeof(struct lightcone_neutrino_data);
  default:
    error("Unhandled particle type");
    return 0;
  }
}

void lightcone_io_make_output_fields(void);

#endif
