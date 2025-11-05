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

#ifndef SWIFT_LIGHTCONE_PARTICLE_IO_H
#define SWIFT_LIGHTCONE_PARTICLE_IO_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <hdf5.h>
#include <string.h>

/* Local headers. */
#include "chemistry.h"
#include "common_io.h"
#include "error.h"
#include "io_compression.h"
#include "part_type.h"
#include "stars.h"
#include "units.h"

/* Forward declarations */
struct gpart;
struct part;
struct xpart;
struct spart;
struct bpart;
struct lightcone_props;
struct engine;

/*
 * Struct to describe an output field in the lightcone
 */
struct lightcone_io_field {

  /* Name */
  char name[FIELD_BUFFER_SIZE];

  /* Type of the field */
  enum IO_DATA_TYPE type;

  /* Dimension (1D, 3D, ...) */
  int dimension;

  /* Offset to this field in the data struct  */
  size_t offset;

  /* Units of this quantity */
  enum unit_conversion_factor units;

  /* Scale-factor exponent to apply for unit conversion to physical */
  float scale_factor_exponent;

  /* Lossy compression to use for this field */
  enum lossy_compression_schemes compression;

  /* Pointer to the next field */
  struct lightcone_io_field *next;
};

/*
 * Struct to store a linked list of lightcone_io_props
 */
struct lightcone_io_field_list {

  /* Pointer to the first field */
  struct lightcone_io_field *first;

  /* Pointer to the last field */
  struct lightcone_io_field *last;

  /* Number of fields */
  int num_fields;
};

/**
 * @brief Gas particle data for lightcone output
 */
struct lightcone_gas_data {
  long long id;
  double x[3];
  float vel[3];
  float mass;
  float a;
  float h;
  float rho;
  float temperature;
#ifdef CHEMISTRY_EAGLE
  float smoothed_metal_mass_fraction[chemistry_element_count];
  float metal_mass_fraction_total;
  float smoothed_metal_mass_fraction_total;
#endif
#ifdef COOLING_PS2020
  double electron_density;
  double ycompton;
#endif
#ifdef WITH_FOF
  long long group_id;
#endif
#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  float last_AGN_injection_scale_factor;
#endif
#ifdef STAR_FORMATION_EAGLE
  float sfr;
#endif
};

int lightcone_store_gas(const struct engine *e, struct lightcone_props *props,
                        const struct gpart *gp, const struct part *p,
                        const struct xpart *xp, const double a_cross,
                        const double x_cross[3],
                        struct lightcone_gas_data *data);

/**
 * @brief Dark matter particle data for lightcone output
 */
struct lightcone_dark_matter_data {
  long long id;
  double x[3];
  float vel[3];
  float mass;
  float a;
};

int lightcone_store_dark_matter(const struct engine *e,
                                struct lightcone_props *props,
                                const struct gpart *gp, const double a_cross,
                                const double x_cross[3],
                                struct lightcone_dark_matter_data *data);

/**
 * @brief Star particle data for lightcone output
 */
struct lightcone_stars_data {
  long long id;
  double x[3];
  float vel[3];
  float mass;
  float a;
#ifdef WITH_FOF
  long long group_id;
#endif
#ifdef STARS_EAGLE
  float mass_init;
  float birth_scale_factor;
  float birth_density;
  float luminosities[luminosity_bands_count];
#endif
#ifdef CHEMISTRY_EAGLE
  float smoothed_metal_mass_fraction[chemistry_element_count];
  float metal_mass_fraction_total;
  float smoothed_metal_mass_fraction_total;
#endif
#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  float last_AGN_injection_scale_factor;
#endif
};

int lightcone_store_stars(const struct engine *e, struct lightcone_props *props,
                          const struct gpart *gp, const struct spart *sp,
                          const double a_cross, const double x_cross[3],
                          struct lightcone_stars_data *data);

/**
 * @brief Black hole particle data for lightcone output
 */
struct lightcone_black_hole_data {
  long long id;
  double x[3];
  float vel[3];
  float mass;
  float a;
#ifdef BLACK_HOLES_EAGLE
  float subgrid_mass;
  float formation_scale_factor;
  float accretion_rate;
  float total_accreted_mass;
  float last_minor_merger_scale_factor;
  float last_major_merger_scale_factor;
  int number_of_mergers;
  float last_AGN_event_scale_factor;
  int AGN_number_of_AGN_events;
  int AGN_number_of_energy_injections;
  float last_high_Eddington_fraction_scale_factor;
  int cumulative_number_seeds;
#ifdef WITH_FOF
  long long group_id;
#endif
#endif
};

int lightcone_store_black_hole(const struct engine *e,
                               struct lightcone_props *props,
                               const struct gpart *gp, const struct bpart *bp,
                               const double a_cross, const double x_cross[3],
                               struct lightcone_black_hole_data *data);

/**
 * @brief Neutrino particle data for lightcone output
 */
struct lightcone_neutrino_data {
  long long id;
  double x[3];
  float vel[3];
  float mass;
  float weight;
  float a;
};

int lightcone_store_neutrino(const struct engine *e,
                             struct lightcone_props *props,
                             const struct gpart *gp, const double a_cross,
                             const double x_cross[3],
                             struct lightcone_neutrino_data *data);

void lightcone_write_particles(struct lightcone_props *props,
                               const struct unit_system *internal_units,
                               const struct unit_system *snapshot_units,
                               int ptype, hid_t file_id);

inline static size_t lightcone_io_struct_size(int ptype) {
  switch (ptype) {
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

void lightcone_io_field_list_init(struct lightcone_io_field_list *list);
void lightcone_io_field_list_clean(struct lightcone_io_field_list *list);
void lightcone_io_field_list_append(struct lightcone_io_field_list *list,
                                    char *name, enum IO_DATA_TYPE type,
                                    int dimension, size_t offset,
                                    enum unit_conversion_factor units,
                                    float scale_factor_exponent,
                                    char *compression);

void lightcone_io_append_gas_output_fields(
    struct lightcone_io_field_list *list);
void lightcone_io_append_dark_matter_output_fields(
    struct lightcone_io_field_list *list);
void lightcone_io_append_dark_matter_background_output_fields(
    struct lightcone_io_field_list *list);
void lightcone_io_append_stars_output_fields(
    struct lightcone_io_field_list *list);
void lightcone_io_append_black_hole_output_fields(
    struct lightcone_io_field_list *list);
void lightcone_io_append_neutrino_output_fields(
    struct lightcone_io_field_list *list);

#endif
