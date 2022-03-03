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

#ifndef SWIFT_LIGHTCONE_MAP_TYPES_H
#define SWIFT_LIGHTCONE_MAP_TYPES_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "io_compression.h"
#include "parser.h"
#include "part_type.h"
#include "units.h"

/* Avoid cyclic inclusions */
struct cosmology;
struct engine;
struct lightcone_map;
struct lightcone_props;
struct gpart;

/* Type to store pointer to function for updating a healpix map */
typedef double (*map_update_function_t)(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/* Type to store pointer to function for providing baseline map value */
typedef double (*map_baseline_function_t)(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map);

/* Type to store pointer to function to check which types contribute to a map */
typedef int (*map_contrib_function_t)(int ptype);

enum lightcone_map_smoothing { map_unsmoothed, map_smoothed };

/**
 * @brief Struct to store information on one type of lightcone map
 */
struct lightcone_map_type {
  char name[PARSER_MAX_LINE_SIZE];
  map_update_function_t update_map;
  map_contrib_function_t ptype_contributes;
  map_baseline_function_t baseline_func;
  enum unit_conversion_factor units;
  enum lightcone_map_smoothing smoothing;
  enum lossy_compression_schemes compression;
  double buffer_scale_factor;
};

/*
  Function used for defining maps which only include gas (e.g. EAGLE x-ray
  outputs)
*/
int lightcone_map_gas_only(int ptype);

/*
   Healpix map of total mass
*/
int lightcone_map_total_mass_type_contributes(int ptype);

double lightcone_map_total_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

double lightcone_map_total_mass_baseline_value(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map);

/*
   Healpix map of gas mass
*/
int lightcone_map_gas_mass_type_contributes(int ptype);

double lightcone_map_gas_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of dark matter mass
*/
int lightcone_map_dark_matter_mass_type_contributes(int ptype);

double lightcone_map_dark_matter_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of stellar mass
*/
int lightcone_map_stellar_mass_type_contributes(int ptype);

double lightcone_map_stellar_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of black hole mass
*/
int lightcone_map_black_hole_mass_type_contributes(int ptype);

double lightcone_map_black_hole_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of star formation rate
*/
int lightcone_map_sfr_type_contributes(int ptype);

double lightcone_map_sfr_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/* This associates map names to the appropriate update function and unit info */
static const struct lightcone_map_type lightcone_map_types[] = {
    {
        .name = "TotalMass",
        .update_map = lightcone_map_total_mass_get_value,
        .ptype_contributes = lightcone_map_total_mass_type_contributes,
        .baseline_func = lightcone_map_total_mass_baseline_value,
        .units = UNIT_CONV_MASS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "SmoothedGasMass",
        .update_map = lightcone_map_gas_mass_get_value,
        .ptype_contributes = lightcone_map_gas_mass_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_MASS,
        .smoothing = map_smoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "UnsmoothedGasMass",
        .update_map = lightcone_map_gas_mass_get_value,
        .ptype_contributes = lightcone_map_gas_mass_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_MASS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "DarkMatterMass",
        .update_map = lightcone_map_dark_matter_mass_get_value,
        .ptype_contributes = lightcone_map_dark_matter_mass_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_MASS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "StellarMass",
        .update_map = lightcone_map_stellar_mass_get_value,
        .ptype_contributes = lightcone_map_stellar_mass_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_MASS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "BlackHoleMass",
        .update_map = lightcone_map_black_hole_mass_get_value,
        .ptype_contributes = lightcone_map_black_hole_mass_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_MASS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        .name = "StarFormationRate",
        .update_map = lightcone_map_sfr_get_value,
        .ptype_contributes = lightcone_map_sfr_type_contributes,
        .baseline_func = NULL,
        .units = UNIT_CONV_SFR,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
    {
        /* NULL functions indicate end of array */
        .name = "",
        .update_map = NULL,
        .ptype_contributes = NULL,
        .baseline_func = NULL,
        .units = UNIT_CONV_NO_UNITS,
        .smoothing = map_unsmoothed,
        .compression = compression_write_lossless,
        .buffer_scale_factor = 1.0,
    },
};

#endif
