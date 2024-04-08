/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Marcel van Daalen (daalen@strw.leidenuniv.nl)
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
#include <config.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* Standard headers */
#include <stdio.h>
#include <string.h>

/* This object's header. */
#include "power_spectrum.h"

/* Local includes. */
#include "cooling.h"
#include "engine.h"
#include "minmax.h"
#include "neutrino.h"
#include "random.h"
#include "row_major_id.h"
#include "space.h"
#include "threadpool.h"
#include "tools.h"
#include "units.h"

#define power_data_default_grid_side_length 256
#define power_data_default_fold_factor 4
#define power_data_default_window_order 3

/* The way to calculate these shifts is to consider a 3D cube of (kx,ky,kz)
 * cells and check which cells fall inside a spherical shell with boundaries
 * (i+0.5,i+1.5), then calculate the average k=sqrt(kx^2+ky^2+kz^2). So for i=0
 * you'd find 6 cells k=1 and 12 cells k=sqrt(2), so the weighted k becomes
 * (6 * 1 + 12 * sqrt(2)) / 18 = 1.2761424 – etc.
 * Note that beyond the 7th term, the correction is < 1%. */
#define number_of_corrected_bins 128
static const float correction_shift_k_values[number_of_corrected_bins] = {
    1.2761424f, 1.1154015f, 1.0447197f, 1.0151449f, 1.0195166f, 1.0203214f,
    1.0102490f, 1.0031348f, 1.0063766f, 1.0093355f, 1.0055681f, 1.0024279f,
    1.0034435f, 1.0038386f, 1.0011069f, 1.0002888f, 1.0018693f, 1.0029172f,
    1.0019128f, 1.0009282f, 1.0015312f, 1.0016361f, 1.0009436f, 1.0003777f,
    1.0005931f, 1.0010948f, 1.0010581f, 1.0009779f, 1.0010282f, 1.0008224f,
    1.0006637f, 1.0004002f, 1.0002419f, 1.0005172f, 1.0005523f, 1.0004342f,
    1.0005183f, 1.0005357f, 1.0003162f, 1.0001836f, 1.0003737f, 1.0004792f,
    1.0004169f, 1.0003660f, 1.0004468f, 1.0004218f, 1.0001436f, 1.0000479f,
    1.0002012f, 1.0003710f, 1.0003234f, 1.0002661f, 1.0003446f, 1.0003313f,
    1.0001844f, 1.0000630f, 1.0001714f, 1.0002382f, 1.0001507f, 1.0001663f,
    1.0002199f, 1.0002403f, 1.0000911f, 0.9999714f, 1.0001136f, 1.0001907f,
    1.0001917f, 1.0001684f, 1.0001875f, 1.0002158f, 1.0000941f, 1.0000646f,
    1.0000930f, 1.0001497f, 1.0001589f, 1.0001215f, 1.0001563f, 1.0001254f,
    1.0000557f, 1.0000220f, 1.0000517f, 1.0001039f, 1.0001185f, 1.0000778f,
    1.0000848f, 1.0001415f, 1.0001108f, 1.0000709f, 1.0000724f, 1.0001201f,
    1.0001480f, 1.0001204f, 1.0001185f, 1.0000844f, 1.0000224f, 0.9999752f,
    0.9999997f, 1.0000969f, 1.0001076f, 1.0000756f, 1.0000700f, 1.0000854f,
    1.0001067f, 1.0000390f, 1.0000443f, 1.0000863f, 1.0000585f, 1.0000352f,
    1.0000677f, 1.0001081f, 1.0000537f, 1.0000199f, 1.0000308f, 1.0000585f,
    1.0000479f, 1.0000304f, 1.0000751f, 1.0000710f, 1.0000152f, 1.0000083f,
    1.0000342f, 1.0000530f, 1.0000543f, 1.0000442f, 1.0000680f, 1.0000753f,
    1.0000369f, 1.0000117f};

#ifdef HAVE_FFTW

/**
 * @brief Return the #power_type corresponding to a given string.
 */
INLINE static enum power_type power_spectrum_get_type(const char* name) {
  if (strcasecmp(name, "matter") == 0)
    return pow_type_matter;
  else if (strcasecmp(name, "cdm") == 0)
    return pow_type_cdm;
  else if (strcasecmp(name, "gas") == 0)
    return pow_type_gas;
  else if (strcasecmp(name, "starBH") == 0)
    return pow_type_starBH;
  else if (strcasecmp(name, "pressure") == 0)
    return pow_type_pressure;
  else if (strcasecmp(name, "neutrino") == 0)
    return pow_type_neutrino;
  else if (strcasecmp(name, "neutrino0") == 0)
    return pow_type_neutrino_0;
  else if (strcasecmp(name, "neutrino1") == 0)
    return pow_type_neutrino_1;
  else
    error("Do not recognize the power spectrum type '%s'.", name);
  return pow_type_count;
}

/**
 * @brief Get the name of the type of power component.
 *
 * @param type The #power_type for which we want the name.
 */
INLINE static const char* get_powtype_name(const enum power_type type) {

  static const char* powtype_names[pow_type_count] = {"Matter",
                                                      "CDM",
                                                      "gas",
                                                      "stars/BHs",
                                                      "electron pressure",
                                                      "neutrino",
                                                      "neutrino (even)",
                                                      "neutrino (odd)"};

  return powtype_names[type];
}

INLINE static const char* get_powtype_filename(const enum power_type type) {

  static const char* powtype_filenames[pow_type_count] = {
      "matter",   "cdm",      "gas",       "starBH",
      "pressure", "neutrino", "neutrino0", "neutrino1"};

  return powtype_filenames[type];
}

/**
 * @brief Shared information for shot noise to be used by all the threads in the
 * pool.
 */
struct shot_mapper_data {
  const struct cell* cells;
  double tot12;
  enum power_type type1;
  enum power_type type2;
  const struct engine* e;
  struct neutrino_model* nu_model;
};

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct grid_mapper_data {
  const struct cell* cells;
  double* dens;
  int N;
  enum power_type type;
  int windoworder;
  double dim[3];
  double fac;
  const struct engine* e;
  struct neutrino_model* nu_model;
};

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct conv_mapper_data {
  double* grid;
  int Ngrid;
  double invcellmean;
};

/**
 * @brief Shared information needed for calculating power from a Fourier grid.
 */
struct pow_mapper_data {
  fftw_complex* powgridft;
  fftw_complex* powgridft2;
  int Ngrid;
  int windoworder;
  int* kbin;
  int* modecounts;
  double* powersum;
  double jfac;
};

/**
 * @brief Decided whether or not a given particle should be considered or
 * not for the specific #power_type we are computing.
 */
int should_collect_mass(const enum power_type type, const struct gpart* gp,
                        const integertime_t ti_current) {

  if (type == pow_type_matter) {
    /* Select all particles */
    return 1;
  } else if (type == pow_type_cdm) {
    if (gp->type == swift_type_dark_matter ||
        gp->type == swift_type_dark_matter_background)
      return 1;
  } else if (type == pow_type_gas) {
    if (gp->type == swift_type_gas) return 1;
  } else if (type == pow_type_starBH) {
    if (gp->type == swift_type_stars || gp->type == swift_type_black_hole)
      return 1;
  } else if (type == pow_type_neutrino) {
    if (gp->type == swift_type_neutrino) return 1;
  } else if (type == pow_type_neutrino_0) {
    if (gp->type == swift_type_neutrino) {
      /* Randomly divide the neutrino ensemble in half */
      return random_unit_interval(gp->id_or_neg_offset, ti_current,
                                  random_number_powerspectrum_split) < 0.5;
    }
  } else if (type == pow_type_neutrino_1) {
    if (gp->type == swift_type_neutrino)
      /* Randomly divide the neutrino ensemble in half */
      return random_unit_interval(gp->id_or_neg_offset, ti_current,
                                  random_number_powerspectrum_split) >= 0.5;
  } else {
#ifdef SWIFT_DEBUG_CHECKS
    error("Invalid type!");
#endif
  }

  return 0;
}

/**
 * @brief Calculates the necessary mass terms for shot noise.
 *
 * @param c The #cell.
 * @param tot12 The shot noise contribution returned.
 * @param type1 The component type of field 1.
 * @param type2 The component type of field 2.
 * @param e The #engine.
 */
void shotnoiseterms(const struct cell* c, double* tot12,
                    const enum power_type type1, const enum power_type type2,
                    const struct engine* e, struct neutrino_model* nu_model) {

  const int gcount = c->grav.count;
  const struct gpart* gparts = c->grav.parts;

  /* Handle on the other particle types */
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  /* Local accumulator for this cell */
  double local_tot12 = 0.;

  /* Calculate the value each particle adds to the grids */
  for (int i = 0; i < gcount; ++i) {

    /* Skip invalid particles */
    if (gparts[i].time_bin == time_bin_inhibited) continue;

    double quantity1;

    /* Special case first for the electron pressure */
    if (type1 == pow_type_pressure) {

      /* Skip non-gas particles */
      if (gparts[i].type != swift_type_gas) continue;

      const struct part* p = &parts[-gparts[i].id_or_neg_offset];
      const struct xpart* xp = &xparts[-gparts[i].id_or_neg_offset];
      quantity1 = cooling_get_electron_pressure(phys_const, hydro_props, us,
                                                cosmo, cool_func, p, xp);
    } else {

      /* We are collecting a mass of some kind.
       * We skip any particle not matching the PS type we want */
      if (!should_collect_mass(type1, &gparts[i], e->ti_current)) continue;

      /* Compute weight (for neutrino delta-f weighting) */
      double weight1 = 1.0;
      if (gparts[i].type == swift_type_neutrino)
        gpart_neutrino_weight_mesh_only(&gparts[i], nu_model, &weight1);

      /* And eventually... collect */
      quantity1 = gparts[i].mass * weight1;
    }

    /* Can we assign already? (i.e. this is an auto-spectrum) */
    if (type1 == type2) {

      /* Now assign to the shot noise collection */
      local_tot12 += quantity1 * quantity1;
    }

    /* This is a cross-spectrum */
    else {

      double quantity2;

      /* Special case first for the electron pressure */
      if (type2 == pow_type_pressure) {

        /* Skip non-gas particles */
        if (gparts[i].type != swift_type_gas) continue;

        const struct part* p = &parts[-gparts[i].id_or_neg_offset];
        const struct xpart* xp = &xparts[-gparts[i].id_or_neg_offset];
        quantity2 = cooling_get_electron_pressure(phys_const, hydro_props, us,
                                                  cosmo, cool_func, p, xp);
      } else {

        /* We are collecting a mass of some kind.
         * We skip any particle not matching the PS type we want */
        if (!should_collect_mass(type2, &gparts[i], e->ti_current)) continue;

        /* Compute weight (for neutrino delta-f weighting) */
        double weight2 = 1.0;
        if (gparts[i].type == swift_type_neutrino)
          gpart_neutrino_weight_mesh_only(&gparts[i], nu_model, &weight2);

        /* And eventually... collect */
        quantity2 = gparts[i].mass * weight2;
      }

      /* Now assign to the shot noise collection */
      local_tot12 += quantity1 * quantity2;
    }

  } /* Loop over particles */

  /* Now that we are done with this cell, write back to the global accumulator
   */
  atomic_add_d(tot12, local_tot12);
}

/**
 * @brief Threadpool mapper function for the shot noise calculation.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the cells.
 */
void shotnoise_mapper(void* map_data, int num, void* extra) {

  /* Unpack the shared information */
  struct shot_mapper_data* data = (struct shot_mapper_data*)extra;
  const struct cell* cells = data->cells;
  const struct engine* e = data->e;
  struct neutrino_model* nu_model = data->nu_model;

  /* Pointer to the chunk to be processed */
  int* local_cells = (int*)map_data;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {
    /* Pointer to local cell */
    const struct cell* c = &cells[local_cells[i]];

    /* Calculate the necessary mass terms */
    shotnoiseterms(c, &data->tot12, data->type1, data->type2, e, nu_model);
  }
}

__attribute__((always_inline)) INLINE static void TSC_set(
    double* mesh, const int N, const int i, const int j, const int k,
    const double dx, const double dy, const double dz, const double value) {

  const double lx = 0.5 * (0.5 - dx) * (0.5 - dx); /* left side, dist 1 + dx  */
  const double mx = 0.75 - dx * dx;                /* center, dist |dx|  */
  const double rx = 0.5 * (0.5 + dx) * (0.5 + dx); /* right side, dist 1 - dx */

  const double ly = 0.5 * (0.5 - dy) * (0.5 - dy); /* left side, dist 1 + dy  */
  const double my = 0.75 - dy * dy;                /* center, dist |dy|  */
  const double ry = 0.5 * (0.5 + dy) * (0.5 + dy); /* right side, dist 1 - dy */

  const double lz = 0.5 * (0.5 - dz) * (0.5 - dz); /* left side, dist 1 + dz  */
  const double mz = 0.75 - dz * dz;                /* center, dist |dz|  */
  const double rz = 0.5 * (0.5 + dz) * (0.5 + dz); /* right side, dist 1 - dz */

  /* TSC interpolation */
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j - 1, k - 1, N, 2)],
      value * lx * ly * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j - 1, k + 0, N, 2)],
      value * lx * ly * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j - 1, k + 1, N, 2)],
      value * lx * ly * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 0, k - 1, N, 2)],
      value * lx * my * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 0, k + 0, N, 2)],
      value * lx * my * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 0, k + 1, N, 2)],
      value * lx * my * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 1, k - 1, N, 2)],
      value * lx * ry * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 1, k + 0, N, 2)],
      value * lx * ry * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i - 1, j + 1, k + 1, N, 2)],
      value * lx * ry * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j - 1, k - 1, N, 2)],
      value * mx * ly * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j - 1, k + 0, N, 2)],
      value * mx * ly * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j - 1, k + 1, N, 2)],
      value * mx * ly * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 0, k - 1, N, 2)],
      value * mx * my * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 0, k + 0, N, 2)],
      value * mx * my * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 0, k + 1, N, 2)],
      value * mx * my * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 1, k - 1, N, 2)],
      value * mx * ry * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 1, k + 0, N, 2)],
      value * mx * ry * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 1, k + 1, N, 2)],
      value * mx * ry * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j - 1, k - 1, N, 2)],
      value * rx * ly * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j - 1, k + 0, N, 2)],
      value * rx * ly * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j - 1, k + 1, N, 2)],
      value * rx * ly * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 0, k - 1, N, 2)],
      value * rx * my * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 0, k + 0, N, 2)],
      value * rx * my * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 0, k + 1, N, 2)],
      value * rx * my * rz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 1, k - 1, N, 2)],
      value * rx * ry * lz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 1, k + 0, N, 2)],
      value * rx * ry * mz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 1, k + 1, N, 2)],
      value * rx * ry * rz);
}

INLINE static void gpart_to_grid_TSC(const struct gpart* gp, double* rho,
                                     const int N, const double fac,
                                     const double dim[3], const double value) {

  /* Fold the particle position position */
  const double pos_x = box_wrap_multiple(gp->x[0], 0., dim[0]) * fac;
  const double pos_y = box_wrap_multiple(gp->x[1], 0., dim[1]) * fac;
  const double pos_z = box_wrap_multiple(gp->x[2], 0., dim[2]) * fac;

  /* Workout the TSC coefficients */
  const int i = (int)(pos_x + 0.5);
  const double dx = pos_x - i;

  const int j = (int)(pos_y + 0.5);
  const double dy = pos_y - j;

  const int k = (int)(pos_z + 0.5);
  const double dz = pos_z - k;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i > N) error("Invalid gpart position in x");
  if (j < 0 || j > N) error("Invalid gpart position in y");
  if (k < 0 || k > N) error("Invalid gpart position in z");
#endif

  TSC_set(rho, N, i, j, k, dx, dy, dz, value);
}

__attribute__((always_inline)) INLINE static void CIC_set(
    double* mesh, const int N, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz, const double value) {

  /* Classic CIC interpolation */
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 0, k + 0, N, 2)],
      value * tx * ty * tz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 0, k + 1, N, 2)],
      value * tx * ty * dz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 1, k + 0, N, 2)],
      value * tx * dy * tz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 0, j + 1, k + 1, N, 2)],
      value * tx * dy * dz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 0, k + 0, N, 2)],
      value * dx * ty * tz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 0, k + 1, N, 2)],
      value * dx * ty * dz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 1, k + 0, N, 2)],
      value * dx * dy * tz);
  atomic_add_d(
      &mesh[row_major_id_periodic_with_padding(i + 1, j + 1, k + 1, N, 2)],
      value * dx * dy * dz);
}

INLINE static void gpart_to_grid_CIC(const struct gpart* gp, double* rho,
                                     const int N, const double fac,
                                     const double dim[3], const double value) {

  /* Fold the particle position position */
  const double pos_x = box_wrap_multiple(gp->x[0], 0., dim[0]) * fac;
  const double pos_y = box_wrap_multiple(gp->x[1], 0., dim[1]) * fac;
  const double pos_z = box_wrap_multiple(gp->x[2], 0., dim[2]) * fac;

  /* Workout the CIC coefficients */
  const int i = (int)pos_x;
  const double dx = pos_x - i;
  const double tx = 1. - dx;

  const int j = (int)pos_y;
  const double dy = pos_y - j;
  const double ty = 1. - dy;

  const int k = (int)pos_z;
  const double dz = pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i > N) error("Invalid gpart position in x");
  if (j < 0 || j > N) error("Invalid gpart position in y");
  if (k < 0 || k > N) error("Invalid gpart position in z");
#endif

  CIC_set(rho, N, i, j, k, tx, ty, tz, dx, dy, dz, value);
}

INLINE static void gpart_to_grid_NGP(const struct gpart* gp, double* rho,
                                     const int N, const double fac,
                                     const double dim[3], const double value) {

  /* Fold the particle position position */
  const double pos_x = box_wrap_multiple(gp->x[0], 0., dim[0]) * fac;
  const double pos_y = box_wrap_multiple(gp->x[1], 0., dim[1]) * fac;
  const double pos_z = box_wrap_multiple(gp->x[2], 0., dim[2]) * fac;

  const int xi = (int)(pos_x + 0.5) % N;
  const int yi = (int)(pos_y + 0.5) % N;
  const int zi = (int)(pos_z + 0.5) % N;
  atomic_add_d(&rho[(xi * N + yi) * (N + 2) + zi], value);
}

/**
 * @brief Assigns all the #gpart of a #cell to the power grid using the
 * chosen mass assignment method.
 *
 * @param c The #cell.
 * @param rho The density grid.
 * @param N the size of the grid along one axis.
 * @param fac Conversion factor of wrapped position to grid.
 * @param type The #power_type we want to assign to the grid.
 * @param windoworder The window to use for grid assignment.
 * @param e The #engine.
 */
void cell_to_powgrid(const struct cell* c, double* rho, const int N,
                     const double fac, const enum power_type type,
                     const int windoworder, const double dim[3],
                     const struct engine* e, struct neutrino_model* nu_model) {

  const int gcount = c->grav.count;
  const struct gpart* gparts = c->grav.parts;

  /* Handle on the other particle types */
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  /* Assign all the gpart of that cell to the mesh */
  for (int i = 0; i < gcount; ++i) {

    /* Skip invalid particles */
    if (gparts[i].time_bin == time_bin_inhibited) continue;

    /* Collect the quantity to assign to the mesh */
    double quantity;

    /* Special case first for the electron pressure */
    if (type == pow_type_pressure) {

      /* Skip non-gas particles */
      if (gparts[i].type != swift_type_gas) continue;

      const struct part* p = &parts[-gparts[i].id_or_neg_offset];
      const struct xpart* xp = &xparts[-gparts[i].id_or_neg_offset];
      quantity = cooling_get_electron_pressure(phys_const, hydro_props, us,
                                               cosmo, cool_func, p, xp);
    } else {

      /* We are collecting a mass of some kind.
       * We skip any particle not matching the PS type we want */
      if (!should_collect_mass(type, &gparts[i], e->ti_current)) continue;

      /* Compute weight (for neutrino delta-f weighting) */
      double weight = 1.0;
      if (gparts[i].type == swift_type_neutrino)
        gpart_neutrino_weight_mesh_only(&gparts[i], nu_model, &weight);

      /* And eventually... collect */
      quantity = gparts[i].mass * weight;
    }

    /* Assign the quantity to the grid */
    switch (windoworder) {
      case 1:
        gpart_to_grid_NGP(&gparts[i], rho, N, fac, dim, quantity);
        break;
      case 2:
        gpart_to_grid_CIC(&gparts[i], rho, N, fac, dim, quantity);
        break;
      case 3:
        gpart_to_grid_TSC(&gparts[i], rho, N, fac, dim, quantity);
        break;
      default:
#ifdef SWIFT_DEBUG_CHECKS
        error("Not implemented!");
#endif
        break;
    }

  } /* Loop over particles */
}

/**
 * @brief Threadpool mapper function for the power grid assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the grid and cells.
 */
void cell_to_powgrid_mapper(void* map_data, int num, void* extra) {

  /* Unpack the shared information */
  const struct grid_mapper_data* data = (struct grid_mapper_data*)extra;
  const struct cell* cells = data->cells;
  double* grid = data->dens;
  const int Ngrid = data->N;
  const enum power_type type = data->type;
  const int order = data->windoworder;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};
  const double gridfac = data->fac;
  const struct engine* e = data->e;
  struct neutrino_model* nu_model = data->nu_model;

  /* Pointer to the chunk to be processed */
  int* local_cells = (int*)map_data;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {
    /* Pointer to local cell */
    const struct cell* c = &cells[local_cells[i]];

    /* Assign this cell's content to the grid */
    cell_to_powgrid(c, grid, Ngrid, gridfac, type, order, dim, e, nu_model);
  }
}

/**
 * @brief Mapper function for the conversion of mass to density contrast.
 *
 * @param map_data The density field array.
 * @param num The number of elements to iterate on (along the x-axis).
 * @param extra The information about the grid and conversion.
 */
void mass_to_contrast_mapper(void* map_data, int num, void* extra) {

  /* Unpack the shared information */
  const struct conv_mapper_data* data = (struct conv_mapper_data*)extra;
  double* grid = data->grid;
  const int Ngrid = data->Ngrid;
  const double invcellmean = data->invcellmean;

  /* Range handled by this call */
  const int xi_start = (double*)map_data - grid;
  const int xi_end = xi_start + num;

  /* Loop over the assigned cells, convert to density contrast */
  for (int xi = xi_start; xi < xi_end; ++xi) {
    for (int yi = 0; yi < Ngrid; ++yi) {
      for (int zi = 0; zi < Ngrid; ++zi) {
        const int index = (xi * Ngrid + yi) * (Ngrid + 2) +
                          zi; /* Ngrid+2 because of padding */
        grid[index] = grid[index] * invcellmean; /* -1.0; -1 only affects k=0 */
      }
    }
  }
}

/**
 * @brief Mapper function for calculating the power from a Fourier grid.
 *
 * @param map_data The array of the density field Fourier transform.
 * @param num The number of elements to iterate on (along the x-axis).
 * @param extra Arrays to store the results/helper variables.
 */
void pow_from_grid_mapper(void* map_data, const int num, void* extra) {

  struct pow_mapper_data* data = (struct pow_mapper_data*)extra;

  /* Unpack the data struct */
  fftw_complex* restrict powgridft = data->powgridft;
  fftw_complex* restrict powgridft2 = data->powgridft2;
  const int Ngrid = data->Ngrid;
  const int Nhalf = Ngrid / 2;
  const int nyq2 = Nhalf * Nhalf;
  const int windoworder = data->windoworder;
  const int* restrict kbin = data->kbin;
  const double jfac = data->jfac;

  /* Output data */
  int* restrict modecounts = data->modecounts;
  double* restrict powersum = data->powersum;

  /* Range handled by this call */
  const int xi_start = (fftw_complex*)map_data - powgridft;
  const int xi_end = xi_start + num;

  /* Loop over the assigned FT'd cells, get deconvolved power from them */
  for (int xi = xi_start; xi < xi_end; ++xi) {

    int kx = xi;
    if (kx > Nhalf) kx -= Ngrid;
    const double fx = kx * jfac;
    const double invsincx = (xi == 0) ? 1.0 : fx / sin(fx);

    for (int yi = 0; yi < Ngrid; ++yi) {

      int ky = yi;
      if (ky > Nhalf) ky -= Ngrid;
      const double fy = ky * jfac;
      const double invsincy = (yi == 0) ? 1.0 : fy / sin(fy);

      /* Real-data FFT, so there is no second z half (symmetric) */
      for (int zi = 0; zi < (Nhalf + 1); ++zi) {

        const int kz = zi;
        const double fz = kz * jfac;
        const double invsincz = (zi == 0) ? 1.0 : fz / sin(fz);

        const int kk = kx * kx + ky * ky + kz * kz;

        /* Avoid singularity and don't go beyond Nyquist (otherwise we spend
           lots of time for little additional info */
        if (kk == 0 || kk > nyq2) continue;

        /* De-aliasing/deconvolution of mass assignment
           (Jing 2005 but without iterations) */
        double W = invsincx * invsincy * invsincz;

        /* W = W^windoworder */
        W = integer_pow(W, windoworder);

        /* Avoid sqrt with a lookup table kbin[kk] (does cost more memory)
         * bin = (int)(sqrt(kk) + 0.5); */
        const int bin = kbin[kk];

        const int mult = (zi == 0) ? 1 : 2;
        atomic_add(&modecounts[bin], mult);

        const int index = (xi * Ngrid + yi) * (Nhalf + 1) + zi;
        atomic_add_d(&powersum[bin],
                     mult * W * W *
                         (powgridft[index][0] * powgridft2[index][0] +
                          powgridft[index][1] * powgridft2[index][1]));
      } /* Loop over z */
    }   /* Loop over y */
  }     /* Loop over z */
}

/**
 * @brief Initialize a power spectrum output file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void power_init_output_file(FILE* fp, const enum power_type type1,
                                          const enum power_type type2,
                                          const struct unit_system* restrict us,
                                          const struct phys_const* phys_const) {

  /* Write a header to the output file */
  if (type1 != type2)
    fprintf(fp, "# %s-%s cross-spectrum\n", get_powtype_name(type1),
            get_powtype_name(type2));
  else
    fprintf(fp, "# %s power spectrum\n", get_powtype_name(type1));
  fprintf(fp, "######################################################\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0)  Redshift\n");
  fprintf(fp, "#      Unit = dimensionless\n");
  fprintf(fp, "# (1)  Fourier scale k\n");
  fprintf(fp, "#      Unit = %e cm**-1\n", pow(us->UnitLength_in_cgs, -1.));
  fprintf(fp, "#      Unit = %e pc**-1\n",
          pow(1. / phys_const->const_parsec, -1.));
  fprintf(fp, "#      Unit = %e Mpc**-1\n",
          pow(1. / phys_const->const_parsec / 1e6, -1.));
  fprintf(fp, "# (2)  The shot-noise subtracted power\n");
  if (type1 != pow_type_pressure && type2 != pow_type_pressure) {
    fprintf(fp, "#      Unit = %e cm**3\n", pow(us->UnitLength_in_cgs, 3.));
    fprintf(fp, "#      Unit = %e pc**3\n",
            pow(1. / phys_const->const_parsec, 3.));
    fprintf(fp, "#      Unit = %e Mpc**3\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  } else if (type1 == pow_type_pressure && type2 == pow_type_pressure) {
    fprintf(fp, "#     Unit = %e Mpc**3 * eV**2 * cm**-6\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  } else {
    fprintf(fp, "#     Unit = %e Mpc**3 * eV * cm**-3\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  }
  fprintf(fp, "# (3)  The shot-noise power (scale-independant)\n");
  if (type1 != pow_type_pressure && type2 != pow_type_pressure) {
    fprintf(fp, "#      Unit = %e cm**3\n", pow(us->UnitLength_in_cgs, 3.));
    fprintf(fp, "#      Unit = %e pc**3\n",
            pow(1. / phys_const->const_parsec, 3.));
    fprintf(fp, "#      Unit = %e Mpc**3\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  } else if (type1 == pow_type_pressure && type2 == pow_type_pressure) {
    fprintf(fp, "#     Unit = %e Mpc**3 * eV**2 * cm**-6\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  } else {
    fprintf(fp, "#     Unit = %e Mpc**3 * eV * cm**-3\n",
            pow(1. / phys_const->const_parsec / 1e6, 3.));
  }
  fprintf(fp, "#\n");
  fprintf(fp, "# %13s %15s %15s %15s\n", "(0)    ", "(1)    ", "(2)    ",
          "(3)    ");
  fprintf(fp, "# %13s %15s %15s %15s\n", "z", "k", "P(k)", "P_noise");
}

/**
 * @brief Compute the power spectrum between type1 and type2, including
 * foldings and dealiasing. Only the real part of the power is returned.
 *
 * Distributes the cells into slabs, performs the FFT, calculates the dealiased
 * power and writes this to file. Then repeats this procedure for a given
 * number of foldings.
 *
 * @param type1 The #enum power_type indicating the first quantity to correlate.
 * @param type2 The #enum power_type indicating the second quantity to
 * correlate.
 * @param powdata The #power_spectrum_data containing power spectrum
 * parameters, FFT plans and pointers to the grids.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
void power_spectrum(const enum power_type type1, const enum power_type type2,
                    struct power_spectrum_data* pow_data, const struct space* s,
                    struct threadpool* tp, const int verbose) {

  const int* local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const struct engine* e = s->e;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const int snapnum = e->ps_output_count; /* -1 if after snapshot dump */

  /* Extract some useful constants */
  const int Ngrid = pow_data->Ngrid;
  const int Ngrid2 = Ngrid * Ngrid;
  const int Ngrid3 = Ngrid2 * Ngrid;
  const int Nhalf = Ngrid / 2;
  const int Nfold = pow_data->Nfold;
  const int foldfac = pow_data->foldfac;
  const double jfac = M_PI / Ngrid;

  /* Gather some neutrino constants if using delta-f weighting on the mesh */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);

  if (verbose)
    message("Preparing to calculate the %s-%s power spectrum.",
            get_powtype_name(type1), get_powtype_name(type2));

  /* could loop over particles but for now just abort */
  if (nr_local_cells == 0)
    error("Cell infrastructure is not in place for power spectra.");

  /* Allocate the grids based on whether this is an auto- or cross-spectrum*/
  pow_data->powgrid = fftw_alloc_real(Ngrid2 * (Ngrid + 2));
  memuse_log_allocation("fftw_grid.grid", pow_data->powgrid, 1,
                        sizeof(double) * Ngrid2 * (Ngrid + 2));
  pow_data->powgridft = (fftw_complex*)pow_data->powgrid;
  if (type1 != type2) {
    pow_data->powgrid2 = fftw_alloc_real(Ngrid2 * (Ngrid + 2));
    memuse_log_allocation("fftw_grid.grid2", pow_data->powgrid2, 1,
                          sizeof(double) * Ngrid2 * (Ngrid + 2));
    pow_data->powgridft2 = (fftw_complex*)pow_data->powgrid2;
  } else {
    pow_data->powgrid2 = pow_data->powgrid;
    pow_data->powgridft2 = pow_data->powgridft;
  }

  /* Constants used for the normalization */
  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double volume = dim[0] * dim[1] * dim[2]; /* units Mpc^3 */
  const double Omega_m = s->e->cosmology->Omega_cdm + s->e->cosmology->Omega_b +
                         s->e->cosmology->Omega_nu_0;
  const double meanrho =
      Omega_m * s->e->cosmology->critical_density_0; /* 1e10 Msun/Mpc^3 */
  const double conv_EV = units_cgs_conversion_factor(us, UNIT_CONV_INV_VOLUME) /
                         phys_const->const_electron_volt;

  /* Inverse of the cosmic mean mass per grid cell in code units */
  double invcellmean, invcellmean2;

  if (type1 != pow_type_pressure)
    invcellmean = Ngrid3 / (meanrho * volume);
  else
    invcellmean = Ngrid3 / volume * conv_EV;

  if (type2 != pow_type_pressure)
    invcellmean2 = Ngrid3 / (meanrho * volume);
  else
    invcellmean2 = Ngrid3 / volume * conv_EV;

  /* When splitting the neutrino ensemble in half, double the inverse mean */
  if (type1 == pow_type_neutrino_0 || type1 == pow_type_neutrino_1)
    invcellmean *= 2.0;
  if (type2 == pow_type_neutrino_0 || type2 == pow_type_neutrino_1)
    invcellmean2 *= 2.0;

  /* Gather the shared information to be used by the threads
     for density computation */
  struct grid_mapper_data densdata;
  struct grid_mapper_data densdata2;

  densdata.cells = s->cells_top;
  densdata.dens = pow_data->powgrid;
  densdata.N = Ngrid;
  densdata.type = type1;
  densdata.windoworder = pow_data->windoworder;
  densdata.e = s->e;
  densdata.nu_model = &nu_model;
  if (type1 != type2) {
    densdata2.cells = s->cells_top;
    densdata2.dens = pow_data->powgrid2;
    densdata2.N = Ngrid;
    densdata2.type = type2;
    densdata2.windoworder = pow_data->windoworder;
    densdata2.e = s->e;
    densdata2.nu_model = &nu_model;
  }

  if (verbose) message("Calculating the shot noise.");

  /* Calculate mass terms for shot noise */
  double shot = 0.;
  if (type1 == pow_type_matter || type2 == pow_type_matter || type1 == type2 ||
      (type1 == pow_type_gas && type2 == pow_type_pressure) ||
      (type2 == pow_type_gas && type1 == pow_type_pressure)) {

    /* Note that for cross-power, there is only shot noise for particles
       that occur in both fields */
    struct shot_mapper_data shotdata;
    shotdata.cells = s->cells_top;
    shotdata.tot12 = 0;
    shotdata.type1 = type1;
    shotdata.type2 = type2;
    shotdata.e = s->e;
    shotdata.nu_model = &nu_model;
    threadpool_map(tp, shotnoise_mapper, (void*)local_cells, nr_local_cells,
                   sizeof(int), threadpool_auto_chunk_size, (void*)&shotdata);
#ifdef WITH_MPI
    /* Add up everybody's shot noise term */
    MPI_Allreduce(MPI_IN_PLACE, &shotdata.tot12, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    /* Store shot noise */
    shot = shotdata.tot12 / volume;
    if (type1 != pow_type_pressure)
      shot /= meanrho;
    else
      shot *= conv_EV;
    if (type2 != pow_type_pressure)
      shot /= meanrho;
    else
      shot *= conv_EV;
  }

  /* Gather the shared information to be used by the threads
     for density conversion */
  struct conv_mapper_data convdata;
  convdata.Ngrid = Ngrid;

  /* Create a lookup table for k (could also do this when initializing) */
  int* kbin = (int*)malloc((Nhalf * Nhalf + 1) * sizeof(int));
  for (int i = 0; i < Nhalf; ++i) {
    for (int j = 0; j <= i; ++j) kbin[i * i + j] = i;
    for (int j = i + 1; j <= 2 * i; ++j) kbin[i * i + j] = i + 1;
  }
  kbin[Nhalf * Nhalf] = Nhalf;

  /* Allocate arrays for power computation from FFT */
  int* modecounts = (int*)malloc((Nhalf + 1) * sizeof(int));
  double* powersum = (double*)malloc((Nhalf + 1) * sizeof(double));

  struct pow_mapper_data powmapdata;
  powmapdata.powgridft = pow_data->powgridft;
  powmapdata.powgridft2 = pow_data->powgridft2;
  powmapdata.Ngrid = Ngrid;
  powmapdata.windoworder = pow_data->windoworder;
  powmapdata.modecounts = modecounts;
  powmapdata.powersum = powersum;
  powmapdata.kbin = kbin;
  powmapdata.jfac = jfac;

  /* Allocate arrays for combined power spectrum */
  const int kcutn = (pow_data->windoworder >= 3) ? 90 : 70;
  const int kcutleft = (int)(Ngrid / 256.0 * kcutn);
  const int kcutright = (int)(Ngrid / 256.0 * (double)kcutn / foldfac);
  /* numtot = 76 * Nfold + 14;
   * assumes a 256 grid, foldfac=6 and  windoworder=2 */
  const int numtot = kcutleft + (Nfold - 1) * (kcutleft - kcutright + 1);
  int numstart = 0;

  double* kcomb = (double*)malloc(numtot * sizeof(double));
  double* pcomb = (double*)malloc(numtot * sizeof(double));

  /* Determine output file name */
  char outputfileBase[200] = "";
  char outputfileName[256] = "";

  sprintf(outputfileBase, "power_%s", get_powtype_filename(type1));
  if (type1 != type2) {
    const int length = strlen(outputfileBase);
    sprintf(outputfileBase + length, "-%s", get_powtype_filename(type2));
  }

  /* Loop over foldings */
  for (int i = 0; i < Nfold; ++i) {

    if (verbose) message("Calculating the power for folding num. %d.", i);

    /* Note:  implicitly assuming a cubic box here */
    densdata.fac = Ngrid / dim[0];
    densdata2.fac = Ngrid / dim[0];
    densdata.dim[0] = dim[0];
    densdata.dim[1] = dim[1];
    densdata.dim[2] = dim[2];
    densdata2.dim[0] = dim[0];
    densdata2.dim[1] = dim[1];
    densdata2.dim[2] = dim[2];
    const double kfac = 2 * M_PI / dim[0];

    /* Empty the grid(s) */
    bzero(pow_data->powgrid, Ngrid2 * (Ngrid + 2) * sizeof(double));
    if (type1 != type2)
      bzero(pow_data->powgrid2, Ngrid2 * (Ngrid + 2) * sizeof(double));

    /* Fill out the folded grid(s) */
    threadpool_map(tp, cell_to_powgrid_mapper, (void*)local_cells,
                   nr_local_cells, sizeof(int), threadpool_auto_chunk_size,
                   (void*)&densdata);
    if (type1 != type2)
      threadpool_map(tp, cell_to_powgrid_mapper, (void*)local_cells,
                     nr_local_cells, sizeof(int), threadpool_auto_chunk_size,
                     (void*)&densdata2);
#ifdef WITH_MPI
    /* Merge everybody's share of the grid onto rank 0 */
    if (e->nodeID == 0)
      MPI_Reduce(MPI_IN_PLACE, pow_data->powgrid, Ngrid2 * (Ngrid + 2),
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce(pow_data->powgrid, NULL, Ngrid2 * (Ngrid + 2), MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);

    /* Same for the secondary grid */
    if (type1 != type2) {

      if (e->nodeID == 0)
        MPI_Reduce(MPI_IN_PLACE, pow_data->powgrid2, Ngrid2 * (Ngrid + 2),
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      else
        MPI_Reduce(pow_data->powgrid2, NULL, Ngrid2 * (Ngrid + 2), MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    /* Only rank 0 needs to perform all the remaining work */
    if (e->nodeID == 0) {

      /* Convert mass to density contrast or pressure to eV/cm^3 */
      convdata.grid = pow_data->powgrid;
      convdata.invcellmean = invcellmean;
      if (Ngrid < 32) {
        mass_to_contrast_mapper(pow_data->powgrid, Ngrid, &convdata);
      } else {
        threadpool_map(tp, mass_to_contrast_mapper, pow_data->powgrid, Ngrid,
                       sizeof(double), threadpool_auto_chunk_size, &convdata);
      }

      if (type1 != type2) {
        convdata.grid = pow_data->powgrid2;
        convdata.invcellmean = invcellmean2;
        if (Ngrid < 32) {
          mass_to_contrast_mapper(pow_data->powgrid2, Ngrid, &convdata);
        } else {
          threadpool_map(tp, mass_to_contrast_mapper, pow_data->powgrid2, Ngrid,
                         sizeof(double), threadpool_auto_chunk_size, &convdata);
        }
      }

      /* Perform FFT(s) */
      fftw_execute_dft_r2c(pow_data->fftplanpow, pow_data->powgrid,
                           pow_data->powgridft);
      if (type1 != type2)
        fftw_execute_dft_r2c(pow_data->fftplanpow2, pow_data->powgrid2,
                             pow_data->powgridft2);

      powmapdata.powgridft = pow_data->powgridft;
      powmapdata.powgridft2 = pow_data->powgridft2;

      /* Zero the mode arrays */
      bzero(modecounts, (Nhalf + 1) * sizeof(int));
      bzero(powersum, (Nhalf + 1) * sizeof(double));

      /* Calculate compensated mode contributions */
      if (Ngrid < 32) {
        pow_from_grid_mapper(pow_data->powgridft, Ngrid, &powmapdata);
      } else {
        threadpool_map(tp, pow_from_grid_mapper, pow_data->powgridft, Ngrid,
                       sizeof(fftw_complex), threadpool_auto_chunk_size,
                       &powmapdata);
      }

      /* Write this folding to the detail file */
      const double volfac = (volume / Ngrid3) / Ngrid3;
      sprintf(outputfileName, "%s/%s_%04d_%d.txt", "power_spectra/foldings",
              outputfileBase, snapnum, i);
      FILE* outputfile = fopen(outputfileName, "w");

      /* Determine units of power */
      char powunits[32] = "";
      if (type1 != pow_type_pressure && type2 != pow_type_pressure)
        sprintf(powunits, "Mpc^3");
      else if (type1 == pow_type_pressure && type2 == pow_type_pressure)
        sprintf(powunits, "Mpc^3 (eV cm^(-3))^2");
      else
        sprintf(powunits, "Mpc^3 eV cm^(-3)");

      fprintf(outputfile,
              "# Folding %d, all lengths/volumes are comoving. k-bin centres "
              "are not corrected for the weights of the modes.\n",
              i);
      fprintf(outputfile, "# Shotnoise [%s]\n", powunits);
      fprintf(outputfile, "%g\n", shot);
      fprintf(outputfile, "# Redshift [dimensionless]\n");
      fprintf(outputfile, "%g\n", s->e->cosmology->z);
      fprintf(outputfile, "# k [Mpc^(-1)]   p [%s]\n", powunits);

      for (int j = 1; j <= Nhalf; ++j) {
        fprintf(outputfile, "%g %g\n", j * kfac,
                powersum[j] / modecounts[j] * volfac);
      }
      fclose(outputfile);

      /* Combine most accurate measurements from foldings */
      if (i == 0) {

        for (int j = 0; j < kcutleft; ++j) {
          kcomb[j] = (j + 1) * kfac;
          pcomb[j] = powersum[j + 1] / modecounts[j + 1] * volfac;
        }

        numstart += kcutleft;

      } else {

        const int off = kcutright + 1;
        for (int j = 0; j < (kcutleft - kcutright + 1); ++j) {
          kcomb[j + numstart] = (j + off) * kfac;
          pcomb[j + numstart] =
              powersum[j + off] / modecounts[j + off] * volfac;
        }
        numstart += (kcutleft - kcutright + 1);
      }

    } /* Work of rank 0 */

    /* Fold the box */
    for (int j = 0; j < 3; ++j) dim[j] /= foldfac;

  } /* Loop over the foldings */

  if (e->nodeID == 0) {

    /* Output attempt at combined measurement */
    sprintf(outputfileName, "%s/%s_%04d.txt", "power_spectra", outputfileBase,
            snapnum);

    FILE* outputfile = fopen(outputfileName, "w");

    /* Header and units */
    power_init_output_file(outputfile, type1, type2, us, phys_const);

    for (int j = 0; j < numtot; ++j) {

      float k = kcomb[j];

      /* Shall we correct the position of the k-space bin
       * to account for the different weights of the modes entering the bin? */
      if (pow_data->shift_centre_small_k_bins && j < number_of_corrected_bins) {
        k *= correction_shift_k_values[j];
      }

      fprintf(outputfile, "%15.8f %15.8e %15.8e %15.8e\n", s->e->cosmology->z,
              k, (pcomb[j] - shot), shot);
    }
    fclose(outputfile);
  }

  /* Done. Just clean up memory */
  if (verbose) message("Calculated the power! Cleaning up and leaving.");

  free(pcomb);
  free(kcomb);
  free(powersum);
  free(modecounts);
  free(kbin);
  if (type1 != type2) {
    memuse_log_allocation("fftw_grid.grid2", pow_data->powgrid2, 0, 0);
    fftw_free(pow_data->powgrid2);
  }
  pow_data->powgrid2 = NULL;
  pow_data->powgridft2 = NULL;
  memuse_log_allocation("fftw_grid.grid", pow_data->powgrid, 0, 0);
  fftw_free(pow_data->powgrid);
  pow_data->powgrid = NULL;
  pow_data->powgridft = NULL;
}

#endif /* HAVE_FFTW */

/**
 * @brief Initialize power spectra calculation.
 *
 * Reads in the power spectrum parameters, sets up FFTW
 * (likely already done for PM), then sets up an FFT plan
 * that can be reused every time.
 *
 * @param nr_threads The number of threads used.
 */
void power_init(struct power_spectrum_data* p, struct swift_params* params,
                int nr_threads) {

#ifdef HAVE_FFTW

  /* Get power spectrum parameters */
  p->Ngrid = parser_get_opt_param_int(params, "PowerSpectrum:grid_side_length",
                                      power_data_default_grid_side_length);
  p->Nfold = parser_get_param_int(params, "PowerSpectrum:num_folds");
  p->foldfac = parser_get_opt_param_int(params, "PowerSpectrum:fold_factor",
                                        power_data_default_fold_factor);
  p->windoworder = parser_get_opt_param_int(
      params, "PowerSpectrum:window_order", power_data_default_window_order);

  if (p->windoworder > 3 || p->windoworder < 1)
    error("Power spectrum calculation is not implemented for %dth order!",
          p->windoworder);
  if (p->windoworder == 1)
    message("WARNING: window order is recommended to be at least 2 (CIC).");
  if (p->windoworder <= 2 && p->foldfac > 4)
    message(
        "WARNING: fold factor is recommended not to exceed 4 for a "
        "mass assignment order of 2 (CIC) or below.");
  if (p->windoworder == 3 && p->foldfac > 6)
    message(
        "WARNING: fold factor is recommended not to exceed 6 for a "
        "mass assignment order of 3 (TSC) or below.");

  p->shift_centre_small_k_bins = parser_get_opt_param_int(
      params, "PowerSpectrum:shift_centre_small_k_bins", 1);

  /* Make sensible choices for the k-cuts */
  const int kcutn = (p->windoworder >= 3) ? 90 : 70;
  const int kcutleft = (int)(p->Ngrid / 256.0 * kcutn);
  const int kcutright = (int)(p->Ngrid / 256.0 * (double)kcutn / p->foldfac);
  if (kcutright < 10 || (kcutleft - kcutright) < 30)
    error(
        "Combination of power grid size and fold factor do not allow for "
        "enough overlap between foldings!");

  p->nr_threads = nr_threads;

#ifdef HAVE_THREADED_FFTW
  /* Initialise the thread-parallel FFTW version
     (probably already done for the PM, but does not matter) */
  if (p->Ngrid >= 64) {
    fftw_init_threads();
    fftw_plan_with_nthreads(nr_threads);
  }
#else
  message("Note that FFTW is not threaded!");
#endif

  char** requested_spectra = NULL;
  parser_get_param_string_array(params, "PowerSpectrum:requested_spectra",
                                &p->spectrumcount, &requested_spectra);

  p->types1 =
      (enum power_type*)malloc(p->spectrumcount * sizeof(enum power_type));
  p->types2 =
      (enum power_type*)malloc(p->spectrumcount * sizeof(enum power_type));

  /* Parse which spectra are being requested */
  for (int i = 0; i < p->spectrumcount; ++i) {

    char* pstr = strtok(requested_spectra[i], "-");
    if (pstr == NULL)
      error("Requested power spectra are not in the format type1-type2!");
    char type1[32];
    strcpy(type1, pstr);

    pstr = strtok(NULL, "-");
    if (pstr == NULL)
      error("Requested power spectra are not in the format type1-type2!");
    char type2[32];
    strcpy(type2, pstr);

    p->types1[i] = power_spectrum_get_type(type1);
    p->types2[i] = power_spectrum_get_type(type2);
  }

  /* Initialize the plan only once -- much faster for FFTs run often!
   * Does require us to allocate the grids, but we delete them right away.
   * Plan can only be used for the same FFTW call */
  const int Ngrid = p->Ngrid;

  /* Grid is padded to allow for in-place FFT */
  p->powgrid = fftw_alloc_real(Ngrid * Ngrid * (Ngrid + 2));
  /* Pointer to grid to interpret it as complex data */
  p->powgridft = (fftw_complex*)p->powgrid;

  p->fftplanpow = fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, p->powgrid,
                                       p->powgridft, FFTW_MEASURE);

  fftw_free(p->powgrid);
  p->powgrid = NULL;
  p->powgridft = NULL;

  /* Do the same for a second grid/plan to allow for cross power */

  /* Grid is padded to allow for in-place FFT */
  p->powgrid2 = fftw_alloc_real(Ngrid * Ngrid * (Ngrid + 2));
  /* Pointer to grid to interpret it as complex data */
  p->powgridft2 = (fftw_complex*)p->powgrid2;

  p->fftplanpow2 = fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, p->powgrid2,
                                        p->powgridft2, FFTW_MEASURE);

  fftw_free(p->powgrid2);
  p->powgrid2 = NULL;
  p->powgridft2 = NULL;

  /* Create directories for power spectra and foldings */
  if (engine_rank == 0) {
    safe_checkdir("power_spectra", /*create=*/1);
    safe_checkdir("power_spectra/foldings", /*create=*/1);
  }

#else
  error("Trying to initialize the PS code without FFTW present!");
#endif
}

void calc_all_power_spectra(struct power_spectrum_data* pow_data,
                            const struct space* s, struct threadpool* tp,
                            const int verbose) {
#ifdef HAVE_FFTW

  const ticks tic = getticks();

  /* Loop over all type combinations the user requested */
  for (int i = 0; i < pow_data->spectrumcount; ++i)
    power_spectrum(pow_data->types1[i], pow_data->types2[i], pow_data, s, tp,
                   verbose);

  /* Increment the PS output counter */
  s->e->ps_output_count++;

  /* Flag we did something special this step */
  s->e->step_props |= engine_step_prop_power_spectra;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("Can't use the PS code without FFTW present!");
#endif /* HAVE_FFTW */
}

void power_clean(struct power_spectrum_data* pow_data) {
#ifdef HAVE_FFTW
  fftw_destroy_plan(pow_data->fftplanpow);
  fftw_destroy_plan(pow_data->fftplanpow2);
  free(pow_data->types2);
  free(pow_data->types1);
#ifdef HAVE_THREADED_FFTW
  // Probably already done for PM at this point
  fftw_cleanup_threads();
#endif
#else
  error("Can't use the PS code without FFTW present!");
#endif /* HAVE_FFTW */
}

/**
 * @brief Write a power_spectrum_data struct to the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void power_spectrum_struct_dump(const struct power_spectrum_data* p,
                                FILE* stream) {
#ifdef HAVE_FFTW
  restart_write_blocks((void*)p, sizeof(struct power_spectrum_data), 1, stream,
                       "power spectrum data", "power spectrum data");
  restart_write_blocks(p->types1, p->spectrumcount, sizeof(enum power_type),
                       stream, "power types 1", "power types 1");
  restart_write_blocks(p->types2, p->spectrumcount, sizeof(enum power_type),
                       stream, "power types 2", "power types 2");
#endif
}

/**
 * @brief Restore a power_spectrum_data struct from the given FILE as a stream
 * of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void power_spectrum_struct_restore(struct power_spectrum_data* p,
                                   FILE* stream) {
#ifdef HAVE_FFTW
  restart_read_blocks((void*)p, sizeof(struct power_spectrum_data), 1, stream,
                      NULL, "power spectrum data");
  p->types1 =
      (enum power_type*)malloc(p->spectrumcount * sizeof(enum power_type));
  restart_read_blocks(p->types1, p->spectrumcount, sizeof(enum power_type),
                      stream, NULL, "power types 1");
  p->types2 =
      (enum power_type*)malloc(p->spectrumcount * sizeof(enum power_type));
  restart_read_blocks(p->types2, p->spectrumcount, sizeof(enum power_type),
                      stream, NULL, "power types 2");

#ifdef HAVE_THREADED_FFTW
  /* Initialise the thread-parallel FFTW version
     (probably already done for the PM, but does not matter) */
  if (p->Ngrid >= 64) {
    fftw_init_threads();
    fftw_plan_with_nthreads(p->nr_threads);
  }  // if
#else
  message("Note that FFTW is not threaded!");
#endif

  /* Initialize the plan only once -- much faster for FFTs run often!
     Does require us to allocate the grids, but we delete them right away */
  int Ngrid = p->Ngrid;

  /* Grid is padded to allow for in-place FFT */
  p->powgrid = fftw_alloc_real(Ngrid * Ngrid * (Ngrid + 2));
  /* Pointer to grid to interpret it as complex data */
  p->powgridft = (fftw_complex*)p->powgrid;

  p->fftplanpow = fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, p->powgrid,
                                       p->powgridft, FFTW_MEASURE);

  fftw_free(p->powgrid);
  p->powgrid = NULL;
  p->powgridft = NULL;

  /* Do the same for a second grid/plan to allow for cross power */

  /* Grid is padded to allow for in-place FFT */
  p->powgrid2 = fftw_alloc_real(Ngrid * Ngrid * (Ngrid + 2));
  /* Pointer to grid to interpret it as complex data */
  p->powgridft2 = (fftw_complex*)p->powgrid2;

  p->fftplanpow2 = fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, p->powgrid2,
                                        p->powgridft2, FFTW_MEASURE);

  fftw_free(p->powgrid2);
  p->powgrid2 = NULL;
  p->powgridft2 = NULL;
#endif /* HAVE_FFTW */
}
