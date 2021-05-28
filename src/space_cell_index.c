/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "threadpool.h"

/* Some standard headers. */
#include <float.h>

/**
 * @brief Information required to compute the particle cell indices.
 */
struct index_data {
  struct space *s;
  int *ind;
  int *cell_counts;
  size_t count_inhibited_part;
  size_t count_inhibited_gpart;
  size_t count_inhibited_spart;
  size_t count_inhibited_bpart;
  size_t count_inhibited_dmpart;
  size_t count_inhibited_sink;
  size_t count_extra_part;
  size_t count_extra_gpart;
  size_t count_extra_spart;
  size_t count_extra_bpart;
  size_t count_extra_dmpart;
  size_t count_extra_sink;
};

/**
 * @brief #threadpool mapper function to compute the particle cell indices.
 *
 * @param map_data Pointer towards the particles.
 * @param nr_parts The number of particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_parts_get_cell_index_mapper(void *map_data, int nr_parts,
                                       void *extra_data) {

  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(parts - s->parts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;
  size_t count_inhibited_part = 0;
  size_t count_extra_part = 0;

  /* Loop over the parts. */
  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    struct part *restrict p = &parts[k];

    double old_pos_x = p->x[0];
    double old_pos_y = p->x[1];
    double old_pos_z = p->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic && p->time_bin != time_bin_inhibited) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

    /* Put it back into the simulation volume */
    double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Treat the case where a particle was wrapped back exactly onto
     * the edge because of rounding issues (more accuracy around 0
     * than around dim) */
    if (pos_x == dim_x) pos_x = 0.0;
    if (pos_y == dim_y) pos_y = 0.0;
    if (pos_z == dim_z) pos_z = 0.0;

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    if (p->time_bin == time_bin_inhibited) {
      /* Is this particle to be removed? */
      ind[k] = -1;
      ++count_inhibited_part;
    } else if (p->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_part;

    } else {
      /* Normal case: list its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      /* Compute minimal mass */
      min_mass = min(min_mass, hydro_get_mass(p));

      /* Compute sum of velocity norm */
      sum_vel_norm += p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2];

      /* Update the position */
      p->x[0] = pos_x;
      p->x[1] = pos_y;
      p->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra parts */
  if (count_inhibited_part)
    atomic_add(&data->count_inhibited_part, count_inhibited_part);
  if (count_extra_part) atomic_add(&data->count_extra_part, count_extra_part);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_part_mass, min_mass);
  atomic_add_f(&s->sum_part_vel_norm, sum_vel_norm);
}

/**
 * @brief #threadpool mapper function to compute the g-particle cell indices.
 *
 * @param map_data Pointer towards the g-particles.
 * @param nr_gparts The number of g-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_gparts_get_cell_index_mapper(void *map_data, int nr_gparts,
                                        void *extra_data) {

  /* Unpack the data */
  struct gpart *restrict gparts = (struct gpart *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(gparts - s->gparts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;
  size_t count_inhibited_gpart = 0;
  size_t count_extra_gpart = 0;

  for (int k = 0; k < nr_gparts; k++) {

    /* Get the particle */
    struct gpart *restrict gp = &gparts[k];

    double old_pos_x = gp->x[0];
    double old_pos_y = gp->x[1];
    double old_pos_z = gp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic && gp->time_bin != time_bin_inhibited) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

    /* Put it back into the simulation volume */
    double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Treat the case where a particle was wrapped back exactly onto
     * the edge because of rounding issues (more accuracy around 0
     * than around dim) */
    if (pos_x == dim_x) pos_x = 0.0;
    if (pos_y == dim_y) pos_y = 0.0;
    if (pos_z == dim_z) pos_z = 0.0;

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    if (gp->time_bin == time_bin_inhibited) {
      /* Is this particle to be removed? */
      ind[k] = -1;
      ++count_inhibited_gpart;
    } else if (gp->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_gpart;

    } else {
      /* List its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      if (gp->type == swift_type_dark_matter) {

        /* Compute minimal mass */
        min_mass = min(min_mass, gp->mass);

        /* Compute sum of velocity norm */
        sum_vel_norm += gp->v_full[0] * gp->v_full[0] +
                        gp->v_full[1] * gp->v_full[1] +
                        gp->v_full[2] * gp->v_full[2];
      }

      /* Update the position */
      gp->x[0] = pos_x;
      gp->x[1] = pos_y;
      gp->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra gparts */
  if (count_inhibited_gpart)
    atomic_add(&data->count_inhibited_gpart, count_inhibited_gpart);
  if (count_extra_gpart)
    atomic_add(&data->count_extra_gpart, count_extra_gpart);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_gpart_mass, min_mass);
  atomic_add_f(&s->sum_gpart_vel_norm, sum_vel_norm);
}

/**
 * @brief #threadpool mapper function to compute the s-particle cell indices.
 *
 * @param map_data Pointer towards the s-particles.
 * @param nr_sparts The number of s-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_sparts_get_cell_index_mapper(void *map_data, int nr_sparts,
                                        void *extra_data) {

  /* Unpack the data */
  struct spart *restrict sparts = (struct spart *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(sparts - s->sparts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;
  size_t count_inhibited_spart = 0;
  size_t count_extra_spart = 0;

  for (int k = 0; k < nr_sparts; k++) {

    /* Get the particle */
    struct spart *restrict sp = &sparts[k];

    double old_pos_x = sp->x[0];
    double old_pos_y = sp->x[1];
    double old_pos_z = sp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic && sp->time_bin != time_bin_inhibited) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

    /* Put it back into the simulation volume */
    double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Treat the case where a particle was wrapped back exactly onto
     * the edge because of rounding issues (more accuracy around 0
     * than around dim) */
    if (pos_x == dim_x) pos_x = 0.0;
    if (pos_y == dim_y) pos_y = 0.0;
    if (pos_z == dim_z) pos_z = 0.0;

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    /* Is this particle to be removed? */
    if (sp->time_bin == time_bin_inhibited) {
      ind[k] = -1;
      ++count_inhibited_spart;
    } else if (sp->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_spart;

    } else {
      /* List its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      /* Compute minimal mass */
      min_mass = min(min_mass, sp->mass);

      /* Compute sum of velocity norm */
      sum_vel_norm +=
          sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2];

      /* Update the position */
      sp->x[0] = pos_x;
      sp->x[1] = pos_y;
      sp->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra sparts */
  if (count_inhibited_spart)
    atomic_add(&data->count_inhibited_spart, count_inhibited_spart);
  if (count_extra_spart)
    atomic_add(&data->count_extra_spart, count_extra_spart);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_spart_mass, min_mass);
  atomic_add_f(&s->sum_spart_vel_norm, sum_vel_norm);
}

/**
 * @brief #threadpool mapper function to compute the b-particle cell indices.
 *
 * @param map_data Pointer towards the b-particles.
 * @param nr_bparts The number of b-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_bparts_get_cell_index_mapper(void *map_data, int nr_bparts,
                                        void *extra_data) {

  /* Unpack the data */
  struct bpart *restrict bparts = (struct bpart *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(bparts - s->bparts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;
  size_t count_inhibited_bpart = 0;
  size_t count_extra_bpart = 0;

  for (int k = 0; k < nr_bparts; k++) {

    /* Get the particle */
    struct bpart *restrict bp = &bparts[k];

    double old_pos_x = bp->x[0];
    double old_pos_y = bp->x[1];
    double old_pos_z = bp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic && bp->time_bin != time_bin_inhibited) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

    /* Put it back into the simulation volume */
    double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Treat the case where a particle was wrapped back exactly onto
     * the edge because of rounding issues (more accuracy around 0
     * than around dim) */
    if (pos_x == dim_x) pos_x = 0.0;
    if (pos_y == dim_y) pos_y = 0.0;
    if (pos_z == dim_z) pos_z = 0.0;

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    /* Is this particle to be removed? */
    if (bp->time_bin == time_bin_inhibited) {
      ind[k] = -1;
      ++count_inhibited_bpart;
    } else if (bp->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_bpart;

    } else {
      /* List its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      /* Compute minimal mass */
      min_mass = min(min_mass, bp->mass);

      /* Compute sum of velocity norm */
      sum_vel_norm +=
          bp->v[0] * bp->v[0] + bp->v[1] * bp->v[1] + bp->v[2] * bp->v[2];

      /* Update the position */
      bp->x[0] = pos_x;
      bp->x[1] = pos_y;
      bp->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra bparts */
  if (count_inhibited_bpart)
    atomic_add(&data->count_inhibited_bpart, count_inhibited_bpart);
  if (count_extra_bpart)
    atomic_add(&data->count_extra_bpart, count_extra_bpart);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_bpart_mass, min_mass);
  atomic_add_f(&s->sum_bpart_vel_norm, sum_vel_norm);
}

/**
 * @brief #threadpool mapper function to compute the DM-particle cell indices.
 *
 * @param map_data Pointer towards the DM-particles.
 * @param nr_dmparts The number of DM-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_dmparts_get_cell_index_mapper(void *map_data, int nr_dmparts,
                                         void *extra_data) {

    /* Unpack the data */
    struct dmpart *restrict dmparts = (struct dmpart *)map_data;
    struct index_data *data = (struct index_data *)extra_data;
    struct space *s = data->s;
    int *const ind = data->ind + (ptrdiff_t)(dmparts - s->dmparts);

    /* Get some constants */
    const double dim_x = s->dim[0];
    const double dim_y = s->dim[1];
    const double dim_z = s->dim[2];
    const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
    const double ih_x = s->iwidth[0];
    const double ih_y = s->iwidth[1];
    const double ih_z = s->iwidth[2];

    /* Init the local count buffer. */
    int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
    if (cell_counts == NULL)
        error("Failed to allocate temporary cell count buffer.");

    /* Init the local collectors */
    float min_mass = FLT_MAX;
    float sum_vel_norm = 0.f;
    size_t count_inhibited_dmpart = 0;
    size_t count_extra_dmpart = 0;

    for (int k = 0; k < nr_dmparts; k++) {

        /* Get the particle */
        struct dmpart *restrict dmp = &dmparts[k];

        double old_pos_x = dmp->x[0];
        double old_pos_y = dmp->x[1];
        double old_pos_z = dmp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
        if (!s->periodic && dmp->time_bin != time_bin_inhibited) {
            if (old_pos_x < 0. || old_pos_x > dim_x)
                error("Particle outside of volume along X.");
            if (old_pos_y < 0. || old_pos_y > dim_y)
                error("Particle outside of volume along Y.");
            if (old_pos_z < 0. || old_pos_z > dim_z)
                error("Particle outside of volume along Z.");
        }
#endif

        /* Put it back into the simulation volume */
        double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
        double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
        double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

        /* Treat the case where a particle was wrapped back exactly onto
         * the edge because of rounding issues (more accuracy around 0
         * than around dim) */
        if (pos_x == dim_x) pos_x = 0.0;
        if (pos_y == dim_y) pos_y = 0.0;
        if (pos_z == dim_z) pos_z = 0.0;

        /* Get its cell index */
        const int index =
                cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
        if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
            error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
                  cdim[1], cdim[2], pos_x, pos_y, pos_z);

        if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
            pos_y < 0. || pos_z < 0.)
            error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
                  pos_z);
#endif

        /* Is this particle to be removed? */
        if (dmp->time_bin == time_bin_inhibited) {
            ind[k] = -1;
            ++count_inhibited_dmpart;
        } else if (dmp->time_bin == time_bin_not_created) {
            /* Is this a place-holder for on-the-fly creation? */
            ind[k] = index;
            cell_counts[index]++;
            ++count_extra_dmpart;

        } else {
            /* List its top-level cell index */
            ind[k] = index;
            cell_counts[index]++;

            /* Compute minimal mass */
            min_mass = min(min_mass, dmp->mass);

            /* Compute sum of velocity norm */
            sum_vel_norm +=
                    dmp->v_full[0] * dmp->v_full[0] + dmp->v_full[1] * dmp->v_full[1] + dmp->v_full[2] * dmp->v_full[2];

            /* Update the position */
            dmp->x[0] = pos_x;
            dmp->x[1] = pos_y;
            dmp->x[2] = pos_z;
        }
    }

    /* Write the counts back to the global array. */
    for (int k = 0; k < s->nr_cells; k++)
        if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
    free(cell_counts);

    /* Write the count of inhibited and extra dmparts */
    if (count_inhibited_dmpart)
        atomic_add(&data->count_inhibited_dmpart, count_inhibited_dmpart);
    if (count_extra_dmpart)
        atomic_add(&data->count_extra_dmpart, count_extra_dmpart);

    /* Write back the minimal part mass and velocity sum */
    atomic_min_f(&s->min_dmpart_mass, min_mass);
    atomic_add_f(&s->sum_dmpart_vel_norm, sum_vel_norm);
}


/**
 * @brief #threadpool mapper function to compute the sink-particle cell indices.
 *
 * @param map_data Pointer towards the sink-particles.
 * @param nr_sinks The number of sink-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_sinks_get_cell_index_mapper(void *map_data, int nr_sinks,
                                       void *extra_data) {

  /* Unpack the data */
  struct sink *restrict sinks = (struct sink *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(sinks - s->sinks);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;
  size_t count_inhibited_sink = 0;
  size_t count_extra_sink = 0;

  for (int k = 0; k < nr_sinks; k++) {

    /* Get the particle */
    struct sink *restrict sink = &sinks[k];

    double old_pos_x = sink->x[0];
    double old_pos_y = sink->x[1];
    double old_pos_z = sink->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic && sink->time_bin != time_bin_inhibited) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

    /* Put it back into the simulation volume */
    double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Treat the case where a particle was wrapped back exactly onto
     * the edge because of rounding issues (more accuracy around 0
     * than around dim) */
    if (pos_x == dim_x) pos_x = 0.0;
    if (pos_y == dim_y) pos_y = 0.0;
    if (pos_z == dim_z) pos_z = 0.0;

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    /* Is this particle to be removed? */
    if (sink->time_bin == time_bin_inhibited) {
      ind[k] = -1;
      ++count_inhibited_sink;
    } else if (sink->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_sink;

    } else {
      /* List its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      /* Compute minimal mass */
      min_mass = min(min_mass, sink->mass);

      /* Compute sum of velocity norm */
      sum_vel_norm += sink->v[0] * sink->v[0] + sink->v[1] * sink->v[1] +
                      sink->v[2] * sink->v[2];

      /* Update the position */
      sink->x[0] = pos_x;
      sink->x[1] = pos_y;
      sink->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra sinks */
  if (count_inhibited_sink)
    atomic_add(&data->count_inhibited_sink, count_inhibited_sink);
  if (count_extra_sink) atomic_add(&data->count_extra_sink, count_extra_sink);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_spart_mass, min_mass);
  atomic_add_f(&s->sum_spart_vel_norm, sum_vel_norm);
}

/**
 * @brief Computes the cell index of all the particles.
 *
 * Also computes the minimal mass of all #part.
 *
 * @param s The #space.
 * @param ind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_parts (return) The number of #part to remove.
 * @param count_extra_parts (return) The number of #part for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_parts_get_cell_index(struct space *s, int *ind, int *cell_counts,
                                size_t *count_inhibited_parts,
                                size_t *count_extra_parts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_part_mass = FLT_MAX;
  s->sum_part_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = ind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_inhibited_bpart = 0;
  data.count_inhibited_sink = 0;
  data.count_inhibited_dmpart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;
  data.count_extra_bpart = 0;
  data.count_extra_sink = 0;
  data.count_extra_dmpart = 0;

  threadpool_map(&s->e->threadpool, space_parts_get_cell_index_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                 &data);

  *count_inhibited_parts = data.count_inhibited_part;
  *count_extra_parts = data.count_extra_part;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the g-particles.
 *
 * Also computes the minimal mass of all dark-matter #gpart.
 *
 * @param s The #space.
 * @param gind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_gparts (return) The number of #gpart to remove.
 * @param count_extra_gparts (return) The number of #gpart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_gparts_get_cell_index(struct space *s, int *gind, int *cell_counts,
                                 size_t *count_inhibited_gparts,
                                 size_t *count_extra_gparts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_gpart_mass = FLT_MAX;
  s->sum_gpart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = gind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_inhibited_bpart = 0;
  data.count_inhibited_sink = 0;
  data.count_inhibited_dmpart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;
  data.count_extra_bpart = 0;
  data.count_extra_sink = 0;
  data.count_extra_dmpart = 0;

  threadpool_map(&s->e->threadpool, space_gparts_get_cell_index_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart),
                 threadpool_auto_chunk_size, &data);

  *count_inhibited_gparts = data.count_inhibited_gpart;
  *count_extra_gparts = data.count_extra_gpart;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the s-particles.
 *
 * Also computes the minimal mass of all #spart.
 *
 * @param s The #space.
 * @param sind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_sparts (return) The number of #spart to remove.
 * @param count_extra_sparts (return) The number of #spart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_sparts_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                 size_t *count_inhibited_sparts,
                                 size_t *count_extra_sparts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_spart_mass = FLT_MAX;
  s->sum_spart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = sind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_inhibited_sink = 0;
  data.count_inhibited_bpart = 0;
  data.count_inhibited_dmpart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;
  data.count_extra_bpart = 0;
  data.count_extra_sink = 0;
  data.count_extra_dmpart = 0;

  threadpool_map(&s->e->threadpool, space_sparts_get_cell_index_mapper,
                 s->sparts, s->nr_sparts, sizeof(struct spart),
                 threadpool_auto_chunk_size, &data);

  *count_inhibited_sparts = data.count_inhibited_spart;
  *count_extra_sparts = data.count_extra_spart;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the sink-particles.
 *
 * @param s The #space.
 * @param sink_ind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_sinks (return) The number of #sink to remove.
 * @param count_extra_sinks (return) The number of #sink for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_sinks_get_cell_index(struct space *s, int *sink_ind,
                                int *cell_counts, size_t *count_inhibited_sinks,
                                size_t *count_extra_sinks, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_sink_mass = FLT_MAX;
  s->sum_sink_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = sink_ind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_inhibited_bpart = 0;
  data.count_inhibited_sink = 0;
  data.count_inhibited_dmpart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;
  data.count_extra_bpart = 0;
  data.count_extra_sink = 0;
  data.count_extra_dmpart = 0;

  threadpool_map(&s->e->threadpool, space_sinks_get_cell_index_mapper, s->sinks,
                 s->nr_sinks, sizeof(struct sink), threadpool_auto_chunk_size,
                 &data);

  *count_inhibited_sinks = data.count_inhibited_sink;
  *count_extra_sinks = data.count_extra_sink;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the b-particles.
 *
 * Also computes the minimal mass of all #bpart.
 *
 * @param s The #space.
 * @param bind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_bparts (return) The number of #bpart to remove.
 * @param count_extra_bparts (return) The number of #bpart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_bparts_get_cell_index(struct space *s, int *bind, int *cell_counts,
                                 size_t *count_inhibited_bparts,
                                 size_t *count_extra_bparts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_bpart_mass = FLT_MAX;
  s->sum_bpart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = bind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_inhibited_bpart = 0;
  data.count_inhibited_sink = 0;
  data.count_inhibited_dmpart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;
  data.count_extra_bpart = 0;
  data.count_extra_sink = 0;
  data.count_extra_dmpart = 0;

  threadpool_map(&s->e->threadpool, space_bparts_get_cell_index_mapper,
                 s->bparts, s->nr_bparts, sizeof(struct bpart),
                 threadpool_auto_chunk_size, &data);

  *count_inhibited_bparts = data.count_inhibited_bpart;
  *count_extra_bparts = data.count_extra_bpart;

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the DM-particles.
 *
 * Also computes the minimal mass of all #dmpart.
 *
 * @param s The #space.
 * @param bind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param count_inhibited_dmparts (return) The number of #dmpart to remove.
 * @param count_extra_dmparts (return) The number of #dmpart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_dmparts_get_cell_index(struct space *s, int *dmind, int *cell_counts,
                                  size_t *count_inhibited_dmparts,
                                  size_t *count_extra_dmparts, int verbose) {

    const ticks tic = getticks();

    /* Re-set the counters */
    s->min_dmpart_mass = FLT_MAX;
    s->sum_dmpart_vel_norm = 0.f;

    /* Pack the extra information */
    struct index_data data;
    data.s = s;
    data.ind = dmind;
    data.cell_counts = cell_counts;
    data.count_inhibited_part = 0;
    data.count_inhibited_gpart = 0;
    data.count_inhibited_spart = 0;
    data.count_inhibited_bpart = 0;
    data.count_inhibited_dmpart = 0;
    data.count_inhibited_sink = 0;
    data.count_extra_part = 0;
    data.count_extra_gpart = 0;
    data.count_extra_spart = 0;
    data.count_extra_bpart = 0;
    data.count_extra_dmpart = 0;
    data.count_extra_sink = 0;

    threadpool_map(&s->e->threadpool, space_dmparts_get_cell_index_mapper,
                   s->dmparts, s->nr_dmparts, sizeof(struct dmpart),
                   threadpool_auto_chunk_size, &data);

    *count_inhibited_dmparts = data.count_inhibited_dmpart;
    *count_extra_dmparts = data.count_extra_dmpart;

    if (verbose)
        message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
                clocks_getunit());
}


