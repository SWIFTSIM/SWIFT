/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "cell.h"

/**
 * @brief Pack the data of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param pc Pointer to an array of packed cells in which the
 *      cells will be packed.
 * @param with_gravity Are we running with gravity and hence need
 *      to exchange multipoles?
 *
 * @return The number of packed cells.
 */
int cell_pack(struct cell *restrict c, struct pcell *restrict pc,
              const int with_gravity) {
#ifdef WITH_MPI

  /* Start by packing the data of the current cell. */
  pc->hydro.h_max = c->hydro.h_max;
  pc->stars.h_max = c->stars.h_max;
  pc->black_holes.h_max = c->black_holes.h_max;
  pc->sinks.r_cut_max = c->sinks.r_cut_max;

  pc->hydro.ti_end_min = c->hydro.ti_end_min;
  pc->grav.ti_end_min = c->grav.ti_end_min;
  pc->stars.ti_end_min = c->stars.ti_end_min;
  pc->sinks.ti_end_min = c->sinks.ti_end_min;
  pc->black_holes.ti_end_min = c->black_holes.ti_end_min;

  pc->hydro.ti_old_part = c->hydro.ti_old_part;
  pc->grav.ti_old_part = c->grav.ti_old_part;
  pc->grav.ti_old_multipole = c->grav.ti_old_multipole;
  pc->stars.ti_old_part = c->stars.ti_old_part;
  pc->black_holes.ti_old_part = c->black_holes.ti_old_part;
  pc->sinks.ti_old_part = c->sinks.ti_old_part;

  pc->hydro.count = c->hydro.count;
  pc->grav.count = c->grav.count;
  pc->stars.count = c->stars.count;
  pc->sinks.count = c->sinks.count;
  pc->black_holes.count = c->black_holes.count;
  pc->maxdepth = c->maxdepth;

  /* Copy the Multipole related information */
  if (with_gravity) {
    const struct gravity_tensors *mp = c->grav.multipole;

    pc->grav.m_pole = mp->m_pole;
    pc->grav.CoM[0] = mp->CoM[0];
    pc->grav.CoM[1] = mp->CoM[1];
    pc->grav.CoM[2] = mp->CoM[2];
    pc->grav.CoM_rebuild[0] = mp->CoM_rebuild[0];
    pc->grav.CoM_rebuild[1] = mp->CoM_rebuild[1];
    pc->grav.CoM_rebuild[2] = mp->CoM_rebuild[2];
    pc->grav.r_max = mp->r_max;
    pc->grav.r_max_rebuild = mp->r_max_rebuild;
  }

#ifdef SWIFT_DEBUG_CHECKS
  pc->cellID = c->cellID;
#endif

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      pc->progeny[k] = count;
      count += cell_pack(c->progeny[k], &pc[count], with_gravity);
    } else {
      pc->progeny[k] = -1;
    }

  /* Return the number of packed cells used. */
  c->mpi.pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the tag of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param tags Pointer to an array of packed tags.
 *
 * @return The number of packed tags.
 */
int cell_pack_tags(const struct cell *c, int *tags) {
#ifdef WITH_MPI

  /* Start by packing the data of the current cell. */
  tags[0] = c->mpi.tag;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL)
      count += cell_pack_tags(c->progeny[k], &tags[count]);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->mpi.pcell_size != count) error("Inconsistent tag and pcell count!");
#endif  // SWIFT_DEBUG_CHECKS

  /* Return the number of packed tags used. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

void cell_pack_part_swallow(const struct cell *c,
                            struct black_holes_part_data *data) {

  const size_t count = c->hydro.count;
  const struct part *parts = c->hydro.parts;

  for (size_t i = 0; i < count; ++i) {
    data[i] = parts[i].black_holes_data;
  }
}

void cell_unpack_part_swallow(struct cell *c,
                              const struct black_holes_part_data *data) {

  const size_t count = c->hydro.count;
  struct part *parts = c->hydro.parts;

  for (size_t i = 0; i < count; ++i) {
    parts[i].black_holes_data = data[i];
  }
}

void cell_pack_bpart_swallow(const struct cell *c,
                             struct black_holes_bpart_data *data) {

  const size_t count = c->black_holes.count;
  const struct bpart *bparts = c->black_holes.parts;

  for (size_t i = 0; i < count; ++i) {
    data[i] = bparts[i].merger_data;
  }
}

void cell_unpack_bpart_swallow(struct cell *c,
                               const struct black_holes_bpart_data *data) {

  const size_t count = c->black_holes.count;
  struct bpart *bparts = c->black_holes.parts;

  for (size_t i = 0; i < count; ++i) {
    bparts[i].merger_data = data[i];
  }
}

/**
 * @brief Unpack the data of a given cell and its sub-cells.
 *
 * @param pc An array of packed #pcell.
 * @param c The #cell in which to unpack the #pcell.
 * @param s The #space in which the cells are created.
 * @param with_gravity Are we running with gravity and hence need
 *      to exchange multipoles?
 *
 * @return The number of cells created.
 */
int cell_unpack(struct pcell *restrict pc, struct cell *restrict c,
                struct space *restrict s, const int with_gravity) {
#ifdef WITH_MPI

  /* Unpack the current pcell. */
  c->hydro.h_max = pc->hydro.h_max;
  c->stars.h_max = pc->stars.h_max;
  c->black_holes.h_max = pc->black_holes.h_max;
  c->sinks.r_cut_max = pc->sinks.r_cut_max;

  c->hydro.ti_end_min = pc->hydro.ti_end_min;
  c->grav.ti_end_min = pc->grav.ti_end_min;
  c->stars.ti_end_min = pc->stars.ti_end_min;
  c->black_holes.ti_end_min = pc->black_holes.ti_end_min;
  c->sinks.ti_end_min = pc->sinks.ti_end_min;

  c->hydro.ti_old_part = pc->hydro.ti_old_part;
  c->grav.ti_old_part = pc->grav.ti_old_part;
  c->grav.ti_old_multipole = pc->grav.ti_old_multipole;
  c->stars.ti_old_part = pc->stars.ti_old_part;
  c->black_holes.ti_old_part = pc->black_holes.ti_old_part;
  c->sinks.ti_old_part = pc->sinks.ti_old_part;

  c->hydro.count = pc->hydro.count;
  c->grav.count = pc->grav.count;
  c->stars.count = pc->stars.count;
  c->sinks.count = pc->sinks.count;
  c->black_holes.count = pc->black_holes.count;
  c->maxdepth = pc->maxdepth;

#ifdef SWIFT_DEBUG_CHECKS
  c->cellID = pc->cellID;
#endif

  /* Copy the Multipole related information */
  if (with_gravity) {
    struct gravity_tensors *mp = c->grav.multipole;

    mp->m_pole = pc->grav.m_pole;
    mp->CoM[0] = pc->grav.CoM[0];
    mp->CoM[1] = pc->grav.CoM[1];
    mp->CoM[2] = pc->grav.CoM[2];
    mp->CoM_rebuild[0] = pc->grav.CoM_rebuild[0];
    mp->CoM_rebuild[1] = pc->grav.CoM_rebuild[1];
    mp->CoM_rebuild[2] = pc->grav.CoM_rebuild[2];
    mp->r_max = pc->grav.r_max;
    mp->r_max_rebuild = pc->grav.r_max_rebuild;
  }

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  c->split = 0;
  for (int k = 0; k < 8; k++)
    if (pc->progeny[k] >= 0) {
      struct cell *temp;
      space_getcells(s, 1, &temp);
      temp->hydro.count = 0;
      temp->grav.count = 0;
      temp->stars.count = 0;
      temp->loc[0] = c->loc[0];
      temp->loc[1] = c->loc[1];
      temp->loc[2] = c->loc[2];
      temp->width[0] = c->width[0] / 2;
      temp->width[1] = c->width[1] / 2;
      temp->width[2] = c->width[2] / 2;
      temp->dmin = c->dmin / 2;
      if (k & 4) temp->loc[0] += temp->width[0];
      if (k & 2) temp->loc[1] += temp->width[1];
      if (k & 1) temp->loc[2] += temp->width[2];
      temp->depth = c->depth + 1;
      temp->split = 0;
      temp->hydro.dx_max_part = 0.f;
      temp->hydro.dx_max_sort = 0.f;
      temp->stars.dx_max_part = 0.f;
      temp->stars.dx_max_sort = 0.f;
      temp->black_holes.dx_max_part = 0.f;
      temp->nodeID = c->nodeID;
      temp->parent = c;
      temp->top = c->top;
      c->progeny[k] = temp;
      c->split = 1;
      count += cell_unpack(&pc[pc->progeny[k]], temp, s, with_gravity);
    }

  /* Return the total number of unpacked cells. */
  c->mpi.pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the tags of a given cell and its sub-cells.
 *
 * @param tags An array of tags.
 * @param c The #cell in which to unpack the tags.
 *
 * @return The number of tags created.
 */
int cell_unpack_tags(const int *tags, struct cell *restrict c) {
#ifdef WITH_MPI

  /* Unpack the current pcell. */
  c->mpi.tag = tags[0];

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_tags(&tags[count], c->progeny[k]);
    }

#ifdef SWIFT_DEBUG_CHECKS
  if (c->mpi.pcell_size != count) error("Inconsistent tag and pcell count!");
#endif  // SWIFT_DEBUG_CHECKS

  /* Return the total number of unpacked tags. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the cell information about time-step sizes and displacements
 * of a cell hierarchy.
 *
 * @param c The #cells to pack.
 * @param pcells the packed cell structures to pack into.
 *
 * @return The number of cells that were packed.
 */
int cell_pack_end_step(const struct cell *c, struct pcell_step *pcells) {

#ifdef WITH_MPI

  /* Pack this cell's data. */
  pcells[0].hydro.ti_end_min = c->hydro.ti_end_min;
  pcells[0].hydro.dx_max_part = c->hydro.dx_max_part;

  pcells[0].grav.ti_end_min = c->grav.ti_end_min;

  pcells[0].stars.ti_end_min = c->stars.ti_end_min;
  pcells[0].stars.dx_max_part = c->stars.dx_max_part;

  pcells[0].black_holes.ti_end_min = c->black_holes.ti_end_min;
  pcells[0].black_holes.dx_max_part = c->black_holes.dx_max_part;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_end_step(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the cell information about time-step sizes and displacements
 * of a cell hierarchy.
 *
 * @param c The #cells to unpack into.
 * @param pcells the packed cell structures to unpack from.
 *
 * @return The number of cells that were packed.
 */
int cell_unpack_end_step(struct cell *c, const struct pcell_step *pcells) {

#ifdef WITH_MPI

  /* Unpack this cell's data. */
  c->hydro.ti_end_min = pcells[0].hydro.ti_end_min;
  c->hydro.dx_max_part = pcells[0].hydro.dx_max_part;

  c->grav.ti_end_min = pcells[0].grav.ti_end_min;

  c->stars.ti_end_min = pcells[0].stars.ti_end_min;
  c->stars.dx_max_part = pcells[0].stars.dx_max_part;

  c->black_holes.ti_end_min = pcells[0].black_holes.ti_end_min;
  c->black_holes.dx_max_part = pcells[0].black_holes.dx_max_part;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_end_step(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the hydro timebin information of the given cell.
 *
 * t needs to be aligned on SWIFT_CACHE_ALIGNMENT.
 *
 * @param c The #cell.
 * @param t (output) The array of timebins we pack into.
 */
void cell_pack_timebin(const struct cell *const c, timebin_t *const t) {

#ifdef WITH_MPI

  swift_declare_aligned_ptr(timebin_t, t_align, t, SWIFT_CACHE_ALIGNMENT);

  for (int i = 0; i < c->hydro.count; ++i)
    t_align[i] = c->hydro.parts[i].time_bin;

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Unpack the hydro timebin information of a given cell.
 *
 * t needs to be aligned on SWIFT_CACHE_ALIGNMENT.
 *
 * @param c The #cell
 * @param t The array of timebins we unpack from.
 */
void cell_unpack_timebin(struct cell *const c, timebin_t *const t) {

#ifdef WITH_MPI

  swift_declare_aligned_ptr(timebin_t, t_align, t, SWIFT_CACHE_ALIGNMENT);

  for (int i = 0; i < c->hydro.count; ++i)
    c->hydro.parts[i].time_bin = t_align[i];

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Pack the multipole information of the given cell and all it's
 * sub-cells.
 *
 * @param c The #cell.
 * @param pcells (output) The multipole information we pack into
 *
 * @return The number of packed cells.
 */
int cell_pack_multipoles(struct cell *restrict c,
                         struct gravity_tensors *restrict pcells) {
#ifdef WITH_MPI

  /* Pack this cell's data. */
  pcells[0] = *c->grav.multipole;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_multipoles(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the multipole information of a given cell and its sub-cells.
 *
 * @param c The #cell
 * @param pcells The multipole information to unpack
 *
 * @return The number of cells created.
 */
int cell_unpack_multipoles(struct cell *restrict c,
                           struct gravity_tensors *restrict pcells) {
#ifdef WITH_MPI

  /* Unpack this cell's data. */
  *c->grav.multipole = pcells[0];

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_multipoles(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the counts for star formation of the given cell and all it's
 * sub-cells.
 *
 * @param c The #cell.
 * @param pcells (output) The multipole information we pack into
 *
 * @return The number of packed cells.
 */
int cell_pack_sf_counts(struct cell *restrict c,
                        struct pcell_sf *restrict pcells) {

#ifdef WITH_MPI

  /* Pack this cell's data. */
  pcells[0].stars.delta_from_rebuild = c->stars.parts - c->stars.parts_rebuild;
  pcells[0].stars.count = c->stars.count;
  pcells[0].stars.dx_max_part = c->stars.dx_max_part;

  /* Pack this cell's data. */
  pcells[0].grav.delta_from_rebuild = c->grav.parts - c->grav.parts_rebuild;
  pcells[0].grav.count = c->grav.count;

#ifdef SWIFT_DEBUG_CHECKS
  /* Stars */
  if (c->stars.parts_rebuild == NULL)
    error("Star particles array at rebuild is NULL! c->depth=%d", c->depth);

  if (pcells[0].stars.delta_from_rebuild < 0)
    error("Stars part pointer moved in the wrong direction!");

  if (pcells[0].stars.delta_from_rebuild > 0 && c->depth == 0)
    error("Shifting the top-level pointer is not allowed!");

  /* Grav */
  if (c->grav.parts_rebuild == NULL)
    error("Grav. particles array at rebuild is NULL! c->depth=%d", c->depth);

  if (pcells[0].grav.delta_from_rebuild < 0)
    error("Grav part pointer moved in the wrong direction!");

  if (pcells[0].grav.delta_from_rebuild > 0 && c->depth == 0)
    error("Shifting the top-level pointer is not allowed!");
#endif

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_sf_counts(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the counts for star formation of a given cell and its
 * sub-cells.
 *
 * @param c The #cell
 * @param pcells The multipole information to unpack
 *
 * @return The number of cells created.
 */
int cell_unpack_sf_counts(struct cell *restrict c,
                          struct pcell_sf *restrict pcells) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.parts_rebuild == NULL)
    error("Star particles array at rebuild is NULL!");
  if (c->grav.parts_rebuild == NULL)
    error("Grav particles array at rebuild is NULL!");
#endif

  /* Unpack this cell's data. */
  c->stars.count = pcells[0].stars.count;
  c->stars.parts = c->stars.parts_rebuild + pcells[0].stars.delta_from_rebuild;
  c->stars.dx_max_part = pcells[0].stars.dx_max_part;

  c->grav.count = pcells[0].grav.count;
  c->grav.parts = c->grav.parts_rebuild + pcells[0].grav.delta_from_rebuild;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_sf_counts(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}
