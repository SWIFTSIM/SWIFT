/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2022 Will Roper (w.roper@sussex.ac.uk)
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

/* Includes. */
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/* Local includes. */
#include "cell.h"
#include "fof.h"
#include "halo_finder/halo.h"

/**
 * @brief Perform a FOF search using union-find on a given leaf-cell
 *
 * @param props The properties fof the FOF scheme.
 * @param l_x2 The square of the FOF linking length.
 * @param halo_level The type of halo we are finding (FOF group = 0,
 *                   6D Host = 1, 6D subhalo = 2)
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param c The #cell in which to perform FOF.
 */
void halo_finder_search_self_cell_gpart(const struct fof_props *props,
                                        const double l_x2,
                                        const enum halo_types halo_level,
                                        const struct gpart *const space_gparts,
                                        struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->split) error("Performing the Halo search at a non-leaf level!");
  if (halo_level == 0) error("Somehow we're in a halo finder function "
                             "and running a FOF!")
#endif

  /* Get particle counts and pointers to the particles. */
  const size_t count = c->grav.count;
  struct gpart *gparts = c->grav.parts;

  /* Index of particles in the global group list */
  size_t *group_index = props->group_index;
  const size_t group_id_offset = props->group_id_offset;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset = group_index + (ptrdiff_t)(gparts - space_gparts);

  if (c->nodeID != engine_rank)
    error("Performing self FOF search on foreign cell.");

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Ignore neutrinos */
    if (pi->type == swift_type_neutrino) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    /* Position of particle i. */
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Velocity of particle i. */
    const double pivx = pi->v_full[0];
    const double pivy = pi->v_full[1];
    const double pivz = pi->v_full[2];

    /* Get the mass of the FOF group/host halo of this particle. */
    double halo_mass;
    if (halo_level == 1) {
      halo_mass = pi->fof_data.group_mass;
    } else if (halo_level == 2) {
      halo_mass = pi->fof_data.host_mass;
    } else {
      error("Trying to find halos at a non-existent overdensity level.");
    }

    /* Define the velocity space linking length for the halo this particle
     * is in. */
    double l_v = props->ini_l_v_coeff * pow(halo_mass, 1.0 / 3.0);
    if (halo_level == host_halo)
      l_v *= props->const_l_v;
    else if (halo_level == sub_halo)
      l_v *= props->sub_const_l_v;
    double l_v2 = l_v * l_v;

    /* Find the root of pi. */
    size_t root_i = fof_find(offset[i], group_index);

    for (size_t j = i + 1; j < count; j++) {

      struct gpart *pj = &gparts[j];

      /* Ignore particles in different FOF groups */
      if (pi->hf_data.group_id != pj->hf_data.group_id) continue;

      /* If we're at the subhalo level ignore particles in different hosts */
      if (halo_level == 2 && (pi->hf_data.host_id != pj->hf_data.host_id))
        continue;

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Ignore neutrinos */
      if (pj->type == swift_type_neutrino) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      /* Position of particle j. */
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Velocity of particle j. */
      const double pjvx = pj->v_full[0];
      const double pjvy = pj->v_full[1];
      const double pjvz = pj->v_full[2];

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Compute the pairwise spatial distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      /* Compute the pairwise velocity "distance" */
      float dv[3], v2 = 0.0f;
      dv[0] = pivx - pjvx;
      dv[1] = pivy - pjvy;
      dv[2] = pivz - pjvz;

      for (int k = 0; k < 3; k++) {
        r2 += (dx[k] * dx[k]) / l_x2;
        v2 += (dv[k] * dv[k]) / l_v2;
      }

      /* Hit or miss? */
      if ((r2 + v2) < 2) {

        /* Merge the groups */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/**
 * @brief Perform a FOF search using union-find between two cells
 *
 * @param props The properties of the FOF scheme.
 * @param dim The dimension of the simulation volume.
 * @param l_x2 The square of the FOF linking length.
  * @param halo_level The type of halo we are finding (FOF group = 0,
 *                   6D Host = 1, 6D subhalo = 2)
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void halo_finder_search_pair_cells_gpart(const struct fof_props *props,
                                         const double dim[3],
                                         const double l_x2,
                                         const enum halo_types halo_level,
                                         const int periodic,
                                         const struct gpart *const space_gparts,
                                         struct cell *restrict ci,
                                         struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (halo_level == 0) error("Somehow we're in a halo finder function "
                             "and running a FOF!")
#endif

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = ci->grav.parts;
  struct gpart *gparts_j = cj->grav.parts;

  /* Index of particles in the global group list */
  size_t *group_index = props->group_index;
  const size_t group_id_offset = props->group_id_offset;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - space_gparts);
  size_t *const offset_j = group_index + (ptrdiff_t)(gparts_j - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (offset_j > offset_i && (offset_j < offset_i + count_i))
    error("Overlapping cells");
  if (offset_i > offset_j && (offset_i < offset_j + count_j))
    error("Overlapping cells");
#endif

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    struct gpart *pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Ignore neutrinos */
    if (pi->type == swift_type_neutrino) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Velocity of particle i. */
    const double pivx = pi->v_full[0];
    const double pivy = pi->v_full[1];
    const double pivz = pi->v_full[2];

    /* Get the mass of the FOF group/host halo of this particle. */
    double halo_mass;
    if (halo_level == host_halo) {
      halo_mass = pi->fof_data.group_mass;
    } else if (halo_level == sub_halo) {
      halo_mass = pi->fof_data.host_mass;
    } else {
      error("Trying to find halos at a non-existent overdensity level.");
    }

    /* Define the velocity space linking length for the halo this particle
     * is in. */
    double l_v = props->ini_l_v_coeff * pow(halo_mass, 1.0 / 3.0);
    if (halo_level == 1)
      l_v *= props->const_l_v;
    else if (halo_level == 2)
      l_v *= props->sub_const_l_v;
    double l_v2 = l_v * l_v;

    /* Find the root of pi. */
    size_t root_i = fof_find(offset_i[i], group_index);

    for (size_t j = 0; j < count_j; j++) {

      struct gpart *pj = &gparts_j[j];

      /* Ignore particles in different FOF groups */
      if (pi->hf_data.group_id != pj->hf_data.group_id) continue;

      /* If we're at the subhalo level ignore particles in different hosts */
      if (halo_level == 2 && (pi->hf_data.host_id != pj->hf_data.host_id))
        continue;

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Ignore neutrinos */
      if (pj->type == swift_type_neutrino) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset_j[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Velocity of particle j. */
      const double pjvx = pj->v_full[0];
      const double pjvy = pj->v_full[1];
      const double pjvz = pj->v_full[2];

      /* Compute the pairwise spatial distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      /* Compute the pairwise velocity "distance" */
      float dv[3], v2 = 0.0f;
      dv[0] = pivx - pjvx;
      dv[1] = pivy - pjvy;
      dv[2] = pivz - pjvz;

      for (int k = 0; k < 3; k++) {
        r2 += (dx[k] * dx[k]) / l_x2;
        v2 += (dv[k] * dv[k]) / l_v2;
      }

      /* Hit or miss? */
      if ((r2 + v2) < 2) {

        /* Merge the groups */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}
