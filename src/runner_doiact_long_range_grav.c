/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will Roper (w.roper@sussex.ac.uk)
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

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "runner.h"
#include "runner_doiact_grav.h"
#include "space.h"
#include "timers.h"

/**
 * @brief Performs M-M interactions between a given top-level cell and
 *        all other top level cells not interacted with via pair tasks.
 *
 * This is the non-periodic case where there is no mesh so all cells not
 * handled by a pair task are interacted with here in this long range
 * gravity function.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_non_periodic(struct runner *r, struct cell *ci,
                                            struct cell *top) {

  struct engine *e = r->e;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Recover the list of top-level cells */
  struct cell *cells = e->s->cells_top;
  int *cells_with_particles = e->s->cells_with_particles_top;
  const int nr_cells_with_particles = e->s->nr_cells_with_particles;

  /* Loop over all the top-level cells and go for a M-M interaction if
   * well-separated */
  for (int n = 0; n < nr_cells_with_particles; ++n) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[cells_with_particles[n]];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0)) {

      /* Call the PM interaction function on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  } /* Loop over top-level cells */
}

/**
 * @brief Performs M-M interactions between a given top-level cell and all other
 * top level cells not interacted with via pair tasks or the mesh.
 *
 * This is used when the space is periodic and there is a mesh, therefore we
 * only interact with cells that are closer than the mesh interaction distance
 * but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_periodic(struct runner *r, struct cell *ci,
                                        struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Get the (i,j,k) location of the top-level cell in the grid. */
  int top_i = top->loc[0] * s->iwidth[0];
  int top_j = top->loc[1] * s->iwidth[1];
  int top_k = top->loc[2] * s->iwidth[2];

  /* Maximal distance any interaction can take place before the mesh kicks in,
   * rounded up to the next integer */
  int d =
      ceil(max_distance * max3(s->iwidth[0], s->iwidth[1], s->iwidth[2])) + 1;

  /* Loop over plausibly useful cells */
  for (int ii = top_i - d; ii <= top_i + d; ++ii) {
    for (int jj = top_j - d; jj <= top_j + d; ++jj) {
      for (int kk = top_k - d; kk <= top_k + d; ++kk) {

        /* Box wrap */
        const int iii = (ii + s->cdim[0]) % s->cdim[0];
        const int jjj = (jj + s->cdim[1]) % s->cdim[1];
        const int kkk = (kk + s->cdim[2]) % s->cdim[2];

        /* Get the cell */
        const int cell_index = cell_getid(s->cdim, iii, jjj, kkk);

        /* Handle on the top-level cell */
        struct cell *cj = &cells[cell_index];

        /* Avoid self contributions  */
        if (top == cj) continue;

        /* Handle on the top-level cell's gravity business*/
        const struct gravity_tensors *multi_j = cj->grav.multipole;

        /* Skip empty cells */
        if (multi_j->m_pole.M_000 == 0.f) continue;

        /* Minimal distance between any pair of particles */
        const double min_radius2 =
            cell_min_dist2_same_size(top, cj, periodic, dim);

        /* Are we beyond the distance where the truncated forces are 0 ?*/
        if (min_radius2 > max_distance2) {

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

          /* We are done here. */
          continue;
        }

        /* Shall we interact with this cell? */
        if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/0)) {

          /* Call the PM interaction function on the active sub-cells of ci
           */
          runner_dopair_grav_mm_nonsym(r, ci, cj);
          // runner_dopair_recursive_grav_pm(r, ci, cj);

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

        } /* We can interact with this cell */
      } /* Loop over relevant top-level cells (k) */
    } /* Loop over relevant top-level cells (j) */
  } /* Loop over relevant top-level cells (i) */
}

/**
 * @brief Recurse accumulating mesh interactions.
 *
 * @param multi_i The #multipole whose counter we are updating.
 * @param ci The zoom #cell.
 * @param cj The background #cell.
 * @param s The #space.
 */
void runner_count_mesh_interactions_recursive(struct gravity_tensors *multi_i,
                                              struct cell *ci, struct cell *cj,
                                              struct space *s) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  struct engine *e = s->e;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Handle on the pair's gravity business. */
  struct gravity_tensors *multi_j = cj->grav.multipole;

  /* Are we beyond the mesh distance? */
  const double min_radius2 =
      cell_min_dist2_same_size(ci, cj, s->periodic, s->dim);
  if (min_radius2 > max_distance2) {
#ifdef SWIFT_DEBUG_CHECKS
    /* Need to account for the interactions we missed */
    accumulate_add_ll(&multi_i->pot.num_interacted, multi_j->m_pole.num_gpart);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
    /* Need to account for the interactions we missed */
    accumulate_add_ll(&multi_i->pot.num_interacted_pm,
                      multi_j->m_pole.num_gpart);
#endif

    /* Record that this multipole received a contribution */
    multi_i->pot.interacted = 1;

  }

  /* Ok, recurse down but don't go further than where the tasks are. */
  else if (ci->grav.super != ci) {
    for (int i = 0; i < 8; i++) {
      if (ci->progeny[i] == NULL) continue;
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] == NULL) continue;
        runner_count_mesh_interactions_recursive(multi_i, ci->progeny[i],
                                                 cj->progeny[j], s);
      }
    }
  }

#else
  error("This function should not be called without debugging checks enabled!");
#endif
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The current top-level cell.
 */
void runner_count_mesh_interactions(struct runner *r, struct cell *ci,
                                    struct cell *top) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Loop over all cells. */
  for (int n = 0; n < s->nr_cells; n++) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Avoid self contributions */
    if (top == cj) continue;

    /* Did we interact via the mesh with this pair or any of their progeny? */
    runner_count_mesh_interactions_recursive(multi_i, top, cj, s);
  }
#else
  error(
      "This function should not be called without debugging checks or "
      "force checks enabled!");
#endif
}

/**
 * @brief Performs all M-M interactions between a given top-level cell and
 * all the other top-levels that are far enough.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param timer Are we timing this ?
 */
void runner_do_grav_long_range(struct runner *r, struct cell *ci,
                               const int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;

  /* Is the space periodic? */
  const int periodic = s->periodic;

  /* Anything to do here? */
  if (!cell_is_active_gravity(ci, r->e)) return;

  if (ci->nodeID != engine_rank)
    error("Non-local cell in long-range gravity task!");

  /* Check multipole has been drifted */
  if (ci->grav.ti_old_multipole < r->e->ti_current)
    cell_drift_multipole(ci, r->e);

  /* Find this cell's top-level (great-)parent */
  struct cell *top = ci;
  while (top->parent != NULL) top = top->parent;

  /* Call the appropriate interaction function based on the type of the
   * cell in question. */
  if (periodic) {
    runner_do_grav_long_range_periodic(r, ci, top);
  } else {
    runner_do_grav_long_range_non_periodic(r, ci, top);
  }

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  /* Count the number of mesh interactions if using the mesh. */
  if (periodic) {
    runner_count_mesh_interactions(r, ci, top);
  }
#endif

  if (timer) TIMER_TOC(timer_dograv_long_range);
}
