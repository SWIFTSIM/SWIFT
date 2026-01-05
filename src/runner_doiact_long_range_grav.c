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

  /* Ensure we don't go out of bounds */
  if (d > s->cdim[0] / 2) d = s->cdim[0] / 2;
  if (d > s->cdim[1] / 2) d = s->cdim[1] / 2;
  if (d > s->cdim[2] / 2) d = s->cdim[2] / 2;

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
 * @brief Increment the mesh interaction counters.
 *
 * This is a helper function for incrementing the mesh interaction counters
 * for debugging purposes.
 *
 * @param multi_i The multipole receiving the interaction.
 * @param multi_j The multipole giving the interaction.
 */
static void runner_count_mesh_interaction(
    struct gravity_tensors *restrict multi_i,
    struct gravity_tensors *restrict multi_j) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Need to account for the mesh interactions we missed */
  accumulate_add_ll(&multi_i->pot.num_interacted, multi_j->m_pole.num_gpart);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Need to account for the mesh interactions we missed */
  accumulate_add_ll(&multi_i->pot.num_interacted_pm, multi_j->m_pole.num_gpart);
#endif
  /* Record that this multipole received a contribution */
  multi_i->pot.interacted = 1;
}

/**
 * @brief Recursively accumulate mesh interactions for pair interactions.
 *
 * This function mirrors the logic in scheduler_splittask_gravity for pair
 * tasks, recursing down the cell hierarchy and counting mesh interactions.
 *
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param cpi The current #cell from ci's hierarchy being processed.
 * @param cpj The current #cell from cj's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_pair_recursive(struct cell *ci,
                                                          struct cell *cpi,
                                                          struct cell *cpj,
                                                          struct space *s) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  /* No self interactions here */
  if (cpi == cpj) {
    error("Self interactions should not be handled in this function!");
  }

  struct engine *e = s->e;

  /* Handle on ci's gravity business. */
  struct gravity_tensors *multi_i = cpi->grav.multipole;
  struct gravity_tensors *multi_j = cpj->grav.multipole;

  /* Can we use the mesh for this pair? */
  if (engine_gravity_can_use_mesh(e, cpi, cpj)) {
    /* Record the mesh interaction */
    runner_count_mesh_interaction(multi_i, multi_j);
    return;
  }

  /* Should this pair be split? */
  if (cell_can_split_pair_gravity_task(cpi) &&
      cell_can_split_pair_gravity_task(cpj)) {

    /* Can we use M-M for this pair? */
    if (cell_can_use_pair_mm(cpi, cpj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/1)) {
      /* This would be handled by a M-M task, nothing to count */
      return;
    }

    /* We would create real tasks, so recurse to find mesh interactions */
    for (int i = 0; i < 8; i++) {
      if (cpi->progeny[i] == NULL) continue;
      /* Only recurse if this progeny is in a branch related to ci:
       * either the progeny contains ci (we're above ci) or ci contains the
       * progeny (we're at/below ci) */
      if (!cell_contains_progeny(cpi->progeny[i], ci) &&
          !cell_contains_progeny(ci, cpi->progeny[i]))
        continue;
      for (int j = 0; j < 8; j++) {
        if (cpj->progeny[j] == NULL) continue;
        runner_count_mesh_interactions_pair_recursive(ci, cpi->progeny[i],
                                                      cpj->progeny[j], s);
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */

#else
  error("This function should not be called without debugging checks enabled!");
#endif
}

/**
 * @brief Recursively accumulate mesh interactions for self interactions.
 *
 * This function mirrors the logic in scheduler_splittask_gravity for self
 * tasks, recursing down the cell hierarchy and counting mesh interactions.
 *
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param cpi The current #cell from ci's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_self_recursive(struct cell *ci,
                                                          struct cell *cpi,
                                                          struct space *s) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  struct engine *e = s->e;

  /* Handle on ci's gravity business. */
  struct gravity_tensors *multi_i = cpi->grav.multipole;

  /* Should this self task be split? */
  if (cell_can_split_self_gravity_task(cpi)) {

    /* Recurse on self interactions for each progeny */
    for (int k = 0; k < 8; k++) {
      if (cpi->progeny[k] == NULL) continue;
      if (!cell_contains_progeny(cpi->progeny[k], ci) &&
          !cell_contains_progeny(ci, cpi->progeny[k]))
        continue;
      runner_count_mesh_interactions_self_recursive(ci, cpi->progeny[k], s);
    }

    /* Now handle pair interactions between progeny */
    for (int j = 0; j < 8; j++) {
      if (cpi->progeny[j] == NULL) continue;
      if (!cell_contains_progeny(cpi->progeny[j], ci) &&
          !cell_contains_progeny(ci, cpi->progeny[j]))
        continue;
      struct cell *cpj = cpi->progeny[j];
      for (int k = j + 1; k < 8; k++) {
        if (cpi->progeny[k] == NULL) continue;
        if (!cell_contains_progeny(cpi->progeny[k], ci) &&
            !cell_contains_progeny(ci, cpi->progeny[k]))
          continue;
        struct cell *cpk = cpi->progeny[k];

        /* Can we use the mesh for this pair? */
        if (engine_gravity_can_use_mesh(e, cpj, cpk)) {
          /* Record the mesh interaction */
          runner_count_mesh_interaction(multi_i, cpk->grav.multipole);
          continue;
        }

        /* Otherwise recurse as a pair interaction */
        runner_count_mesh_interactions_pair_recursive(ci, cpj, cpk, s);
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */

#else
  error("This function should not be called without debugging checks enabled!");
#endif
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes.
 *
 * This function mirrors the task creation and splitting logic to count
 * mesh interactions that ci would receive from all top-level cells.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param top The current top-level cell (ci's top-level parent).
 */
void runner_count_mesh_interactions(struct runner *r, struct cell *ci,
                                    struct cell *top) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* First, handle self interactions from the top-level cell.
   * This mirrors the self task created at the top level. */
  runner_count_mesh_interactions_self_recursive(ci, top, s);

  /* Now loop over all other top-level cells for pair interactions.
   * This mirrors the pair tasks created between top-level cells. */
  for (int n = 0; n < s->nr_cells; n++) {

    /* Handle on the top-level cell and its gravity business */
    struct cell *cj = &cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions (already handled above) */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Can we use the mesh for this top-level pair? */
    if (engine_gravity_can_use_mesh(e, top, cj)) {

      /* If so, record the mesh interaction */
      runner_count_mesh_interaction(multi_i, multi_j);
      continue;
    }

    /* Can we use M-M for this top-level pair? */
    if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0)) {

      /* M-M task handles this, nothing to count */
      continue;
    }

    /* We would create a pair task here, so recurse to count mesh interactions
     * that arise from task splitting */
    runner_count_mesh_interactions_pair_recursive(ci, top, cj, s);
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
