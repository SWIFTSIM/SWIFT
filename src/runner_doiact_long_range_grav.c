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
 * This function is used when running a non-periodic uniform box (i.e. non-zoom
 * simulation) and just loops over every other cell with particles and interacts
 * with them.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_uniform_non_periodic(struct runner *r,
                                                    struct cell *ci,
                                                    struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;

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
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* Call the PM interaction function on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  } /* Loop over top-level cells */
}

/**
 * @brief Performs M-M interactions between a given top-level cell and
 *        all other top level cells not interacted with via pair tasks.
 *
 * This is the non-periodic case where there is no mesh so all cells not
 * handled by a pair task are interacted with here in this long range
 * gravity function.
 *
 * This function is used when running a non-periodic zoom simulation and
 * will loop over all non-zoom cells but use the void cell hierarchy to
 * interact with all zoom cells.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_zoom_non_periodic(struct runner *r,
                                                 struct cell *ci,
                                                 struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Recover the list of top-level cells */
  struct cell *bkg_cells = e->s->zoom_props->bkg_cells_top;

#ifdef SWIFT_DEBUG_CHECKS
  /* Define counters used to count gparts. */
  size_t tested_gparts = 0;
#endif

  /* Since the zoom cells will be handled by the void cell hierarchy we can
   * just loop over all other cells which are not zoom cells. This is
   * trivial since the zoom cells are first in cells_top. */
  for (int cjd = 0; cjd < s->zoom_props->nr_bkg_cells; cjd++) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &bkg_cells[cjd];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

#ifdef SWIFT_DEBUG_CHECKS
    tested_gparts += multi_j->m_pole.num_gpart;
#endif

    /* Avoid self contributions */
    if (top == cj) continue;

    if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {
      /* Call the PM interaction function on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  } /* Loop over top-level cells */

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we at tested against all possible gparts. */
  if (tested_gparts != e->s->nr_gparts) {
    error(
        "Not all gparts were tested in long range gravity task! (tested: %ld, "
        "total: %ld)",
        tested_gparts, e->s->nr_gparts);
  }
#endif
}

/**
 * @brief Performs M-M interactions between a given top-level cell and all other
 * top level cells not interacted with via pair tasks or the mesh.
 *
 * This is used when running a zoom simulation, the space is periodic and there
 * is a mesh, therefore we only interact with cells that are closer than the
 * mesh interaction distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * This function handles the following long range interactions:
 * - zoom -> bkg
 * - zoom -> buffer (if their are buffer cells)
 * - bkg -> bkg
 * - bkg -> buffer (if their are buffer cells)
 *
 * Since we define the original pair tasks at the void level we only need to
 * consider combinations of background cells. However, If a zoom cell has a long
 * range task it means there were no tasks in the void level but we still only
 * need to consider interactions between the zoom cell and the background cells.
 * Note that the top cell pointer here is already going to be the top level void
 * parent cell if ci is a zoom cell.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_zoom_periodic(struct runner *r, struct cell *ci,
                                             struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *bkg_cells = s->zoom_props->bkg_cells_top;
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
        struct cell *cj = &bkg_cells[cell_index];

        /* Avoid self contributions  */
        if (top == cj) continue;

        /* Handle on the top-level cell's gravity business*/
        const struct gravity_tensors *multi_j = cj->grav.multipole;

        /* Skip empty cells */
        if (multi_j->m_pole.M_000 == 0.f) continue;

        /* Minimal distance between any pair of particles */
        const double min_radius2 = cell_min_dist2(top, cj, periodic, dim);

        /* Are we beyond the distance where the truncated forces are 0 ?*/
        if (min_radius2 > max_distance2) {

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

          /* We are done here. */
          continue;
        }

        /* Shall we interact with this cell? */
        if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/0,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {

          /* Call the PM interaction function on the active sub-cells of ci */
          runner_dopair_grav_mm_nonsym(r, ci, cj);

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

        } /* We can interact with this cell */
      } /* Loop over relevant top-level cells (k) */
    } /* Loop over relevant top-level cells (j) */
  } /* Loop over relevant top-level cells (i) */
}

/**
 * @brief Performs M-M interactions between a given top-level cell and all other
 * top level cells not interacted with via pair tasks or the mesh.
 *
 * This is used when running a uniform box, the space is periodic and there is a
 * mesh, therefore we only interact with cells that are closer than the mesh
 * interaction distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_uniform_periodic(struct runner *r,
                                                struct cell *ci,
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
        const double min_radius2 = cell_min_dist2(top, cj, periodic, dim);

        /* Are we beyond the distance where the truncated forces are 0 ?*/
        if (min_radius2 > max_distance2) {

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

          /* We are done here. */
          continue;
        }

        /* Shall we interact with this cell? */
        if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/0,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {

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

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

/**
 * @brief Increment the mesh interaction counters.
 *
 * This is a helper function for incrementing the mesh interaction counters
 * for debugging purposes.
 *
 * @param multi_i The multipole receiving the interaction.
 * @param multi_j The multipole giving the interaction.
 */
static void runner_accumulate_interaction(
    struct gravity_tensors *restrict multi_i,
    struct gravity_tensors *restrict multi_j) {

  /* Ensure we aren't self-interacting */
  if (multi_i == multi_j) {
    error("Self interactions should not be handled in this function!");
  }

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
 * @brief Count a mesh interaction between two related cells.
 *
 * Since the counts are accumulated downwards from the super level in the
 * grav/down task we need to update different cells based on the cells we have
 * been passed. This function will select what cell should be updated based on
 * the nesting:
 *   - If the super cell is nested within ci, then the super cell is updated.
 *   - If the super cell is ci, then they are both the same and it does not
 *     matter which is updated.
 *   - If ci is nested within the super cell, then ci is updated.
 *   - If we have a self interaction (pair nested within the same top cell):
 *     - If the super cell is nested within cj, then the super cell is updated.
 *     - If the super cell is cj, then they are both the same and it does not
 *       matter which is updated.
 *     - If cj is nested within the super cell, then cj is updated.
 *
 * @param super The super-cell being updated.
 * @param ci The first #cell in the pair.
 * @param cj The second #cell in the pair.
 */
static void runner_count_mesh_interaction(struct cell *super, struct cell *ci,
                                          struct cell *cj) {

  /* Do we share the same top level cell? i.e. are we self-interacting? */
  int is_self = ci->top == cj->top;

  /* Decide which cell we are updating. */
  if (super == ci) {
    runner_accumulate_interaction(super->grav.multipole, cj->grav.multipole);
  } else if (cell_contains_progeny(ci, super)) {
    runner_accumulate_interaction(super->grav.multipole, cj->grav.multipole);
  } else if (cell_contains_progeny(super, ci)) {
    runner_accumulate_interaction(ci->grav.multipole, cj->grav.multipole);
  }

  /* Handle the symmetric case for self interactions */
  if (is_self) {
    if (super == cj) {
      runner_accumulate_interaction(super->grav.multipole, ci->grav.multipole);
    } else if (cell_contains_progeny(cj, super)) {
      runner_accumulate_interaction(super->grav.multipole, ci->grav.multipole);
    } else if (cell_contains_progeny(super, cj)) {
      runner_accumulate_interaction(cj->grav.multipole, ci->grav.multipole);
    }
  }
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
static void runner_count_mesh_interactions_pair_recursive(struct cell *c,
                                                          struct cell *ci,
                                                          struct cell *cj,
                                                          struct space *s) {

  /* No self interactions here */
  if (ci == cj) {
    return;
  }
  if (c == cj) {
    return;
  }

  struct engine *e = s->e;

  /* Should this pair be split? */
  if (cell_can_split_pair_gravity_task(ci, cj)) {

    /* Check particle count threshold - mirrors scheduler_splittask_gravity */
    const long long gcount_i = ci->grav.count;
    const long long gcount_j = cj->grav.count;
    if (gcount_i * gcount_j < ((long long)space_subsize_pair_grav)) {
      return;
    }

    /* Recurse on all progeny pairs */
    for (int i = 0; i < 8; i++) {
      if (ci->progeny[i] == NULL) continue;
      struct cell *cpi = ci->progeny[i];
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] == NULL) continue;
        struct cell *cpj = cj->progeny[j];

        /* Can we use the mesh for this pair? */
        if (cell_can_use_mesh(e, cpi, cpj)) {
          /* Record the mesh interaction */
          runner_count_mesh_interaction(c, cpi, cpj);
          continue;
        }

        /* Can we use M-M for this pair? */
        if (cell_can_use_pair_mm(cpi, cpj, e, s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/1,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {
          /* This would be handled by a M-M task, nothing to count */
          continue;
        }

        /* We would create real tasks, so recurse to find mesh interactions */
        runner_count_mesh_interactions_pair_recursive(c, cpi, cpj, s);
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */
}

/**
 * @brief Recursively accumulate mesh interactions for self interactions.
 *
 * This function mirrors the logic in scheduler_splittask_gravity for self
 * tasks, recursing down the cell hierarchy and counting mesh interactions.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_self_recursive(struct cell *c,
                                                          struct cell *ci,
                                                          struct space *s) {

  struct engine *e = s->e;

  /* Should this self task be split? */
  if (cell_can_split_self_gravity_task(ci)) {

    /* Check particle count threshold - mirrors scheduler_splittask_gravity
     */
    if (ci->grav.count < space_subsize_self_grav) {
      return;
    }

    /* Recurse on self interactions for each progeny */
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] == NULL) continue;
      runner_count_mesh_interactions_self_recursive(c, ci->progeny[k], s);
    }

    /* Now handle pair interactions between progeny */
    for (int j = 0; j < 8; j++) {
      if (ci->progeny[j] == NULL) continue;
      struct cell *cpj = ci->progeny[j];
      for (int k = j + 1; k < 8; k++) {
        if (ci->progeny[k] == NULL) continue;
        struct cell *cpk = ci->progeny[k];

        /* Can we use the mesh for this pair? */
        if (cell_can_use_mesh(e, cpj, cpk)) {
          /* Record the mesh interaction */
          runner_count_mesh_interaction(c, cpj, cpk);
          continue;
        }

        /* Otherwise recurse as a pair interaction */
        runner_count_mesh_interactions_pair_recursive(c, cpj, cpk, s);
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes.
 *
 * This function mirrors the task creation and splitting logic to count
 * mesh interactions that ci would receive from all top-level cells.
 *
 * NOTE: This will recurse over cells that are not directly realted to c (the
 * super cell of ci). It will not add their contribution though so is "safe".
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param top The current top-level cell (ci's top-level parent).
 */
static void runner_count_mesh_interactions_uniform(struct runner *r,
                                                   struct cell *ci,
                                                   struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;

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
    if (cell_can_use_mesh(e, top, cj)) {

      /* If so, record the mesh interaction */
      runner_count_mesh_interaction(ci, top, cj);
      continue;
    }

    /* Can we use M-M for this top-level pair? */
    if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* M-M task handles this, nothing to count */
      continue;
    }

    /* We would create a pair task here, so recurse to count mesh interactions
     * that arise from task splitting */
    runner_count_mesh_interactions_pair_recursive(ci, top, cj, s);
  }
}

/**
 * @brief Recursively accumulate mesh interactions for zoom pair interactions.
 *
 * This function mirrors the logic in
 * zoom_scheduler_splittask_gravity_void_pair, recursing through the void cell
 * hierarchy and then using the normal pair recursive function once we reach
 * non-void cells.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param cj The current #cell being paired with.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_zoom_pair_recursive(
    struct cell *c, struct cell *ci, struct cell *cj, struct space *s) {

  struct engine *e = s->e;

  /* If neither cell is a void cell, use the normal pair recursive function */
  if (ci->subtype != cell_subtype_void && cj->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_pair_recursive(c, ci, cj, s);
    return;
  }

  /* Loop over progeny pairs, mirroring
   * zoom_scheduler_splittask_gravity_void_pair */
  for (int i = 0; i < 8; i++) {
    struct cell *cpi = ci->progeny[i];

    /* Skip NULL progeny */
    if (cpi == NULL) continue;

    /* Skip empty non-void progeny */
    if (cpi->grav.count == 0 && cpi->subtype != cell_subtype_void) continue;

    for (int j = 0; j < 8; j++) {
      struct cell *cpj = cj->progeny[j];

      /* Skip NULL progeny */
      if (cpj == NULL) continue;

      /* Skip empty non-void progeny */
      if (cpj->grav.count == 0 && cpj->subtype != cell_subtype_void) continue;

      /* Can we use the mesh for this pair? */
      if (cell_can_use_mesh(e, cpi, cpj)) {
        /* Record the mesh interaction */
        runner_count_mesh_interaction(c, cpi, cpj);
        continue;
      }

      /* Can we use M-M for this pair? */
      if (cell_can_use_pair_mm(cpi, cpj, e, s, /*use_rebuild_data=*/1,
                               /*is_tree_walk=*/0,
                               /*periodic boundaries*/ s->periodic,
                               /*use_mesh*/ s->periodic)) {
        /* M-M task handles this, nothing to count */
        continue;
      }

      /* Recurse to find more mesh interactions */
      runner_count_mesh_interactions_zoom_pair_recursive(c, cpi, cpj, s);
    }
  }
}

/**
 * @brief Recursively accumulate mesh interactions for zoom self interactions.
 *
 * This function mirrors the logic in
 * zoom_scheduler_splittask_gravity_void_self, recursing through the void cell
 * hierarchy and then using the normal self recursive function once we reach
 * non-void cells.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_zoom_self_recursive(
    struct cell *c, struct cell *ci, struct space *s) {

  struct engine *e = s->e;

  /* If not a void cell, use the normal self recursive function */
  if (ci->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_self_recursive(c, ci, s);
    return;
  }

  /* Loop over progeny for self interactions */
  for (int k = 0; k < 8; k++) {
    if (ci->progeny[k] == NULL) continue;

    /* Skip empty non-void progeny */
    if (ci->progeny[k]->subtype != cell_subtype_void &&
        ci->progeny[k]->grav.count == 0)
      continue;

    runner_count_mesh_interactions_zoom_self_recursive(c, ci->progeny[k], s);
  }

  /* Now handle pair interactions between progeny */
  for (int j = 0; j < 8; j++) {
    if (ci->progeny[j] == NULL) continue;

    /* Skip empty non-void progeny */
    if (ci->progeny[j]->subtype != cell_subtype_void &&
        ci->progeny[j]->grav.count == 0)
      continue;

    struct cell *cpj = ci->progeny[j];

    for (int k = j + 1; k < 8; k++) {
      if (ci->progeny[k] == NULL) continue;

      /* Skip empty non-void progeny */
      if (ci->progeny[k]->subtype != cell_subtype_void &&
          ci->progeny[k]->grav.count == 0)
        continue;

      struct cell *cpk = ci->progeny[k];

      /* Can we use the mesh for this pair? */
      if (cell_can_use_mesh(e, cpj, cpk)) {
        /* Record the mesh interaction */
        runner_count_mesh_interaction(c, cpj, cpk);
        continue;
      }

      /* Otherwise recurse as a pair interaction */
      runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpk, s);
    }
  }
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes in zoom simulations.
 *
 * This function mirrors the task creation and splitting logic for zoom
 * simulations to count mesh interactions that ci would receive from all
 * top-level cells.
 *
 * NOTE: This will recurse over cells that are not directly related to c (the
 * super cell of ci). It will not add their contribution though so is "safe".
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param top The current top-level cell (ci's top-level parent).
 */
static void runner_count_mesh_interactions_zoom(struct runner *r,
                                                struct cell *ci,
                                                struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *bkg_cells = s->zoom_props->bkg_cells_top;

  /* First, handle self interactions from the top-level cell.
   * This mirrors the self task created at the void level. */
  runner_count_mesh_interactions_zoom_self_recursive(ci, top, s);

  /* Now loop over all other background/void top-level cells for pair
   * interactions. This mirrors the pair tasks created between void cells. */
  for (int n = 0; n < s->zoom_props->nr_bkg_cells; n++) {

    /* Handle on the top-level cell and its gravity business */
    struct cell *cj = &bkg_cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions (already handled above) */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Can we use the mesh for this top-level pair? */
    if (cell_can_use_mesh(e, top, cj)) {

      /* If so, record the mesh interaction */
      runner_count_mesh_interaction(ci, top, cj);
      continue;
    }

    /* Can we use M-M for this top-level pair? */
    if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* M-M task handles this, nothing to count */
      continue;
    }

    /* We would create a pair task here, so recurse to count mesh interactions
     * that arise from task splitting through the void hierarchy */
    runner_count_mesh_interactions_zoom_pair_recursive(ci, top, cj, s);
  }
}

#endif /* SWIFT_DEBUG_CHECKS || SWIFT_GRAVITY_FORCE_CHECKS */

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

  /* If we have a nested cell the true top level cell where we defined the
   * interactions is the background void parent cell. */
  /* TODO: This is a bit of a hack, we should actually have the top pointers set
   * properly during void tree splitting. */
  while (top->void_parent != NULL) top = top->void_parent->top;

  /* Call the appropriate interaction function based on the type of the
   * cell in question. */
  if (periodic) {
    switch (ci->type) {

      case cell_type_regular:
        runner_do_grav_long_range_uniform_periodic(r, ci, top);
        break;
      case cell_type_zoom:
        runner_do_grav_long_range_zoom_periodic(r, ci, top);
        break;
      case cell_type_buffer:
        runner_do_grav_long_range_zoom_periodic(r, ci, top);
        break;
      case cell_type_bkg:
        runner_do_grav_long_range_zoom_periodic(r, ci, top);
        break;
      default:
        error("Unknown cell type in long-range gravity task!");
    }
  } else {

    switch (ci->type) {

      case cell_type_regular:
        runner_do_grav_long_range_uniform_non_periodic(r, ci, top);
        break;
      case cell_type_zoom:
        runner_do_grav_long_range_zoom_non_periodic(r, ci, top);
        break;
      case cell_type_buffer:
        runner_do_grav_long_range_zoom_non_periodic(r, ci, top);
        break;
      case cell_type_bkg:
        runner_do_grav_long_range_zoom_non_periodic(r, ci, top);
        break;
      default:
        error("Unknown cell type in long-range gravity task!");
    }
  }

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  /* Count the number of mesh interactions if using the mesh. */
  if (periodic && s->with_zoom_region) {
    runner_count_mesh_interactions_zoom(r, ci, top);
  } else if (periodic) {
    runner_count_mesh_interactions_uniform(r, ci, top);
  }
#endif

  if (timer) TIMER_TOC(timer_dograv_long_range);
}
