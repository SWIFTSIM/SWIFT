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
 * This function is used when running a non-periodic uniform box and just loops
 * over every other cell with particles and interacts with them.
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

  /* Get the mutlipole of the cell we are interacting. */
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

  /* Get the mutlipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Recover the list of top-level cells */
  struct cell *cells = e->s->cells_top;

#ifdef SWIFT_DEBUG_CHECKS
  if (top.type == cell_type_zoom) {
    error(
        "Zoom top cell found in long range gravity task! These should be "
        "handled "
        "by the void cell hierarchy.");
  }
#endif

  /* Since the zoom cells will be handled by the void cell hierarchy we can
   * just loop over all other cells which are not zoom cells. This is
   * trivial since the zoom cells are first in cells_top. */
  for (int cjd = s->zoom_props->bkg_cell_offset; cjd < s->nr_cells; cjd++) {

#ifdef SWIFT_DEBUG_CHECKS
    if (cells[cjd].type == cell_type_zoom) {
      error(
          "Zoom cell found in long range gravity task! These should be handled "
          "by the void cell hierarchy.");
    }
#endif

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[cjd];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (cj->subtype == cell_subtype_empty) {
      error("Empty cell found in long range gravity task!");
    }
#endif

    /* Avoid self contributions */
    if (top == cj) continue;
    runner_dopair_grav_mm_nonsym(r, ci, cj);
    continue;

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
 * @brief Performs M-M interactions between a given top-level cell
 *        (regular/background) and all other top level cells not interacted
 *        with via pair tasks or the mesh.
 *
 * This is used when the space is periodic and there is a mesh, therefore we
 * only interact with cells that are closer than the mesh interaction
 * distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * In the zoom case this function handles the following long range
 * interactions:
 * - bkg -> bkg
 * - bkg -> buffer (if their are buffer cells)
 *
 * Interactions with the zoom cells are not necessary since we defined the
 * interactions between zoom cells and any other cell at the top level of the
 * void cells.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_periodic(struct runner *r, struct cell *ci,
                                        struct cell *top, struct cell *cells) {

  struct engine *e = r->e;
  struct space *s = e->s;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Get the mutlipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Get the (i,j,k) location of the top-level cell in the grid. */
  int top_i = top->loc[0] * s->iwidth[0];
  int top_j = top->loc[1] * s->iwidth[1];
  int top_k = top->loc[2] * s->iwidth[2];

  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
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

  /* When running with a zoom region we also need to loop over useful buffer
   * cells. We can play the same game here walking out only a certain number
   * of buffer cells to interact with those we need to. */
  if (s->with_zoom_region && s->zoom_props->with_buffer_cells) {

    /* Get the buffer cell pointers. */
    struct cell *buffer_cells = s->zoom_props->buffer_cells_top;

    /* Get useful values. */
    const double *buffer_bounds = s->zoom_props->buffer_lower_bounds;
    const int *buffer_cdim = s->zoom_props->buffer_cdim;
    const double *buffer_iwidth = s->zoom_props->buffer_iwidth;

    /* Get the (i,j,k) location of the top-level buffer cell in the grid. */
    top_i = (top->loc[0] - buffer_bounds[0]) * buffer_iwidth[0];
    top_j = (top->loc[1] - buffer_bounds[1]) * buffer_iwidth[1];
    top_k = (top->loc[2] - buffer_bounds[2]) * buffer_iwidth[2];

    /* Maximal distance any interaction can take place
     * before the mesh kicks in, rounded up to the next integer */
    d = 1 + ceil(max_distance * max3(s->zoom_props->buffer_iwidth[0],
                                     s->zoom_props->buffer_iwidth[1],
                                     s->zoom_props->buffer_iwidth[2]));

    /* Loop over plausibly useful cells, exiting if beyond the buffer cells
     * which are not periodic. */
    for (int ii = top_i - d; ii <= top_i + d; ++ii) {
      if (ii < 0 || ii >= buffer_cdim[0]) continue;
      for (int jj = top_j - d; jj <= top_j + d; ++jj) {
        if (jj < 0 || jj >= buffer_cdim[1]) continue;
        for (int kk = top_k - d; kk <= top_k + d; ++kk) {
          if (kk < 0 || kk >= buffer_cdim[2]) continue;

          /* Get the cell */
          const int cjd = cell_getid(buffer_cdim, ii, jj, kk);

          /* Handle on the top-level cell */
          struct cell *cj = &buffer_cells[cjd];

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

            /* Call the PM interaction function on the active sub-cells of
             * ci */
            runner_dopair_grav_mm_nonsym(r, ci, cj);
            // runner_dopair_recursive_grav_pm(r, ci, cj);

            /* Record that this multipole received a contribution */
            multi_i->pot.interacted = 1;

          } /* We can interact with this cell */
        } /* Buffer cell i loop. */
      } /* Buffer cell j loop. */
    } /* Buffer cell k loop. */
  } /* Buffer cells enabled. */
}

/**
 * @brief Performs M-M interactions between a given top-level zoom cell and
 *        all other top level cells not interacted with via pair tasks or
 * the mesh.
 *
 * This is used when the space is periodic and there is a mesh, therefore we
 * only interact with cells that are closer than the mesh interaction
 * distance but further than the direct interaction distance.
 *
 * This function is only used when running a zoom simulation and handles the
 * following long range interactions:
 * - zoom -> nieghbour (where a neighbour can be either a background or
 * buffer)
 *
 * We only consider iteractions at the void->bkg/buffer level since we defined
 * all interactions at this level. Any interactions that were not pair tasks
 * below this level will have been handled by grav_mm tasks.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_long_range_zoom_periodic(struct runner *r, struct cell *ci,
                                        struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};

  /* Get the maximum distance at which we can have a non-mesh interaction.
   */
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Get the mutlipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* We need to test whether we need a long range interaction at the same
   * level we used to define pair interactions (i.e. the top level void
   * cells). Therefore, get the top level void cell. */
  struct cell *void_top = ci->top->void_parent->top;

  /* Because neighbours are by definition the background cells (or buffer
   * cells) within the mesh distance criterion we only need to check them
   * for long range interactions. */
  const int nr_neighbours = s->zoom_props->nr_neighbour_cells;
  const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

  /* Now loop over the neighbouring background cells.  */
  for (int k = 0; k < nr_neighbours; k++) {

    /* Handle on the neighbouring background cell. */
    struct cell *bkg_cj = &cells[neighbour_cells[k]];

    /* Handle on the top-level cell's gravity business*/
    const struct gravity_tensors *multi_j = bkg_cj->grav.multipole;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Minimal distance between any pair of particles */
    const double min_radius2 = cell_min_dist2(void_top, bkg_cj, periodic, dim);

    /* Are we beyond the distance where the truncated forces are 0 ?*/
    if (min_radius2 > max_distance2) {

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

      /* We are done here. */
      continue;
    }

    /* Shall we interact with this cell? */
    if (cell_can_use_pair_mm(void_top, bkg_cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* Call the PM interaction function on the active sub-cells of ci
       */
      runner_dopair_grav_mm_nonsym(r, ci, bkg_cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We can interact with this cell. */
  } /* Neighbour cell loop. */
}

/**
 * @brief Performs M-M interactions between a given top-level buffer cell
 * and all other top level cells not interacted with via pair tasks or the
 *        mesh.
 *
 * This is used when the space is periodic and there is a mesh, therefore we
 * only interact with cells that are closer than the mesh interaction
 * distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grids that are not covered by the mesh and some padding for
 * safety.
 *
 * This function is only used when running a zoom simulation and handles the
 * following long range interactions:
 * - buffer -> buffer
 * - buffer -> bkg
 *
 * Interactions with the zoom cells are not necessary since we defined the
 * interactions between zoom cells and any other cell at the top level of the
 * void cells.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_buffer_periodic(struct runner *r,
                                               struct cell *ci,
                                               struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};

  /* Get the maximum distance at which we can have a non-mesh interaction.
   */
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Get the mutlipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* For these we need to interact with all other buffer cells, all
   * background cells in range, and all zoom cells in range. */

  /* Get useful values. */
  const double *buffer_bounds = s->zoom_props->buffer_lower_bounds;
  const int *buffer_cdim = s->zoom_props->buffer_cdim;
  const double *buffer_iwidth = s->zoom_props->buffer_iwidth;
  struct cell *buffer_cells = s->zoom_props->buffer_cells_top;

  /* Get the (i,j,k) location of the top-level buffer cell in the grid. */
  int top_i = (top->loc[0] - buffer_bounds[0]) * buffer_iwidth[0];
  int top_j = (top->loc[1] - buffer_bounds[1]) * buffer_iwidth[1];
  int top_k = (top->loc[2] - buffer_bounds[2]) * buffer_iwidth[2];

  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
  int d = 1 + ceil(max_distance * max3(s->zoom_props->buffer_iwidth[0],
                                       s->zoom_props->buffer_iwidth[1],
                                       s->zoom_props->buffer_iwidth[2]));

  /* Loop over plausibly useful cells, exiting if beyond the buffer cells
   * which are not periodic. */
  for (int ii = top_i - d; ii <= top_i + d; ++ii) {
    if (ii < 0 || ii >= buffer_cdim[0]) continue;
    for (int jj = top_j - d; jj <= top_j + d; ++jj) {
      if (jj < 0 || jj >= buffer_cdim[1]) continue;
      for (int kk = top_k - d; kk <= top_k + d; ++kk) {
        if (kk < 0 || kk >= buffer_cdim[2]) continue;

        /* Get the cell */
        const int cjd = cell_getid(buffer_cdim, ii, jj, kk);

        /* Handle on the top-level cell */
        struct cell *cj = &buffer_cells[cjd];

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
        if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
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
      } /* Buffer cell i loop. */
    } /* Buffer cell j loop. */
  } /* Buffer cell k loop. */

  /* Finally we can interact with the background cells we need to by walking
   * out only as far as we need to from the empty cells above
   * this buffer cell. */

  /* Get the (i,j,k) location of the top-level cell in the grid. */
  top_i = top->loc[0] * s->iwidth[0];
  top_j = top->loc[1] * s->iwidth[1];
  top_k = top->loc[2] * s->iwidth[2];

  /* Get the background cells. */
  struct cell *bkg_cells = s->zoom_props->bkg_cells_top;

  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
  d = ceil(max_distance * max3(s->iwidth[0], s->iwidth[1], s->iwidth[2])) + 1;

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
        if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
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

/**
 * @brief Accumalate the number of particle mesh interactions for debugging
 * purposes.
 *
 * @param s The #space.
 * @param cells The top-level cells.
 * @param top The current top-level cell.
 * @param periodic Is the space periodic?
 * @param dim The dimensions of the space.
 * @param max_distance2 The maximum distance for a pair or mm interaction.
 */
void runner_count_mesh_interactions(struct runner *r, struct cell *ci,
                                    struct cell *top) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Get the mutlipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* We need to treat zoom cells and background/buffer cells differently
   * when counting interactions since we've made all interactions choices
   * at the void cell top level. Therefore, for zoom cells we need to use
   * the top level void cell for considerations. */
  if (ci->type == cell_type_zoom) {
    top = ci->top->void_parent->top;
  }

  /* Loop over all other non-zoom cells and account for the
   * mesh contribution. */
  for (int n = s->zoom_props->bkg_cell_offset; n < s->nr_cells; n++) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Minimal distance between any pair of particles */
    const double min_radius2 = cell_min_dist2(top, cj, periodic, dim);

    /* Are we beyond the distance where the truncated forces are 0 ?*/
    if (min_radius2 > max_distance2) {
#ifdef SWIFT_DEBUG_CHECKS
      /* Need to account for the interactions we missed */
      accumulate_add_ll(&multi_i->pot.num_interacted,
                        multi_j->m_pole.num_gpart);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
      /* Need to account for the interactions we missed */
      accumulate_add_ll(&multi_i->pot.num_interacted_pm,
                        multi_j->m_pole.num_gpart);
#endif
      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;
    }
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
 * This function is the main long-range gravity function. It is called by
 * the runner and will call the appropriate interaction function for the
 * given cell type/space periodicity.
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

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->subtype == cell_subtype_empty) {
    error("Empty cell found in long range gravity task!");
  }
#endif

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

  /* If we have a zoom cell the true top cell is a void cell. */
  if (top->type == cell_type_zoom) {
    top = top->void_parent->top;
  }

  /* Call the appropriate interaction function. */
  if (periodic) {
    switch (ci->type) {

      case cell_type_regular:
        runner_do_grav_long_range_periodic(r, ci, top, s->cells_top);
        break;
      case cell_type_zoom:
        runner_do_long_range_zoom_periodic(r, ci, top);
        break;
      case cell_type_buffer:
        runner_do_grav_long_range_buffer_periodic(r, ci, top);
        break;
      case cell_type_bkg:
        runner_do_grav_long_range_periodic(r, ci, top,
                                           s->zoom_props->bkg_cells_top);
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
  if (periodic) {
    runner_count_mesh_interactions(r, ci, top);
  }
#endif

  if (timer) TIMER_TOC(timer_dograv_long_range);
}
