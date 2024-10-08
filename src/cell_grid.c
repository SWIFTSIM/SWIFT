/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *                             Yolan Uyttenhove (Yolan.Uyttenhove@UGent.be)
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

/* Corresponding header */
#include "cell_grid.h"

/* Local headers */
#include "engine.h"
#include "error.h"
#include "space_getsid.h"

/**
 * @brief Recursively free grid memory for cell.
 *
 * @param c The #cell.
 */
void cell_free_grid_rec(struct cell *c) {

#ifndef MOVING_MESH
  /* Nothing to do as we have no tessellations */
#else
#ifdef SWIFT_DEBUG_CHECKS
  if (c->grid.construction_level != c && c->grid.voronoi != NULL)
    error("Grid allocated, but not on grid construction level!");
#endif
  if (c->grid.construction_level == NULL) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_free_grid_rec(c->progeny[k]);

  } else if (c->grid.construction_level == c) {
    cell_free_grid(c);
  } else if (c->grid.construction_level != c) {
    error("Somehow ended up below grid construction level!");
  }
#endif
}

/**
 * @brief Updates the grid.self_completeness flag on this cell and its
 * sub-cells.
 *
 * This cell satisfies the local completeness criterion for the Voronoi grid.
 *
 * A cell is defined as locally complete if, when we would split that cell in
 * thirds along each dimension (i.e. in 27 smaller cells), every small cube
 * would contain at least one particle.
 *
 * If a given cell and its direct neighbours on the same level in the AMR tree
 * are self-complete, the Voronoi grid of that cell can completely be
 * constructed using only particles of this cell and its direct neighbours,
 * i.e. by doing a normal SWIFT neighbour loop.
 *
 * @param c The #cell to be checked
 * @param force Whether to forcefully recompute the completeness if it was not
 * invalidated.
 * */
void cell_grid_update_self_completeness(struct cell *c, int force) {
  if (c == NULL) return;
  if (!force && (c->grid.self_completeness != grid_invalidated_completeness))
    return;

  if (c->split) {
    int all_complete = 1;

    /* recurse */
    for (int i = 0; all_complete && i < 8; i++) {
      if (c->progeny[i] != NULL) {
        cell_grid_update_self_completeness(c->progeny[i], force);
        /* As long as all progeny is complete, this cell can safely be split for
         * the grid construction (when not considering neighbouring cells) */
        all_complete &=
            (c->progeny[i]->grid.self_completeness == grid_complete);
      }
    }

    /* If all sub-cells are complete, this cell is also complete. */
    if (all_complete) {
      c->grid.self_completeness = grid_complete;
      /* We set complete to true for now */
      c->grid.complete = 1;
      /* We are done here */
      return;
    }
  }

  /* If this cell is not split, or not all subcells are complete, we need to
   * check if this cell is complete by looping over all the particles. */

  /* criterion = 0b111_111_111_111_111_111_111_111_111*/
#ifdef HYDRO_DIMENSION_1D
  const int criterion = (1 << 3) - 1;
#elif defined(HYDRO_DIMENSION_2D)
  const int criterion = (1 << 9) - 1;
#elif defined(HYDRO_DIMENSION_3D)
  const int criterion = (1 << 27) - 1;
#else
#error "Unknown hydro dimension"
#endif
  int flags = 0;
  for (int i = 0; flags != criterion && i < c->hydro.count; i++) {
    struct part *p = &c->hydro.parts[i];
    int x_bin = (int)(3. * (p->x[0] - c->loc[0]) / c->width[0]);
    int y_bin = (int)(3. * (p->x[1] - c->loc[1]) / c->width[1]);
    int z_bin = (int)(3. * (p->x[2] - c->loc[2]) / c->width[2]);
    if (x_bin >= 0 && x_bin < 3 && y_bin >= 0 && y_bin < 3 && z_bin >= 0 &&
        z_bin < 3) {
      flags |= 1 << (x_bin + 3 * y_bin + 9 * z_bin);
    }
  }

  /* Set completeness flags accordingly */
  if (flags == criterion) {
    c->grid.self_completeness = grid_complete;
    c->grid.complete = 1;
  } else {
    c->grid.self_completeness = grid_incomplete;
    c->grid.complete = 0;
  }
}

/**
 * @brief (mapper) Updates the grid.self_completeness flag on this cell and its
 * sub-cells.
 *
 * This function recomputes the self_completeness flag for all cells containing
 * particles.
 */
void cell_grid_set_self_completeness_mapper(void *map_data, int num_elements,
                                            void *extra_data) {
  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {
    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *c = &cells[cid];

    /* A top level cell can be empty in 1D and 2D simulations, just skip it */
    if (c->hydro.count == 0) {
      continue;
#ifdef SWIFT_DEBUG_CHECKS
      if (hydro_dimension == 3)
        error("Found empty top-level cell while running in 3D!");
#endif
    }
    if (c->nodeID != nodeID) continue;

    /* Set the splittable attribute for the moving mesh */
    cell_grid_update_self_completeness(c, /*force*/ 1);
  }
}

/**
 * @brief Updates the grid.completeness flag for the given pair of cells and all
 * recursive pairs of sub-cells.
 *
 * If one of the cells of the pair is not self-complete, or requires a rebuild,
 * we mark the other cell in the pair as incomplete (it cannot construct its
 * Voronoi grid on that level).
 *
 *
 *
 * @param ci The first #cell of the pair to be checked.
 * @param cj The second #cell of the pair to be checked. If NULL, only check
 * pairs of subcells from #ci.
 * @param sid The sort_id (direction) of the pair.
 * @param e The #engine.
 * */
void cell_grid_set_pair_completeness(struct cell *restrict ci,
                                     struct cell *restrict cj, int sid,
                                     const struct engine *e) {

  int ci_local = ci->nodeID == e->nodeID;
  int cj_local = cj != NULL ? cj->nodeID == e->nodeID : 0;
  /* Anything to do here? */
  if (!ci_local && !cj_local) return;

  /* Self or pair? */
  if (cj == NULL) {
    /* Self: Here we just need to recurse to hit all the pairs of sub-cells */
    if (ci->split) {
      /* recurse */
      for (int k = 0; k < 7; k++) {
        if (ci->progeny[k] != NULL) {
          /* Self: Recurse for pairs of sub-cells of this sub-cell */
          cell_grid_set_pair_completeness(ci->progeny[k], NULL, 0, e);

          /* Recurse for pairs of sub-cells */
          for (int l = k + 1; l < 8; l++) {
            if (ci->progeny[l] != NULL) {
              /* Get sid for pair */
              int sid_sub = sub_sid_flag[k][l];
              /* Pair: Recurse for pairs of sub-cells of this pair of sub-cells
               */
              cell_grid_set_pair_completeness(ci->progeny[k], ci->progeny[l],
                                              sid_sub, e);
            }
          }
        }
      }
    }
  } else {
    /* pair: Here we need to recurse further AND check whether one of the
     * neighbouring cells invalidates the completeness of the other. */
    struct cell_split_pair pairs = cell_split_pairs[sid];
    if (ci->split && cj->split) {
      /* recurse */
      for (int i = 0; i < pairs.count; i++) {
        struct cell *ci_sub = ci->progeny[pairs.pairs[i].pid];
        struct cell *cj_sub = cj->progeny[pairs.pairs[i].pjd];
        if (ci_sub == NULL || cj_sub == NULL) continue;
        double shift[3];
        int sid_sub =
            space_getsid_and_swap_cells(e->s, &ci_sub, &cj_sub, shift);
#ifdef SWIFT_DEBUG_CHECKS
        assert(sid_sub == pairs.pairs[i].sid);
#endif
        cell_grid_set_pair_completeness(ci_sub, cj_sub, sid_sub, e);
      }
    } else if (!ci->split && cj->split) {
      /* Set the completeness for the sub-cells of cj for this sid to 0 (they
       * have no neighbouring cell on the same level for this SID) */
      for (int i = 0; i < pairs.count; i++) {
        int l = pairs.pairs[i].pjd;
        if (cj->progeny[l] != NULL) cj->progeny[l]->grid.complete = 0;
      }
    } else if (!cj->split && ci->split) {
      /* Set the completeness for the sub-cells of ci for this sid to 0 (they
       * have no neighbouring cell on the same level for this SID) */
      for (int i = 0; i < pairs.count; i++) {
        int k = pairs.pairs[i].pid;
        if (ci->progeny[k] != NULL) ci->progeny[k]->grid.complete = 0;
      }
    }

    /* Update these cells' completeness flags (i.e. check whether the
     * neighbouring cell invalidates completeness)
     * We need to use atomics here, since multiple threads may change this at
     * the same time. */
    if (ci_local) {
      atomic_and(&ci->grid.complete,
                 !cell_grid_pair_invalidates_completeness(ci, cj));
    }
    if (cj_local) {
      atomic_and(&cj->grid.complete,
                 !cell_grid_pair_invalidates_completeness(cj, ci));
    }
  }
}

/**
 * @brief (mapper) Sets the grid.completeness flag for all cells, by looping
 * aggregating the self_completeness flags of all neighbours of each cell.
 *
 * A cell is considered complete if it and all its neighbours are
 * self_complete. The Voronoi grid may be constructed at any level in the AMR
 * tree as long as the cell at that level is complete.
 * */
void cell_set_grid_completeness_mapper(void *map_data, int num_elements,
                                       void *extra_data) {
  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;
  const int periodic = e->s->periodic;

  struct space *s = e->s;
  const int nodeID = e->nodeID;
  const int *cdim = s->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL, to set
   * the neighbour flags. */
  for (int ind = 0; ind < num_elements; ind++) {
    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Anything to do here? */
    if (ci->hydro.count == 0) continue;

    const int ci_local = ci->nodeID == nodeID;

    /* Update completeness for all the pairs of sub cells of this cell */
    if (ci_local) cell_grid_set_pair_completeness(ci, NULL, 0, e);

    /* Now loop over all the neighbours of this cell to also update the
     * completeness for pairs with this cell and all pairs of sub-cells */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Treat pairs only once. */
          const int cj_local = cj->nodeID == nodeID;
          if (cid >= cjd || cj->hydro.count == 0 || (!ci_local && !cj_local))
            continue;

          /* Update the completeness flag for this pair of cells and all pair of
           * sub-cells */
          int sid = (kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1));
          const int flip = runner_flip[sid];
          sid = sortlistID[sid];
          if (flip) {
            cell_grid_set_pair_completeness(cj, ci, sid, e);
          } else {
            cell_grid_set_pair_completeness(ci, cj, sid, e);
          }
        }
      }
    } /* Now loop over all the neighbours of this cell */
  } /* Loop through the elements, which are just byte offsets from NULL. */
}

/**
 * @brief Select a suitable construction level for the Voronoi grid.
 *
 * The Voronoi grid is constructed at the lowest level for which the cell
 * has more than #space_grid_split_threshold hydro particles *and* is marked
 * as complete for the Voronoi construction.
 *
 * Cells below the construction level store a pointer to its higher level cell
 * at the construction level.
 *
 * @param c The #cell
 * @param construction_level NULL, if we are yet to encounter the suitable
 * construction level, or a pointer to the parent-cell of #c at the construction
 * level.
 */
void cell_set_grid_construction_level(struct cell *c,
                                      struct cell *construction_level) {

  /* Above construction level? */
  if (construction_level == NULL) {
    /* Check if we can split this cell (i.e. all sub-cells are complete) */
    int splittable = c->split && c->hydro.count > space_grid_split_threshold;
    if (!c->grid.complete) {
      /* Are we on the top level? */
      if (c->top == c) {
        warning("Found incomplete top level cell!");
        splittable = 0;
      } else {
        error("Found incomplete cell above construction level!");
      }
    }
    for (int k = 0; splittable && k < 8; k++) {
      if (c->progeny[k] != NULL) splittable &= c->progeny[k]->grid.complete;
    }

    if (!splittable) {
      /* This is the first time we encounter an unsplittable cell, meaning that
       * it has too few particles to be split further or one of its progenitors
       * is not complete. I.e. we have arrived at the construction level! */
      construction_level = c;
    }
  }

  /* Set the construction level of this cell */
  c->grid.construction_level = construction_level;

  /* Recurse. */
  if (c->split)
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        cell_set_grid_construction_level(c->progeny[k], construction_level);
    }
}

void cell_set_grid_construction_level_mapper(void *map_data, int num_elements,
                                             void *extra_data) {
  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL, to set
   * the neighbour flags. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Anything to do here? */
    if (ci->hydro.count == 0) continue;
    const int ci_local = ci->nodeID == nodeID;

    /* This cell's completeness flags are now set all the way down the cell
     * hierarchy. We can now set the construction level. */
    if (ci_local) {
      cell_set_grid_construction_level(ci, NULL);
    }
  }
}