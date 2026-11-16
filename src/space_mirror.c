/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will Roper (w.roper@sussex.ac.uk)
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

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "cell_boundary.h"
#include "engine.h"
#include "kernel_hydro.h"
#include "memuse.h"

/**
 * @brief Make a mirror cell the reflection of a real cell hierarchy.
 *
 * The mirror cell gets the same tree as the cell it reflects, with the
 * progeny swapped along the reflected axes, and as many particles. The
 * particles of a mirror leaf are in the order of the leaf it reflects; those
 * of a split mirror cell are in the order of its own progeny, as for any
 * other cell.
 *
 * @param s The #space.
 * @param cm The mirror #cell. Its location, size and depth are already set.
 * @param c The real #cell it reflects.
 * @param parts The mirror particles of this cell.
 */
void space_mirror_clone(struct space *s, struct cell *cm, struct cell *c,
                        struct part *parts) {

  cm->mirror_of = c;
  cm->hydro.parts = parts;
  cm->hydro.xparts = NULL;
  cm->hydro.count = c->hydro.count;
  cm->hydro.count_total = c->hydro.count;

  /* Mirror particles are never updated, so mirror cells are never active.
   * The rest of the cell's state comes with the first refresh. */
  cm->hydro.ti_end_min = max_nr_timesteps;
  cm->hydro.ti_beg_max = 0;
  cm->maxdepth = c->maxdepth;
  cm->split = c->split;
  if (!c->split) return;

  /* The progeny index has one bit per axis: flip those we reflect in. */
  const unsigned int walls = cm->mirror_walls;
  const int flip = ((walls & 3u) ? 4 : 0) | ((walls & 12u) ? 2 : 0) |
                   ((walls & 48u) ? 1 : 0);

  /* ponytail: serial, so always thread 0's cell buffer. */
  space_getcells(s, 8, cm->progeny, /*tpid=*/0);
  for (int k = 0; k < 8; k++) {
    struct cell *cmp = cm->progeny[k];
    struct cell *cp = c->progeny[k ^ flip];
    if (cp == NULL) {
      space_recycle(s, cmp);
      cm->progeny[k] = NULL;
      continue;
    }

    cmp->width[0] = cm->width[0] / 2;
    cmp->width[1] = cm->width[1] / 2;
    cmp->width[2] = cm->width[2] / 2;
    cmp->loc[0] = cm->loc[0] + ((k & 4) ? cmp->width[0] : 0.);
    cmp->loc[1] = cm->loc[1] + ((k & 2) ? cmp->width[1] : 0.);
    cmp->loc[2] = cm->loc[2] + ((k & 1) ? cmp->width[2] : 0.);
    cmp->dmin = cm->dmin / 2;
    cmp->h_min_allowed = cmp->dmin * 0.5 * (1. / kernel_gamma);
    cmp->h_max_allowed = cmp->dmin * (1. / kernel_gamma);
    cmp->depth = cm->depth + 1;
    cmp->nodeID = cm->nodeID;
    cmp->type = cm->type;
    cmp->subtype = cm->subtype;
    cmp->mirror_walls = cm->mirror_walls;
    cmp->parent = cm;
    cmp->top = cm->top;
#ifdef WITH_MPI
    cmp->mpi.tag = -1;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
    cell_assign_cell_index(cmp, cm);
#endif

    space_mirror_clone(s, cmp, cp, parts);
    parts += cp->hydro.count;
  }
}

/**
 * @brief Fill the mirror cells with reflections of the cells next to the
 * mirror walls.
 *
 * Must be called once the real cells have been split. Does nothing in runs
 * without mirror walls.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_mirror_populate(struct space *s, int verbose) {

  if (s->boundary_conditions.wall_model != boundary_wall_model_mirror) return;

  const ticks tic = getticks();

  /* ponytail: only mirrors of local cells, and all of each wall cell; images
   * of foreign cells come with the MPI rank rule, trimming to the wall-side
   * progeny when the memory matters. */
  size_t nr_mirrors = 0;
  for (int k = 0; k < s->nr_cells; k++) {
    const struct cell *c = &s->cells_top[k];
    if (c->subtype == cell_subtype_mirror &&
        c->mirror_of->nodeID == engine_rank)
      nr_mirrors += c->mirror_of->hydro.count;
  }

  /* Make room for them */
  if (nr_mirrors > s->size_parts_mirror) {
    swift_free("parts_mirror", s->parts_mirror);
    s->size_parts_mirror = nr_mirrors * engine_redistribute_alloc_margin;
    if (swift_memalign("parts_mirror", (void **)&s->parts_mirror, part_align,
                       sizeof(struct part) * s->size_parts_mirror) != 0)
      error("Failed to allocate the mirror particles.");
  }
  s->nr_parts_mirror = nr_mirrors;

  /* Reflect each wall cell into its mirror cells */
  struct part *parts = s->parts_mirror;
  for (int k = 0; k < s->nr_cells; k++) {
    struct cell *cm = &s->cells_top[k];
    if (cm->subtype != cell_subtype_mirror) continue;
    struct cell *c = cm->mirror_of;
    if (c->nodeID != engine_rank || c->hydro.count == 0) continue;

    space_mirror_clone(s, cm, c, parts);
    cell_mirror_refresh(cm, s->dim, /*moved=*/1);
    parts += c->hydro.count;
  }

  if (verbose)
    message("Made %zd mirror particles, took %.3f %s.", nr_mirrors,
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}
