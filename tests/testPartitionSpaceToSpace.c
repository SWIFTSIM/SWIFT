/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit.
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

#include <config.h>

/* Standard headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

#define NR_NODES 2

/**
 * @brief Set up a space of n^3 top-level cells in a unit box.
 */
static void make_space(struct space *s, struct engine *e, const int n) {

  bzero(e, sizeof(struct engine));
  bzero(s, sizeof(struct space));
  e->nr_nodes = NR_NODES;
  s->e = e;
  for (int d = 0; d < 3; d++) {
    s->dim[d] = 1.0;
    s->cdim[d] = n;
    s->width[d] = s->dim[d] / n;
    s->iwidth[d] = 1.0 / s->width[d];
  }
  s->nr_cells = n * n * n;
  s->cells_top = (struct cell *)calloc(s->nr_cells, sizeof(struct cell));
  if (s->cells_top == NULL) error("Failed to allocate the top-level cells");
}

/**
 * @brief Map an o^3 grid of owners onto an n^3 grid and check that every new
 * cell gets the owner of the old cell that holds its centre.
 *
 * The old owner is the old x index modulo NR_NODES, so a wrong old x index
 * gives a wrong owner.
 */
static void check_mapping(const int o, const int n) {

  struct engine e;
  struct space s;
  make_space(&s, &e, n);

  double oldh[3], oldcdim[3];
  for (int d = 0; d < 3; d++) {
    oldh[d] = 1.0 / o;
    oldcdim[d] = o;
  }

  int *old_owners = (int *)malloc(sizeof(int) * o * o * o);
  if (old_owners == NULL) error("Failed to allocate the old owners");
  for (int i = 0; i < o; i++)
    for (int j = 0; j < o; j++)
      for (int k = 0; k < o; k++)
        old_owners[(int)cell_getid(oldcdim, i, j, k)] = i % NR_NODES;

  const int complete = partition_space_to_space(oldh, oldcdim, old_owners, &s);

  for (int i = 0; i < n; i++) {
    const int expected = ((int)floor((i + 0.5) * o / (double)n)) % NR_NODES;
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        const int owner = s.cells_top[cell_getid(s.cdim, i, j, k)].nodeID;
        if (owner != expected)
          error("old %d new %d: cell (%d %d %d) has owner %d, expected %d", o,
                n, i, j, k, owner, expected);
      }
    }
  }
  if (!complete)
    error("old %d new %d: a node is reported missing but all are present", o,
          n);

  free(old_owners);
  free(s.cells_top);
}

/**
 * @brief Give one old cell an owner outside [0, NR_NODES) and check that the
 * mapping reports the partition as incomplete.
 */
static void check_bad_owner(const int bad_owner) {

  const int n = 4;
  struct engine e;
  struct space s;
  make_space(&s, &e, n);

  double oldh[3], oldcdim[3];
  for (int d = 0; d < 3; d++) {
    oldh[d] = 1.0 / n;
    oldcdim[d] = n;
  }

  int *old_owners = (int *)malloc(sizeof(int) * n * n * n);
  if (old_owners == NULL) error("Failed to allocate the old owners");
  for (int i = 0; i < n * n * n; i++) old_owners[i] = i % NR_NODES;
  old_owners[n * n * n - 1] = bad_owner;

  if (partition_space_to_space(oldh, oldcdim, old_owners, &s))
    error("owner %d outside [0, %d) was accepted", bad_owner, NR_NODES);

  free(old_owners);
  free(s.cells_top);
}

int main(int argc, char *argv[]) {

  /* Unchanged grid. */
  check_mapping(6, 6);

  /* Coarser grids. */
  check_mapping(8, 6);
  check_mapping(6, 5);

  /* Finer grids. Old 4 and new 9 also rules out the plain ratio i * o / n. */
  check_mapping(6, 7);
  check_mapping(6, 8);
  check_mapping(4, 9);

  /* Owners outside the valid range. */
  check_bad_owner(NR_NODES);
  check_bad_owner(-1);

  return 0;
}
