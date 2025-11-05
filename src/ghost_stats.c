/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include "ghost_stats.h"

#include "cell.h"
#include "space.h"

#ifdef SWIFT_GHOST_STATS

/**
 * @brief Resets all the individual cell ghost counters to 0.
 *
 * @param c The #cell to reset.
 */
void cell_reset_ghost_histograms(struct cell *c) {

  ghost_stats_reset_entries(&c->ghost_statistics);

  /* recurse */
  for (int k = 0; k < 8; ++k)
    if (c->progeny[k] != NULL) cell_reset_ghost_histograms(c->progeny[k]);
}

/**
 * @brief Write the ghost statistics for the given cell to the given file,
 * using the given cell ID.
 *
 * This function recurses into the progenies. They all get different cell IDs
 * by shifting the parent cell ID by 3 bits. Note that this does not
 * guarantee unique cell IDs (but it is good enough to distinguish parents
 * from their progeny).
 *
 * @param f File to write to.
 * @param c Cell to write.
 * @param cellID ID used to distinguish this cell from its progeny and
 * immediate neigbours (not guaranteed to be unique!).
 */
void cell_write_ghost_stats(FILE *f, const struct cell *c,
                            const long long cellID) {
  if (c == NULL) return;

  ghost_stats_write_cell_stats(f, &c->ghost_statistics, cellID);

  /* Write children */
  const long long pID = cellID << 3;
  for (int i = 0; i < 8; i++) {
    cell_write_ghost_stats(f, c->progeny[i], pID + i);
  }
}

/**
 * @brief Reset the ghost histograms for all top level cells in the space.
 *
 * @param s Space.
 */
void space_reset_ghost_histograms(struct space *s) {
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_reset_ghost_histograms(&s->cells_top[i]);
  }
}

/**
 * @brief Write the ghost statistics for all the top level cells in the space
 * to a file.
 *
 * The file name is composed as follows:
 *   ghost_stats_%04i_%04i.txt
 * where the first counter is an argument to this function (intended to be the
 * step counter), and the second counter is the rank that does the write.
 *
 * @param s Space.
 * @param j First counter in the output file name.
 */
void space_write_ghost_stats(const struct space *s, int j) {
  /* Open file */
  char filename[200];
  sprintf(filename, "ghost_stats_%04i_%04i.txt", j, engine_rank);
  FILE *f = fopen(filename, "w");
  if (f == NULL) error("Error opening ghost stats file.");

  /* Write header */
  ghost_stats_write_header(f);

  /* Write all the top level cells (and their children) */
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->nodeID == engine_rank) cell_write_ghost_stats(f, c, i);
  }

  /* Cleanup */
  fclose(f);
}

#endif
