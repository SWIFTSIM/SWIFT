/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include "swift.h"

#include <stdlib.h>
#include <string.h>

/**
 * @brief Constructs a cell with N SPH particles
 */
struct cell *make_cell(size_t N, float cellSize, int offset[3], int id_offset) {
  size_t count = N * N * N;
  struct cell *cell = malloc(sizeof(struct cell));
  struct part *part;
  struct xpart *xpart;
  float h;
  size_t x, y, z, size;

  size = count * sizeof(struct part);
  if (posix_memalign((void **)&cell->parts, 32, size) != 0) {
    error("couldn't allocate particles");
  }

  size = count * sizeof(struct xpart);
  if (posix_memalign((void **)&cell->xparts, 32, size) != 0) {
    error("couldn't allocate extended particles");
  }

  h = 1.127 * cellSize / N;

  part = cell->parts;
  xpart = cell->xparts;
  memset(part, 0, count * sizeof(struct part));
  memset(xpart, 0, count * sizeof(struct xpart));
  for (x = 0; x < N; ++x) {
    for (y = 0; y < N; ++y) {
      for (z = 0; z < N; ++z) {
        part->x[0] =
            offset[0] * cellSize + x * cellSize / N + cellSize / (2 * N);
        part->x[1] =
            offset[1] * cellSize + y * cellSize / N + cellSize / (2 * N);
        part->x[2] =
            offset[2] * cellSize + z * cellSize / N + cellSize / (2 * N);
        part->h = h;
        part->id = x * N * N + y * N + z + id_offset;
        ++part;
      }
    }
  }

  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->h[0] = cellSize;
  cell->h[1] = cellSize;
  cell->h[2] = cellSize;

  return cell;
}

#ifdef DEFAULT_SPH

/* Just a forward declaration... */
void runner_doself1_density(struct runner *r, struct cell *ci);
void runner_doself2_force(struct runner *r, struct cell *ci);

/* Run a full time step integration for one cell */
int main() {

  int i, j, k, offset[3];
  struct part *p;

  int N = 10;
  float dim = 1.;
  float rho = 2.;
  float P = 1.;

  /* Create cells */
  struct cell *cells[27];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++) {
        offset[0] = i;
        offset[1] = j;
        offset[2] = k;
        cells[i * 9 + j * 3 + k] =
            make_cell(N, dim, offset, (i * 9 + j * 3 + k) * N * N * N);
      }

  /* Set particle properties */
  for (j = 0; j < 27; ++j)
    for (i = 0; i < cells[j]->count; ++i) {
      cells[j]->parts[i].mass = dim * dim * dim * rho / (N * N * N);
      cells[j]->parts[i].u = P / ((const_hydro_gamma - 1.) * rho);
    }

  message("m=%f", dim * dim * dim * rho / (N * N * N));

  /* Pick the central cell */
  struct cell *ci = cells[13];

  /* Create the infrastructure */
  struct engine e;
  struct runner r;
  r.e = &e;

  /* Simulation properties */
  e.timeBegin = 0.;
  e.timeEnd = 1.;
  e.timeOld = 0.;
  e.time = 0.;
  e.dt_min = 0.;
  e.dt_max = 1e10;

  /* The tracked particle */
  p = &(ci->parts[N * N * N / 2 + N * N / 2 + N / 2]);

  message("Studying particle p->id=%lld", p->id);

  /* Initialise the particles */
  for (j = 0; j < 27; ++j) {
    runner_doinit(&r, cells[j], 0);
  }

  /* Compute density */
  runner_doself1_density(&r, ci);
  runner_doghost(&r, ci);

  message("h=%f rho=%f N_ngb=%f", p->h, p->rho, p->density.wcount);
  message("c=%f", p->force.c);

  runner_doself2_force(&r, ci);
  runner_dokick(&r, ci, 1);

  message("ti_end=%d", p->ti_end);

  free(ci->parts);
  free(ci->xparts);

  return 0;
}

#else

int main() { return 0; }

#endif
