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
  bzero(cell, sizeof(struct cell));
  struct part *part;
  struct xpart *xpart;
  float h;
  size_t x, y, z, size;

  size = count * sizeof(struct part);
  if (posix_memalign((void **)&cell->parts, part_align, size) != 0) {
    error("couldn't allocate particles");
  }

  size = count * sizeof(struct xpart);
  if (posix_memalign((void **)&cell->xparts, xpart_align, size) != 0) {
    error("couldn't allocate extended particles");
  }

  h = 1.2348 * cellSize / N;

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
        part->time_bin = 1;
        ++part;
      }
    }
  }

  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->gcount = 0;
  cell->dx_max_part = 0.;
  cell->dx_max_sort = 0.;
  cell->width[0] = cellSize;
  cell->width[1] = cellSize;
  cell->width[2] = cellSize;

  cell->ti_end_min = 1;
  cell->ti_end_max = 1;

  cell->sorted = 0;
  for (int k = 0; k < 13; k++) cell->sort[k] = NULL;

  return cell;
}

/* Just a forward declaration... */
void runner_doself1_density(struct runner *r, struct cell *ci);
void runner_doself2_force(struct runner *r, struct cell *ci);
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_dopair2_force(struct runner *r, struct cell *ci, struct cell *cj);

/* Run a full time step integration for one cell */
int main() {

#ifndef DEFAULT_SPH
  return 0;
#else

  int i, j, k, offset[3];
  struct part *p;

  int N = 10;
  float dim = 1.;
  float rho = 2.;
  float P = 1.;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Get some randomness going */
  srand(0);

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
      cells[j]->parts[i].u = P / (hydro_gamma_minus_one * rho);
    }

  message("m=%f", dim * dim * dim * rho / (N * N * N));

  /* Pick the central cell */
  struct cell *ci = cells[13];

  /* Create the infrastructure */
  struct space space;
  space.periodic = 0;
  space.cell_min = 1.;

  struct phys_const prog_const;
  prog_const.const_newton_G = 1.f;

  struct hydro_props hp;
  hp.target_neighbours = 48.f;
  hp.delta_neighbours = 2.;
  hp.max_smoothing_iterations = 1;
  hp.CFL_condition = 0.1;

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.hydro_properties = &hp;
  e.physical_constants = &prog_const;
  e.s = &space;
  e.time = 0.1f;
  e.ti_current = 1;

  struct runner r;
  r.e = &e;

  /* Simulation properties */
  e.timeBegin = 0.;
  e.timeEnd = 1.;
  e.timeOld = 0.;
  e.time = 0.1f;
  e.ti_current = 1;

  /* The tracked particle */
  p = &(ci->parts[N * N * N / 2 + N * N / 2 + N / 2]);

  message("Studying particle p->id=%lld", p->id);

  /* Sort the particles */
  for (j = 0; j < 27; ++j) {
    runner_do_sort(&r, cells[j], 0x1FFF, 0);
  }

  message("Sorting done");

  /* Initialise the particles */
  for (j = 0; j < 27; ++j) {
    runner_do_init(&r, cells[j], 0);
  }

  message("Init done");

  /* Compute density */
  runner_doself1_density(&r, ci);
  message("Self done");
  for (int j = 0; j < 27; ++j)
    if (cells[j] != ci) runner_dopair1_density(&r, ci, cells[j]);

  message("Density done");

  /* Ghost task */
  runner_do_ghost(&r, ci);

  message("h=%f rho=%f N_ngb=%f", p->h, p->rho, p->density.wcount);
  message("soundspeed=%f", p->force.soundspeed);

  runner_doself2_force(&r, ci);
  runner_do_kick(&r, ci, 1);

  message("ti_end=%d", p->ti_end);

  for (int j = 0; j < 27; ++j) {
    free(cells[j]->parts);
    free(cells[j]->xparts);
    for (int k = 0; k < 13; k++)
      if (cells[j]->sort[k] != NULL) free(cells[j]->sort[k]);
    free(cells[j]);
  }

  return 0;
#endif
}
