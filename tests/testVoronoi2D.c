/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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
#include "../config.h"

/* Local headers. */
#include "hydro/Shadowswift/voronoi2d_algorithm.h"
#include "tools.h"

/* Number of cells used to test the 2D interaction algorithm */
#define TESTVORONOI2D_NUMCELL 100

int main(int argc, char *argv[]) {

  /* initialize simulation box */
  double anchor[3] = {-0.5f, -0.5f, -0.5f};
  double side[3] = {2.0f, 2.0f, 2.0f};

  /* test initialization and finalization algorithms */
  {
    message("Testing initialization and finalization algorithm...");

    struct voronoi_cell cell;
    double x[3] = {0.5, 0.5, 0.5};

    voronoi_cell_init(&cell, x, anchor, side);

    float maxradius = voronoi_cell_finalize(&cell);

    assert(maxradius == 2.0f * sqrtf(2.0f));

    assert(cell.volume == 4.0f);

    assert(cell.centroid[0] == 0.5f);
    assert(cell.centroid[1] == 0.5f);

    message("Done.");
  }

  /* test interaction algorithm: normal case */
  {
    message("Testing %i cell grid with random positions...",
            TESTVORONOI2D_NUMCELL);

    /* create 100 cells with random positions in [0,1]x[0,1] */
    struct voronoi_cell cells[TESTVORONOI2D_NUMCELL];
    double x[2];
    float dx[2];
    int i, j;
    float Atot;
    struct voronoi_cell *cell_i, *cell_j;

    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
      x[0] = random_uniform(0., 1.);
      x[1] = random_uniform(0., 1.);
      voronoi_cell_init(&cells[i], x, anchor, side);
#ifdef VORONOI_VERBOSE
      message("cell[%i]: %g %g", i, x[0], x[1]);
#endif
    }

    /* interact the cells (with periodic boundaries) */
    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
      cell_i = &cells[i];
      for (j = 0; j < TESTVORONOI2D_NUMCELL; ++j) {
        if (i != j) {
          cell_j = &cells[j];
          dx[0] = cell_i->x[0] - cell_j->x[0];
          dx[1] = cell_i->x[1] - cell_j->x[1];
          /* periodic boundaries */
          if (dx[0] >= 0.5f) {
            dx[0] -= 1.0f;
          }
          if (dx[0] < -0.5f) {
            dx[0] += 1.0f;
          }
          if (dx[1] >= 0.5f) {
            dx[1] -= 1.0f;
          }
          if (dx[1] < -0.5f) {
            dx[1] += 1.0f;
          }
#ifdef VORONOI_VERBOSE
          message("Cell %i before:", i);
          voronoi_print_cell(&cells[i]);
          message("Interacting cell %i with cell %i (%g %g, %g %g", i, j,
                  cells[i].x[0], cells[i].x[1], cells[j].x[0], cells[j].x[1]);
#endif
          voronoi_cell_interact(cell_i, dx, j);
        }
      }
    }

    Atot = 0.0f;
    /* print the cells to the stdout */
    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
#ifdef VORONOI_VERBOSE
      printf("Cell %i:\n", i);
      voronoi_print_cell(&cells[i]);
#endif
      voronoi_cell_finalize(&cells[i]);
      Atot += cells[i].volume;
    }

    /* Check the total surface area */
    assert(fabs(Atot - 1.0f) < 1.e-6);

    /* Check the neighbour relations for an arbitrary cell: cell 44
       We plotted the grid and manually found the correct neighbours and their
       order. */
    assert(cells[44].nvert == 7);
    assert(cells[44].ngbs[0] == 26);
    assert(cells[44].ngbs[1] == 38);
    assert(cells[44].ngbs[2] == 3);
    assert(cells[44].ngbs[3] == 33);
    assert(cells[44].ngbs[4] == 5);
    assert(cells[44].ngbs[5] == 90);
    assert(cells[44].ngbs[6] == 4);

    message("Done.");
  }

  /* test interaction algorithm */
  {
    message("Testing 100 cell grid with Cartesian mesh positions...");

    struct voronoi_cell cells[100];
    double x[2];
    float dx[2];
    int i, j;
    float Atot;
    struct voronoi_cell *cell_i, *cell_j;

    for (i = 0; i < 10; ++i) {
      for (j = 0; j < 10; ++j) {
        x[0] = (i + 0.5f) * 0.1;
        x[1] = (j + 0.5f) * 0.1;
        voronoi_cell_init(&cells[10 * i + j], x, anchor, side);
      }
    }

    /* interact the cells (with periodic boundaries) */
    for (i = 0; i < 100; ++i) {
      cell_i = &cells[i];
      for (j = 0; j < 100; ++j) {
        if (i != j) {
          cell_j = &cells[j];
          dx[0] = cell_i->x[0] - cell_j->x[0];
          dx[1] = cell_i->x[1] - cell_j->x[1];
          /* periodic boundaries */
          if (dx[0] >= 0.5f) {
            dx[0] -= 1.0f;
          }
          if (dx[0] < -0.5f) {
            dx[0] += 1.0f;
          }
          if (dx[1] >= 0.5f) {
            dx[1] -= 1.0f;
          }
          if (dx[1] < -0.5f) {
            dx[1] += 1.0f;
          }
#ifdef VORONOI_VERBOSE
          message("Cell %i before:", i);
          voronoi_print_cell(&cells[i]);
          message("Interacting cell %i with cell %i (%g %g, %g %g", i, j,
                  cells[i].x[0], cells[i].x[1], cells[j].x[0], cells[j].x[1]);
#endif
          voronoi_cell_interact(cell_i, dx, j);
        }
      }
    }

    Atot = 0.0f;
    /* print the cells to the stdout */
    for (i = 0; i < 100; ++i) {
#ifdef VORONOI_VERBOSE
      printf("Cell %i:\n", i);
      voronoi_print_cell(&cells[i]);
#endif
      voronoi_cell_finalize(&cells[i]);
      Atot += cells[i].volume;
    }

    /* Check the total surface area */
    assert(fabs(Atot - 1.0f) < 1.e-6);

    /* Check the neighbour relations for an arbitrary cell: cell 44
       We plotted the grid and manually found the correct neighbours and their
       order. */
    assert(cells[44].nvert == 4);
    assert(cells[44].ngbs[0] == 34);
    assert(cells[44].ngbs[1] == 45);
    assert(cells[44].ngbs[2] == 54);
    assert(cells[44].ngbs[3] == 43);

    message("Done.");
  }

  return 0;
}
