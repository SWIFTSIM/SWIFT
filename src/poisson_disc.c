/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/**
 *  @file poisson_disc.c
 *  @brief Generate 3D random points using a poisson disc distribution.
 *
 *  Poisson disc requires that each position is at least a given
 *  radius from all other points, so is a good way to sample a grid
 *  without introducing aliasing.
 *
 *  See: Robert Bridson. 2007. Fast Poisson disk sampling in
 *       arbitrary dimensions. In ACM SIGGRAPH 2007 sketches
 *       (SIGGRAPH '07). ACM, New York, NY, USA, , Article 22.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <values.h>

/* Local headers. */
#include "error.h"

/* Useful defines. */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define CHUNK 512

/**
 * @brief 3d position.
 */
struct point {
  float x[3];
};

/**
 * @brief collection of data for the grid.
 */
struct grid_data {
  struct point *cells;
  float cell_size;
  float radius;
  int dims[3];
  int ncells;

  struct point *activelist;
  int actlen;
  int actsize;

  struct point *samplelist;
  int samplelen;
  int samplesize;
};

/**
 * @brief Get random number in range 0 to 1.
 */
static float randone() { return (float)(rand() / (double)RAND_MAX); }

/**
 *  @brief Add a position to the active list.
 */
static void add_active(struct grid_data *grid, struct point *p) {

  for (int i = 0; i < 3; i++)
    grid->activelist[grid->actlen].x[i] = p->x[i];
  grid->actlen++;

  /*  Add more space to activelist, if needed. */
  if (grid->actlen >= grid->actsize) {
    grid->actsize = grid->actlen + CHUNK;
    grid->activelist =
        realloc(grid->activelist, grid->actsize * sizeof(struct point));
    if ( grid->activelist == NULL ) {
        error( "Failed to realloc active list" );
    }
  }
}

/**
 *  @brief Remove a position from the active list.
 */
static void remove_active(struct grid_data *grid, int index) {

  /*  Shuffle activelist removing index. */
  for (int i = index; i < grid->actlen - 1; i++)
    for (int j = 0; j < 3; j++)
      grid->activelist[i].x[j] = grid->activelist[i + 1].x[j];

  grid->actlen--;
}

/**
 * @brief Mark a position on the grid and add to the active list.
 */
static void mark_position(struct grid_data *grid, struct point *p) {

  /* Index of point on grid. */
  int index = (int)grid->dims[1] * grid->dims[0] * floor(p->x[2] / grid->cell_size) +
              grid->dims[0] * floor(p->x[1] / grid->cell_size) +
              floor(p->x[0] / grid->cell_size);

  /* Check if already seen, nothing to do. */
  if (grid->cells[index].x[0] == -1.0f) {
    for (int i = 0; i < 3; i++)
      grid->cells[index].x[i] = p->x[i];

    add_active(grid, p);
  }
}

/**
 *  @brief Test if a position is available on the grid.
 *
 *  That is further than the radius away from positions already marked on the
 *  grid and not marked.
 */
static int not_marked(struct grid_data *grid, struct point *p) {

  int i = (int)p->x[1] / grid->cell_size;
  int j = (int)p->x[0] / grid->cell_size;
  int l = (int)p->x[2] / grid->cell_size;
  int i0 = MAX(i - 2, 0);
  int j0 = MAX(j - 2, 0);
  int l0 = MAX(l - 2, 0);
  int j1 = MIN(j + 3, grid->dims[0]);
  int i1 = MIN(i + 3, grid->dims[1]);
  int l1 = MIN(l + 3, grid->dims[2]);

  for (int j = j0; j < j1; j++) {
    for (int i = i0; i < i1; i++) {
      for (int l = l0; l < l1; l++) {
        int index = l * grid->dims[1] * grid->dims[0] + i * grid->dims[0] + j;
        if (grid->cells[index].x[0] != -1.0) {
          double dsq = 0.0;
          for (int j = 0; j < 3; j++) {
            float dx = grid->cells[index].x[j] - p->x[j];
            dsq += dx * dx;
          }
          if (dsq < (grid->radius * grid->radius)) {
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

/**
 * @brief Add a position to final sample.
 */
static void add_sample(struct grid_data *grid, struct point *p) {

  for (int i = 0; i < 3; i++)
    grid->samplelist[grid->samplelen].x[i] = p->x[i];
  grid->samplelen++;

  /*  Add more space to samples, if needed. */
  if (grid->samplelen >= grid->samplesize) {
    grid->samplesize = grid->samplesize + CHUNK;
    grid->samplelist =
        realloc(grid->samplelist, grid->samplesize * sizeof(struct point));
    if ( grid->samplelist == NULL ) {
      error("Failed to realloc sample list");
    }
  }
}

/**
 * @brief free up any resources allocated by the grid.
 */
void poisson_free(struct grid_data *grid) {

  if (grid->cells != NULL) free(grid->cells);
  if (grid->activelist != NULL) free(grid->activelist);
  if (grid->samplelist != NULL) free(grid->samplelist);
}

/**
 * @brief Generate a sample of positions.
 *
 * On return the grid_data struct is populated with a samplelist containing
 * the selected positions.
 *
 * @param grid the struct to contain all the data about the sample selection.
 * @param dims the dimensions of the space to sample.
 * @param radius minimum radius between selected positions.
 * @param k number of attempts to pick an unmarked position (30 is a good
 * choice).
 */
void poisson_disc(struct grid_data *grid, int dims[3], float radius, int k) {

  grid->radius = radius;
  grid->cell_size = radius / sqrtf(3.0f);
  grid->ncells = 1;
  for (int i = 0; i < 3; i++) {
    grid->dims[i] = (int)ceilf(dims[i] / grid->cell_size);
    grid->ncells *= grid->dims[i];
  }

  grid->cells = (struct point *)malloc(sizeof(struct point) * grid->ncells );
  if (grid->cells == NULL) {
    error("Failed to allocate grid cells");
  }

  for (int i = 0; i < grid->ncells; i++) {
    grid->cells[i].x[0] = -1.0;
  }

  /*  Queue for active list. */
  grid->activelist = (struct point *)malloc(CHUNK * sizeof(struct point));
  if (grid->activelist == NULL) {
    error("Failed to allocate an active list");
  }
  grid->actsize = CHUNK;
  grid->actlen = 0;

  /*  Space for results. */
  grid->samplelist = (struct point *)malloc(CHUNK * sizeof(struct point));
  if ( grid->samplelist == NULL ) {
    error("Failed to allocate a sample list");
  }
  grid->samplesize = CHUNK;
  grid->samplelen = 0;

  /* Initialise a seed point. */
  struct point ip;
  for (int i = 0; i < 3; i++)
      ip.x[i] = randone() * dims[i];
  mark_position(grid, &ip);

  while (grid->actlen > 0) {

    /*  Grab a position from the activelist. */
    int index = (int)(randone() * grid->actlen);
    struct point *p = &grid->activelist[index];
    int havenew = 0;

    /*  Look for a random position that is far enough away and not already
     *  in use. */
    for (int j = 0; j < k; j++) {
      float theta = 2.0f * M_PI * randone();
      float phi = 2.0f * M_PI * randone();

      /*  Radius to this position is outside the expected radius. */
      float r = randone() * radius + radius;
      struct point np;
      np.x[0] = p->x[0] + r * sin(theta) * cos(phi);
      np.x[1] = p->x[1] + r * sin(theta) * sin(phi);
      np.x[2] = p->x[2] + r * cos(theta);

      if (np.x[0] >= 0.0f && np.x[0] < dims[0] &&
          np.x[1] >= 0.0f && np.x[1] < dims[1] &&
          np.x[2] >= 0.0f && np.x[2] < dims[2] &&
          not_marked(grid, &np)) {

        /* Accepted, so mark and add to the active list. */
        mark_position(grid, &np);
        havenew = 1;
        break;
      }
    }

    /*  Remove point from active list and keep as no other position is
     *  near. */
    if (!havenew) {
      add_sample(grid, p);
      remove_active(grid, index);
    }
  }
}

int main(int argc, char *argv[]) {
  int dims[3];
  float radius = 0.0f;
  int k = 30;
  struct grid_data grid;

  /* Expected sample size. */
  int N = 90;

  /* Dimensions. */
  dims[0] = 300;
  dims[1] = 300;
  dims[2] = 300;

  srand(time(NULL));

  /*  Pick radius for expected sample size. */
  radius = pow((dims[0] * dims[1] * dims[2]) / (N), 0.3333);
  printf("# Radius = %f\n", radius);

  /*  Sample is stocastic, so we may need to ask more than one to get the
   *  number of samples we require as a minimum. */
  grid.samplelen = 0;
  while (grid.samplelen < N) {
    printf("# Sampling...\n");
    poisson_disc(&grid, dims, radius, k);
    printf("# Samples = %d\n", grid.samplelen);
  }

  for (int i = 0; i < grid.samplelen; i++) {
    printf("# %f %f %f\n", grid.samplelist[i].x[0], grid.samplelist[i].x[1],
           grid.samplelist[i].x[2]);
  }

  /*  Partition the space. Slow .... */

  unsigned long int counts[N];
  for (int i = 0; i < N; i++) {
    counts[i] = 0;
  }

  for (int i = 0; i < dims[0]; i++) {
    for (int j = 0; j < dims[1]; j++) {
      for (int k = 0; k < dims[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        for (int l = 0; l < N; l++) {
          float dx = grid.samplelist[l].x[0] - (i + 0.5);
          float dy = grid.samplelist[l].x[1] - (j + 0.5);
          float dz = grid.samplelist[l].x[2] - (k + 0.5);
          float rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        counts[select]++;
        printf("%f %f %f %d\n", i + 0.5, j + 0.5, k + 0.5, select);
      }
    }
  }

  printf("# Counts:\n");
  unsigned long int total = 0;
  for (int i = 0; i < N; i++) {
    printf("#  %d %ld\n", i, counts[i]);
    total += counts[i];
  }
  printf("# total = %ld\n", total);

  poisson_free( &grid );

  return 0;
}
