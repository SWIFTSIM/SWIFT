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
#include "poisson_disc.h"

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
 * @brief collection of data for the various elements in constructing the grid.
 */
struct poisson_data {
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
static void add_active(struct poisson_data *data, struct point *p) {

  for (int i = 0; i < 3; i++)
    data->activelist[data->actlen].x[i] = p->x[i];
  data->actlen++;

  /*  Add more space to activelist, if needed. */
  if (data->actlen >= data->actsize) {
    data->actsize = data->actlen + CHUNK;
    data->activelist =
        realloc(data->activelist, data->actsize * sizeof(struct point));
    if ( data->activelist == NULL ) {
        error( "Failed to realloc active list" );
    }
  }
}

/**
 *  @brief Remove a position from the active list.
 */
static void remove_active(struct poisson_data *data, int index) {

  /*  Shuffle activelist removing index. */
  for (int i = index; i < data->actlen - 1; i++)
    for (int j = 0; j < 3; j++)
      data->activelist[i].x[j] = data->activelist[i + 1].x[j];

  data->actlen--;
}

/**
 * @brief Mark a position on the grid and add to the active list.
 */
static void mark_position(struct poisson_data *data, struct point *p) {

  /* Index of point on grid. */
  int index = (int)data->dims[1] * data->dims[0] * floor(p->x[2] / data->cell_size) +
              data->dims[0] * floor(p->x[1] / data->cell_size) +
              floor(p->x[0] / data->cell_size);

  /* Check if already seen, nothing to do. */
  if (data->cells[index].x[0] == -1.0f) {
    for (int i = 0; i < 3; i++)
      data->cells[index].x[i] = p->x[i];

    add_active(data, p);
  }
}

/**
 *  @brief Test if a position is available on the grid.
 *
 *  That is further than the radius away from positions already marked on the
 *  grid and not marked.
 */
static int not_marked(struct poisson_data *data, struct point *p) {

  int i = (int)p->x[1] / data->cell_size;
  int j = (int)p->x[0] / data->cell_size;
  int l = (int)p->x[2] / data->cell_size;
  int i0 = MAX(i - 2, 0);
  int j0 = MAX(j - 2, 0);
  int l0 = MAX(l - 2, 0);
  int j1 = MIN(j + 3, data->dims[0]);
  int i1 = MIN(i + 3, data->dims[1]);
  int l1 = MIN(l + 3, data->dims[2]);

  for (int j = j0; j < j1; j++) {
    for (int i = i0; i < i1; i++) {
      for (int l = l0; l < l1; l++) {
        int index = l * data->dims[1] * data->dims[0] + i * data->dims[0] + j;
        if (data->cells[index].x[0] != -1.0) {
          double dsq = 0.0;
          for (int j = 0; j < 3; j++) {
            float dx = data->cells[index].x[j] - p->x[j];
            dsq += dx * dx;
          }
          if (dsq < (data->radius * data->radius)) {
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
static void add_sample(struct poisson_data *data, struct point *p) {

  for (int i = 0; i < 3; i++)
    data->samplelist[data->samplelen].x[i] = p->x[i];
  data->samplelen++;

  /*  Add more space to samples, if needed. */
  if (data->samplelen >= data->samplesize) {
    data->samplesize = data->samplesize + CHUNK;
    data->samplelist =
        realloc(data->samplelist, data->samplesize * sizeof(struct point));
    if ( data->samplelist == NULL ) {
      error("Failed to realloc sample list");
    }
  }
}

/**
 * @brief Generate a sample of positions.
 *
 * On return the poisson_data struct is populated with a samplelist containing
 * the selected positions.
 *
 * @param data the struct to contain all the data about the sample selection.
 * @param dims the dimensions of the space to sample.
 * @param radius minimum radius between selected positions.
 * @param k number of attempts to pick an unmarked position (30 is a good
 * choice).
 */
static void poisson_disc(struct poisson_data *data, int dims[3], float radius,
                         int k) {

  /* Initialise a seed point. */
  struct point ip;
  for (int i = 0; i < 3; i++)
      ip.x[i] = randone() * dims[i];
  mark_position(data, &ip);

  while (data->actlen > 0) {

    /*  Grab a position from the activelist. */
    int index = (int)(randone() * data->actlen);
    struct point *p = &data->activelist[index];
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
          not_marked(data, &np)) {

        /* Accepted, so mark and add to the active list. */
        mark_position(data, &np);
        havenew = 1;
        break;
      }
    }

    /*  Remove point from active list and keep as no other position is
     *  near. */
    if (!havenew) {
      add_sample(data, p);
      remove_active(data, index);
    }
  }
}

/**
 * @brief Generate a poisson disc sample for a 3D space.
 *
 * Generates a poisson disc sample for the 3D space of the cells returning
 * a structure that can be applied to the cells of a space to partition
 * the space into contiguous regions.
 *
 * @param s the space to partitions.
 * @param nparts number of parititions to generate.
 * @param samplelist memory for a list of nparts*3 float coordinates that are
 *                   the sample points. 
 * @return 1 for success, 0 for failure.
 */
int poisson_generate(struct space *s, int nparts, float *samplelist) {

  int result = 1;

  /* Randomize randoms. */
  srand(time(NULL));

  /* Pick radius for expected sample size. This part is tricky, we want a
   * radius that will sample the space well, but the random selection stops a
   * space filling solution from ever working, so we guess and need to
   * possibly seek more than one solution.
   */
  float radius = pow((s->cdim[0]*s->cdim[1]*s->cdim[2])/(1.5*nparts), 0.3333);
  message("# Radius = %f", radius);

  /* Initialise the data struct. */
  struct poisson_data data;
  data.radius = radius;

  /* Cell size for this resolution. */
  data.cell_size = radius / sqrtf(3.0f);
  data.ncells = 1;
  for (int i = 0; i < 3; i++) {
    data.dims[i] = (int)ceilf(s->cdim[i] / data.cell_size);
    data.ncells *= data.dims[i];
  }

  data.cells = (struct point *)malloc(sizeof(struct point) * data.ncells );
  if (data.cells == NULL) {
    error("Failed to allocate data cells");
  }
  for (int i = 0; i < data.ncells; i++) {
    data.cells[i].x[0] = -1.0;
  }

  /*  Queue for active list. */
  data.activelist = (struct point *)malloc(CHUNK * sizeof(struct point));
  if (data.activelist == NULL) {
    error("Failed to allocate an active list");
  }
  data.actsize = CHUNK;
  data.actlen = 0;

  /*  Space for results. */
  data.samplelist = (struct point *)malloc(CHUNK * sizeof(struct point));
  if (data.samplelist == NULL) {
    error("Failed to allocate a sample list");
  }
  data.samplesize = CHUNK;
  data.samplelen = 0;

  /* Get the sample, try a number of times, but give up so that we never get a
   * loop if our guess ever fails completely. */
  int ntries = 0;
  while (data.samplelen < nparts && ntries < 10) {
    poisson_disc(&data, s->cdim, radius, 30);
    message("# Samples = %d", data.samplelen);
    ntries++;
  }
  if ( data.samplelen < nparts ) {
    message("Failed to partition space (last sample size: %d, expected: %d)",
            data.samplelen, nparts);
    result = 0;
  } else {

     /* Extract the samplelist. */
     for (int i = 0, k = 0; i < nparts; i++)
       for (int j = 0; j < 3; j++)
         samplelist[k++] = data.samplelist[i].x[j];

     for ( int i = 0; i < data.samplelen; i++) {
       message("# %f %f %f",
               data.samplelist[i].x[0],
               data.samplelist[i].x[1],
               data.samplelist[i].x[2]);
     }
  }

  if (data.cells != NULL) free(data.cells);
  if (data.activelist != NULL) free(data.activelist);
  if (data.samplelist != NULL) free(data.samplelist);

  return result;
}

/**
 * @brief Apply poisson disc partitioning to a cell structure.
 *
 * Uses the results of poisson_generate to assign each cell's nodeID to the
 * nearest sample point index, thus partitioning the space into contiguous
 * regions.
 *
 * @param s the space containing the cells to split into partitions.
 * @param nparts number of parititions.
 * @param samplelist pointer to sample coordinates (from poisson_generate).
 */
void poisson_split(struct space *s, int nparts, float *samplelist) {

  for (int l = 0, m = 0; l < nparts; l++, m += 3)
    message( "# %f %f %f", samplelist[m], samplelist[m+1], samplelist[m+2]);

  /*  Partition the space. Slow ??? */
  unsigned long int counts[nparts];
  for (int i = 0; i < nparts; i++) {
    counts[i] = 0;
  }

  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nparts; l++) {
          float dx = samplelist[m++] - (i + 0.5);
          float dy = samplelist[m++] - (j + 0.5);
          float dz = samplelist[m++] - (k + 0.5);
          float rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells[n++].nodeID = select;
        counts[select]++;
        //message("%f %f %f %d", i + 0.5, j + 0.5, k + 0.5, select);
      }
    }
  }

  message("# Counts:");
  unsigned long int total = 0;
  for (int i = 0; i < nparts; i++) {
    message("#  %d %ld", i, counts[i]);
    total += counts[i];
  }
  message("# total = %ld", total);
}

/*
int main( int argc, char *argv[] )
{

    float *samplelist = NULL;
    struct space s;
    s.cdim[0] = 100;
    s.cdim[1] = 12;
    s.cdim[2] = 12;

    samplelist = poisson_generate(&s, 2);

    poisson_split(&s, 2, samplelist);

}
*/
