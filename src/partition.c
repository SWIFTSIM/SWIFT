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
 *  @file partition.c
 *  @brief file of various techniques for partitioning a grid of cells
 *         into geometrically connected regions.
 *
 *  Currently supported types, grid, random, vectorise and metis.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS headers only used when MPI is also available. */
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "partition.h"
#include "const.h"
#include "error.h"
#include "debug.h"

/* Useful defines. */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define CHUNK 512

/* Simple descriptions of initial partition types for reports. */
const char *initpart_name[] = {
  "gridded cells",
  "random point associated cells",
  "vectorized point associated cells",
  "METIS particle weighted cells",
  "METIS unweighted cells"
};

/* Simple descriptions of repartition types for reports. */
const char *repart_name[] = {
  "no",
  "METIS edge and vertex weighted cells",
  "METIS vertex weighted cells",
  "METIS edge weights"
};

/*  Vectorisation support */
/*  ===================== */

/* Vectorise test routine
int main(int argc, char *argv[]) {

  int N = 0;
  int D1 = 0;
  int D2 = 0;
  int D3 = 0;
  if ( argc > 4 ) {
    D1 = atoi( argv[1] );
    D2 = atoi( argv[2] );
    D3 = atoi( argv[3] );
    N = atoi( argv[4] );

    struct space s;
    int samplecells[3*N];
    unsigned int counts[N];

    s.cdim[0] = D1;
    s.cdim[1] = D2;
    s.cdim[2] = D3;

    part_pick_vector(&s, N, samplecells);
    part_apply_vector(&s, N, samplecells, 0, counts);

  } else {
    message( "no parts supplied" );
  }
}
*/


/**
 *  @brief Pick a number of cell positions from a vectorised list.
 *
 *  Vectorise the cell space and pick positions in it for the number of
 *  expected regions using a single step. Vectorisation is guaranteed
 *  to work, providing there are more cells than regions.
 *
 *  @param s the space.
 *  @param nregions the number of regions
 *  @param samplecells the list of sample cell positions, size of 3*nregions
 */
void part_pick_vector(struct space *s, int nregions, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }

  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;

  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        if (n == 0 && l < nregions) {
          samplecells[m++] = i;
          samplecells[m++] = j;
          samplecells[m++] = k;
          l++;
        }
        n++;
        if (n == step) n = 0;
      }
    }
  }
}

/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometry closest
 * to each and apply the partition to the space.
 */
void part_split_vector(struct space *s, int nregions, int *samplecells) {


  /* XXX While testing keep the counts per component. */
  int counts[nregions];
  for (int i = 0; i < nregions; i++) {
    counts[i] = 0;
  }

  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells[n++].nodeID = select;
        //message("@ %d %d %d %d", i, j, k, select);
        counts[select]++;
      }
    }
  }

  /* XXX Test section */
  message("# counts:");
  unsigned int total = 0;
  for (int i = 0; i < nregions; i++) {
    message("#  %d %d", i, counts[i]);
    if ( counts[i] == 0 ) {
      message( "sampling failed" );
    }
    total += counts[i];
  }
  message("# total = %d", total);
}

/* Random points support.
 * ======================
 *
 *  Generates 3D random points using a poisson disc distribution.
 *
 *  Poisson disc requires that each position is at least a given radius from
 *  all other points, so is a good way to sample a grid without introducing
 *  aliasing and keeps a strong geometric relationship between cells. The
 *  downside is that it is not possible to completely guarantee a solution as
 *  the required radius is difficult to estimate and the match between
 *  spherical regions and the cell fixed grid leaves gaps that can result in a
 *  position not being close to a cell. Should only happen in extreme cases.
 *
 *  See: Robert Bridson. 2007. Fast Poisson disk sampling in
 *       arbitrary dimensions. In ACM SIGGRAPH 2007 sketches
 *       (SIGGRAPH '07). ACM, New York, NY, USA, , Article 22.
 */

/* Test routine for random sampling.
int main( int argc, char *argv[] )
{
    int N = 5;
    if ( argc > 1 ) {
        N = atoi( argv[1] );
    }

    float samplelist[N*4];
    struct space s;
    s.cdim[0] = 6;
    s.cdim[1] = 6;
    s.cdim[2] = 6;

    if ( part_pick_random(&s, N, samplelist))
        part_split_random(&s, N, samplelist);

    return 0;
}
*/

/**
 * @brief 3d position.
 */
struct random_point {
  float x[3];
};

/**
 * @brief collection of data for the various elements used.
 */
struct part_random_data {
  struct random_point *cells;
  float cell_size;
  float radius;
  int dims[3];
  int ncells;

  struct random_point *activelist;
  int actlen;
  int actsize;

  struct random_point *samplelist;
  int samplelen;
  int samplesize;
};

/**
 * @brief Get random number in range 0 to 1.
 */
static float random_one() {
  return (float)(rand() / (double)RAND_MAX);
}

/**
 *  @brief Add a position to the active list.
 */
static void random_add_active(struct part_random_data *data,
                              struct random_point *p) {

  for (int i = 0; i < 3; i++)
    data->activelist[data->actlen].x[i] = p->x[i];
  data->actlen++;

  /*  Add more space to activelist, if needed. */
  if (data->actlen >= data->actsize) {
    data->actsize = data->actlen + CHUNK;
    data->activelist =
        realloc(data->activelist, data->actsize * sizeof(struct random_point));
    if ( data->activelist == NULL ) {
        error( "Failed to realloc active list" );
    }
  }
}

/**
 *  @brief Remove a position from the active list.
 */
static void random_remove_active(struct part_random_data *data, int index) {

  /*  Shuffle activelist removing index. */
  for (int i = index; i < data->actlen - 1; i++)
    for (int j = 0; j < 3; j++)
      data->activelist[i].x[j] = data->activelist[i + 1].x[j];

  data->actlen--;
}

/**
 * @brief Mark a position on the grid and add to the active list.
 */
static void random_mark_position(struct part_random_data *data,
                                 struct random_point *p) {

  /* Index of point on grid. */
    int index =
        data->dims[1] * data->dims[0] * (int)(p->x[2] / data->cell_size) +
        data->dims[0] * (int)(p->x[1] / data->cell_size) +
        (int)(p->x[0] / data->cell_size);

  /* Check if already seen, nothing to do. */
  if (data->cells[index].x[0] == -1.0f) {
    for (int i = 0; i < 3; i++)
      data->cells[index].x[i] = p->x[i];

    random_add_active(data, p);
  }
}

/**
 *  @brief Test if a position is available on the grid.
 *
 *  That is further than the radius away from positions already marked on the
 *  grid and not marked.
 */
static int random_not_marked(struct part_random_data *data,
                             struct random_point *p) {

  int i = (int)(p->x[1] / data->cell_size);
  int j = (int)(p->x[0] / data->cell_size);
  int l = (int)(p->x[2] / data->cell_size);
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
static void random_add_sample(struct part_random_data *data,
                              struct random_point *p) {

  for (int i = 0; i < 3; i++)
    data->samplelist[data->samplelen].x[i] = p->x[i];
  data->samplelen++;

  /*  Add more space to samples, if needed. */
  if (data->samplelen >= data->samplesize) {
    data->samplesize = data->samplesize + CHUNK;
    data->samplelist = realloc(data->samplelist,
                               data->samplesize * sizeof(struct random_point));
    if ( data->samplelist == NULL ) {
      error("Failed to realloc sample list");
    }
  }
}

/**
 * @brief Partition the space returning the counts per cells and optionally
 * applying the partition to the space.
 */
static void random_partition(struct space *s, int nregions, float *samplelist,
                             int apply, unsigned long int *counts) {

  for (int i = 0; i < nregions; i++) {
    counts[i] = 0;
  }

  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplelist[m++] - (i + 0.5);
          float dy = samplelist[m++] - (j + 0.5);
          float dz = samplelist[m++] - (k + 0.5);
          float w = samplelist[m++];
          float rsq = w * (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        if ( apply )
          s->cells[n++].nodeID = select;
        counts[select]++;
      }
    }
  }
}

/**
 * @brief Generate a sample of positions.
 *
 * On return the part_random_data struct is populated with a samplelist
 * containing the selected positions.
 *
 * @param data the struct to contain all the data about the sample selection.
 * @param dims the dimensions of the space to sample.
 * @param radius minimum radius between selected positions.
 * @param k number of attempts to pick an unmarked position (30 is a good
 * choice).
 */
static void random_disc(struct part_random_data *data, int dims[3],
                        float radius, int k) {

  /* Initialise a seed point. */
  struct random_point ip;
  for (int i = 0; i < 3; i++)
      ip.x[i] = random_one() * dims[i];
  random_mark_position(data, &ip);
  message( "# initial position: %f, %f, %f", ip.x[0], ip.x[1], ip.x[2]);

  while (data->actlen > 0) {

    /*  Grab a position from the activelist. */
    int index = (int)(random_one() * data->actlen);
    struct random_point *p = &data->activelist[index];
    int havenew = 0;

    /*  Look for a random position that is far enough away and not already
     *  in use. */
    for (int j = 0; j < k; j++) {
      float theta = 2.0f * M_PI * random_one();
      float phi = 2.0f * M_PI * random_one();

      /*  Radius to this position is outside the expected radius. */
      float r = random_one() * radius + 2.0f * radius;
      struct random_point np;
      np.x[0] = p->x[0] + r * sin(theta) * cos(phi);
      np.x[1] = p->x[1] + r * sin(theta) * sin(phi);
      np.x[2] = p->x[2] + r * cos(theta);

      if (np.x[0] >= 0.0f && np.x[0] < dims[0] &&
          np.x[1] >= 0.0f && np.x[1] < dims[1] &&
          np.x[2] >= 0.0f && np.x[2] < dims[2] &&
          random_not_marked(data, &np)) {

        /* Accepted, so mark and add to the active list. */
        random_mark_position(data, &np);
        havenew = 1;
        break;
      }
    }

    /*  Remove point from active list and keep as no other position is
     *  near. */
    if (!havenew) {
      random_add_sample(data, p);
      random_remove_active(data, index);
    }
  }
}

/**
 * @brief Generate a poisson disc sample for a 3D space.
 *
 * Generates a poisson disc sample for the 3D space of the cells returning
 * a structure that can be applied to the cells of a space to partition
 * the space into regions.
 *
 * @param s the space to partitions.
 * @param nregions number of regions to generate.
 * @param samplelist memory for a list of nregions*3 float coordinates that are
 *                   the sample points, plus one weight.
 * @return 1 for success, 0 for failure.
 */
int part_pick_random(struct space *s, int nregions, float *samplelist) {

  int result = 1;
  struct part_random_data data;

  /* Randomize randoms. */
  srand(time(NULL));

  /* Pick radius for expected sample size. This part is tricky, we want a
   * radius that will sample the space well, but the random selection and,
   * potentially, non-cube nature of the space stops a space filling solution
   * from ever working, so we guess and need to possibly seek more than one
   * solution. The 4.18 is space filling, but we allow centres to be closer to
   * the edges than radius, giving the slop and allowing a dimension to be
   * smaller than 2*radius.
   */
  float radius = pow((s->cdim[0]*s->cdim[1]*s->cdim[2])/(4.18*nregions), 0.3333);
  data.radius = radius;
  message("# radius = %f", radius);

  /* Cell size for this resolution. */
  data.cell_size = radius / sqrtf(3.0f);

  /* Count cells in the sample space. */
  data.ncells = 1;
  for (int i = 0; i < 3; i++) {
    data.dims[i] = (int)ceilf(s->cdim[i] / data.cell_size);
    data.ncells *= data.dims[i];
  }
  message("# partition space has %d cells (%dx%dx%d, size:%f)", data.ncells,
          data.dims[0], data.dims[1], data.dims[2], data.cell_size);

  /* If we want more samples than cells, then things are not going to work. */
  if (nregions > data.ncells) {
    message("Cannot sub-sample space");
    return 0;
  }
  data.cells =
      (struct random_point *)malloc(sizeof(struct random_point) * data.ncells);
  if (data.cells == NULL) {
    error("Failed to allocate data cells");
  }

  /*  Queue for active list. */
  data.activelist =
      (struct random_point *)malloc(CHUNK * sizeof(struct random_point));
  if (data.activelist == NULL) {
    error("Failed to allocate an active list");
  }
  data.actsize = CHUNK;

  /*  Space for results. */
  data.samplelist =
      (struct random_point *)malloc(CHUNK * sizeof(struct random_point));
  if (data.samplelist == NULL) {
    error("Failed to allocate a sample list");
  }
  data.samplesize = CHUNK;

  /* Get the sample, only try a number of times. */
  int ntries = 0;
  data.samplelen = 0;
  while (data.samplelen < nregions && ntries < 30) {

    /* Re-initialize data struct making sure it is clean of previous
     * attempts. */
    for (int i = 0; i < data.ncells; i++) {
      data.cells[i].x[0] = -1.0;
    }
    data.samplelen = 0;
    data.actlen = 0;

    random_disc(&data, s->cdim, radius, 30);
    message("# samples = %d", data.samplelen);
    ntries++;

    /* If samplelen is greater than 2 * nregions, shrink the radius. */
    if (data.samplelen > 2 * nregions && ntries < 10) {
        radius = radius * 0.75f;
        data.radius = radius;
        message("# radius becomes = %f", radius);

        /* Cell size for this resolution. */
        data.cell_size = radius / sqrtf(3.0f);

        /* Count cells in the sample space. */
        data.ncells = 1;
        for (int i = 0; i < 3; i++) {
            data.dims[i] = (int)ceilf(s->cdim[i] / data.cell_size);
            data.ncells *= data.dims[i];
        }
        message("# partition space has %d cells (%dx%dx%d, size:%f)",
                data.ncells, data.dims[0], data.dims[1], data.dims[2],
                data.cell_size);
        data.samplelen = 0;
    }
  }

  if (data.samplelen < nregions) {
    message("failed to partition space (last sample size: %d, expected: %d)",
            data.samplelen, nregions);
    result = 0;
  } else {

    /* Extract the samplelist. Weights initially set to 1.*/
    for (int i = 0, k = 0; i < nregions; i++) {
      for (int j = 0; j < 3; j++)
        samplelist[k++] = data.samplelist[i].x[j];
      samplelist[k++] = 1.0f;
    }

    for (int i = 0; i < data.samplelen; i++) {
      message("# %d: %f %f %f", i,
              data.samplelist[i].x[0],
              data.samplelist[i].x[1],
              data.samplelist[i].x[2]);
    }

    /* Check the radii are as expected. */
    for (int i = 0; i < data.samplelen; i++) {
      float ref[3];
      ref[0] = data.samplelist[i].x[0];
      ref[1] = data.samplelist[i].x[1];
      ref[2] = data.samplelist[i].x[2];
      for (int j = i + 1; j < data.samplelen; j++) {
        float dx = ref[0] - data.samplelist[j].x[0];
        float dy = ref[1] - data.samplelist[j].x[1];
        float dz = ref[2] - data.samplelist[j].x[2];
        if (sqrtf(dx*dx+dy*dy+dz*dz) < radius) {
          message( "point too close: %d,%d (%f)", i, j,
                   sqrtf(dx*dx+dy*dy+dz*dz));
        }
      }
    }

    /* Attempt to share cells around more evenly by adding weight factors
     * based on the cells numbers per component. First need to get the
     * partition. */
    unsigned long int *counts = malloc( sizeof(unsigned long int) * nregions);
    if (counts == NULL) {
      error("Failed to allocate counts array");
    }
    random_partition(s, nregions, samplelist, 0, counts);

    /* Use counts to populate the weights. */
    double total = 0.0;
    double cmin = DBL_MAX;
    double cmax = 0.0;
    for (int i = 0; i < nregions; i++) {
        message("#  %d %ld", i, counts[i]);
        if (counts[i] == 0) {
          message( "sampling failed" );
        }
        total += counts[i];
        if ( counts[i] < cmin ) cmin = counts[i];
        if ( counts[i] > cmax ) cmax = counts[i];
    }

    for (int i = 0, k = -1; i < nregions; i++) {
        k = k + 4;
        if (counts[i] != 0) {
            samplelist[k] = (float)1.0 + (counts[i]-cmin)/(cmax-cmin);

        } else {
            samplelist[k] = 0.75f;
        }
        message("weight %d: %f", i, samplelist[k]);
    }
    free(counts);
  }

  if (data.cells != NULL) free(data.cells);
  if (data.activelist != NULL) free(data.activelist);
  if (data.samplelist != NULL) free(data.samplelist);

  return result;
}

/**
 * @brief Apply poisson disc partitioning to a cell structure.
 *
 * Uses the results of random_generate to assign each cell's nodeID to the
 * nearest sample point index, thus partitioning the space into regions.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param samplelist pointer to sample coordinates and weights (from
 *                   random_generate).
 */
void part_split_random(struct space *s, int nregions, float *samplelist) {

  for (int l = 0, m = 0; l < nregions; l++, m += 4)
    message( "# %f %f %f (%f)", samplelist[m], samplelist[m+1],
             samplelist[m+2], samplelist[m+3]);

  /* Partition the space. */
  unsigned long int counts[nregions];
  random_partition(s, nregions, samplelist, 1, counts);

  /* Test section */
  message("# counts:");
  unsigned long int total = 0;
  for (int i = 0; i < nregions; i++) {
    message("#  %d %ld", i, counts[i]);
    if ( counts[i] == 0 ) {
        message( "sampling failed" );
    }
    total += counts[i];
  }
  message("# total = %ld", total);
}

/* METIS support
 * =============
 *
 * METIS partitions using a multi-level k-way scheme. We support using this in
 * a unweighted scheme, which works well and seems to be guaranteed, and a
 * weighted by the number of particles scheme. Note METIS is optional.
 */
/* Test routine for METIS
int main(int argc, char *argv[]) {

  int N = 0;
  int D1 = 0;
  int D2 = 0;
  int D3 = 0;
  if ( argc > 4 ) {
    D1 = atoi( argv[1] );
    D2 = atoi( argv[2] );
    D3 = atoi( argv[3] );
    N = atoi( argv[4] );

    struct space s;

    s.cdim[0] = D1;
    s.cdim[1] = D2;
    s.cdim[2] = D3;

    message("# Partition space of %d,%d,%d cells into %d", s.cdim[0], s.cdim[1],
            s.cdim[2], N);

    part_metis(&s, N);

  } else {
    message( "no parts supplied" );
  }
}
*/

/* Convert cell location to sequence number. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the
 * cell particle counts as weights.
 *
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param vertexw weights for the cells, sizeof number of cells if used, 
 *        NULL for unit weights.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjcny array.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 */
void part_pick_metis(struct space *s, int nregions, int *vertexw,
                     int *edgew, int *celllist) {

#if defined(HAVE_METIS)

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if ( nregions == 1 ) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* Allocate weights and adjacency arrays . */
  idx_t *xadj;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
    error("Failed to allocate xadj buffer.");
  idx_t *adjncy;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
    error("Failed to allocate adjncy array.");
  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate vertex weights array");
  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)malloc( 26 *sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate edge weights array");
  idx_t *regionid;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
    error("Failed to allocate regionid array");

  /* Fill the xadj and adjncy array to define the graph of cells. */
  /* Loop over all cells. */
  int cid = 0;
  for (int l = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {

        /* Visit all neighbours of this cell, wrapping space at edges. */
        int p = 0;
        for (int i = -1; i <= 1; i++) {
          int ii = l + i;
          if (ii < 0)
            ii += s->cdim[0];
          else if (ii >= s->cdim[0])
            ii -= s->cdim[0];
          for (int j = -1; j <= 1; j++) {
            int jj = m + j;
            if (jj < 0)
              jj += s->cdim[1];
            else if (jj >= s->cdim[1])
              jj -= s->cdim[1];
            for (int k = -1; k <= 1; k++) {
              int kk = n + k;
              if (kk < 0)
                kk += s->cdim[2];
              else if (kk >= s->cdim[2])
                kk -= s->cdim[2];

              /* If not self, record id of neighbour. */
              if (i || j || k) {
                adjncy[cid * 26 + p] = cell_getid(s->cdim, ii, jj, kk);
                p++;
              }
            }
          }
        }

        /* Next cell. */
        cid++;
      }
    }
  }
  xadj[0] = 0;
  for (int k = 0; k < ncells; k++) xadj[k + 1] = xadj[k] + 26;

  /* Init the vertex weights array. */
  if (vertexw != NULL) {
    for (int k = 0; k < ncells; k++) {
      if ( vertexw[k] > 0 ) {
        weights_v[k] = vertexw[k];
      } else {
        weights_v[k] = 1;
      }
    }
  }

  /* Init the edges weights array. */
  if (edgew != NULL) {
    for (int k = 0; k < 26 * ncells; k++) {
      if (edgew[k] > 0 ) {
        weights_e[k] = edgew[k];
      } else {
        weights_e[k] = 1;
      }
    }
  }

  /* Set the METIS options. */
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_NCUTS] = 10;
  options[METIS_OPTION_NITER] = 20;

  /* Call METIS. */
  idx_t one = 1;
  idx_t idx_ncells = ncells;
  idx_t idx_nregions = nregions;
  idx_t objval;

  /* Dump graph in METIS format */
  dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy,
                 weights_v, weights_e, NULL);

  if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v,
                          weights_e, NULL, &idx_nregions, NULL, NULL,
                          options, &objval, regionid) != METIS_OK)
    error("Call to METIS_PartGraphKway failed.");

  /* Check that the regionid are ok. */
  for (int k = 0; k < ncells; k++)
    if (regionid[k] < 0 || regionid[k] >= nregions)
      error("Got bad nodeID %"PRIDX" for cell %i.", regionid[k], k);


  /* Check that the partition is complete and all regions have some
   * cells. */
  int *present = malloc(sizeof(int) * nregions);
  if (present == NULL)
    error("Failed to allocate present buffer.");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < ncells; i++) present[regionid[i]]++;
  for (int i = 0; i < nregions; i++) {
    message("count %d %d", i, present[i]);
    if (! present[i]) {
      failed = 1;
      message("Region %d is not present in partition", i);
    }
  }
  if ( failed ) {
    message("Partitioning failed");
  }

  /* Set the cell list to the region index. */
  for (int k = 0; k < ncells; k++) {
    celllist[k] = regionid[k];
  }

  for (int l = 0, k = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {
        message("node %d %d %d -> %d", l, m, n, regionid[k]);
        k++;
      }
    }
  }

  /* Clean up. */
  free(adjncy);
  free(weights_v);
  free(regionid);
  free(present);

#else
  error("SWIFT was not compiled with METIS support.");
#endif
}

/**
 * @brief Apply METIS cell list partitioning to a cell structure.
 *
 * Uses the results of part_metis_pick to assign each cell's nodeID to the
 * picked region index, thus partitioning the space into regions.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
void part_split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++)
    s->cells[i].nodeID = celllist[i];
}

