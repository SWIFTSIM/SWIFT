/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "../config.h"

#include <float.h>

#include "cell.h"
#include "gravity_properties.h"
#include "engine.h"
#include "proxy.h"
#include "partition.h"
#include "space.h"
#include "zoom_region.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/*  Vectorisation support */
/*  ===================== */

#if defined(WITH_MPI)
#ifdef WITH_ZOOM_REGION
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
void pick_vector_zoom(struct space *s, int nregions, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  length += s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }
	/* Set up variables */
  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;

  /* Loop over zoom grid */
	for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
	  for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
	    for (int k = 0; k < s->zoom_props->cdim[2]; k++) {
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

  /* Loop over natural grid */
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
#endif
#endif

#if defined(WITH_MPI)
#ifdef WITH_ZOOM_REGION
/**
 * @brief Partition the space including the zoom cells.
 *
 * Using the sample positions as seeds pick cells that are geometrically
 * closest and apply the partition to the space.
 */
void split_vector_zoom(struct space *s, int nregions, int *samplecells) {

	/* Define variables for selection */
  int cid = 0;

  /* Loop over zoom cells*/
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {
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
        s->cells_top[cid++].nodeID = select;
      }
    }
  }

  /* Loop over natural cells*/
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
        s->cells_top[cid++].nodeID = select;
      }
    }
  }
}
#endif
#endif

