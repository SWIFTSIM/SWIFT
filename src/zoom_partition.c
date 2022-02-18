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

/**
 * @brief Check if all regions have been assigned a node in the
 *        cells of a space.
 *
 * @param s the space containing the cells to check.
 * @param nregions number of regions expected.
 * @param verbose if true report the missing regions.
 * @return true if all regions have been found, false otherwise.
 */
static int check_complete(struct space *s, int verbose, int nregions) {

	int *present = NULL;
	if ((present = (int *)malloc(sizeof(int) * nregions)) == NULL)
		error("Failed to allocate present array");

	int failed = 0;
	for (int i = 0; i < nregions; i++) present[i] = 0;
	for (int i = 0; i < s->nr_cells; i++) {
		if (s->cells_top[i].nodeID <= nregions)
			present[s->cells_top[i].nodeID]++;
		else
			message("Bad nodeID: s->cells_top[%d].nodeID = %d", i,
			        s->cells_top[i].nodeID);
	}
	for (int i = 0; i < nregions; i++) {
		if (!present[i]) {
			failed = 1;
			if (verbose) message("Region %d is not present in partition", i);
		}
	}
	free(present);
	return (!failed);
}

/**
 * @brief Partition a space of cells based on another space of cells.
 *
 * The two spaces are expected to be at different cell sizes, so what we'd
 * like to do is assign the second space to geometrically closest nodes
 * of the first, with the effect of minimizing particle movement when
 * rebuilding the second space from the first.
 *
 * Since two spaces cannot exist simultaneously the old space is actually
 * required in a decomposed state. These are the old cells sizes and counts
 * per dimension, along with a list of the old nodeIDs. The old nodeIDs are
 * indexed by the cellid (see cell_getid()), so should be stored that way.
 *
 * On exit the new space cells will have their nodeIDs assigned.
 *
 * @param oldh the cell dimensions of old space.
 * @param oldcdim number of cells per dimension in old space.
 * @param oldnodeIDs the nodeIDs of cells in the old space, indexed by old
 *cellid.
 * @param s the space to be partitioned.
 *
 * @return 1 if the new space contains nodeIDs from all nodes, 0 otherwise.
 */
int partition_space_to_space_zoom(double *oldh, double *oldcdim, double *oldzoomh,
		                              double *oldzoomcdim, int *oldnodeIDs, struct space *s) {

	/* Define the old tl_cell_offset */
	const int old_bkg_cell_offset = oldzoomcdim[0] * oldzoomcdim[1] * oldzoomcdim[2];

	/* Loop over all the new zoom cells. */
	for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
		for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
			for (int k = 0; k < s->zoom_props->cdim[2]; k++) {

				/* Scale indices to old cell space. */
				const int ii = rint(i * s->zoom_props->iwidth[0] * oldzoomh[0]);
				const int jj = rint(j * s->zoom_props->iwidth[1] * oldzoomh[1]);
				const int kk = rint(k * s->zoom_props->iwidth[2] * oldzoomh[2]);

				const int cid = cell_getid(s->zoom_props->cdim, i, j, k);
				const int oldcid = cell_getid(oldzoomcdim, ii, jj, kk);
				s->cells_top[cid].nodeID = oldnodeIDs[oldcid];
			}
		}
	}

	/* Loop over all the new cells. */
	for (int i = 0; i < s->cdim[0]; i++) {
		for (int j = 0; j < s->cdim[1]; j++) {
			for (int k = 0; k < s->cdim[2]; k++) {

				/* Scale indices to old cell space. */
				const int ii = rint(i * s->iwidth[0] * oldh[0]);
				const int jj = rint(j * s->iwidth[1] * oldh[1]);
				const int kk = rint(k * s->iwidth[2] * oldh[2]);

				const int cid = cell_getid(s->cdim, i, j, k) + s->zoom_props->tl_cell_offset;
				const int oldcid = cell_getid(oldcdim, ii, jj, kk) + old_bkg_cell_offset;
				s->cells_top[cid].nodeID = oldnodeIDs[oldcid];
			}
		}
	}

	/* Check we have all nodeIDs present in the resample. */
	return check_complete(s, 1, s->e->nr_nodes);
}

/*  Vectorisation support */
/*  ===================== */

#ifdef WITH_MPI
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
  int length = (s->cdim[0] * s->cdim[1] * s->cdim[2]) + (s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2]);
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

#ifdef WITH_MPI
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

