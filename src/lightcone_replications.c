/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "align.h"
#include "cell.h"
#include "error.h"
#include "memuse.h"
#include "lightcone_replications.h"


/**
 * @brief Comparison function for sorting replications
 *
 * a The first replication
 * b The second replication
 *
 */
static int compare_replication_rmin(const void *a, const void *b) {
  const struct replication *rep_a = (struct replication *) a;
  const struct replication *rep_b = (struct replication *) b;
  if(rep_a->rmin2 < rep_b->rmin2)
    return -1;
  else if(rep_a->rmin2 > rep_b->rmin2)
    return 1;
  else
    return 0;
}


/**
 * @brief Make a list of periodic box replications which overlap 
 *        the specified distance range from an observer.
 *
 * @param boxsize Size of the cubic simulation box.
 * @param observer_position Location of the observer.
 * @param lightcone_rmin Minimum distance from the observer.
 * @param lightcone_rmax Maximum distance from the observer.
 *        of particles.
 * @param replication_list Pointer to the struct to initialise.
 */

void replication_list_init(struct replication_list *replication_list,
                           double boxsize, double observer_position[3],
                           double lightcone_rmin, double lightcone_rmax) {
  
  /* Find range of replications to examine in each dimension */
  int rep_min[3];
  int rep_max[3];
  for(int i=0; i<3; i+=1) {
    rep_min[i] = (int) floor((observer_position[i] - lightcone_rmax) / boxsize);
    rep_max[i] = (int) floor((observer_position[i] + lightcone_rmax) / boxsize);
  }

  /* On first pass just count replications */
  for(int ipass=0; ipass<2; ipass+=1) {

    replication_list->nrep = 0;

    /* Loop over periodic replications */
    for(int i=rep_min[0]; i<=rep_max[0]; i+=1) {
      for(int j=rep_min[1]; j<=rep_max[1]; j+=1) {
        for(int k=rep_min[2]; k<=rep_max[2]; k+=1) {
          
          /* Find centre of this replication relative to observer */
          double cx = boxsize*i + 0.5*boxsize - observer_position[0];
          double cy = boxsize*j + 0.5*boxsize - observer_position[1];
          double cz = boxsize*k + 0.5*boxsize - observer_position[2];

          /* Find distance to closest point in this replication  */
          double dx, dy, dz;
          dx = abs(cx) - 0.5*boxsize;
          if(dx < 0) dx = 0;
          dy = abs(cy) - 0.5*boxsize;
          if(dy < 0) dy = 0;
          dz = abs(cz) - 0.5*boxsize;
          if(dz < 0) dz = 0;
          double rep_rmin = sqrt(dx*dx+dy*dy+dz*dz);

          /* Find distance to most distant point in this replication  */
          dx = abs(cx) + 0.5*boxsize;
          dy = abs(cy) + 0.5*boxsize;
          dz = abs(cz) + 0.5*boxsize;
          double rep_rmax = sqrt(dx*dx+dy*dy+dz*dz);

          /* Flag if any point in this replication could be in the lightcone */
          int in_lightcone = 1;

          /* Check distance limits */
          if(rep_rmax < lightcone_rmin || rep_rmin > lightcone_rmax) in_lightcone = 0;

          if(in_lightcone) {
            /* Store replications on second pass */
            if(ipass==1) {
              /* Get a pointer to the next replication */
              const int nrep = replication_list->nrep;
              struct replication *rep = replication_list->replication+nrep;
              /* Store info about this replication */
              rep->rmin2 = pow(rep_rmin, 2.0);
              rep->rmax2 = pow(rep_rmax, 2.0);
              rep->coord[0] = i*boxsize;
              rep->coord[1] = j*boxsize;
              rep->coord[2] = k*boxsize;
            }
            replication_list->nrep += 1; 
          }
        } /* Next replication in z */
      } /* Next replication in y */
    } /* Next replication in x */

    /* Allocate storage after first pass */
    if(ipass==0) {
      const int nrep = replication_list->nrep;
      if(swift_memalign("lightcone_replications", (void **) &replication_list->replication,
                        SWIFT_STRUCT_ALIGNMENT, sizeof(struct replication)*nrep) != 0) {
        error("Failed to allocate lightcone replication list");
      }
    }
  } /* Next pass */

  /* Now sort replications by minimum distance */
  qsort(replication_list->replication, 
        (size_t) replication_list->nrep, 
        sizeof(struct replication),
        compare_replication_rmin);
}


/**
 * @brief Deallocate a replication list
 *
 * @param replication_list Pointer to the struct to deallocate.
 */
void replication_list_clean(struct replication_list *replication_list) {
  swift_free("lightcone_replications", replication_list->replication);
  replication_list->replication = NULL;
  replication_list->nrep = 0;
}

/**
 * @brief Write a replication list to a file as text
 *
 * @param replication_list The replication list
 * @param fd The file to write to
 */

void replication_list_write(struct replication_list *replication_list, FILE *fd) {

  for(int i=0; i<replication_list->nrep; i+=1) {
    fprintf(fd, "%e, %e, %e, %e, %e\n",
            replication_list->replication[i].coord[0],
            replication_list->replication[i].coord[1],
            replication_list->replication[i].coord[2],
            sqrt(replication_list->replication[i].rmin2),
            sqrt(replication_list->replication[i].rmax2));
  }
}


/**
 * Determine subset of replications which overlap a #cell
 *
 * @param rep_in The input replication list
 * @param cell The input cell
 * @param rep_out The output replication list
 *
 * Initializes rep_out, which must then be freed with
 * replication_list_clean().
 *
 */
void replication_list_subset_for_cell(const struct replication_list *rep_in,
                                      const struct cell *cell,
                                      const double observer_position[3],
                                      struct replication_list *rep_out) {
  
  /* Find centre coordinates of this cell */
  const double cell_centre[] = {cell->loc[0]+0.5*cell->width[0],
                                cell->loc[1]+0.5*cell->width[1],
                                cell->loc[2]+0.5*cell->width[2]};
  
  /* Find 'effective' width of this cell - particles can wander out of the cell */
  const double cell_eff_width[] = {2.0*cell->width[0],
                                   2.0*cell->width[1],
                                   2.0*cell->width[2]};

  /* Allocate array of replications for the new list */
  const int nrep_max = rep_in->nrep;
  if(swift_memalign("lightcone_replications", (void **) &rep_out->replication,
                    SWIFT_STRUCT_ALIGNMENT, sizeof(struct replication)*nrep_max) != 0) {
    error("Failed to allocate pruned lightcone replication list");
  }

  /* Loop over all replications */
  rep_out->nrep = 0;
  for(int i=0; i<nrep_max; i+=1) {
    
    /* Get a pointer to this input replication */
    const struct replication *rep = rep_in->replication+i;

    /* Find coordinates of centre of this replication of the cell relative to the observer */
    double cell_rep_centre[3];
    for(int j=0; j<3; j+=1) {
      cell_rep_centre[j] = rep->coord[j] + cell_centre[j] - observer_position[j];
    }

    /* Compute minimum possible distance squared from observer to this replication of this cell */
    double cell_rmin2 = 0.0;
    for(int j=0; j<3; j+=1) {
      double dx = abs(cell_rep_centre[j]) - 0.5*cell_eff_width[j];
      if(dx < 0.0)dx = 0.0;
      cell_rmin2 += dx*dx;
    }

    /* Compute maximum possible distance squared from observer to this replication of this cell */
    double cell_rmax2 = 0.0;
    for(int j=0; j<3; j+=1) {
      double dx = abs(cell_rep_centre[j]) + 0.5*cell_eff_width[j];
      cell_rmax2 += dx*dx;
    }

    /* Decide whether this cell could contribute to this replication */
    if(cell_rmax2 >= rep->rmin2 && cell_rmin2 <= rep->rmax2) {
      memcpy(rep_out->replication+rep_out->nrep, rep, sizeof(struct replication));
      rep_out->nrep +=1;
    }
    /* Next input replication */
  }
}
