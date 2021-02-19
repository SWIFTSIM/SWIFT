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

#include "periodic_replications.h"


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

/* TODO: ensure lightcone_rmin is reduced to allow for drifting */

void replication_list_init(struct replication_list *replication_list,
                           double boxsize,
                           double observer_position[3],
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
          
          /* Find centre of this replication */
          double cx = boxsize*i + 0.5*boxsize;
          double cy = boxsize*j + 0.5*boxsize;
          double cz = boxsize*k + 0.5*boxsize;

          /* Find distance to closest point in this replication  */
          double dx, dy, dz;
          dx = abs(observer_position[0] - cx) - 0.5*boxsize;
          if(dx < 0) dx = 0;
          dy = abs(observer_position[1] - cy) - 0.5*boxsize;
          if(dy < 0) dy = 0;
          dz = abs(observer_position[2] - cz) - 0.5*boxsize;
          if(dz < 0) dz = 0;
          double rep_rmin = sqrt(dx*dx+dy*dy+dz*dz);

          /* Find distance to most distant point in this replication  */
          dx = abs(observer_position[0] - cx) + 0.5*boxsize;
          dy = abs(observer_position[1] - cy) + 0.5*boxsize;
          dz = abs(observer_position[2] - cz) + 0.5*boxsize;
          double rep_rmax = sqrt(dx*dx+dy*dy+dz*dz);

          /* Check if any point in this replication could be in the lightcone */
          if(rep_rmax > lightcone_rmin && rep_rmin < lightcone_rmax) {
    
            /* Store replications on second pass */
            if(ipass==1) {
              /* Get a pointer to the next replication */
              const int nrep = replication_list->nrep;
              struct replication *rep = replication_list->replication+nrep;
              /* Store info about this replication */
              rep->rmin2 = pow(rep_rmin, 2.0);
              rep->rmax2 = pow(rep_rmax, 2.0);
              rep->coord[0] = i;
              rep->coord[1] = j;
              rep->coord[2] = k;
            }
            replication_list->nrep += 1;
            
          }
        } /* Next replication in z */
      } /* Next replication in y */
    } /* Next replication in x */

    /* Allocate storage after first pass */
    if(ipass==0) {
      const int nrep = replication_list->nrep;
      replication_list->replication = malloc(sizeof(struct replication)*nrep); 
    }
  } /* Next pass */

  /* Now sort replications by minimum distance */
  qsort(replication_list->replication, 
        (size_t) replication_list->replication, 
        sizeof(struct replication),
        compare_replication_rmin);
}


/**
 * @brief Deallocate a replication list
 *
 * @param replication_list Pointer to the struct to deallocate.
 */
void replication_list_clean(struct replication_list *replication_list) {
  free(replication_list->replication);
  replication_list->replication = NULL;
  replication_list->nrep = 0;
}
