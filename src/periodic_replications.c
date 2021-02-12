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
 * @brief Make a list of periodic box replications which overlap 
 *        the specified distance range from an observer.
 *
 * @param boxsize Size of the cubic simulation box.
 * @param observer_position Location of the observer.
 * @param lightcone_rmin Minimum distance from the observer.
 * @param lightcone_rmax Maximum distance from the observer.
 * @param lightcone_boundary Size of buffer zone to allow for movement
 *        of particles.
 * @param replication_list Pointer to the struct to initialise.
 */
void replication_list_init(double boxsize,
                           double observer_position[3],
                           double lightcone_rmin, double lightcone_rmax,
                           double lightcone_boundary,
                           struct replication_list *replication_list) {
  
  /* Find range of replications to examine in each dimension */
  int rep_min[3];
  int rep_max[3];
  for(i=0; i<3; i+=1) {
    rep_min[i] = (int) floor((observer_position[i] - lightcone_rmax) / boxsize);
    rep_min[i] = (int) floor((observer_position[i] + lightcone_rmax) / boxsize);
  }

  /* On first pass just count replications */
  for(ipass=0; ipass<2; ipass+=1) {

    replication_list->nrep = 0;

    /* Loop over periodic replications */
    for(i=rep_min[0]; i<=rep_max[0]; i+=1) {
      for(j=rep_min[1]; j<=rep_max[1]; j+=1) {
        for(k=rep_min[2]; k<=rep_max[2]; k+=1) {
          
          /* Find centre of this replication */
          double cx = boxsize*i + 0.5*boxsize;
          double cy = boxsize*j + 0.5*boxsize;
          double cz = boxsize*k + 0.5*boxsize;

          /* Find distance to closest point in this replication  */
          double dx = abs(observer_position[0] - cx) - 0.5*boxsize;
          if(dx < 0) dx = 0;
          double dy = abs(observer_position[1] - cy) - 0.5*boxsize;
          if(dy < 0) dy = 0;
          double dz = abs(observer_position[2] - cz) - 0.5*boxsize;
          if(dz < 0) dz = 0;
          double rep_rmin = sqrt(dx*dx+dy*dy+dz*dz);

          /* Find distance to most distant point in this replication  */
          double dx = abs(observer_position[0] - cx) + 0.5*boxsize;
          double dy = abs(observer_position[1] - cy) + 0.5*boxsize;
          double dz = abs(observer_position[2] - cz) + 0.5*boxsize;
          double rep_rmax = sqrt(dx*dx+dy*dy+dz*dz);

          /* Check if any point in this replication could be in the lightcone,
             allowing a boundary layer to account for the distance particles 
             can be drifted */
          if(rep_rmax > lightcone_rmin-lightcone_boundary &&
             rep_rmin < lightcone_rmax+lightcone_boundary) {
    
            /* Store replications on second pass */
            if(ipass==1) {
              /* Get a pointer to the next replication */
              const int nrep = replication_list->nrep;
              struct replication *rep = replication_list->replication+nrep;
              /* Store info about this replication */
              rep.rmin = rep_rmin;
              rep.rmax = rep_rmax;
              rep.coord[0] = i;
              rep.coord[1] = j;
              rep.coord[2] = k;
            }
            replication_list->nrep += 1;
            
          }
        } /* Next replication in z */
      } /* Next replication in y */
    } /* Next replication in x */

    /* Allocate storage after first pass */
    if(ipass==0) {
      const int nrep = replication_list->nrep;
      replication_list->replication = malloc(sizeof(struct replication)*(*nrep)); 
    }
  } /* Next pass */
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


/* int main(int argc, char *argv[]) { */

/*   double boxsize = 100.0; */
/*   double observer_position = {50.0, 50.0, 50.0}; */
/*   double lightcone_rmin = 250.0; */
/*   double lightcone_rmax = 500.0; */
/*   double lightcone_boundary = 0.0; */
/*   int nrep; */
/*   struct replication *replication; */

/*   make_replication_list(boxsize, observer_position, */
/*                         lightcone_rmin, lightcone_rmax, */
/*                            double lightcone_boundary, */
/*                            int *nrep, struct replication *replication) */
  
  
/*   return 0; */
/* } */
