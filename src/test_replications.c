#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "periodic_replications.h"

int main(int argc, char *argv[]) {

  /* Input parameters */
  double boxsize = strtod(argv[1], NULL);
  double observer_position[3];
  for(int i=0; i<3; i+=1)
    observer_position[i] = strtod(argv[2+i], NULL);
  double lightcone_rmin = strtod(argv[5], NULL);
  double lightcone_rmax = strtod(argv[6], NULL);

  /* Make the list */
  struct replication_list replication_list;
  replication_list_init(&replication_list,
                        boxsize, observer_position,
                        lightcone_rmin, lightcone_rmax);
  
  /* Output the list */
  //printf("nrep=%d\n", replication_list.nrep);
  for(int i=0; i<replication_list.nrep; i+=1)
    {
      printf("%d, %d, %d, %e, %e\n",
             replication_list.replication[i].coord[0],
             replication_list.replication[i].coord[1],
             replication_list.replication[i].coord[2],
             sqrt(replication_list.replication[i].rmin2),
             sqrt(replication_list.replication[i].rmax2));
    }

  /* Free the list */
  replication_list_clean(&replication_list);

  return 0;
}
