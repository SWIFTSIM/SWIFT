/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 STFC (Author email aidan.chalk@stfc.ac.uk)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <stdint.h>

#ifdef HAVE_LIBNUMA
#include <unistd.h>
#include <numa.h>
#include <numaif.h>
#endif

/* Include own header. */
#include "numa_helper.h"

#include "engine.h"
#include "space.h"
#include "cell.h"



/**
 * @brief Move a cell structure's hydro parts into a single NUMA region.
 *
 * This can be done in parallel. Any pages not completely inside this cell will only be moved if
 * the over half of the page is in this cell.
 * This function is THREAD SAFE.
 * @param c The #cell.
 * @param node The numa node to move to.
 * @param verbose Are we talkative/do we check things.
 */
int swiftnuma_cell_move_hydro_parts(struct cell *c, int32_t node, int32_t verbose){
#ifdef HAVE_LIBNUMA
  struct part *start_part = c->hydro.parts;
  struct part *end_part = &c->hydro.parts[c->hydro.count];
  const uintptr_t page_size = getpagesize();
  uint8_t *start_page = (uint8_t*)(((uintptr_t)start_part) & ~(page_size-1));
  uint8_t *last_page = (uint8_t*)(((uintptr_t)end_part) & ~(page_size-1));
  /* We only move a page if over half of it belongs to this cell. */
  if( ((uintptr_t)(((uint8_t*) start_part) - start_page)) <= (page_size/2) ){
    start_page += page_size;
  }
  if( ((uintptr_t)(((uint8_t*)(end_part+1)) - last_page)) <= (page_size/2) ){
    end_page -= page_size;
  }
  uintptr_t nr_pages = ((uintptr_t)(last_page - start_page)) / page_size;
  if(nr_pages > 0){
    void** pages = malloc(sizeof(void*) * nr_pages);
    int* domain = malloc(sizeof(int) * nr_pages); // NUMA call requires int inputs
    int* status = malloc(sizeof(int) * nr_pages);
    if(pages == NULL || domain == NULL || status == NULL){
      error("Failed to allocate arrays to move pages");
    }

#ifdef SWIFT_DEBUG_CHECKS
  //TODO check we are allowed/it is sane to put things on this NUMA domain.
#ifdef HAVE_SETAFFINITY
  cpu_set_t *global_affinity = engine_entry_affinity();
  const int nr_affinity_cores = CPU_COUNT(entry_affinity);
  int* cpuid = (int *)malloc(nr_affinity_cores * sizeof(int));

  int skip = 0;
  for (int k = 0; k < nr_affinity_cores; k++) {
    int c;
    for (c = skip; c < CPU_SETSIZE && !CPU_ISSET(c, entry_affinity); ++c)
      ;
    cpuid[k] = c;
    skip = c + 1;
  }
  int found = 0;
  for(int i = 0; i < nr_affinity_cores; i++) {
    int l = numa_node_of_cpu(cpuid[i]);
    if(l == node){
      found = 1;
      break;
    }
  }
  if(!found){
    error("We are attempting to allocate a cell's hydro particles to a NUMA node that is not in the affinity mask of the program");
  }
#endif
#endif
    for(uintptr_t z = 0; z < nr_pages; z++){
      pages[z] = (void*) (start_page + z*page_size);
      domain[z] = node;
    }
    long returnval = numa_move_pages(0,(unsigned long)nr_pages, pages, domain, status, MPOL_MF_MOVE);
    if( verbose && returnval < 0 ){
      error("numa_move_pages returned an error code %li", returnval);
    }
    if( verbose ){
        for(int32_t k = 0; k < nr_pages; k++){
          if(status[k] != domain[k]){
             message("Not all pages moved successfully");
             break;
          }
        }
    }
    free(pages);
    free(domain);
    free(status);
  }
#else
  error("SWIFT not compiled with NUMA support")
#endif
}
