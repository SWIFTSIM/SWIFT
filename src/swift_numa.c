/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <config.h>

#include <numa.h>
#include <numaif.h>
#include <stdio.h>
#include <errno.h>

#include "error.h"
#include "swift_numa.h"


// XXX need to get this from the machine...
#define NR_NUMA_NODES 8


/**
 *  @file numa.c
 *  @brief file of various techniques for using NUMA awareness.
 */


/**
 * @param binds a memory region to a NUMA node and move them if necessary.
 */
static void memory_to_numa_node(const char *label, void *mem, size_t size,
                                int node, int pol) {
  struct bitmask *nodes;
  nodes = numa_allocate_nodemask();

  /* Move memory if we have already touched it. */
  int mbind_flags = MPOL_MF_MOVE;
  
  if (pol == MPOL_DEFAULT) {
    /* Need the empty set for this policy. */
    nodes->maskp = NULL;
    nodes->size = 0;
  } else {
    numa_bitmask_setbit(nodes, node);
  }
  int res = mbind(mem, size, pol, nodes->maskp, nodes->size, mbind_flags);
  if (res > 0) {
    message("%s: error in mbind: %d", label, res);
  }

  numa_bitmask_free(nodes);
}

//  Spread an allocation over a number of NUMA regions.
void swift_numa_spread(const char *label, void *ptr, size_t size, int unbind) {

  // Small allocations stay where they are.
  if (size < 4096 * NR_NUMA_NODES * 64) {
    //message("small allocation %zd not spread", size);
    return;
  }

  // Spread this over the available nodes.
  //size_t chunk_size = size / NR_NUMA_NODES;

  size_t chunk_size = size / (NR_NUMA_NODES * 64);

  //message("%s: chunk_size = %zd size = %zd @ %p", label, chunk_size, size, ptr);

  // chunk_size has to be a multiple of the page size.
  chunk_size = (chunk_size / 4096 + 1) * 4096;
  //message("%s: page chunk_size = %zd size = %zd", label, chunk_size, size);

  char *lptr = (char *)ptr;
  size_t eaten = 0;

  /* The policy, either bind to the nodes, or set back to default, an
   * effective unbind, I hope. */
  int policy = MPOL_PREFERRED;
  if (unbind) policy = MPOL_DEFAULT;

  // Start at our current NUMA node to avoid adding all memory to 0 first,
  // important for small allocations, which should remain local.
  int ind = numa_preferred();

  //for (int k = 0; k < NR_NUMA_NODES; k++) {
  while (chunk_size > 0) {
    //message("%s: %d %zd %p", label, ind, chunk_size, lptr);
    memory_to_numa_node(label, lptr, chunk_size, ind, policy);
    lptr += chunk_size;
    eaten += chunk_size;

    //  Ragged limit.
    if (eaten + chunk_size > size) chunk_size = size - eaten;

    //  Roll over when more nodes to consider.
    ind++;
    if (ind >= NR_NUMA_NODES) ind = 0;
  }
  //fflush(stdout);
}

// Bind memory to a NUMA node.
void swift_numa_bind(const char *label, void *mem, size_t size, int node) {
  memory_to_numa_node(label, mem, size, node, MPOL_PREFERRED);
}

//  Get the NUMA node that holds the given memory address.
//  Use with care this will touch a page at that address if not already done,
//  so will be created with the current memory policy. Returns -1
//  if the address isn't valid, usually some alias or NULL.
int swift_get_numa_node(void *ptr) {
  int numa_node = -1;
  get_mempolicy(&numa_node, NULL, 0, (void*)ptr, MPOL_F_NODE | MPOL_F_ADDR);
  return numa_node;
}

//  Bind the current thread and all following memory allocations to a NUMA
//  node. Only when node is valid. If both is false only bind the thread
//  not the memory. Passing a -1 will clear any previous binding of the
//  thread to the CPU, but that may not unbind the memory, not clear about that.
void swift_set_numa_node(int node, int both) {

  // Clear any existing affinity. XXX is this necessary, it shouldn't be.
  numa_run_on_node_mask(numa_all_nodes_ptr);
  numa_run_on_node(node);
  if (both) {
    struct bitmask *nodemask = numa_allocate_nodemask();
    if (node != -1) numa_bitmask_setbit(nodemask, node);
    numa_set_membind(nodemask);
    numa_free_nodemask(nodemask);
  }
}
