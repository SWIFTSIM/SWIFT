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

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
#include <numa.h>
#include <numaif.h>
#endif
#include <stdio.h>
#include <errno.h>
#include <unistd.h>

#include "error.h"
#include "swift_numa.h"

/* If these functions are active. If they are all noops. */
static int ACTIVE;

/* Number of NUMA nodes. */
static int NR_NUMA_NODES;

/* Page size of machine. */
static int PAGESIZE;

/**
 *  @file numa.c
 *  @brief file of various techniques for using NUMA awareness.
 */

/**
 * @brief Must be called before any other routines in this file and before any
 *        swift_ memory allocations.
 * @param active whether functions should be active, otherwise they do nothing.
 */
void swift_numa_init(int active) {

  ACTIVE = active;
  if (!ACTIVE) return;

  /* Must be true on any machine. */
  NR_NUMA_NODES = 1;
  PAGESIZE = sysconf(_SC_PAGESIZE);

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  if (numa_available() >= 0) NR_NUMA_NODES = numa_num_task_nodes();
#endif

}

/**
 * @brief mbind a memory region to a NUMA node and moving if necessary.
 *
 * @param label descriptive context for memory region (see swift_malloc etc).
 * @param ptr pointer to the memory.
 * @param size size in bytes of the memory region.
 * @param node the NUMA node to bind memory to.
 * @param pol the NUMA memory policy to use. Use MPOL_DEFAULT to clear an
 *            existing policy.
 */
static void memory_to_numa_node(const char *label, void *ptr, size_t size,
                                int node, int pol) {
#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
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
  int res = mbind(ptr, size, pol, nodes->maskp, nodes->size, mbind_flags);
  if (res > 0) {
    message("%s: error in mbind: %d", label, res);
  }
  numa_bitmask_free(nodes);
#endif
}

/**
 * @brief Spread an memory allocation over all the NUMA nodes.
 *
 * @param label descriptive context for memory region (see swift_malloc etc).
 * @param ptr pointer to the memory.
 * @param size size in bytes of the memory region.
 * @param unbind if 1 undo the previous memory bind.
 */
void swift_numa_spread(const char *label, void *ptr, size_t size, int unbind) {

  if (!ACTIVE) return;

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)

  /* Smaller allocations stay where they are. These should probably be
   * tunable. */
  if (size < (size_t)PAGESIZE * NR_NUMA_NODES * 128) return;

  /* Spread this over the available nodes in some larger chunking than the
   * default, which is the pagesize. */
  size_t chunk_size = size / (NR_NUMA_NODES * 128);

  //message("%s: chunk_size = %zd size = %zd @ %p", label, chunk_size, size, ptr);

  /* chunk_size has to be a multiple of the page size. */
  chunk_size = (chunk_size / PAGESIZE + 1) * PAGESIZE;
  //message("%s: page chunk_size = %zd size = %zd", label, chunk_size, size);

  char *lptr = (char *)ptr;
  size_t eaten = 0;

  /* The policy, either bind to the nodes, or set back to default, an
   * effective unbind, I hope. */
  int policy = MPOL_PREFERRED;
  if (unbind) policy = MPOL_DEFAULT;

  /* Start at our current NUMA node to avoid adding all memory to 0 first,
   * important for small allocations, which should start locally. */
  int ind = numa_preferred();

  while (chunk_size > 0) {
    //message("%s: %d %zd %p", label, ind, chunk_size, lptr);
    memory_to_numa_node(label, lptr, chunk_size, ind, policy);
    lptr += chunk_size;
    eaten += chunk_size;

    /*  Ragged limit. */
    if (eaten + chunk_size > size) chunk_size = size - eaten;

    /*  Roll over when more nodes to consider. */
    ind++;
    if (ind >= NR_NUMA_NODES) ind = 0;
  }
#endif
}

/**
 * @brief Bind memory to a NUMA node.
 *
 * @param label descriptive context for memory region (see swift_malloc etc).
 * @param ptr pointer to the memory.
 * @param size size in bytes of the memory region.
 * @param node the NUMA node to bind memory to.
 */
void swift_numa_bind(const char *label, void *ptr, size_t size, int node) {

  if (!ACTIVE) return;

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  memory_to_numa_node(label, ptr, size, node, MPOL_PREFERRED);
#endif
}

/**
 * @brief Get the NUMA node that holds the given memory address.
 *
 * Use with care this will touch a page at that address if not already done,
 * so will be created with the current memory policy.
 *
 * @param ptr the memory address to check.
 *
 * @result -1 if the address isn't valid, usually some alias or NULL.
 */
int swift_get_numa_node(void *ptr) {

  if (!ACTIVE) return -1;

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  int numa_node = -1;
  get_mempolicy(&numa_node, NULL, 0, (void*)ptr, MPOL_F_NODE | MPOL_F_ADDR);
  return numa_node;
#else
  return -1;
#endif
}

/**
 * @brief Bind the current thread and, optionally, all following memory
 * allocations to a NUMA node.
 *
 * Only works when the node is valid (i.e. <=0). If both is false only binds
 * the thread not the memory. Passing a -1 will clear any previous binding of
 * the thread to the CPU, but that may not unbind the memory, not clear about
 * that.
 *
 * @param node the NUMA node to bind the currrent thread to.
 * @param both if true also arrange to bind local memory allocations.
 */
void swift_set_numa_node(int node, int both) {

  if (!ACTIVE) return;

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  /* Clear any existing affinity. */
  numa_run_on_node_mask(numa_all_nodes_ptr);

  /* And bind the thread. */
  numa_run_on_node(node);

  /* Also memory allocations. */
  if (both) {
    struct bitmask *nodemask = numa_allocate_nodemask();
    if (node != -1) numa_bitmask_setbit(nodemask, node);
    numa_set_membind(nodemask);
    numa_free_nodemask(nodemask);
  }
#endif
}
