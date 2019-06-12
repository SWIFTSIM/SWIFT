/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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
 *  @file memuse.c
 *  @brief file of routines to report about memory use in SWIFT.
 *  Note reports are in KB.
 */

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "memuse.h"

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "engine.h"

#ifdef SWIFT_MEMUSE_REPORTS

/* The initial size and increment of the log entries buffer. */
#define MEMUSE_INITLOG 1000000

/* A megabyte for conversions. */
#define MEGABYTE 1048576.0

/* Maximum length of label in log entry. */
#define MEMUSE_MAXLABLEN 32

/* Also recorded in logger. */
extern int engine_rank;
extern int engine_current_step;

/* Entry for logger of memory allocations and deallocations in a step. */
struct memuse_log_entry {

  /* Step of action. */
  int step;

  /* Whether allocated or deallocated. */
  int allocated;

  /* Memory allocated in bytes. */
  size_t size;

  /* Address of memory. Use union as easy way to convert into an array of
   * bytes. */
  union {
    void *ptr;
    uint8_t vptr[sizeof(uintptr_t)];
  };

  /* Relative time of this action. */
  ticks dtic;

  /* Label associated with the memory. */
  char label[MEMUSE_MAXLABLEN + 1];

  /* Whether log is still active, i.e. not matched with a free or allocation. */
  int active;
};

/* The log of allocations and frees. All volatile as accessed from threads
 * that use the value to synchronise. */
static struct memuse_log_entry *volatile memuse_log = NULL;
static volatile size_t memuse_log_size = 0;
static volatile size_t memuse_log_count = 0;
static volatile size_t memuse_old_count = 0;
static volatile size_t memuse_log_done = 0;

/* Current sum of memory in use. Only used in dumping. */
static size_t memuse_current = 0;

/* Label usage gathering struct. Only used in dumping. */
struct memuse_labelled_item {
  size_t sum;
  size_t count;
};

/* A radix node, this has a single byte key and a pointer to some related
 * resource. It also holds a sorted list of children, if any. */
struct memuse_rnode {

  /* Byte key of this node. */
  uint8_t keypart;

  /* Value of this node, if set. */
  void *ptr;

  /* Sorted pointers to children of this node. */
  struct memuse_rnode **children;
  unsigned int count;
};

/* Persistent radix trie root node. Holds active logs between dumps. */
static struct memuse_rnode *memuse_rnode_root;
static int memuse_rnode_root_init = 1;

#ifdef MEMUSE_RNODE_DUMP
/**
 * @brief Dump a representation of the radix tree rooted at a node to stdout.
 *
 * @param depth the depth of the node in the tree, root is 0.
 * @param node the node at which to start dumping.
 * @param full if not zero then nodes that are not storing a value
 *              are also reported.
 */
static void memuse_rnode_dump(int depth, struct memuse_rnode *node, int full) {

  /* Value of the full key, to this depth. Assumes full key is a pointer,
   * so uncomment when using strings. */
  static union {
    //uint8_t key[MEMUSE_MAXLABLEN];
    //char ptr[MEMUSE_MAXLABLEN];
    uint8_t key[sizeof(uintptr_t)];
    void *ptr;
  } keyparts = {0};

  /* Record keypart at this depth. Root has no keypart. */
  if (depth != 0) keyparts.key[depth - 1] = node->keypart;

  //if (node->ptr != NULL || full) {
  //  keyparts.key[depth] = '\0';
  //
  //    /* Gather children's keys if full. */
  //    char fullkey[MEMUSE_MAXLABLEN];
  //    if (full) {
  //      for (size_t k = 0; k < node->count; k++) {
  //        fullkey[k] = node->children[k]->keypart;
  //      }
  //      fullkey[node->count] = '\0';
  //      printf("dump @ depth: %d keypart: %d key: %s value: %p fullkey: %s\n",
  //             depth, node->keypart, keyparts.ptr, node->ptr, fullkey);
  //    } else {
  //      printf("dump @ depth: %d keypart: %d key: %s value: %p\n", depth,
  //             node->keypart, keyparts.ptr, node->ptr);
  //    }
  //}

  if (node->ptr != NULL || full) {
    printf("dump @ depth: %d keypart: %d key: %p value: %p\n", depth,
           node->keypart, keyparts.ptr, node->ptr);
  }

  /* Recurse to all children. */
  for (size_t k = 0; k < node->count; k++) {
    memuse_rnode_dump(depth + 1, node->children[k], full);
  }
}
#endif

/**
 * @brief Return the position of a keypart for a list of children.
 *        If not found returns where it would be inserted.
 *
 * @param keypart the keypart to locate.
 * @param children the list of sorted children.
 * @param count the number of children
 *
 * @return the index of key or where it should be inserted.
 */
static unsigned int memuse_rnode_bsearch(uint8_t keypart,
                                         struct memuse_rnode **children,
                                         unsigned int count) {

  /* Search for lower bound. */
  unsigned int lower = 0;
  unsigned int upper = count;
  while (lower < upper) {
    unsigned int middle = (upper + lower) / 2;
    if (keypart > children[middle]->keypart)
      lower = middle + 1;
    else
      upper = middle;
  }
  return lower;
}

/**
 * @brief Insert a child, if needed, into a list of children. Assumes
 *        we have sufficient room.
 *
 * @param child the child to insert, if needed.
 * @param children the list of sorted children.
 * @param count the number of children
 */
static void memuse_rnode_binsert_child(struct memuse_rnode *child,
                                       struct memuse_rnode **children,
                                       unsigned int *count) {
  unsigned int pos = 0;
  if (*count > 0) {

    /* Find the child or insertion point. */
    pos = memuse_rnode_bsearch(child->keypart, children, *count);

    /* If not found move all children to make a space, unless we're inserting
     * after the end. */
    if (pos < *count && children[pos]->keypart != child->keypart) {
      memmove(&children[pos + 1], &children[pos],
              (*count - pos) * sizeof(struct memuse_rnode *));
    }
  }

  /* Insert new child */
  children[pos] = child;
  *count += 1;
}

/**
 * @brief Add a child rnode to an rnode. Making sure we have room and keeping
 *        the sort order.
 *
 * @param node the parent node.
 * @param child the node to add to the parent,
 */
static void memuse_rnode_add_child(struct memuse_rnode *node,
                                   struct memuse_rnode *child) {

  /* XXX use a cache of allocations. */
  /* Extend the children list to include a new entry .*/
  void *mem = realloc(node->children,
                      (node->count + 1) * sizeof(struct memuse_rnode *));
  if (mem == NULL) error("Failed to reallocate rnodes\n");
  node->children = mem;

  /* Insert the new child. */
  memuse_rnode_binsert_child(child, node->children, &node->count);
}

/**
 * @brief Find a child of a node with the given key part.
 *
 * @param node the node to search.
 * @param keypart the key part of the child.
 * @return NULL if not found.
 */
static struct memuse_rnode *memuse_rnode_lookup(const struct memuse_rnode *node,
                                                uint8_t keypart) {

  /* Locate the key, or where it would be inserted. */
  if (node->count > 0) {
    unsigned int index = memuse_rnode_bsearch(keypart, node->children,
                                              node->count);
    if (index < node->count && keypart == node->children[index]->keypart) {
      return node->children[index];
    }
  }
  return NULL;
}

/**
 * @brief insert a child into a node's children list and add a pointer, iff
 *        this is the destination node for the given key.
 *
 * @param node the parent node.
 * @param depth the depth of the parent node.
 * @param key the full key of the eventual leaf node.
 * @param keylen the numbers of bytes in the full key.
 * @param value pointer that will be stored as the value of the leaf node.
 */
static void memuse_rnode_insert_child(struct memuse_rnode *node, uint8_t depth,
                                      uint8_t *key, uint8_t keylen,
                                      void *value) {

  /* Check if keypart this already exists at this level and add new child if
   * not. */
  uint8_t keypart = key[depth];
  struct memuse_rnode *child = memuse_rnode_lookup(node, keypart);
  if (child == NULL) {
    child = calloc(1, sizeof(struct memuse_rnode));
    child->keypart = keypart;
    memuse_rnode_add_child(node, child);
  }

  /* Are we at the lowest level yet? */
  depth++;
  if (depth == keylen) {
    /* Our destination node. */

#if SWIFT_DEBUG_CHECKS
    if (child->ptr != NULL)
      message("Overwriting rnode value: %p with %p", child->ptr, value);
#endif
    child->ptr = value;
    return;
  }

  /* Down we go to the next level. */
  memuse_rnode_insert_child(child, depth, key, keylen, value);
  return;
}

/**
 * @brief Find a child node for the given full key.
 *
 * @param node the current parent node.
 * @param depth the depth of the parent node, 0 for first call.
 * @param key the full key of the expected child node.
 * @param keylen the number of bytes in the key.
 */
static struct memuse_rnode *memuse_rnode_find_child(struct memuse_rnode *node,
                                                    uint8_t depth, uint8_t *key,
                                                    uint8_t keylen) {
  uint8_t keypart = key[depth];
  struct memuse_rnode *child = NULL;
  if (node->count > 0) child = memuse_rnode_lookup(node, keypart);
  if (child != NULL && (depth + 1) < keylen) {
    return memuse_rnode_find_child(child, depth + 1, key, keylen);
  }
  return child;
}

/**
 * @brief Free all resources associated with a node.
 *
 * @param node the rnode.
 */
static void memuse_rnode_cleanup(struct memuse_rnode *node) {

  if (!node) return;

  for (size_t k = 0; k < node->count; k++) {
    memuse_rnode_cleanup(node->children[k]);
    free(node->children[k]);
  }
  if (node->count > 0) free(node->children);
}

/**
 * @brief reallocate the entries log if space is needed.
 */
static void memuse_log_reallocate(size_t ind) {

  if (ind == 0) {

    /* Need to perform initialization. Be generous. */
    if ((memuse_log = (struct memuse_log_entry *)malloc(
             sizeof(struct memuse_log_entry) * MEMUSE_INITLOG)) == NULL)
      error("Failed to allocate memuse log.");

    /* Last action. */
    memuse_log_size = MEMUSE_INITLOG;

  } else {
    struct memuse_log_entry *new_log;
    if ((new_log = (struct memuse_log_entry *)malloc(
             sizeof(struct memuse_log_entry) *
             (memuse_log_size + MEMUSE_INITLOG))) == NULL)
      error("Failed to re-allocate memuse log.");

    /* Wait for all writes to the old buffer to complete. */
    while (memuse_log_done < memuse_log_size)
      ;

    /* Copy to new buffer. */
    memcpy(new_log, memuse_log,
           sizeof(struct memuse_log_entry) * memuse_log_size);
    free(memuse_log);
    memuse_log = new_log;

    /* Last action, releases waiting threads. */
    atomic_add(&memuse_log_size, MEMUSE_INITLOG);
  }
}

/**
 * @brief Log an allocation or deallocation of memory.
 *
 * @param label the label associated with the memory.
 * @param ptr the memory pointer.
 * @param allocated whether this is an allocation or deallocation.
 * @param size the size in byte of memory allocated, set to 0 when
 *             deallocating.
 */
void memuse_log_allocation(const char *label, void *ptr, int allocated,
                           size_t size) {

  size_t ind = atomic_inc(&memuse_log_count);

  /* If we are at the current size we need more space. */
  if (ind == memuse_log_size) memuse_log_reallocate(ind);

  /* Other threads wait for space. */
  while (ind > memuse_log_size)
    ;

  /* Record the log. */
  memuse_log[ind].step = engine_current_step;
  memuse_log[ind].allocated = allocated;
  memuse_log[ind].size = size;
  memuse_log[ind].ptr = ptr;
  strncpy(memuse_log[ind].label, label, MEMUSE_MAXLABLEN);
  memuse_log[ind].label[MEMUSE_MAXLABLEN] = '\0';
  memuse_log[ind].dtic = getticks() - clocks_start_ticks;
  memuse_log[ind].active = 1;
  atomic_inc(&memuse_log_done);
}

/**
 * @brief dump the log to a file and reset, if anything to dump.
 *
 * @param filename name of file for log dump.
 */
void memuse_log_dump(const char *filename) {

  /* Skip if nothing allocated this step. */
  if (memuse_log_count == memuse_old_count) return;

  //ticks tic = getticks();

  /* Create the radix tree. If not already done. */
  if (memuse_rnode_root_init) {
    memuse_rnode_root =
        (struct memuse_rnode *)calloc(1, sizeof(struct memuse_rnode));
    memuse_rnode_root_init = 0;
  }

  /* Stop any new logs from being processed while we are dumping.
   * Remember to not abort with error() in this section, that is recursive
   * with the exit handler. */
  size_t log_count = memuse_log_count;
  size_t old_count = memuse_old_count;

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL) {
    message("Failed to create memuse log file '%s', logs not dumped.",
            filename);
    return;
  }

  /* Write a header. */
  fprintf(fd, "# dtic step label size sum\n");

  size_t memuse_maxmem = memuse_current;
  for (size_t k = old_count; k < log_count; k++) {

    /* Check if this address has already been recorded. */
    struct memuse_rnode *child = memuse_rnode_find_child(
        memuse_rnode_root, 0, memuse_log[k].vptr, sizeof(uintptr_t));

    if (child != NULL && child->ptr != NULL) {

      /* Found the allocation, this should be the free. */
      if (memuse_log[k].allocated) {

      /* Allocated twice, this is an error, but we cannot abort as that will
       * attempt another memory dump, so just complain. */
#if SWIFT_DEBUG_CHECKS
        message("Allocated the same address twice (%s: %zd)",
                memuse_log[k].label, memuse_log[k].size);
#endif
        continue;
      }

      /* Free, update the size to remove the allocation. */
      struct memuse_log_entry *oldlog = (struct memuse_log_entry *)child->ptr;
      memuse_log[k].size = -oldlog->size;

      /* And deactivate this key. */
      child->ptr = NULL;

      /* And mark this as matched. */
      memuse_log[k].active = 0;
      oldlog->active = 0;

    } else if (child == NULL && memuse_log[k].allocated) {

      /* Not found, so new allocation which we store the log against the
       * address. */
      memuse_rnode_insert_child(memuse_rnode_root, 0, memuse_log[k].vptr,
                                sizeof(uintptr_t), &memuse_log[k]);

    } else if (child == NULL && !memuse_log[k].allocated) {

    /* Unmatched free, OK if NULL. */
#if SWIFT_DEBUG_CHECKS
      if (memuse_log[k].ptr != NULL) {
        message("Unmatched non-NULL free: %s", memuse_log[k].label);
      }
#endif
      continue;
    } else if (memuse_log[k].allocated) {

      /* Must be previously released allocation with same address, so we
       * store. */
      memuse_rnode_insert_child(memuse_rnode_root, 0, memuse_log[k].vptr,
                                sizeof(uintptr_t), &memuse_log[k]);

    } else {
      /* Should not happen ... */
      message("weird memory log record for label '%s' skipped",
              memuse_log[k].label);
      continue;
    }

    /* Keep maximum and rolling sum. */
    memuse_current += memuse_log[k].size;
    if (memuse_current > memuse_maxmem) memuse_maxmem = memuse_current;

    /* And output. */
    fprintf(fd, "%lld %d %s %zd %zd\n", memuse_log[k].dtic,
            memuse_log[k].step, memuse_log[k].label,
            memuse_log[k].size, memuse_current);
  }

#ifdef MEMUSE_RNODE_DUMP
  /* Debug dump of tree. */
  //memuse_rnode_dump(0, memuse_rnode_root, 0);
#endif

  /* Now we find all the still active nodes and gather their sizes against the
   * labels. */
  struct memuse_rnode *activernodes =
      (struct memuse_rnode *)calloc(1, sizeof(struct memuse_rnode));
  size_t newcount = 0;
  struct memuse_rnode *labellednodes =
      (struct memuse_rnode *)calloc(1, sizeof(struct memuse_rnode));
  size_t *lindices = (size_t *)calloc(log_count, sizeof(size_t));
  size_t lcount = 0;
  for (size_t k = 0; k < log_count; k++) {

    /* Only allocations are stored also is it active? */
    if (memuse_log[k].allocated && memuse_log[k].active) {

      /* Look for this label in our tree. */
      struct memuse_rnode *labelchild = memuse_rnode_find_child(
          labellednodes, 0, (uint8_t *)memuse_log[k].label,
          strlen(memuse_log[k].label));
      struct memuse_labelled_item *item = NULL;
      if (labelchild == NULL || labelchild->ptr == NULL) {

        /* New, so create an instance to keep the count. */
        item = (struct memuse_labelled_item *)calloc(
            1, sizeof(struct memuse_labelled_item));
        item->sum = 0;
        item->count = 0;
        memuse_rnode_insert_child(labellednodes, 0,
                                  (uint8_t *)memuse_log[k].label,
                                  strlen(memuse_log[k].label), item);

        /* Keep for indexing next time. */
        lindices[lcount] = newcount;
        lcount++;
      } else {
        item = (struct memuse_labelled_item *)labelchild->ptr;
      }

      /* And increment sum. */
      item->sum += memuse_log[k].size;
      item->count++;

      /* Keep this in new log entry tree. Move to head. */
      memcpy(&memuse_log[newcount], &memuse_log[k],
             sizeof(struct memuse_log_entry));
      memuse_rnode_insert_child(activernodes, 0, memuse_log[newcount].vptr,
                                sizeof(uintptr_t), &memuse_log[newcount]);
      newcount++;
    }
  }

  /* And move all active logs to a clean new tree for next time. */
  memuse_log_count = newcount;
  memuse_old_count = newcount;
  memuse_rnode_cleanup(memuse_rnode_root);
  free(memuse_rnode_root);
  memuse_rnode_root = activernodes;

  /* Now dump the labelled counts. */
  fprintf(fd, "# Memory use by label:\n");
  fprintf(fd, "##  %30s %16s %16s\n", "label", "MB", "numactive");
  fprintf(fd, "##\n");

  size_t total_mem = 0;
  for (size_t k = 0; k < lcount; k++) {
    size_t ind = lindices[k];

    /* Find this entry. */
    struct memuse_rnode *labelchild = memuse_rnode_find_child(
        labellednodes, 0, (uint8_t *)memuse_log[ind].label,
        strlen(memuse_log[ind].label));
    struct memuse_labelled_item *item =
        (struct memuse_labelled_item *)labelchild->ptr;
    fprintf(fd, "## %30s %16.3f %16zd\n", memuse_log[ind].label,
            item->sum / MEGABYTE, item->count);
    total_mem += item->sum;

    /* Don't need this again. */
    free(item);
  }
  fprintf(fd, "##\n");
  fprintf(fd, "# Total memory still in use : %.3f (MB)\n",
          total_mem / MEGABYTE);
  fprintf(fd, "# Peak memory usage         : %.3f (MB)\n",
          memuse_maxmem / MEGABYTE);
  fprintf(fd, "#\n");
  fprintf(fd, "# Memory use by process (all/system): %s\n", memuse_process(1));
  fprintf(fd, "# cpufreq: %lld\n", clocks_get_cpufreq());

  /* Clean up tree. */
  memuse_rnode_cleanup(labellednodes);
  free(labellednodes);
  free(lindices);

  /* Close the file. */
  fflush(fd);
  fclose(fd);

  //message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
  //        clocks_getunit());
}

/**
 * @brief dump the log for using the given rank to generate a standard
 *        name for the output. Used when exiting in error.
 *
 * @param rank the rank exiting in error.
 */
void memuse_log_dump_error(int rank) {
  char filename[60];
  sprintf(filename, "memuse-error-report-rank%d.txt", rank);
  memuse_log_dump(filename);
}

#endif /* SWIFT_MEMUSE_REPORTS */

/**
 * @brief parse the process /proc/self/statm file to get the process
 *        memory use (in KB). Top field in ().
 *
 * @param size     total virtual memory (VIRT/VmSize)
 * @param resident resident non-swapped memory (RES/VmRSS)
 * @param shared   shared (mmap'd) memory  (SHR, RssFile+RssShmem)
 * @param text     text (exe) resident set (CODE, note also includes data
 *                 segment, so is considered broken for Linux)
 * @param data     data+stack resident set (DATA, note also includes library,
 *                 so is considered broken for Linux)
 * @param library  library resident set (0 for Linux)
 * @param dirty    dirty pages (nDRT = 0 for Linux)
 */
void memuse_use(long *size, long *resident, long *shared, long *text,
                long *data, long *library, long *dirty) {

  /* Open the file. */
  FILE *file = fopen("/proc/self/statm", "r");
  if (file != NULL) {
    int nscan = fscanf(file, "%ld %ld %ld %ld %ld %ld %ld", size, resident,
                       shared, text, library, data, dirty);

    if (nscan == 7) {
      /* Convert pages into bytes. Usually 4096, but could be 512 on some
       * systems so take care in conversion to KB. */
      long sz = sysconf(_SC_PAGESIZE);
      *size *= sz;
      *resident *= sz;
      *shared *= sz;
      *text *= sz;
      *library *= sz;
      *data *= sz;
      *dirty *= sz;

      *size /= 1024;
      *resident /= 1024;
      *shared /= 1024;
      *text /= 1024;
      *library /= 1024;
      *data /= 1024;
      *dirty /= 1024;
    } else {
      error("Failed to read sufficient fields from /proc/self/statm");
    }
    fclose(file);
  } else {
    error("Failed to open /proc/self/statm");
  }
}

/**
 * @brief Return a string with the current memory use of the process described.
 *
 * Not thread safe.
 *
 * @param inmb if true then report in MB, not KB.
 *
 * @result the memory use of the process, note make a copy if not used
 * immediately.
 */
const char *memuse_process(int inmb) {
  static char buffer[256];
  long size;
  long resident;
  long shared;
  long text;
  long library;
  long data;
  long dirty;
  memuse_use(&size, &resident, &shared, &text, &data, &library, &dirty);

  if (inmb) {
    snprintf(buffer, 256,
             "VIRT = %.3f SHR = %.3f CODE = %.3f DATA = %.3f "
             "RES = %.3f (MB)",
             size / 1024.0, shared / 1024.0, text / 1024.0, data / 1024.0,
             resident / 1024.0);
  } else {
    snprintf(buffer, 256,
             "VIRT = %ld SHR = %ld CODE = %ld DATA = %ld "
             "RES = %ld (KB)",
             size, shared, text, data, resident);
  }
  return buffer;
}
