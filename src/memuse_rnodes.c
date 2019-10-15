/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
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
 *  @file memuse_rnodes.c
 *  @brief file of routines used for radix nodes in memory loggers.
 */

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "memuse_rnodes.h"

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"

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

  /* Extend the children list to include a new entry .*/
  void *mem = realloc(node->children,
                      (node->count + 1) * sizeof(struct memuse_rnode *));
  if (mem == NULL) error("Failed to reallocate rnodes\n");
  node->children = (struct memuse_rnode **)mem;

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
    unsigned int index =
        memuse_rnode_bsearch(keypart, node->children, node->count);
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
void memuse_rnode_insert_child(struct memuse_rnode *node, uint8_t depth,
                               uint8_t *key, uint8_t keylen, void *value) {

  /* Check if keypart this already exists at this level and add new child if
   * not. */
  uint8_t keypart = key[depth];
  struct memuse_rnode *child = memuse_rnode_lookup(node, keypart);
  if (child == NULL) {
    child = (struct memuse_rnode *)calloc(1, sizeof(struct memuse_rnode));
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
struct memuse_rnode *memuse_rnode_find_child(struct memuse_rnode *node,
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
void memuse_rnode_cleanup(struct memuse_rnode *node) {

  if (!node) return;

  for (size_t k = 0; k < node->count; k++) {
    memuse_rnode_cleanup(node->children[k]);
    free(node->children[k]);
  }
  if (node->count > 0) free(node->children);
}

/**
 * @brief Dump a representation of the radix tree rooted at a node to stdout.
 *
 * Debugging code.
 *
 * @param depth the depth of the node in the tree, root is 0.
 * @param node the node at which to start dumping.
 * @param full if not zero then nodes that are not storing a value
 *              are also reported.
 */
void memuse_rnode_dump(int depth, struct memuse_rnode *node, int full) {

  /* Value of the full key, to this depth. Assumes full key is a pointer,
   * so uncomment when using strings. */
  static union {
    // uint8_t key[MEMUSE_MAXLABLEN];
    // char ptr[MEMUSE_MAXLABLEN];
    uint8_t key[sizeof(uintptr_t)];
    void *ptr;
  } keyparts = {{0}};

  /* Record keypart at this depth. Root has no keypart. */
  if (depth != 0) keyparts.key[depth - 1] = node->keypart;

  // if (node->ptr != NULL || full) {
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
