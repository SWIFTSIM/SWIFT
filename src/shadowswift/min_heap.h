#ifndef SHADOWSWIFT_MIN_HEAP_H
#define SHADOWSWIFT_MIN_HEAP_H

#include "bvh.h"

#include <stdlib.h>

typedef struct bboxtree_node {
  struct bboxtree_node* progeny[8];
  bbox_t bbox;
  struct part* parts;
  int count;
} bboxtree_node_t;

typedef struct heap_node {
  double distance2;

  struct {
    int sid;
    int offset;
    bboxtree_node_t* bboxtree_node;
  } data;
} heap_node_t;

inline static int heap_node_is_leaf(heap_node_t* self) {
  return self->data.bboxtree_node == NULL;
}

inline static void heap_node_init(heap_node_t* self, double distance2, int sid,
                                  int offset, bboxtree_node_t* bboxtree_node) {
  self->distance2 = distance2;
  self->data.sid = sid;
  self->data.offset = offset;
  self->data.bboxtree_node = bboxtree_node;
}

inline static heap_node_t* heap_node_new(double distance2, int sid, int offset,
                                         bboxtree_node_t* bboxtree_node) {
  heap_node_t* node = (heap_node_t*)malloc(sizeof(*node));
  heap_node_init(node, distance2, sid, offset, bboxtree_node);
}

inline static void heap_node_destroy(heap_node_t* self) { free(self); }

inline static void swap_nodes(heap_node_t** a, heap_node_t** b) {
  heap_node_t* temp = *a;
  *a = *b;
  *b = temp;
}

typedef struct min_heap {
  heap_node_t** heap;

  int size;

  int capacity;
} min_heap_t;

inline static void min_heap_sift_up(min_heap_t* self);
inline static void min_heap_sift_down(min_heap_t* self);

inline static min_heap_t* min_heap_new(int capacity) {
  min_heap_t* min_heap = (min_heap_t*)malloc(sizeof(*min_heap));
  min_heap->heap = (heap_node_t**)malloc(capacity * sizeof(*min_heap->heap));
  min_heap->size = 0;
  min_heap->capacity = capacity;
  return min_heap;
}

inline static void min_heap_destroy(min_heap_t* self) {
  free(self->heap);
  free(self);
}

/**
 * @brief Push a new node onto the `min_heap` and restore the min_heap property
 * afterwards.
 */
inline static void min_heap_push(min_heap_t* self, heap_node_t* node) {
  // grow if necessary
  if (self->size == self->capacity) {
    self->capacity = self->capacity << 1;
    self->heap = (heap_node_t**)realloc(self->heap,
                                        self->capacity * sizeof(*self->heap));
  }

  // Add new node to end and sift up
  self->heap[self->size] = node;
  self->size++;
  min_heap_sift_up(self);
}

/**
 * @brief Pop the first (minimal) element of the `min_heap` and restore the
 * min_heap property afterwards.
 */
inline static heap_node_t* min_heap_pop(min_heap_t* self) {
  heap_node_t* min_node = self->heap[0];

  // Swap the last element to the front and sift down
  swap_nodes(&self->heap[0], &self->heap[--self->size]);
  min_heap_sift_down(self);

  return min_node;
}

/**
 * @brief Peek at the first (minimal) element of the `min_heap` without 
 * modifying it.
 */
inline static heap_node_t* min_heap_peek(const min_heap_t* self) {
  return self->heap[0];
}

/**
 * @brief Retrieve the first (minimal) element of the `min_heap` and replace it 
 * with a new element. Restore the min_heap property afterwards.
 * 
 * This is more efficient than first popping and then pushing (which would 
 * restore the min-heap property twice).
 */
inline static heap_node_t* min_heap_replace(min_heap_t* self, heap_node_t* node) {
  heap_node_t* min_node = self->heap[0];
  self->heap[0] = node;
  min_heap_sift_down(self);
  return min_node;
}

/**
 * @brief Sift the last element of this 'min_heap' upwards until it is at the
 * right level.
 */
inline static void min_heap_sift_up(min_heap_t* self) {
  int index = self->size - 1;
  int parent_index = (index - 1) / 2;

  while (index > 0 &&
         self->heap[index]->distance2 < self->heap[parent_index]->distance2) {
    swap_nodes(&self->heap[index], &self->heap[parent_index]);
    index = parent_index;
    parent_index = (index - 1) / 2;
  }
}

/**
 * @brief Sift the first element of this `min_heap` downwards until it is at the
 * right level
 */
inline static void min_heap_sift_down(min_heap_t* self) {
  int index = 0;
  int left_child_index = 1;
  int right_child_index = 2;

  while (index < self->size - 2) {
    if (self->heap[index]->distance2 >
        self->heap[left_child_index]->distance2) {
      swap_nodes(&self->heap[index], &self->heap[left_child_index]);
      index = left_child_index;
    } else if (self->heap[index]->distance2 >
               self->heap[right_child_index]->distance2) {
      swap_nodes(&self->heap[index], &self->heap[right_child_index]);
      index = right_child_index;
    } else {
      // Node is at correct level
      break;
    }
    left_child_index = index * 2 + 1;
    right_child_index = index * 2 + 2;
  }
}

/**
 * A struct used to iterate over the nearest neighbours of a certain site `loc`
 *
 * Iteration happens from closest to furthest and uses a priority queue/min_heap
 * under the hood.
 */
typedef struct nearest_neighbour_iterator {
  const double loc[3];
  min_heap_t* min_heap;
} nearest_neighbour_iterator_t;

/**
 * @brief Construct a new `nearest_neighbour_iterator` struct.
 *
 * @param loc The query site
 * @param tree_nodes Array of `bboxtrees` containing the candidate nearest
 * neighbours.
 * @param sid Array of sid's of the corresponding bboxtrees.
 * @param size Length of the `tree_nodes` and `sid` arrays.
 * @param capacity Capacity to reserve for the internal min_heap.
 * @return A pointer to a new `nearest_neighbour_iterator` struct.
 */
inline static nearest_neighbour_iterator_t* nearest_neighbour_iterator_new(
    const double loc[3], bboxtree_node_t* tree_nodes, int* sid, int size,
    int capacity) {

  nearest_neighbour_iterator_t* nearest_neighbour_iterator =
      (nearest_neighbour_iterator_t*)malloc(
          sizeof(*nearest_neighbour_iterator));
  nearest_neighbour_iterator->min_heap = min_heap_new(capacity);

  for (int i = 0; i < size; i++) {
    bboxtree_node_t* tree_node = &tree_nodes[i];
    double distance2 = bbox_distance2(&tree_node->bbox, loc);
    heap_node_t* heap_node = heap_node_new(distance2, sid[i], 0, tree_node);
    min_heap_push(nearest_neighbour_iterator->min_heap, heap_node);
  }

  return nearest_neighbour_iterator;
}

inline static void nearest_neighbour_iterator_destroy(
    nearest_neighbour_iterator_t* self) {
  // Destroy any nodes left in the min_heap
  for (int i = 0; i < self->min_heap->size; i++) {
    heap_node_destroy(self->min_heap->heap[i]);
  }
  // Destroy the min_heap itself
  min_heap_destroy(self->min_heap);
  free(self);
}

/**
 * @brief Retrieve the index and sid of the next nearest neighbour of the query
 * site.
 *
 * @param self The `nearest_neighbour_iterator` to advance.
 * @param sid (Return) The sid of the next nearest neighbour.
 * @return The index of the next nearest neighbour in it's cell.
 */
inline static int nearest_neighbour_iterator_next(
    nearest_neighbour_iterator_t* self, /*return*/ int* sid) {

  heap_node_t* next = min_heap_peek(self->min_heap);
  while (!heap_node_is_leaf(next)) {
    // Add children of next to heap
    bboxtree_node_t* tree_node = next->data.bboxtree_node;
    int this_sid = next->data.sid;
    int offset = next->data.offset;
    if (tree_node->parts == NULL) {
      // Add subtrees
      for (int i = 0; i < 8; i++) {
        bboxtree_node_t* subtree = tree_node->progeny[i];
        double distance2 = bbox_distance2(&subtree->bbox, self->loc);
        heap_node_t* new_node =
            heap_node_new(distance2, this_sid, offset, subtree);
        if (i == 0) {
          // replace the root of the min_heap with the new node and heapify.
          min_heap_replace(self->min_heap, new_node);
        } else {
          // push the new node onto the heap and heapify.
          min_heap_push(self->min_heap, new_node);
        }
        offset += subtree->count;
      }
    } else {
      // Add individual parts
      for (int i = 0; i < tree_node->count; i++) {
        const struct part* part = &tree_node->parts[i];
        double dx[3] = {part->x[0] - self->loc[0], part->x[1] - self->loc[1],
                        part->x[2] - self->loc[2]};
        double distance2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        heap_node_t* new_node =
            heap_node_new(distance2, this_sid, offset, NULL);
        if (i == 0) {
          // replace the root of the min_heap with the new node and heapify.
          min_heap_replace(self->min_heap, new_node);
        } else {
          // push the new node onto the heap and heapify.
          min_heap_push(self->min_heap, new_node);
        }
        offset++;
      }
    }

    // destroy the old root node, we will not need it anymore.
    heap_node_destroy(next);

    // peek at new head
    next = min_heap_peek(self->min_heap);
  }

  int offset = next->data.offset;
  *sid = next->data.sid;

  // Destroy next, we will not need it anymore.
  heap_node_destroy(next);

  return offset;
}

#endif