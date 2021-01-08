#ifndef SWIFT_CELL_SHADOWFAX_H
#define SWIFT_CELL_SHADOWFAX_H

#include "cell.h"
#include "delaunay.h"

__attribute__((always_inline)) INLINE static void
cell_malloc_delaunay_tesselation(struct cell *c, const struct hydro_space *hs) {

  delaunay_init(&c->hydro.deltess, hs, 100, 100);
}

#endif /* SWIFT_CELL_SHADOWFAX_H */
