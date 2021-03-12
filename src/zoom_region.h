#include "cell.h"
#include "space.h"

#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

void zoom_region_init(struct swift_params *params, struct space *s);
void construct_zoom_region(struct space *s, int verbose);

#endif /* SWIFT_ZOOM_H */
