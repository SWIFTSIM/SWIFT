#include "cell.h"
#include "space.h"

#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

void zoom_region_init(struct swift_params *params, struct space *s);
void construct_zoom_region(struct space *s, int verbose);
void construct_tl_cells_with_zoom_region(struct space *s, const int *cdim, const float dmin, 
        const integertime_t ti_current, const int verbose);
void find_neighbouring_cells(struct space *s, const int verbose);

#endif /* SWIFT_ZOOM_H */
