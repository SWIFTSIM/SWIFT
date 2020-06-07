#include "space.h"

void zoom_region_init(struct swift_params *params, struct space *s);
void construct_zoom_region(struct space *s, int verbose);
void check_zoom_region(const struct space *s, const int verbose);
int cell_getid_zoom(const int cdim[3], const double x, const double y, const double z,
    const struct zoom_region_properties *zoom_props,
    const int i, const int j, const int k);
void makeproxies_between_top_levels(struct engine *e);
