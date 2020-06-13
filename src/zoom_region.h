#include "space.h"
#include "cell.h"

void zoom_region_init(struct swift_params *params, struct space *s);
void construct_zoom_region(struct space *s, int verbose);
void check_zoom_region(const struct space *s, const int verbose);
int cell_getid_zoom(const int cdim[3], const double x, const double y, const double z,
    const struct zoom_region_properties *zoom_props,
    const int i, const int j, const int k);
void makeproxies_between_top_levels(struct engine *e);
double cell_min_dist2(
            const struct cell *restrict ci, const struct cell *restrict cj,
                const int periodic, const double dim[3]);
double cell_min_dist2_diff_size(
            const struct cell *restrict ci, const struct cell *restrict cj,
                const int periodic, const double dim[3]);
void engine_make_self_gravity_tasks_mapper_zoom(void *map_data, int num_elements,
                                                   void *extra_data);
void engine_make_self_gravity_tasks_mapper_between_toplevels(void *map_data, int num_elements,
                                                   void *extra_data);
void find_neighbouring_cells(struct space *s, const int verbose);
