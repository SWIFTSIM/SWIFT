#include "cell.h"
#include "gravity_properties.h"
#include "space.h"

#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

void zoom_region_init(struct swift_params *params, struct space *s);
int cell_getid_zoom(const struct space *s, const double x, const double y,
                    const double z);
void construct_zoom_region(struct space *s, int verbose);
void construct_tl_cells_with_zoom_region(struct space *s, const int *cdim, const float dmin, 
        const integertime_t ti_current, struct gravity_props *gravity_properties, int verbose);
void find_neighbouring_cells(struct space *s, struct gravity_props *gravity_properties, const int verbose);
double cell_min_dist2_diff_size(const struct cell *restrict ci,
                                const struct cell *restrict cj,
                                const int periodic, const double dim[3]);
double cell_min_dist2(const struct cell *restrict ci,
                      const struct cell *restrict cj, const int periodic,
                      const double dim[3]);
void engine_makeproxies_with_natural_cells(struct engine *e);
void engine_makeproxies_with_zoom_cells(struct engine *e);
void engine_makeproxies_with_between_grids(struct engine *e);
void engine_makeproxies_with_zoom_region(struct engine *e);
void engine_make_self_gravity_tasks_mapper_natural_cells(void *map_data,
																										     int num_elements,
																										     void *extra_data);
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data);
void engine_make_self_gravity_tasks_mapper_with_zoom_diffsize(void *map_data,
                                                              int num_elements,
                                                              void *extra_data);
void engine_make_hydroloop_tasks_mapper_with_zoom(void *map_data, int num_elements,
                                                  void *extra_data);
void engine_make_fofloop_tasks_mapper_with_zoom(void *map_data, int num_elements,
                                                void *extra_data);
/* Parition prototypes */

int partition_space_to_space_zoom(double *oldh, double *oldcdim, double *oldzoomh,
		                              double *oldzoomcdim, int *oldnodeIDs, struct space *s);
void pick_vector_zoom(struct space *s, int nregions, int *samplecells);
void split_vector_zoom(struct space *s, int nregions, int *samplecells);

/* Regrid prototypes */

void space_regrid_zoom(struct space *s, struct gravity_props *p, int verbose);

#endif /* SWIFT_ZOOM_H */
