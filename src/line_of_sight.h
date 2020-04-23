/* Config parameters. */
#include "../config.h"
#include "engine.h"
#include "io_properties.h"

struct line_of_sight {
  int xaxis, yaxis, zaxis;
  double Xpos, Ypos;
  int particles_in_los_total;
  int particles_in_los_local;
};

struct los_props {
  int num_along_xy;
  int num_along_yz;
  int num_along_xz;

  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;

  int num_tot;
};

double los_periodic(double x, double dim);
void generate_line_of_sights(struct line_of_sight *Los,
                             const struct los_props *params);
void print_los_info(const struct line_of_sight *Los,
                    const struct los_props *params);
void do_line_of_sight(struct engine *e);
void los_init(double dim[3], struct los_props *los_params,
        struct swift_params *params);
void write_los_hdf5_datasets(hid_t grp, int j, int N, struct part* parts,
                struct engine* e);
void write_los_hdf5_dataset(const struct io_props p, int N, int j, struct engine* e, hid_t grp);
