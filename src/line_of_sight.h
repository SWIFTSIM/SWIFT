/* Config parameters. */
#include "../config.h"

struct line_of_sight {
  int xaxis, yaxis, zaxis;
  double Xpos, Ypos;
  int particles_in_los_total;
  int particles_in_los_local;
};

struct line_of_sight_particles {
  double pos[3];
  float h;
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
                             struct los_props *params);
void print_los_info(struct line_of_sight *Los,
                    struct los_props *params);
void do_line_of_sight(struct engine *e);
void los_init(double dim[3], struct los_props *los_params,
        struct swift_params *params);
