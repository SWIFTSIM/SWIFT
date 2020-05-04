/* Config parameters. */
#include "../config.h"
#include "engine.h"
#include "io_properties.h"

struct line_of_sight {
  int xaxis, yaxis, zaxis;
  double Xpos, Ypos;
  size_t particles_in_los_total;
  size_t particles_in_los_local;
  int periodic;
  double dim[3];
};

struct los_props {
  int num_along_xy;
  int num_along_yz;
  int num_along_xz;

  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;

  int num_tot;

  char basename[200];
};

double los_periodic(double x, double dim);
void generate_line_of_sights(struct line_of_sight *Los,
                             const struct los_props *params,
                             const int periodic, const double dim[3]);
void print_los_info(const struct line_of_sight *Los, const int i);
void do_line_of_sight(struct engine *e);
void los_init(double dim[3], struct los_props *los_params,
        struct swift_params *params);
void write_los_hdf5_datasets(hid_t grp, int j, size_t N, const struct part* parts,
                struct engine* e, const struct xpart* xparts);
void write_los_hdf5_dataset(const struct io_props p, size_t N, int j, struct engine* e, hid_t grp);
void write_hdf5_header(hid_t h_file, const struct engine *e, const struct los_props* LOS_params);
void create_line_of_sight(const double Xpos, const double Ypos,   
        const int xaxis, const int yaxis, const int zaxis,
        const int periodic, const double dim[3], struct line_of_sight *los);
void los_struct_dump(const struct los_props *internal_los,
                            FILE *stream);
void los_struct_restore(const struct los_props *internal_los,
                                       FILE *stream);
