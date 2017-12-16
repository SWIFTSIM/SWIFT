/***********************************************************************
/
/ Grackle c wrapper
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include "grackle_wrapper.h"

#ifdef SWIFT_DEBUG_CHECKS
#include <assert.h>
#define GRACKLE_ASSERT(X) assert((X))
#else
#define GRACKLE_ASSERT(X)
#endif

code_units my_units;

// arrays passed to grackle as input and to be filled
#define FIELD_SIZE 1

gr_float HDI_density[FIELD_SIZE];

// Set grid dimension and size.
// grid_start and grid_end are used to ignore ghost zones.
const int grid_rank = 1;

int wrap_init_cooling(char *CloudyTable, int UVbackground, double udensity,
                      double ulength, double utime, int grackle_chemistry) {

#ifdef SWIFT_DEBUG_CHECKS
  grackle_verbose = 1;
#endif
  double velocity_units;

  // First, set up the units system.
  // These are conversions from code units to cgs.

  my_units.a_units = 1.0;  // units for the expansion factor (1/1+zi)

  my_units.comoving_coordinates =
      0; /* so, according to the doc, we assume here all physical quantities to
            be in proper coordinate (not comobile)  */
  my_units.density_units = udensity;
  my_units.length_units = ulength;
  my_units.time_units = utime;
  velocity_units =
      my_units.a_units * my_units.length_units / my_units.time_units;
  my_units.velocity_units = velocity_units;

  // Second, create a chemistry object for parameters and rate data.
  if (set_default_chemistry_parameters() == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry.
  grackle_data.use_grackle = 1;
  grackle_data.with_radiative_cooling = 1;
  grackle_data.primordial_chemistry =
      grackle_chemistry;           // molecular network with H, He, D
  grackle_data.metal_cooling = 1;  // metal cooling on
  grackle_data.UVbackground = UVbackground;
  grackle_data.grackle_data_file = CloudyTable;

  // Finally, initialize the chemistry object.
  // snl: a_value is not the true initial ae!! This should get set during
  // update_UVbackground
  double initial_redshift = 0.;
  double a_value = 1. / (1. + initial_redshift);

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units, a_value) == 0) {
    error("Error in initialize_chemistry_data.");
  }

  return 1;
}

int wrap_set_UVbackground_on() {
  // The UV background rates is enabled
  grackle_data.UVbackground = 1;
  return 1;
}

int wrap_set_UVbackground_off() {
  // The UV background rates is disabled
  grackle_data.UVbackground = 0;
  return 1;
}

int wrap_get_cooling_time(double rho, double u, double Z, double a_now,
                          double *coolingtime) {
  gr_float den_factor = 1.0;
  gr_float u_factor = 1.0;

  gr_float x_velocity[FIELD_SIZE] = {0.0};
  gr_float y_velocity[FIELD_SIZE] = {0.0};
  gr_float z_velocity[FIELD_SIZE] = {0.0};

  gr_float density[FIELD_SIZE] = {rho * den_factor};
  gr_float metal_density[FIELD_SIZE] = {Z * density[0]};
  gr_float energy[FIELD_SIZE] = {u * u_factor};

  gr_float cooling_time[FIELD_SIZE] = {0.0};

  int grid_dimension[3] = {1, 0, 0};
  int grid_start[3] = {0, 0, 0};
  int grid_end[3] = {0, 0, 0};

  if (FIELD_SIZE != 1) {
    error("field_size must currently be set to 1.");
  }

  if (calculate_cooling_time_table(&my_units, a_now, grid_rank, grid_dimension,
                                   grid_start, grid_end, density, energy,
                                   x_velocity, y_velocity, z_velocity,
                                   metal_density, cooling_time) == 0) {
    error("Error in calculate_cooling_time.");
  }

  // return updated chemistry and energy
  for (int i = 0; i < FIELD_SIZE; i++) {
    coolingtime[i] = cooling_time[i];
  }

  return 1;
}

int wrap_do_cooling(double rho, double *u, double dt, double Z, double a_now) {

  GRACKLE_ASSERT(FIELD_SIZE == 1);

  gr_float den_factor = 1.0;
  gr_float u_factor = 1.0;

  gr_float x_velocity[FIELD_SIZE] = {0.0};
  gr_float y_velocity[FIELD_SIZE] = {0.0};
  gr_float z_velocity[FIELD_SIZE] = {0.0};

  gr_float density[FIELD_SIZE] = {rho * den_factor};
  gr_float metal_density[FIELD_SIZE] = {Z * density[0]};
  gr_float energy[FIELD_SIZE] = {(*u) * u_factor};

  int grid_dimension[3] = {1, 0, 0};
  int grid_start[3] = {0, 0, 0};
  int grid_end[3] = {0, 0, 0};

  if (solve_chemistry_table(&my_units, a_now, dt, grid_rank, grid_dimension,
                            grid_start, grid_end, density, energy, x_velocity,
                            y_velocity, z_velocity, metal_density) == 0) {
    error("Error in solve_chemistry.");
    return 0;
  }
  // return updated chemistry and energy
  for (int i = 0; i < FIELD_SIZE; i++) {
    u[i] = energy[i] / u_factor;
  }

  return 1;
}

void grackle_print_data() {
  message("Grackle Data:");
  message("\t Data file: %s", grackle_data.grackle_data_file);
  message("\t With grackle: %i", grackle_data.use_grackle);
  message("\t With radiative cooling: %i", grackle_data.with_radiative_cooling);
  message("\t With UV background: %i", grackle_data.UVbackground);
  message("\t With primordial chemistry: %i",
          grackle_data.primordial_chemistry);
  message("\t Number temperature bins: %i",
          grackle_data.NumberOfTemperatureBins);
  message("\t T = (%g, ..., %g)", grackle_data.TemperatureStart,
          grackle_data.TemperatureEnd);

  message("Primordial Cloudy");
  cloudy_print_data(grackle_data.cloudy_primordial, 1);
  if (grackle_data.metal_cooling) {
    message("Metal Cooling");
    cloudy_print_data(grackle_data.cloudy_metal, 0);
  }

  message("\t Gamma: %g", grackle_data.Gamma);

  /* UVB */
  if (grackle_data.UVbackground && grackle_data.primordial_chemistry != 0) {
    struct UVBtable uvb = grackle_data.UVbackground_table;
    long long N = uvb.Nz;
    message("\t UV Background");
    message("\t\t Redshift from %g to %g with %lli steps", uvb.zmin, uvb.zmax,
            N);
    message("\t\t z = (%g, ..., %g)", uvb.z[0], uvb.z[N - 1]);
  }
}

void cloudy_print_data(const cloudy_data c, const int print_mmw) {
  long long N = c.data_size;
  message("\t Data size: %lli", N);
  message("\t Grid rank: %lli", c.grid_rank);

  char msg[200] = "\t Dimension: (";
  for (long long i = 0; i < c.grid_rank; i++) {
    char tmp[200] = "%lli%s";
    if (i == c.grid_rank - 1)
      sprintf(tmp, tmp, c.grid_dimension[i], ")");
    else
      sprintf(tmp, tmp, c.grid_dimension[i], ", ");

    strcat(msg, tmp);
  }
  message("%s", msg);

  if (c.heating_data)
    message("\t Heating: (%g, ..., %g)", c.heating_data[0],
            c.heating_data[N - 1]);
  if (c.cooling_data)
    message("\t Cooling: (%g, ..., %g)", c.cooling_data[0],
            c.cooling_data[N - 1]);
  if (c.mmw_data && print_mmw)
    message("\t Mean molecular weigth: (%g, ..., %g)", c.mmw_data[0],
            c.mmw_data[N - 1]);
}
