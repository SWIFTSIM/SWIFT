//
// Created by yuyttenh on 11/01/23.
//

#ifndef SWIFTSIM_COOLING_TABLES_H
#define SWIFTSIM_COOLING_TABLES_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "cooling_properties.h"
#include "error.h"
#include "memuse.h"

/* System includes. */
#include <hdf5.h>
#include <string.h>

#define de_rijcke_cooling_N_temperatures 351

static void get_cooling_tables(struct cooling_function_data *restrict cooling) {

  char fname[de_rijcke_table_path_name_length + 16];
  sprintf(fname, "%s/rates.hdf5", cooling->cooling_table_path);
  hid_t tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) error("unable to open file %s\n", fname);

  /* Get the group that contains the tables */
  hid_t group_id = H5Gopen(tempfile_id, "/HeatingRates", H5P_DEFAULT);
  if (group_id < 0) error("Unable to read groop %s\n", "HeatingRates");

  /* First read the length of the arrays */
  int N_bins;
  hid_t attr_id = H5Aopen(group_id, "N_bins", H5P_DEFAULT);
  herr_t status = H5Aread(attr_id, H5T_NATIVE_INT, &N_bins);
  if (status < 0) error("error reading number of bins");
  if (H5Aclose(attr_id) < 0) error("Error closing cooling attribute");
  if (N_bins != de_rijcke_cooling_N_temperatures)
    error("Invalid number of bins in file!");

  /* Allocate memory for the temperature and the cooling arrays */
  if (swift_memalign("cooling", (void **)&cooling->table.temperature,
                     SWIFT_STRUCT_ALIGNMENT, N_bins * sizeof(float)) != 0)
    error("Failed to allocate temperature table");
  if (swift_memalign("cooling", (void **)&cooling->table.cooling_rate,
                     SWIFT_STRUCT_ALIGNMENT, N_bins * sizeof(float)) != 0)
    error("Failed to allocate temperature table");

  /* Now read the datasets themselves */
  hid_t dataset =
      H5Dopen(tempfile_id, "/HeatingRates/Temperature", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.temperature);
  if (status < 0) error("error reading temperature bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing temperature dataset");

  dataset = H5Dopen(tempfile_id, "/HeatingRates/CoolingRate", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.cooling_rate);
  if (status < 0) error("error reading cooling rate bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling rate dataset");

  if (H5Gclose(group_id) < 0) error("Error closing cooling tables group");
  if (H5Fclose(tempfile_id) < 0) error("Error closing cooling tables file");
}

#endif  // SWIFTSIM_COOLING_TABLES_H
