/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "common_io.h"

/* Local includes. */
#include "engine.h"
#include "error.h"
#include "kernel_hydro.h"
#include "part.h"
#include "part_type.h"
#include "threadpool.h"
#include "tools.h"
#include "version.h"

/* I/O functions of each sub-module */
#include "black_holes_io.h"
#include "chemistry_io.h"
#include "cooling_io.h"
#include "extra_io.h"
#include "feedback.h"
#include "fof_io.h"
#include "gravity_io.h"
#include "hydro_io.h"
#include "mhd_io.h"
#include "neutrino_io.h"
#include "particle_splitting.h"
#include "rt_io.h"
#include "sink_io.h"
#include "star_formation_io.h"
#include "stars_io.h"
#include "tracers_io.h"
#include "velociraptor_io.h"

/* Some standard headers. */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_HDF5)

#include <hdf5.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Converts a C data type to the HDF5 equivalent.
 *
 * This function is a trivial wrapper around the HDF5 types but allows
 * to change the exact storage types matching the code types in a transparent
 *way.
 */
hid_t io_hdf5_type(enum IO_DATA_TYPE type) {

  switch (type) {
    case INT:
      return H5T_NATIVE_INT;
    case UINT8:
      return H5T_NATIVE_UINT8;
    case UINT:
      return H5T_NATIVE_UINT;
    case UINT64:
      return H5T_NATIVE_UINT64;
    case LONG:
      return H5T_NATIVE_LONG;
    case ULONG:
      return H5T_NATIVE_ULONG;
    case LONGLONG:
      return H5T_NATIVE_LLONG;
    case ULONGLONG:
      return H5T_NATIVE_ULLONG;
    case FLOAT:
      return H5T_NATIVE_FLOAT;
    case DOUBLE:
      return H5T_NATIVE_DOUBLE;
    case CHAR:
      return H5T_NATIVE_CHAR;
    default:
      error("Unknown type");
      return 0;
  }
}

/**
 * @brief Return 1 if the type has double precision
 *
 * Returns an error if the type is not FLOAT or DOUBLE
 */
int io_is_double_precision(enum IO_DATA_TYPE type) {

  switch (type) {
    case FLOAT:
      return 0;
    case DOUBLE:
      return 1;
    default:
      error("Invalid type");
      return 0;
  }
}

/**
 * @brief Reads an attribute (scalar) from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the attribute to read.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group.
 *
 * Calls #error() if an error occurs.
 */
void io_read_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                       void* data) {

  const hid_t h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  const hid_t h_err = H5Aread(h_attr, io_hdf5_type(type), data);
  if (h_err < 0) error("Error while reading attribute '%s'", name);

  H5Aclose(h_attr);
}

/**
 * @brief Reads an attribute from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the attribute to read.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group.
 *
 * Exits gracefully (i.e. does not read the attribute at all) if
 * it is not present, unless debugging checks are activated. If they are,
 * and the read fails, we print a warning.
 */
void io_read_attribute_graceful(hid_t grp, const char* name,
                                enum IO_DATA_TYPE type, void* data) {

  /* First, we need to check if this attribute exists to avoid raising errors
   * within the HDF5 library if we attempt to access an attribute that does
   * not exist. */
  const htri_t h_exists = H5Aexists(grp, name);

  if (h_exists <= 0) {
    /* Attribute either does not exist (0) or function failed (-ve) */
#ifdef SWIFT_DEBUG_CHECKS
    message("WARNING: attribute '%s' does not exist.", name);
#endif
  } else {
    /* Ok, now we know that it exists we can read it. */
    const hid_t h_attr = H5Aopen(grp, name, H5P_DEFAULT);

    if (h_attr >= 0) {
      const hid_t h_err = H5Aread(h_attr, io_hdf5_type(type), data);
      if (h_err < 0) {
        /* Explicitly do nothing unless debugging checks are activated */
#ifdef SWIFT_DEBUG_CHECKS
        message("WARNING: unable to read attribute '%s'", name);
#endif
      }
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      if (h_attr < 0) {
        message("WARNING: was unable to open attribute '%s'", name);
      }
#endif
    }

    H5Aclose(h_attr);
  }
}

/**
 * @brief Asserts that the redshift in the initial conditions and the one
 *        specified by the parameter file match.
 *
 * @param h_grp The Header group from the ICs
 * @param a Current scale factor as specified by parameter file
 */
void io_assert_valid_header_cosmology(hid_t h_grp, double a) {

  double redshift_from_snapshot = -1.0;
  io_read_attribute_graceful(h_grp, "Redshift", DOUBLE,
                             &redshift_from_snapshot);

  /* If the Header/Redshift value is not present, then we skip this check */
  if (redshift_from_snapshot == -1.0) return;

  const double current_redshift = 1.0 / a - 1.0;

  /* Escape if non-cosmological */
  if (current_redshift == 0.) return;

  const double redshift_fractional_difference =
      fabs(redshift_from_snapshot - current_redshift) / current_redshift;

  if (redshift_fractional_difference >= io_redshift_tolerance) {
    error(
        "Initial redshift specified in parameter file (%lf) and redshift "
        "read from initial conditions (%lf) are inconsistent.",
        current_redshift, redshift_from_snapshot);
  }
}

/**
 * @brief Reads the number of elements in a HDF5 attribute.
 *
 * @param attr The attribute from which to read.
 *
 * @return The number of elements.
 *
 * Calls #error() if an error occurs.
 */
hsize_t io_get_number_element_in_attribute(hid_t attr) {
  /* Get the dataspace */
  hid_t space = H5Aget_space(attr);
  if (space < 0) error("Failed to get data space");

  /* Read the number of dimensions */
  const int ndims = H5Sget_simple_extent_ndims(space);

  /* Read the dimensions */
  hsize_t* dims = (hsize_t*)malloc(sizeof(hsize_t) * ndims);
  H5Sget_simple_extent_dims(space, dims, NULL);

  /* Compute number of elements */
  hsize_t count = 1;
  for (int i = 0; i < ndims; i++) {
    count *= dims[i];
  }

  /* Cleanup */
  free(dims);
  H5Sclose(space);
  return count;
};

/**
 * @brief Reads an attribute (array) from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the dataset to read.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group (need to be
 * already allocated).
 * @param number_element Number of elements in the attribute.
 *
 * Calls #error() if an error occurs.
 */
void io_read_array_attribute(hid_t grp, const char* name,
                             enum IO_DATA_TYPE type, void* data,
                             hsize_t number_element) {

  /* Open attribute */
  const hid_t h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  /* Get the number of elements */
  hsize_t count = io_get_number_element_in_attribute(h_attr);

  /* Check if correct number of element */
  if (count != number_element) {
    error(
        "Error found a different number of elements than expected (%llu != "
        "%llu) in attribute %s",
        (unsigned long long)count, (unsigned long long)number_element, name);
  }

  /* Read attribute */
  const hid_t h_err = H5Aread(h_attr, io_hdf5_type(type), data);
  if (h_err < 0) error("Error while reading attribute '%s'", name);

  /* Cleanup */
  H5Aclose(h_attr);
}

/**
 * @brief Reads the number of elements in a HDF5 dataset.
 *
 * @param dataset The dataset from which to read.
 *
 * @return The number of elements.
 *
 * Calls #error() if an error occurs.
 */
hsize_t io_get_number_element_in_dataset(hid_t dataset) {
  /* Get the dataspace */
  hid_t space = H5Dget_space(dataset);
  if (space < 0) error("Failed to get data space");

  /* Read the number of dimensions */
  const int ndims = H5Sget_simple_extent_ndims(space);

  /* Read the dimensions */
  hsize_t* dims = (hsize_t*)malloc(sizeof(hsize_t) * ndims);
  H5Sget_simple_extent_dims(space, dims, NULL);

  /* Compute number of elements */
  hsize_t count = 1;
  for (int i = 0; i < ndims; i++) {
    count *= dims[i];
  }

  /* Cleanup */
  free(dims);
  H5Sclose(space);
  return count;
};

/**
 * @brief Reads a dataset (array) from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the dataset to read.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group (need to be
 * already allocated).
 * @param number_element Number of elements in the attribute.
 *
 * Calls #error() if an error occurs.
 */
void io_read_array_dataset(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                           void* data, hsize_t number_element) {

  /* Open dataset */
  const hid_t h_dataset = H5Dopen(grp, name, H5P_DEFAULT);
  if (h_dataset < 0) error("Error while opening attribute '%s'", name);

  /* Get the number of elements */
  hsize_t count = io_get_number_element_in_dataset(h_dataset);

  /* Check if correct number of element */
  if (count != number_element) {
    error(
        "Error found a different number of elements than expected (%llu != "
        "%llu) in dataset %s",
        (unsigned long long)count, (unsigned long long)number_element, name);
  }

  /* Read dataset */
  const hid_t h_err = H5Dread(h_dataset, io_hdf5_type(type), H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, data);
  if (h_err < 0) error("Error while reading dataset '%s'", name);

  /* Cleanup */
  H5Dclose(h_dataset);
}

/**
 * @brief Write an attribute to a given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param name The name of the attribute to write.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data The attribute to write.
 * @param num The number of elements to write
 *
 * Calls #error() if an error occurs.
 */
void io_write_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                        const void* data, int num) {

  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating dataspace for attribute '%s'.", name);

  hsize_t dim[1] = {(hsize_t)num};
  const hid_t h_err = H5Sset_extent_simple(h_space, 1, dim, NULL);
  if (h_err < 0)
    error("Error while changing dataspace shape for attribute '%s'.", name);

  const hid_t h_attr =
      H5Acreate1(grp, name, io_hdf5_type(type), h_space, H5P_DEFAULT);
  if (h_attr < 0) error("Error while creating attribute '%s'.", name);

  const hid_t h_err2 = H5Awrite(h_attr, io_hdf5_type(type), data);
  if (h_err2 < 0) error("Error while reading attribute '%s'.", name);

  H5Sclose(h_space);
  H5Aclose(h_attr);
}

/**
 * @brief Write a string as an attribute to a given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param name The name of the attribute to write.
 * @param str The string to write.
 * @param length The length of the string
 *
 * Calls #error() if an error occurs.
 */
void io_writeStringAttribute(hid_t grp, const char* name, const char* str,
                             int length) {

  const hid_t h_space = H5Screate(H5S_SCALAR);
  if (h_space < 0)
    error("Error while creating dataspace for attribute '%s'.", name);

  const hid_t h_type = H5Tcopy(H5T_C_S1);
  if (h_type < 0) error("Error while copying datatype 'H5T_C_S1'.");

  const hid_t h_err = H5Tset_size(h_type, length);
  if (h_err < 0) error("Error while resizing attribute type to '%i'.", length);

  const hid_t h_attr = H5Acreate1(grp, name, h_type, h_space, H5P_DEFAULT);
  if (h_attr < 0) error("Error while creating attribute '%s'.", name);

  const hid_t h_err2 = H5Awrite(h_attr, h_type, str);
  if (h_err2 < 0) error("Error while reading attribute '%s'.", name);

  H5Tclose(h_type);
  H5Sclose(h_space);
  H5Aclose(h_attr);
}

/**
 * @brief Writes a double value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void io_write_attribute_d(hid_t grp, const char* name, double data) {
  io_write_attribute(grp, name, DOUBLE, &data, 1);
}

/**
 * @brief Writes a float value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void io_write_attribute_f(hid_t grp, const char* name, float data) {
  io_write_attribute(grp, name, FLOAT, &data, 1);
}

/**
 * @brief Writes an int value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void io_write_attribute_i(hid_t grp, const char* name, int data) {
  io_write_attribute(grp, name, INT, &data, 1);
}

/**
 * @brief Writes a long value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void io_write_attribute_l(hid_t grp, const char* name, long data) {
  io_write_attribute(grp, name, LONG, &data, 1);
}

/**
 * @brief Writes a long long value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void io_write_attribute_ll(hid_t grp, const char* name, long long data) {
  io_write_attribute(grp, name, LONGLONG, &data, 1);
}

/**
 * @brief Writes a string value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param str The string to write
 */
void io_write_attribute_s(hid_t grp, const char* name, const char* str) {
  io_writeStringAttribute(grp, name, str, strlen(str));
}

/**
 * @brief Writes the meta-data of the run to an open hdf5 snapshot file.
 *
 * @param h_file The opened hdf5 file.
 * @param e The #engine containing the meta-data.
 * @param internal_units The system of units used internally.
 * @param snapshot_units The system of units used in snapshots.
 * @param fof Is this a FOF output? If so don't write subgrid info.
 */
void io_write_meta_data(hid_t h_file, const struct engine* e,
                        const struct unit_system* internal_units,
                        const struct unit_system* snapshot_units,
                        const int fof) {

  hid_t h_grp;

  /* Print the code version */
  io_write_code_description(h_file);

  /* Print the run's policy */
  io_write_engine_policy(h_file, e);

  /* Print the physical constants */
  phys_const_print_snapshot(h_file, e->physical_constants);

  if (!fof) {

    /* Print the SPH parameters */
    if (e->policy & engine_policy_hydro) {
      h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating SPH group");
      hydro_props_print_snapshot(h_grp, e->hydro_properties);
      hydro_write_flavour(h_grp);
      mhd_write_flavour(h_grp);
      H5Gclose(h_grp);
    }

    /* Print the subgrid parameters */
    h_grp = H5Gcreate(h_file, "/SubgridScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating subgrid group");
    hid_t h_grp_columns =
        H5Gcreate(h_grp, "NamedColumns", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_columns < 0) error("Error while creating named columns group");
    entropy_floor_write_flavour(h_grp);
    extra_io_write_flavour(h_grp, h_grp_columns);
    cooling_write_flavour(h_grp, h_grp_columns, e->cooling_func);
    chemistry_write_flavour(h_grp, h_grp_columns, e);
    tracers_write_flavour(h_grp);
    feedback_write_flavour(e->feedback_props, h_grp);
    rt_write_flavour(h_grp, h_grp_columns, e, internal_units, snapshot_units,
                     e->rt_props);
    H5Gclose(h_grp);

    /* Print the gravity parameters */
    if (e->policy & engine_policy_self_gravity) {
      h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating gravity group");
      gravity_props_print_snapshot(h_grp, e->gravity_properties);
      H5Gclose(h_grp);
    }

    /* Print the stellar parameters */
    if (e->policy & engine_policy_stars) {
      h_grp = H5Gcreate(h_file, "/StarsScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating stars group");
      stars_props_print_snapshot(h_grp, h_grp_columns, e->stars_properties);
      H5Gclose(h_grp);
    }

    H5Gclose(h_grp_columns);
  }

  /* Print the cosmological model  */
  h_grp =
      H5Gcreate(h_file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cosmology group");
  if (e->policy & engine_policy_cosmology)
    io_write_attribute_i(h_grp, "Cosmological run", 1);
  else
    io_write_attribute_i(h_grp, "Cosmological run", 0);
  cosmology_write_model(h_grp, e->cosmology);
  H5Gclose(h_grp);

  /* Print the runtime parameters */
  h_grp =
      H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, /*write_used=*/1);
  H5Gclose(h_grp);

  /* Print the runtime unused parameters */
  h_grp = H5Gcreate(h_file, "/UnusedParameters", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, /*write_used=*/0);
  H5Gclose(h_grp);

  /* Print the recording triggers */
  h_grp = H5Gcreate(h_file, "/RecordingTriggers", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating recording triggers group");
  if (num_snapshot_triggers_part) {
    io_write_attribute(h_grp, "DesiredRecordingTimesGas", DOUBLE,
                       e->snapshot_recording_triggers_desired_part,
                       num_snapshot_triggers_part);
    io_write_attribute(h_grp, "ActualRecordingTimesGas", DOUBLE,
                       e->snapshot_recording_triggers_part,
                       num_snapshot_triggers_part);
  }
  if (num_snapshot_triggers_spart) {
    io_write_attribute(h_grp, "DesiredRecordingTimesStars", DOUBLE,
                       e->snapshot_recording_triggers_desired_spart,
                       num_snapshot_triggers_spart);
    io_write_attribute(h_grp, "ActualRecordingTimesStars", DOUBLE,
                       e->snapshot_recording_triggers_spart,
                       num_snapshot_triggers_spart);
  }
  if (num_snapshot_triggers_bpart) {
    io_write_attribute(h_grp, "DesiredRecordingTimesBlackHoles", DOUBLE,
                       e->snapshot_recording_triggers_desired_bpart,
                       num_snapshot_triggers_bpart);
    io_write_attribute(h_grp, "ActualRecordingTimesBlackHoles", DOUBLE,
                       e->snapshot_recording_triggers_bpart,
                       num_snapshot_triggers_bpart);
  }
  H5Gclose(h_grp);

  /* Print the system of Units used in the spashot */
  io_write_unit_system(h_file, snapshot_units, "Units");

  /* Print the system of Units used internally */
  io_write_unit_system(h_file, internal_units, "InternalCodeUnits");

  /* Tell the user if a conversion will be needed */
  if (e->verbose) {
    if (units_are_equal(snapshot_units, internal_units)) {

      message("Snapshot and internal units match. No conversion needed.");

    } else {

      message("Conversion needed from:");
      message("(Snapshot) Unit system: U_M =      %e g.",
              snapshot_units->UnitMass_in_cgs);
      message("(Snapshot) Unit system: U_L =      %e cm.",
              snapshot_units->UnitLength_in_cgs);
      message("(Snapshot) Unit system: U_t =      %e s.",
              snapshot_units->UnitTime_in_cgs);
      message("(Snapshot) Unit system: U_I =      %e A.",
              snapshot_units->UnitCurrent_in_cgs);
      message("(Snapshot) Unit system: U_T =      %e K.",
              snapshot_units->UnitTemperature_in_cgs);
      message("to:");
      message("(internal) Unit system: U_M = %e g.",
              internal_units->UnitMass_in_cgs);
      message("(internal) Unit system: U_L = %e cm.",
              internal_units->UnitLength_in_cgs);
      message("(internal) Unit system: U_t = %e s.",
              internal_units->UnitTime_in_cgs);
      message("(internal) Unit system: U_I = %e A.",
              internal_units->UnitCurrent_in_cgs);
      message("(internal) Unit system: U_T = %e K.",
              internal_units->UnitTemperature_in_cgs);
    }
  }
}

/**
 * @brief Reads the Unit System from an IC file.
 *
 * If the 'Units' group does not exist in the ICs, we will use the internal
 * system of units.
 *
 * @param h_file The (opened) HDF5 file from which to read.
 * @param ic_units The unit_system to fill.
 * @param internal_units The internal system of units to copy if needed.
 * @param mpi_rank The MPI rank we are on.
 */
void io_read_unit_system(hid_t h_file, struct unit_system* ic_units,
                         const struct unit_system* internal_units,
                         int mpi_rank) {

  /* First check if it exists as this is *not* required. */
  const htri_t exists = H5Lexists(h_file, "/Units", H5P_DEFAULT);

  if (exists == 0) {

    if (mpi_rank == 0)
      message("'Units' group not found in ICs. Assuming internal unit system.");

    units_copy(ic_units, internal_units);

    return;

  } else if (exists < 0) {
    error("Serious problem with 'Units' group in ICs. H5Lexists gives %d",
          exists);
  }

  if (mpi_rank == 0) message("Reading IC units from ICs.");
  hid_t h_grp = H5Gopen(h_file, "/Units", H5P_DEFAULT);

  /* Ok, Read the damn thing */
  io_read_attribute(h_grp, "Unit length in cgs (U_L)", DOUBLE,
                    &ic_units->UnitLength_in_cgs);
  io_read_attribute(h_grp, "Unit mass in cgs (U_M)", DOUBLE,
                    &ic_units->UnitMass_in_cgs);
  io_read_attribute(h_grp, "Unit time in cgs (U_t)", DOUBLE,
                    &ic_units->UnitTime_in_cgs);
  io_read_attribute(h_grp, "Unit current in cgs (U_I)", DOUBLE,
                    &ic_units->UnitCurrent_in_cgs);
  io_read_attribute(h_grp, "Unit temperature in cgs (U_T)", DOUBLE,
                    &ic_units->UnitTemperature_in_cgs);

  /* Clean up */
  H5Gclose(h_grp);
}

/**
 * @brief Writes the current Unit System
 * @param h_file The (opened) HDF5 file in which to write
 * @param us The unit_system to dump
 * @param groupName The name of the HDF5 group to write to
 */
void io_write_unit_system(hid_t h_file, const struct unit_system* us,
                          const char* groupName) {

  const hid_t h_grpunit = H5Gcreate1(h_file, groupName, 0);
  if (h_grpunit < 0) error("Error while creating Unit System group");

  io_write_attribute_d(h_grpunit, "Unit mass in cgs (U_M)",
                       units_get_base_unit(us, UNIT_MASS));
  io_write_attribute_d(h_grpunit, "Unit length in cgs (U_L)",
                       units_get_base_unit(us, UNIT_LENGTH));
  io_write_attribute_d(h_grpunit, "Unit time in cgs (U_t)",
                       units_get_base_unit(us, UNIT_TIME));
  io_write_attribute_d(h_grpunit, "Unit current in cgs (U_I)",
                       units_get_base_unit(us, UNIT_CURRENT));
  io_write_attribute_d(h_grpunit, "Unit temperature in cgs (U_T)",
                       units_get_base_unit(us, UNIT_TEMPERATURE));

  H5Gclose(h_grpunit);
}

/**
 * @brief Writes the code version to the file
 * @param h_file The (opened) HDF5 file in which to write
 */
void io_write_code_description(hid_t h_file) {

  const hid_t h_grpcode = H5Gcreate1(h_file, "/Code", 0);
  if (h_grpcode < 0) error("Error while creating code group");

  io_write_attribute_s(h_grpcode, "Code", "SWIFT");
  io_write_attribute_s(h_grpcode, "Code Version", package_version());
  io_write_attribute_s(h_grpcode, "Compiler Name", compiler_name());
  io_write_attribute_s(h_grpcode, "Compiler Version", compiler_version());
  io_write_attribute_s(h_grpcode, "Git Branch", git_branch());
  io_write_attribute_s(h_grpcode, "Git Revision", git_revision());
  io_write_attribute_s(h_grpcode, "Git Date", git_date());
  io_write_attribute_s(h_grpcode, "Configuration options",
                       configuration_options());
  io_write_attribute_s(h_grpcode, "CFLAGS", compilation_cflags());
  io_write_attribute_s(h_grpcode, "HDF5 library version", hdf5_version());
  io_write_attribute_s(h_grpcode, "Thread barriers", thread_barrier_version());
  io_write_attribute_s(h_grpcode, "Allocators", allocator_version());
#ifdef HAVE_FFTW
  io_write_attribute_s(h_grpcode, "FFTW library version", fftw3_version());
#endif
#ifdef HAVE_LIBGSL
  io_write_attribute_s(h_grpcode, "GSL library version", libgsl_version());
#endif
#ifdef HAVE_SUNDIALS
  io_write_attribute_s(h_grpcode, "SUNDIALS library version",
                       sundials_version());
#endif
#ifdef WITH_MPI
  io_write_attribute_s(h_grpcode, "MPI library", mpi_version());
#ifdef HAVE_METIS
  io_write_attribute_s(h_grpcode, "METIS library version", metis_version());
#endif
#ifdef HAVE_PARMETIS
  io_write_attribute_s(h_grpcode, "ParMETIS library version",
                       parmetis_version());
#endif
#else
  io_write_attribute_s(h_grpcode, "MPI library", "Non-MPI version of SWIFT");
#endif
  io_write_attribute_i(h_grpcode, "RandomSeed", SWIFT_RANDOM_SEED_XOR);
  H5Gclose(h_grpcode);
}

/**
 * @brief Write the #engine policy to the file.
 * @param h_file File to write to.
 * @param e The #engine to read the policy from.
 */
void io_write_engine_policy(hid_t h_file, const struct engine* e) {

  const hid_t h_grp = H5Gcreate1(h_file, "/Policy", 0);
  if (h_grp < 0) error("Error while creating policy group");

  for (int i = 1; i < engine_maxpolicy; ++i)
    if (e->policy & (1 << i))
      io_write_attribute_i(h_grp, engine_policy_names[i + 1], 1);
    else
      io_write_attribute_i(h_grp, engine_policy_names[i + 1], 0);

  H5Gclose(h_grp);
}

void io_write_part_type_names(hid_t h_grp) {

  io_write_attribute_i(h_grp, "NumPartTypes", swift_type_count);

  /* Create an array of partcle type names */
  const int name_length = 128;
  char names[swift_type_count][name_length];
  for (int i = 0; i < swift_type_count; ++i)
    strcpy(names[i], part_type_names[i]);

  hsize_t dims[1] = {swift_type_count};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, name_length);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp, "PartTypeNames", type, space, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names[0]);
  H5Dclose(dset);
  H5Tclose(type);
  H5Sclose(space);
}

#endif /* HAVE_HDF5 */

/**
 * @brief Returns the memory size of the data type
 */
size_t io_sizeof_type(enum IO_DATA_TYPE type) {

  switch (type) {
    case INT:
      return sizeof(int);
    case UINT8:
      return sizeof(uint8_t);
    case UINT:
      return sizeof(unsigned int);
    case UINT64:
      return sizeof(uint64_t);
    case LONG:
      return sizeof(long);
    case ULONG:
      return sizeof(unsigned long);
    case LONGLONG:
      return sizeof(long long);
    case ULONGLONG:
      return sizeof(unsigned long long);
    case FLOAT:
      return sizeof(float);
    case DOUBLE:
      return sizeof(double);
    case CHAR:
      return sizeof(char);
    case SIZE_T:
      return sizeof(size_t);
    default:
      error("Unknown type");
      return 0;
  }
}

void io_prepare_dm_gparts_mapper(void* restrict data, int Ndm, void* dummy) {

  struct gpart* restrict gparts = (struct gpart*)data;

  /* Let's give all these gparts a negative id */
  for (int i = 0; i < Ndm; ++i) {

    /* Negative ids are not allowed */
    if (gparts[i].id_or_neg_offset < 0)
      error("Negative ID for DM particle %i: ID=%lld", i,
            gparts[i].id_or_neg_offset);

    /* Set gpart type */
    gparts[i].type = swift_type_dark_matter;
  }
}

/**
 * @brief Prepare the DM particles (in gparts) read in for the addition of the
 * other particle types
 *
 * This function assumes that the DM particles are all at the start of the
 * gparts array
 *
 * @param tp The current #threadpool.
 * @param gparts The array of #gpart freshly read in.
 * @param Ndm The number of DM particles read in.
 */
void io_prepare_dm_gparts(struct threadpool* tp, struct gpart* const gparts,
                          size_t Ndm) {

  threadpool_map(tp, io_prepare_dm_gparts_mapper, gparts, Ndm,
                 sizeof(struct gpart), threadpool_auto_chunk_size, NULL);
}

void io_prepare_dm_background_gparts_mapper(void* restrict data, int Ndm,
                                            void* dummy) {

  struct gpart* restrict gparts = (struct gpart*)data;

  /* Let's give all these gparts a negative id */
  for (int i = 0; i < Ndm; ++i) {

    /* Negative ids are not allowed */
    if (gparts[i].id_or_neg_offset < 0)
      error("Negative ID for DM particle %i: ID=%lld", i,
            gparts[i].id_or_neg_offset);

    /* Set gpart type */
    gparts[i].type = swift_type_dark_matter_background;
  }
}

/**
 * @brief Prepare the DM backgorund particles (in gparts) read in
 * for the addition of the other particle types
 *
 * This function assumes that the DM particles are all at the start of the
 * gparts array and that the background particles directly follow them.
 *
 * @param tp The current #threadpool.
 * @param gparts The array of #gpart freshly read in.
 * @param Ndm The number of DM particles read in.
 */
void io_prepare_dm_background_gparts(struct threadpool* tp,
                                     struct gpart* const gparts, size_t Ndm) {

  threadpool_map(tp, io_prepare_dm_background_gparts_mapper, gparts, Ndm,
                 sizeof(struct gpart), threadpool_auto_chunk_size, NULL);
}

void io_prepare_dm_neutrino_gparts_mapper(void* restrict data, int Ndm,
                                          void* dummy) {

  struct gpart* restrict gparts = (struct gpart*)data;

  /* Let's give all these gparts a negative id */
  for (int i = 0; i < Ndm; ++i) {

    /* Negative ids are not allowed */
    if (gparts[i].id_or_neg_offset < 0)
      error("Negative ID for DM particle %i: ID=%lld", i,
            gparts[i].id_or_neg_offset);

    /* Set gpart type */
    gparts[i].type = swift_type_neutrino;
  }
}

/**
 * @brief Prepare the neutrino dark matter particles (in gparts) read in
 * for the addition of the other particle types
 *
 * This function assumes that the DM & background DM particles are all at the
 * start of the gparts array and that the neutrinos directly follow them.
 *
 * @param tp The current #threadpool.
 * @param gparts The array of #gpart freshly read in.
 * @param Ndm The number of DM particles read in.
 */
void io_prepare_dm_neutrino_gparts(struct threadpool* tp,
                                   struct gpart* const gparts, size_t Ndm) {

  threadpool_map(tp, io_prepare_dm_neutrino_gparts_mapper, gparts, Ndm,
                 sizeof(struct gpart), threadpool_auto_chunk_size, NULL);
}

struct duplication_data {

  struct part* parts;
  struct gpart* gparts;
  struct spart* sparts;
  struct bpart* bparts;
  struct sink* sinks;
  int Ndm;
  int Ngas;
  int Nstars;
  int Nblackholes;
  int Nsinks;
};

void io_duplicate_hydro_gparts_mapper(void* restrict data, int Ngas,
                                      void* restrict extra_data) {

  struct duplication_data* temp = (struct duplication_data*)extra_data;
  const int Ndm = temp->Ndm;
  struct part* parts = (struct part*)data;
  const ptrdiff_t offset = parts - temp->parts;
  struct gpart* gparts = temp->gparts + offset;

  for (int i = 0; i < Ngas; ++i) {

    /* Duplicate the crucial information */
    gparts[i + Ndm].x[0] = parts[i].x[0];
    gparts[i + Ndm].x[1] = parts[i].x[1];
    gparts[i + Ndm].x[2] = parts[i].x[2];

    gparts[i + Ndm].v_full[0] = parts[i].v[0];
    gparts[i + Ndm].v_full[1] = parts[i].v[1];
    gparts[i + Ndm].v_full[2] = parts[i].v[2];

    gparts[i + Ndm].mass = hydro_get_mass(&parts[i]);

    /* Set gpart type */
    gparts[i + Ndm].type = swift_type_gas;

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -(long long)(offset + i);
    parts[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every #part into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles are all at the start of the
 * gparts array and adds the hydro particles afterwards
 *
 * @param tp The current #threadpool.
 * @param parts The array of #part freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM particles
 * at the start
 * @param Ngas The number of gas particles read in.
 * @param Ndm The number of DM particles read in.
 */
void io_duplicate_hydro_gparts(struct threadpool* tp, struct part* const parts,
                               struct gpart* const gparts, size_t Ngas,
                               size_t Ndm) {
  struct duplication_data data;
  data.parts = parts;
  data.gparts = gparts;
  data.Ndm = Ndm;

  threadpool_map(tp, io_duplicate_hydro_gparts_mapper, parts, Ngas,
                 sizeof(struct part), threadpool_auto_chunk_size, &data);
}

void io_duplicate_stars_gparts_mapper(void* restrict data, int Nstars,
                                      void* restrict extra_data) {

  struct duplication_data* temp = (struct duplication_data*)extra_data;
  const int Ndm = temp->Ndm;
  struct spart* sparts = (struct spart*)data;
  const ptrdiff_t offset = sparts - temp->sparts;
  struct gpart* gparts = temp->gparts + offset;

  for (int i = 0; i < Nstars; ++i) {

    /* Duplicate the crucial information */
    gparts[i + Ndm].x[0] = sparts[i].x[0];
    gparts[i + Ndm].x[1] = sparts[i].x[1];
    gparts[i + Ndm].x[2] = sparts[i].x[2];

    gparts[i + Ndm].v_full[0] = sparts[i].v[0];
    gparts[i + Ndm].v_full[1] = sparts[i].v[1];
    gparts[i + Ndm].v_full[2] = sparts[i].v[2];

    gparts[i + Ndm].mass = sparts[i].mass;

    /* Set gpart type */
    gparts[i + Ndm].type = swift_type_stars;

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -(long long)(offset + i);
    sparts[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every #spart into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles and gas particles are all at
 * the start of the gparts array and adds the star particles afterwards
 *
 * @param tp The current #threadpool.
 * @param sparts The array of #spart freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM and gas
 * particles at the start.
 * @param Nstars The number of stars particles read in.
 * @param Ndm The number of DM and gas particles read in.
 */
void io_duplicate_stars_gparts(struct threadpool* tp,
                               struct spart* const sparts,
                               struct gpart* const gparts, size_t Nstars,
                               size_t Ndm) {

  struct duplication_data data;
  data.gparts = gparts;
  data.sparts = sparts;
  data.Ndm = Ndm;

  threadpool_map(tp, io_duplicate_stars_gparts_mapper, sparts, Nstars,
                 sizeof(struct spart), threadpool_auto_chunk_size, &data);
}

void io_duplicate_sinks_gparts_mapper(void* restrict data, int Nsinks,
                                      void* restrict extra_data) {

  struct duplication_data* temp = (struct duplication_data*)extra_data;
  const int Ndm = temp->Ndm;
  struct sink* sinks = (struct sink*)data;
  const ptrdiff_t offset = sinks - temp->sinks;
  struct gpart* gparts = temp->gparts + offset;

  for (int i = 0; i < Nsinks; ++i) {

    /* Duplicate the crucial information */
    gparts[i + Ndm].x[0] = sinks[i].x[0];
    gparts[i + Ndm].x[1] = sinks[i].x[1];
    gparts[i + Ndm].x[2] = sinks[i].x[2];

    gparts[i + Ndm].v_full[0] = sinks[i].v[0];
    gparts[i + Ndm].v_full[1] = sinks[i].v[1];
    gparts[i + Ndm].v_full[2] = sinks[i].v[2];

    gparts[i + Ndm].mass = sinks[i].mass;

    /* Set gpart type */
    gparts[i + Ndm].type = swift_type_sink;

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -(long long)(offset + i);
    sinks[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every #sink into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles, gas particles and star particles
 * are all at the start of the gparts array and adds the sinks particles
 * afterwards
 *
 * @param tp The current #threadpool.
 * @param sinks The array of #sink freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM, gas
 * and star particles at the start.
 * @param Nsinks The number of sink particles read in.
 * @param Ndm The number of DM, gas and star particles read in.
 */
void io_duplicate_sinks_gparts(struct threadpool* tp, struct sink* const sinks,
                               struct gpart* const gparts, size_t Nsinks,
                               size_t Ndm) {

  struct duplication_data data;
  data.gparts = gparts;
  data.sinks = sinks;
  data.Ndm = Ndm;

  threadpool_map(tp, io_duplicate_sinks_gparts_mapper, sinks, Nsinks,
                 sizeof(struct sink), threadpool_auto_chunk_size, &data);
}

void io_duplicate_black_holes_gparts_mapper(void* restrict data,
                                            int Nblackholes,
                                            void* restrict extra_data) {

  struct duplication_data* temp = (struct duplication_data*)extra_data;
  const int Ndm = temp->Ndm;
  struct bpart* bparts = (struct bpart*)data;
  const ptrdiff_t offset = bparts - temp->bparts;
  struct gpart* gparts = temp->gparts + offset;

  for (int i = 0; i < Nblackholes; ++i) {

    /* Duplicate the crucial information */
    gparts[i + Ndm].x[0] = bparts[i].x[0];
    gparts[i + Ndm].x[1] = bparts[i].x[1];
    gparts[i + Ndm].x[2] = bparts[i].x[2];

    gparts[i + Ndm].v_full[0] = bparts[i].v[0];
    gparts[i + Ndm].v_full[1] = bparts[i].v[1];
    gparts[i + Ndm].v_full[2] = bparts[i].v[2];

    gparts[i + Ndm].mass = bparts[i].mass;

    /* Set gpart type */
    gparts[i + Ndm].type = swift_type_black_hole;

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -(long long)(offset + i);
    bparts[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every #bpart into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles, gas particles and star particles
 * are all at the start of the gparts array and adds the black hole particles
 * afterwards
 *
 * @param tp The current #threadpool.
 * @param bparts The array of #bpart freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM, gas
 * and star particles at the start.
 * @param Nblackholes The number of blackholes particles read in.
 * @param Ndm The number of DM, gas and star particles read in.
 */
void io_duplicate_black_holes_gparts(struct threadpool* tp,
                                     struct bpart* const bparts,
                                     struct gpart* const gparts,
                                     size_t Nblackholes, size_t Ndm) {

  struct duplication_data data;
  data.gparts = gparts;
  data.bparts = bparts;
  data.Ndm = Ndm;

  threadpool_map(tp, io_duplicate_black_holes_gparts_mapper, bparts,
                 Nblackholes, sizeof(struct bpart), threadpool_auto_chunk_size,
                 &data);
}

/**
 * @brief Copy every non-inhibited #part into the parts_written array.
 *
 * Also takes into account possible downsampling.
 *
 * @param parts The array of #part containing all particles.
 * @param xparts The array of #xpart containing all particles.
 * @param parts_written The array of #part to fill with particles we want to
 * write.
 * @param xparts_written The array of #xpart  to fill with particles we want to
 * write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Nparts The total number of #part.
 * @param Nparts_written The total number of #part to write.
 */
void io_collect_parts_to_write(const struct part* restrict parts,
                               const struct xpart* restrict xparts,
                               struct part* restrict parts_written,
                               struct xpart* restrict xparts_written,
                               const int subsample, const float subsample_ratio,
                               const int snap_num, const size_t Nparts,
                               const size_t Nparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Nparts; ++i) {

    /* And collect the ones that have not been removed */
    if (parts[i].time_bin != time_bin_inhibited &&
        parts[i].time_bin != time_bin_not_created) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r = random_unit_interval(parts[i].id, snap_num,
                                             random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      parts_written[count] = parts[i];
      xparts_written[count] = xparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Nparts_written)
    error("Collected the wrong number of particles (%zu vs. %zu expected)",
          count, Nparts_written);
}

/**
 * @brief Copy every non-inhibited #spart into the sparts_written array.
 *
 * Also takes into account possible downsampling.
 *
 * @param sparts The array of #spart containing all particles.
 * @param sparts_written The array of #spart to fill with particles we want to
 * write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Nsparts The total number of #part.
 * @param Nsparts_written The total number of #part to write.
 */
void io_collect_sparts_to_write(const struct spart* restrict sparts,
                                struct spart* restrict sparts_written,
                                const int subsample,
                                const float subsample_ratio, const int snap_num,
                                const size_t Nsparts,
                                const size_t Nsparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Nsparts; ++i) {

    /* And collect the ones that have not been removed */
    if (sparts[i].time_bin != time_bin_inhibited &&
        sparts[i].time_bin != time_bin_not_created) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r = random_unit_interval(sparts[i].id, snap_num,
                                             random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      sparts_written[count] = sparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Nsparts_written)
    error("Collected the wrong number of s-particles (%zu vs. %zu expected)",
          count, Nsparts_written);
}

/**
 * @brief Copy every non-inhibited #sink into the sinks_written array.
 *
 * Also takes into account possible downsampling.
 *
 * @param sinks The array of #sink containing all particles.
 * @param sinks_written The array of #sink to fill with particles we want to
 * write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Nsinks The total number of #sink.
 * @param Nsinks_written The total number of #sink to write.
 */
void io_collect_sinks_to_write(const struct sink* restrict sinks,
                               struct sink* restrict sinks_written,
                               const int subsample, const float subsample_ratio,
                               const int snap_num, const size_t Nsinks,
                               const size_t Nsinks_written) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Nsinks; ++i) {

    /* And collect the ones that have not been removed */
    if (sinks[i].time_bin != time_bin_inhibited &&
        sinks[i].time_bin != time_bin_not_created) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r = random_unit_interval(sinks[i].id, snap_num,
                                             random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      sinks_written[count] = sinks[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Nsinks_written)
    error("Collected the wrong number of sink-particles (%zu vs. %zu expected)",
          count, Nsinks_written);
}

/**
 * @brief Copy every non-inhibited #bpart into the bparts_written array.
 *
 * Also takes into account possible downsampling.
 *
 * @param bparts The array of #bpart containing all particles.
 * @param bparts_written The array of #bpart to fill with particles we want to
 * write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Nbparts The total number of #part.
 * @param Nbparts_written The total number of #part to write.
 */
void io_collect_bparts_to_write(const struct bpart* restrict bparts,
                                struct bpart* restrict bparts_written,
                                const int subsample,
                                const float subsample_ratio, const int snap_num,
                                const size_t Nbparts,
                                const size_t Nbparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Nbparts; ++i) {

    /* And collect the ones that have not been removed */
    if (bparts[i].time_bin != time_bin_inhibited &&
        bparts[i].time_bin != time_bin_not_created) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r = random_unit_interval(bparts[i].id, snap_num,
                                             random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      bparts_written[count] = bparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Nbparts_written)
    error("Collected the wrong number of s-particles (%zu vs. %zu expected)",
          count, Nbparts_written);
}

/**
 * @brief Copy every non-inhibited regulat DM #gpart into the gparts_written
 * array.
 *
 * Also takes into account possible downsampling.
 *
 * @param gparts The array of #gpart containing all particles.
 * @param vr_data The array of gpart-related VELOCIraptor output.
 * @param gparts_written The array of #gpart to fill with particles we want to
 * write.
 * @param vr_data_written The array of gpart-related VELOCIraptor with particles
 * we want to write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Ngparts The total number of #part.
 * @param Ngparts_written The total number of #part to write.
 * @param with_stf Are we running with STF? i.e. do we want to collect vr data?
 */
void io_collect_gparts_to_write(
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* restrict vr_data,
    struct gpart* restrict gparts_written,
    struct velociraptor_gpart_data* restrict vr_data_written,
    const int subsample, const float subsample_ratio, const int snap_num,
    const size_t Ngparts, const size_t Ngparts_written, const int with_stf) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Ngparts; ++i) {

    /* And collect the ones that have not been removed */
    if ((gparts[i].time_bin != time_bin_inhibited) &&
        (gparts[i].time_bin != time_bin_not_created) &&
        (gparts[i].type == swift_type_dark_matter)) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r =
            random_unit_interval(gparts[i].id_or_neg_offset, snap_num,
                                 random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      if (with_stf) vr_data_written[count] = vr_data[i];

      gparts_written[count] = gparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Ngparts_written)
    error("Collected the wrong number of g-particles (%zu vs. %zu expected)",
          count, Ngparts_written);
}

/**
 * @brief Copy every non-inhibited background DM #gpart into the gparts_written
 * array.
 *
 * Also takes into account possible downsampling.
 *
 * @param gparts The array of #gpart containing all particles.
 * @param vr_data The array of gpart-related VELOCIraptor output.
 * @param gparts_written The array of #gpart to fill with particles we want to
 * write.
 * @param vr_data_written The array of gpart-related VELOCIraptor with particles
 * we want to write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Ngparts The total number of #part.
 * @param Ngparts_written The total number of #part to write.
 * @param with_stf Are we running with STF? i.e. do we want to collect vr data?
 */
void io_collect_gparts_background_to_write(
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* restrict vr_data,
    struct gpart* restrict gparts_written,
    struct velociraptor_gpart_data* restrict vr_data_written,
    const int subsample, const float subsample_ratio, const int snap_num,
    const size_t Ngparts, const size_t Ngparts_written, const int with_stf) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Ngparts; ++i) {

    /* And collect the ones that have not been removed */
    if ((gparts[i].time_bin != time_bin_inhibited) &&
        (gparts[i].time_bin != time_bin_not_created) &&
        (gparts[i].type == swift_type_dark_matter_background)) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r =
            random_unit_interval(gparts[i].id_or_neg_offset, snap_num,
                                 random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      if (with_stf) vr_data_written[count] = vr_data[i];

      gparts_written[count] = gparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Ngparts_written)
    error("Collected the wrong number of g-particles (%zu vs. %zu expected)",
          count, Ngparts_written);
}

/**
 * @brief Copy every non-inhibited neutrino DM #gpart into the gparts_written
 * array.
 *
 * Also takes into account possible downsampling.
 *
 * @param gparts The array of #gpart containing all particles.
 * @param vr_data The array of gpart-related VELOCIraptor output.
 * @param gparts_written The array of #gpart to fill with particles we want to
 * write.
 * @param vr_data_written The array of gpart-related VELOCIraptor with particles
 * we want to write.
 * @param subsample Are we subsampling the particles?
 * @param subsample_ratio The fraction of particles to write if subsampling.
 * @param snap_num The snapshot ID (used to seed the RNG when sub-sampling).
 * @param Ngparts The total number of #part.
 * @param Ngparts_written The total number of #part to write.
 * @param with_stf Are we running with STF? i.e. do we want to collect vr data?
 */
void io_collect_gparts_neutrino_to_write(
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* restrict vr_data,
    struct gpart* restrict gparts_written,
    struct velociraptor_gpart_data* restrict vr_data_written,
    const int subsample, const float subsample_ratio, const int snap_num,
    const size_t Ngparts, const size_t Ngparts_written, const int with_stf) {

  size_t count = 0;

  /* Loop over all parts */
  for (size_t i = 0; i < Ngparts; ++i) {

    /* And collect the ones that have not been removed */
    if ((gparts[i].time_bin != time_bin_inhibited) &&
        (gparts[i].time_bin != time_bin_not_created) &&
        (gparts[i].type == swift_type_neutrino)) {

      /* When subsampling, select particles at random */
      if (subsample) {
        const float r =
            random_unit_interval(gparts[i].id_or_neg_offset, snap_num,
                                 random_number_snapshot_sampling);

        if (r > subsample_ratio) continue;
      }

      if (with_stf) vr_data_written[count] = vr_data[i];

      gparts_written[count] = gparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Ngparts_written)
    error("Collected the wrong number of g-particles (%zu vs. %zu expected)",
          count, Ngparts_written);
}

/**
 * @brief Create the subdirectory for snapshots if the user demanded one.
 *
 * Does nothing if the directory is '.'
 *
 * @param dirname The name of the directory.
 */
void io_make_snapshot_subdir(const char* dirname) {

  if (strcmp(dirname, ".") != 0 && strnlen(dirname, PARSER_MAX_LINE_SIZE) > 0) {
    safe_checkdir(dirname, /*create=*/1);
  }
}

/**
 * @brief Construct the file names for a single-file hdf5 snapshots and
 * corresponding XMF descriptor file.
 *
 * The XMF file always uses the default basename.
 *
 * @param filename (return) The file name of the hdf5 snapshot.
 * @param xmf_filename (return) The file name of the associated XMF file.
 * @param use_time_label Are we using time labels for the snapshot indices?
 * @param snapshots_invoke_stf Are we calling STF when dumping a snapshot?
 * @param time The current simulation time.
 * @param stf_count The counter of STF outputs.
 * @param snap_count The counter of snapshot outputs.
 * @param default_subdir The common part of the default sub-directory names.
 * @param subdir The sub-directory in which the snapshots are written.
 * @param default_basename The common part of the default snapshot names.
 * @param basename The common part of the snapshot names.
 */
void io_get_snapshot_filename(char filename[FILENAME_BUFFER_SIZE],
                              char xmf_filename[FILENAME_BUFFER_SIZE],
                              const struct output_list* output_list,
                              const int snapshots_invoke_stf,
                              const int stf_count, const int snap_count,
                              const char* default_subdir, const char* subdir,
                              const char* default_basename,
                              const char* basename) {

  int snap_number = -1;
  int number_digits = -1;
  if (output_list && output_list->alternative_labels_on) {
    snap_number = output_list->snapshot_labels[snap_count];
    number_digits = 0;
  } else if (snapshots_invoke_stf) {
    snap_number = stf_count;
    number_digits = 4;
  } else {
    snap_number = snap_count;
    number_digits = 4;
  }

  /* Are we using a sub-dir? */
  if (strlen(subdir) > 0) {
    sprintf(filename, "%s/%s_%0*d.hdf5", subdir, basename, number_digits,
            snap_number);
    sprintf(xmf_filename, "%s/%s.xmf", default_subdir, default_basename);
  } else {
    sprintf(filename, "%s_%0*d.hdf5", basename, number_digits, snap_number);
    sprintf(xmf_filename, "%s.xmf", default_basename);
  }
}
/**
 * @brief Set all ParticleIDs for each gpart to 1.
 *
 * Function is called when remap_ids is 1.
 *
 * Note only the gparts IDs have to be set to 1, as other parttypes can survive
 * as ParticleIDs=0 until the remapping routine.
 *
 * @param gparts The array of loaded gparts.
 * @param Ngparts Number of loaded gparts.
 */
void io_set_ids_to_one(struct gpart* gparts, const size_t Ngparts) {
  for (size_t i = 0; i < Ngparts; i++) gparts[i].id_or_neg_offset = 1;
}

/**
 * @brief Select the fields to write to snapshots for the gas particles.
 *
 * @param parts The #part's
 * @param xparts The #xpart's
 * @param with_cosmology Are we running with cosmology switched on?
 * @param with_cooling Are we running with cooling switched on?
 * @param with_temperature Are we running with temperature switched on?
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param with_rt Are we running with radiative transfer?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_hydro_fields(const struct part* const parts,
                            const struct xpart* const xparts,
                            const int with_cosmology, const int with_cooling,
                            const int with_temperature, const int with_fof,
                            const int with_stf, const int with_rt,
                            const struct engine* const e, int* const num_fields,
                            struct io_props* const list) {

  hydro_write_particles(parts, xparts, list, num_fields);

  *num_fields += mhd_write_particles(parts, xparts, list + *num_fields);

  *num_fields += particle_splitting_write_particles(
      parts, xparts, list + *num_fields, with_cosmology);
  *num_fields += chemistry_write_particles(parts, xparts, list + *num_fields,
                                           with_cosmology);
  if (with_cooling || with_temperature) {
    *num_fields += cooling_write_particles(parts, xparts, list + *num_fields);
  }
  if (with_fof) {
    *num_fields += fof_write_parts(parts, xparts, list + *num_fields);
  }
  if (with_stf) {
    *num_fields += velociraptor_write_parts(parts, xparts, list + *num_fields);
  }
  *num_fields += tracers_write_particles(parts, xparts, list + *num_fields,
                                         with_cosmology);
  *num_fields +=
      star_formation_write_particles(parts, xparts, list + *num_fields);
  if (with_rt) {
    *num_fields += rt_write_particles(parts, list + *num_fields);
  }
  *num_fields += extra_io_write_particles(parts, xparts, list + *num_fields,
                                          with_cosmology);
}

/**
 * @brief Select the fields to write to snapshots for the DM particles.
 *
 * @param gparts The #gpart's
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_dm_fields(const struct gpart* const gparts,
                         const struct velociraptor_gpart_data* gpart_group_data,
                         const int with_fof, const int with_stf,
                         const struct engine* const e, int* const num_fields,
                         struct io_props* const list) {

  darkmatter_write_particles(gparts, list, num_fields);
  if (with_fof) {
    *num_fields += fof_write_gparts(gparts, list + *num_fields);
  }
  if (with_stf) {
    *num_fields +=
        velociraptor_write_gparts(gpart_group_data, list + *num_fields);
  }
}

/**
 * @brief Select the fields to write to snapshots for the neutrino particles.
 *
 * @param gparts The #gpart's
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_neutrino_fields(
    const struct gpart* const gparts,
    const struct velociraptor_gpart_data* gpart_group_data, const int with_fof,
    const int with_stf, const struct engine* const e, int* const num_fields,
    struct io_props* const list) {

  darkmatter_write_particles(gparts, list, num_fields);

  *num_fields += neutrino_write_particles(gparts, list + *num_fields);
}

/**
 * @brief Select the fields to write to snapshots for the sink particles.
 *
 * @param sinks The #sink's
 * @param with_cosmology Are we running with cosmology switched on?
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_sink_fields(const struct sink* const sinks,
                           const int with_cosmology, const int with_fof,
                           const int with_stf, const struct engine* const e,
                           int* const num_fields, struct io_props* const list) {

  sink_write_particles(sinks, list, num_fields, with_cosmology);
}

/**
 * @brief Select the fields to write to snapshots for the star particles.
 *
 * @param sparts The #spart's
 * @param with_cosmology Are we running with cosmology switched on?
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param with_rt Are we running with radiative transfer?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_star_fields(const struct spart* const sparts,
                           const int with_cosmology, const int with_fof,
                           const int with_stf, const int with_rt,
                           const struct engine* const e, int* const num_fields,
                           struct io_props* const list) {

  stars_write_particles(sparts, list, num_fields, with_cosmology);
  *num_fields +=
      particle_splitting_write_sparticles(sparts, list + *num_fields);
  *num_fields += chemistry_write_sparticles(sparts, list + *num_fields);
  *num_fields +=
      tracers_write_sparticles(sparts, list + *num_fields, with_cosmology);
  *num_fields += star_formation_write_sparticles(sparts, list + *num_fields);
  if (with_fof) {
    *num_fields += fof_write_sparts(sparts, list + *num_fields);
  }
  if (with_stf) {
    *num_fields += velociraptor_write_sparts(sparts, list + *num_fields);
  }
  if (with_rt) {
    *num_fields += rt_write_stars(sparts, list + *num_fields);
  }
  *num_fields +=
      extra_io_write_sparticles(sparts, list + *num_fields, with_cosmology);
}

/**
 * @brief Select the fields to write to snapshots for the BH particles.
 *
 * @param bparts The #bpart's
 * @param with_cosmology Are we running with cosmology switched on?
 * @param with_fof Are we running FoF?
 * @param with_stf Are we running with structure finding?
 * @param e The #engine (to access scheme properties).
 * @param num_fields (return) The number of fields to write.
 * @param list (return) The list of fields to write.
 */
void io_select_bh_fields(const struct bpart* const bparts,
                         const int with_cosmology, const int with_fof,
                         const int with_stf, const struct engine* const e,
                         int* const num_fields, struct io_props* const list) {

  black_holes_write_particles(bparts, list, num_fields, with_cosmology);
  *num_fields +=
      particle_splitting_write_bparticles(bparts, list + *num_fields);
  *num_fields += chemistry_write_bparticles(bparts, list + *num_fields);
  *num_fields +=
      tracers_write_bparticles(bparts, list + *num_fields, with_cosmology);
  if (with_fof) {
    *num_fields += fof_write_bparts(bparts, list + *num_fields);
  }
  if (with_stf) {
    *num_fields += velociraptor_write_bparts(bparts, list + *num_fields);
  }
  *num_fields +=
      extra_io_write_bparticles(bparts, list + *num_fields, with_cosmology);
}
