/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

#if defined(HAVE_HDF5)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "common_io.h"

/* Local includes. */
#include "error.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "part.h"
#include "units.h"
#include "version.h"

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
    case UINT:
      return H5T_NATIVE_UINT;
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
      return H5T_C_S1;
    default:
      error("Unknown type");
      return 0;
  }
}

/**
 * @brief Returns the memory size of the data type
 */
size_t io_sizeof_type(enum IO_DATA_TYPE type) {

  switch (type) {
    case INT:
      return sizeof(int);
    case UINT:
      return sizeof(unsigned int);
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
 * @brief Reads an attribute from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the attribute to read.
 * @param type The #IO_DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group.
 *
 * Calls #error() if an error occurs.
 */
void io_read_attribute(hid_t grp, char* name, enum IO_DATA_TYPE type,
                       void* data) {
  hid_t h_attr = 0, h_err = 0;

  h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if (h_attr < 0) {
    error("Error while opening attribute '%s'", name);
  }

  h_err = H5Aread(h_attr, io_hdf5_type(type), data);
  if (h_err < 0) {
    error("Error while reading attribute '%s'", name);
  }

  H5Aclose(h_attr);
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
                        void* data, int num) {
  hid_t h_space = 0, h_attr = 0, h_err = 0;
  hsize_t dim[1] = {num};

  h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0) {
    error("Error while creating dataspace for attribute '%s'.", name);
  }

  h_err = H5Sset_extent_simple(h_space, 1, dim, NULL);
  if (h_err < 0) {
    error("Error while changing dataspace shape for attribute '%s'.", name);
  }

  h_attr = H5Acreate1(grp, name, io_hdf5_type(type), h_space, H5P_DEFAULT);
  if (h_attr < 0) {
    error("Error while creating attribute '%s'.", name);
  }

  h_err = H5Awrite(h_attr, io_hdf5_type(type), data);
  if (h_err < 0) {
    error("Error while reading attribute '%s'.", name);
  }

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
  hid_t h_space = 0, h_attr = 0, h_err = 0, h_type = 0;

  h_space = H5Screate(H5S_SCALAR);
  if (h_space < 0) {
    error("Error while creating dataspace for attribute '%s'.", name);
  }

  h_type = H5Tcopy(H5T_C_S1);
  if (h_type < 0) {
    error("Error while copying datatype 'H5T_C_S1'.");
  }

  h_err = H5Tset_size(h_type, length);
  if (h_err < 0) {
    error("Error while resizing attribute type to '%i'.", length);
  }

  h_attr = H5Acreate1(grp, name, h_type, h_space, H5P_DEFAULT);
  if (h_attr < 0) {
    error("Error while creating attribute '%s'.", name);
  }

  h_err = H5Awrite(h_attr, h_type, str);
  if (h_err < 0) {
    error("Error while reading attribute '%s'.", name);
  }

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
 * @brief Writes a string value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param str The string to write
 */
void io_write_attribute_s(hid_t grp, const char* name, const char* str) {
  io_writeStringAttribute(grp, name, str, strlen(str));
}

/**
 * @brief Reads the Unit System from an IC file.
 * @param h_file The (opened) HDF5 file from which to read.
 * @param us The UnitSystem to fill.
 *
 * If the 'Units' group does not exist in the ICs, cgs units will be assumed
 */
void io_read_UnitSystem(hid_t h_file, struct UnitSystem* us) {

  hid_t h_grp = H5Gopen(h_file, "/Units", H5P_DEFAULT);

  if (h_grp < 0) {
    message("'Units' group not found in ICs. Assuming CGS unit system.");

    /* Default to CGS */
    us->UnitMass_in_cgs = 1.;
    us->UnitLength_in_cgs = 1.;
    us->UnitTime_in_cgs = 1.;
    us->UnitCurrent_in_cgs = 1.;
    us->UnitTemperature_in_cgs = 1.;

    return;
  }

  /* Ok, Read the damn thing */
  io_read_attribute(h_grp, "Unit length in cgs (U_L)", DOUBLE,
                    &us->UnitLength_in_cgs);
  io_read_attribute(h_grp, "Unit mass in cgs (U_M)", DOUBLE,
                    &us->UnitMass_in_cgs);
  io_read_attribute(h_grp, "Unit time in cgs (U_t)", DOUBLE,
                    &us->UnitTime_in_cgs);
  io_read_attribute(h_grp, "Unit current in cgs (U_I)", DOUBLE,
                    &us->UnitCurrent_in_cgs);
  io_read_attribute(h_grp, "Unit temperature in cgs (U_T)", DOUBLE,
                    &us->UnitTemperature_in_cgs);

  /* Clean up */
  H5Gclose(h_grp);
}

/**
 * @brief Writes the current Unit System
 * @param h_file The (opened) HDF5 file in which to write
 * @param us The UnitSystem to dump
 * @param groupName The name of the HDF5 group to write to
 */
void io_write_UnitSystem(hid_t h_file, const struct UnitSystem* us,
                         const char* groupName) {

  hid_t h_grpunit = 0;
  h_grpunit = H5Gcreate1(h_file, groupName, 0);
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
  hid_t h_grpcode = 0;

  h_grpcode = H5Gcreate1(h_file, "/Code", 0);
  if (h_grpcode < 0) error("Error while creating code group");

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
#ifdef HAVE_FFTW
  io_write_attribute_s(h_grpcode, "FFTW library version", fftw3_version());
#endif
#ifdef WITH_MPI
  io_write_attribute_s(h_grpcode, "MPI library", mpi_version());
#ifdef HAVE_METIS
  io_write_attribute_s(h_grpcode, "METIS library version", metis_version());
#endif
#else
  io_write_attribute_s(h_grpcode, "MPI library", "Non-MPI version of SWIFT");
#endif
  H5Gclose(h_grpcode);
}

#endif /* HAVE_HDF5 */

/* ------------------------------------------------------------------------------------------------
 * This part writes the XMF file descriptor enabling a visualisation through
 * ParaView
 * ------------------------------------------------------------------------------------------------
 */

/**
 * @brief Prepares the XMF file for the new entry
 *
 * Creates a temporary file on the disk in order to copy the right lines.
 *
 * @todo Use a proper XML library to avoid stupid copies.
 */

/**
 * @brief Prepare the DM particles (in gparts) read in for the addition of the
 * other particle types
 *
 * This function assumes that the DM particles are all at the start of the
 * gparts array
 *
 * @param gparts The array of #gpart freshly read in.
 * @param Ndm The number of DM particles read in.
 */
void io_prepare_dm_gparts(struct gpart* const gparts, size_t Ndm) {

  /* Let's give all these gparts a negative id */
  for (size_t i = 0; i < Ndm; ++i) {
    /* 0 or negative ids are not allowed */
    if (gparts[i].id_or_neg_offset <= 0)
      error("0 or negative ID for DM particle %zu: ID=%lld", i,
            gparts[i].id_or_neg_offset);

    /* Set gpart type */
    gparts[i].type = swift_type_dark_matter;
  }
}

/**
 * @brief Copy every #part into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles are all at the start of the
 * gparts array and adds the hydro particles afterwards
 *
 * @param parts The array of #part freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM particles
 * at the start
 * @param Ngas The number of gas particles read in.
 * @param Ndm The number of DM particles read in.
 */
void io_duplicate_hydro_gparts(struct part* const parts,
                               struct gpart* const gparts, size_t Ngas,
                               size_t Ndm) {

  for (size_t i = 0; i < Ngas; ++i) {

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
    gparts[i + Ndm].id_or_neg_offset = -i;
    parts[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every #spart into the corresponding #gpart and link them.
 *
 * This function assumes that the DM particles and gas particles are all at
 * the start of the gparts array and adds the star particles afterwards
 *
 * @param sparts The array of #spart freshly read in.
 * @param gparts The array of #gpart freshly read in with all the DM and gas
 * particles at the start.
 * @param Nstars The number of stars particles read in.
 * @param Ndm The number of DM and gas particles read in.
 */
void io_duplicate_star_gparts(struct spart* const sparts,
                              struct gpart* const gparts, size_t Nstars,
                              size_t Ndm) {

  for (size_t i = 0; i < Nstars; ++i) {

    /* Duplicate the crucial information */
    gparts[i + Ndm].x[0] = sparts[i].x[0];
    gparts[i + Ndm].x[1] = sparts[i].x[1];
    gparts[i + Ndm].x[2] = sparts[i].x[2];

    gparts[i + Ndm].v_full[0] = sparts[i].v[0];
    gparts[i + Ndm].v_full[1] = sparts[i].v[1];
    gparts[i + Ndm].v_full[2] = sparts[i].v[2];

    gparts[i + Ndm].mass = sparts[i].mass;

    /* Set gpart type */
    gparts[i + Ndm].type = swift_type_star;

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -i;
    sparts[i].gpart = &gparts[i + Ndm];
  }
}

/**
 * @brief Copy every DM #gpart into the dmparts array.
 *
 * @param gparts The array of #gpart containing all particles.
 * @param Ntot The number of #gpart.
 * @param dmparts The array of #gpart containg DM particles to be filled.
 * @param Ndm The number of DM particles.
 */
void io_collect_dm_gparts(const struct gpart* const gparts, size_t Ntot,
                          struct gpart* const dmparts, size_t Ndm) {

  size_t count = 0;

  /* Loop over all gparts */
  for (size_t i = 0; i < Ntot; ++i) {

    /* message("i=%zd count=%zd id=%lld part=%p", i, count, gparts[i].id,
     * gparts[i].part); */

    /* And collect the DM ones */
    if (gparts[i].type == swift_type_dark_matter) {
      dmparts[count] = gparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Ndm)
    error("Collected the wrong number of dm particles (%zu vs. %zu expected)",
          count, Ndm);
}
