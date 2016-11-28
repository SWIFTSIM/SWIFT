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

const char* particle_type_names[NUM_PARTICLE_TYPES] = {
    "Gas", "DM", "Boundary", "Dummy", "Star", "BH"};

/**
 * @brief Converts a C data type to the HDF5 equivalent.
 *
 * This function is a trivial wrapper around the HDF5 types but allows
 * to change the exact storage types matching the code types in a transparent
 *way.
 */
hid_t hdf5Type(enum DATA_TYPE type) {

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
size_t sizeOfType(enum DATA_TYPE type) {

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
int isDoublePrecision(enum DATA_TYPE type) {

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
 * @param type The #DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group.
 *
 * Calls #error() if an error occurs.
 */
void readAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data) {
  hid_t h_attr = 0, h_err = 0;

  h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if (h_attr < 0) {
    error("Error while opening attribute '%s'", name);
  }

  h_err = H5Aread(h_attr, hdf5Type(type), data);
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
 * @param type The #DATA_TYPE of the attribute.
 * @param data The attribute to write.
 * @param num The number of elements to write
 *
 * Calls #error() if an error occurs.
 */
void writeAttribute(hid_t grp, const char* name, enum DATA_TYPE type,
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

  h_attr = H5Acreate1(grp, name, hdf5Type(type), h_space, H5P_DEFAULT);
  if (h_attr < 0) {
    error("Error while creating attribute '%s'.", name);
  }

  h_err = H5Awrite(h_attr, hdf5Type(type), data);
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
void writeStringAttribute(hid_t grp, const char* name, const char* str,
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
void writeAttribute_d(hid_t grp, const char* name, double data) {
  writeAttribute(grp, name, DOUBLE, &data, 1);
}

/**
 * @brief Writes a float value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_f(hid_t grp, const char* name, float data) {
  writeAttribute(grp, name, FLOAT, &data, 1);
}

/**
 * @brief Writes an int value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */

void writeAttribute_i(hid_t grp, const char* name, int data) {
  writeAttribute(grp, name, INT, &data, 1);
}

/**
 * @brief Writes a long value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_l(hid_t grp, const char* name, long data) {
  writeAttribute(grp, name, LONG, &data, 1);
}

/**
 * @brief Writes a string value as an attribute
 * @param grp The group in which to write
 * @param name The name of the attribute
 * @param str The string to write
 */
void writeAttribute_s(hid_t grp, const char* name, const char* str) {
  writeStringAttribute(grp, name, str, strlen(str));
}

/**
 * @brief Reads the Unit System from an IC file.
 * @param h_file The (opened) HDF5 file from which to read.
 * @param us The UnitSystem to fill.
 *
 * If the 'Units' group does not exist in the ICs, cgs units will be assumed
 */
void readUnitSystem(hid_t h_file, struct UnitSystem* us) {

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
  readAttribute(h_grp, "Unit length in cgs (U_L)", DOUBLE,
                &us->UnitLength_in_cgs);
  readAttribute(h_grp, "Unit mass in cgs (U_M)", DOUBLE, &us->UnitMass_in_cgs);
  readAttribute(h_grp, "Unit time in cgs (U_t)", DOUBLE, &us->UnitTime_in_cgs);
  readAttribute(h_grp, "Unit current in cgs (U_I)", DOUBLE,
                &us->UnitCurrent_in_cgs);
  readAttribute(h_grp, "Unit temperature in cgs (U_T)", DOUBLE,
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
void writeUnitSystem(hid_t h_file, const struct UnitSystem* us,
                     const char* groupName) {

  hid_t h_grpunit = 0;
  h_grpunit = H5Gcreate1(h_file, groupName, 0);
  if (h_grpunit < 0) error("Error while creating Unit System group");

  writeAttribute_d(h_grpunit, "Unit mass in cgs (U_M)",
                   units_get_base_unit(us, UNIT_MASS));
  writeAttribute_d(h_grpunit, "Unit length in cgs (U_L)",
                   units_get_base_unit(us, UNIT_LENGTH));
  writeAttribute_d(h_grpunit, "Unit time in cgs (U_t)",
                   units_get_base_unit(us, UNIT_TIME));
  writeAttribute_d(h_grpunit, "Unit current in cgs (U_I)",
                   units_get_base_unit(us, UNIT_CURRENT));
  writeAttribute_d(h_grpunit, "Unit temperature in cgs (U_T)",
                   units_get_base_unit(us, UNIT_TEMPERATURE));

  H5Gclose(h_grpunit);
}

/**
 * @brief Writes the code version to the file
 * @param h_file The (opened) HDF5 file in which to write
 */
void writeCodeDescription(hid_t h_file) {
  hid_t h_grpcode = 0;

  h_grpcode = H5Gcreate1(h_file, "/Code", 0);
  if (h_grpcode < 0) error("Error while creating code group");

  writeAttribute_s(h_grpcode, "Code Version", package_version());
  writeAttribute_s(h_grpcode, "Compiler Name", compiler_name());
  writeAttribute_s(h_grpcode, "Compiler Version", compiler_version());
  writeAttribute_s(h_grpcode, "Git Branch", git_branch());
  writeAttribute_s(h_grpcode, "Git Revision", git_revision());
  writeAttribute_s(h_grpcode, "Configuration options", configuration_options());
  writeAttribute_s(h_grpcode, "CFLAGS", compilation_cflags());
  writeAttribute_s(h_grpcode, "HDF5 library version", hdf5_version());
#ifdef HAVE_FFTW
  writeAttribute_s(h_grpcode, "FFTW library version", fftw3_version());
#endif
#ifdef WITH_MPI
  writeAttribute_s(h_grpcode, "MPI library", mpi_version());
#ifdef HAVE_METIS
  writeAttribute_s(h_grpcode, "METIS library version", metis_version());
#endif
#else
  writeAttribute_s(h_grpcode, "MPI library", "Non-MPI version of SWIFT");
#endif
  H5Gclose(h_grpcode);
}

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
FILE* prepareXMFfile(const char* baseName) {
  char buffer[1024];

  char fileName[FILENAME_BUFFER_SIZE];
  char tempFileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s.xmf", baseName);
  snprintf(tempFileName, FILENAME_BUFFER_SIZE, "%s_temp.xmf", baseName);
  FILE* xmfFile = fopen(fileName, "r");
  FILE* tempFile = fopen(tempFileName, "w");

  if (xmfFile == NULL) error("Unable to open current XMF file.");

  if (tempFile == NULL) error("Unable to open temporary file.");

  /* First we make a temporary copy of the XMF file and count the lines */
  int counter = 0;
  while (fgets(buffer, 1024, xmfFile) != NULL) {
    counter++;
    fprintf(tempFile, "%s", buffer);
  }
  fclose(tempFile);
  fclose(xmfFile);

  /* We then copy the XMF file back up to the closing lines */
  xmfFile = fopen(fileName, "w");
  tempFile = fopen(tempFileName, "r");

  if (xmfFile == NULL) error("Unable to open current XMF file.");

  if (tempFile == NULL) error("Unable to open temporary file.");

  int i = 0;
  while (fgets(buffer, 1024, tempFile) != NULL && i < counter - 3) {
    i++;
    fprintf(xmfFile, "%s", buffer);
  }
  fprintf(xmfFile, "\n");
  fclose(tempFile);
  remove(tempFileName);

  return xmfFile;
}

/**
 * @brief Writes the begin of the XMF file
 *
 * @todo Exploit the XML nature of the XMF format to write a proper XML writer
 *and simplify all the XMF-related stuff.
 */
void createXMFfile(const char* baseName) {

  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s.xmf", baseName);
  FILE* xmfFile = fopen(fileName, "w");

  fprintf(xmfFile, "<?xml version=\"1.0\" ?> \n");
  fprintf(xmfFile, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
  fprintf(
      xmfFile,
      "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.1\">\n");
  fprintf(xmfFile, "<Domain>\n");
  fprintf(xmfFile,
          "<Grid Name=\"TimeSeries\" GridType=\"Collection\" "
          "CollectionType=\"Temporal\">\n\n");

  fprintf(xmfFile, "</Grid>\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");

  fclose(xmfFile);
}

/**
 * @brief Writes the part of the XMF entry presenting the geometry of the
 *snapshot
 *
 * @param xmfFile The file to write in.
 * @param hdfFileName The name of the HDF5 file corresponding to this output.
 * @param time The current simulation time.
 */
void writeXMFoutputheader(FILE* xmfFile, char* hdfFileName, float time) {
  /* Write end of file */

  fprintf(xmfFile, "<!-- XMF description for file: %s -->\n", hdfFileName);
  fprintf(xmfFile,
          "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(xmfFile, "<Time Type=\"Single\" Value=\"%f\"/>\n", time);
}

/**
 * @brief Writes the end of the XMF file (closes all open markups)
 *
 * @param xmfFile The file to write in.
 * @param output The number of this output.
 * @param time The current simulation time.
 */
void writeXMFoutputfooter(FILE* xmfFile, int output, float time) {
  /* Write end of the section of this time step */

  fprintf(xmfFile,
          "\n</Grid> <!-- End of meta-data for output=%03i, time=%f -->\n",
          output, time);
  fprintf(xmfFile, "\n</Grid> <!-- timeSeries -->\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");

  fclose(xmfFile);
}

void writeXMFgroupheader(FILE* xmfFile, char* hdfFileName, size_t N,
                         enum PARTICLE_TYPE ptype) {
  fprintf(xmfFile, "\n<Grid Name=\"%s\" GridType=\"Uniform\">\n",
          particle_type_names[ptype]);
  fprintf(xmfFile,
          "<Topology TopologyType=\"Polyvertex\" Dimensions=\"%zu\"/>\n", N);
  fprintf(xmfFile, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf(xmfFile,
          "<DataItem Dimensions=\"%zu 3\" NumberType=\"Double\" "
          "Precision=\"8\" "
          "Format=\"HDF\">%s:/PartType%d/Coordinates</DataItem>\n",
          N, hdfFileName, (int)ptype);
  fprintf(xmfFile,
          "</Geometry>\n <!-- Done geometry for %s, start of particle fields "
          "list -->\n",
          particle_type_names[ptype]);
}

void writeXMFgroupfooter(FILE* xmfFile, enum PARTICLE_TYPE ptype) {
  fprintf(xmfFile, "</Grid> <!-- End of meta-data for parttype=%s -->\n",
          particle_type_names[ptype]);
}

/**
 * @brief Writes the lines corresponding to an array of the HDF5 output
 *
 * @param xmfFile The file in which to write
 * @param fileName The name of the HDF5 file associated to this XMF descriptor.
 * @param partTypeGroupName The name of the group containing the particles in
 *the HDF5 file.
 * @param name The name of the array in the HDF5 file.
 * @param N The number of particles.
 * @param dim The dimension of the quantity (1 for scalars, 3 for vectors).
 * @param type The type of the data to write.
 *
 * @todo Treat the types in a better way.
 */
void writeXMFline(FILE* xmfFile, const char* fileName,
                  const char* partTypeGroupName, const char* name, size_t N,
                  int dim, enum DATA_TYPE type) {
  fprintf(xmfFile,
          "<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
          name, dim == 1 ? "Scalar" : "Vector");
  if (dim == 1)
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%zu\" NumberType=\"Double\" "
            "Precision=\"%d\" Format=\"HDF\">%s:%s/%s</DataItem>\n",
            N, type == FLOAT ? 4 : 8, fileName, partTypeGroupName, name);
  else
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%zu %d\" NumberType=\"Double\" "
            "Precision=\"%d\" Format=\"HDF\">%s:%s/%s</DataItem>\n",
            N, dim, type == FLOAT ? 4 : 8, fileName, partTypeGroupName, name);
  fprintf(xmfFile, "</Attribute>\n");
}

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
void prepare_dm_gparts(struct gpart* const gparts, size_t Ndm) {

  /* Let's give all these gparts a negative id */
  for (size_t i = 0; i < Ndm; ++i) {
    /* 0 or negative ids are not allowed */
    if (gparts[i].id_or_neg_offset <= 0)
      error("0 or negative ID for DM particle %zu: ID=%lld", i,
            gparts[i].id_or_neg_offset);
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
 *at the start
 * @param Ngas The number of gas particles read in.
 * @param Ndm The number of DM particles read in.
 */
void duplicate_hydro_gparts(struct part* const parts,
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

    /* Link the particles */
    gparts[i + Ndm].id_or_neg_offset = -i;
    parts[i].gpart = &gparts[i + Ndm];
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
void collect_dm_gparts(const struct gpart* const gparts, size_t Ntot,
                       struct gpart* const dmparts, size_t Ndm) {

  size_t count = 0;

  /* Loop over all gparts */
  for (size_t i = 0; i < Ntot; ++i) {

    /* message("i=%zd count=%zd id=%lld part=%p", i, count, gparts[i].id,
     * gparts[i].part); */

    /* And collect the DM ones */
    if (gparts[i].id_or_neg_offset > 0) {
      dmparts[count] = gparts[i];
      count++;
    }
  }

  /* Check that everything is fine */
  if (count != Ndm)
    error("Collected the wrong number of dm particles (%zu vs. %zu expected)",
          count, Ndm);
}

#endif
