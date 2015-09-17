/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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
#include "const.h"
#include "error.h"
#include "kernel.h"
#include "version.h"

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
void writeAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data,
                    int num) {
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
void writeStringAttribute(hid_t grp, char* name, const char* str, int length) {
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
    error("Error while resizing attribute tyep to '%i'.", length);
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
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_d(hid_t grp, char* name, double data) {
  writeAttribute(grp, name, DOUBLE, &data, 1);
}

/**
 * @brief Writes a float value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_f(hid_t grp, char* name, float data) {
  writeAttribute(grp, name, FLOAT, &data, 1);
}

/**
 * @brief Writes an int value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */

void writeAttribute_i(hid_t grp, char* name, int data) {
  writeAttribute(grp, name, INT, &data, 1);
}

/**
 * @brief Writes a long value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_l(hid_t grp, char* name, long data) {
  writeAttribute(grp, name, LONG, &data, 1);
}

/**
 * @brief Writes a string value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param str The string to write
 */
void writeAttribute_s(hid_t grp, char* name, const char* str) {
  writeStringAttribute(grp, name, str, strlen(str));
}

/* ------------------------------------------------------------------------------------------------
 * This part writes the XMF file descriptor enabling a visualisation through
 * ParaView
 * ------------------------------------------------------------------------------------------------
 */
/**
 * @brief Writes the current model of SPH to the file
 * @param h_file The (opened) HDF5 file in which to write
 */
void writeSPHflavour(hid_t h_file) {
  hid_t h_grpsph = 0;

  h_grpsph = H5Gcreate1(h_file, "/SPH", 0);
  if (h_grpsph < 0) error("Error while creating SPH group");

  writeAttribute_f(h_grpsph, "Kernel eta", const_eta_kernel);
  writeAttribute_f(h_grpsph, "Weighted N_ngb", kernel_nwneigh);
  writeAttribute_f(h_grpsph, "Delta N_ngb", const_delta_nwneigh);
  writeAttribute_f(h_grpsph, "Hydro gamma", const_hydro_gamma);

#ifdef LEGACY_GADGET2_SPH
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model",
                   "(No treatment) Legacy Gadget-2 as in Springel (2005)");
  writeAttribute_s(h_grpsph, "Viscosity Model",
                   "Legacy Gadget-2 as in Springel (2005)");
  writeAttribute_f(h_grpsph, "Viscosity alpha", const_viscosity_alpha);
  writeAttribute_f(h_grpsph, "Viscosity beta", 3.f);
#else
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model",
                   "Price (2008) without switch");
  writeAttribute_f(h_grpsph, "Thermal Conductivity alpha",
                   const_conductivity_alpha);
  writeAttribute_s(h_grpsph, "Viscosity Model",
                   "Morris & Monaghan (1997), Rosswog, Davies, Thielemann & "
                   "Piran (2000) with additional Balsara (1995) switch");
  writeAttribute_f(h_grpsph, "Viscosity alpha_min", const_viscosity_alpha_min);
  writeAttribute_f(h_grpsph, "Viscosity alpha_max", const_viscosity_alpha_max);
  writeAttribute_f(h_grpsph, "Viscosity beta", 2.f);
  writeAttribute_f(h_grpsph, "Viscosity decay length", const_viscosity_length);
#endif

  writeAttribute_f(h_grpsph, "CFL parameter", const_cfl);
  writeAttribute_f(h_grpsph, "Maximal ln(Delta h) change over dt",
                   const_ln_max_h_change);
  writeAttribute_f(h_grpsph, "Maximal Delta h change over dt",
                   exp(const_ln_max_h_change));
  writeAttribute_f(h_grpsph, "Maximal Delta u change over dt",
                   const_max_u_change);
  writeAttribute_s(h_grpsph, "Kernel", kernel_name);

  H5Gclose(h_grpsph);
}

/**
 * @brief Writes the current Unit System
 * @param h_file The (opened) HDF5 file in which to write
 * @param us The UnitSystem used in the run
 */
void writeUnitSystem(hid_t h_file, struct UnitSystem* us) {
  hid_t h_grpunit = 0;

  h_grpunit = H5Gcreate1(h_file, "/Units", 0);
  if (h_grpunit < 0) error("Error while creating Unit System group");

  writeAttribute_d(h_grpunit, "Unit mass in cgs (U_M)",
                   getBaseUnit(us, UNIT_MASS));
  writeAttribute_d(h_grpunit, "Unit length in cgs (U_L)",
                   getBaseUnit(us, UNIT_LENGTH));
  writeAttribute_d(h_grpunit, "Unit time in cgs (U_t)",
                   getBaseUnit(us, UNIT_TIME));
  writeAttribute_d(h_grpunit, "Unit current in cgs (U_I)",
                   getBaseUnit(us, UNIT_CURRENT));
  writeAttribute_d(h_grpunit, "Unit temperature in cgs (U_T)",
                   getBaseUnit(us, UNIT_TEMPERATURE));

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
  writeAttribute_s(h_grpcode, "Git Branch", git_branch());
  writeAttribute_s(h_grpcode, "Git Revision", git_revision());

  H5Gclose(h_grpcode);
}

/**
 * @brief Prepares the XMF file for the new entry
 *
 * Creates a temporary file on the disk in order to copy the right lines.
 *
 * @todo Use a proper XML library to avoid stupid copies.
 */
FILE* prepareXMFfile() {
  char buffer[1024];

  FILE* xmfFile = fopen("output.xmf", "r");
  FILE* tempFile = fopen("output_temp.xmf", "w");

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
  xmfFile = fopen("output.xmf", "w");
  tempFile = fopen("output_temp.xmf", "r");

  if (xmfFile == NULL) error("Unable to open current XMF file.");

  if (tempFile == NULL) error("Unable to open temporary file.");

  int i = 0;
  while (fgets(buffer, 1024, tempFile) != NULL && i < counter - 3) {
    i++;
    fprintf(xmfFile, "%s", buffer);
  }
  fprintf(xmfFile, "\n");
  fclose(tempFile);
  remove("output_temp.xmf");

  return xmfFile;
}

/**
 * @brief Writes the begin of the XMF file
 *
 * @todo Exploit the XML nature of the XMF format to write a proper XML writer
 *and simplify all the XMF-related stuff.
 */
void createXMFfile() {
  FILE* xmfFile = fopen("output.xmf", "w");

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
 * @param Nparts The number of particles.
 * @param hdfFileName The name of the HDF5 file corresponding to this output.
 * @param time The current simulation time.
 */
void writeXMFheader(FILE* xmfFile, long long Nparts, char* hdfFileName,
                    float time) {
  /* Write end of file */

  fprintf(xmfFile,
          "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(xmfFile, "<Time Type=\"Single\" Value=\"%f\"/>\n", time);
  fprintf(xmfFile, "<Grid Name=\"Gas\" GridType=\"Uniform\">\n");
  fprintf(xmfFile,
          "<Topology TopologyType=\"Polyvertex\" Dimensions=\"%lld\"/>\n",
          Nparts);
  fprintf(xmfFile, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf(xmfFile,
          "<DataItem Dimensions=\"%lld 3\" NumberType=\"Double\" "
          "Precision=\"8\" "
          "Format=\"HDF\">%s:/PartType0/Coordinates</DataItem>\n",
          Nparts, hdfFileName);
  fprintf(xmfFile, "</Geometry>");
}

/**
 * @brief Writes the end of the XMF file (closes all open markups)
 *
 * @param xmfFile The file to write in.
 */
void writeXMFfooter(FILE* xmfFile) {
  /* Write end of the section of this time step */

  fprintf(xmfFile, "\n</Grid>\n");
  fprintf(xmfFile, "</Grid>\n");
  fprintf(xmfFile, "\n</Grid>\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");

  fclose(xmfFile);
}

/**
 * @brief Writes the lines corresponding to an array of the HDF5 output
 *
 * @param xmfFile The file in which to write
 * @param fileName The name of the HDF5 file associated to this XMF descriptor.
 * @param name The name of the array in the HDF5 file.
 * @param N The number of particles.
 * @param dim The dimension of the quantity (1 for scalars, 3 for vectors).
 * @param type The type of the data to write.
 *
 * @todo Treat the types in a better way.
 */
void writeXMFline(FILE* xmfFile, char* fileName, char* name, long long N,
                  int dim, enum DATA_TYPE type) {
  fprintf(xmfFile,
          "<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
          name, dim == 1 ? "Scalar" : "Vector");
  if (dim == 1)
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%lld\" NumberType=\"Double\" "
            "Precision=\"%d\" Format=\"HDF\">%s:/PartType0/%s</DataItem>\n",
            N, type == FLOAT ? 4 : 8, fileName, name);
  else
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%lld %d\" NumberType=\"Double\" "
            "Precision=\"%d\" Format=\"HDF\">%s:/PartType0/%s</DataItem>\n",
            N, dim, type == FLOAT ? 4 : 8, fileName, name);
  fprintf(xmfFile, "</Attribute>\n");
}

#endif
