/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *                    Peter W. Draper   (p.w.draper@durham.ac.uk)
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

/* Some standard headers. */
#include <libgen.h>
#include <stdio.h>
#include <string.h>

/* This object's header. */
#include "xmf.h"

/* Local headers. */
#include "common_io.h"
#include "error.h"

/**
 * @brief Return the basename of an HDF5 path.
 *
 * Need basename as XML paths are relative to the container, and XMF file is
 * written with the same baseName as the HDF5 snapshots.
 *
 * @param hdfFileName
 * @return the basename part of hdfFileName.
 */
static const char* xmf_basename(const char* hdfFileName) {
  static char buffer[FILENAME_BUFFER_SIZE];
  strcpy(buffer, hdfFileName);
  return basename(buffer);
}

/**
 * @brief Prepare the XMF file corresponding to a snapshot.
 *
 * @param baseName The common part of the file name.
 */
FILE* xmf_prepare_file(const char* baseName) {
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
 * and simplify all the XMF-related stuff.
 */
void xmf_create_file(const char* baseName) {

  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s.xmf", baseName);
  FILE* xmfFile = fopen(fileName, "w");
  if (xmfFile == NULL) error("Unable to create XMF file.");

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
 * snapshot
 *
 * @param xmfFile The file to write in.
 * @param hdfFileName The name of the HDF5 file corresponding to this output.
 * @param time The current simulation time.
 */
void xmf_write_outputheader(FILE* xmfFile, char* hdfFileName, float time) {
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
void xmf_write_outputfooter(FILE* xmfFile, int output, float time) {
  /* Write end of the section of this time step */

  fprintf(xmfFile,
          "\n</Grid> <!-- End of meta-data for output=%04i, time=%f -->\n",
          output, time);
  fprintf(xmfFile, "\n</Grid> <!-- timeSeries -->\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");

  fclose(xmfFile);
}

/**
 * @brief Writes the header of an XMF group for a given particle type.
 *
 * @param xmfFile The file to write to.
 * @param hdfFileName The name of the corresponding HDF5 file.
 * @param N The number of particles to write.
 * @param ptype The particle type we are writing.
 */
void xmf_write_groupheader(FILE* xmfFile, char* hdfFileName, size_t N,
                           enum part_type ptype) {

  fprintf(xmfFile, "\n<Grid Name=\"%s\" GridType=\"Uniform\">\n",
          part_type_names[ptype]);
  fprintf(xmfFile,
          "<Topology TopologyType=\"Polyvertex\" Dimensions=\"%zu\"/>\n", N);
  fprintf(xmfFile, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf(xmfFile,
          "<DataItem Dimensions=\"%zu 3\" NumberType=\"Double\" "
          "Precision=\"8\" "
          "Format=\"HDF\">%s:/PartType%d/Coordinates</DataItem>\n",
          N, xmf_basename(hdfFileName), (int)ptype);
  fprintf(xmfFile,
          "</Geometry>\n <!-- Done geometry for %s, start of particle fields "
          "list -->\n",
          part_type_names[ptype]);
}

/**
 * @brief Writes the footer of an XMF group for a given particle type.
 *
 * @param xmfFile The file to write to.
 * @param ptype The particle type we are writing.
 */
void xmf_write_groupfooter(FILE* xmfFile, enum part_type ptype) {
  fprintf(xmfFile, "</Grid> <!-- End of meta-data for parttype=%s -->\n",
          part_type_names[ptype]);
}

/**
 * @brief Returns the precision of a given dataset type
 */
int xmf_precision(enum IO_DATA_TYPE type) {
  switch (type) {
    case INT:
    case FLOAT:
      return 4;
      break;
    case DOUBLE:
      return 8;
      break;
    case ULONGLONG:
    case LONGLONG:
      return 8;
      break;
    case CHAR:
      return 1;
      break;
    default:
      error("Unsupported type");
  }
  return 0;
}

/**
 * @brief Returns the Xdmf type name of a given dataset type
 */
const char* xmf_type(enum IO_DATA_TYPE type) {
  switch (type) {
    case FLOAT:
    case DOUBLE:
      return "Float";
      break;
    case INT:
    case ULONGLONG:
    case LONGLONG:
      return "Int";
      break;
    case CHAR:
      return "Char";
      break;
    default:
      error("Unsupported type");
  }
  return "";
}

/**
 * @brief Writes the lines corresponding to an array of the HDF5 output
 *
 * @param xmfFile The file in which to write
 * @param fileName The name of the HDF5 file associated to this XMF descriptor.
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param name The name of the array in the HDF5 file.
 * @param N The number of particles.
 * @param dim The dimension of the quantity (1 for scalars, 3 for vectors).
 * @param type The type of the data to write.
 *
 * @todo Treat the types in a better way.
 */
void xmf_write_line(FILE* xmfFile, const char* fileName,
                    const char* partTypeGroupName, const char* name, size_t N,
                    int dim, enum IO_DATA_TYPE type) {
  fprintf(xmfFile,
          "<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
          name, dim == 1 ? "Scalar" : "Vector");
  if (dim == 1)
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%zu\" NumberType=\"%s\" "
            "Precision=\"%d\" Format=\"HDF\">%s:%s/%s</DataItem>\n",
            N, xmf_type(type), xmf_precision(type), xmf_basename(fileName),
            partTypeGroupName, name);
  else
    fprintf(xmfFile,
            "<DataItem Dimensions=\"%zu %d\" NumberType=\"%s\" "
            "Precision=\"%d\" Format=\"HDF\">%s:%s/%s</DataItem>\n",
            N, dim, xmf_type(type), xmf_precision(type), xmf_basename(fileName),
            partTypeGroupName, name);
  fprintf(xmfFile, "</Attribute>\n");
}
