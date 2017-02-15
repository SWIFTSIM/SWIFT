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
#ifndef SWIFT_COMMON_IO_H
#define SWIFT_COMMON_IO_H

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)

#include "part.h"
#include "units.h"

/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum DATA_TYPE {
  INT,
  LONG,
  LONGLONG,
  UINT,
  ULONG,
  ULONGLONG,
  FLOAT,
  DOUBLE,
  CHAR
};

/**
 * @brief The different particle types present in a GADGET IC file
 *
 */
enum PARTICLE_TYPE {
  GAS = 0,
  DM = 1,
  BOUNDARY = 2,
  DUMMY = 3,
  STAR = 4,
  BH = 5,
  NUM_PARTICLE_TYPES
};

extern const char* particle_type_names[];

#define FILENAME_BUFFER_SIZE 150
#define FIELD_BUFFER_SIZE 200
#define PARTICLE_GROUP_BUFFER_SIZE 50

hid_t hdf5Type(enum DATA_TYPE type);
size_t sizeOfType(enum DATA_TYPE type);
int isDoublePrecision(enum DATA_TYPE type);

void collect_dm_gparts(const struct gpart* const gparts, size_t Ntot,
                       struct gpart* const dmparts, size_t Ndm);
void prepare_dm_gparts(struct gpart* const gparts, size_t Ndm);
void duplicate_hydro_gparts(struct part* const parts,
                            struct gpart* const gparts, size_t Ngas,
                            size_t Ndm);
void duplicate_star_gparts(struct spart* const sparts,
                           struct gpart* const gparts, size_t Nstars,
                           size_t Ndm);

void readAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data);

void writeAttribute(hid_t grp, const char* name, enum DATA_TYPE type,
                    void* data, int num);

void writeAttribute_d(hid_t grp, const char* name, double data);
void writeAttribute_f(hid_t grp, const char* name, float data);
void writeAttribute_i(hid_t grp, const char* name, int data);
void writeAttribute_l(hid_t grp, const char* name, long data);
void writeAttribute_s(hid_t grp, const char* name, const char* str);

void createXMFfile(const char* baseName);
FILE* prepareXMFfile(const char* baseName);
void writeXMFoutputheader(FILE* xmfFile, char* hdfFileName, float time);
void writeXMFoutputfooter(FILE* xmfFile, int outputCount, float time);
void writeXMFgroupheader(FILE* xmfFile, char* hdfFileName, size_t N,
                         enum PARTICLE_TYPE ptype);
void writeXMFgroupfooter(FILE* xmfFile, enum PARTICLE_TYPE ptype);
void writeXMFline(FILE* xmfFile, const char* fileName,
                  const char* partTypeGroupName, const char* name, size_t N,
                  int dim, enum DATA_TYPE type);

void writeCodeDescription(hid_t h_file);

void readUnitSystem(hid_t h_file, struct UnitSystem* us);
void writeUnitSystem(hid_t h_grp, const struct UnitSystem* us,
                     const char* groupName);

#endif /* defined HDF5 */

#endif /* SWIFT_COMMON_IO_H */
