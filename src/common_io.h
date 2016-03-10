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

/* Includes. */
#include "part.h"
#include "units.h"

#if defined(HAVE_HDF5)

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
 * @brief The two sorts of data present in the GADGET IC files: compulsory to
 *start a run or optional.
 *
 */
enum DATA_IMPORTANCE {
  COMPULSORY = 1,
  OPTIONAL = 0
};

/**
 * @brief The different particle types present in a GADGET IC file
 *
 */
enum PARTICLE_TYPE {
  GAS = 0,
  DM = 1,
  STAR = 4,
  BH = 5
};

hid_t hdf5Type(enum DATA_TYPE type);
size_t sizeOfType(enum DATA_TYPE type);

void collect_dm_gparts(struct gpart* gparts, int Ntot, struct gpart* dmparts,
                       int Ndm);
void prepare_dm_gparts(struct gpart* gparts, int Ndm);
void duplicate_hydro_gparts(struct part* parts, struct gpart* gparts, int Ngas,
                            int Ndm);

void readAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data);

void writeAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data,
                    int num);

void writeAttribute_d(hid_t grp, char* name, double data);
void writeAttribute_f(hid_t grp, char* name, float data);
void writeAttribute_i(hid_t grp, char* name, int data);
void writeAttribute_l(hid_t grp, char* name, long data);
void writeAttribute_s(hid_t grp, char* name, const char* str);

void createXMFfile();
FILE* prepareXMFfile();
void writeXMFfooter(FILE* xmfFile);
void writeXMFheader(FILE* xmfFile, long long N, char* hdfFileName, float time);
void writeXMFline(FILE* xmfFile, char* fileName, char* name, long long N,
                  int dim, enum DATA_TYPE type);

void writeCodeDescription(hid_t h_file);
void writeSPHflavour(hid_t h_file);
void writeUnitSystem(hid_t h_file, struct UnitSystem* us);

#endif

#endif /* SWIFT_COMMON_IO_H */
