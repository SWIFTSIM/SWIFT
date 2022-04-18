/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017  Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_XMF_H
#define SWIFT_XMF_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "common_io.h"
#include "part_type.h"

void xmf_create_file(const char* fileName);
FILE* xmf_prepare_file(const char* fileName);
void xmf_write_outputheader(FILE* xmfFile, const char* hdfFileName, float time);
void xmf_write_outputfooter(FILE* xmfFile, int outputCount, float time);
void xmf_write_groupheader(FILE* xmfFile, const char* hdfFileName,
                           const int distributed, size_t N,
                           enum part_type ptype);
void xmf_write_groupfooter(FILE* xmfFile, enum part_type ptype);
void xmf_write_line(FILE* xmfFile, const char* fileName, const int distributed,
                    const char* partTypeGroupName, const char* name, size_t N,
                    int dim, enum IO_DATA_TYPE type);

#endif /* SWIFT_XMF_H */
