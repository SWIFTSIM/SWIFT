/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PARSER_H
#define SWIFT_PARSER_H

#include <stdio.h>

#define MAX_LINE_SIZE 128
#define MAX_NO_OF_PARAMS 512 

#define COMMENT "#"
#define VALUE char ":"
#define END_OF_FILE "..."

struct parameter {
    char name [MAX_LINE_SIZE];
    char value [MAX_LINE_SIZE];
};

struct swift_params {
    struct parameter data [MAX_NO_OF_PARAMS];
    int count;
};

/* Public API. */
void parseFile(const char *file_name, struct swift_params *params);
void printParameters(struct swift_params *params);
void getParamInt(struct swift_params *params, char *name, int *retParam);
void getParamFloat(struct swift_params *params, char *name, float *retParam);
void getParamString(struct swift_params *params, char *name, char *retParam);

/* Private functions. */
static void readParameter(FILE *fp, struct swift_params *params);

#endif /* SWIFT_PARSER_H */
