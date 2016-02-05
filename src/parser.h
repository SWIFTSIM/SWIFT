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
#define VALUE ':'
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
void parser_read_file(const char *file_name, struct swift_params *params);
void parser_print_params(struct swift_params *params);
void parser_get_param_int(struct swift_params *params, char *name, int *retParam);
void parser_get_param_float(struct swift_params *params, char *name, float *retParam);
void parser_get_param_string(struct swift_params *params, char *name, char *retParam);

/* Private functions. */
static void read_param(FILE *fp, struct swift_params *params);
static int count_char(char *str, char val); 

#endif /* SWIFT_PARSER_H */
