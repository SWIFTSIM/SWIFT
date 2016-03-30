/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
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

#define PARSER_MAX_LINE_SIZE 128
#define PARSER_MAX_NO_OF_PARAMS 512

#define PARSER_COMMENT_CHAR "#"
#define PARSER_VALUE_CHAR ':'
#define PARSER_VALUE_STRING ":"
#define PARSER_END_OF_FILE "..."

struct parameter {
  char name[PARSER_MAX_LINE_SIZE];
  char value[PARSER_MAX_LINE_SIZE];
};

struct swift_params {
  struct parameter data[PARSER_MAX_NO_OF_PARAMS];
  int count;
};

/* Public API. */
void parser_read_file(const char *file_name, struct swift_params *params);
void parser_print_params(struct swift_params *params);
void parser_get_param_int(struct swift_params *params, char *name,
                          int *retParam);
void parser_get_param_float(struct swift_params *params, char *name,
                            float *retParam);
void parser_get_param_double(struct swift_params *params, char *name,
                             double *retParam);
void parser_get_param_string(struct swift_params *params, char *name,
                             char *retParam);

#endif /* SWIFT_PARSER_H */
