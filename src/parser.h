/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
 *               2017 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Standard headers */
#include <stdio.h>

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Some constants. */
#define PARSER_MAX_LINE_SIZE 256
#define PARSER_MAX_NO_OF_PARAMS 600
#define PARSER_MAX_NO_OF_SECTIONS 64

/* A parameter in the input file */
struct parameter {
  char name[PARSER_MAX_LINE_SIZE];
  char value[PARSER_MAX_LINE_SIZE];
  int used;
  int is_default;
};

struct section {
  char name[PARSER_MAX_LINE_SIZE];
};

/* The array of parameters read from a file */
struct swift_params {
  struct section section[PARSER_MAX_NO_OF_SECTIONS];
  struct parameter data[PARSER_MAX_NO_OF_PARAMS];
  int sectionCount;
  int paramCount;
  char fileName[PARSER_MAX_LINE_SIZE];
};

/* Public API. */
void parser_init(const char *file_name, struct swift_params *params);
void parser_read_file(const char *file_name, struct swift_params *params);
void parser_print_params(const struct swift_params *params);
void parser_write_params_to_file(const struct swift_params *params,
                                 const char *file_name, int write_all);
void parser_set_param(struct swift_params *params, const char *desc);

char parser_get_param_char(struct swift_params *params, const char *name);
int parser_get_param_int(struct swift_params *params, const char *name);
float parser_get_param_float(struct swift_params *params, const char *name);
double parser_get_param_double(struct swift_params *params, const char *name);
long long parser_get_param_longlong(struct swift_params *params,
                                    const char *name);
void parser_get_param_string(struct swift_params *params, const char *name,
                             char *retParam);

int parser_does_param_exist(struct swift_params *params, const char *name);

char parser_get_opt_param_char(struct swift_params *params, const char *name,
                               char def);
int parser_get_opt_param_int(struct swift_params *params, const char *name,
                             int def);
float parser_get_opt_param_float(struct swift_params *params, const char *name,
                                 float def);
double parser_get_opt_param_double(struct swift_params *params,
                                   const char *name, double def);
long long parser_get_opt_param_longlong(struct swift_params *params,
                                        const char *name, long long def);
void parser_get_opt_param_string(struct swift_params *params, const char *name,
                                 char *retParam, const char *def);
void parser_get_param_char_array(struct swift_params *params, const char *name,
                                 int nval, char *values);
void parser_get_param_int_array(struct swift_params *params, const char *name,
                                int nval, int *values);
void parser_get_param_float_array(struct swift_params *params, const char *name,
                                  int nval, float *values);
void parser_get_param_double_array(struct swift_params *params,
                                   const char *name, int nval, double *values);
void parser_get_param_longlong_array(struct swift_params *params,
                                     const char *name, int nval,
                                     long long *values);
int parser_get_opt_param_char_array(struct swift_params *params,
                                    const char *name, int nval, char *values);
int parser_get_opt_param_int_array(struct swift_params *params,
                                   const char *name, int nval, int *values);
int parser_get_opt_param_float_array(struct swift_params *params,
                                     const char *name, int nval, float *values);
int parser_get_opt_param_double_array(struct swift_params *params,
                                      const char *name, int nval,
                                      double *values);
int parser_get_opt_param_longlong_array(struct swift_params *params,
                                        const char *name, int nval,
                                        long long *values);
void parser_get_param_string_array(struct swift_params *params,
                                   const char *name, int *nval, char ***values);
int parser_get_opt_param_string_array(struct swift_params *params,
                                      const char *name, int *nval,
                                      char ***values, int ndef,
                                      const char *def[]);
void parser_free_param_string_array(int nval, char **values);

#if defined(HAVE_HDF5)
void parser_write_params_to_hdf5(const struct swift_params *params, hid_t grp,
                                 int write_all);
#endif

/* Dump/restore. */
void parser_struct_dump(const struct swift_params *params, FILE *stream);
void parser_struct_restore(const struct swift_params *params, FILE *stream);

/* Lookup functions */
int parser_get_section_id(const struct swift_params *params, const char *name);

/* Compare two param structs for changed values. */
int parser_compare_params(const struct swift_params *refparams,
                          struct swift_params *params);
#endif /* SWIFT_PARSER_H */
