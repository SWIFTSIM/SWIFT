/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 James Willis (james.s.willis@durham.ac.uk).
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

#include "parser.h"
#include <assert.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]) {

  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Create variables that will be set from the parameter file. */
  int no_of_threads = 0;
  int no_of_time_steps = 0;
  float max_h = 0.0f;
  double start_time = 0.0;
  char ic_file[PARSER_MAX_LINE_SIZE];

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Print the contents of the structure. */
  parser_print_params(&param_file);

  /* Retrieve parameters and store them in variables defined above.
   * Have to specify the name of the parameter as it appears in the
   * input file: testParserInput.yaml.*/
  parser_get_param_int(&param_file, "no_of_threads", &no_of_threads);
  parser_get_param_int(&param_file, "no_of_time_steps", &no_of_time_steps);
  parser_get_param_float(&param_file, "max_h", &max_h);
  parser_get_param_double(&param_file, "start_time", &start_time);
  parser_get_param_string(&param_file, "ic_file", ic_file);

  /* Print the variables to check their values are correct. */
  printf(
      "no_of_threads: %d, no_of_time_steps: %d, max_h: %f, start_time: %lf, "
      "ic_file: %s\n",
      no_of_threads, no_of_time_steps, max_h, start_time, ic_file);

  assert(no_of_threads == 16);
  assert(no_of_time_steps == 10);
  assert(fabs(max_h - 1.1255) < 0.00001);
  assert(fabs(start_time - 1.23456789) < 0.00001);
  assert(strcmp(ic_file, "ic_file.ini") == 0); /*strcmp returns 0 if correct.*/

  return 0;
}
