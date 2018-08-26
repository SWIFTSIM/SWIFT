/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 James Willis (james.s.willis@durham.ac.uk).
 *               2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "parser.h"

int main(int argc, char *argv[]) {
  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Print the contents of the structure to stdout. */
  printf("\n--- Values read from file:\n");
  parser_print_params(&param_file);

  /* Retrieve parameters and store them in variables defined above.
   * Have to specify the name of the parameter as it appears in the
   * input file: testParserInput.yaml.*/
  const int no_of_threads =
      parser_get_param_int(&param_file, "Scheduler:no_of_threads");
  const int no_of_time_steps =
      parser_get_param_int(&param_file, "Simulation:no_of_time_steps");
  const float max_h = parser_get_param_float(&param_file, "Simulation:max_h");
  const double start_time =
      parser_get_param_double(&param_file, "Simulation:start_time");
  const int kernel = parser_get_param_int(&param_file, "kernel");

  const int optional1 =
      parser_get_opt_param_int(&param_file, "Simulation:optional", 1);

  /* Check if we can read it again */
  const int optional2 =
      parser_get_opt_param_int(&param_file, "Simulation:optional", 1);

  /* Optional things not mentioned in parameter file. Should be in output
   * files.*/
  parser_get_opt_param_int(&param_file, "Virtual:param1", 1);
  parser_get_opt_param_int(&param_file, "Virtual:param2", 2);
  parser_get_opt_param_int(&param_file, "Virtual:param3", 3);
  parser_get_opt_param_int(&param_file, "Virtual:param4", 4);

  /* Check if we can read it again */
  parser_get_opt_param_int(&param_file, "Virtual:param1", 1);

  char ic_file[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(&param_file, "IO:ic_file", ic_file);

  char csides[3];
  parser_get_param_char_array(&param_file, "Box:csides", 3, csides);

  int isides[3];
  parser_get_param_int_array(&param_file, "Box:isides", 3, isides);

  float fsides[3];
  parser_get_param_float_array(&param_file, "Box:fsides", 3, fsides);

  double dsides[3];
  parser_get_param_double_array(&param_file, "Box:dsides", 3, dsides);

  int optsides[5] = {1, 2, 3, 4, 5};
  int haveopt1 =
      parser_get_opt_param_int_array(&param_file, "Box:moresides", 5, optsides);
  /* Check if we can read it again */
  int haveopt2 =
      parser_get_opt_param_int_array(&param_file, "Box:moresides", 5, optsides);

  char **var_result;
  int nvar_result;
  parser_get_param_string_array(&param_file, "Words:list", &nvar_result,
                                &var_result);

  printf("\nList from Words:list parameter\n");
  for (int i = 0; i < nvar_result; i++) printf("   %d: %s\n", i, var_result[i]);

  /* Get same list without []. */
  char **var_result2;
  int nvar_result2;
  parser_get_param_string_array(&param_file, "Words:nonstdlist", &nvar_result2,
                                &var_result2);

  assert(nvar_result == nvar_result2);
  for (int i = 0; i < nvar_result; i++)
    assert(strcmp(var_result[i], var_result2[i]) == 0);

  parser_free_param_string_array(nvar_result, var_result);
  parser_free_param_string_array(nvar_result2, var_result2);

  const char *optwords[4] = {"Word1", "Word2", "Word3", "Word4"};
  int noptwords = 4;
  int haveoptwords1 = parser_get_opt_param_string_array(
      &param_file, "Simulation:optwords", &nvar_result, &var_result, noptwords,
      optwords);
  /* Check if we can read it again */
  int haveoptwords2 = parser_get_opt_param_string_array(
      &param_file, "Simulation:optwords", &nvar_result, &var_result, noptwords,
      optwords);
  printf("\nList from Simulation:optwords parameter (%d of %d)\n", nvar_result,
         noptwords);
  assert(noptwords == nvar_result);
  for (int i = 0; i < nvar_result; i++) {
    assert(strcmp(optwords[i], var_result[i]) == 0);
    printf("   %d: %s\n", i, var_result[i]);
  }
  parser_free_param_string_array(nvar_result, var_result);

  /* Long list of values. */
  parser_get_param_string_array(&param_file, "Words:long", &nvar_result,
                                &var_result);
  printf("No of letters in alphabet = %d\n", nvar_result);
  assert(nvar_result == 26);
  char alphabet[27];
  for (int i = 0; i < nvar_result; i++) {
    alphabet[i] = var_result[i][0];
  }
  alphabet[26] = '\0';
  printf("The alphabet: %s\n", alphabet);
  assert(strcmp("abcdefghijklmnopqrstuvwxyz", alphabet) == 0);
  parser_free_param_string_array(nvar_result, var_result);

  /* Print the contents of the structure to stdout now used. */
  printf("\n--- Values after being used:\n");
  parser_print_params(&param_file);

  /* Print the variables to check their values are correct. */
  printf(
      "\nValues read from file:\n"
      "no_of_threads: %d, no_of_time_steps: %d, max_h: %f, start_time: %lf,"
      " ic_file: %s, kernel: %d optional: %d\n",
      no_of_threads, no_of_time_steps, max_h, start_time, ic_file, kernel,
      optional1);
  printf("Box: [%d, %d, %d]\n", isides[0], isides[1], isides[2]);

  /* Print the contents of the structure to a file in YAML format. */
  parser_write_params_to_file(&param_file, "used_parser_output.yml", 1);
  parser_write_params_to_file(&param_file, "unused_parser_output.yml", 0);

  assert(no_of_threads == 16);
  assert(no_of_time_steps == 10);
  assert(fabs(max_h - 1.1255) < 0.00001);
  assert(fabs(start_time - 1.23456789) < 0.00001);
  assert(strcmp(ic_file, "ic_file.ini") == 0); /*strcmp returns 0 if correct.*/
  assert(kernel == 4);
  assert(optional1 == 1);
  assert(optional2 == 1);

  assert(csides[0] == 'a');
  assert(csides[1] == 'b');
  assert(csides[2] == 'c');

  assert(isides[0] == 2);
  assert(isides[1] == 3);
  assert(isides[2] == 4);

  assert(fsides[0] == 2.0);
  assert(fsides[1] == 3.0);
  assert(fsides[2] == 4.0);

  assert(dsides[0] == 2.0);
  assert(dsides[1] == 3.0);
  assert(dsides[2] == 4.0);

  assert(haveopt1 == 0);
  assert(haveopt2 == 1);
  assert(optsides[0] == 1);
  assert(optsides[1] == 2);
  assert(optsides[2] == 3);
  assert(optsides[3] == 4);
  assert(optsides[4] == 5);

  assert(haveoptwords1 == 0);
  assert(haveoptwords2 == 1);

  return 0;
}
