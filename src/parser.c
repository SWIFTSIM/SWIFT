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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
/* Needs to be included so that strtok returns char * instead of a int *. */
#include <string.h>
#include <stdlib.h>

/* This object's header. */
#include "parser.h"

/* Local headers. */
#include "error.h"

/* Private functions. */
static int count_char(char *str, char val);
static void parse_line(FILE *fp, struct swift_params *params);

/**
 * @brief Reads an input file and stores each parameter in a structure.
 *
 * @param file_name Name of file to be read
 * @param params Structure to be populated from file
 */

void parser_read_file(const char *file_name, struct swift_params *params) {

  FILE *fp;

  params->count = 0;

  /* Open file for reading */
  fp = fopen(file_name, "r");

  if (fp == NULL) {
    error("Error opening parameter file: %s", file_name);
  }

  /* Read until the end of the file is reached.*/
  while (!feof(fp)) {
    parse_line(fp, params);
  }

  fclose(fp);
}

/**
 * @brief Counts the number of times a specific character appears in a string.
 *
 * @param str String to be checked
 * @param val Character to be counted
 */

static int count_char(char *str, char val) {

  int count = 0;

  /* Check if the line contains the character */
  while (*str) {
    if (*str++ == val) ++count;
  }

  return count;
}

/**
 * @brief Parses a line from a file and stores any parameters in a structure.
 *
 * @param fp File pointer to file to be read
 * @param params Structure to be populated from file
 *
 */

static void parse_line(FILE *fp, struct swift_params *params) {

  char line[PARSER_MAX_LINE_SIZE];
  char trim_line[PARSER_MAX_LINE_SIZE];

  /* Read a line of the file */
  if (fgets(line, PARSER_MAX_LINE_SIZE, fp) != NULL) {

    char *token;
    /* Remove comments */
    token = strtok(line, PARSER_COMMENT_CHAR);
    strcpy(trim_line, token);

    /* Check if the line contains a value */
    if (strchr(trim_line, PARSER_VALUE_CHAR)) {
      /* Check for more than one parameter on the same line. */
      if (count_char(trim_line, PARSER_VALUE_CHAR) > 1) {
        error("Found more than one parameter in '%s', only one allowed.", line);
      } else {
        /* Take first token as the parameter name. */
        token = strtok(trim_line, PARSER_VALUE_STRING);
        strcpy(params->data[params->count].name, token);

        /* Take second token as the parameter value. */
        token = strtok(NULL, " #\n");
        strcpy(params->data[params->count++].value, token);
      }
    }
  }
}

/**
 * @brief Retrieve integer parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam Value of the parameter found
 *
 */

void parser_get_param_int(struct swift_params *params, char *name,
                          int *retParam) {

  char str[128];

  for (int i = 0; i < params->count; i++) {

    /*strcmp returns 0 if both strings are the same.*/
    if (!strcmp(name, params->data[i].name)) {

      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%d%s", retParam, str) != 1) {
        error(
            "Tried parsing int '%s' but found '%s' with illegal integer "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      return;
    }
  }

  message("Cannot find '%s' in the structure.", name);
}

/**
 * @brief Retrieve float parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam Value of the parameter found
 *
 */

void parser_get_param_float(struct swift_params *params, char *name,
                            float *retParam) {

  char str[128];

  for (int i = 0; i < params->count; i++) {

    /*strcmp returns 0 if both strings are the same.*/
    if (!strcmp(name, params->data[i].name)) {

      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%f%s", retParam, str) != 1) {
        error(
            "Tried parsing float '%s' but found '%s' with illegal float "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      return;
    }
  }

  message("Cannot find '%s' in the structure.", name);
}

/**
 * @brief Retrieve double parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam Value of the parameter found
 *
 */

void parser_get_param_double(struct swift_params *params, char *name,
                             double *retParam) {

  char str[128];

  for (int i = 0; i < params->count; i++) {

    /*strcmp returns 0 if both strings are the same.*/
    if (!strcmp(name, params->data[i].name)) {

      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%lf", retParam) != 1) {
        error(
            "Tried parsing double '%s' but found '%s' with illegal double "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      return;
    }
  }

  message("Cannot find '%s' in the structure.", name);
}

/**
 * @brief Retrieve string parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam Value of the parameter found
 *
 */

void parser_get_param_string(struct swift_params *params, char *name,
                             char *retParam) {

  for (int i = 0; i < params->count; i++) {

    /*strcmp returns 0 if both strings are the same.*/
    if (!strcmp(name, params->data[i].name)) {
      strcpy(retParam, params->data[i].value);
      return;
    }
  }
}

/**
 * @brief Prints the contents of the parameter structure.
 *
 * @param params Structure that holds the parameters
 *
 */

void parser_print_params(struct swift_params *params) {

  printf("\n--------------------------\n");
  printf("|  SWIFT Parameter File  |\n");
  printf("--------------------------\n");

  for (int i = 0; i < params->count; i++) {
    printf("Parameter name: %s\n", params->data[i].name);
    printf("Parameter value: %s\n", params->data[i].value);
  }
}
