/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
 *               2017-2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "parser.h"

/* Local headers. */
#include "common_io.h"
#include "error.h"
#include "restart.h"

#define PARSER_COMMENT_STRING "#"
#define PARSER_COMMENT_CHAR '#'
#define PARSER_VALUE_CHAR ':'
#define PARSER_VALUE_STRING ":"
#define PARSER_START_OF_FILE "---"
#define PARSER_END_OF_FILE "..."

/* Private functions. */
static int is_empty(const char *str);
static int count_indentation(const char *str);
static void parse_line(char *line, struct swift_params *params);
static void parse_value(char *line, struct swift_params *params);
static void parse_section_param(char *line, int *isFirstParam,
                                char *sectionName, struct swift_params *params);
static void find_duplicate_params(const struct swift_params *params,
                                  const char *param_name);
static void find_duplicate_section(const struct swift_params *params,
                                   const char *section_name);
static int lineNumber = 0;

/**
 * @brief trim leading white space from a string.
 *
 * @param s the string.
 * @result the result.
 */
static char *trim_leading(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  while(isspace(*s)) s++;
  return s;
}

/**
 * @brief trim trailing white space from a string.
 *
 * @param s the string.
 * @result the result.
 */
static char *trim_trailing(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  char *end = s + strlen(s) - 1;
  while(isspace(*end)) end--;
  *(end+1) = '\0';
  return s;
}

/**
 * @brief trim leading and trailing white space from a string.
 *
 * @param s the string.
 * @result the result.
 */
static char *trim_both(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  return trim_trailing(trim_leading(s));
}

/**
 * @brief parse a YAML list of strings returning a set of pointers to
 *        the strings.
 *
 * It is assumed that the [] have been removed (also no lists in lists)
 * words are separated by commas and the strings may or may not be quoted.
 * So lines like:
 *
 *    'xyz', 'ABC', "ab'c", "de:f", "g,hi", "zzz", Hello World, again
 *
 * Are supported as expected.
 *
 * @param line the line to parse.
 * @param result array of pointers to the strings.
 * @return the number of strings
 */
static int parse_quoted_strings(const char *line, char ***result) {

    char word[PARSER_MAX_LINE_SIZE];
    int nchar = 0;
    int nwords = 0;
    char quote = '\0';

    /* Preallocate a number of pointers. */
    char **strings;
    strings = (char **)malloc( 10 * sizeof(char *));
    int count = 10;

    word[0] = '\0';
    for (unsigned int i = 0; i < strlen(line); i++) {
      char c = line[i];
      if (c == '"' || c == '\'') {
        if (c == quote) {
          quote = '\0';
        } else if (!quote) {
          quote = c;
        } else {
          word[nchar++] = c;
        }
      } else if (c == ',') {
        if (!quote) {

          /* Save word. */
          word[nchar++] = '\0';
          if (count <= nwords) {
            count += 10;
            strings = (char **)realloc(strings, count * sizeof(char));
          }
          strings[nwords] = (char *)malloc((strlen(word) + 1) * sizeof(char));
          strcpy(strings[nwords], trim_both(word));
          nwords++;

          /* Ready for next. */
          nchar = 0;
          word[0] = '\0';

        } else {
          word[nchar++] = c;
        }
      } else {
        word[nchar++] = c;
      }
    }

    /* Keep unfinished words. */
    if (nchar > 0) {
      word[nchar] = '\0';
      if (count <= nwords) {
          count += 1;
          strings = (char **)realloc(strings, count * sizeof(char));
      }
      strings[nwords] = (char *)malloc((strlen(word) + 1) * sizeof(char));
      strcpy(strings[nwords], trim_both(word));
      nwords++;
    }

    *result = strings;
    return nwords;
}

/**
 * @brief Initialize the parser structure.
 *
 * @param file_name Name of file to be read
 * @param params Structure to be populated from file
 */
void parser_init(const char *file_name, struct swift_params *params) {
  params->paramCount = 0;
  params->sectionCount = 0;
  strcpy(params->fileName, file_name);
}

/**
 * @brief Reads an input file and stores each parameter in a structure.
 *
 * @param file_name Name of file to be read
 * @param params Structure to be populated from file
 */
void parser_read_file(const char *file_name, struct swift_params *params) {
  /* Open file for reading */
  FILE *file = fopen(file_name, "r");

  /* Line to parsed. */
  char line[PARSER_MAX_LINE_SIZE];

  /* Initialise parameter count. */
  parser_init(file_name, params);

  /* Check if parameter file exits. */
  if (file == NULL) {
    error("Error opening parameter file: %s", file_name);
  }

  /* Read until the end of the file is reached.*/
  while (!feof(file)) {
    if (fgets(line, PARSER_MAX_LINE_SIZE, file) != NULL) {
      lineNumber++;
      parse_line(line, params);
    }
  }

  fclose(file);
}

/**
 * @brief Set or update a parameter using a compressed format.
 *
 * The compressed format allows a value to be given as a single
 * string and has the format "section:parameter:value", with all
 * names as would be given in the parameter file.
 *
 * @param params Structure that holds the parameters.
 * @param namevalue the parameter name and value as described.
 */
void parser_set_param(struct swift_params *params, const char *namevalue) {

  /* Get the various parts. */
  char name[PARSER_MAX_LINE_SIZE];
  char value[PARSER_MAX_LINE_SIZE];
  name[0] = '\0';
  value[0] = '\0';

  /* Name is part until second colon. */
  const char *p1 = strchr(namevalue, ':');
  if (p1 != NULL) {
    const char *p2 = strchr(p1 + 1, ':');
    if (p2 != NULL) {
      memcpy(name, namevalue, p2 - namevalue);
      name[p2 - namevalue] = '\0';

      /* Value is rest after second colon. */
      p2++;
      strcpy(value, p2);
    }
  }

  /* Sanity check. */
  if (strlen(name) == 0 || strlen(value) == 0 || strchr(value, ':') != NULL)
    error(
        "Cannot parse compressed parameter string: '%s', check syntax "
        "should be section:parameter:value",
        namevalue);

  /* And update or set. */
  int updated = 0;
  for (int i = 0; i < params->paramCount; i++) {
    if (strcmp(name, params->data[i].name) == 0) {
      message("Value of '%s' changed from '%s' to '%s'", params->data[i].name,
              params->data[i].value, value);
      strcpy(params->data[i].value, value);
      updated = 1;
    }
  }
  if (!updated) {
    strcpy(params->data[params->paramCount].name, name);
    strcpy(params->data[params->paramCount].value, value);
    params->data[params->paramCount].used = 0;
    params->paramCount++;
    if (params->paramCount == PARSER_MAX_NO_OF_PARAMS)
      error("Too many parameters, current maximum is %d.", params->paramCount);
  }
}

/**
 * @brief Counts the number of white spaces that prefix a string.
 *
 * @param str String to be checked
 *
 * @return Number of white spaces prefixing str
 */

static int count_indentation(const char *str) {
  int count = 0;

  /* Check if the line contains the character */
  while (*(++str) == ' ') {
    count++;
  }
  return count;
}

/**
 * @brief Checks if a string is empty.
 *
 * @param str String to be checked
 *
 * @return Returns 1 if str is empty, 0 otherwise
 */

static int is_empty(const char *str) {
  int retParam = 1;
  while (*str != '\0') {
    if (!isspace(*str)) {
      retParam = 0;
      break;
    }
    str++;
  }

  return retParam;
}

/**
 * @brief Look for duplicate parameters.
 *
 * @param params Structure that holds the parameters
 * @param param_name Name of parameter to be searched for
 */

static void find_duplicate_params(const struct swift_params *params,
                                  const char *param_name) {
  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(param_name, params->data[i].name)) {
      error("Invalid line:%d '%s', parameter is a duplicate.", lineNumber,
            param_name);
    }
  }
}

/**
 * @brief Look for duplicate sections.
 *
 * @param params Structure that holds the parameters
 * @param section_name Name of section to be searched for
 */

static void find_duplicate_section(const struct swift_params *params,
                                   const char *section_name) {
  for (int i = 0; i < params->sectionCount; i++) {
    if (!strcmp(section_name, params->section[i].name)) {
      error("Invalid line:%d '%s', section is a duplicate.", lineNumber,
            section_name);
    }
  }
}
/**
 * @brief Parses a line from a file and stores any parameters in a structure.
 *
 * @param line Line to be parsed.
 * @param params Structure to be populated from file.
 */

static void parse_line(char *line, struct swift_params *params) {
  /* Parse line if it doesn't begin with a comment. */
  if (*line != PARSER_COMMENT_CHAR) {
    char trim_line[PARSER_MAX_LINE_SIZE];
    char tmp_str[PARSER_MAX_LINE_SIZE];
    char *token;

    /* Remove comments at the end of a line. */
    token = strtok(line, PARSER_COMMENT_STRING);
    strcpy(tmp_str, token);

    /* Check if the line is just white space. */
    if (!is_empty(tmp_str)) {
      /* Trim '\n' characters from string. */
      token = strtok(tmp_str, "\n");
      strcpy(trim_line, token);

      /* Check if the line contains a value and parse it. */
      if (strchr(trim_line, PARSER_VALUE_CHAR)) {
        parse_value(trim_line, params);
      }
      /* Check for invalid lines,not including the start and end of file. */
      else if (strcmp(trim_line, PARSER_START_OF_FILE) &&
               strcmp(trim_line, PARSER_END_OF_FILE)) {
        error("Invalid line:%d '%s'.", lineNumber, trim_line);
      }
    }
  }
}

/**
 * @brief Performs error checking and stores a parameter in a structure.
 *
 * @param line Line containing the parameter
 * @param params Structure to be written to
 *
 */

static void parse_value(char *line, struct swift_params *params) {
  static int inSection = 0;
  static char section[PARSER_MAX_LINE_SIZE]; /* Keeps track of current section
                                                name. */
  static int isFirstParam = 1;
  char tmpStr[PARSER_MAX_LINE_SIZE];
  char tmpSectionName[PARSER_MAX_LINE_SIZE];

  char *token;

  /* Check that standalone parameters have correct indentation. */
  if (!inSection && *line == ' ') {
    error(
        "Invalid line:%d '%s', standalone parameter defined with incorrect "
        "indentation.",
        lineNumber, line);
  }

  /* Check that it is a parameter inside a section.*/
  if (*line == ' ' || *line == '\t') {
    parse_section_param(line, &isFirstParam, section, params);
  } else { /*Else it is the start of a new section or standalone parameter. */
    /* Take first token as the parameter name. */
    token = strtok(line, ":\t");
    strcpy(tmpStr, trim_trailing(token));

    /* Take second token as the parameter value. */
    token = trim_both(strtok(NULL, "#\n"));

    /* If second token is NULL or empty then the line must be a section
     * heading. */
    if (token == NULL || strlen(token) == 0) {
      strcpy(tmpSectionName, tmpStr);
      strcat(tmpSectionName, PARSER_VALUE_STRING);

      /* Check for duplicate section name. */
      find_duplicate_section(params, tmpSectionName);

      /* Check for duplicate standalone parameter name used as a section name.
       */
      find_duplicate_params(params, tmpStr);

      strcpy(section, tmpSectionName);
      strcpy(params->section[params->sectionCount].name, tmpSectionName);
      if (params->sectionCount == PARSER_MAX_NO_OF_SECTIONS - 1) {
        error(
            "Maximal number of sections in parameter file reached. Aborting !");
      } else {
        params->sectionCount++;
      }
      inSection = 1;
      isFirstParam = 1;
    } else {
      /* Create string with standalone parameter name appended with ":" to aid
       * duplicate search as section names are stored with ":" at the end.*/
      strcpy(tmpSectionName, tmpStr);
      strcat(tmpSectionName, PARSER_VALUE_STRING);

      /* Check for duplicate parameter name. */
      find_duplicate_params(params, tmpStr);

      /* Check for duplicate section name used as standalone parameter name. */
      find_duplicate_section(params, tmpSectionName);

      /* Must be a standalone parameter so no need to prefix name with a
       * section. */
      strcpy(params->data[params->paramCount].name, tmpStr);
      strcpy(params->data[params->paramCount].value, token);
      params->data[params->paramCount].used = 0;
      if (params->paramCount == PARSER_MAX_NO_OF_PARAMS - 1) {
        error(
            "Maximal number of parameters in parameter file reached. Aborting "
            "!");
      } else {
        params->paramCount++;
      }
      inSection = 0;
      isFirstParam = 1;
    }
  }
}

/**
 * @brief Parses a parameter that appears in a section and stores it in a
 *structure.
 *
 * @param line Line containing the parameter
 * @param isFirstParam Shows if the first parameter of a section has been found
 * @param sectionName String containing the current section name
 * @param params Structure to be written to
 *
 */

static void parse_section_param(char *line, int *isFirstParam,
                                char *sectionName,
                                struct swift_params *params) {
  static int sectionIndent = 0;
  char tmpStr[PARSER_MAX_LINE_SIZE];
  char paramName[PARSER_MAX_LINE_SIZE];
  char *token;

  /* Count indentation of each parameter and check that it
   * is consistent with the first parameter in the section. */
  if (*isFirstParam) {
    sectionIndent = count_indentation(line);
    *isFirstParam = 0;
  } else if (count_indentation(line) != sectionIndent) {
    error("Invalid line:%d '%s', parameter has incorrect indentation.",
          lineNumber, line);
  }

  /* Take first token as the parameter name and trim leading white space. */
  token = trim_both(strtok(line, ":\t"));
  strcpy(tmpStr, token);

  /* Take second token as the parameter value. */
  token = trim_both(strtok(NULL, "#\n"));

  /* Prefix the parameter name with its section name and
   * copy it into the parameter structure. */
  strcpy(paramName, sectionName);
  strcat(paramName, tmpStr);

  /* Check for duplicate parameter name. */
  find_duplicate_params(params, paramName);

  strcpy(params->data[params->paramCount].name, paramName);
  strcpy(params->data[params->paramCount].value, token);
  params->data[params->paramCount].used = 0;
  if (params->paramCount == PARSER_MAX_NO_OF_PARAMS - 1) {
    error("Maximal number of parameters in parameter file reached. Aborting !");
  } else {
    params->paramCount++;
  }
}

/**
 * @brief Retrieve integer parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
int parser_get_param_int(struct swift_params *params, const char *name) {

  char str[PARSER_MAX_LINE_SIZE];
  int retParam = 0;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%d%s", &retParam, str) != 1) {
        error(
            "Tried parsing int '%s' but found '%s' with illegal integer "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  error("Cannot find '%s' in the structure, in file '%s'.", name,
        params->fileName);
  return 0;
}

/**
 * @brief Retrieve char parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
char parser_get_param_char(struct swift_params *params,
                           const char *name) {

  char str[PARSER_MAX_LINE_SIZE];
  char retParam = 0;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%c%s", &retParam, str) != 1) {
        error(
            "Tried parsing char '%s' but found '%s' with illegal char "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  error("Cannot find '%s' in the structure, in file '%s'.", name,
        params->fileName);
  return 0;
}

/**
 * @brief Retrieve float parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
float parser_get_param_float(struct swift_params *params,
                             const char *name) {

  char str[PARSER_MAX_LINE_SIZE];
  float retParam = 0.f;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%f%s", &retParam, str) != 1) {
        error(
            "Tried parsing float '%s' but found '%s' with illegal float "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  error("Cannot find '%s' in the structure, in file '%s'.", name,
        params->fileName);
  return 0.f;
}

/**
 * @brief Retrieve double parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
double parser_get_param_double(struct swift_params *params,
                               const char *name) {

  char str[PARSER_MAX_LINE_SIZE];
  double retParam = 0.;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%lf%s", &retParam, str) != 1) {
        error(
            "Tried parsing double '%s' but found '%s' with illegal double "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  error("Cannot find '%s' in the structure, in file '%s'.", name,
        params->fileName);
  return 0.;
}

/**
 * @brief Retrieve string parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam (return) Value of the parameter found
 */
void parser_get_param_string(struct swift_params *params,
                             const char *name, char *retParam) {

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      strcpy(retParam, params->data[i].value);
      /* this parameter has been used */
      params->data[i].used = 1;
      return;
    }
  }

  error("Cannot find '%s' in the structure.", name);
}

/**
 * @brief Retrieve optional integer parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param def Default value of the parameter of not found.
 * @return Value of the parameter found
 */
int parser_get_opt_param_int(struct swift_params *params,
                             const char *name, int def) {

  char str[PARSER_MAX_LINE_SIZE];
  int retParam = 0;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%d%s", &retParam, str) != 1) {
        error(
            "Tried parsing int '%s' but found '%s' with illegal integer "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  /* Generate string for new parameter */
  sprintf(str, "%s: %i", name, def);

  /* Add it to params */
  parser_set_param(params, str);

  /* Set parameter as used */
  params->data[params->paramCount-1].used = 1;

  return def;
}

/**
 * @brief Retrieve optional char parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param def Default value of the parameter of not found.
 * @return Value of the parameter found
 */
char parser_get_opt_param_char(struct swift_params *params,
                               const char *name, char def) {

  char str[PARSER_MAX_LINE_SIZE];
  char retParam = 0;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%c%s", &retParam, str) != 1) {
        error(
            "Tried parsing char '%s' but found '%s' with illegal char "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  /* Generate string for new parameter */
  sprintf(str, "%s: %c", name, def);

  /* Add it to params */
  parser_set_param(params, str);

  /* Set parameter as used */
  params->data[params->paramCount-1].used = 1;

  return def;
}

/**
 * @brief Retrieve optional float parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param def Default value of the parameter of not found.
 * @return Value of the parameter found
 */
float parser_get_opt_param_float(struct swift_params *params,
                                 const char *name, float def) {

  char str[PARSER_MAX_LINE_SIZE];
  float retParam = 0.f;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%f%s", &retParam, str) != 1) {
        error(
            "Tried parsing float '%s' but found '%s' with illegal float "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  /* Generate string for new parameter */
  sprintf(str, "%s: %f", name, def);

  /* Add it to params */
  parser_set_param(params, str);

  /* Set parameter as used */
  params->data[params->paramCount-1].used = 1;

  return def;
}

/**
 * @brief Retrieve optional double parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param def Default value of the parameter of not found.
 * @return Value of the parameter found
 */
double parser_get_opt_param_double(struct swift_params *params,
                                   const char *name, double def) {

  char str[PARSER_MAX_LINE_SIZE];
  double retParam = 0.;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      /* Check that exactly one number is parsed. */
      if (sscanf(params->data[i].value, "%lf%s", &retParam, str) != 1) {
        error(
            "Tried parsing double '%s' but found '%s' with illegal double "
            "characters '%s'.",
            params->data[i].name, params->data[i].value, str);
      }

      /* this parameter has been used */
      params->data[i].used = 1;

      return retParam;
    }
  }

  /* Generate string for new parameter */
  sprintf(str, "%s: %lf", name, def);

  /* Add it to params */
  parser_set_param(params, str);

  /* Set parameter as used */
  params->data[params->paramCount-1].used = 1;

  return def;
}

/**
 * @brief Retrieve string parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param def Default value of the parameter of not found.
 * @param retParam (return) Value of the parameter found
 */
void parser_get_opt_param_string(struct swift_params *params,
                                 const char *name, char *retParam,
                                 const char *def) {

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      strcpy(retParam, params->data[i].value);

      /* this parameter has been used */
      params->data[i].used = 1;

      return;
    }
  }

  /* Generate string for new parameter */
  char str[PARSER_MAX_LINE_SIZE];
  sprintf(str, "%s: %s", name, def);

  /* Add it to params */
  parser_set_param(params, str);

  /* Set parameter as used */
  params->data[params->paramCount-1].used = 1;

  strcpy(retParam, def);
}

/* Macro defining functions that get primitive types as simple one-line YAML
 * arrays, that is SEC: [v1,v2,v3...] format. TYPE is the data type, float
 * etc. FMT a format to parse a single value with the possibility of
 * extraneous characters after the value, so " %f%s " for a float (spaces can
 * be eaten so safest to have one at beginning and end) and DESC the type
 * description i.e. "float".
 */
#define PARSER_GET_ARRAY(TYPE, FMT, DESC)                               \
  int parser_get_param_##TYPE##_array(const struct swift_params *params,\
                                   const char *name, int required,      \
                                   int nval, TYPE *values) {            \
    char str[PARSER_MAX_LINE_SIZE];                                     \
    char cpy[PARSER_MAX_LINE_SIZE];                                     \
                                                                        \
    for (int i = 0; i < params->paramCount; i++) {                      \
      if (!strcmp(name, params->data[i].name)) {                        \
        char *cp = cpy;                                                 \
        strcpy(cp, params->data[i].value);                              \
        cp = trim_both(cp);                                             \
                                                                        \
        /* Strip off []. */                                             \
        if (cp[0] != '[')                                               \
          error("Array '%s' does not start with '['", name);            \
        cp++;                                                           \
        int l = strlen(cp);                                             \
        if (cp[l-1] != ']')                                             \
          error("Array '%s' does not end with ']'", name);              \
        cp[l-1] = '\0';                                                 \
        cp = trim_both(cp);                                             \
                                                                        \
        /* Parse out values which should now be "v, v, v" with          \
         * internal     whitespace variations. */                       \
        char *p = strtok(cp, ",");                                      \
        for (int k = 0; k < nval; k++) {                                \
          if (p != NULL) {                                              \
            if (sscanf(p, FMT, &values[k], str) != 1) {                 \
              error("Tried parsing " DESC " '%s' but found '%s' with "  \
                    "illegal " DESC " characters '%s'.", name, p, str); \
            }                                                           \
          } else {                                                      \
            error("Array '%s' with value '%s' has too few values, "     \
                  "expected %d", name, params->data[i].value, nval);    \
          }                                                             \
          if (k < nval - 1) p = strtok(NULL, ",");                      \
        }                                                               \
        return 1;                                                       \
      }                                                                 \
    }                                                                   \
    if (required)                                                       \
      error("Cannot find '%s' in the structure, in file '%s'.", name,   \
            params->fileName);                                          \
    return 0;                                                           \
  }

/**
 * @brief Retrieve int char parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether it is an error if the parameter has not been set
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 * @return whether the parameter has been found.
 */
PARSER_GET_ARRAY(char, " %c%s ", "int");

/**
 * @brief Retrieve int array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether it is an error if the parameter has not been set
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 * @return whether the parameter has been found.
 */
PARSER_GET_ARRAY(int, " %d%s ", "int");

/**
 * @brief Retrieve float array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether it is an error if the parameter has not been set
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 * @return whether the parameter has been found.
 */
PARSER_GET_ARRAY(float, " %f%s ", "float");

/**
 * @brief Retrieve double array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether it is an error if the parameter has not been set
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 * @return whether the parameter has been found.
 */
PARSER_GET_ARRAY(double, " %lf%s ", "double");

/**
 * @brief Retrieve string array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether it is an error if the parameter has not been set
 * @param nval number of values located.
 * @param values pointer to an array of [nval] pointers to the strings.
 *        Note this must be freed by a call to
 *        parser_free_param_string_array.
 * @return whether the parameter has been found.
 */
int parser_get_param_string_array(const struct swift_params *params,
                                  const char *name, int required,
                                  int *nval, char ***values) {

    char cpy[PARSER_MAX_LINE_SIZE];

    for (int i = 0; i < params->paramCount; i++) {
      if (!strcmp(name, params->data[i].name)) {
        char *cp = cpy;
        strcpy(cp, params->data[i].value);
        cp = trim_both(cp);

        /* Strip off []. */
        if (cp[0] != '[')
          error("Array '%s' does not start with '['", name);
        cp++;
        int l = strlen(cp);
        if (cp[l-1] != ']')
          error("Array '%s' does not end with ']'", name);
        cp[l-1] = '\0';
        cp = trim_both(cp);

        *nval = parse_quoted_strings(cp, values);

        return 1;
      }
    }
    if (required)
      error("Cannot find '%s' in the structure, in file '%s'.", name,
            params->fileName);
    return 0;
}

/**
 * @brief Free string array allocated by parser_get_param_string_array.
 *
 * @param nval number of strings returned.
 * @param values pointer to the returned values.
 */
void parser_free_param_string_array(int nval, char **values) {
  for (int i = 0; i < nval; i++) {
    free(values[i]);
  }
  free(values);
  return;
}

/**
 * @brief Prints the contents of the parameter structure.
 *
 * @param params Structure that holds the parameters
 */
void parser_print_params(const struct swift_params *params) {
  printf("\n--------------------------\n");
  printf("|  SWIFT Parameter File  |\n");
  printf("--------------------------\n");

  for (int i = 0; i < params->paramCount; i++) {
    printf("Parameter name: %s\n", params->data[i].name);
    printf("Parameter value: %s\n", params->data[i].value);
    printf("Parameter used: %i\n", params->data[i].used);
  }
}

/**
 * @brief Write the contents of the parameter structure to a file in YAML
 *format.
 *
 * @param params Structure that holds the parameters
 * @param file_name Name of file to be written
 * @param write_used Write used fields or unused fields.
 */
void parser_write_params_to_file(const struct swift_params *params,
                                 const char *file_name, int write_used) {
  FILE *file = fopen(file_name, "w");
  char section[PARSER_MAX_LINE_SIZE] = {0};
  char param_name[PARSER_MAX_LINE_SIZE] = {0};
  char *token;

  /* Start of file identifier in YAML. */
  fprintf(file, "%s\n", PARSER_START_OF_FILE);

  for (int i = 0; i < params->paramCount; i++) {
    if (write_used && !params->data[i].used) {
#ifdef SWIFT_DEBUG_CHECKS
      message("Parameter `%s` was not used. "
              "Only the parameter used are written.",
              params->data[i].name);
#endif
      continue;
    }
    else if (!write_used && params->data[i].used)
      continue;
    /* Check that the parameter name contains a section name. */
    if (strchr(params->data[i].name, PARSER_VALUE_CHAR)) {
      /* Copy the parameter name into a temporary string and find the section
       * name. */
      strcpy(param_name, params->data[i].name);
      token = strtok(param_name, PARSER_VALUE_STRING);

      /* If a new section name is found print it to the file. */
      if (strcmp(token, section)) {
        strcpy(section, token);
        fprintf(file, "\n%s%c\n", section, PARSER_VALUE_CHAR);
      }

      /* Remove white space from parameter name and write it to the file. */
      token = trim_both(strtok(NULL, "#\n"));

      fprintf(file, "  %s%c %s\n", token, PARSER_VALUE_CHAR,
              params->data[i].value);
    } else {
      fprintf(file, "\n%s%c %s\n", params->data[i].name, PARSER_VALUE_CHAR,
              params->data[i].value);
    }
  }

  /* End of file identifier in YAML. */
  fprintf(file, "%s\n", PARSER_END_OF_FILE);

  fclose(file);
}

#if defined(HAVE_HDF5)

/**
 * @brief Write the contents of the parameter structure to a hdf5 file
 *
 * @param params Structure that holds the parameters
 * @param grp HDF5 group
 * @param write_used Write used fields or unused fields.
 */
void parser_write_params_to_hdf5(const struct swift_params *params, hid_t grp,
                                 int write_used) {

  for (int i = 0; i < params->paramCount; i++) {
    if (write_used && !params->data[i].used)
      continue;
    else if (!write_used && params->data[i].used)
      continue;
    io_write_attribute_s(grp, params->data[i].name, params->data[i].value);
  }
}
#endif

/**
 * @brief Write a swift_params struct to the given FILE as a stream of bytes.
 *
 * @param params the struct
 * @param stream the file stream
 */
void parser_struct_dump(const struct swift_params *params, FILE *stream) {
  restart_write_blocks((void *)params, sizeof(struct swift_params), 1, stream,
                       "parameters", "parameters");
}

/**
 * @brief Restore a swift_params struct from the given FILE as a stream of
 * bytes.
 *
 * @param params the struct
 * @param stream the file stream
 */
void parser_struct_restore(const struct swift_params *params, FILE *stream) {
  restart_read_blocks((void *)params, sizeof(struct swift_params), 1, stream,
                      NULL, "parameters");
}
