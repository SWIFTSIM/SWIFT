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
#include "tools.h"
#include "version.h"

#define PARSER_COMMENT_STRING "#"
#define PARSER_COMMENT_CHAR '#'
#define PARSER_VALUE_CHAR ':'
#define PARSER_VALUE_STRING ":"
#define PARSER_START_OF_FILE "---"
#define PARSER_END_OF_FILE "..."

#define CHUNK 10

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
  int count = CHUNK;
  strings = (char **)malloc(count * sizeof(char *));

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
          count += CHUNK;
          strings = (char **)realloc(strings, count * sizeof(char *));
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
      strings = (char **)realloc(strings, count * sizeof(char *));
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
  char section[PARSER_MAX_LINE_SIZE];
  name[0] = '\0';
  value[0] = '\0';

  /* Name is part until second colon. */
  const char *p1 = strchr(namevalue, ':');
  if (p1 != NULL) {

    /* Section is first part until a colon. */
    memcpy(section, namevalue, p1 - namevalue);
    section[p1 - namevalue] = ':';
    section[p1 - namevalue + 1] = '\0';

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
      strcpy(params->data[i].value, trim_both(value));
      updated = 1;
    }
  }
  if (!updated) {
    /* Is this a new section? */
    int newsection = 1;
    for (int i = 0; i < params->sectionCount; i++) {
      if (strcmp(section, params->section[i].name) == 0) {
        newsection = 0;
        break;
      }
    }
    if (newsection) {
      strcpy(params->section[params->sectionCount].name, section);
      params->sectionCount++;
      if (params->sectionCount == PARSER_MAX_NO_OF_SECTIONS)
        error("Too many sections, current maximum is %d.",
              params->sectionCount);
    }

    strcpy(params->data[params->paramCount].name, name);
    strcpy(params->data[params->paramCount].value, value);
    params->data[params->paramCount].used = 0;
    params->data[params->paramCount].is_default = 0;
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
  if (line[0] != PARSER_COMMENT_CHAR) {
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

        /* Trim trailing space before parsing line for a value. */
        char no_space_line[PARSER_MAX_LINE_SIZE];
        strcpy(no_space_line, trim_trailing(trim_line));

        parse_value(no_space_line, params);
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
  if (!inSection && line[0] == ' ') {
    error(
        "Invalid line:%d '%s', standalone parameter defined with incorrect "
        "indentation.",
        lineNumber, line);
  }

  /* Check that it is a parameter inside a section.*/
  if (line[0] == ' ' || line[0] == '\t') {
    parse_section_param(line, &isFirstParam, section, params);
  } else {
    /* It is the start of a new section or standalone parameter.
     * Take first token as the parameter name. */
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
      params->data[params->paramCount].is_default = 0;
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
  params->data[params->paramCount].is_default = 0;
  if (params->paramCount == PARSER_MAX_NO_OF_PARAMS - 1) {
    error("Maximal number of parameters in parameter file reached. Aborting !");
  } else {
    params->paramCount++;
  }
}

// Retrieve parameter value from structure. TYPE is the data type, float, int
// etc. FMT the format required for that data type, i.e. %f, %d etc. and DESC
// a one word description of the type, "float", "int" etc.
#define PARSER_GET_VALUE(TYPE, FMT, DESC)                                      \
  static int get_param_##TYPE(struct swift_params *params, const char *name,   \
                              TYPE *def, TYPE *result) {                       \
    char str[PARSER_MAX_LINE_SIZE];                                            \
    for (int i = 0; i < params->paramCount; i++) {                             \
      if (strcmp(name, params->data[i].name) == 0) {                           \
        /* Check that exactly one number is parsed, capture junk. */           \
        if (sscanf(params->data[i].value, " " FMT "%s ", result, str) != 1) {  \
          error("Tried parsing " DESC                                          \
                " '%s' but found '%s' with "                                   \
                "illegal trailing characters '%s'.",                           \
                params->data[i].name, params->data[i].value, str);             \
        }                                                                      \
        /* Ensure same behavior if called multiple times for same parameter */ \
        if (params->data[i].is_default && def == NULL)                         \
          error(                                                               \
              "Tried parsing %s again but cannot parse a default "             \
              "parameter as mandatory",                                        \
              name);                                                           \
        if (params->data[i].is_default && *def != *result)                     \
          error(                                                               \
              "Tried parsing %s again but cannot parse a parameter with "      \
              "two different default value (" FMT "!=" FMT ")",                \
              name, *def, *result);                                            \
        /* This parameter has been used */                                     \
        params->data[i].used = 1;                                              \
        return 1;                                                              \
      }                                                                        \
    }                                                                          \
    if (def == NULL)                                                           \
      error("Cannot find '%s' in the structure, in file '%s'.", name,          \
            params->fileName);                                                 \
    return 0;                                                                  \
  }

// Set a parameter to a value and save for dumping.
#define PARSER_SAVE_VALUE(PREFIX, TYPE, FMT)                      \
  static void save_param_##PREFIX(struct swift_params *params,    \
                                  const char *name, TYPE value) { \
    char str[PARSER_MAX_LINE_SIZE];                               \
    sprintf(str, "%s:" FMT, name, value);                        \
    parser_set_param(params, str);                                \
    params->data[params->paramCount - 1].used = 1;                \
    params->data[params->paramCount - 1].is_default = 0;          \
  }

/* Instantiations. */
PARSER_GET_VALUE(char, "%c", "char");
PARSER_GET_VALUE(int, "%d", "int");
PARSER_GET_VALUE(float, "%f", "float");
PARSER_GET_VALUE(double, "%lf", "double");
PARSER_SAVE_VALUE(char, char, "%c");
PARSER_SAVE_VALUE(int, int, "%d");
PARSER_SAVE_VALUE(float, float, "%g");
PARSER_SAVE_VALUE(double, double, "%g");
PARSER_SAVE_VALUE(string, const char *, "%s");

/**
 * @brief Retrieve integer parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
int parser_get_param_int(struct swift_params *params, const char *name) {
  int result = 0;
  get_param_int(params, name, NULL, &result);
  return result;
}

/**
 * @brief Retrieve char parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
char parser_get_param_char(struct swift_params *params, const char *name) {
  char result = 0;
  get_param_char(params, name, NULL, &result);
  return result;
}

/**
 * @brief Retrieve float parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
float parser_get_param_float(struct swift_params *params, const char *name) {
  float result = 0;
  get_param_float(params, name, NULL, &result);
  return result;
}

/**
 * @brief Retrieve double parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @return Value of the parameter found
 */
double parser_get_param_double(struct swift_params *params, const char *name) {
  double result = 0;
  get_param_double(params, name, NULL, &result);
  return result;
}

/**
 * @brief Retrieve string parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param retParam (return) Value of the parameter found
 */
void parser_get_param_string(struct swift_params *params, const char *name,
                             char *retParam) {

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      if (params->data[i].is_default)
        error(
            "Tried parsing %s again but cannot parse a "
            "default parameter as mandatory",
            name);
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
int parser_get_opt_param_int(struct swift_params *params, const char *name,
                             int def) {
  int result = 0;
  if (get_param_int(params, name, &def, &result)) return result;
  save_param_int(params, name, def);
  params->data[params->paramCount - 1].is_default = 1;
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
char parser_get_opt_param_char(struct swift_params *params, const char *name,
                               char def) {
  char result = 0;
  if (get_param_char(params, name, &def, &result)) return result;
  save_param_char(params, name, def);
  params->data[params->paramCount - 1].is_default = 1;
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
float parser_get_opt_param_float(struct swift_params *params, const char *name,
                                 float def) {
  float result = 0;
  if (get_param_float(params, name, &def, &result)) return result;
  save_param_float(params, name, def);
  params->data[params->paramCount - 1].is_default = 1;
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
  double result = 0;
  if (get_param_double(params, name, &def, &result)) return result;
  save_param_double(params, name, def);
  params->data[params->paramCount - 1].is_default = 1;
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
void parser_get_opt_param_string(struct swift_params *params, const char *name,
                                 char *retParam, const char *def) {

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      strcpy(retParam, params->data[i].value);

      /* Ensure same behavior if called multiple times for same parameter */
      if (params->data[i].is_default && !strcmp(def, retParam))
        error(
            "Tried parsing %s again but cannot parse a parameter with "
            "two different default value ('%s' != '%s')",
            name, def, retParam);
      /* this parameter has been used */
      params->data[i].used = 1;
      return;
    }
  }
  save_param_string(params, name, def);
  params->data[params->paramCount - 1].is_default = 1;
  strcpy(retParam, def);
}

/* Macro defining functions that get primitive types as simple one-line YAML
 * arrays, that is SEC: [v1,v2,v3...] format, with the extension that the []
 * are optional. TYPE is the data type, float etc. FMT a format to parse a
 * single value, so "%f" for a float and DESC the type description
 * i.e. "float".
 */
#define PARSER_GET_ARRAY(TYPE, FMT, DESC)                             \
  static int get_param_##TYPE##_array(struct swift_params *params,    \
                                      const char *name, int required, \
                                      int nval, TYPE *values) {       \
    char str[PARSER_MAX_LINE_SIZE];                                   \
    char cpy[PARSER_MAX_LINE_SIZE];                                   \
                                                                      \
    for (int i = 0; i < params->paramCount; i++) {                    \
      if (!strcmp(name, params->data[i].name)) {                      \
        if (params->data[i].is_default && required)                   \
          error(                                                      \
              "Tried parsing %s again but cannot parse a default "    \
              "parameter as mandatory",                               \
              name);                                                  \
        char *cp = cpy;                                               \
        strcpy(cp, params->data[i].value);                            \
        cp = trim_both(cp);                                           \
                                                                      \
        /* Strip off [], if present. */                               \
        if (cp[0] == '[') cp++;                                       \
        int l = strlen(cp);                                           \
        if (cp[l - 1] == ']') cp[l - 1] = '\0';                       \
        cp = trim_both(cp);                                           \
                                                                      \
        /* Format that captures spaces and trailing junk. */          \
        char fmt[20];                                                 \
        sprintf(fmt, " %s%%s ", FMT);                                 \
                                                                      \
        /* Parse out values which should now be "v, v, v" with        \
         * internal     whitespace variations. */                     \
        char *p = strtok(cp, ",");                                    \
        for (int k = 0; k < nval; k++) {                              \
          if (p != NULL) {                                            \
            TYPE tmp_value;                                           \
            if (sscanf(p, fmt, &tmp_value, str) != 1) {               \
              error("Tried parsing " DESC                             \
                    " '%s' but found '%s' with "                      \
                    "illegal " DESC " characters '%s'.",              \
                    name, p, str);                                    \
            }                                                         \
            if (params->data[i].is_default && tmp_value != values[k]) \
              error(                                                  \
                  "Tried parsing %s again but cannot parse a "        \
                  "parameter with two different default value "       \
                  "(" FMT "!=" FMT ")",                               \
                  name, tmp_value, values[k]);                        \
            values[k] = tmp_value;                                    \
          } else {                                                    \
            error(                                                    \
                "Array '%s' with value '%s' has too few values, "     \
                "expected %d",                                        \
                name, params->data[i].value, nval);                   \
          }                                                           \
          if (k < nval - 1) p = strtok(NULL, ",");                    \
        }                                                             \
        params->data[i].used = 1;                                     \
        return 1;                                                     \
      }                                                               \
    }                                                                 \
    if (required)                                                     \
      error("Cannot find '%s' in the structure, in file '%s'.", name, \
            params->fileName);                                        \
    return 0;                                                         \
  }

// Set values of a default parameter so they will be saved correctly.
#define PARSER_SAVE_ARRAY(TYPE, FMT)                                           \
  static int save_param_##TYPE##_array(                                        \
      struct swift_params *params, const char *name, int nval, TYPE *values) { \
    /* Save values against the parameter. */                                   \
    char str[PARSER_MAX_LINE_SIZE];                                            \
    int k = sprintf(str, "%s: [", name);                                       \
    for (int i = 0; i < nval - 1; i++)                                         \
      k += sprintf(&str[k], FMT ", ", values[i]);                              \
    sprintf(&str[k], FMT "]", values[nval - 1]);                               \
    parser_set_param(params, str);                                             \
    params->data[params->paramCount - 1].used = 1;                             \
    params->data[params->paramCount - 1].is_default = 0;                       \
    return 0;                                                                  \
  }

/* Instantiations. */
PARSER_GET_ARRAY(char, "%c", "char");
PARSER_GET_ARRAY(int, "%d", "int");
PARSER_GET_ARRAY(float, "%f", "float");
PARSER_GET_ARRAY(double, "%lf", "double");
PARSER_SAVE_ARRAY(char, "%c");
PARSER_SAVE_ARRAY(int, "%d");
PARSER_SAVE_ARRAY(float, "%g");
PARSER_SAVE_ARRAY(double, "%g");

/**
 * @brief Retrieve char array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 */
void parser_get_param_char_array(struct swift_params *params, const char *name,
                                 int nval, char *values) {
  get_param_char_array(params, name, 1, nval, values);
}

/**
 * @brief Retrieve optional char array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals. If the
 *               parameter is not found these values will be returned
 *               unmodified, so should be set to the default values.
 * @return whether the parameter has been found.
 */
int parser_get_opt_param_char_array(struct swift_params *params,
                                    const char *name, int nval, char *values) {
  if (get_param_char_array(params, name, 0, nval, values) != 1) {
    save_param_char_array(params, name, nval, values);
    params->data[params->paramCount - 1].is_default = 1;
    return 0;
  }
  return 1;
}

/**
 * @brief Retrieve int array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 */
void parser_get_param_int_array(struct swift_params *params, const char *name,
                                int nval, int *values) {
  get_param_int_array(params, name, 1, nval, values);
}

/**
 * @brief Retrieve optional int array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals. If the
 *               parameter is not found these values will be returned
 *               unmodified, so should be set to the default values.
 * @return whether the parameter has been found.
 */
int parser_get_opt_param_int_array(struct swift_params *params,
                                   const char *name, int nval, int *values) {
  if (get_param_int_array(params, name, 0, nval, values) != 1) {
    save_param_int_array(params, name, nval, values);
    params->data[params->paramCount - 1].is_default = 1;
    return 0;
  }
  return 1;
}

/**
 * @brief Retrieve float array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 */
void parser_get_param_float_array(struct swift_params *params, const char *name,
                                  int nval, float *values) {
  get_param_float_array(params, name, 1, nval, values);
}

/**
 * @brief Retrieve optional float array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals. If the
 *               parameter is not found these values will be returned
 *               unmodified, so should be set to the default values.
 * @return whether the parameter has been found.
 */
int parser_get_opt_param_float_array(struct swift_params *params,
                                     const char *name, int nval,
                                     float *values) {
  if (get_param_float_array(params, name, 0, nval, values) != 1) {
    save_param_float_array(params, name, nval, values);
    params->data[params->paramCount - 1].is_default = 1;
    return 0;
  }
  return 1;
}

/**
 * @brief Retrieve double array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals.
 */
void parser_get_param_double_array(struct swift_params *params,
                                   const char *name, int nval, double *values) {
  get_param_double_array(params, name, 1, nval, values);
}

/**
 * @brief Retrieve optional double array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values expected.
 * @param values Values of the parameter found, of size at least nvals. If the
 *               parameter is not found these values will be returned
 *               unmodified, so should be set to the default values.
 * @return whether the parameter has been found.
 */
int parser_get_opt_param_double_array(struct swift_params *params,
                                      const char *name, int nval,
                                      double *values) {
  if (get_param_double_array(params, name, 0, nval, values) != 1) {
    save_param_double_array(params, name, nval, values);
    params->data[params->paramCount - 1].is_default = 1;
    return 0;
  }
  return 1;
}

/**
 * @brief Retrieve string array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param required whether the parameter is required or not.
 * @param nval number of values located.
 * @param values pointer to an array of [nval] pointers to the strings.
 *        Note this must be freed by a call to
 *        parser_free_param_string_array when no longer required.
 * @result whether the parameter was found or not. Note if required
 *        an error will be thrown.
 */
static int get_string_array(struct swift_params *params, const char *name,
                            int required, int *nval, char ***values) {

  char cpy[PARSER_MAX_LINE_SIZE];
  *nval = 0;

  for (int i = 0; i < params->paramCount; i++) {
    if (!strcmp(name, params->data[i].name)) {
      char *cp = cpy;
      strcpy(cp, params->data[i].value);
      cp = trim_both(cp);

      /* Strip off [], if present. */
      if (cp[0] == '[') cp++;
      int l = strlen(cp);
      if (cp[l - 1] == ']') cp[l - 1] = '\0';
      cp = trim_both(cp);

      *nval = parse_quoted_strings(cp, values);

      params->data[i].used = 1;
      return 1;
    }
  }
  if (required)
    error("Cannot find '%s' in the structure, in file '%s'.", name,
          params->fileName);
  return 0;
}

/**
 * @brief Retrieve string array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values located.
 * @param values pointer to an array of [nval] pointers to the strings.
 *        Note this must be freed by a call to
 *        parser_free_param_string_array when no longer required.
 */
void parser_get_param_string_array(struct swift_params *params,
                                   const char *name, int *nval,
                                   char ***values) {
  get_string_array(params, name, 1, nval, values);
}

/**
 * @brief Retrieve optional string array parameter from structure.
 *
 * @param params Structure that holds the parameters
 * @param name Name of the parameter to be found
 * @param nval number of values located.
 * @param values pointer to an array of [nval] pointers to the strings.
 *        Note this must be freed by a call to
 *        parser_free_param_string_array when no longer required.
 * @param ndef the number of default values.
 * @param def the default values as an array of pointers to strings.
 *        Note copied to values if used.
 * @result whether the parameter was found or not.
 */
int parser_get_opt_param_string_array(struct swift_params *params,
                                      const char *name, int *nval,
                                      char ***values, int ndef,
                                      const char *def[]) {
  if (get_string_array(params, name, 0, nval, values) == 1) return 1;

  /* Not found, so save the default values against the parameter. Look for
   * single quotes in value and use if not found, otherwise use double
   * quotes. We don't support having both in a string. */
  char cpy[PARSER_MAX_LINE_SIZE];
  int k = sprintf(cpy, "%s: [", name);
  int i = 0;
  for (i = 0; i < ndef - 1; i++) {
    if (strchr(def[i], '\'') == 0)
      k += sprintf(&cpy[k], "'%s', ", def[i]);
    else
      k += sprintf(&cpy[k], "\"%s\", ", def[i]);
  }
  if (strchr(def[i], '\'') == 0)
    sprintf(&cpy[k], "'%s']", def[i]);
  else
    sprintf(&cpy[k], "\"%s\"]", def[i]);
  parser_set_param(params, cpy);
  params->data[params->paramCount - 1].is_default = 1;
  params->data[params->paramCount - 1].used = 1;

  /* Now copy to output space. */
  char **strings;
  strings = (char **)malloc(ndef * sizeof(char *));
  for (int j = 0; j < ndef; j++) {
    strings[j] = (char *)malloc((strlen(def[j]) + 1) * sizeof(char));
    strcpy(strings[j], def[j]);
  }
  *values = strings;
  *nval = ndef;

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
  char param_name[PARSER_MAX_LINE_SIZE] = {0};
  char section[PARSER_MAX_LINE_SIZE] = {0};
  char *token;

  /* Start of file identifier in YAML. */
  fprintf(file, "%s\n\n", PARSER_START_OF_FILE);

  fprintf(file, "# SWIFT used parameter file\n");
  fprintf(file, "# Code version: %s\n", package_version());
  fprintf(file, "# git revision: %s\n", git_revision());
  fprintf(file, "# git branch: %s\n", git_branch());
  fprintf(file, "# git date: %s\n", git_date());

  /* Flags to track which parameters are written. */
  int *written = (int *)calloc(params->paramCount, sizeof(int));
  int nwritten = 0;

  /* Loop over all sections. These are not contiguous when storing optional
   * values. */
  for (int k = 0; k < params->sectionCount; k++) {
    int first = 1;

    /* Locate parameters in this section. */
    for (int i = 0; i < params->paramCount; i++) {
      if (!written[i] && ((write_used && params->data[i].used) ||
                          (!write_used && !params->data[i].used))) {

        /* Find section part of name, if have one. */
        if (strchr(params->data[i].name, PARSER_VALUE_CHAR)) {
          strcpy(param_name, params->data[i].name);
          token = strtok(param_name, PARSER_VALUE_STRING);
          strcpy(section, token);
          strcat(section, PARSER_VALUE_STRING);

          /* If in our section name print it to the file. */
          if (strcmp(section, params->section[k].name) == 0) {
            if (first) {
              fprintf(file, "\n%s\n", section);
              first = 0;
            }

            /* Remove white space from parameter name and write it to the
             * file. */
            token = trim_both(strtok(NULL, "#\n"));
            fprintf(file, "  %s%c %s\n", token, PARSER_VALUE_CHAR,
                    params->data[i].value);
            written[i] = 1;
            nwritten++;
          }
        }
      }
    }
  }

  /* Write out any parameters outside of sections. */
  if (nwritten < params->paramCount) {
    for (int i = 0; i < params->paramCount; i++) {
      if (!written[i] && ((write_used && params->data[i].used) ||
                          (!write_used && !params->data[i].used))) {
        fprintf(file, "\n%s%c %s\n", params->data[i].name, PARSER_VALUE_CHAR,
                params->data[i].value);
      }
    }
  }

  /* End of file identifier in YAML. */
  fprintf(file, "%s\n", PARSER_END_OF_FILE);

  free(written);
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
