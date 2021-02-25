/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_python_tools.h"
#include "logger_reader.h"
#include "logger_time.h"

#ifdef HAVE_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/* Include PyMethodDef */
#include <structmember.h>

static PyTypeObject PyParticleReader_Type;

/**
 * @brief Reader for a single particle type.
 */
typedef struct {
  PyObject_HEAD;

  /* Link to the reader */
  void *reader;

  /* Particle type */
  enum part_type type;
} PyParticleReader;

/**
 * @brief Main class for the reader.
 */
typedef struct {
  PyObject_HEAD;

  /* Add logger stuff here */
  struct logger_reader reader;

  /* Is the logger ready? */
  int ready;

  /* basename of the logfile */
  char *basename;

  /* Verbose level */
  int verbose;

  /* Reader for each type of particles */
  PyParticleReader *part_reader[swift_type_count];

} PyObjectReader;

/**
 * @brief Deallocate the memory of a particle reader.
 */
__attribute__((always_inline)) INLINE static void particle_reader_dealloc(
    PyParticleReader *self) {
  Py_TYPE(self)->tp_free((PyObject *)self);
}

/**
 * @brief Deallocate the memory.
 */
static void Reader_dealloc(PyObjectReader *self) {
  if (self->basename != NULL) {
    free(self->basename);
    self->basename = NULL;
  }
  Py_TYPE(self)->tp_free((PyObject *)self);
}

/**
 * @brief Allocate the memory of a particle reader.
 */
__attribute__((always_inline)) INLINE static PyObject *particle_reader_new(
    PyTypeObject *type, PyObject *args, PyObject *kwds) {
  PyParticleReader *self = (PyParticleReader *)type->tp_alloc(type, 0);
  return (PyObject *)self;
}

/**
 * @brief Allocate the memory.
 */
static PyObject *Reader_new(PyTypeObject *type, PyObject *args,
                            PyObject *kwds) {
  PyObjectReader *self;

  self = (PyObjectReader *)type->tp_alloc(type, 0);
  self->ready = 0;
  self->basename = NULL;

  for (int i = 0; i < swift_type_count; i++) {
    self->part_reader[i] = NULL;
  }

  return (PyObject *)self;
}
/**
 * @brief Initializer of the Reader.
 *
 * @param basename The basename of the logger file.
 * @param verbose The verbosity level.
 *
 * @return The Reader object.
 */
static int Reader_init(PyObjectReader *self, PyObject *args, PyObject *kwds) {

  /* input variables. */
  char *basename = NULL;
  int verbose = 0;

  /* List of keyword arguments. */
  static char *kwlist[] = {"basename", "verbose", NULL};

  /* parse the arguments. */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", kwlist, &basename,
                                   &verbose))
    return -1;

  /* Copy the arguments */
  self->verbose = verbose;

  size_t n_plus_null = strlen(basename) + 1;
  self->basename = (char *)malloc(n_plus_null * sizeof(char));
  if (self->basename == NULL) {
    error_python("Failed to allocate memory");
  }
  strcpy(self->basename, basename);

  return 0;
}

/**
 * @brief Read the minimal and maximal time.
 *
 * <b>returns</b> tuple containing min and max time.
 */
static PyObject *getTimeLimits(PyObject *self, PyObject *Py_UNUSED(ignored)) {
  if (!((PyObjectReader *)self)->ready) {
    error_python(
        "The logger is not ready yet."
        "Did you forget to open it with \"with\"?");
  }

  /* initialize the reader. */
  struct logger_reader *reader = &((PyObjectReader *)self)->reader;

  if (reader->verbose > 1) message("Reading time limits.");

  /* Get the time limits */
  double time_min = logger_reader_get_time_begin(reader);
  double time_max = logger_reader_get_time_end(reader);

  /* Create the output */
  PyObject *out = PyTuple_New(2);
  PyTuple_SetItem(out, 0, PyFloat_FromDouble(time_min));
  PyTuple_SetItem(out, 1, PyFloat_FromDouble(time_max));

  return (PyObject *)out;
}

#define find_field_in_module_internal(MODULE, PART)                           \
  for (int local = 0; local < MODULE##_logger_field##PART##_count; local++) { \
    const int global = MODULE##_logger_local_to_global##PART[local];          \
    const int local_shifted = local + total_number_fields;                    \
    if (field_indices[i] == global) {                                         \
      /* Check if we have the same fields for the different modules */        \
      if (current_field != NULL) {                                            \
        if (current_field->dimension !=                                       \
                python_fields[local_shifted].dimension ||                     \
            current_field->typenum != python_fields[local_shifted].typenum) { \
          error_python(                                                       \
              "The python definition of the field %s does not correspond "    \
              "between the modules.",                                         \
              MODULE##_logger_field_names##PART[local]);                      \
        }                                                                     \
      }                                                                       \
      strcpy(field_name, MODULE##_logger_field_names##PART[local]);           \
      current_field = &python_fields[local_shifted];                          \
      break;                                                                  \
    }                                                                         \
  }                                                                           \
  total_number_fields += MODULE##_logger_field##PART##_count;

/**
 * Same function as find_field_in_module_internal before but with a single
 * argument.
 */
#define find_field_in_module_single_particle_type(MODULE) \
  find_field_in_module_internal(MODULE, )

/**
 * Same function as find_field_in_module_internal but with a cleaner argument.
 */
#define find_field_in_module(MODULE, PART) \
  find_field_in_module_internal(MODULE, _##PART)

/**
 * @brief Create a list of numpy array containing the fields.
 *
 * @param output A list of array of fields.
 * @param field_indices The indices (header ordering) of the requested fields.
 * @param n_fields The number of fields requested.
 * @param n_part The number of particles of each type.
 * @param n_tot The total number of particles.
 *
 * @return The python list of numpy array.
 */
__attribute__((always_inline)) INLINE static PyObject *
logger_loader_create_output(void **output, const int *field_indices,
                            const int n_fields, const uint64_t *n_part,
                            uint64_t n_tot) {

  struct logger_python_field python_fields[100];

  /* Create the python list */
  PyObject *dict = PyDict_New();
  struct logger_python_field *current_field;

  /* Get the hydro fields */
  hydro_logger_generate_python(python_fields);
  int total_number_fields = hydro_logger_field_count;
  /* Get the gravity fields */
  gravity_logger_generate_python(python_fields + total_number_fields);
  total_number_fields += gravity_logger_field_count;
  /* Get the stars fields */
  stars_logger_generate_python(python_fields + total_number_fields);
  total_number_fields += stars_logger_field_count;
  /* Get the chemistry (part) fields */
  chemistry_logger_generate_python_part(python_fields + total_number_fields);
  total_number_fields += chemistry_logger_field_part_count;
  /* Get the chemistry (spart) fields */
  chemistry_logger_generate_python_spart(python_fields + total_number_fields);
  total_number_fields += chemistry_logger_field_spart_count;
  /* Get the star formation fields */
  star_formation_logger_generate_python(python_fields + total_number_fields);
  total_number_fields += star_formation_logger_field_count;

  /* Get all the requested fields */
  for (int i = 0; i < n_fields; i++) {
    /* Reset the variables. */
    current_field = NULL;
    total_number_fields = 0;
    char field_name[STRING_SIZE];

    /* Find the fields in the different modules. */
    find_field_in_module_single_particle_type(hydro);
    find_field_in_module_single_particle_type(gravity);
    find_field_in_module_single_particle_type(stars);
    find_field_in_module(chemistry, part);
    find_field_in_module(chemistry, spart);
    find_field_in_module_single_particle_type(star_formation);

    /* Check if we got a field */
    if (current_field == NULL) {
      error_python("Failed to find the required field");
    }
    PyObject *array = NULL;
    /* Simple case where we only have normal type */
    if (current_field->typenum != CUSTOM_NPY_TYPE) {
      if (current_field->dimension > 1) {
        npy_intp dims[2] = {n_tot, current_field->dimension};
        array = PyArray_SimpleNewFromData(2, dims, current_field->typenum,
                                          output[i]);
      } else {
        npy_intp dims = n_tot;
        array = PyArray_SimpleNewFromData(1, &dims, current_field->typenum,
                                          output[i]);
      }
    }
    /* Now the complex types */
    else {
      /* Ensure that the field is properly defined */
      if (current_field->subfields_registred != current_field->dimension) {
        error_python(
            "It seems that you forgot to register a field (found %i and "
            "expecting %i)",
            current_field->subfields_registred, current_field->dimension);
      }

      /* Creates the dtype */
      PyObject *names = PyList_New(current_field->dimension);
      PyObject *formats = PyList_New(current_field->dimension);
      PyObject *offsets = PyList_New(current_field->dimension);

      /* Gather the data into the lists */
      for (int k = 0; k < current_field->dimension; k++) {
        struct logger_python_subfield *subfield = &current_field->subfields[k];

        /* Transform the information into python */
        PyObject *name = PyUnicode_FromString(subfield->name);
        PyObject *offset = PyLong_FromSize_t(subfield->offset);
        PyObject *format = logger_python_tools_get_format_string(
            subfield->typenum, subfield->dimension);

        /* Add everything into the lists */
        PyList_SetItem(names, k, name);
        PyList_SetItem(formats, k, format);
        PyList_SetItem(offsets, k, offset);
      }

      /* Now create the dtype dictionary */
      PyObject *dtype = PyDict_New();
      PyDict_SetItemString(dtype, "names", names);
      PyDict_SetItemString(dtype, "formats", formats);
      PyDict_SetItemString(dtype, "offsets", offsets);

      /* Cleanup a bit */
      Py_DECREF(names);
      Py_DECREF(formats);
      Py_DECREF(offsets);

      /* Now get the pyarray_descr */
      PyArray_Descr *descr = NULL;
      PyArray_DescrConverter(dtype, &descr);

      /* Create the list of flags */

      /* Create the array */
      npy_intp dims = n_tot;
      array = PyArray_NewFromDescr(&PyArray_Type, descr, /* nd */ 1, &dims,
                                   /* strides */ NULL, output[i],
                                   NPY_ARRAY_CARRAY, NULL);
    }

    logger_loader_python_field_free(current_field);
    PyDict_SetItemString(dict, field_name, array);
  }

  return dict;
}

static PyObject *pyEnter(__attribute__((unused)) PyObject *self,
                         PyObject *args) {

  PyObjectReader *self_reader = (PyObjectReader *)self;
  logger_reader_init(&self_reader->reader, self_reader->basename,
                     self_reader->verbose);
  self_reader->ready = 1;

  /* Allocate the particle readers */
  PyObject *args_init = PyTuple_New(0);
  for (int i = 0; i < swift_type_count; i++) {
    PyParticleReader *part_reader = (PyParticleReader *)PyObject_CallObject(
        (PyObject *)&PyParticleReader_Type, args_init);

    /* Ensure that the object was created */
    if (part_reader == NULL) {
      return NULL;
    }
    part_reader->reader = (void *)self_reader;
    self_reader->part_reader[i] = part_reader;
    self_reader->part_reader[i]->type = i;
  }
  Py_DecRef(args_init);

  Py_INCREF(self);
  return self;
}

static PyObject *pyExit(__attribute__((unused)) PyObject *self,
                        PyObject *args) {
  PyObjectReader *self_reader = (PyObjectReader *)self;
  if (!self_reader->ready) {
    error_python("It seems that the reader was not initialized");
  }
  logger_reader_free(&self_reader->reader);
  self_reader->ready = 0;

  /* Do the same for the particle readers */
  for (int i = 0; i < swift_type_count; i++) {
    self_reader->part_reader[i]->reader = NULL;
    Py_DecRef((PyObject *)self_reader->part_reader[i]);
  }

  Py_RETURN_NONE;
}

static PyObject *pyGetListFields(__attribute__((unused)) PyObject *self,
                                 PyObject *args, PyObject *kwds) {
  PyObjectReader *self_reader = (PyObjectReader *)self;
  if (!self_reader->ready) {
    error_python(
        "The logger is not ready yet."
        "Did you forget to open it with \"with\"?");
  }

  /* input variables. */
  PyObject *types = Py_None;

  /* List of keyword arguments. */
  static char *kwlist[] = {"part_type", NULL};

  /* parse the arguments. */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &types))
    return NULL;

  /* Get the type of particles to read. */
  int read_types[swift_type_count] = {0};
  /* By default, we read everything */
  if (types == Py_None) {
    for (int i = 0; i < swift_type_count; i++) {
      read_types[i] = 1;
    }
  }
  /* Deal with the case of a single int. */
  else if (PyLong_Check(types)) {
    const size_t type = PyLong_AsSize_t(types);
    if (type >= swift_type_count) {
      error_python("Unexpected particle type %zi", type);
    }
    read_types[type] = 1;
  }
  /* Deal with the case of a list */
  else if (PyList_Check(types)) {
    const size_t size = PyList_Size(types);
    for (size_t i = 0; i < size; i++) {
      PyObject *cur = PyList_GetItem(types, i);
      const size_t type = PyLong_AsSize_t(cur);
      if (type >= swift_type_count) {
        error_python("Unexpected particle type %zi", type);
      }
      read_types[type] = 1;
    }
  }

  /* initialize the reader. */
  struct logger_reader *reader = &self_reader->reader;
  const struct header *h = &reader->log.header;

  /* Create the array to check if a field is present. */
  int *field_present = (int *)malloc(h->masks_count * sizeof(int));
  if (field_present == NULL) {
    error("Failed to allocate the memory for the fields present.");
  }

  /* Initialize the array */
  for (int i = 0; i < h->masks_count; i++) {
    field_present[i] = 1;
  }

  /* Check all the fields */
  for (int i = 0; i < swift_type_count; i++) {
    /* Skip the types that are not required */
    if (read_types[i] == 0) continue;

    /* Get the number of fields for the current particle type. */
    int number_fields = tools_get_number_fields(i);
    if (number_fields == 0) {
      continue;
    }

    /* Get the list of fields for the current particle type. */
    struct field_information *fields = (struct field_information *)malloc(
        number_fields * sizeof(struct field_information));
    if (fields == NULL) error("Failed to initialize the field information");
    tools_get_list_fields(fields, i, h);

    for (int j = 0; j < h->masks_count; j++) {
      /* Skip the fields not found in previous type. */
      if (field_present[j] == 0) continue;

      /* Check if the field is present */
      int found = 0;
      for (int k = 0; k < number_fields; k++) {
        if (strcmp(h->masks[j].name, fields[k].name) == 0) {
          found = 1;
          break;
        }
      }

      /* Set the field as not found */
      if (!found) field_present[j] = 0;
    }

    /* Free the memory */
    free(fields);
  }

  /* Count the number of fields found */
  int number_fields = 0;
  for (int i = 0; i < h->masks_count; i++) {
    number_fields += field_present[i];
  }

  /* Create the python list for the output*/
  PyObject *list = PyList_New(number_fields);
  int current = 0;
  for (int i = 0; i < h->masks_count; i++) {
    /* Keep only the field present. */
    if (field_present[i] == 0) continue;

    PyObject *name = PyUnicode_FromString(h->masks[i].name);
    PyList_SetItem(list, current, name);
    current += 1;
  }

  /* Free the memory. */
  free(field_present);

  return list;
}

/**
 * @brief shortcut to PyGetListFields for a single particle type
 */
static PyObject *pyParticleGetListFields(__attribute__((unused)) PyObject *self,
                                         PyObject *args, PyObject *kwds) {
  PyParticleReader *part_reader = (PyParticleReader *)self;
  PyObject *reader = (PyObject *)part_reader->reader;

  /* Check if we need to create the dictionary */
  const int dict_created = kwds == NULL;
  if (dict_created) {
    kwds = PyDict_New();
  }

  /* Set the particle type  */
  PyObject *part_type = PyLong_FromLong(part_reader->type);

  /* Add the list to the k-arguments */
  if (PyDict_SetItemString(kwds, "part_type", part_type) != 0) return NULL;

  /* Get the fields */
  PyObject *ret = pyGetListFields(reader, args, kwds);

  /* Cleanup */
  if (dict_created) {
    Py_DecRef(kwds);
    kwds = NULL;
  } else {
    PyDict_DelItemString(kwds, "part_type");
  }

  return ret;
}
/**
 * @brief Check the input data for #pyGetData.
 *
 * @param fields The list of fields to read.
 * @param types The type of particles to read.
 * @param part_ids The list of particle ids to read.
 * @param read_types (output) The list of particles' type to
 *   read (needs to be already allocated).
 */
void pyGetData_CheckInput(PyObject *fields, PyObject *types, PyObject *part_ids,
                          int *read_types) {

  /* Check the inputs. */
  if (!PyList_Check(fields)) {
    error_python("Expecting a list of fields");
  }

  /* Ensure that only the type or the ids are provided. */
  if (types != Py_None && part_ids != Py_None) {
    error_python(
        "Cannot request for both a type of particles and "
        "a list of ids.");
  }

  /* Ensure that a list of ids are provided. */
  if (part_ids != Py_None) {
    if (!PyList_Check(part_ids)) {
      error_python(
          "Expecting a list of numpy arrays containing the ids of the "
          "particles");
    }

    /* Check each element */
    const size_t size = PyList_Size(part_ids);
    if (size != swift_type_count) {
      error_python(
          "The list of ids should contain %i elements (one for each particle "
          "type)",
          swift_type_count);
    }
    for (size_t i = 0; i < size; i++) {
      PyArrayObject *el = (PyArrayObject *)PyList_GetItem(part_ids, i);

      /* Early abort */
      if ((PyObject *)el == Py_None) continue;

      /* Check object type */
      if (!PyArray_Check(el))
        error_python(
            "Expecting either None or a numpy array for each element in "
            "the particle ids");
      /* Check data type */
      if (PyArray_TYPE(el) != NPY_LONGLONG) {
        error_python("Expecting an array of long long integer for the ids.");
      }
    }
  }

  /* Get the type of particles to read. */
  /* By default, we read everything */
  if (types == Py_None) {
    for (int i = 0; i < swift_type_count; i++) {
      read_types[i] = 1;
    }
  }
  /* Deal with the case of a single int. */
  else if (PyLong_Check(types)) {
    const size_t type = PyLong_AsSize_t(types);
    if (type >= swift_type_count) {
      error_python("Unexpected particle type %zi", type);
    }
    read_types[type] = 1;
  }
  /* Deal with the case of a list */
  else if (PyList_Check(types)) {
    const size_t size = PyList_Size(types);
    for (size_t i = 0; i < size; i++) {
      PyObject *cur = PyList_GetItem(types, i);
      const size_t type = PyLong_AsSize_t(cur);
      if (type >= swift_type_count) {
        error_python("Unexpected particle type %zi", type);
      }
      read_types[type] = 1;
    }
  }
}

/**
 * @brief Read some fields at a given time.
 *
 * @param basename The basename of the logger files.
 * @param fields Python list containing the name of the fields (e.g.
 * Coordinates).
 * @param time The time of the fields.
 * @param verbose (Optional) The verbose level of the reader.
 *
 * @return Dictionary of numpy array containing the fields requested.
 */
static PyObject *pyGetData(__attribute__((unused)) PyObject *self,
                           PyObject *args, PyObject *kwds) {

  PyObjectReader *self_reader = (PyObjectReader *)self;
  if (!self_reader->ready) {
    error_python(
        "The logger is not ready yet."
        "Did you forget to open it with \"with\"?");
  }

  /* input variables. */
  PyObject *fields = NULL;
  double time = 0;
  PyObject *types = Py_None;
  PyObject *part_ids = Py_None;

  /* List of keyword arguments. */
  static char *kwlist[] = {"fields", "time", "filter_by_ids", "part_type",
                           NULL};

  /* parse the arguments. */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Od|OO", kwlist, &fields, &time,
                                   &part_ids, &types))
    return NULL;

  int read_types[swift_type_count] = {0};
  pyGetData_CheckInput(fields, types, part_ids, read_types);

  /* initialize the reader. */
  struct logger_reader *reader = &self_reader->reader;
  const struct header *h = &reader->log.header;

  /* Get the fields indexes from the header. */
  const int n_fields = PyList_Size(fields);
  int *field_indices = (int *)malloc(n_fields * sizeof(int));
  for (int i = 0; i < n_fields; i++) {
    field_indices[i] = -1;

    /* Get the an item in the list. */
    PyObject *field = PyList_GetItem(fields, i);
    if (!PyUnicode_Check(field)) {
      error_python("Expecting list of string for the fields");
    }

    /* Convert into C string. */
    Py_ssize_t size = 0;
    const char *field_name = PyUnicode_AsUTF8AndSize(field, &size);

    /* Find the equivalent field inside the header. */
    for (int j = 0; j < h->masks_count; j++) {
      if (strcmp(field_name, h->masks[j].name) == 0) {
        field_indices[i] = j;
        break;
      }
    }

    /* Check if we found a field. */
    if (field_indices[i] == -1) {
      error_python("Failed to find the field %s", field_name);
    }
  }

  /* Set the time. */
  logger_reader_set_time(reader, time);

  /* Get the number of particles. */
  uint64_t n_part[swift_type_count];
  if (part_ids == Py_None) {
    /* If no ids provided, read everything from the index files */
    logger_reader_get_number_particles(reader, n_part, read_types);
  }
  /* If the ids are provided, use them to get the number of particles */
  else {
    for (int i = 0; i < swift_type_count; i++) {
      PyObject *el = PyList_GetItem(part_ids, i);
      n_part[i] = el == Py_None ? 0 : PyArray_Size(el);
    }
  }

  /* Count the total number of particles. */
  uint64_t n_tot = 0;
  for (int i = 0; i < swift_type_count; i++) {
    n_tot += n_part[i];
  }

  /* Allocate the output memory. */
  void **output = malloc(n_fields * sizeof(void *));
  for (int i = 0; i < n_fields; i++) {
    output[i] = malloc(n_tot * h->masks[field_indices[i]].size);
  }

  /* Read the particles. */
  if (part_ids == Py_None) {
    logger_reader_read_all_particles(reader, time, logger_reader_lin,
                                     field_indices, n_fields, output, n_part);
  } else {
    /* Get the arrays in a C compatible way */
    PyObject *list[swift_type_count] = {NULL};
    long long *c_part_ids[swift_type_count];
    for (int i = 0; i < swift_type_count; i++) {
      PyObject *el = PyList_GetItem(part_ids, i);

      /* Check if a list is provided */
      if (el == Py_None) {
        c_part_ids[i] = NULL;
        continue;
      }

      /* Copy the ids in order to modify them */
      list[i] = PyArray_NewCopy((PyArrayObject *)el, NPY_CORDER);
      c_part_ids[i] = (long long *)PyArray_DATA((PyArrayObject *)list[i]);
    }

    /* Read the particles. */
    logger_reader_read_particles_from_ids(reader, time, logger_reader_lin,
                                          field_indices, n_fields, output,
                                          n_part, c_part_ids);

    /* Free the memory and recompute n_tot */
    n_tot = 0;
    for (int i = 0; i < swift_type_count; i++) {
      if (list[i] != NULL) Py_DECREF(list[i]);

      n_tot += n_part[i];
    }
  }

  /* Create the python output. */
  PyObject *array = logger_loader_create_output(output, field_indices, n_fields,
                                                n_part, n_tot);

  /* Free the reader. */
  free(field_indices);

  return array;
}

/**
 * @brief Shortcut to #pyGetData for a single particle type.
 */
static PyObject *pyParticleGetData(__attribute__((unused)) PyObject *self,
                                   PyObject *args, PyObject *kwds) {

  PyParticleReader *part_reader = (PyParticleReader *)self;
  PyObject *reader = (PyObject *)part_reader->reader;

  /* input variables. */
  PyObject *fields = NULL;
  double time = 0;
  PyObject *part_ids = Py_None;

  /* List of keyword arguments. */
  static char *kwlist[] = {"fields", "time", "filter_by_ids", NULL};

  /* parse the arguments. */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Od|O", kwlist, &fields, &time,
                                   &part_ids))
    return NULL;

  /* Check if we need to create the dictionary */
  const int dict_created = kwds == NULL;
  if (dict_created) {
    kwds = PyDict_New();
  }

  /* Adapt the inputs in case we are filtering */
  PyObject *part_ids_list = part_ids;
  if (part_ids != Py_None) {
    if (!PyArray_Check(part_ids)) {
      error("You need to provide a numpy array of ids");
    }

    part_ids_list = PyList_New(swift_type_count);
    for (int i = 0; i < swift_type_count; i++) {
      PyList_SetItem(part_ids_list, i, Py_None);
    }
    PyList_SetItem(part_ids_list, part_reader->type, part_ids);

    if (PyDict_SetItemString(kwds, "filter_by_ids", part_ids_list) != 0)
      return NULL;
  }
  /* Now the case without filtering */
  else {
    /* Set the particle type  */
    PyObject *part_type = PyLong_FromLong(part_reader->type);

    /* Add the list to the k-arguments */
    if (PyDict_SetItemString(kwds, "part_type", part_type) != 0) return NULL;
  }

  /* Get the fields */
  PyObject *ret = pyGetData(reader, args, kwds);

  /* Cleanup */
  if (dict_created) {
    Py_DecRef(kwds);
    kwds = NULL;
  } else {
    if (part_ids == Py_None) PyDict_DelItemString(kwds, "part_type");
  }

  return ret;
}

/* definition of the method table. */
static PyMethodDef libloggerMethods[] = {
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/* Definition of the Reader methods */
static PyMethodDef libloggerReaderMethods[] = {
    {"get_time_limits", getTimeLimits, METH_NOARGS,
     "Read the time limits of the simulation.\n\n"
     "Returns\n"
     "-------\n\n"
     "times: tuple\n"
     "  time min, time max\n"},
    {"get_data", (PyCFunction)pyGetData, METH_VARARGS | METH_KEYWORDS,
     "Read some fields from the logfile at a given time.\n\n"
     "Parameters\n"
     "----------\n\n"
     "fields: list\n"
     "  The list of fields (e.g. 'Coordinates', 'Entropies', ...)\n\n"
     "time: float\n"
     "  The time at which the fields must be read.\n\n"
     "Returns\n"
     "-------\n\n"
     "list_of_fields: list\n"
     "  Each element is a numpy array containing the corresponding field.\n"},
    {"get_list_fields", (PyCFunction)pyGetListFields,
     METH_VARARGS | METH_KEYWORDS,
     "Read the list of available fields in the logfile.\n\n"
     "Parameters\n"
     "----------\n\n"
     "type: int, list\n"
     "  The particle type for the list of fields\n\n"
     "Returns\n"
     "-------\n"
     "fields: tuple\n"
     "  The list of fields present in the logfile.\n"},
    {"__enter__", pyEnter, METH_VARARGS, ""},
    {"__exit__", pyExit, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef libloggermodule = {
    PyModuleDef_HEAD_INIT,
    "liblogger",
    "Module reading a SWIFTsim logger snapshot",
    -1,
    libloggerMethods,
    NULL, /* m_slots */
    NULL, /* m_traverse */
    NULL, /* m_clear */
    NULL  /* m_free */
};

static PyMemberDef PyObjectReader_members[] = {
    {"gas", T_OBJECT_EX, offsetof(PyObjectReader, part_reader), READONLY},
    {"dark_matter", T_OBJECT_EX,
     offsetof(PyObjectReader, part_reader) + sizeof(PyParticleReader *),
     READONLY},
    {"dark_matter_background", T_OBJECT_EX,
     offsetof(PyObjectReader, part_reader) + 2 * sizeof(PyParticleReader *),
     READONLY},
    {"sinks", T_OBJECT_EX,
     offsetof(PyObjectReader, part_reader) + 3 * sizeof(PyParticleReader *),
     READONLY},
    {"stars", T_OBJECT_EX,
     offsetof(PyObjectReader, part_reader) + 4 * sizeof(PyParticleReader *),
     READONLY},
    {"black_holes", T_OBJECT_EX,
     offsetof(PyObjectReader, part_reader) + 5 * sizeof(PyParticleReader *),
     READONLY},
    {NULL}};

static PyTypeObject PyObjectReader_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "liblogger.Reader",
    .tp_basicsize = sizeof(PyObjectReader),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc =
        "Deals with a logger file and provides the data through"
        " its different methods. The reader is really open with the method "
        "__enter__"
        "(called by `with ...``). Outside the area delimited by with, "
        "the logger is not accessible.\n\n"
        "Parameters\n"
        "----------\n\n"
        "basename: str\n"
        "    Basename of the logfile.\n"
        "verbose: int (optional)\n"
        "    Verbose level\n\n"
        "Methods\n"
        "-------\n\n"
        "get_time_limits\n"
        "    Provides the initial and final time of the simulation\n"
        "get_particle_data\n"
        "    Reads some particle's data from the logfile\n\n"
        "Examples\n"
        "--------\n\n"
        ">>> import liblogger as logger\n"
        ">>> with logger.Reader(\"index_0000\") as reader:\n"
        ">>>    t0, t1 = reader.get_time_limits()\n"
        ">>>    fields = reader.get_list_fields()\n"
        ">>>    if \"Coordinates\" not in fields:\n"
        ">>>        raise Exception(\"Field Coordinates not present in the "
        "logfile.\")\n"
        ">>>    pos, ent = reader.get_particle_data([\"Coordinates\", "
        "\"Entropies\"]"
        ", 0.5 * (t0 + t1))\n",
    .tp_init = (initproc)Reader_init,
    .tp_dealloc = (destructor)Reader_dealloc,
    .tp_new = Reader_new,
    .tp_methods = libloggerReaderMethods,
    .tp_members = PyObjectReader_members,
};

/* Definition of the Reader methods */
static PyMethodDef libloggerParticleReaderMethods[] = {
    {"get_data", (PyCFunction)pyParticleGetData, METH_VARARGS | METH_KEYWORDS,
     "Read some fields from the logfile at a given time.\n\n"
     "Parameters\n"
     "----------\n\n"
     "fields: list\n"
     "  The list of fields (e.g. 'Coordinates', 'Entropies', ...)\n\n"
     "time: float\n"
     "  The time at which the fields must be read.\n\n"
     "Returns\n"
     "-------\n\n"
     "list_of_fields: list\n"
     "  Each element is a numpy array containing the corresponding field.\n"},
    {"get_list_fields", (PyCFunction)pyParticleGetListFields,
     METH_VARARGS | METH_KEYWORDS,
     "Read the list of available fields in the logfile.\n\n"
     "Parameters\n"
     "----------\n\n"
     "type: int, list\n"
     "  The particle type for the list of fields\n\n"
     "Returns\n"
     "-------\n"
     "fields: tuple\n"
     "  The list of fields present in the logfile.\n"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static PyTypeObject PyParticleReader_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "liblogger.ParticleReader",
    .tp_basicsize = sizeof(PyParticleReader),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc =
        "Implements liblogger.Reader for a single particle type.\n"
        "See the documentation of liblogger.Reader for more information.\n",
    .tp_init = NULL,
    .tp_dealloc = (destructor)particle_reader_dealloc,
    .tp_new = particle_reader_new,
    .tp_methods = libloggerParticleReaderMethods,
};

PyMODINIT_FUNC PyInit_liblogger(void) {

  /* Create the module. */
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  /* Deal with SWIFT clock */
  clocks_set_cpufreq(0);

  /* Prepare the classes */
  /* Object reader */
  if (PyType_Ready(&PyObjectReader_Type) < 0) {
    return NULL;
  }
  Py_INCREF(&PyObjectReader_Type);

  /* Particle reader */
  if (PyType_Ready(&PyParticleReader_Type) < 0) {
    return NULL;
  }
  Py_INCREF(&PyParticleReader_Type);

  /* Add the reader object */
  PyModule_AddObject(m, "Reader", (PyObject *)&PyObjectReader_Type);
  PyModule_AddIntConstant(m, "gas", swift_type_gas);
  PyModule_AddIntConstant(m, "dark_matter", swift_type_dark_matter);
  PyModule_AddIntConstant(m, "dark_matter_background",
                          swift_type_dark_matter_background);
  PyModule_AddIntConstant(m, "sinks", swift_type_sink);
  PyModule_AddIntConstant(m, "stars", swift_type_stars);
  PyModule_AddIntConstant(m, "black_holes", swift_type_black_hole);

  /* Import numpy. */
  import_array();

  return m;
}

#endif  // HAVE_PYTHON
