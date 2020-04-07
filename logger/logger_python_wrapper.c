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
#include "logger_reader.h"
#include "logger_time.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <errno.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  PyObject_HEAD struct logger_particle part;
} PyLoggerParticle;

static PyTypeObject PyLoggerParticle_Type;
const char *particle_name = "Particle";

PyArray_Descr *logger_particle_descr;

/**
 * @brief load data from the index files.
 *
 * <b>basename</b> Base name of the logger files.
 *
 * <b>time</b> The time requested.
 *
 * <b>verbose</b> Verbose level.
 *
 * <b>returns</b> dictionnary containing the data read.
 */
static PyObject *loadSnapshotAtTime(__attribute__((unused)) PyObject *self,
                                    PyObject *args) {

  /* declare variables. */
  char *basename = NULL;

  double time = 0;
  int verbose = 1;

  /* parse arguments. */
  if (!PyArg_ParseTuple(args, "sd|i", &basename, &time, &verbose)) return NULL;

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, basename, verbose);

  if (verbose > 1) message("Reading particles.");

  /* Number of particles in the index files */
  npy_intp n_tot = 0;

  /* Set the reading time */
  logger_reader_set_time(&reader, time);

  /* Get the number of particles */
  int n_type = 0;
  const uint64_t *n_parts =
      logger_reader_get_number_particles(&reader, &n_type);
  for (int i = 0; i < n_type; i++) {
    n_tot += n_parts[i];
  }

  if (verbose > 0) {
    message("Found %lu particles", n_tot);
  }

  /* Allocate the output memory */
  PyArrayObject *out = (PyArrayObject *)PyArray_SimpleNewFromDescr(
      1, &n_tot, logger_particle_descr);

  void *data = PyArray_DATA(out);
  /* Allows to use threads */
  Py_BEGIN_ALLOW_THREADS;

  /* Read the particle. */
  logger_reader_read_all_particles(&reader, time, logger_reader_const, data,
                                   n_tot);

  /* No need of threads anymore */
  Py_END_ALLOW_THREADS;

  /* Free the memory. */
  logger_reader_free(&reader);

  return (PyObject *)out;
}

/**
 * @brief Read the minimal and maximal time.
 *
 * <b>basename</b> Base name of the logger files.
 *
 * <b>verbose</b> Verbose level.
 *
 * <b>returns</b> tuple containing min and max time.
 */
static PyObject *getTimeLimits(__attribute__((unused)) PyObject *self,
                               PyObject *args) {

  /* declare variables. */
  char *basename = NULL;

  int verbose = 1;

  /* parse arguments. */
  if (!PyArg_ParseTuple(args, "s|i", &basename, &verbose)) return NULL;

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, basename, verbose);

  if (verbose > 1) message("Reading time limits.");

  /* Get the time limits */
  double time_min = logger_reader_get_time_begin(&reader);
  double time_max = logger_reader_get_time_end(&reader);

  /* Free the memory. */
  logger_reader_free(&reader);

  /* Create the output */
  PyObject *out = PyTuple_New(2);
  PyTuple_SetItem(out, 0, PyFloat_FromDouble(time_min));
  PyTuple_SetItem(out, 1, PyFloat_FromDouble(time_max));

  return (PyObject *)out;
}

/**
 * @brief Reverse offset in log file
 *
 * <b>filename</b> string filename of the log file
 * <b>verbose</b> Verbose level
 */
static PyObject *pyReverseOffset(__attribute__((unused)) PyObject *self,
                                 PyObject *args) {
  /* input variables. */
  char *filename = NULL;

  int verbose = 0;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "s|i", &filename, &verbose)) return NULL;

  /* initialize the reader which reverse the offset if necessary. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);

  /* Free the reader. */
  logger_reader_free(&reader);

  return Py_BuildValue("");
}

/**
 * @brief Move forward in time an array of particles.
 *
 * <b>filename</b> string filename of the log file.
 * <b>parts</b> Numpy array containing the particles to evolve.
 * <b>time</b> Time requested for the particles.
 * <b>verbose</b> Verbose level
 *
 * <b>returns</b> The evolved array of particles.
 */
static PyObject *pyMoveForwardInTime(__attribute__((unused)) PyObject *self,
                                     PyObject *args) {
  /* input variables. */
  char *filename = NULL;

  int verbose = 0;
  PyArrayObject *parts = NULL;
  double time = 0;
  int new_array = 1;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "sOd|ii", &filename, &parts, &time, &verbose,
                        &new_array))
    return NULL;

  /* Check parts */
  if (!PyArray_Check(parts)) {
    error("Expecting a numpy array of particles.");
  }

  if (PyArray_NDIM(parts) != 1) {
    error("Expecting a 1D array of particles.");
  }

  if (PyArray_TYPE(parts) != logger_particle_descr->type_num) {
    error("Expecting an array of particles.");
  }

  /* Create the interpolated array. */
  PyArrayObject *interp = NULL;
  if (new_array) {
    interp =
        (PyArrayObject *)PyArray_NewLikeArray(parts, NPY_ANYORDER, NULL, 0);

    /* Check if the allocation was fine */
    if (interp == NULL) {
      return NULL;
    }

    /* Reference stolen in PyArray_NewLikeArray => incref */
    Py_INCREF(PyArray_DESCR(parts));
  } else {
    interp = parts;
    // We return it, therefore one more reference exists.
    Py_INCREF(interp);
  }

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);

  /* Get the offset of the requested time. */
  size_t offset = logger_reader_get_next_offset_from_time(&reader, time);

  /* Loop over all the particles */
  size_t N = PyArray_DIM(parts, 0);
  for (size_t i = 0; i < N; i++) {

    /* Obtain the required particle records. */
    struct logger_particle *p = PyArray_GETPTR1(parts, i);

    /* Check that we are really going forward in time. */
    if (time < p->time) {
      error("Requesting to go backward in time (%g < %g)", time, p->time);
    }
    struct logger_particle new;
    logger_reader_get_next_particle(&reader, p, &new, offset);

    /* Interpolate the particle. */
    struct logger_particle *p_ret = PyArray_GETPTR1(interp, i);
    *p_ret = *p;

    logger_particle_interpolate(p_ret, &new, time);
  }

  /* Free the reader. */
  logger_reader_free(&reader);

  return PyArray_Return(interp);
}

/* definition of the method table. */

static PyMethodDef libloggerMethods[] = {
    {"loadSnapshotAtTime", loadSnapshotAtTime, METH_VARARGS,
     "Load a snapshot directly from the logger using the index files.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "time: double\n"
     "  The (double) time of the snapshot.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "Returns\n"
     "-------\n\n"
     "snapshot: dict\n"
     "  The full output generated for the whole file.\n"},
    {"reverseOffset", pyReverseOffset, METH_VARARGS,
     "Reverse the offset (from pointing backward to forward).\n\n"
     "Parameters\n"
     "----------\n\n"
     "filename: str\n"
     "  The filename of the log file.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n"},
    {"getTimeLimits", getTimeLimits, METH_VARARGS,
     "Read the time limits of the simulation.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "Returns\n"
     "-------\n\n"
     "times: tuple\n"
     "  time min, time max\n"},
    {"moveForwardInTime", pyMoveForwardInTime, METH_VARARGS,
     "Move the particles forward in time.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "parts: np.array\n"
     "  The array of particles.\n\n"
     "time: double\n"
     "  The requested time for the particles.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "new_array: bool, optional\n"
     "  Does the function return a new array (or use the provided one)?\n\n"
     "Returns\n"
     "-------\n\n"
     "parts: np.array\n"
     "  The particles at the requested time.\n"},

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

#define CREATE_FIELD(fields, name, field_name, type)                      \
  ({                                                                      \
    PyObject *tuple = PyTuple_New(2);                                     \
    PyTuple_SetItem(tuple, 0, (PyObject *)PyArray_DescrFromType(type));   \
    PyTuple_SetItem(                                                      \
        tuple, 1,                                                         \
        PyLong_FromSize_t(offsetof(struct logger_particle, field_name))); \
    PyDict_SetItem(fields, PyUnicode_FromString(name), tuple);            \
  })

#define CREATE_FIELD_3D(fields, name, field_name, type)                   \
  ({                                                                      \
    /* Create the 3D descriptor */                                        \
    PyArray_Descr *vec = PyArray_DescrNewFromType(type);                  \
    vec->subarray = malloc(sizeof(PyArray_ArrayDescr));                   \
    vec->subarray->base = PyArray_DescrFromType(type);                    \
    vec->subarray->shape = PyTuple_New(1);                                \
    PyTuple_SetItem(vec->subarray->shape, 0, PyLong_FromSize_t(3));       \
                                                                          \
    /* Create the field */                                                \
    PyObject *tuple = PyTuple_New(2);                                     \
    PyTuple_SetItem(tuple, 0, (PyObject *)vec);                           \
    PyTuple_SetItem(                                                      \
        tuple, 1,                                                         \
        PyLong_FromSize_t(offsetof(struct logger_particle, field_name))); \
    PyDict_SetItem(fields, PyUnicode_FromString(name), tuple);            \
  })

void pylogger_particle_define_typeobject(void) {

  PyLoggerParticle_Type.tp_name = particle_name;
  PyLoggerParticle_Type.tp_print = NULL;
  PyType_Ready(&PyLoggerParticle_Type);
}

void pylogger_particle_define_descr(void) {
  /* Generate list of field names */
  PyObject *names = PyTuple_New(9);
  PyTuple_SetItem(names, 0, PyUnicode_FromString("positions"));
  PyTuple_SetItem(names, 1, PyUnicode_FromString("velocities"));
  PyTuple_SetItem(names, 2, PyUnicode_FromString("accelerations"));
  PyTuple_SetItem(names, 3, PyUnicode_FromString("entropies"));
  PyTuple_SetItem(names, 4, PyUnicode_FromString("smoothing_lengths"));
  PyTuple_SetItem(names, 5, PyUnicode_FromString("densities"));
  PyTuple_SetItem(names, 6, PyUnicode_FromString("masses"));
  PyTuple_SetItem(names, 7, PyUnicode_FromString("ids"));
  PyTuple_SetItem(names, 8, PyUnicode_FromString("times"));

  /* Generate list of fields */
  PyObject *fields = PyDict_New();
  CREATE_FIELD_3D(fields, "positions", pos, NPY_DOUBLE);
  CREATE_FIELD_3D(fields, "velocities", vel, NPY_FLOAT32);
  CREATE_FIELD_3D(fields, "accelerations", acc, NPY_FLOAT32);
  CREATE_FIELD(fields, "entropies", entropy, NPY_FLOAT32);
  CREATE_FIELD(fields, "smoothing_lengths", h, NPY_FLOAT32);
  CREATE_FIELD(fields, "densities", density, NPY_FLOAT32);
  CREATE_FIELD(fields, "masses", mass, NPY_FLOAT32);
  CREATE_FIELD(fields, "ids", id, NPY_LONGLONG);
  CREATE_FIELD(fields, "times", time, NPY_DOUBLE);

  /* Generate descriptor */
  logger_particle_descr = PyObject_New(PyArray_Descr, &PyArrayDescr_Type);
  logger_particle_descr->typeobj = &PyLoggerParticle_Type;
  // V if for an arbitrary kind of array
  logger_particle_descr->kind = 'V';
  // Not well documented (seems any value is fine)
  logger_particle_descr->type = 'v';
  // Native byte ordering
  logger_particle_descr->byteorder = '=';
  // Flags
  logger_particle_descr->flags = NPY_USE_GETITEM | NPY_USE_SETITEM;
  // id of the data type (assigned automatically)
  logger_particle_descr->type_num = 0;
  // Size of an element (using more size than required in order to log
  // everything)
  logger_particle_descr->elsize = sizeof(struct logger_particle);
  // alignment (doc magic)
  logger_particle_descr->alignment = offsetof(
      struct {
        char c;
        struct logger_particle v;
      },
      v);
  // no subarray
  logger_particle_descr->subarray = NULL;
  // functions
  logger_particle_descr->f = NULL;
  // Meta data
  logger_particle_descr->metadata = NULL;
  logger_particle_descr->c_metadata = NULL;
  logger_particle_descr->names = names;
  logger_particle_descr->fields = fields;
}

PyMODINIT_FUNC PyInit_liblogger(void) {
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  /* Deal with SWIFT clock */
  clocks_set_cpufreq(0);

  import_array();
  /* Define the type object */
  pylogger_particle_define_typeobject();

  /* Define the descr of the logger_particle */
  pylogger_particle_define_descr();

  return m;
}
