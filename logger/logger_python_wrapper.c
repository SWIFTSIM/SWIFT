#include "logger_header.h"
#include "logger_io.h"
#include "logger_particle.h"
#include "logger_time.h"
#include "logger_reader.h"

#include <Python.h>
#include <errno.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief load data from the offset without any interpolation
 *
 * @param offset PyArrayObject list of offset for each particle
 * @param filename string filename of the dump file
 * @return dictionnary containing the data read
 */
static PyObject *loadFromIndex(__attribute__((unused)) PyObject *self,
                               PyObject *args) {

  struct header h;

  /* input */
  PyArrayObject *offset = NULL;
  char *filename = NULL;

  /* output */
  PyArrayObject *pos = NULL;
  PyArrayObject *vel = NULL;
  PyArrayObject *acc = NULL;
  PyArrayObject *entropy = NULL;
  PyArrayObject *h_sph = NULL;
  PyArrayObject *rho = NULL;
  PyArrayObject *mass = NULL;
  PyArrayObject *id = NULL;

  size_t time_offset;

  /* parse arguments */

  if (!PyArg_ParseTuple(args, "OsL", &offset, &filename, &time_offset))
    return NULL;

  if (!PyArray_Check(offset)) {
    error_no_return(ENOTSUP, "Offset is not a numpy array");
    return NULL;
  }
  if (PyArray_NDIM(offset) != 1) {
    error_no_return(ENOTSUP, "Offset is not a 1 dimensional array");
    return NULL;
  }
  if (PyArray_TYPE(offset) != NPY_UINT64) {
    error_no_return(ENOTSUP, "Offset does not contain unsigned int");
    return NULL;
  }

  /* open file */
  int fd;
  void *map;
  if (io_open_file(filename, &fd, &map) != 0) return NULL;

  /* read header */
  if (header_read(&h, map) != 0) return NULL;

  /* reverse offset if needed */
  if (!h.forward_offset) {
    if (io_close_file(&fd, &map) != 0) return NULL;

    if (reverse_offset(filename) != 0) return NULL;

    if (io_open_file(filename, &fd, &map) != 0) return NULL;

    /* Reset header */
    header_free(&h);

    if (header_read(&h, map) != 0) return NULL;
  }

  /* read timestamps */
  struct time_array times;

  if (time_array_init(&times, &h, map, fd) != 0) return NULL;

  time_array_print(&times);
  /* get required time */
  double time = time_array_get_time(&times, time_offset);

  /* init array */
  npy_intp dim[2];
  dim[0] = PyArray_DIMS(offset)[0];
  dim[1] = DIM;

  /* init output */
  if (header_is_present(&h, "position")) {
    pos = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }

  if (header_is_present(&h, "velocity")) {
    vel = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_FLOAT);
  }

  if (header_is_present(&h, "acceleration")) {
    acc = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_FLOAT);
  }

  if (header_is_present(&h, "entropy")) {
    entropy =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_is_present(&h, "cutoff radius")) {
    h_sph =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_is_present(&h, "density")) {
    rho =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
  }

  if (header_is_present(&h, "consts")) {
    mass =
        (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_FLOAT);
    id = (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(offset), NPY_ULONG);
  }

  int error_code = 0;

  /* check data type in particles */
  if (!particle_check_data_type(&h)) {
    error_no_return(ENOTSUP, "Particle data type are not compatible");
    return NULL;
  }

  /* loop over all particles */
  for (npy_intp i = 0; i < PyArray_DIMS(offset)[0]; i++) {
    struct particle part;

    size_t *offset_particle = (size_t *)PyArray_GETPTR1(offset, i);

    error_code = particle_read(&part, &h, map, offset_particle, time,
                               reader_lin, &times);
    if (error_code != 0) return NULL;

    double *dtmp;
    float *ftmp;
    size_t *stmp;

    /* copy data */
    for (size_t k = 0; k < DIM; k++) {
      if (pos) {
        dtmp = PyArray_GETPTR2(pos, i, k);
        *dtmp = part.pos[k];
      }

      if (vel) {
        ftmp = PyArray_GETPTR2(vel, i, k);
        *ftmp = part.vel[k];
      }

      if (acc) {
        ftmp = PyArray_GETPTR2(acc, i, k);
        *ftmp = part.acc[k];
      }
    }

    if (entropy) {
      ftmp = PyArray_GETPTR1(entropy, i);
      *ftmp = part.entropy;
    }

    if (rho) {
      ftmp = PyArray_GETPTR1(rho, i);
      *ftmp = part.density;
    }

    if (h_sph) {
      ftmp = PyArray_GETPTR1(h_sph, i);
      *ftmp = part.h;
    }

    if (mass) {
      ftmp = PyArray_GETPTR1(mass, i);
      *ftmp = part.mass;
    }

    if (id) {
      stmp = PyArray_GETPTR1(id, i);
      *stmp = part.id;
    }
  }

  header_free(&h);

  /* construct return */
  PyObject *dict = PyDict_New();
  PyObject *key = PyUnicode_FromString("position");
  PyDict_SetItem(dict, key, PyArray_Return(pos));

  if (vel) {
    key = PyUnicode_FromString("velocity");
    PyDict_SetItem(dict, key, PyArray_Return(vel));
  }

  if (acc) {
    key = PyUnicode_FromString("acceleration");
    PyDict_SetItem(dict, key, PyArray_Return(acc));
  }

  if (entropy) {
    key = PyUnicode_FromString("entropy");
    PyDict_SetItem(dict, key, PyArray_Return(entropy));
  }

  if (rho) {
    key = PyUnicode_FromString("rho");
    PyDict_SetItem(dict, key, PyArray_Return(rho));
  }

  if (h_sph) {
    key = PyUnicode_FromString("h_sph");
    PyDict_SetItem(dict, key, PyArray_Return(h_sph));
  }

  if (mass) {
    key = PyUnicode_FromString("mass");
    PyDict_SetItem(dict, key, PyArray_Return(mass));
  }

  if (id) {
    key = PyUnicode_FromString("id");
    PyDict_SetItem(dict, key, PyArray_Return(id));
  }

  io_close_file(&fd, &map);

  return dict;
}

/**
 * @brief Reverse offset in dump file
 *
 * @param filename string filename of the dump file
 */
static PyObject *pyReverseOffset(__attribute__((unused)) PyObject *self,
                                 PyObject *args) {
  /* input */
  char *filename = NULL;

  /* parse arguments */

  if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

  if (reverse_offset(filename) != 0) return NULL;

  return Py_BuildValue("");
}

/* definition of the method table */

static PyMethodDef libloggerMethods[] = {
    {"loadFromIndex", loadFromIndex, METH_VARARGS,
     "Load snapshot directly from the offset in an index file."},
    {"reverseOffset", pyReverseOffset, METH_VARARGS,
     "Reverse the offset (from pointing backward to forward)."},

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

PyMODINIT_FUNC PyInit_liblogger(void) {
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  import_array();

  return m;
}
