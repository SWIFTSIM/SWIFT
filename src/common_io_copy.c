/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* This object's header. */
#include "common_io.h"

/* Local includes. */
#include "engine.h"
#include "io_properties.h"
#include "threadpool.h"

/**
 * @brief Mapper function to copy #part or #gpart fields into a buffer.
 */
void io_copy_mapper(void* restrict temp, int N, void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;

  /* How far are we with this chunk? */
  char* restrict temp_c = (char*)temp;
  const ptrdiff_t delta = (temp_c - props.start_temp_c) / copySize;

  for (int k = 0; k < N; k++) {
    memcpy(&temp_c[k * copySize], props.field + (delta + k) * props.partSize,
           copySize);
  }
}

/**
 * @brief Mapper function to copy #part into a buffer of floats using a
 * conversion function.
 */
void io_convert_part_f_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct part* restrict parts = props.parts;
  const struct xpart* restrict xparts = props.xparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  float* restrict temp_f = (float*)temp;
  const ptrdiff_t delta = (temp_f - props.start_temp_f) / dim;

  for (int i = 0; i < N; i++)
    props.convert_part_f(e, parts + delta + i, xparts + delta + i,
                         &temp_f[i * dim]);
}

/**
 * @brief Mapper function to copy #part into a buffer of ints using a
 * conversion function.
 */
void io_convert_part_i_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct part* restrict parts = props.parts;
  const struct xpart* restrict xparts = props.xparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  int* restrict temp_i = (int*)temp;
  const ptrdiff_t delta = (temp_i - props.start_temp_i) / dim;

  for (int i = 0; i < N; i++)
    props.convert_part_i(e, parts + delta + i, xparts + delta + i,
                         &temp_i[i * dim]);
}

/**
 * @brief Mapper function to copy #part into a buffer of doubles using a
 * conversion function.
 */
void io_convert_part_d_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct part* restrict parts = props.parts;
  const struct xpart* restrict xparts = props.xparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  double* restrict temp_d = (double*)temp;
  const ptrdiff_t delta = (temp_d - props.start_temp_d) / dim;

  for (int i = 0; i < N; i++)
    props.convert_part_d(e, parts + delta + i, xparts + delta + i,
                         &temp_d[i * dim]);
}

/**
 * @brief Mapper function to copy #part into a buffer of doubles using a
 * conversion function.
 */
void io_convert_part_l_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct part* restrict parts = props.parts;
  const struct xpart* restrict xparts = props.xparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  long long* restrict temp_l = (long long*)temp;
  const ptrdiff_t delta = (temp_l - props.start_temp_l) / dim;

  for (int i = 0; i < N; i++)
    props.convert_part_l(e, parts + delta + i, xparts + delta + i,
                         &temp_l[i * dim]);
}

/**
 * @brief Mapper function to copy #gpart into a buffer of floats using a
 * conversion function.
 */
void io_convert_gpart_f_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct gpart* restrict gparts = props.gparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  float* restrict temp_f = (float*)temp;
  const ptrdiff_t delta = (temp_f - props.start_temp_f) / dim;

  for (int i = 0; i < N; i++)
    props.convert_gpart_f(e, gparts + delta + i, &temp_f[i * dim]);
}

/**
 * @brief Mapper function to copy #gpart into a buffer of ints using a
 * conversion function.
 */
void io_convert_gpart_i_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct gpart* restrict gparts = props.gparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  int* restrict temp_i = (int*)temp;
  const ptrdiff_t delta = (temp_i - props.start_temp_i) / dim;

  for (int i = 0; i < N; i++)
    props.convert_gpart_i(e, gparts + delta + i, &temp_i[i * dim]);
}

/**
 * @brief Mapper function to copy #gpart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_gpart_d_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct gpart* restrict gparts = props.gparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  double* restrict temp_d = (double*)temp;
  const ptrdiff_t delta = (temp_d - props.start_temp_d) / dim;

  for (int i = 0; i < N; i++)
    props.convert_gpart_d(e, gparts + delta + i, &temp_d[i * dim]);
}

/**
 * @brief Mapper function to copy #gpart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_gpart_l_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct gpart* restrict gparts = props.gparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  long long* restrict temp_l = (long long*)temp;
  const ptrdiff_t delta = (temp_l - props.start_temp_l) / dim;

  for (int i = 0; i < N; i++)
    props.convert_gpart_l(e, gparts + delta + i, &temp_l[i * dim]);
}

/**
 * @brief Mapper function to copy #spart into a buffer of floats using a
 * conversion function.
 */
void io_convert_spart_f_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct spart* restrict sparts = props.sparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  float* restrict temp_f = (float*)temp;
  const ptrdiff_t delta = (temp_f - props.start_temp_f) / dim;

  for (int i = 0; i < N; i++)
    props.convert_spart_f(e, sparts + delta + i, &temp_f[i * dim]);
}

/**
 * @brief Mapper function to copy #spart into a buffer of ints using a
 * conversion function.
 */
void io_convert_spart_i_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct spart* restrict sparts = props.sparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  int* restrict temp_i = (int*)temp;
  const ptrdiff_t delta = (temp_i - props.start_temp_i) / dim;

  for (int i = 0; i < N; i++)
    props.convert_spart_i(e, sparts + delta + i, &temp_i[i * dim]);
}

/**
 * @brief Mapper function to copy #spart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_spart_d_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct spart* restrict sparts = props.sparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  double* restrict temp_d = (double*)temp;
  const ptrdiff_t delta = (temp_d - props.start_temp_d) / dim;

  for (int i = 0; i < N; i++)
    props.convert_spart_d(e, sparts + delta + i, &temp_d[i * dim]);
}

/**
 * @brief Mapper function to copy #spart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_spart_l_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct spart* restrict sparts = props.sparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  long long* restrict temp_l = (long long*)temp;
  const ptrdiff_t delta = (temp_l - props.start_temp_l) / dim;

  for (int i = 0; i < N; i++)
    props.convert_spart_l(e, sparts + delta + i, &temp_l[i * dim]);
}

/**
 * @brief Mapper function to copy #bpart into a buffer of floats using a
 * conversion function.
 */
void io_convert_bpart_f_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct bpart* restrict bparts = props.bparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  float* restrict temp_f = (float*)temp;
  const ptrdiff_t delta = (temp_f - props.start_temp_f) / dim;

  for (int i = 0; i < N; i++)
    props.convert_bpart_f(e, bparts + delta + i, &temp_f[i * dim]);
}

/**
 * @brief Mapper function to copy #bpart into a buffer of ints using a
 * conversion function.
 */
void io_convert_bpart_i_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct bpart* restrict bparts = props.bparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  int* restrict temp_i = (int*)temp;
  const ptrdiff_t delta = (temp_i - props.start_temp_i) / dim;

  for (int i = 0; i < N; i++)
    props.convert_bpart_i(e, bparts + delta + i, &temp_i[i * dim]);
}

/**
 * @brief Mapper function to copy #bpart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_bpart_d_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct bpart* restrict bparts = props.bparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  double* restrict temp_d = (double*)temp;
  const ptrdiff_t delta = (temp_d - props.start_temp_d) / dim;

  for (int i = 0; i < N; i++)
    props.convert_bpart_d(e, bparts + delta + i, &temp_d[i * dim]);
}

/**
 * @brief Mapper function to copy #bpart into a buffer of doubles using a
 * conversion function.
 */
void io_convert_bpart_l_mapper(void* restrict temp, int N,
                               void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct bpart* restrict bparts = props.bparts;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  long long* restrict temp_l = (long long*)temp;
  const ptrdiff_t delta = (temp_l - props.start_temp_l) / dim;

  for (int i = 0; i < N; i++)
    props.convert_bpart_l(e, bparts + delta + i, &temp_l[i * dim]);
}

/**
 * @brief Mapper function to copy #sink into a buffer of floats using a
 * conversion function.
 */
void io_convert_sink_f_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct sink* restrict sinks = props.sinks;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  float* restrict temp_f = (float*)temp;
  const ptrdiff_t delta = (temp_f - props.start_temp_f) / dim;

  for (int i = 0; i < N; i++)
    props.convert_sink_f(e, sinks + delta + i, &temp_f[i * dim]);
}

/**
 * @brief Mapper function to copy #sink into a buffer of ints using a
 * conversion function.
 */
void io_convert_sink_i_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct sink* restrict sinks = props.sinks;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  int* restrict temp_i = (int*)temp;
  const ptrdiff_t delta = (temp_i - props.start_temp_i) / dim;

  for (int i = 0; i < N; i++)
    props.convert_sink_i(e, sinks + delta + i, &temp_i[i * dim]);
}

/**
 * @brief Mapper function to copy #sink into a buffer of doubles using a
 * conversion function.
 */
void io_convert_sink_d_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct sink* restrict sinks = props.sinks;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  double* restrict temp_d = (double*)temp;
  const ptrdiff_t delta = (temp_d - props.start_temp_d) / dim;

  for (int i = 0; i < N; i++)
    props.convert_sink_d(e, sinks + delta + i, &temp_d[i * dim]);
}

/**
 * @brief Mapper function to copy #sink into a buffer of doubles using a
 * conversion function.
 */
void io_convert_sink_l_mapper(void* restrict temp, int N,
                              void* restrict extra_data) {

  const struct io_props props = *((const struct io_props*)extra_data);
  const struct sink* restrict sinks = props.sinks;
  const struct engine* e = props.e;
  const size_t dim = props.dimension;

  /* How far are we with this chunk? */
  long long* restrict temp_l = (long long*)temp;
  const ptrdiff_t delta = (temp_l - props.start_temp_l) / dim;

  for (int i = 0; i < N; i++)
    props.convert_sink_l(e, sinks + delta + i, &temp_l[i * dim]);
}

/**
 * @brief Copy the particle data into a temporary buffer ready for i/o.
 *
 * @param temp The buffer to be filled. Must be allocated and aligned properly.
 * @param e The #engine.
 * @param props The #io_props corresponding to the particle field we are
 * copying.
 * @param N The number of particles to copy
 * @param internal_units The system of units used internally.
 * @param snapshot_units The system of units used for the snapshots.
 */
void io_copy_temp_buffer(void* temp, const struct engine* e,
                         struct io_props props, size_t N,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;
  const size_t num_elements = N * props.dimension;

  /* Copy particle data to temporary buffer */
  if (props.conversion == 0) { /* No conversion */

    /* Prepare some parameters */
    char* temp_c = (char*)temp;
    props.start_temp_c = temp_c;

    /* Copy the whole thing into a buffer */
    threadpool_map((struct threadpool*)&e->threadpool, io_copy_mapper, temp_c,
                   N, copySize, threadpool_auto_chunk_size, (void*)&props);

  } else { /* Converting particle to data */

    if (props.convert_part_f != NULL) {

      /* Prepare some parameters */
      float* temp_f = (float*)temp;
      props.start_temp_f = (float*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_part_f_mapper, temp_f, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_part_i != NULL) {

      /* Prepare some parameters */
      int* temp_i = (int*)temp;
      props.start_temp_i = (int*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_part_i_mapper, temp_i, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_part_d != NULL) {

      /* Prepare some parameters */
      double* temp_d = (double*)temp;
      props.start_temp_d = (double*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_part_d_mapper, temp_d, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_part_l != NULL) {

      /* Prepare some parameters */
      long long* temp_l = (long long*)temp;
      props.start_temp_l = (long long*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_part_l_mapper, temp_l, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_gpart_f != NULL) {

      /* Prepare some parameters */
      float* temp_f = (float*)temp;
      props.start_temp_f = (float*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_gpart_f_mapper, temp_f, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_gpart_i != NULL) {

      /* Prepare some parameters */
      int* temp_i = (int*)temp;
      props.start_temp_i = (int*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_gpart_i_mapper, temp_i, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_gpart_d != NULL) {

      /* Prepare some parameters */
      double* temp_d = (double*)temp;
      props.start_temp_d = (double*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_gpart_d_mapper, temp_d, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_gpart_l != NULL) {

      /* Prepare some parameters */
      long long* temp_l = (long long*)temp;
      props.start_temp_l = (long long*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_gpart_l_mapper, temp_l, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_spart_f != NULL) {

      /* Prepare some parameters */
      float* temp_f = (float*)temp;
      props.start_temp_f = (float*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_spart_f_mapper, temp_f, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_spart_i != NULL) {

      /* Prepare some parameters */
      int* temp_i = (int*)temp;
      props.start_temp_i = (int*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_spart_i_mapper, temp_i, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_spart_d != NULL) {

      /* Prepare some parameters */
      double* temp_d = (double*)temp;
      props.start_temp_d = (double*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_spart_d_mapper, temp_d, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_spart_l != NULL) {

      /* Prepare some parameters */
      long long* temp_l = (long long*)temp;
      props.start_temp_l = (long long*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_spart_l_mapper, temp_l, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_sink_f != NULL) {

      /* Prepare some parameters */
      float* temp_f = (float*)temp;
      props.start_temp_f = (float*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_sink_f_mapper, temp_f, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_sink_i != NULL) {

      /* Prepare some parameters */
      int* temp_i = (int*)temp;
      props.start_temp_i = (int*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_sink_i_mapper, temp_i, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_sink_d != NULL) {

      /* Prepare some parameters */
      double* temp_d = (double*)temp;
      props.start_temp_d = (double*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_sink_d_mapper, temp_d, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_sink_l != NULL) {

      /* Prepare some parameters */
      long long* temp_l = (long long*)temp;
      props.start_temp_l = (long long*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_sink_l_mapper, temp_l, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_bpart_f != NULL) {

      /* Prepare some parameters */
      float* temp_f = (float*)temp;
      props.start_temp_f = (float*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_bpart_f_mapper, temp_f, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_bpart_i != NULL) {

      /* Prepare some parameters */
      int* temp_i = (int*)temp;
      props.start_temp_i = (int*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_bpart_i_mapper, temp_i, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_bpart_d != NULL) {

      /* Prepare some parameters */
      double* temp_d = (double*)temp;
      props.start_temp_d = (double*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_bpart_d_mapper, temp_d, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else if (props.convert_bpart_l != NULL) {

      /* Prepare some parameters */
      long long* temp_l = (long long*)temp;
      props.start_temp_l = (long long*)temp;
      props.e = e;

      /* Copy the whole thing into a buffer */
      threadpool_map((struct threadpool*)&e->threadpool,
                     io_convert_bpart_l_mapper, temp_l, N, copySize,
                     threadpool_auto_chunk_size, (void*)&props);

    } else {
      error("Missing conversion function");
    }
  }

  /* Unit conversion if necessary */
  const double factor =
      units_conversion_factor(internal_units, snapshot_units, props.units);
  if (factor != 1.) {

    /* message("Converting ! factor=%e", factor); */

    if (io_is_double_precision(props.type)) {
      swift_declare_aligned_ptr(double, temp_d, (double*)temp,
                                IO_BUFFER_ALIGNMENT);
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= factor;
    } else {
      swift_declare_aligned_ptr(float, temp_f, (float*)temp,
                                IO_BUFFER_ALIGNMENT);
      for (size_t i = 0; i < num_elements; ++i) temp_f[i] *= factor;
    }
  }
}
