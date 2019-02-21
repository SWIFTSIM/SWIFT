/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_INTERPOL_EAGLE_H
#define SWIFT_INTERPOL_EAGLE_H

/**
 * @file src/cooling/EAGLE/interpolate.h
 * @brief Interpolation functions for EAGLE tables
 */

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "align.h"
#include "error.h"
#include "inline.h"

/**
 * @brief Returns the 1d index of element with 2d indices x,y
 * from a flattened 2d array in row major order
 *
 * @param x, y Indices of element of interest
 * @param Nx, Ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_2d(const int x,
                                                             const int y,
                                                             const int Nx,
                                                             const int Ny) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
#endif
  return x * Ny + y;
}

/**
 * @brief Returns the 1d index of element with 3d indices x,y,z
 * from a flattened 3d array in row major order
 *
 * @param x, y, z Indices of element of interest
 * @param Nx, Ny, Nz Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_3d(
    const int x, const int y, const int z, const int Nx, const int Ny,
    const int Nz) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
#endif
  return x * Ny * Nz + y * Nz + z;
}

/**
 * @brief Returns the 1d index of element with 4d indices x,y,z,w
 * from a flattened 4d array in row major order
 *
 * @param x, y, z, w Indices of element of interest
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_4d(
    const int x, const int y, const int z, const int w, const int Nx,
    const int Ny, const int Nz, const int Nw) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
  assert(w < Nw);
#endif
  return x * Ny * Nz * Nw + y * Nz * Nw + z * Nw + w;
}

/**
 * @brief Finds the index of a value in a table and compute delta to nearest
 * element.
 *
 * This function assumes the table is monotonically increasing with a constant
 * difference between adjacent values.
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * We use a small epsilon of 1e-4 to avoid out-of-range accesses due to
 * rounding errors.
 *
 * @param table The table to search in.
 * @param size The number of elements in the table.
 * @param x The value to search for.
 * @param i (return) The index in the table of the element.
 * @param *dx (return) The difference between x and table[i]
 */
__attribute__((always_inline)) INLINE void get_index_1d(
    const float *restrict table, const int size, const float x, int *i,
    float *restrict dx) {

  /* Small epsilon to avoid rounding issues leading to out-of-bound
   * access when using the indices later to read data from the tables. */
  const float epsilon = 1e-4f;

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Distance between elements in the array */
  const float delta = (size - 1) / (table[size - 1] - table[0]);

  if (x < table[0] + epsilon) {
    /* We are below the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x < table[size - 1] - epsilon) {
    /* Normal case */
    *i = (x - table[0]) * delta;

#ifdef SWIFT_DEBUG_CHECKS
    if (*i > size || *i < 0) {
      error(
          "trying to get index for value outside table range. Table size: %d, "
          "calculated index: %d, value: %.5e, table[0]: %.5e, grid size: %.5e",
          size, *i, x, table[0], delta);
    }
#endif

    *dx = (x - table[*i]) * delta;
  } else {
    /* We are after the last element */
    *i = size - 2;
    *dx = 1.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (*dx < -0.001f || *dx > 1.001f) error("Invalid distance found dx=%e", *dx);
#endif
}

/**
 * @brief Interpolate a flattened 2D table at a given position.
 *
 * This function uses linear interpolation along each axis. It also
 * assumes that the table is aligned on SWIFT_STRUCT_ALIGNMENT.
 *
 * @param table The 2D table to interpolate.
 * @param xi, yi Indices of element of interest.
 * @param Nx, Ny Sizes of array dimensions.
 * @param dx, dy Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_2d(
    const float *table, const int xi, const int yi, const float dx,
    const float dy, const int Nx, const int Ny) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx < -0.001f || dx > 1.001f) error("Invalid dx=%e", dx);
  if (dy < -0.001f || dy > 1.001f) error("Invalid dy=%e", dy);
#endif

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Linear interpolation along each axis. We read the table 2^2=4 times */
  float result = tx * ty * table[row_major_index_2d(xi + 0, yi + 0, Nx, Ny)];

  result += tx * dy * table[row_major_index_2d(xi + 0, yi + 1, Nx, Ny)];
  result += dx * ty * table[row_major_index_2d(xi + 1, yi + 0, Nx, Ny)];

  result += dx * dy * table[row_major_index_2d(xi + 1, yi + 1, Nx, Ny)];

  return result;
}

/**
 * @brief Interpolate a flattened 3D table at a given position.
 *
 * This function uses linear interpolation along each axis. It also
 * assumes that the table is aligned on SWIFT_STRUCT_ALIGNMENT.
 *
 * @param table The 3D table to interpolate.
 * @param xi, yi, zi Indices of element of interest.
 * @param Nx, Ny, Nz Sizes of array dimensions.
 * @param dx, dy, dz Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_3d(
    const float *table, const int xi, const int yi, const int zi,
    const float dx, const float dy, const float dz, const int Nx, const int Ny,
    const int Nz) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx < -0.001f || dx > 1.001f) error("Invalid dx=%e", dx);
  if (dy < -0.001f || dy > 1.001f) error("Invalid dy=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dz=%e", dz);
#endif

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;

  /* Indicate that the whole array is aligned on page boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Linear interpolation along each axis. We read the table 2^3=8 times */
  float result = tx * ty * tz *
                 table[row_major_index_3d(xi + 0, yi + 0, zi + 0, Nx, Ny, Nz)];

  result += tx * ty * dz *
            table[row_major_index_3d(xi + 0, yi + 0, zi + 1, Nx, Ny, Nz)];
  result += tx * dy * tz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 0, Nx, Ny, Nz)];
  result += dx * ty * tz *
            table[row_major_index_3d(xi + 1, yi + 0, zi + 0, Nx, Ny, Nz)];

  result += tx * dy * dz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 1, Nx, Ny, Nz)];
  result += dx * ty * dz *
            table[row_major_index_3d(xi + 1, yi + 0, zi + 1, Nx, Ny, Nz)];
  result += dx * dy * tz *
            table[row_major_index_3d(xi + 1, yi + 1, zi + 0, Nx, Ny, Nz)];

  result += dx * dy * dz *
            table[row_major_index_3d(xi + 1, yi + 1, zi + 1, Nx, Ny, Nz)];

  return result;
}

/**
 * @brief Interpolate a flattened 3D table at a given position but avoid the
 * x-dimension.
 *
 * This function uses linear interpolation along each axis.
 * We look at the xi coordoniate but do not interpolate around it. We just
 * interpolate the remaining 2 dimensions.
 * The function also assumes that the table is aligned on
 * SWIFT_STRUCT_ALIGNMENT.
 *
 * @param table The 3D table to interpolate.
 * @param xi, yi, zi Indices of element of interest.
 * @param Nx, Ny, Nz Sizes of array dimensions.
 * @param dx, dy, dz Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_3d_no_x(
    const float *table, const int xi, const int yi, const int zi,
    const float dx, const float dy, const float dz, const int Nx, const int Ny,
    const int Nz) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx != 0.f) error("Attempting to interpolate along x!");
  if (dy < -0.001f || dy > 1.001f) error("Invalid dy=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dz=%e", dz);
#endif

  const float tx = 1.f;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;

  /* Indicate that the whole array is aligned on page boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Linear interpolation along each axis. We read the table 2^2=4 times */
  /* Note that we intentionally kept the table access along the axis where */
  /* we do not interpolate as comments in the code to allow readers to */
  /* understand what is going on. */
  float result = tx * ty * tz *
                 table[row_major_index_3d(xi + 0, yi + 0, zi + 0, Nx, Ny, Nz)];

  result += tx * ty * dz *
            table[row_major_index_3d(xi + 0, yi + 0, zi + 1, Nx, Ny, Nz)];
  result += tx * dy * tz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 0, Nx, Ny, Nz)];
  /* result += dx * ty * tz * */
  /*           table[row_major_index_3d(xi + 1, yi + 0, zi + 0, Nx, Ny, Nz)]; */

  result += tx * dy * dz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 1, Nx, Ny, Nz)];
  /* result += dx * ty * dz * */
  /*           table[row_major_index_3d(xi + 1, yi + 0, zi + 1, Nx, Ny, Nz)]; */
  /* result += dx * dy * tz * */
  /*           table[row_major_index_3d(xi + 1, yi + 1, zi + 0, Nx, Ny, Nz)]; */

  /* result += dx * dy * dz * */
  /*           table[row_major_index_3d(xi + 1, yi + 1, zi + 1, Nx, Ny, Nz)]; */

  return result;
}

/**
 * @brief Interpolate a flattened 4D table at a given position.
 *
 * This function uses linear interpolation along each axis. It also
 * assumes that the table is aligned on SWIFT_STRUCT_ALIGNMENT.
 *
 * @param table The 4D table to interpolate.
 * @param xi, yi, zi, wi Indices of element of interest.
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions.
 * @param dx, dy, dz, dw Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_4d(
    const float *table, const int xi, const int yi, const int zi, const int wi,
    const float dx, const float dy, const float dz, const float dw,
    const int Nx, const int Ny, const int Nz, const int Nw) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx < -0.001f || dx > 1.001f) error("Invalid dx=%e", dx);
  if (dy < -0.001f || dy > 1.001f) error("Invalid dy=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dz=%e", dz);
  if (dw < -0.001f || dw > 1.001f) error("Invalid dw=%e", dw);
#endif

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;
  const float tw = 1.f - dw;

  /* Indicate that the whole array is aligned on page boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Linear interpolation along each axis. We read the table 2^4=16 times */
  float result =
      tx * ty * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      tx * ty * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * ty * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * tz * tw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      tx * ty * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * tz * dw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * dz * tw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * dy * tz * tw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      dx * dy * dz * tw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * dy * tz * dw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * dz * dw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];

  result +=
      dx * dy * dz * dw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];

  return result;
}

/**
 * @brief Interpolate a flattened 4D table at a given position but avoid the
 * x-dimension.
 *
 * This function uses linear interpolation along each axis.
 * We look at the xi coordoniate but do not interpolate around it. We just
 * interpolate the remaining 3 dimensions.
 * The function also assumes that the table is aligned on
 * SWIFT_STRUCT_ALIGNMENT.
 *
 * @param table The 4D table to interpolate.
 * @param xi, yi, zi, wi Indices of element of interest.
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions.
 * @param dx, dy, dz, dw Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_4d_no_x(
    const float *table, const int xi, const int yi, const int zi, const int wi,
    const float dx, const float dy, const float dz, const float dw,
    const int Nx, const int Ny, const int Nz, const int Nw) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx != 0.f) error("Attempting to interpolate along x!");
  if (dy < -0.001f || dy > 1.001f) error("Invalid dy=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dz=%e", dz);
  if (dw < -0.001f || dw > 1.001f) error("Invalid dw=%e", dw);
#endif

  const float tx = 1.f;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;
  const float tw = 1.f - dw;

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Linear interpolation along each axis. We read the table 2^3=8 times */
  /* Note that we intentionally kept the table access along the axis where */
  /* we do not interpolate as comments in the code to allow readers to */
  /* understand what is going on. */
  float result =
      tx * ty * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      tx * ty * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * ty * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];
  /* result += */
  /*     dx * ty * tz * tw * */
  /*     table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz,
   * Nw)]; */

  result +=
      tx * ty * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  /* result += */
  /*     dx * ty * tz * dw * */
  /*     table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz,
   * Nw)]; */
  result +=
      tx * dy * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  /* result += */
  /*     dx * ty * dz * tw * */
  /*     table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz,
   * Nw)]; */
  /* result += */
  /*     dx * dy * tz * tw * */
  /*     table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, */
  /* Nw)]; */

  /* result += */
  /*     dx * dy * dz * tw * */
  /*     table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, */
  /* Nw)]; */
  /* result += */
  /*     dx * dy * tz * dw * */
  /*     table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, */
  /* Nw)]; */
  /* result += */
  /*     dx * ty * dz * dw * */
  /*     table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz,
   * Nw)]; */
  result +=
      tx * dy * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];

  /* result += */
  /*     dx * dy * dz * dw * */
  /*     table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, */
  /* Nw)]; */

  return result;
}

#endif
