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
#ifndef SWIFT_GEAR_INTERPOLATION_H
#define SWIFT_GEAR_INTERPOLATION_H

#include "error.h"
#include "exp10.h"
#include "inline.h"
#include "minmax.h"

#include <stddef.h>

/**
 * @brief Type of boundary condition available.
 */
enum interpolate_boundary_condition {
  /* No extrapolation => raise errors */
  boundary_condition_error,

  /* Zero as boundary conditions */
  boundary_condition_zero,

  /* Zero (left boundary) and constant (right boundary) boundary conditions */
  boundary_condition_zero_const,

  /* constant boundary conditions */
  boundary_condition_const,
};

/*****************************************************************************/
/* Interpolation for float data */
/*****************************************************************************/

/**
 * @brief Structure for the interpolation, float version.
 */
struct interpolation_1d {
  /* Data to interpolate */
  float *data;

  /* Minimal x */
  float xmin;

  /* Step size between x points */
  float dx;

  /* Number of element in the data */
  int N;

  /* Type of boundary conditions. */
  enum interpolate_boundary_condition boundary_condition;
};

/**
 * @brief Initialize the #interpolation_1d.
 *
 * @param interp The #interpolation_1d.
 * @param xmin Minimal value of x (in log).
 * @param xmax Maximal value of x (in log).
 * @param N Requested number of values.
 * @param log_data_xmin The minimal value of the data (in log).
 * @param step_size The size of the x steps (in log).
 * @param N_data The number of element in the data.
 * @param data The data to interpolate (y).
 * @param N The number of element in data.
 * @param boundary_condition The type of #interpolate_boundary_condition.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_init(
    struct interpolation_1d *interp, float xmin, float xmax, int N,
    float log_data_xmin, float step_size, int N_data, const float *data,
    enum interpolate_boundary_condition boundary_condition) {

  /* Save the variables */
  interp->N = N;
  interp->xmin = xmin;
  interp->dx = (xmax - xmin) / (N - 1.f);
  interp->boundary_condition = boundary_condition;

  /* Allocate the memory */
  interp->data = malloc(sizeof(float) * N);
  if (interp->data == NULL)
    error("Failed to allocate memory for the interpolation");

  /* Interpolate the data */
  for (int i = 0; i < N; i++) {
    const float log_x = xmin + i * interp->dx;
    const float x_j = (log_x - log_data_xmin) / step_size;

    /* Check boundaries */
    if (x_j < 0) {
      switch (boundary_condition) {
        case boundary_condition_error:
          error("Cannot extrapolate");
          break;
        case boundary_condition_zero:
          interp->data[i] = 0;
          break;
        case boundary_condition_zero_const:
          interp->data[i] = 0;
          break;
        case boundary_condition_const:
          interp->data[i] = data[0];
          break;
        default:
          error("Interpolation type not implemented");
      }
      continue;
    } else if (x_j >= N_data) {
      switch (boundary_condition) {
        case boundary_condition_error:
          error("Cannot extrapolate");
          break;
        case boundary_condition_zero:
          interp->data[i] = 0;
          break;
        case boundary_condition_zero_const:
        case boundary_condition_const:
          interp->data[i] = interp->data[i - 1];
          break;
        default:
          error("Interpolation type not implemented");
      }
      continue;
    }

    /* Interpolate i */
    int j = (int)x_j;

    /* Handle the edge case where x_j is exactly at or very close to N_data - 1
     */
    if (j >= N_data - 1) {
      interp->data[i] = data[N_data - 1];
    } else {
      const float f = x_j - (float)j;
      interp->data[i] = (1.f - f) * data[j] + f * data[j + 1];
    }
  }
}

/**
 * @brief Interpolate the data.
 *
 * @param interp The #interpolation_1d.
 * @param x The x value where to interpolate.
 *
 * @return The interpolated value y.
 */
__attribute__((always_inline)) static INLINE float interpolate_1d(
    const struct interpolation_1d *interp, float x) {

  /* Find indice */
  const float i = (x - interp->xmin) / interp->dx;
  const int idx = i;
  const float dx = i - idx;

  /* Should we extrapolate? */
  if (i < 0) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
      case boundary_condition_zero_const:
        return 0;
      case boundary_condition_const:
        return interp->data[0];
      default:
        error("Interpolation type not implemented");
    }
  } else if (i >= interp->N - 1) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
        return 0;
      case boundary_condition_zero_const:
      case boundary_condition_const:
        return interp->data[interp->N - 1];
      default:
        error("Interpolation type not implemented");
    }
  }

  /* interpolate */
  return interp->data[idx] * (1. - dx) + interp->data[idx + 1] * dx;
}

/**
 * @brief Print the data.
 *
 * @param interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_print(
    const struct interpolation_1d *interp) {

  message("Interpolation between %g and %g", interp->xmin,
          interp->xmin + interp->dx * interp->N);

  message("Contains %i values and use the boundary condition %i", interp->N,
          interp->boundary_condition);

  /* Print values */
  for (int i = 0; i < interp->N; i++) {
    float x = interp->xmin + i * interp->dx;
    message("%.2g: %g", x, interp->data[i]);
  }
}

/**
 * @brief Cleanup the #interpolation_1d structure.
 *
 * @param interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_free(
    struct interpolation_1d *interp) {

  /* Free the allocated memory */
  free(interp->data);
  interp->data = NULL;
}

/**
 * @brief Zero pointers in an #interpolation_1d struct, so a struct that was
 * never (or is not yet) initialized can still be safely passed to
 * #interpolate_1d_free.
 *
 * @param interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_zero_pointers(
    struct interpolation_1d *interp) {
  interp->data = NULL;
  interp->xmin = 0.0;
  interp->dx = 0.0;
  interp->N = 0;
  interp->boundary_condition = boundary_condition_error;
}

/**************************** Interpolation 2D ********************************/

/**
 * @brief Structure for the 2D interpolation.
 */
struct interpolation_2d {
  /* Data to interpolate */
  float *data;

  /* The log-x nodes this table is stored on: #Nx strictly increasing
     values, owned by this struct. #interpolate_2d brackets these at query
     time when #is_uniform_x is 0, so nothing there depends on their
     spacing being regular. When #is_uniform_x is 1 the nodes are a
     materialised convenience for printing and bracketing helpers: the
     index arithmetic then uses #xmin and #dx, which differ from
     differencing two adjacent nodes by of order |xmin| / dx * 2^-24. */
  float *xvals;

  /* First log-x node. */
  float xmin;

  /* Step between log-x nodes. Meaningful only when #is_uniform_x is 1; 0
     otherwise, since a node list with irregular spacing has no step. */
  float dx;

  /* Minimal y */
  float ymin;

  /* Step size between y points */
  float dy;

  /* Number of element in the x direction of the data */
  int Nx;

  /* Number of element in the y direction of the data */
  int Ny;

  /* Is the x axis exactly #xmin + i * #dx? Selects the index arithmetic. */
  int is_uniform_x;

  /* Type of boundary condition applied when x is out of the data range. */
  enum interpolate_boundary_condition boundary_condition_x;

  /* Type of boundary condition applied when y is out of the data range. */
  enum interpolate_boundary_condition boundary_condition_y;
};

/* #interpolation_2d.xvals falls inside an #interpolation_1d's footprint, so
   wherever the two share a union #interpolate_2d_free() would free() a
   pointer synthesised from float bit patterns if it ever ran on a live 1D
   table: every free site must select the member through its own
   dimensionality flag (see radiation_clean()). */
_Static_assert(offsetof(struct interpolation_2d, xvals) <
                   sizeof(struct interpolation_1d),
               "interpolation_2d.xvals no longer overlaps interpolation_1d; "
               "re-derive the free() dispatch rule.");

/**
 * @brief Does this #interpolate_boundary_condition mean "return zero" on the
 * lower-bound (index < 0) side?
 *
 * boundary_condition_zero_const is asymmetric: zero below, constant above
 * (see #interpolate_boundary_condition), so it groups with
 * boundary_condition_zero here but not in
 * #interpolate_boundary_is_zero_above().
 */
__attribute__((always_inline)) static INLINE int
interpolate_boundary_is_zero_below(enum interpolate_boundary_condition bc) {
  return bc == boundary_condition_zero || bc == boundary_condition_zero_const;
}

/**
 * @brief Does this #interpolate_boundary_condition mean "return zero" on the
 * upper-bound (index >= N) side?
 */
__attribute__((always_inline)) static INLINE int
interpolate_boundary_is_zero_above(enum interpolate_boundary_condition bc) {
  return bc == boundary_condition_zero;
}

/**
 * @brief Raise an error if this #interpolate_boundary_condition is
 * boundary_condition_error; a no-op otherwise.
 */
__attribute__((always_inline)) static INLINE void
interpolate_boundary_check_error(enum interpolate_boundary_condition bc) {
  if (bc == boundary_condition_error) error("Cannot extrapolate");
}

/**
 * @brief Decide whether an out-of-range 2D interpolation query evaluates to
 * zero, and raise an error for any out-of-range axis whose boundary
 * condition is boundary_condition_error.
 *
 * Must be called only once the caller has established that x_raw or y_raw
 * is out of range. x_raw/y_raw must be the raw (untruncated) index-space
 * coordinates: truncating first would let a value in (-1, 0) wrongly read
 * as "in range". A zero policy on the out-of-range axis wins over a const
 * policy on the other axis.
 *
 * @param boundary_condition_x The #interpolate_boundary_condition applied
 * when x is out of range.
 * @param boundary_condition_y The #interpolate_boundary_condition applied
 * when y is out of range.
 * @param x_raw The raw (untruncated) x coordinate in index space.
 * @param y_raw The raw (untruncated) y coordinate in index space.
 * @param nx The number of indices along x.
 * @param ny The number of indices along y.
 *
 * @return 1 if the query should evaluate to zero, 0 otherwise.
 */
__attribute__((always_inline)) static INLINE int
interpolate_2d_boundary_is_zero(
    enum interpolate_boundary_condition boundary_condition_x,
    enum interpolate_boundary_condition boundary_condition_y, float x_raw,
    float y_raw, int nx, int ny) {

  const int x_below = x_raw < 0;
  const int x_above = !x_below && x_raw >= nx - 1;
  const int y_below = y_raw < 0;
  const int y_above = !y_below && y_raw >= ny - 1;

  if (x_below || x_above)
    interpolate_boundary_check_error(boundary_condition_x);
  if (y_below || y_above)
    interpolate_boundary_check_error(boundary_condition_y);

  return (x_below &&
          interpolate_boundary_is_zero_below(boundary_condition_x)) ||
         (x_above &&
          interpolate_boundary_is_zero_above(boundary_condition_x)) ||
         (y_below &&
          interpolate_boundary_is_zero_below(boundary_condition_y)) ||
         (y_above && interpolate_boundary_is_zero_above(boundary_condition_y));
}

/**
 * @brief Convert a log-x query into the index space of a list of log-x
 * nodes.
 *
 * Brackets the nodes and returns node index + the linear-in-log-x fraction
 * inside the bracketing interval. A value bit-equal to a node returns that
 * node's integer index with a zero fraction, which is what makes the
 * build-time resample copy a source row unchanged. A runtime query only
 * lands on a node to within the rounding of its own log: the build uses
 * log10f(Z) and the caller (#radiation_get_log_metallicity) a narrowed
 * log10((double)Z), which can differ by one ULP, so the query returns that
 * node's row to within one blend step of a neighbouring row rather than
 * exactly. The spacing of the nodes is never assumed to be regular.
 *
 * Out of range, the returned value stays below 0 or at/above N - 1 so the
 * caller's boundary branch (see #interpolate_2d_boundary_is_zero) fires.
 * The end intervals set the scale there; nothing is extrapolated by this
 * function itself.
 *
 * @param xv The log-x nodes, strictly increasing.
 * @param Nx The number of nodes in @p xv (at least 2).
 * @param log_x The x value where to interpolate in log.
 *
 * @return The (possibly out-of-range) fractional index along x.
 */
__attribute__((always_inline)) static INLINE float interpolate_index_in_nodes(
    const float *xv, int Nx, float log_x) {

  if (log_x < xv[0]) {
    const float step = xv[1] - xv[0];
    if (step <= 0.f) return -1.f;
    return (log_x - xv[0]) / step;
  }

  if (log_x >= xv[Nx - 1]) {
    const float step = xv[Nx - 1] - xv[Nx - 2];
    if (step <= 0.f) return (float)(Nx - 1);
    return (float)(Nx - 1) + (log_x - xv[Nx - 1]) / step;
  }

  /* xv[lo] <= log_x < xv[lo + 1], by bisection over the native nodes. */
  int lo = 0;
  int hi = Nx - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (log_x >= xv[mid])
      lo = mid;
    else
      hi = mid;
  }

  const float step = xv[lo + 1] - xv[lo];
  if (step <= 0.f) return (float)lo;
  return (float)lo + (log_x - xv[lo]) / step;
}

/**
 * @brief Convert a log-x query into the #interpolation_2d's own x index
 * space.
 *
 * A uniform axis divides by the step the table was built with. Recovering
 * that step by differencing two adjacent nodes instead would cost a
 * relative error of order |xmin| / dx * 2^-24, which reaches 9e-06 on the
 * stellar-wind metallicity axis (|xmin| / dx = 133).
 *
 * @param interp The #interpolation_2d.
 * @param log_x The x value where to interpolate in log.
 *
 * @return The (possibly out-of-range) fractional index along x.
 */
__attribute__((always_inline)) static INLINE float interpolate_2d_index_x(
    const struct interpolation_2d *interp, float log_x) {

  if (interp->is_uniform_x) return (log_x - interp->xmin) / interp->dx;

  if (interp->xvals == NULL || interp->Nx < 2)
    error("Cannot extrapolate: this interpolation table was never built");

  return interpolate_index_in_nodes(interp->xvals, interp->Nx, log_x);
}

/**
 * @brief Return the two x indices an #interpolation_2d query at @p log_x
 * blends between, clamped into the table.
 *
 * For a caller that needs the bracketing rows themselves rather than an
 * interpolated value. Both indices collapse to the nearest end row when the
 * query falls outside the table.
 *
 * @param interp The #interpolation_2d.
 * @param log_x The x value where to interpolate in log.
 * @param idx_lo (output) Lower bracketing x index.
 * @param idx_hi (output) Upper bracketing x index.
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_bracket_x(
    const struct interpolation_2d *interp, float log_x, int *idx_lo,
    int *idx_hi) {

  /* Nx = 0 would clamp both indices to -1 and hand the caller an
     out-of-bounds row. */
  if (interp->Nx < 1)
    error("Cannot bracket an interpolation table that was never built");

  const float x_raw = interpolate_2d_index_x(interp, log_x);
  const int clamped_low = max((int)x_raw, 0);
  const int lo = min(clamped_low, interp->Nx - 1);

  *idx_lo = lo;
  *idx_hi = min(lo + 1, interp->Nx - 1);
}

/**
 * @brief Build an #interpolation_2d from x source indices its caller has
 * already computed.
 *
 * The shared body of #interpolate_2d_init() and
 * #interpolate_2d_init_uniform_x(). Each of those owns the arithmetic that
 * turns an output node into a fractional index into @p data, so neither
 * has to express its x axis in the other's terms.
 *
 * @param interp The #interpolation_2d result, stored in swift.
 * @param x_index Fractional index into the source x axis for each of the
 * @p Nx output nodes.
 * @param log_x_out The output log-x nodes, strictly increasing.
 * @param Nx The number of nodes in @p log_x_out (at least 2).
 * @param is_uniform_x Are the @p log_x_out nodes exactly @p log_xmin +
 * i * @p log_dx?
 * @param log_xmin The first output log-x node.
 * @param log_dx The output log-x step, 0 when @p is_uniform_x is 0.
 * @param log_ymin Minimal value of y (in log).   Interpolation limits
 * @param log_ymax Maximal value of y (in log).   Interpolation limits
 * @param Ny Requested number of values in y axes.  Interpolation limits
 * @param log_data_ymin The minimal value of the data in y (in log).  Data
 * limits
 * @param log_step_size_y The size of the y steps (in log).   Data limits
 * @param N_data_x The number of element in the data x axis. Data limits
 * @param N_data_y The number of element in the data y axis. Data limits
 * @param data The data coming from hdf5 table to interpolate.
 * @param boundary_condition_x The #interpolate_boundary_condition applied
 * when x is out of the data range.
 * @param boundary_condition_y The #interpolate_boundary_condition applied
 * when y is out of the data range.
 */
__attribute__((always_inline)) static INLINE void
interpolate_2d_init_from_x_indices(
    struct interpolation_2d *interp, const float *x_index,
    const float *log_x_out, int Nx, int is_uniform_x, float log_xmin,
    float log_dx, float log_ymin, float log_ymax, int Ny, float log_data_ymin,
    float log_step_size_y, int N_data_x, int N_data_y, const double *data,
    enum interpolate_boundary_condition boundary_condition_x,
    enum interpolate_boundary_condition boundary_condition_y) {

  /* Save the variables */
  interp->Nx = Nx;
  interp->is_uniform_x = is_uniform_x;
  interp->xmin = log_xmin;
  interp->dx = log_dx;
  interp->boundary_condition_x = boundary_condition_x;
  interp->boundary_condition_y = boundary_condition_y;

  interp->Ny = Ny;
  interp->ymin = log_ymin;
  interp->dy = (log_ymax - log_ymin) / (Ny - 1.f);

  /* Allocate the memory */
  interp->data = malloc(sizeof(float) * Nx * Ny);
  if (interp->data == NULL) {
    error("Failed to allocate memory for the interpolation");
  }

  interp->xvals = malloc(sizeof(float) * Nx);
  if (interp->xvals == NULL) {
    error("Failed to allocate memory for the interpolation x axis");
  }
  for (int i = 0; i < Nx; i++) interp->xvals[i] = log_x_out[i];

  /* Interpolate the data */
  for (int i = 0; i < Nx; i++) {
    const float x_k = x_index[i];

    for (int j = 0; j < Ny; j++) {
      const float log_y = log_ymin + j * interp->dy;
      const float y_k = (log_y - log_data_ymin) / log_step_size_y;

      /* Data indexes */
      const int idx = x_k;
      const float fx = x_k - idx;
      const int idy = y_k;
      const float fy = y_k - idy;

      int current_cell = i * Ny + j;
      if (current_cell >= Nx * Ny) {
        error("Index %d out of boundaries for interp->data", current_cell);
      }
      /* Extrapolate? See #interpolate_2d_boundary_is_zero(). */
      if (x_k < 0 || x_k >= N_data_x - 1 || y_k < 0 || y_k >= N_data_y - 1) {
        const int is_zero = interpolate_2d_boundary_is_zero(
            boundary_condition_x, boundary_condition_y, x_k, y_k, N_data_x,
            N_data_y);

        if (is_zero) {
          interp->data[current_cell] = 0;
        } else {
          const int midx = max(idx, 0);
          const int midy = max(idy, 0);
          const int row = min(midx, N_data_x - 1);
          const int col = min(midy, N_data_y - 1);
          const int cell_to_get = row * N_data_y + col;
          if (cell_to_get >= N_data_x * N_data_y) {
            error(
                "Index row=%d col=%d is out of boundary for the target data "
                "which has dimension row=%d col=%d",
                row, col, N_data_x, N_data_y);
          }
          interp->data[current_cell] = data[cell_to_get];
        }
        continue;
      }

      /* Interpolate data[i][j] <=> data[i * Ny + j] */
      const float fx1 = data[idx * N_data_y + idy] * (1. - fx) +
                        data[(idx + 1) * N_data_y + idy] * fx;
      const float fx2 = data[idx * N_data_y + idy + 1] * (1. - fx) +
                        data[(idx + 1) * N_data_y + idy + 1] * fx;
      interp->data[current_cell] = fx1 * (1. - fy) + fx2 * fy;
    }
  }
}

/**
 * @brief Initialize the #interpolation_2d.
 *
 * Resamples @p data onto the @p Nx output x nodes @p log_x_out and onto
 * @p Ny points uniform in log-y, and keeps @p log_x_out so #interpolate_2d
 * brackets it at query time. The x axis is a list of nodes both on input
 * and on output, so how the source file describes its own x spacing never
 * enters the arithmetic.
 *
 * Passing @p log_data_x itself as @p log_x_out keeps the x axis exact:
 * every output row then lands on one source row with a zero blend
 * fraction. That is the only choice available to a source axis whose node
 * spacings share no common divisor, since the narrowest of them can be
 * finer than any practical uniform step.
 *
 * @param interp The #interpolation_2d result, stored in swift.
 * @param log_data_x The data's own log-x nodes, strictly increasing.
 * @param N_data_x The number of nodes in @p log_data_x (at least 2).
 * @param log_x_out The output log-x nodes, strictly increasing.
 * @param Nx The number of nodes in @p log_x_out (at least 2).
 * @param log_ymin Minimal value of y (in log).   Interpolation limits
 * @param log_ymax Maximal value of y (in log).   Interpolation limits
 * @param Ny Requested number of values in y axes.  Interpolation limits
 * @param log_data_ymin The minimal value of the data in y (in log).  Data
 * limits
 * @param log_step_size_y The size of the y steps (in log).   Data limits
 * @param N_data_y The number of element in the data y axis. Data limits
 * @param data The data coming from hdf5 table to interpolate.
 * @param boundary_condition_x The #interpolate_boundary_condition applied
 * when x is out of the data range.
 * @param boundary_condition_y The #interpolate_boundary_condition applied
 * when y is out of the data range.
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_init(
    struct interpolation_2d *interp, const float *log_data_x, int N_data_x,
    const float *log_x_out, int Nx, float log_ymin, float log_ymax, int Ny,
    float log_data_ymin, float log_step_size_y, int N_data_y,
    const double *data,
    enum interpolate_boundary_condition boundary_condition_x,
    enum interpolate_boundary_condition boundary_condition_y) {

  if (N_data_x < 2)
    error("An interpolation source x axis needs at least 2 nodes, got %d",
          N_data_x);
  if (Nx < 2)
    error("An interpolation output x axis needs at least 2 nodes, got %d", Nx);

  for (int i = 1; i < N_data_x; i++) {
    if (!(log_data_x[i] > log_data_x[i - 1]))
      error(
          "Interpolation source x nodes are not strictly increasing at "
          "index %d",
          i);
  }
  /* The query-time bisection depends on this one too. */
  for (int i = 1; i < Nx; i++) {
    if (!(log_x_out[i] > log_x_out[i - 1]))
      error(
          "Interpolation output x nodes are not strictly increasing at "
          "index %d",
          i);
  }

  float *x_index = (float *)malloc(sizeof(float) * Nx);
  if (x_index == NULL)
    error("Failed to allocate memory for the interpolation x indices");

  for (int i = 0; i < Nx; i++)
    x_index[i] = interpolate_index_in_nodes(log_data_x, N_data_x, log_x_out[i]);

  interpolate_2d_init_from_x_indices(
      interp, x_index, log_x_out, Nx, /*is_uniform_x=*/0, log_x_out[0],
      /*log_dx=*/0.f, log_ymin, log_ymax, Ny, log_data_ymin, log_step_size_y,
      N_data_x, N_data_y, data, boundary_condition_x, boundary_condition_y);

  free(x_index);
}

/**
 * @brief Initialize an #interpolation_2d from a source x axis and an output
 * x axis both described as a minimum/step/count triple.
 *
 * Both axes divide by the step they were given, so the resampled table and
 * every later query reproduce the arithmetic a regularly spaced axis has
 * always used. Going through the node list instead would recover each step
 * by differencing two adjacent nodes and lose of order
 * |log_xmin| / log_step_size_x * 2^-24 of relative precision.
 *
 * @param interp The #interpolation_2d result, stored in swift.
 * @param log_xmin Minimal value of x (in log).  Interpolation limits
 * @param log_xmax Maximal value of x (in log).  Interpolation limits
 * @param Nx Requested number of values in x axes.  Interpolation limits
 * @param log_ymin Minimal value of y (in log).   Interpolation limits
 * @param log_ymax Maximal value of y (in log).   Interpolation limits
 * @param Ny Requested number of values in y axes.  Interpolation limits
 * @param log_data_xmin The minimal value of the data in x (in log).  Data
 * limits
 * @param log_data_ymin The minimal value of the data in y (in log).  Data
 * limits
 * @param log_step_size_x The size of the x steps (in log).   Data limits
 * @param log_step_size_y The size of the y steps (in log).   Data limits
 * @param N_data_x The number of element in the data x axis. Data limits
 * @param N_data_y The number of element in the data y axis. Data limits
 * @param data The data coming from hdf5 table to interpolate.
 * @param boundary_condition_x The #interpolate_boundary_condition applied
 * when x is out of the data range.
 * @param boundary_condition_y The #interpolate_boundary_condition applied
 * when y is out of the data range.
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_init_uniform_x(
    struct interpolation_2d *interp, float log_xmin, float log_xmax, int Nx,
    float log_ymin, float log_ymax, int Ny, float log_data_xmin,
    float log_data_ymin, float log_step_size_x, float log_step_size_y,
    int N_data_x, int N_data_y, const double *data,
    enum interpolate_boundary_condition boundary_condition_x,
    enum interpolate_boundary_condition boundary_condition_y) {

  if (N_data_x < 2 || Nx < 2)
    error(
        "A uniform interpolation x axis needs at least 2 nodes, got %d "
        "source and %d output",
        N_data_x, Nx);

  float *log_x_out = (float *)malloc(sizeof(float) * Nx);
  float *x_index = (float *)malloc(sizeof(float) * Nx);
  if (log_x_out == NULL || x_index == NULL)
    error("Failed to allocate memory for the interpolation x axis");

  const float dx = (log_xmax - log_xmin) / (Nx - 1.f);
  for (int i = 0; i < Nx; i++) {
    log_x_out[i] = log_xmin + i * dx;
    x_index[i] = (log_x_out[i] - log_data_xmin) / log_step_size_x;
  }

  interpolate_2d_init_from_x_indices(
      interp, x_index, log_x_out, Nx, /*is_uniform_x=*/1, log_xmin, dx,
      log_ymin, log_ymax, Ny, log_data_ymin, log_step_size_y, N_data_x,
      N_data_y, data, boundary_condition_x, boundary_condition_y);

  free(log_x_out);
  free(x_index);
}

/**
 * @brief Interpolate the data.
 *
 * @param interp The #interpolation_2d.
 * @param x The x value where to interpolate in log.
 * @param y The y value where to interpolate in log.
 *
 * @return The interpolated value.
 */
__attribute__((always_inline)) static INLINE double interpolate_2d(
    const struct interpolation_2d *interp, float log_x, float log_y) {

  /* Find indices */
  const float i = interpolate_2d_index_x(interp, log_x);
  const int idx = i;
  const float dx = i - idx;

  const int Nx = interp->Nx;
  const int Ny = interp->Ny;
  const int array_size = Nx * Ny;

  const float j = (log_y - interp->ymin) / interp->dy;
  const int idy = j;
  const float dy = j - idy;

  /* Extrapolate? See #interpolate_2d_boundary_is_zero(). */
  if (i < 0 || i >= Nx - 1 || j < 0 || j >= Ny - 1) {
    const int is_zero = interpolate_2d_boundary_is_zero(
        interp->boundary_condition_x, interp->boundary_condition_y, i, j, Nx,
        Ny);

    if (is_zero) {
#if defined(SWIFT_TEST_STELLAR_WIND)
      message(
          "interp->Nx=%d interp->Ny=%d interp->xmin=%g interp->ymin=%g "
          "interp->dx=%g interp->dy=%g idx=%d idy=%d "
          "out_of_boundary_type=zero",
          Nx, Ny, interp->xmin, interp->ymin, interp->dx, interp->dy, idx, idy);
#endif /* !defined SWIFT_TEST_STELLAR_WIND */
      return 0;
    }

    const int midx = max(idx, 0);
    const int midy = max(idy, 0);
    const int row = min(midx, Nx - 1);
    const int col = min(midy, Ny - 1);
    const int cell_to_get = row * Ny + col;
    if (cell_to_get >= array_size) {
      error("Index %d is out of boundary for the target data", cell_to_get);
    }
#if defined(SWIFT_TEST_STELLAR_WIND)
    message(
        "interp->Nx=%d interp->Ny=%d interp->xmin=%g interp->ymin=%g "
        "interp->dx=%g interp->dy=%g idx=%d idy=%d "
        "out_of_boundary_type=const cell_to_get=%d E[%d][%d]=%g",
        Nx, Ny, interp->xmin, interp->ymin, interp->dx, interp->dy, idx, idy,
        cell_to_get, row, col, exp10(interp->data[cell_to_get]));
#endif /* !defined SWIFT_TEST_STELLAR_WIND */
    return interp->data[cell_to_get];
  }

  /* interpolate */
  if ((idx * Ny + idy) >= array_size || ((idx + 1) * Ny + idy) >= array_size ||
      (idx * Ny + idy + 1) >= array_size ||
      ((idx + 1) * Ny + idy + 1) >= array_size) {
    error("Index is out of boundaries for the interpolation");
  }

#if defined(SWIFT_TEST_STELLAR_WIND)
  message(
      "interp->Nx=%d interp->Ny=%d interp->xmin=%g interp->ymin=%g "
      "interp->dx=%g interp->dy=%g idx=%d idy=%d out_of_boundary_type=none "
      "E[idx][idy]=%g E[idx][idy+1]=%g E[idx+1][idy]=%g E[idx+1][idy+1]=%g",
      Nx, Ny, interp->xmin, interp->ymin, interp->dx, interp->dy, idx, idy,
      exp10(interp->data[idx * Ny + idy]),
      exp10(interp->data[idx * Ny + idy + 1]),
      exp10(interp->data[(idx + 1) * Ny + idy]),
      exp10(interp->data[(idx + 1) * Ny + idy + 1]));
#endif /* !defined SWIFT_TEST_STELLAR_WIND */

  const float fx1 = interp->data[idx * Ny + idy] * (1. - dx) +
                    interp->data[(idx + 1) * Ny + idy] * dx;
  const float fx2 = interp->data[idx * Ny + idy + 1] * (1. - dx) +
                    interp->data[(idx + 1) * Ny + idy + 1] * dx;
  return fx1 * (1. - dy) + fx2 * dy;
}

/**
 * @brief Cleanup the #interpolation_2d structure.
 *
 * @param interp The #interpolation_2d.
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_free(
    struct interpolation_2d *interp) {

  /* Free the allocated memory */
  free(interp->data);
  interp->data = NULL;
  free(interp->xvals);
  interp->xvals = NULL;
}

/**
 * @brief zero pointers in interpolation_2d struct
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_zero_pointers(
    struct interpolation_2d *interp) {
  interp->data = NULL;
  interp->xvals = NULL;
  interp->Nx = 0;
  interp->Ny = 0;
  interp->is_uniform_x = 0;
  interp->xmin = 0.0;
  interp->dx = 0.0;
  interp->ymin = 0.0;
  interp->dy = 0.0;
  interp->boundary_condition_x = boundary_condition_error;
  interp->boundary_condition_y = boundary_condition_error;
}

#endif  // SWIFT_GEAR_INTERPOLATION_H
