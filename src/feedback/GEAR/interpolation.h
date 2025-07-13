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

/**
 * @brief Structure for the interpolation.
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
    const int j = x_j;
    const float f = x_j - j;
    interp->data[i] = (1. - f) * data[j] + f * data[j + 1];
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


////////////////////////////// Interpolation 2D /////////////////////////////////

/**
 * @brief Structure for the 2D interpolation.
 */
struct interpolation_2d {
  /* Data to interpolate */
  float *data;

  /* Minimal x */
  float xmin;

  /* Step size between x points */
  float dx;

  /* Minimal y */
  float ymin;

  /* Step size between y points */
  float dy;

  /* Number of element in the x direction of the data */
  int Nx;

  /* Number of element in the y direction of the data */
  int Ny;

  /* Type of boundary conditions. */
  enum interpolate_boundary_condition boundary_condition;
};


/**
 * @brief Initialize the #interpolation_2d. Stock the data (with "Data limits" proportion) into a flattened 1D array (with "Interpolation limits" proportion).
 *
 * @param interp The #interpolation_2d.
 * @param log_xmin Minimal value of x to (in log).  Interpolation limits
 * @param log_xmax Maximal value of x (in log).   Interpolation limts
 * @param Nx Requested number of values in x axes.  Interpolation limits
 * @param log_ymin Minimal value of y (in log).   Interpolation limits
 * @param log_ymax Maximal value of y (in log).   Interpolation limits
 * @param Ny Requested number of values in y axes.  Interpolation limits
 * @param log_data_xmin The minimal value of the data in x (in log).  Data limits
 * @param log_data_ymin The minimal value of the data in y (in log).  Data limits
 * @param log_step_size_x The size of the x steps (in log).   Data limits
 * @param log_step_size_y The size of the y steps (in log).   Data limits
 * @param N_data_x The number of element in the data x axis. Data limits
 * @param N_data_y The number of element in the data y axis. Data limits
 * @param data The data to interpolate (y).
 * @param boundary_condition The type of #interpolate_boundary_condition.
 */
__attribute__((always_inline)) static INLINE void interpolate_2d_init(
  struct interpolation_2d *interp, float log_xmin, float log_xmax, int Nx, 
  float log_ymin, float log_ymax, int Ny, float log_data_xmin, float log_data_ymin, 
  float log_step_size_x, float log_step_size_y, int N_data_x, int N_data_y, const float *data,
  enum interpolate_boundary_condition boundary_condition) {

  /* Save the variables */
  interp->Nx = Nx;
  interp->xmin = log_xmin;
  interp->dx = (log_xmax - log_xmin) / (Nx - 1.f);
  interp->boundary_condition = boundary_condition;

  interp->Ny = Ny;
  interp->ymin = log_ymin;
  interp->dy = (log_ymax - log_ymin) / (Ny - 1.f);

  /* Allocate the memory */
  interp->data = malloc(sizeof(double) * Nx * Ny);
  if (interp->data == NULL){
    error("Failed to allocate memory for the interpolation");
  }

  /* Interpolate the data */
  for (int i = 0; i < Nx; i++) {
    const float log_x = log_xmin + i * interp->dx;               
    const float x_k = (log_x - log_data_xmin) / log_step_size_x; 

    for (int j = 0; j < Ny; j++) {
      const float log_y = log_ymin + j * interp->dy;
      const float y_k = (log_y - log_data_ymin) / log_step_size_y;

      /* Data indexes */
      const int idx = x_k;
      const float fx = x_k - idx;
      const int idy = y_k;
      const float fy = y_k - idy;

      /* Check i boundaries */
      if (x_k < 0) {
        switch (boundary_condition) {
          case boundary_condition_error:
            error("Cannot extrapolate");
            break;
          case boundary_condition_zero:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_zero_const:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_const:
            if (idy >= 0) 
              interp->data[i * Ny + j] = data[idy];
            else
              interp->data[i * Ny + j] = data[0];
            break;
          default:
            error("Interpolation type not implemented");
        }
        continue;
      } else if (x_k >= N_data_x) {
        switch (boundary_condition) {
          case boundary_condition_error:
            error("Cannot extrapolate");
            break;
          case boundary_condition_zero:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_zero_const:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_const:
            if (i != 0) 
              interp->data[i * Ny + j] = interp->data[(i - 1) * Ny + j];
            else
              interp->data[i * Ny + j] = 0;
            break;
          default:
            error("Interpolation type not implemented");
        }
        continue;
      }

      /* Check boundaries */
      if (y_k < 0) {
        switch (boundary_condition) {
          case boundary_condition_error:
            error("Cannot extrapolate");
            break;
          case boundary_condition_zero:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_zero_const:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_const:
            if (idx >= 0)
              interp->data[i * Ny + j] = data[idx * N_data_y];
            else
              interp->data[i * Ny + j] = data[0];
            break;
          default:
            error("Interpolation type not implemented");
        }
        continue;
      } else if (y_k >= N_data_y) {
        switch (boundary_condition) {
          case boundary_condition_error:
            error("Cannot extrapolate");
            break;
          case boundary_condition_zero:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_zero_const:
            interp->data[i * Ny + j] = 0;
            break;
          case boundary_condition_const:
            if (j != 0)
              interp->data[i * Ny + j] = interp->data[i * Ny + (j - 1)];
            else
              interp->data[i * Ny + j] = 0;
            break;

          default:
            error("Interpolation type not implemented");
        }
        continue;
      }

      /* Interpolate data[i][j] <=> data[i * Ny + j] */
      const float fx1 = data[idx * N_data_y + idy] * (1. - fx) + data[(idx + 1) * N_data_y + idy] * fx;
      const float fx2 = data[idx * N_data_y + idy + 1] * (1. - fx) + data[(idx + 1) * N_data_y + idy + 1] * fx;
      interp->data[i * Ny + j] = fx1 * (1. - fy) + fx2 * fy;
    }
  }

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
  const float i = (log_x - interp->xmin) / interp->dx;
  const int idx = i;
  const float dx = i - idx;

  const int Nx = interp->Nx;
  const int Ny = interp->Ny;

  const float j = (log_y - interp->ymin) / interp->dy;
  const int idy = j;
  const float dy = j - idy;

  /* Extrapolate */
  if (i < 0) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
        return 0;
      case boundary_condition_zero_const:
        return 0;
      case boundary_condition_const:
        if (j >= 0 && j <= Ny - 1)
          return interp->data[idy];
        else if (j < 0)
          return interp->data[0];
        else
          return interp->data[Ny - 1];
      default:
        error("Interpolation type not implemented");
    }
  } else if (i >= Nx - 1) {
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
        return 0;
      case boundary_condition_zero_const:
        return 0;
      case boundary_condition_const:
        if (j >= 0 && j <= Ny - 1)
          return interp->data[(Nx-1) * Ny + idy];
        else if (j < 0)
          return interp->data[(Nx-1) * Ny];
        else
          return interp->data[(Nx-1) * Ny + Ny - 1];
      default:
        error("Interpolation type not implemented");
    }
  } else if (j < 0){
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
        return 0;
      case boundary_condition_zero_const:
        return 0;
      case boundary_condition_const:
        return interp->data[idx * Ny];
      default:
        error("Interpolation type not implemented");
    }
  } else if (j >= Ny - 1){
    switch (interp->boundary_condition) {
      case boundary_condition_error:
        error("Cannot extrapolate");
        break;
      case boundary_condition_zero:
        return 0;
      case boundary_condition_zero_const:
        return 0;
      case boundary_condition_const:
        return interp->data[idx * Ny + (Ny - 1)];
      default:
        error("Interpolation type not implemented");
    }
  }
  
  /* interpolate */
  const float fx1 = interp->data[idx * Ny + idy] * (1. - dx) + interp->data[(idx + 1) * Ny + idy] * dx;
  const float fx2 = interp->data[idx * Ny + idy + 1] * (1. - dx) + interp->data[(idx + 1) * Ny + idy + 1] * dx;
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
}


#endif  // SWIFT_GEAR_INTERPOLATION_H
