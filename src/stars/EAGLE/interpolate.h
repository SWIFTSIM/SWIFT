/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_STARS_INTERPOLATE_H
#define SWIFT_EAGLE_STARS_INTERPOLATE_H

/**
 * @brief linear interpolation of 1d table at bin i with offset dx
 *
 * @param table array to interpolate
 * @param i index of cell to interpolate
 * @param dx offset within cell to interpolate
 */
inline static double interpolate_1d(double* table, int i, float dx) {
  double result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/**
 * @brief linear interpolation of 2d table at bin i,j with offset dx, dy
 *
 * @param table array to interpolate
 * @param i row index of cell to interpolate
 * @param j column index of cell to interpolate
 * @param dx row offset within cell to interpolate
 * @param dy column offset within cell to interpolate
 */
inline static double interpolate_2d(double** table, int i, int j, float dx,
                                    float dy) {
  double result;

  result = (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] +
           dx * (1 - dy) * table[i + 1][j] + dx * dy * table[i + 1][j + 1];

  return result;
}

/**
 * @brief linear interpolation of non-uniformly spaced 1d array, array_y, whose
 * positions are specified in array_x. The function takes an input value in the
 * range of array_x and returns a value interpolated from array_y with the same
 * offset in the corresponding bin.
 *
 * @param array_x array of values indicating positions of the array to be
 * interpolated
 * @param array_y array to interpolate
 * @param size length of array_x and array_y
 * @param x value within range of array_x indicating bin and offset within
 * array_y to interpolate
 */
inline static double interpolate_1D_non_uniform(double* array_x,
                                                double* array_y, int size,
                                                double x) {

  double result;

  /* Check that x within range of array_x */
  if (x < array_x[0])
    error("interpolating value less than array min. value %.5e array min %.5e",
          x, array_x[0]);
  else if (x > array_x[size - 1])
    error(
        "interpolating value greater than array max. value %.5e array max %.5e",
        x, array_x[size - 1]);
  else {
    /* Find bin index and offset of x within array_x */
    int index = 0;
    while (array_x[index] <= x) index++;
    double offset =
        (array_x[index] - x) / (array_x[index] - array_x[index - 1]);

    /* Interpolate array_y */
    result = offset * array_y[index - 1] + (1 - offset) * array_y[index];
  }

  return result;
}

#endif /* SWIFT_EAGLE_STARS_INTERPOLATE_H */
