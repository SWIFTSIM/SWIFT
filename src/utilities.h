/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

#ifndef SWIFT_UTILITIES_H
#define SWIFT_UTILITIES_H

/**
 * @brief Search for a value in a monotonic array to find the index such that
 *      array[index] < value < array[index + 1]
 *
 * @param x The value to find
 * @param arr The array to search
 * @param n The length of the array
 * @param i The found index
 *
 * Set -1 and n for x below and above the array edge values respectively.
 */
INLINE static void find_value_in_monotonic_array(float x, float *arr, int n,
                                                 int *i) {

    int is_incr = (arr[n-1] > arr[0]);  // Increasing or decreasing?
    int i_mid, i_low = 0, i_high = n;

    // Until arr[i_low] < x < arr[i_high=i_low+1]
    while (i_high - i_low > 1) {
        i_mid = (i_high + i_low) >> 1;  // Middle index

        // If mid < x and increasing or x < mid and decreasing
        if ((arr[i_mid] <= x) == is_incr)
            i_low = i_mid;
        else
            i_high = i_mid;
    }

    // Set i with the found i_low, or an error value if outside the array
    if (x < arr[0])
        *i = -1;
    else if (arr[n-1] <= x)
        *i = n;
    else
        *i = i_low;
}

#endif /* SWIFT_UTILITIES_H */

































