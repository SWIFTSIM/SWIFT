/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

#include "swift.h"
#include "utilities.h"

/**
 * @brief Test generic utility functions
 */
int main() {
    /// Test find_value_in_monotonic_array()
    int n = 100;
    float arr[n];
    int i;
    float x;

    // Initialise test array
    for (int j = 0; j < n; j++) {
        arr[j] = j;
    }

    // Typical value
    x = 42.42f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != 42) {
        printf("Failed with a typical value \n");
        return 1;
    }

    // Value on array element
    x = 33.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != 33) {
        printf("Failed with an array element \n");
        return 1;
    }

    // Value below array
    x = -123.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != -1) {
        printf("Failed with a value below the array \n");
        return 1;
    }

    // Value above array
    x = 123.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != n) {
        printf("Failed with a value above the array \n");
        return 1;
    }

    // Array slice with typical value
    x = 9.81f;
    n = 10;
    find_value_in_monotonic_array(x, arr + 5, n, &i);
    if (i != 4) {
        printf("Failed with an array slice \n");
        return 1;
    }

    return 0;
}
