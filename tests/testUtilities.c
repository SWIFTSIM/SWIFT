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
    const int n = 100;
    float arr[n];
    int i;
    float x;

    // Initialise test array
    for (int j = 0; j < n; j++) {
        arr[j] = j;
    }

    // Typical example
    x = 42.42f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != 42) {
        printf("Failed with a normal value \n");
        return 1;
    }

    // On array element
    x = 33.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != 33) {
        printf("Failed with an array element \n");
        return 1;
    }

    // Below array
    x = -123.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != -1) {
        printf("Failed with a value below the array \n");
        return 1;
    }

    // Above array
    x = 123.f;
    find_value_in_monotonic_array(x, arr, n, &i);
    if (i != n) {
        printf("Failed with a value above the array \n");
        return 1;
    }

    return 0;
}
