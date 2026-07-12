#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

from pathlib import Path

import numpy as np
import swiftsimio

from makeIC import ELEMENT_COUNT


def test(name):
    def test_decorator(test_func):
        def test_runner(*args, **kwargs):
            try:
                test_func(*args, **kwargs)
            except Exception as e:
                print(f"\tTest {name} failed!")
                raise e
            else:
                print(f"\tTest {name} \u2705")

        return test_runner

    return test_decorator


def approx_equal(a, b, threshold=0.001):
    return 2 * abs(a - b) / (abs(a) + abs(b)) < threshold


@test("total metal mass conservation")
def test_total_metal_mass_conservation(data_start, data_end):
    def metal_mass(data):
        columns = getattr(data.gas.metal_mass_fractions, "named_columns", None)
        if columns is not None:
            return sum(
                data.gas.masses * getattr(data.gas.metal_mass_fractions, c)
                for c in columns
            )
        else:
            return data.gas.masses.reshape(-1, 1) * data.gas.metal_mass_fractions

    assert approx_equal(np.sum(metal_mass(data_start)), np.sum(metal_mass(data_end)))


def element_mass(data, element_idx):
    columns = getattr(data.gas.metal_mass_fractions, "named_columns", None)
    if columns is not None:
        return data.gas.masses * getattr(
            data.gas.metal_mass_fractions, columns[element_idx]
        )
    else:
        return data.gas.masses * data.gas.metal_mass_fractions[:, element_idx]


@test("element-wise mass conservation")
def test_element_wise_mass_conservation(data_start, data_end):
    for i in range(ELEMENT_COUNT):
        assert approx_equal(
            np.sum(element_mass(data_start, i)), np.sum(element_mass(data_end, i))
        )


if __name__ == "__main__":
    print("Running sanity checks...")

    cwd = Path(__file__).parent
    start = swiftsimio.load(cwd / "output_0000.hdf5")
    end = swiftsimio.load(cwd / "output_0001.hdf5")

    test_total_metal_mass_conservation(start, end)
    test_element_wise_mass_conservation(start, end)

    print("Done!")
