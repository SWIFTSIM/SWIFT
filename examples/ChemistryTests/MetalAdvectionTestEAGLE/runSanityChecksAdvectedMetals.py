#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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


@test("mass fractions sum to 1")
def test_sum_mass_fractions(data):
    metal_fractions = data.gas.metal_mass_fractions
    h_fraction = data.gas.element_mass_fractions.hydrogen
    he_fraction = data.gas.element_mass_fractions.helium
    total = metal_fractions + he_fraction + h_fraction

    assert np.all(approx_equal(total, 1))


@test("total metal mass conservation")
def test_total_metal_mass_conservation(data_start, data_end):
    def metal_mass(data):
        return data.gas.masses * data.gas.metal_mass_fractions

    assert approx_equal(np.sum(metal_mass(data_start)), np.sum(metal_mass(data_end)))


def element_mass(data, element_name):
    return data.gas.masses * getattr(data.gas.element_mass_fractions, element_name)


@test("hydrogen mass conservation")
def test_h_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "hydrogen")),
        np.sum(element_mass(data_end, "hydrogen")),
    )


@test("helium mass conservation")
def test_he_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "helium")),
        np.sum(element_mass(data_end, "helium")),
    )


@test("carbon mass conservation")
def test_c_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "carbon")),
        np.sum(element_mass(data_end, "carbon")),
    )


@test("nitrogen mass conservation")
def test_n_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "nitrogen")),
        np.sum(element_mass(data_end, "nitrogen")),
    )


@test("oxygen mass conservation")
def test_o_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "oxygen")),
        np.sum(element_mass(data_end, "oxygen")),
    )


@test("neon mass conservation")
def test_ne_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "neon")), np.sum(element_mass(data_end, "neon"))
    )


@test("magnesium mass conservation")
def test_mg_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "magnesium")),
        np.sum(element_mass(data_end, "magnesium")),
    )


@test("silicon mass conservation")
def test_si_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "silicon")),
        np.sum(element_mass(data_end, "silicon")),
    )


@test("iron mass conservation")
def test_fe_mass_conservation(data_start, data_end):
    assert approx_equal(
        np.sum(element_mass(data_start, "iron")), np.sum(element_mass(data_end, "iron"))
    )


if __name__ == "__main__":
    print("Running sanity checks...")

    cwd = Path(__file__).parent
    start = swiftsimio.load(cwd / "output_0000.hdf5")
    end = swiftsimio.load(cwd / "output_0001.hdf5")

    test_sum_mass_fractions(end)
    test_total_metal_mass_conservation(start, end)
    test_h_mass_conservation(start, end)
    test_he_mass_conservation(start, end)
    test_c_mass_conservation(start, end)
    test_n_mass_conservation(start, end)
    test_o_mass_conservation(start, end)
    test_ne_mass_conservation(start, end)
    test_mg_mass_conservation(start, end)
    test_si_mass_conservation(start, end)
    test_fe_mass_conservation(start, end)

    print("Done!")
