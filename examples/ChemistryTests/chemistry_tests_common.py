#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
# This file is part of SWIFT.
# Copyright (c)  2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
################################################################################
"""Helpers shared by the ChemistryTests metal_profile.py analysis scripts."""


def get_fe_metal_mass(data):
    """Get the Fe metal mass (in code units) from swiftsimio data."""
    if hasattr(data.gas.metal_mass_fractions, "fe"):
        m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:  # This case happens when we run the simulation without --feedback
        m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses
    return m_fe
