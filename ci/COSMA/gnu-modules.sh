#!/bin/bash

# This file is part of SWIFT.
# Copyright (C) 2026 p.w.draper@durham.ac.uk.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#+
#  GNU/GCC toolchain modules with serial HDF5.
#  Source this file.
#
#  Peter W. Draper 20-JAN-2026.
#-
echo "Loading the GNU/GCC toolchain for SWIFT."

module purge
module load gnu_comp/14.1.0
module load openmpi/5.0.3
module load gsl/2.8
module load hdf5/1.14.4
module load parmetis/4.0.3
module load fftw/3.3.10
module load utils
module load llvm/20.1.0
export CLANG_FORMAT_CMD="clang-format"
