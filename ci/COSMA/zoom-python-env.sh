#!/bin/bash

# This file is part of SWIFT.
# Copyright (C) 2026 p.w.draper@durham.ac.uk.
#               2026 w.roper@sussex.ac.uk
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
#  Python environment for zoom CI examples on COSMA.
#  Source this file from a CI job script.
#-

if [ -z "${WORKDIR}" ]; then
	echo "ERROR: WORKDIR is not set. Source ci/setup.sh first."
	return 1
fi

VENV_DIR="${WORKDIR}/.venv-zoom-ci"

echo "Setting up Python environment for zoom CI examples."
python3 -m venv "${VENV_DIR}" || return 1
source "${VENV_DIR}/bin/activate" || return 1

python3 -m pip install --upgrade pip || return 1
python3 -m pip install numpy h5py unyt swiftsimio || return 1

python3 -c "import numpy, h5py, unyt, swiftsimio" || return 1
