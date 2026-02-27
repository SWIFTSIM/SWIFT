#!/bin/bash -l

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
#  Runtime integration checks for zoom functionality on COSMA.
#
#  This script builds SWIFT with debugging checks and runs zoom integration
#  workflows from examples/ZoomSimulations.
#
#  Peter W. Draper 25-FEB-2026.
#-

#  Build toolchain.
source ci/COSMA/intel-modules.sh
source ci/setup.sh

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

#  Python environment for zoom IC generation scripts.
source ci/COSMA/zoom-python-env.sh

echo
echo "------------------------------------------------------"
echo "Configure group: DMO + hydro zoom integration examples"
echo "------------------------------------------------------"
do_configure --with-parmetis --enable-debugging-checks --disable-vec --enable-debug
do_make

echo
echo "-------------------"
echo "DMO integration runs"
echo "-------------------"

echo
echo "---------------------------------------------------"
echo "Zoom integration check: UniformDMGravity, 16 steps"
echo "---------------------------------------------------"
cd examples/ZoomSimulations/UniformDMGravity
do_run bash run.sh

echo
echo "------------------------------------------------------------"
echo "Zoom integration check: UniformDMGravityWithHoles, 16 steps"
echo "------------------------------------------------------------"
cd ../UniformDMGravityWithHoles
do_run bash run.sh

echo
echo "------------------------------------------------------"
echo "Zoom integration check: OffsetUniDMGravity, 16 steps"
echo "------------------------------------------------------"
cd ../OffsetUniDMGravity
do_run bash run.sh

echo
echo "------------------------------------------------------------"
echo "Zoom integration check: BoundaryOffsetUniDMGrav, 16 steps"
echo "------------------------------------------------------------"
cd ../BoundaryOffsetUniDMGrav
do_run bash run.sh

echo
echo "---------------------"
echo "Hydro integration runs"
echo "---------------------"

cd "${WORKDIR}"
do_make clean

echo
echo "-------------------------------------------"
echo "Configure group: hydro zoom integration runs"
echo "-------------------------------------------"
do_configure --with-parmetis --with-hydro=sphenix --with-kernel=wendland-C2 --enable-debugging-checks --disable-vec --enable-debug
do_make

echo
echo "-------------------------------------------------"
echo "Zoom integration check: UniformDMHydro, 16 steps"
echo "-------------------------------------------------"
cd ../UniformDMHydro
do_run bash run.sh

cd "${WORKDIR}"
do_make clean

echo
echo "-----------------------------------------------------"
echo "Configure group: full physics zoom integration builds"
echo "-----------------------------------------------------"

do_configure --with-parmetis --with-subgrid=EAGLE-XL --with-hydro=sphenix --with-kernel=wendland-C2 --enable-debugging-checks --disable-vec --enable-debug
do_make

echo "------------------------------------------------"
echo "Full physics runtime checks currently not enabled"
echo "------------------------------------------------"
echo "UniformDMEAGLE run remains disabled until tuned."

exit
