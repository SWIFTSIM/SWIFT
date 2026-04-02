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

#  Print the runtime of an integration run.
run_with_timer() {
	local label="$1"
	shift

	local start_time=$SECONDS
	"$@"
	local duration=$((SECONDS - start_time))

	echo
	echo "Runtime: ${label} completed in ${duration}s"
}

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
cd "${WORKDIR}/examples/ZoomSimulations/UniformDMGravity"
run_with_timer "UniformDMGravity, 16 steps" do_run bash run.sh

echo
echo "------------------------------------------------------------"
echo "Zoom integration check: UniformDMGravityWithHoles, 16 steps"
echo "------------------------------------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/UniformDMGravityWithHoles"
run_with_timer "UniformDMGravityWithHoles, 16 steps" do_run bash run.sh

echo
echo "------------------------------------------------------"
echo "Zoom integration check: OffsetUniDMGravity, 16 steps"
echo "------------------------------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/OffsetUniDMGravity"
run_with_timer "OffsetUniDMGravity, 16 steps" do_run bash run.sh

echo
echo "------------------------------------------------------------"
echo "Zoom integration check: BoundaryOffsetUniDMGrav, 16 steps"
echo "------------------------------------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/BoundaryOffsetUniDMGrav"
run_with_timer "BoundaryOffsetUniDMGrav, 16 steps" do_run bash run.sh

echo
echo "------------------------------------------------------"
echo "Zoom MPI integration runs, 4 ranks, 3 steps"
echo "------------------------------------------------------"

echo "-------------------------------"
echo "UniformDMGravity MPI integration run"
echo "-------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/UniformDMGravity"
run_with_timer "UniformDMGravity MPI grid/none, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:grid" \
	--param="DomainDecomposition:repartition_type:none" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_gravity.yml

echo "------------------------------------------"
echo "UniformDMGravity wedge MPI integration run"
echo "------------------------------------------"
run_with_timer "UniformDMGravity MPI wedge/none, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:wedge" \
	--param="DomainDecomposition:repartition_type:none" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_gravity.yml

echo "-----------------------------------------------------"
echo "UniformDMGravity memory/fullcosts integration run"
echo "-----------------------------------------------------"
run_with_timer "UniformDMGravity MPI memory/fullcosts, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:memory" \
	--param="DomainDecomposition:repartition_type:fullcosts" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_gravity.yml

echo "--------------------------------------------------------"
echo "UniformDMGravity edgememory/edgecosts integration run"
echo "--------------------------------------------------------"
run_with_timer "UniformDMGravity MPI edgememory/edgecosts, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:edgememory" \
	--param="DomainDecomposition:repartition_type:edgecosts" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_gravity.yml

echo "------------------------------------------------"
echo "UniformDMGravity region/memory integration run"
echo "------------------------------------------------"
run_with_timer "UniformDMGravity MPI region/memory, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:region" \
	--param="DomainDecomposition:repartition_type:memory" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_gravity.yml

echo "--------------------------------"
echo "OffsetUniDMGravity MPI integration run"
echo "--------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/OffsetUniDMGravity"
run_with_timer "OffsetUniDMGravity MPI grid/none, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:grid" \
	--param="DomainDecomposition:repartition_type:none" \
	--param="DomainDecomposition:trigger:0.1" \
	offset_uni_dm_gravity.yml

echo "--------------------------------------"
echo "BoundaryOffsetUniDMGrav MPI integration run"
echo "--------------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/BoundaryOffsetUniDMGrav"
run_with_timer "BoundaryOffsetUniDMGrav MPI grid/none, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:grid" \
	--param="DomainDecomposition:repartition_type:none" \
	--param="DomainDecomposition:trigger:0.1" \
	boundary_offset_uni_dm_grav.yml

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
cd "${WORKDIR}/examples/ZoomSimulations/UniformDMHydro"
run_with_timer "UniformDMHydro, 16 steps" do_run bash run.sh

echo
echo "------------------------------------------------------"
echo "UniformDMHydro MPI integration run, 4 ranks, 3 steps"
echo "------------------------------------------------------"
run_with_timer "UniformDMHydro MPI grid/none, 4 ranks, 3 steps" do_run mpirun -np 4 ../../../swift_mpi --cosmology --hydro --self-gravity --zoom \
	--threads=1 -n 3 --no-io \
	--param="DomainDecomposition:initial_type:grid" \
	--param="DomainDecomposition:repartition_type:none" \
	--param="DomainDecomposition:trigger:0.1" \
	zoom_uniform_dm_hydro.yml

cd "${WORKDIR}"
do_make clean

echo
echo "-----------------------------------------------------"
echo "Configure group: full physics zoom integration builds"
echo "-----------------------------------------------------"

do_configure --with-parmetis --with-subgrid=EAGLE-XL --with-hydro=sphenix --with-kernel=wendland-C2 --enable-debugging-checks --disable-vec --enable-debug
do_make

echo
echo "-------------------------------------------------"
echo "Zoom integration check: UniformDMEAGLE, 16 steps"
echo "-------------------------------------------------"
cd "${WORKDIR}/examples/ZoomSimulations/UniformDMEAGLE"

# Link tables from shared HOME location when available; run.sh downloads
# any missing dependencies.
if [ ! -e yieldtables ] && [ -e "$HOME/yieldtables" ]; then
	link_data yieldtables
fi
if [ ! -e coolingtables ] && [ -e "$HOME/coolingtables" ]; then
	link_data coolingtables
fi
if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ] && [ -e "$HOME/UV_dust1_CR1_G1_shield1.hdf5" ]; then
	link_data UV_dust1_CR1_G1_shield1.hdf5
fi
if [ ! -e photometry ] && [ -e "$HOME/photometry" ]; then
	link_data photometry
fi

run_with_timer "UniformDMEAGLE, 16 steps" do_run bash run.sh

exit
