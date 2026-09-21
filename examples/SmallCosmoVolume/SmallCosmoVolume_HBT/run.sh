#!/bin/bash

### Options for SWIFT

# Number of threads to use
num_swift_threads=8

### Options for HBT-HERONS

# Path to the HBT executable, found in the build directory (see README)
hbt_path=/path/to/HBT-HERONS/build/HBT
# Number of MPI ranks to use
num_mpi_ranks=4
# Number of OpenMP threads per MPI rank
num_omp_threads=2

### Fetch the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

### Create output directories
mkdir -p ./outputs/SWIFT ./outputs/HBT-HERONS

### Run SWIFT with FOF
# This box is small enough to run on a single node, so we use the plain (non-MPI)
# SWIFT binary and parallelise using threads alone. Larger runs would instead use
# swift_mpi launched with mpirun. HBT-HERONS is a hybrid MPI+OpenMP code and is
# always launched with mpirun, even when running on a single node.
../../../swift --cosmology --self-gravity --fof --threads=${num_swift_threads} small_cosmo_volume_hbt.yml 2>&1 | tee output.log

### Run HBT-HERONS
export OMP_NUM_THREADS=${num_omp_threads}
mpirun -np ${num_mpi_ranks} ${hbt_path} HBT.conf

### Plot the halo mass function
python plotHMF.py
