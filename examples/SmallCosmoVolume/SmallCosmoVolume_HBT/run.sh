#!/bin/bash

set -e

# TODO: Remove
module purge
module load old-modules
module load intel_comp/2024.2.0 compiler-rt tbb compiler
module load openmpi/5.0.3
module load ucx/1.13.0rc2
module load parallel_hdf5/1.14.4
module load fftw/3.3.10
module load parmetis
module load gsl/2.5



### Options for HBT-HERONS
# Path to the executable
hbt_path=/path/to/HBT-HERONS
#TODO:
hbt_path=/cosma/home/dp004/dc-mcgi1/swift/hbt_example/HBT-HERONS/build/HBT
# Number of MPI ranks to use
num_mpi_ranks=4
# Number of OpenMP threads per MPI rank
num_omp_threads=2

# Fetch the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

# Create output directories
mkdir -p ./outputs/SWIFT ./outputs/HBT-HERONS

# Run SWIFT with FOF
# TODO: Set threads
../../../swift --cosmology --self-gravity --fof --threads=16 small_cosmo_volume_hbt.yml 2>&1 | tee output.log

# Run HBT-HERONS
export OMP_NUM_THREADS=${num_omp_threads}
mpirun -np ${num_mpi_ranks} ${hbt_path} HBT.conf

# Plot the halo mass function
python plotHMF.py
