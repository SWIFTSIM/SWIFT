#!/bin/bash

### Options for HBT-HERONS

# Path to the executable
hbt_path=/path/to/HBT-HERONS
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
../../../swift --cosmology --self-gravity --fof --threads=8 small_cosmo_volume_hbt.yml 2>&1 | tee output.log

### Run HBT-HERONS
export OMP_NUM_THREADS=${num_omp_threads}
mpirun -np ${num_mpi_ranks} ${hbt_path} HBT.conf

### Plot the halo mass function
python plotHMF.py
