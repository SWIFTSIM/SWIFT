#!/bin/bash -l
# Example with $ntasks MPI tasks and 16 cpus per task
# Targeting $ntasks NUMA regions
# Project/Account (use your own)
#SBATCH -A dp004

#SBATCH -p cosma8

# Number of MPI tasks
#SBATCH --ntasks=32 

#SBATCH --cpus-per-task=16

# Runtime of this jobs is less then 12 hours.
#SBATCH --time=00:40:00

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

module load cosma/2018 python/3.6.5 intel_comp/2022.1.2 compiler openmpi/4.1.1 fftw/3.3.9 parallel_hdf5/1.12.0 parmetis/4.0.3-64bit metis/5.1.0-64bit gsl/2.5

# And finally run the job
mpirun --map-by numa /cosma8/data/dp004/dc-gile1/swiftsim-scotch/swift_mpi --threads=16  --cosmology --hydro --self-gravity --stars eagle_6.yml | tee output.log
# End of submit file