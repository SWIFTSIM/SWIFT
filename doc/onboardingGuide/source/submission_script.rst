.. Submission Script

Submission Script
=================


Below is an example submission script for the SLURM batch system. This runs 
SWIFT with MPI, thread pinning, hydrodynamics, and self-gravity.

.. code-block:: bash

    #SBATCH --partition=<queue>
    #SBATCH --account-name=<groupName>
    #SBATCH --job-name=<jobName>
    #SBATCH --nodes=<nNodes>
    #SBATCH --ntasks-per-node=<nMPIRank>
    #SBATCH --cpus-per-task=<nThreadsPerMPIRank>
    #SBATCH --time=<hh>:<mm>:<ss>

    srun -n $SLURM_NTASKS ./swift_mpi \ 
        --threads=$SLURM_CPUS_PER_TASK --pin \
        --hydro --self-gravity parameter_file.yml

