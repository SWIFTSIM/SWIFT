.. What about MPI? Running SWIFT on more than one node
   Josh Borrow, 5th April 2018

What about MPI? Running SWIFT on more than one node
===================================================

After compilation, you will be left with two binaries. One is called ``swift``,
and the other ``swift_mpi``. Current wisdom is to run ``swift`` if you are only
using one node (i.e. without any interconnect), and one MPI rank per NUMA
region using ``swift_mpi`` for anything larger. You will need some initial 
conditions in the GADGET-2 HDF5 format (see :ref:`Initial Conditions <Initial_Conditions_label>`) 
to run SWIFT, as well as a compatible :ref:`yaml parameter file <Parameter_File_label>`.



SLURM Submission Script Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Below is an example submission script for the SLURM
batch system. This runs SWIFT with thread pinning, SPH,
and self-gravity.

.. code-block:: bash

    #SBATCH --partition=<queue>
    #SBATCH --account-name=<groupName>
    #SBATCH --job-name=<jobName>
    #SBATCH --nodes=<nNodes>
    #SBATCH --ntasks-per-node=<nMPIRank>
    #SBATCH --cpus-per-task=<nThreadsPerMPIRank>
    #SBATCH --output=outFile.out
    #SBATCH --error=errFile.err

    ## expected runtime
    #SBATCH --time-<hh>:<mm>:<ss>

    ./swift_mpi --pin --threads=<nThreadsPerMPIRank> \
        --hydro --self-gravity \
        parameter_file.yml

