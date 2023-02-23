.. Running on Large Systems
   Josh Borrow, 5th April 2018

Running on Large Systems
========================

There are a few extra things to keep in mind when running SWIFT on a large
system (i.e. over MPI on several nodes). Here are some recommendations:

+ Compile and run with
  `tbbmalloc <https://www.threadingbuildingblocks.org>`_.  You can add this
  to the configuration of SWIFT by running configure with the
  ``--with-tbbmalloc`` flag. Using this allocator, over the one included in the
  standard library, is particularly important on systems with large core counts
  per node. Alternatives include
  `jemalloc <https://github.com/jemalloc/jemalloc>`_ and
  `tcmalloc <https://github.com/gperftools/gperftools>`_, and using these
  other allocation tools also improves performance on single-node jobs.
+ Run with one MPI rank per NUMA region, usually a socket, rather than per node.
  Typical HPC clusters now use two chips per node. Consult with your local system
  manager if you are unsure about your system configuration. This can be done
  by invoking ``mpirun -np <NUMBER OF CHIPS> swift_mpi -t <NUMBER OF CORES PER CHIP>``.
  You should also be careful to include this in your batch script, for example
  with the `SLURM <https://slurm.schedmd.com>`_ batch system you will need to
  include ``#SBATCH --tasks-per-node=2``.
+ Run with threads pinned. You can do this by passing the ``-a`` flag to the
  SWIFT binary. This ensures that processes stay on the same core that spawned
  them, ensuring that cache is accessed more efficiently.
+ Ensure that you compile with ParMETIS or METIS. These are required if
  want to load balance between MPI ranks.

Your batch script should look something like the following (to run on 16 nodes
each with 2x16 core processors for a total of 512 cores):

.. code-block:: bash
  
   #SBATCH -N 16  # Number of nodes to run on
   #SBATCH --tasks-per-node=2  # This system has 2 chips per node
   
   mpirun -np 32 swift_mpi --threads=16 --pin parameter.yml

