#!/bin/bash -l

# ==============
# Batch script for GIHR simulations with HOT
# ==============
#BSUB -L /bin/bash

# ============== Number of processors (and one MPI rank per node)
#BSUB -n 1
#BSUB -R "span[ptile=1]"

# ============== Job name
#BSUB -J SWIFT

# ============== Output files
#BSUB -oo /cosma5/data/dp004/dc-kege1/gihr_data5/swift/outs_and_errs/swift_%J.out
#BSUB -eo /cosma5/data/dp004/dc-kege1/gihr_data5/swift/outs_and_errs/swift_%J.err

# ============== Queue
# #BSUB -q bench1
#BSUB -q cosma
# #BSUB -q cosma5
# #BSUB -q cosma5-prince

# ============== Project and user
# #BSUB -P dp004
#BSUB -P durham
#BSUB -u jacob.kegerreis@durham.ac.uk
#BSUB -N

# ============== Wall clock time limit
# #BSUB -W 504:00
#BSUB -W 12:00
# #BSUB -W 0:30

# ============== Ensure that the right MPI module is loaded -- i.e. the same
#                module with which the program was compiled.
module purge
module load swift
module load swift/c4/gcc/intelmpi/5.1.3

# ============== Run the program
../swift -G -s -t 12 ./uranus_m2L2v5_1e5.yml
# mpirun -np $LSB_DJOB_NUMPROC ../swift_mpi -G -s -t 12 ./uranus_m2L2v5_1e5.yml




