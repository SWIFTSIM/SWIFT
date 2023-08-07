#!/bin/bash -l

#SBATCH -n 2
# #SBATCH --nodes=4
#SBATCH --ntasks=32
# #SBATCH --ntasks=64
#SBATCH --cpus-per-task=64
# ---- #SBATCH --cpus-per-task=32
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#SBATCH --mem-per-cpu=4096
#SBATCH --ntasks-per-core=1

# NOTE: In nodes with hyper-threading enabled, a task not requesting full cores 
#  may be distributed across sockets. 
#  This can be avoided by specifying --ntasks-per-core=1,
#  which forces tasks to allocate full cores.

#SBATCH --tasks-per-node 64
# #SBATCH --sockets-per-node=<sockets>
# #SBATCH --cores-per-socket=<cores>
# #SBATCH --threads-per-core=<threads>
# #SBATCH --ntasks-per-core=1




./create-target-file.sh $0

# We can specify options even below the call of create-target-file.sh 
#SBATCH --tasks-per-node 64

# mpirun -np $SLURM_NTASKS your_program your_inputs $SLURM_ARRAY_TASK_ID


# swift_mpi -np $SLURM_NTASKS --with-arch=cosma8 your_program your_inputs $SLURM_ARRAY_TASK_ID
# swift_mpi --with-arch=cosma8 -np $SLURM_NTASKS your_program your_inputs $SLURM_ARRAY_TASK_ID
# swift_mpi --with-arch=cosma8 -np $SLURM_NTASKS your_program your_inputs $SLURM_ARRAY_TASK_ID

    # swift-mpi --with-arch=cosma8 -np $SLURM_NTASKS your_program your_inputs $SLURM_ARRAY_TASK_ID
mpirun /cosma/home/dp004/dc-kots1/swiftsim-ucl-dp004/swift_mpi --with-arch=cosma8 -np $SLURM_NTASKS --threads=16  --cosmology --hydro --self-gravity --stars eagle_6.yml | tee output.log

