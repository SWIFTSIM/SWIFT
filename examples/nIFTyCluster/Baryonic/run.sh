#!/bin/bash -l

#SBATCH -J nIFTyClusterSWIFT
#SBATCH -N 1
#SBATCH --tasks-per-node=2
#SBATCH -o nifty_%j.out
#SBATCH -e nifty_%j.err
#SBATCH -p <CHANGEME>
#SBATCH -A <CHANGEME>
#SBATCH --exclusive

#SBATCH -t 72:00:00

mpirun -np 2 ../../swift_mpi --cosmology --hydro --self-gravity -v 1 --pin --threads=14  nifty.yml
