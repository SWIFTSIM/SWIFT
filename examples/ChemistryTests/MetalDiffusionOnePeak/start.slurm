#!/bin/bash
#SBATCH --job-name=L6_chemistry_test
#SBATCH --ntasks=1
#SBATCH --nodes=1 #-N
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=darwin.roduit@epfl.ch
#SBATCH --mail-type=ALL

#srun python3 makeIC.py --level 6 --mass 1e7 -o ICs_homogeneous_box.hdf5 --rho 1 --epsilon 0.01 --metal_mass_fraction 1e-4
srun ~/program/swiftsim/swift --hydro --external-gravity --stars --feedback --threads=${SLURM_CPUS_PER_TASK} params.yml
