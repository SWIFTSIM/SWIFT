#!/bin/bash -l

#SBATCH -J nIFTyClusterSWIFT
#SBATCH -N 1
#SBATCH -o nifty_%j.out
#SBATCH -e nifty_%j.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

#SBATCH -t 72:00:00

../../swift --cosmology --hydro --self-gravity -v 1 --pin --threads=56  nifty.yml
