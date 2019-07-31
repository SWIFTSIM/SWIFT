#!/bin/bash -l

#SBATCH -J SWIFTDiffusionCalibration
#SBATCH -N 1
#SBATCH -o swift_diffusion.out
#SBATCH -e swift_diffusion.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

#SBATCH -t 1:00:00


for diffusion_alpha_max in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0;
do
    mkdir default_diffmax_$diffusion_alpha_max

    cd default_diffmax_$diffusion_alpha_max

    ../../../swift --hydro --cooling --limiter --threads=16 --param="SPH:diffusion_alpha_max:${diffusion_alpha_max}" ../feedback.yml 2>&1 | tee output.log 

    cd ..
    
    mkdir nocool_diffmax_$diffusion_alpha_max

    cd nocool_diffmax_$diffusion_alpha_max

    ../../../swift --hydro --temperature --limiter --threads=16 --param="SPH:diffusion_alpha_max:${diffusion_alpha_max}" ../feedback.yml 2>&1 | tee output.log 

    cd ..
done
