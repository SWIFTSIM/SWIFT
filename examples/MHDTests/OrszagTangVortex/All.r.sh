#!/bin/bash

#SBATCH --job-name=Vortex
#SBATCH --partition=short
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=2
#SBATCH --nodes=1

### Cargar el entorno del usuario incluyendo la funcionalidad de modules
### No tocar
#. /etc/profile

### Cargar los módulos para la tarea
# FALTA: Agregar los módulos necesarios

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_N=128

#module purge

#module load gcc
#module load gsl
#module load openmpi
#module load hdf5 parmetis

#ln -sf out.$SLURM_JOB_ID.log last.log
echo "GO GO GO"

IDs=(VeP FDI ODI)
IDs=(VeP)

for J in ${IDs[*]}
do
   echo $J
   cd $J
   export PAR_FILE="../OT.yml"
   #srun ./sw_$J"_mpi" --hydro --threads=$OMP_NUM_THREADS $PAR_FILE 2>&1 | tee $J/out.$SLURM_JOB_ID.log
   cp ../../../../sw_$J .
   ./sw_$J --hydro --threads=$OMP_N $PAR_FILE 2>&1 1> out.log
   cd ..
done


