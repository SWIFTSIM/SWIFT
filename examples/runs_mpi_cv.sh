#!/bin/bash

# Set the number of runs
FINALTIME=1.

# Cores per node
CPN=12

# The queue on which to run
QUEUE=cosma
PROJECT=durham
PREFIX=CosmoVolume
INPUT=$PREFIX/CosmoVolume.hdf5

# Make sure the OMP threads don't go wild
export OMP_WAIT_POLICY=PASSIVE

# Set the library path so that libmetis is found
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cosma/home/nnrw56/lib

# Single-node runs
for cpu in $(seq 1 $CPN)
do
    if [ ! -e ${PREFIX}_${QUEUE}_1x${cpu}.dump ]
    then
        bsub -oo ${PREFIX}_${QUEUE}_1x${cpu}.dump -q ${QUEUE} -P ${PROJECT} -x -n 1 -R "span[ptile=1]" ./swift -c $FINALTIME -t $cpu -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
    fi
done

# Multi-node runs
if [ ! -e ${PREFIX}_${QUEUE}_2x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_2x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 2 -R "span[ptile=1]" mpirun -np 2 ./swift_mpi -c $FINALTIME -t $CPN -g "2 1 1" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_4x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_4x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 4 -R "span[ptile=1]" mpirun -np 4 ./swift_mpi -c $FINALTIME -t $CPN -g "2 2 1" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_8x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_8x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 8 -R "span[ptile=1]" mpirun -np 8 ./swift_mpi -c $FINALTIME -t $CPN -g "2 2 2" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_16x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_16x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 16 -R "span[ptile=1]" mpirun -np 16 ./swift_mpi -c $FINALTIME -t $CPN -g "4 2 2" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_32x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_32x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 32 -R "span[ptile=1]" mpirun -np 32 ./swift_mpi -c $FINALTIME -t $CPN -g "4 4 2" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_64x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_64x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 64 -R "span[ptile=1]" mpirun -np 64 ./swift_mpi -c $FINALTIME -t $CPN -g "4 4 4" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

if [ ! -e ${PREFIX}_${QUEUE}_128x${cpu}.dump ]
then
    bsub -oo ${PREFIX}_${QUEUE}_128x${CPN}.dump -q ${QUEUE} -P ${PROJECT} -x -W 02:00 -n 128 -R "span[ptile=1]" mpirun -np 128 ./swift_mpi -c $FINALTIME -t $CPN -g "8 4 4" -f ${INPUT} -m 0.705 -w 6000 -z 300 -d 1e-7 -e 0.01
fi

