#/bin/bash

export HDF5_INC=/usr/include/hdf5/serial
export HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/

rm testcuda.o

nvcc -c -g -D_FORCE_INLINES -O3 -lineinfo -I$HDF5_INC -src-in-ptx --maxrregcount=32 -ftz=true -m64 -ccbin=gcc testcuda.cu -o testcuda.o

echo "############################"
echo "linking"

#gcc -g -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -ffast-math -funroll-loops -march=core-avx2 -mavx2 -fno-tree-vectorize -m64 -Wall -Wextra -Wno-unused-parameter -Werror -o testcuda testcuda.o  -lcudart -lcuda -lnuma -lm -O0 -gdwarf -fvar-tracking-assignments -ldl -lz -lsz
gcc -g -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -ffast-math -funroll-loops -march=core-avx2 -mavx2 -fno-tree-vectorize -m64 -Wall -Wextra -Wno-unused-parameter -Werror -o testcuda testcuda.o -L/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/lib64 -lcudart  -lnuma -lm -O0 -gdwarf -fvar-tracking-assignments -ldl -lz