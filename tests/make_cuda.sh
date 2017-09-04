#/bin/bash

CUDA_CFLAGS="-I/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/include -I/opt/cray/pe/hdf5/1.10.0/GNU/5.1/include"
LIBS="-lnuma  -lm"

nvcc -c -I../src -g -D_FORCE_INLINES -O3 -lineinfo $CUDA_CFLAGS $LIBS -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_CUDA  -ccbin=gcc-4.8 -m64 testcuda.cu -o testcuda.o

echo "############################"
echo "linking"

cc -I.. -g -DWITH_CUDA -I/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/include -I/opt/cray/pe/hdf5/1.10.0/GNU/5.1/include -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -ffast-math -funroll-loops -march=haswell -mavx2 -fno-tree-vectorize -m64 -Wall -Wextra -Wno-unused-parameter -Werror -I/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/include -I/opt/cray/pe/hdf5/1.10.0/GNU/5.1/include-L/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/lib -L/opt/nvidia/cudatoolkit8.0/8.0.54_2.2.8_ga620558-2.1/lib64 -lhdf5 -lcudart -lnuma -lm  -o testcuda testcuda.o  ../src/CUDA/.libs/libswiftCUDA.a
