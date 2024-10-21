SOURCES_CUDA = GPU_runner_functions.cu tester.cu ../files_for_new_functions/arrays_malloc.cu ../files_for_new_functions/host_device_data_transfer.cu #../runner_main.cu
include_HEADERS = GPU_runner_functions.h device_functions.h BLOCK_SIZE.h tester.h ../files_for_new_functions/arrays_malloc.h ../files_for_new_functions/host_device_data_transfer.h
EXTRA_DIST = $(SOURCES_CUDA) $(include_HEADERS)

if HAVECUDA

AM_CFLAGS = -I.. $(HDF5_CPPFLAGS)
CUDA_MYFLAGS = -D_FORCE_INLINES -O4 -lineinfo -src-in-ptx --maxrregcount=64 -ftz=true -DWITH_CUDA --default-stream per-thread --use_fast_math -lcudadevrt #-dlink -ccbin=gcc 
CUDA_MYFLAGS += -arch=sm_70
CUDA_MYFLAGS += --extra-device-vectorization

#CUDA_MYFLAGS = -D_FORCE_INLINES -O3 -lineinfo -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_CUDA -ccbin=gcc -m64 --default-stream per-thread #-dlink
#CUDA_MYFLAGS += -arch=sm_80 \
#-gencode=arch=compute_80,code=sm_80 \
#-gencode=arch=compute_86,code=sm_86 \
#-gencode=arch=compute_87,code=sm_87 \
#-gencode=arch=compute_86,code=compute_86
#CUDA_MYFLAGS += --extra-device-vectorization

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS) -version-info 0:0:0

#bin_PROGRAMS = test_27_cells test_125_cells

# Rules to compile CUDA code.
.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $< -o $@
.cu.lo:
	PATH=$(top_srcdir):$(PATH) && cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $<

# The library. Dummy C library so that we get libtool linking setup.
lib_LTLIBRARIES = libswiftCUDA.la libswiftdummy.la

# Special link command to avoid including CFLAGS which are not understood.
libswiftCUDA_la_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CC \
        $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=link $(CCLD) \
        $(libswiftCUDA_la_LDFLAGS) \
        $(LDFLAGS) -o $@

libswiftCUDA_la_SOURCES = $(SOURCES_CUDA)
libswiftCUDA_la_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) ../libswiftsim_cuda.la -I../
libswiftCUDA_la_CXXFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) ../libswiftsim_cuda.la -I../
libswiftCUDA_la_LIBADD = ../.libs/libswiftsim_cuda.la
libswiftCUDA_la_LDFLAGS = $(AM_LDFLAGS)

if HAVEMPI
libswiftCUDA_la_CFLAGS += ../libswiftsim_mpicuda.la
libswiftCUDA_la_CXXFLAGS += ../libswiftsim_mpicuda.la
libswiftCUDA_la_LIBADD += ../.libs/libswiftsim_mpicuda.la
endif

libswiftdummy_la_SOURCES = dummy.c
libswiftdummy_la_CFLAGS = $(AM_CFLAGS)
libswiftdummy_la_LDFLAGS = $(AM_LDFLAGS)

#test_27_cells_SOURCES=test27cells.c
#test_27_cells_CFLAGS=$(AM_CFLAGS) -DWITH_CUDA $(CUDA_CFLAGS) 
#test_27_cells_LDADD= ../.libs/libswiftsim_cuda.la ../.libs/libswiftsim_mpicuda.la libswiftCUDA.la $(MPI_LIBS) $(EXTRA_LIBS) $(CUDA_LIBS) -L/home/aidan/cuda_7.5/lib64/ -lcudart
#test_27_cells_LDFLAGS = $(AM_LDFLAGS) $(CUDA_CFLAGS)

#test_125_cells_SOURCES=test125cells.c
#test_125_cells_CFLAGS=$(AM_CFLAGS) -DWITH_CUDA $(CUDA_CFLAGS)
#test_125_cells_LDADD= ../libswiftsim_cuda.la ../libswiftsim_mpicuda.la libswiftCUDA.la $(MPI_LIBS) $(EXTRA_LIBS) $(CUDA_LIBS) -L/home/aidan/cuda_7.5/lib64/ -lcudart
#test_125_cells_LDFLAGS = $(AM_LDFLAGS) $(CUDA_CFLAGS) 

endif
