SOURCES_CUDA = runner_cuda_main.cu
include_HEADERS = cell_cuda.h queue_cuda.h runner_cuda_main.h task_cuda.h
EXTRA_DIST = $(SOURCES_CUDA) $(include_HEADERS)

if HAVECUDA

AM_CFLAGS = -I.. $(HDF5_CPPFLAGS) -g
CUDA_MYFLAGS = -D_FORCE_INLINES -O0 -g -lineinfo -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_CUDA -G -ccbin=gcc-4.8
#-dc

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS) -version-info 0:0:0

bin_PROGRAMS = test_27_cells test_125_cells

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
libswiftCUDA_la_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS)
libswiftCUDA_la_LDFLAGS = $(AM_LDFLAGS)

libswiftdummy_la_SOURCES = dummy.c
libswiftdummy_la_CFLAGS = $(AM_CFLAGS)
libswiftdummy_la_LDFLAGS = $(AM_LDFLAGS)

test_27_cells_SOURCES=test27cells.c
test_27_cells_CFLAGS=$(AM_CFLAGS) -DWITH_CUDA $(CUDA_CFLAGS) 
test_27_cells_LDADD= libswiftCUDA.la ../libswiftsim_cuda.la ../libswiftsim_mpicuda.la $(MPI_LIBS) $(EXTRA_LIBS) $(CUDA_LIBS) -lcudart
test_27_cells_LDFLAGS = $(AM_LDFLAGS) $(CUDA_CFLAGS)

test_125_cells_SOURCES=test125cells.c
test_125_cells_CFLAGS=$(AM_CFLAGS) -DWITH_CUDA $(CUDA_CFLAGS)
test_125_cells_LDADD= libswiftCUDA.la ../libswiftsim_cuda.la ../libswiftsim_mpicuda.la $(MPI_LIBS) $(EXTRA_LIBS) $(CUDA_LIBS) -L/home/aidan/cuda_7.5/lib64/ -lcudart
test_125_cells_LDFLAGS = $(AM_LDFLAGS) $(CUDA_CFLAGS) 


endif
