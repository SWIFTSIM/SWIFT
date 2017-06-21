
# Add the debug flag to the whole thing
AM_CFLAGS = -I.. $(HDF5_CPPFLAGS)

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS) -version-info 0:0:0

# CUDA sources
if HAVECUDA
    SOURCES_CUDA = runner_cuda_main.cu queue_cuda.cu
    CUDA_MYFLAGS = -O3 -g -lineinfo -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_CUDA
#.cu: cell_cuda.h queue_cuda.h runner_cuda_main.h task_cuda.h
.cu.o: 
	$(NVCC) -c $(NVCCFLAGS) $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $< -o $@
.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $<
else
    SOURCES_CUDA =
endif

if HAVECUDA
lib_LTLIBRARIES = libswiftCUDA.la
endif

include_HEADERS = cell_cuda.h queue_cuda.h runner_cuda_main.h task_cuda.h

libswiftCUDA_la_SOURCES = $(SOURCES_CUDA)
libswiftCUDA_la_CFLAGS = $(AM_CFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS)
libswiftCUDA_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_MYFLAGS)
