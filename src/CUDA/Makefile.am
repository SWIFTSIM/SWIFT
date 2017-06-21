# CUDA sources
if HAVECUDA
    SOURCES_CUDA = runner_cuda_main.cu queue_cuda.cu
    CUDA_MYFLAGS = -O3 -g -lineinfo -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_CUDA
.cu: cell_cuda.h queue_cuda.h runner_cuda_main.h task_cuda.h
.cu.o: $(NVCC) -c $(NVCCFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $< -o $@
.cu.lo: $(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $(CUDA_CFLAGS) $(CUDA_MYFLAGS) $<
else
    SOURCES_CUDA =
endif

if HAVECUDA
lib_LTLIBRARIES = libswiftCUDA.la
endif

include_HEADERS = ../space.h ../runner.h ../queue.h ../task.h ../cell.h ../part.h ../const.h \
  ../engine.h ../swift.h ../serial_io.h ../timers.h ../map.h ../tools.h ../physical_constraints.h \
  ../physical_constants_cgs.h ../potential.h ../version.h ../hydro_properties.h ../dump.h \
  ../hydro_space.h ../sort_part.h

libswiftCUDA_la_SOURCES = $(SOURCES_CUDA)
libswiftCUDA_la_CFLAGS = $(CUDA_CFLAGS) $(CUDA_MYFLAGS)
libswiftCUDA_la_LDFLAGS = $(CUDA_MYFLAGS)
