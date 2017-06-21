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
