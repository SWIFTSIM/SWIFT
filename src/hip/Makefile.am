SOURCES_HIP = HIP_runner_functions.hip
include_HEADERS = HIP_runner_functions.h device_functions.h BLOCK_SIZE.h tester.h
EXTRA_DIST = $(SOURCES_HIP) $(include_HEADERS)

if HAVEHIP

AM_CFLAGS = -I.. $(HDF5_CPPFLAGS)
HIP_MYFLAGS = -D_FORCE_INLINES -O3 -g -DWITH_HIP --offload-arch=gfx90a
#HIP_MYFLAGS = -D_FORCE_INLINES -O3 -g -v -lineinfo -src-in-ptx --maxrregcount=32 -ftz=true -DWITH_HIP -ccbin=gcc -m64 --default-stream per-thread#-dlink

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS) -version-info 0:0:0

#bin_PROGRAMS = test_27_cells test_125_cells

# Rules to compile HIP code.
.hip.o:
	$(HIPCC) -c $(HIPFLAGS) $(AM_CFLAGS) $(HIP_CFLAGS) $(HIP_MYFLAGS) $< -o $@
.hip.lo:
	PATH=$(top_srcdir):$(PATH) && cudalt.py $@ $(HIPCC) -c $(HIPFLAGS) $(AM_CFLAGS) $(HIP_CFLAGS) $(HIP_MYFLAGS) $<

# The library. Dummy C library so that we get libtool linking setup.
lib_LTLIBRARIES = libswiftHIP.la libswiftdummy.la

# Special link command to avoid including CFLAGS which are not understood.
libswiftHIP_la_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CC \
        $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=link $(CCLD) \
        $(libswiftHIP_la_LDFLAGS) \
        $(LDFLAGS) -o $@

libswiftHIP_la_SOURCES = $(SOURCES_HIP)
libswiftHIP_la_CFLAGS = $(AM_CFLAGS) $(HIP_CFLAGS) $(HIP_MYFLAGS) ../libswiftsim_hip.la -I../
libswiftHIP_la_LIBADD = ../.libs/libswiftsim_hip.la
libswiftHIP_la_LDFLAGS = $(AM_LDFLAGS)

if HAVEMPI
libswiftHIP_la_CFLAGS += ../libswiftsim_mpihip.la
libswiftHIP_la_LIBADD += ../.libs/libswiftsim_mpihip.la
endif

libswiftdummy_la_SOURCES = dummy.c
libswiftdummy_la_CFLAGS = $(AM_CFLAGS)
libswiftdummy_la_LDFLAGS = $(AM_LDFLAGS)

#test_27_cells_SOURCES=test27cells.c
#test_27_cells_CFLAGS=$(AM_CFLAGS) -DWITH_HIP $(HIP_CFLAGS) 
#test_27_cells_LDADD= ../.libs/libswiftsim_hip.la ../.libs/libswiftsim_mpihip.la libswiftHIP.la $(MPI_LIBS) $(EXTRA_LIBS) $(HIP_LIBS)
#test_27_cells_LDFLAGS = $(AM_LDFLAGS) $(HIP_CFLAGS)

#test_125_cells_SOURCES=test125cells.c
#test_125_cells_CFLAGS=$(AM_CFLAGS) -DWITH_HIP $(HIP_CFLAGS)
#test_125_cells_LDADD= ../libswiftsim_hip.la ../libswiftsim_mpihip.la libswiftHIP.la $(MPI_LIBS) $(EXTRA_LIBS) $(HIP_LIBS)
#test_125_cells_LDFLAGS = $(AM_LDFLAGS) $(HIP_CFLAGS) 

endif
