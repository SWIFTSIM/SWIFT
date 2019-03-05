# tHIS FIle is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Add the source directory and the non-standard paths to the included library headers to CFLAGS
AM_CFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/examples $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS)

AM_LDFLAGS = $(HDF5_LDFLAGS) $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) $(GSL_LIBS) $(PROFILER_LIBS)

# Extra libraries.
EXTRA_LIBS = $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS) $(GSL_LIBS)

# Programs.
bin_PROGRAMS = cooling_rates

# Sources
cooling_rates_SOURCES = cooling_rates.c
cooling_rates_CFLAGS = $(AM_CFLAGS)
cooling_rates_LDADD =  ../../../src/.libs/libswiftsim.a $(EXTRA_LIBS)

