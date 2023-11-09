# This file is part of SWIFT.
# Copyright (c) 2012 pedro.gonnet@durham.ac.uk
#               2015 schaller@strw.leidenuniv.nl
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

# Automake stuff
ACLOCAL_AMFLAGS = -I m4

# Show the way...
if HAVECSDS
SUBDIRS = csds
else
SUBDIRS =
endif
SUBDIRS += src argparse examples doc tests tools
if HAVEEAGLECOOLING
SUBDIRS += examples/Cooling/CoolingRates
DIST_SUBDIRS = $(SUBDIRS)
else
DIST_SUBDIRS = $(SUBDIRS) examples/Cooling/CoolingRates
endif

# Common flags
MYFLAGS =

# Add the source directory and the non-standard paths to the included library headers to CFLAGS
AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/argparse $(HDF5_CPPFLAGS) \
	$(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS) $(GRACKLE_INCS) \
	$(CHEALPIX_CFLAGS)

AM_LDFLAGS = $(HDF5_LDFLAGS)

# Extra libraries.
EXTRA_LIBS = $(GSL_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(PROFILER_LIBS) \
	$(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) \
	$(CHEALPIX_LIBS)

# MPI libraries.
MPI_LIBS = $(PARMETIS_LIBS) $(METIS_LIBS) $(MPI_THREAD_LIBS) $(FFTW_MPI_LIBS)
MPI_FLAGS = -DWITH_MPI $(PARMETIS_INCS) $(METIS_INCS) $(FFTW_MPI_INCS)

# Programs.
bin_PROGRAMS = swift

# Also build the FOF tool?
if HAVESTANDALONEFOF
bin_PROGRAMS += fof
endif

# Do we have the CSDS?
if HAVECSDS
LD_CSDS = csds/src/.libs/libcsds_writer.a
else
LD_CSDS =
endif

# Build MPI versions as well?
if HAVEMPI
bin_PROGRAMS += swift_mpi
if HAVESTANDALONEFOF
bin_PROGRAMS += fof_mpi
endif
endif

# engine_policy_setaffinity is available?
if HAVESETAFFINITY
ENGINE_POLICY_SETAFFINITY=| engine_policy_setaffinity
else
ENGINE_POLICY_SETAFFINITY=
endif

# Sources for swift
swift_SOURCES = swift.c
swift_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_LDADD =  src/libswiftsim.la argparse/libargparse.la $(VELOCIRAPTOR_LIBS) $(EXTRA_LIBS) $(LD_CSDS)

# Sources for swift_mpi, do we need an affinity policy for MPI?
swift_mpi_SOURCES = swift.c
swift_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_mpi_LDADD =  src/libswiftsim_mpi.la argparse/libargparse.la $(MPI_LIBS) $(VELOCIRAPTOR_MPI_LIBS) $(EXTRA_LIBS) $(LD_CSDS)

# Sources for fof
fof_SOURCES = swift_fof.c
fof_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
fof_LDADD =  src/.libs/libswiftsim.a argparse/.libs/libargparse.a $(VELOCIRAPTOR_LIBS) $(EXTRA_LIBS) $(LD_CSDS)

# Sources for fof_mpi, do we need an affinity policy for MPI?
fof_mpi_SOURCES = swift_fof.c
fof_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
fof_mpi_LDADD =  src/.libs/libswiftsim_mpi.a argparse/.libs/libargparse.a $(MPI_LIBS) $(VELOCIRAPTOR_MPI_LIBS) $(EXTRA_LIBS) $(LD_CSDS)

# Non-standard files that should be part of the distribution.
EXTRA_DIST = INSTALL.swift .clang-format format.sh
