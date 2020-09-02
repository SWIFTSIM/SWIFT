# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
#                    Loic Hausammann (loic.hausammann@epfl.ch)
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

# Add the non-standard paths to the included library headers
AM_CFLAGS = $(PYTHON_INCS) -I$(top_srcdir)/src $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS) $(GRACKLE_INCS)


AM_LDFLAGS = $(HDF5_LDFLAGS)

# Assign a "safe" version number
BIN_LDFLAGS = -version-info 0:0:0

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Additional dependencies for shared libraries.
EXTRA_LIBS = $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) $(GRACKLE_LIBS) \
	$(VELOCIRAPTOR_LIBS) $(GSL_LIBS) -L../src/.libs -lswiftsim

# MPI libraries.
# MPI_LIBS = $(MPI_THREAD_LIBS)
# MPI_FLAGS = -DWITH_MPI

# Build the liblogger library
lib_LTLIBRARIES = liblogger.la
# Build a MPI-enabled version too?
# if HAVEMPI
# lib_LTLIBRARIES += liblogger_mpi.la
# endif

GRAVITY_SRC = gravity/MultiSoftening/logger_gravity.c
STARS_SRC = stars/Default/logger_stars.c
HYDRO_SRC = hydro/Gadget2/logger_hydro.c

# subdirectories
SUBDIRS = tests

# List required headers
include_HEADERS = logger_header.h logger_loader_io.h logger_particle.h logger_time.h logger_tools.h logger_reader.h \
	logger_logfile.h logger_index.h quick_sort.h logger_python_tools.h

# Common source files
AM_SOURCES = logger_header.c logger_loader_io.c logger_particle.c logger_time.c logger_tools.c logger_reader.c \
	logger_logfile.c logger_index.c quick_sort.c $(GRAVITY_SRC) $(STARS_SRC) $(HYDRO_SRC)
if HAVEPYTHON
AM_SOURCES += logger_python_wrapper.c
endif

# Include files for distribution, not installation.
nobase_noinst_HEADERS = logger_hydro.h hydro/Gadget2/logger_hydro.h \
	logger_stars.h stars/Default/logger_stars.h \
	logger_gravity.h gravity/MultiSoftening/logger_gravity.h

# Sources and flags for regular library
liblogger_la_SOURCES = $(AM_SOURCES)
liblogger_la_CFLAGS = $(AM_CFLAGS)
liblogger_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS) $(BIN_LDFLAGS)

# Sources and flags for MPI library
# liblogger_mpi_la_SOURCES = $(AM_SOURCES)
# liblogger_mpi_la_CFLAGS = $(AM_CFLAGS) $(MPI_FLAGS)
# liblogger_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS)
# liblogger_mpi_la_SHORTNAME = mpi
# liblogger_mpi_la_LIBADD =
