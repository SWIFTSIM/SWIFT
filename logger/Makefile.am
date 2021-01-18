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

GRAVITY_SRC = gravity/MultiSoftening/logger_gravity.c ../src/gravity/MultiSoftening/gravity_logger.c
# Include the stars modules
if HAVE_STARS_DEFAULT
STARS_SRC = stars/Default/logger_stars.c  ../src/stars/Default/stars_logger.c
endif
if HAVE_STARS_GEAR
STARS_SRC = stars/GEAR/logger_stars.c ../src/stars/GEAR/stars_logger.c
endif

# Include the hydro modules
if HAVE_SPHENIX
HYDRO_SRC = hydro/SPHENIX/logger_hydro.c ../src/hydro/SPHENIX/hydro_logger.c
endif
if HAVE_GADGET2
HYDRO_SRC = hydro/Gadget2/logger_hydro.c ../src/hydro/Gadget2/hydro_logger.c
endif

# Include the chemistry modules
if HAVE_CHEMISTRY_NONE
CHEMISTRY_SRC = chemistry/none/logger_chemistry.c ../src/chemistry/none/chemistry_logger.c
endif
if HAVE_CHEMISTRY_GEAR
CHEMISTRY_SRC = chemistry/GEAR/logger_chemistry.c ../src/chemistry/GEAR/chemistry_logger.c
endif

# Include the star formation modules
if HAVE_STAR_FORMATION_DEFAULT
STAR_FORMATION_SRC = star_formation/none/logger_star_formation.c ../src/star_formation/none/star_formation_particle_logger.c
endif
if HAVE_STAR_FORMATION_GEAR
STAR_FORMATION_SRC = star_formation/GEAR/logger_star_formation.c ../src/star_formation/GEAR/star_formation_particle_logger.c
endif

# subdirectories
SUBDIRS = tests

# List required headers
include_HEADERS = logger_header.h logger_loader_io.h logger_particle.h logger_time.h logger_tools.h logger_reader.h \
	logger_logfile.h logger_index.h quick_sort.h logger_python_tools.h logger_interpolation.h

# Common source files
AM_SOURCES = logger_header.c logger_loader_io.c logger_particle.c logger_time.c logger_tools.c logger_reader.c
AM_SOURCES += logger_logfile.c logger_index.c quick_sort.c
AM_SOURCES += $(GRAVITY_SRC) $(STARS_SRC) $(HYDRO_SRC) $(CHEMISTRY_SRC) $(STAR_FORMATION_SRC)

if HAVEPYTHON
AM_SOURCES += logger_python_wrapper.c
endif

if HAVEPYTHON
# The main python structure is dependent on the version, thus we need to
# remove this one in order to be safe.
PYTHON_EXTRA_COMPILER_FLAG = -Wno-missing-field-initializers
endif

# Include files for distribution, not installation.
nobase_noinst_HEADERS = logger_hydro.h hydro/Gadget2/logger_hydro.h hydro/SPHENIX/logger_hydro.h
nobase_noinst_HEADERS += logger_stars.h  stars/Default/logger_stars.h stars/GEAR/logger_stars.h
nobase_noinst_HEADERS += logger_gravity.h gravity/MultiSoftening/logger_gravity.h
nobase_noinst_HEADERS += logger_chemistry.h chemistry/none/logger_chemistry.h chemistry/GEAR/logger_chemistry.h
nobase_noinst_HEADERS += logger_star_formation.h star_formation/none/logger_star_formation.h
nobase_noinst_HEADERS += star_formation/GEAR/logger_star_formation.h


# Sources and flags for regular library
liblogger_la_SOURCES = $(AM_SOURCES)
liblogger_la_CFLAGS = $(AM_CFLAGS) $(PYTHON_EXTRA_COMPILER_FLAG)
liblogger_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS) $(BIN_LDFLAGS)

# Sources and flags for MPI library
# liblogger_mpi_la_SOURCES = $(AM_SOURCES)
# liblogger_mpi_la_CFLAGS = $(AM_CFLAGS) $(MPI_FLAGS) $(PYTHON_EXTRA_COMPILER_FLAG)
# liblogger_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS)
# liblogger_mpi_la_SHORTNAME = mpi
# liblogger_mpi_la_LIBADD =
