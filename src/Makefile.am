
# This file is part of SWIFT.
# Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
AUTOMAKE_OPTIONS=gnu

# Add the debug flag to the whole thing
AM_CFLAGS = -g -O3 -std=gnu99 -Wall -Werror -ffast-math -fstrict-aliasing \
    -ftree-vectorize -funroll-loops $(SIMD_FLAGS) $(OPENMP_CFLAGS) \
    -DTIMER -DCOUNTER -DCPU_TPS=2.30e9 \
    # -fsanitize=address -fno-omit-frame-pointer
# AM_CFLAGS = -Wall -Werror $(OPENMP_CFLAGS) \
#     -DTIMER -DCOUNTER -DCPU_TPS=2.67e9

# Assign a "safe" version number
AM_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LDFLAGS) -version-info 0:0:0 # -fsanitize=address

# Build the libswiftsim library
lib_LTLIBRARIES = libswiftsim.la
# Build a MPI-enabled version too?
if HAVEMPI
lib_LTLIBRARIES += libswiftsim_mpi.la
endif

# List required headers
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h \
    engine.h swift.h serial_io.h timers.h debug.h scheduler.h proxy.h parallel_io.h \
    common_io.h multipole.h

# Common source files
AM_SOURCES = space.c runner.c queue.c task.c cell.c engine.c \
    serial_io.c timers.c debug.c scheduler.c proxy.c parallel_io.c \
    units.c common_io.c multipole.c
    
# Sources and flags for regular library
libswiftsim_la_SOURCES = $(AM_SOURCES)

# Sources and flags for MPI library
libswiftsim_mpi_la_SOURCES = $(AM_SOURCES)
libswiftsim_mpi_la_CFLAGS = $(AM_CFLAGS) -DWITH_MPI
libswiftsim_mpi_la_SHORTNAME = mpi

