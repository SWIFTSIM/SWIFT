
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

# Add the debug flag to the whole thing
AM_CFLAGS = -DTIMER -DCOUNTER -DCPU_TPS=2.30e9

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
    common_io.h multipole.h git_revision.h

# Common source files
AM_SOURCES = space.c runner.c queue.c task.c cell.c engine.c \
    serial_io.c timers.c debug.c scheduler.c proxy.c parallel_io.c \
    units.c common_io.c multipole.c git_revision.c

# Include files for distribution, not installation.
noinst_HEADERS = atomic.h cycle.h error.h inline.h kernel.h vector.h \
	runner_iact_legacy.h runner_iact_grav.h runner_doiact.h \
	runner_doiact_grav.h units.h

# Sources and flags for regular library
libswiftsim_la_SOURCES = $(AM_SOURCES)

# Sources and flags for MPI library
libswiftsim_mpi_la_SOURCES = $(AM_SOURCES)
libswiftsim_mpi_la_CFLAGS = $(AM_CFLAGS) -DWITH_MPI
libswiftsim_mpi_la_SHORTNAME = mpi


# Git revision
# Generate a file with the Git revision chain using the pre-defined command
git_revision.c: GIT_REV_CHAIN
	echo 'const char* git_revision(void) { const char* Git_Revision = "\' > git_revision.c
	@my_git_revision@                  >> git_revision.c
	echo '"; return Git_Revision; }'   >> git_revision.c

# If we are in a git directory save the revision number
if HAVE_GIT_WC
git_revision: GIT_REV_CHAIN
	@my_git_revision@ > git_revision
endif

GIT_REV_CHAIN:

BUILT_SOURCES=git_revision.c 
EXTRA_DIST=git_revision
