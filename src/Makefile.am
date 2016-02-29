
# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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
AM_CFLAGS = -DTIMER -DCOUNTER -DCPU_TPS=2.30e9 $(HDF5_CPPFLAGS)

# Assign a "safe" version number
AM_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LDFLAGS) -version-info 0:0:0 # -fsanitize=address

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Build the libswiftsim library
lib_LTLIBRARIES = libswiftsim.la
# Build a MPI-enabled version too?
if HAVEMPI
lib_LTLIBRARIES += libswiftsim_mpi.la
endif

# List required headers
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h \
    engine.h swift.h serial_io.h timers.h debug.h scheduler.h proxy.h parallel_io.h \
    common_io.h single_io.h multipole.h map.h tools.h

# Common source files
AM_SOURCES = space.c runner.c queue.c task.c cell.c engine.c \
    serial_io.c timers.c debug.c scheduler.c proxy.c parallel_io.c \
    units.c common_io.c single_io.c multipole.c version.c map.c \
    kernel.c tools.c part.c

# Include files for distribution, not installation.
nobase_noinst_HEADERS = approx_math.h atomic.h cycle.h error.h inline.h kernel.h vector.h \
		 runner_doiact.h runner_doiact_grav.h units.h intrinsics.h minmax.h \
		 gravity.h \
		 gravity/Default/gravity.h gravity/Default/runner_iact_grav.h \
		 gravity/Default/gravity_part.h \
	 	 hydro.h hydro_io.h \
		 hydro/Minimal/hydro.h hydro/Minimal/hydro_iact.h hydro/Minimal/hydro_io.h \
                 hydro/Minimal/hydro_debug.h hydro/Minimal/hydro_part.h \
		 hydro/Default/hydro.h hydro/Default/hydro_iact.h hydro/Default/hydro_io.h \
                 hydro/Default/hydro_debug.h hydro/Default/hydro_part.h \
		 hydro/Gadget2/hydro.h hydro/Gadget2/hydro_iact.h hydro/Gadget2/hydro_io.h \
                 hydro/Gadget2/hydro_debug.h hydro/Gadget2/hydro_part.h \
		 hydro/Gizmo/hydro.h hydro/Gizmo/hydro_iact.h hydro/Gizmo/hydro_io.h \
                 hydro/Gizmo/hydro_debug.h hydro/Gizmo/hydro_part.h \
	         riemann.h \
		 riemann/riemann_hllc.h riemann/riemann_trrs.h riemann/riemann_exact.h

# Sources and flags for regular library
libswiftsim_la_SOURCES = $(AM_SOURCES)

# Sources and flags for MPI library
libswiftsim_mpi_la_SOURCES = $(AM_SOURCES)
libswiftsim_mpi_la_CFLAGS = $(AM_CFLAGS) -DWITH_MPI $(METIS_INCS)
libswiftsim_mpi_la_SHORTNAME = mpi


# Versioning. If any sources change then update the version.h file with
# the current git revision and package version.
# May have a checkout without a version.h file and no git command (tar/zip
# download), allow that, but make sure we know it.
version.h: version.h.in $(AM_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	if test "X$(GIT_CMD)" != "X"; then \
	    GIT_REVISION=`$(GIT_CMD) describe --abbrev=8  --always --tags --dirty`; \
	    GIT_BRANCH=`$(GIT_CMD) branch | sed -n 's/^\* \(.*\)/\1/p'`; \
	    sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	        -e "s,@GIT_REVISION\@,$${GIT_REVISION}," \
	        -e "s|@GIT_BRANCH\@|$${GIT_BRANCH}|" version.h.in > version.h; \
	else \
	    if test ! -f version.h; then \
	        sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	            -e "s,@GIT_REVISION\@,unknown," \
		    -e "s,@GIT_BRANCH\@,unknown," version.h.in > version.h; \
	    fi; \
	fi

#  Make sure version.h is built first.
BUILT_SOURCES = version.h

#  And distribute the built files.
EXTRA_DIST = version.h version.h.in
