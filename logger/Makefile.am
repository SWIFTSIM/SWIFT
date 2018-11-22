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

# Add the non-standard paths to the included library headers
AM_CFLAGS = $(PYTHON_INCS)

# Assign a "safe" version number
AM_LDFLAGS = -version-info 0:0:0

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Additional dependencies for shared libraries.
EXTRA_LIBS = $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(PYTHON_LIBS)

# MPI libraries.
# MPI_LIBS = $(MPI_THREAD_LIBS)
# MPI_FLAGS = -DWITH_MPI

# Build the libswiftsim library
lib_LTLIBRARIES = libswiftlogger.la
# Build a MPI-enabled version too?
# if HAVEMPI
# lib_LTLIBRARIES += libswiftlogger_mpi.la
# endif

# List required headers
include_HEADERS = header.h io.h particle.h timeline.h logger_tools.h

# Common source files
AM_SOURCES = header.c io.c logger_loader.c particle.c \
	timeline.c logger_tools.c

# Include files for distribution, not installation.
nobase_noinst_HEADERS = 

# Sources and flags for regular library
libswiftlogger_la_SOURCES = $(AM_SOURCES)
libswiftlogger_la_CFLAGS = $(AM_CFLAGS)
libswiftlogger_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS)

# Sources and flags for MPI library
# libswiftlogger_mpi_la_SOURCES = $(AM_SOURCES)
# libswiftlogger_mpi_la_CFLAGS = $(AM_CFLAGS) $(MPI_FLAGS)
# libswiftlogger_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS)
# libswiftlogger_mpi_la_SHORTNAME = mpi
# libswiftlogger_mpi_la_LIBADD =


# Versioning. If any sources change then update the version_string.h file with
# the current git revision and package version.
# May have a checkout without a version_string.h file and no git command (tar/zip
# download), allow that, but make sure we know it.
version_string.h: ../src/version_string.h.in $(AM_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	if test "X$(GIT_CMD)" != "X"; then \
	    GIT_REVISION=`$(GIT_CMD) describe --abbrev=8  --always --tags --dirty`; \
	    GIT_BRANCH=`$(GIT_CMD) branch | sed -n 's/^\* \(.*\)/\1/p'`; \
            GIT_DATE=`$(GIT_CMD) log -1 --format=%ci`; \
	    sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	        -e "s,@GIT_REVISION\@,$${GIT_REVISION}," \
	        -e "s|@GIT_BRANCH\@|$${GIT_BRANCH}|" \
	        -e "s|@GIT_DATE\@|$${GIT_DATE}|" \
	        -e "s|@SWIFT_CFLAGS\@|$(CFLAGS)|" $< > version_string.h; \
	else \
	    if test ! -f version_string.h; then \
	        sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	            -e "s,@GIT_REVISION\@,unknown," \
		    -e "s,@GIT_BRANCH\@,unknown," \
		    -e "s,@GIT_DATE\@,unknown," \
	            -e "s|@SWIFT_CFLAGS\@|$(CFLAGS)|" $< > version_string.h; \
	    fi; \
	fi

#  Make sure version_string.h is built first.
BUILT_SOURCES = version_string.h

#  And distribute the built files.
EXTRA_DIST = version_string.h ../src/version_string.h.in
