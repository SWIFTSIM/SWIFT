
# This file is part of GadgetSMP.
# Coypright (c) 2012 pedro.gonnet@durham.ac.uk.
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
AM_CFLAGS = -g -O3 -Wall -Werror -ffast-math -fstrict-aliasing $(SIMD_FLAGS) $(CFLAGS) $(OPENMP_CFLAGS) -DTIMER -DCOUNTER -DCPU_TPS=2.67e9

# Assign a "safe" version number
AM_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) -version-info 0:0:0

# Build the libgadgetsmp library
lib_LTLIBRARIES = libgadgetsmp.la
libgadgetsmp_la_SOURCES = space.c runner.c queue.c task.c cell.c engine.c

# List required headers
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h \
    engine.h gadgetsmp.h

