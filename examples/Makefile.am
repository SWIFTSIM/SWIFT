
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

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -g -Wall -Werror -I../src $(OPENMP_CFLAGS) -DCPU_TPS=2.67e9

AM_LDFLAGS =

# Set-up the library
bin_PROGRAMS = test test_fixed

# Sources for test
test_SOURCES = test.c
test_CFLAGS = -DCOUNTER -DTIMER $(AM_CFLAGS)
test_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS)

# Sources for test_fixed
test_fixed_SOURCES = test.c
test_fixed_CFLAGS = -DCOUNTER -DTIMER $(AM_CFLAGS) -DENGINE_FLAGS=engine_policy_fixdt
test_fixed_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS)

