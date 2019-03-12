# This file is part of SWIFT.
# Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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

lib_LTLIBRARIES = libargparse.la

include_HEADERS = argparse.h

AM_SOURCES = argparse.c

# Sources and flags for regular library
libargparse_la_SOURCES = $(AM_SOURCES)
libargparse_la_CFLAGS = $(AM_CFLAGS)
libargparse_la_LDFLAGS = $(AM_LDFLAGS)

EXTRA_DIST = LICENSE README.md
