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

doxyfile.stamp:
if HAVE_DOXYGEN
	@DX_ENV@ @DX_DOXYGEN@ Doxyfile
	echo `date` > doxyfile.stamp
	sed -i "/^#.*breathe/s/^#//" @srcdir@/RTD/source/conf.py
	cp @srcdir@/RTD/source/index_doxygen.rst @srcdir@/RTD/source/index.rst
	breathe-apidoc -f -g struct,file xml/ -o "@srcdir@/RTD/source/Doxygen/" -p SWIFT
else
	cp @srcdir@/RTD/source/index_no_doxygen.rst @srcdir@/RTD/source/index.rst
	sed -i "s/^[^\#]*breathe.*/\#&/" @srcdir@/RTD/source/conf.py
endif

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

clean-local:
	rm -rf html/ latex/ man/ xml/

distclean: clean
