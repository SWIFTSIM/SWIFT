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
	sed -i "/master_doc = .\/*/c\master_doc = 'index_doxygen'" RTD/source/conf.py
	sed -i "/^#.*breathe/s/^#//" source/conf.py
	breathe-apidoc -f -g struct,file xml/ -o "RTD/source/Doxygen/" -p SWIFT
else
	sed -i "/master_doc = .\/*/c\master_doc = 'index'" RTD/source/conf.py
	sed -i "s/^[^\#]*breathe.*/\#&/" source/conf.py
endif

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

clean-local:
	rm -rf html/ latex/ man/

distclean: clean
