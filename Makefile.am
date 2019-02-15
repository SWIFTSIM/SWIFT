# This file is part of SWIFT.
# Copyright (c) 2012 pedro.gonnet@durham.ac.uk
#               2015 matthieu.schaller@durham.ac.uk
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
ACLOCAL_AMFLAGS = -I m4

# Show the way...
SUBDIRS = src argparse examples doc tests tools
if HAVEEAGLECOOLING
SUBDIRS += examples/Cooling/CoolingRates
endif

# Non-standard files that should be part of the distribution.
EXTRA_DIST = INSTALL.swift .clang-format format.sh
