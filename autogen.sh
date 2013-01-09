#! /bin/sh
# autogen.sh

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

# make libtool and stuff
libtoolize --force --copy

# run aclocal
aclocal -I m4/

# run autoconf
autoconf -I m4

# run autoheader to generate config.h.in
autoheader -I m4

# run automake
automake --add-missing

