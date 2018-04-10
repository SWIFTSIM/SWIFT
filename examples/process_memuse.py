#!/usr/bin/env python
"""
Usage:
    process_memuse.py output.dat

Parse the output of a run of SWIFT to convert the memuse output strings
into memory use per rank per step.

This file is part of SWIFT.
Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
from collections import OrderedDict

#  Command-line arguments.
if len(sys.argv) != 2:
    print "usage: ", sys.argv[0], " output.dat"
    sys.exit(1)

memuse = OrderedDict()

with open(sys.argv[1]) as infile:
    for line in infile:
        if ":memuse:" in line:
            words = line.split()
            key = words[3] + ":" + words[4]
            if key not in memuse:
                memuse[key] = 0
            memuse[key] = memuse[key] + int(words[-1])

print "# step label function memuse"
for key, value in memuse.items():
    words = key.split(":", 2)
    print words[0], words[1], words[2], value

sys.exit(0)
