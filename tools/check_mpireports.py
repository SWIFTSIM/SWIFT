#!/usr/bin/env python
"""
Usage:
    check_mpireportx.py [options] mpi-report...

Check any mpi reports from a step for any requests which are not fulfilled,
that any sends or receives that are not matched to a successful completion.
Note that a report with a final sum of memory equal to zero will not have any
of these!

This file is part of SWIFT.

Copyright (C) 2020 Peter W. Draper (p.w.draper@durham.ac.uk)
All Rights Reserved.

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
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Check MPI report")

parser.add_argument("input",
                    nargs="+",
                    metavar="mpi-report",
                    help="MPI report")
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    help="Verbose output",
    default=False,
    action="store_true",
)
args = parser.parse_args()
infiles = args.input

#  Indices for words in a line.
sticcol=0
eticcol=1
dticcol=2
stepcol=3
rankcol=4
otherrankcol=5
typecol=6
itypecol=7
subtypecol=8
isubtypecol=9
activationcol=10
tagcol=11
sizecol=12
sum=13

#  Keyed lines.
isends = {}
irecvs = {}
esends = {}
erecvs = {}

#  Gather keys from input file. We created dicts with matchable keys
#  for the sends and recvs, when they start and complete.
for f in infiles:
    if args.verbose:
        print "Processing: " + f
    with open(f, "r") as fp:
        for line in fp:
            if line[0] == '#':
                continue
            words = line.split()
            if words[activationcol] == "1":
                key = words[otherrankcol] + "/" + \
                      words[rankcol] + "/" + \
                      words[subtypecol] + "/" + \
                      words[tagcol] + "/" + \
                      words[sizecol]
                if words[typecol] == "send":
                    isends[key] = [line[:-1]]
                else:
                    irecvs[key] = [line[:-1]]
                    
            else:
                # Size will be negative.
                key = words[otherrankcol] + "/" + \
                      words[rankcol] + "/" + \
                      words[subtypecol] + "/" + \
                      words[tagcol] + "/" + \
                      words[sizecol][1:]
                if words[typecol] == "send":
                    esends[key] = [line[:-1]]
                else:
                    erecvs[key] = [line[:-1]]

#  Now report any uncompleted sends or receives.
print "# stic etic dtic step rank otherrank type itype subtype isubtype activation tag size sum"
for key in isends:
    if not key in esends:
        print isends[key][0]
for key in irecvs:
    if not key in erecvs:
        print irecvs[key][0]

sys.exit(0)
