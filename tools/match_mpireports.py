#!/usr/bin/env python
"""
Usage:
    match_mpireports.py [options] mpi-reports...

Match the rows that sends start and recvs complete from a set of mpi-reports
of a single step, and output the matched rows to standard output. If captured
the output can be analysed to see how long the send to recvs took to complete.

This file is part of SWIFT.

Copyright (C) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
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
parser = argparse.ArgumentParser(description="Match MPI reports")

parser.add_argument("input",
                    nargs="+",
                    metavar="mpi-reports",
                    help="MPI reports")
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
sends = {}
recvs = {}

#  Gather keys from input files. We created dicts with matchable keys
#  for when sends start and recvs end. Other pairings are possible...
#  Note size of completion recv is negative.
for f in infiles:
    if args.verbose:
        print "Processing: " + f
    with open(f, "r") as fp:
        for line in fp:
            if line[0] == '#':
                continue
            words = line.split()
            if words[activationcol] == "1" and words[typecol] == "send":
                key = words[otherrankcol] + "/" + \
                      words[rankcol] + "/" + \
                      words[subtypecol] + "/" + \
                      words[tagcol] + "/" + \
                      words[sizecol]
                if not key in sends:
                    sends[key] = [line[:-1]]
                else:
                    sends[key].append(line[:-1])

            elif words[activationcol] == "0" and words[typecol] == "recv":
                key = words[rankcol] + "/" + \
                      words[otherrankcol] + "/" + \
                      words[subtypecol] + "/" + \
                      words[tagcol] + "/" + \
                      words[sizecol][1:]

                if not key in recvs:
                    recvs[key] = [line[:-1]]
                else:
                    recvs[key].append(line[:-1])

#  Now output. Note we could have unmatched recv keys, we don't check for that.
print "# send_stic send_etic send_dtic send_step send_rank send_otherrank " + \
        "send_type send_itype send_subtype send_isubtype send_activation " + \
        "send_tag send_size send_sum " + \
        "recv_stic recv_etic recv_dtic recv_step recv_rank recv_otherrank " + \
        "recv_type recv_itype recv_subtype recv_isubtype recv_activation " + \
        "recv_tag recv_size recv_sum "
for key in sends:
    if key in recvs:
        if len(sends[key]) == 1 and len(recvs[key]) == 1:
            print sends[key][0], recvs[key][0]
        else:
            print "# ERROR: found ", len(sends[key]), "/", len(recvs[key]), " matches for key: ", key, " should be 1/1"
    else:
        print "# ERROR: missing recv key: ", key


sys.exit(0)
