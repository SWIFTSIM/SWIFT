#!/usr/bin/env python
"""
Usage:
    analysedumpcells.py nx ny nx cell<1>.dat cell<2>.dat ...

Analyses a number of output files created by calls to the dumpCells() debug
function (presumably called in engine_step()) to output a list of active
top-level cells, identifying those that are on the edges of the volumes being
processed on by various nodes. The point is that these should be minimised to
reduce the MPI communications.

The nx, ny and nz arguments are the number of cells in the complete space,
we need these so that we can wrap around the edge of space correctly.

This file is part of SWIFT.
Copyright (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)

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
import pylab as pl
import numpy as np
import sys
import pandas

xcol = 0
ycol = 1
zcol = 2
xwcol = 3
ywcol = 4
zwcol = 5
supercol = 15
topcol = 16
activecol = 17
localcol = 19
mpicol = 20

#  Command-line arguments.
if len(sys.argv) < 5:
    print("usage: ", sys.argv[0], " nx ny nz cell1.dat cell2.dat ...")
    sys.exit(1)
nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])

print("# x y z onedge")
allactives = []
onedge = 0
tcount = 0
for i in range(4, len(sys.argv)):

    #  Read the file.
    data = pl.loadtxt(sys.argv[i])
    if len(data) == 0 or len(data) == 20:
        continue

    #  Select cells that are on the current rank and are top-level cells.
    rdata = data[data[:, localcol] == 1]
    tdata = rdata[rdata[:, topcol] == 1]

    #  Separation of the cells is in data.
    xwidth = tdata[0, xwcol]
    ywidth = tdata[0, ywcol]
    zwidth = tdata[0, zwcol]

    #  Fill space nx, ny,n nz with all toplevel cells and flag their active
    #  state.
    space = np.zeros((nx, ny, nz))
    actives = []
    for line in tdata:
        ix = int(np.rint(line[xcol] / xwidth))
        iy = int(np.rint(line[ycol] / ywidth))
        iz = int(np.rint(line[zcol] / zwidth))
        active = int(line[activecol])
        space[ix, iy, iz] = 1 + active
        tcount = tcount + 1
        if active == 1:
            actives.append([ix, iy, iz, line])

    #  Report all active cells and flag any without 26 neighbours. These are
    #  on the edge of the partition volume and will have foreign neighbour
    #  cells.
    for active in actives:
        count = 0
        for ii in [-1, 0, 1]:
            i = active[0] + ii
            if i < 0:
                i = i + nx
            elif i >= nx:
                i = i - nx

            for jj in [-1, 0, 1]:
                j = active[1] + jj
                if j < 0:
                    j = j + ny
                elif j >= ny:
                    j = j - ny

                for kk in [-1, 0, 1]:
                    k = active[2] + kk
                    if k < 0:
                        k = k + nz
                    elif k >= nz:
                        k = k - nz
                    if space[i, j, k] > 0:
                        count = count + 1
        if count < 27:
            onedge = onedge + 1
            print(active[3][0], active[3][1], active[3][2], 1)
        else:
            print(active[3][0], active[3][1], active[3][2], 0)

    allactives.extend(actives)

print("# top cells: ", tcount, " active: ", len(allactives), " on edge: ", onedge)

sys.exit(0)
