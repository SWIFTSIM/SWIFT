###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

from numpy import *
import sys

abs_tol = 1e-7
rel_tol = 1e-7

# Compares the content of two ASCII tables of floats line by line and
# reports all differences beyond the given tolerances
# Comparisons are done both in absolute and relative terms

# Individual tolerances for each column can be provided in a file

file1 = sys.argv[1]
file2 = sys.argv[2]
fileTol = ""

if len(sys.argv) == 4:
    fileTol = sys.argv[3]

data1 = loadtxt(file1)
data2 = loadtxt(file2)
if fileTol != "":
    dataTol = loadtxt(fileTol)
    n_linesTol = shape(dataTol)[0]
    n_columnsTol = shape(dataTol)[1]


if shape(data1) != shape(data2):
    print "Non-matching array sizes in the files", file1, "and", file2, "."
    sys.exit(1)

n_lines = shape(data1)[0]
n_columns = shape(data1)[1]

if fileTol != "":
    if n_linesTol != 2:
        print "Incorrect number of lines in tolerance file '%s'."%fileTol
    if n_columnsTol != n_columns:
        print "Incorrect number of columns in tolerance file '%s'."%fileTol

if fileTol == "":
    print "Absolute difference tolerance:", abs_tol
    print "Relative difference tolerance:", rel_tol
    absTol = ones(n_columns) * abs_tol
    relTol = ones(n_columns) * rel_tol
else:
    print "Tolerances read from file"
    absTol = dataTol[0,:]
    relTol = dataTol[1,:]

error = False
for i in range(n_lines):
    for j in range(n_columns):

        abs_diff = abs(data1[i,j] - data2[i,j])

        sum = abs(data1[i,j] + data2[i,j])
        if sum > 0:
            rel_diff = abs(data1[i,j] - data2[i,j]) / sum
        else:
            rel_diff = 0.

        if( abs_diff > absTol[j]):
            print "Absolute difference larger than tolerance (%e) on line %d, column %d:"%(absTol[j], i,j)
            print "%10s:           a = %e"%("File 1", data1[i,j])
            print "%10s:           b = %e"%("File 2", data2[i,j])
            print "%10s:       |a-b| = %e"%("Difference", abs_diff)
            print ""
            error = True

        if( rel_diff > relTol[j]):
            print "Relative difference larger than tolerance (%e) on line %d, column %d:"%(relTol[j], i,j)
            print "%10s:           a = %e"%("File 1", data1[i,j])
            print "%10s:           b = %e"%("File 2", data2[i,j])
            print "%10s: |a-b|/|a+b| = %e"%("Difference", rel_diff)
            print ""
            error = True


if error:
    exit(1)
else:
    print "No differences found"
    exit(0)
