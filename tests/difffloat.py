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
# The (cube root of) the number of lines to check is provided as
# an optional 4th argument

file1 = sys.argv[1]
file2 = sys.argv[2]
number_to_check = -1

fileTol = ""
if len(sys.argv) >= 4:
    fileTol = sys.argv[3]

if len(sys.argv) >= 5:
    number_to_check = int(sys.argv[4])

# Get the particle properties being compared from the header.
with open(file1, 'r') as f:
  line = f.readline()
  if 'ID' in line:
    part_props = line.split()[1:]

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
    if n_linesTol != 3:
        print "Incorrect number of lines in tolerance file '%s'."%fileTol
    if n_columnsTol != n_columns:
        print "Incorrect number of columns in tolerance file '%s'."%fileTol

if fileTol == "":
    print "Absolute difference tolerance:", abs_tol
    print "Relative difference tolerance:", rel_tol
    absTol = ones(n_columns) * abs_tol
    relTol = ones(n_columns) * rel_tol
    limTol = zeros(n_columns)
else:
    print "Tolerances read from file"
    absTol = dataTol[0,:]
    relTol = dataTol[1,:]
    limTol = dataTol[2,:]

n_lines_to_check = 0
if number_to_check > 0:
    n_lines_to_check = number_to_check**3
    n_lines_to_check = min(n_lines_to_check, n_lines)
    print "Checking the first %d particles."%n_lines_to_check
else:
    n_lines_to_check = n_lines
    print "Checking all particles in the file."

error = False
for i in range(n_lines_to_check):
    for j in range(n_columns):

        abs_diff = abs(data1[i,j] - data2[i,j])

        sum = abs(data1[i,j] + data2[i,j])
        if sum > 0:
            rel_diff = abs(data1[i,j] - data2[i,j]) / sum
        else:
            rel_diff = 0.

        if( abs_diff > 1.1*absTol[j]):
            print "Absolute difference larger than tolerance (%e) for particle %d, column %s:"%(absTol[j], data1[i,0], part_props[j])
            print "%10s:           a = %e"%("File 1", data1[i,j])
            print "%10s:           b = %e"%("File 2", data2[i,j])
            print "%10s:       |a-b| = %e"%("Difference", abs_diff)
            print ""
            error = True

        if abs(data1[i,j]) + abs(data2[i,j]) < limTol[j] : continue

        if( rel_diff > 1.1*relTol[j]):
            print "Relative difference larger than tolerance (%e) for particle %d, column %s:"%(relTol[j], data1[i,0], part_props[j])
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
