# This file is part of SWIFT.
# Copyright (c) 2015 matthieu.schaller@durham.ac.uk.
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

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -I../src $(HDF5_CPPFLAGS) -DTIMER


AM_LDFLAGS = ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

# List of programs and scripts to run in the test suite
TESTS = testGreetings testReading.sh testSingle testPair.sh testPairPerturbed.sh \
	test27cells.sh test27cellsPerturbed.sh testParser.sh

# List of test programs to compile
check_PROGRAMS = testGreetings testReading testSingle testTimeIntegration \
		 testSPHStep testPair test27cells testParser

# Sources for the individual programs
testGreetings_SOURCES = testGreetings.c

testReading_SOURCES = testReading.c

testTimeIntegration_SOURCES = testTimeIntegration.c

testSPHStep_SOURCES = testSPHStep.c

testSingle_SOURCES = testSingle.c

testPair_SOURCES = testPair.c

test27cells_SOURCES = test27cells.c

testParser_SOURCES = testParser.c

# Files necessary for distribution
EXTRA_DIST = testReading.sh makeInput.py testPair.sh testPairPerturbed.sh \
	     test27cells.sh test27cellsPerturbed.sh tolerance.dat testParser.sh \ 
	     testParserInput.yaml
