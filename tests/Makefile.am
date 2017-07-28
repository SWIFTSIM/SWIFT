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
AM_CFLAGS = -I$(top_srcdir)/src $(HDF5_CPPFLAGS)

AM_LDFLAGS = ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(FFTW_LIBS)

# List of programs and scripts to run in the test suite
TESTS = testGreetings testMaths testReading.sh testSingle testKernel testSymmetry \
        testPair.sh testPairPerturbed.sh test27cells.sh test27cellsPerturbed.sh  \
        testParser.sh testSPHStep test125cells.sh test125cellsPerturbed.sh testFFT \
        testAdiabaticIndex testRiemannExact testRiemannTRRS testRiemannHLLC \
        testMatrixInversion testThreadpool testDump testLogger \
        testVoronoi1D testVoronoi2D testVoronoi3D \
	testPeriodicBC.sh testPeriodicBCPerturbed.sh

# List of test programs to compile
check_PROGRAMS = testGreetings testReading testSingle testTimeIntegration \
		 testSPHStep testPair test27cells test125cells testParser \
                 testKernel testFFT testInteractions testMaths \
                 testSymmetry testThreadpool benchmarkInteractions \
                 testAdiabaticIndex testRiemannExact testRiemannTRRS \
                 testRiemannHLLC testMatrixInversion testDump testLogger \
		 testVoronoi1D testVoronoi2D testVoronoi3D testPeriodicBC

# Rebuild tests when SWIFT is updated.
$(check_PROGRAMS): ../src/.libs/libswiftsim.a

# Sources for the individual programs
testGreetings_SOURCES = testGreetings.c

testMaths_SOURCES = testMaths.c

testReading_SOURCES = testReading.c

testSymmetry_SOURCES = testSymmetry.c

testTimeIntegration_SOURCES = testTimeIntegration.c

testSPHStep_SOURCES = testSPHStep.c

testSingle_SOURCES = testSingle.c

testPair_SOURCES = testPair.c

test27cells_SOURCES = test27cells.c

testPeriodicBC_SOURCES = testPeriodicBC.c

test125cells_SOURCES = test125cells.c

testParser_SOURCES = testParser.c

testKernel_SOURCES = testKernel.c

testFFT_SOURCES = testFFT.c

testInteractions_SOURCES = testInteractions.c

benchmarkInteractions_SOURCES = benchmarkInteractions.c

testAdiabaticIndex_SOURCES = testAdiabaticIndex.c

testRiemannExact_SOURCES = testRiemannExact.c

testRiemannTRRS_SOURCES = testRiemannTRRS.c

testRiemannHLLC_SOURCES = testRiemannHLLC.c

testMatrixInversion_SOURCES = testMatrixInversion.c

testVoronoi1D_SOURCES = testVoronoi1D.c

testVoronoi2D_SOURCES = testVoronoi2D.c

testVoronoi3D_SOURCES = testVoronoi3D.c

testThreadpool_SOURCES = testThreadpool.c

testDump_SOURCES = testDump.c

testLogger_SOURCES = testLogger.c

# Files necessary for distribution
EXTRA_DIST = testReading.sh makeInput.py testPair.sh testPairPerturbed.sh \
	     test27cells.sh test27cellsPerturbed.sh testParser.sh testPeriodicBC.sh \
	     testPeriodicBCPerturbed.sh test125cells.sh test125cellsPerturbed.sh testParserInput.yaml \
	     difffloat.py tolerance_125_normal.dat tolerance_125_perturbed.dat \
             tolerance_27_normal.dat tolerance_27_perturbed.dat tolerance_27_perturbed_h.dat \
	     tolerance_pair_normal.dat tolerance_pair_perturbed.dat \
	     fft_params.yml tolerance_periodic_BC_normal.dat tolerance_periodic_BC_perturbed.dat
