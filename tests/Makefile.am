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

# Add the source directory and the non-standard paths to the included library headers to CFLAGS
AM_CFLAGS = -I$(top_srcdir)/src $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS)

AM_LDFLAGS = ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) $(GSL_LIBS) $(PROFILER_LIBS)

# List of programs and scripts to run in the test suite
TESTS = testGreetings testMaths testReading.sh testSingle testKernel \
        testActivePair.sh test27cells.sh test27cellsPerturbed.sh  \
        testParser.sh testSPHStep test125cells.sh test125cellsPerturbed.sh testFFT \
        testAdiabaticIndex testRandom \
        testMatrixInversion testThreadpool testDump testLogger testInteractions.sh \
        testVoronoi1D testVoronoi2D testVoronoi3D testGravityDerivatives \
	testPeriodicBC.sh testPeriodicBCPerturbed.sh testPotentialSelf \
	testPotentialPair testEOS testUtilities testSelectOutput.sh \
	testCbrt testCosmology testOutputList testFormat.sh \
	test27cellsStars.sh test27cellsStarsPerturbed.sh

# List of test programs to compile
check_PROGRAMS = testGreetings testReading testSingle testTimeIntegration \
		 testSPHStep testActivePair test27cells test27cells_subset test125cells testParser \
                 testKernel testFFT testInteractions testMaths testRandom \
                 testSymmetry testThreadpool \
                 testAdiabaticIndex testRiemannExact testRiemannTRRS \
                 testRiemannHLLC testMatrixInversion testDump testLogger \
		 testVoronoi1D testVoronoi2D testVoronoi3D testPeriodicBC \
		 testGravityDerivatives testPotentialSelf testPotentialPair testEOS testUtilities \
		 testSelectOutput testCbrt testCosmology testOutputList test27cellsStars \
		 test27cellsStars_subset testCooling testHashmap

# Rebuild tests when SWIFT is updated.
$(check_PROGRAMS): ../src/.libs/libswiftsim.a

# Sources for the individual programs
testGreetings_SOURCES = testGreetings.c

testMaths_SOURCES = testMaths.c

testRandom_SOURCES = testRandom.c

testReading_SOURCES = testReading.c

testSelectOutput_SOURCES = testSelectOutput.c

testCosmology_SOURCES = testCosmology.c

testOutputList_SOURCES = testOutputList.c

testSymmetry_SOURCES = testSymmetry.c

# Added because of issues using memcmp on clang 4.x
testSymmetry_CFLAGS = $(AM_CFLAGS) -fno-builtin-memcmp

testTimeIntegration_SOURCES = testTimeIntegration.c

testSPHStep_SOURCES = testSPHStep.c

testSingle_SOURCES = testSingle.c

testActivePair_SOURCES = testActivePair.c

test27cells_SOURCES = test27cells.c

test27cells_subset_SOURCES = test27cells.c

test27cells_subset_CFLAGS = $(AM_CFLAGS) -DTEST_DOSELF_SUBSET -DTEST_DOPAIR_SUBSET

test27cellsStars_SOURCES = test27cellsStars.c

test27cellsStars_subset_SOURCES = test27cellsStars.c

test27cellsStars_subset_CFLAGS = $(AM_CFLAGS) -DTEST_DOSELF_SUBSET -DTEST_DOPAIR_SUBSET

testPeriodicBC_SOURCES = testPeriodicBC.c

test125cells_SOURCES = test125cells.c

testParser_SOURCES = testParser.c

testKernel_SOURCES = testKernel.c

testFFT_SOURCES = testFFT.c

testInteractions_SOURCES = testInteractions.c

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

testGravityDerivatives_SOURCES = testGravityDerivatives.c

testPotentialSelf_SOURCES = testPotentialSelf.c

testPotentialPair_SOURCES = testPotentialPair.c

testEOS_SOURCES = testEOS.c

testUtilities_SOURCES = testUtilities.c

testCooling_SOURCES = testCooling.c

testHashmap_SOURCES = testHashmap.c

# Files necessary for distribution
EXTRA_DIST = testReading.sh makeInput.py testActivePair.sh \
	     test27cells.sh test27cellsPerturbed.sh testParser.sh testPeriodicBC.sh \
	     testPeriodicBCPerturbed.sh test125cells.sh test125cellsPerturbed.sh testParserInput.yaml \
	     difffloat.py tolerance_125_normal.dat tolerance_125_perturbed.dat \
             tolerance_27_normal.dat tolerance_27_perturbed.dat tolerance_27_perturbed_h.dat tolerance_27_perturbed_h2.dat \
	     tolerance_testInteractions.dat tolerance_pair_active.dat tolerance_pair_force_active.dat \
	     fft_params.yml tolerance_periodic_BC_normal.dat tolerance_periodic_BC_perturbed.dat \
	     testEOS.sh testEOS_plot.sh testSelectOutput.sh selectOutput.yml \
             output_list_params.yml output_list_time.txt output_list_redshift.txt \
             output_list_scale_factor.txt testEOS.sh testEOS_plot.sh \
	     test27cellsStars.sh test27cellsStarsPerturbed.sh star_tolerance_27_normal.dat \
	     star_tolerance_27_perturbed.dat star_tolerance_27_perturbed_h.dat star_tolerance_27_perturbed_h2.dat
