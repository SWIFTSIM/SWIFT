# tHIS FIle is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

# Common flags
MYFLAGS = 

# Add the source directory and the non-standard paths to the included library headers to CFLAGS
AM_CFLAGS = -I$(top_srcdir)/src $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS)

AM_LDFLAGS = $(HDF5_LDFLAGS)

# Extra libraries.
EXTRA_LIBS = $(HDF5_LIBS) $(FFTW_LIBS) $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS) $(GSL_LIBS)

# MPI libraries.
MPI_LIBS = $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(METIS_INCS)

# Programs.
bin_PROGRAMS = swift

# Build MPI versions as well?
if HAVEMPI
bin_PROGRAMS += swift_mpi
endif

# engine_policy_setaffinity is available?
if HAVESETAFFINITY
ENGINE_POLICY_SETAFFINITY=| engine_policy_setaffinity
else
ENGINE_POLICY_SETAFFINITY=
endif

# Sources for swift
swift_SOURCES = main.c
swift_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_LDADD =  ../src/.libs/libswiftsim.a $(EXTRA_LIBS)

# Sources for swift_mpi, do we need an affinity policy for MPI?
swift_mpi_SOURCES = main.c
swift_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(MPI_LIBS) $(EXTRA_LIBS)

# Scripts to generate ICs
EXTRA_DIST = CoolingBox/coolingBox.yml CoolingBox/energy_plot.py CoolingBox/makeIC.py CoolingBox/run.sh \
	     EAGLE_6/eagle_6.yml EAGLE_6/getIC.sh EAGLE_6/README EAGLE_6/run.sh \
	     EAGLE_12/eagle_12.yml EAGLE_12/getIC.sh EAGLE_12/README EAGLE_12/run.sh \
	     EAGLE_25/eagle_25.yml EAGLE_25/getIC.sh EAGLE_25/README EAGLE_25/run.sh \
	     EAGLE_50/eagle_50.yml EAGLE_50/getIC.sh EAGLE_50/README EAGLE_50/run.sh \
	     EAGLE_100/eagle_100.yml EAGLE_100/getIC.sh EAGLE_100/README EAGLE_100/run.sh \
	     EAGLE_DMO_12/eagle_12.yml EAGLE_DMO_12/getIC.sh EAGLE_DMO_12/README EAGLE_DMO_12/run.sh \
	     EAGLE_DMO_25/eagle_25.yml EAGLE_DMO_25/getIC.sh EAGLE_DMO_25/README EAGLE_DMO_25/run.sh \
	     EAGLE_DMO_50/eagle_50.yml EAGLE_DMO_50/getIC.sh EAGLE_DMO_50/README EAGLE_DMO_50/run.sh \
	     EAGLE_DMO_100/eagle_100.yml EAGLE_DMO_100/getIC.sh EAGLE_DMO_100/README EAGLE_DMO_100/run.sh \
	     EvrardCollapse_3D/evrard.yml EvrardCollapse_3D/makeIC.py EvrardCollapse_3D/plotSolution.py EvrardCollapse_3D/run.sh EvrardCollapse_3D/getReference.sh \
	     ExternalPointMass/externalPointMass.yml ExternalPointMass/makeIC.py ExternalPointMass/run.sh ExternalPointMass/energy_plot.py \
	     GreshoVortex_2D/getGlass.sh GreshoVortex_2D/gresho.yml GreshoVortex_2D/makeIC.py GreshoVortex_2D/plotSolution.py GreshoVortex_2D/run.sh \
	     GreshoVortex_3D/getGlass.sh GreshoVortex_3D/gresho.yml GreshoVortex_3D/makeIC.py GreshoVortex_3D/plotSolution.py GreshoVortex_3D/run.sh \
	     HydrostaticHalo/README HydrostaticHalo/hydrostatic.yml HydrostaticHalo/makeIC.py HydrostaticHalo/run.sh \
	     HydrostaticHalo/density_profile.py HydrostaticHalo/velocity_profile.py HydrostaticHalo/internal_energy_profile.py HydrostaticHalo/test_energy_conservation.py \
	     InteractingBlastWaves_1D/run.sh InteractingBlastWaves_1D/makeIC.py InteractingBlastWaves_1D/plotSolution.py InteractingBlastWaves_1D/interactingBlastWaves.yml InteractingBlastWaves_1D/getReference.sh \
	     IsothermalPotential/README IsothermalPotential/run.sh IsothermalPotential/energy_plot.py IsothermalPotential/isothermal.yml IsothermalPotential/makeIC.py \
	     KelvinHelmholtz_2D/kelvinHelmholtz.yml KelvinHelmholtz_2D/makeIC.py KelvinHelmholtz_2D/plotSolution.py KelvinHelmholtz_2D/run.sh \
	     MultiTypes/makeIC.py  MultiTypes/multiTypes.yml MultiTypes/run.sh \
             Noh_1D/makeIC.py Noh_1D/noh.yml Noh_1D/plotSolution.py Noh_1D/run.sh \
             Noh_2D/makeIC.py Noh_2D/noh.yml Noh_2D/plotSolution.py Noh_2D/run.sh Noh_2D/getGlass.sh \
             Noh_3D/makeIC.py Noh_3D/noh.yml Noh_3D/plotSolution.py Noh_3D/run.sh Noh_3D/getGlass.sh \
	     PerturbedBox_2D/makeIC.py PerturbedBox_2D/perturbedPlane.yml \
	     PerturbedBox_3D/makeIC.py PerturbedBox_3D/perturbedBox.yml PerturbedBox_3D/run.sh \
	     SedovBlast_1D/makeIC.py SedovBlast_1D/plotSolution.py SedovBlast_1D/run.sh SedovBlast_1D/sedov.yml \
	     SedovBlast_2D/getGlass.sh SedovBlast_2D/makeIC.py SedovBlast_2D/plotSolution.py SedovBlast_2D/run.sh SedovBlast_2D/sedov.yml \
	     SedovBlast_3D/getGlass.sh SedovBlast_3D/makeIC.py SedovBlast_3D/plotSolution.py SedovBlast_3D/run.sh SedovBlast_3D/sedov.yml \
             SineWavePotential_1D/makeIC.py SineWavePotential_1D/plotSolution.py SineWavePotential_1D/run.sh SineWavePotential_1D/sineWavePotential.yml \
             SineWavePotential_2D/makeIC.py SineWavePotential_2D/plotSolution.py SineWavePotential_2D/run.sh SineWavePotential_2D/sineWavePotential.yml \
             SineWavePotential_3D/makeIC.py SineWavePotential_3D/plotSolution.py SineWavePotential_3D/run.sh SineWavePotential_3D/sineWavePotential.yml \
             SmallCosmoVolume/README SmallCosmoVolume/getIC.sh SmallCosmoVolume/run.sh SmallCosmoVolume/small_cosmo_volume.yml \
	     SodShock_1D/makeIC.py SodShock_1D/plotSolution.py SodShock_1D/run.sh SodShock_1D/sodShock.yml \
	     SodShock_2D/getGlass.sh SodShock_2D/makeIC.py SodShock_2D/plotSolution.py SodShock_2D/run.sh SodShock_2D/sodShock.yml \
	     SodShock_3D/getGlass.sh SodShock_3D/makeIC.py SodShock_3D/plotSolution.py SodShock_3D/run.sh SodShock_3D/sodShock.yml \
	     SquareTest_2D/makeIC.py SquareTest_2D/plotSolution.py SquareTest_2D/run.sh SquareTest_2D/square.yml \
	     UniformBox_2D/makeIC.py UniformBox_2D/run.sh UniformBox_2D/uniformPlane.yml \
	     UniformBox_3D/makeICbig.py UniformBox_3D/makeIC.py UniformBox_3D/run.sh UniformBox_3D/uniformBox.yml \
	     UniformDMBox/makeIC.py \
             ZeldovichPancake_3D/makeIC.py ZeldovichPancake_3D/zeldovichPancake.yml ZeldovichPancake_3D/run.sh ZeldovichPancake_3D/plotSolution.py

# Default parameter file
EXTRA_DIST += parameter_example.yml

# Scripts to plot task graphs
EXTRA_DIST += plot_tasks.py analyse_tasks.py process_plot_tasks_MPI process_plot_tasks

# Scripts to plot threadpool 'task' graphs
EXTRA_DIST += analyse_threadpool_tasks.py \
              plot_threadpool.py \
              process_plot_threadpool

# Script for scaling plot
EXTRA_DIST += plot_scaling_results.py \
              plot_scaling_results_breakdown.py

# Script for gravity accuracy
EXTRA_DIST += plot_gravity_checks.py
