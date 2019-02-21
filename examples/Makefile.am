# This file is part of SWIFT.
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
AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/argparse $(HDF5_CPPFLAGS) \
	$(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS) $(GRACKLE_INCS)

AM_LDFLAGS = $(HDF5_LDFLAGS)

# Extra libraries.
EXTRA_LIBS = $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(PROFILER_LIBS) \
	$(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) \
	$(VELOCIRAPTOR_LIBS) $(GSL_LIBS)

# MPI libraries.
MPI_LIBS = $(PARMETIS_LIBS) $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(PARMETIS_INCS) $(METIS_INCS)

# Programs.
bin_PROGRAMS = swift

bin_PROGRAMS += fof

# Build MPI versions as well?
if HAVEMPI
bin_PROGRAMS += swift_mpi
bin_PROGRAMS += fof_mpi
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
swift_LDADD =  ../src/.libs/libswiftsim.a ../argparse/.libs/libargparse.a $(EXTRA_LIBS)

# Sources for swift_mpi, do we need an affinity policy for MPI?
swift_mpi_SOURCES = main.c
swift_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a ../argparse/.libs/libargparse.a $(MPI_LIBS) $(EXTRA_LIBS)

# Sources for fof
fof_SOURCES = main_fof.c
fof_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
fof_LDADD =  ../src/.libs/libswiftsim.a $(EXTRA_LIBS)

# Sources for fof_mpi, do we need an affinity policy for MPI?
fof_mpi_SOURCES = main_fof.c
fof_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
fof_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(MPI_LIBS) $(EXTRA_LIBS)

# Scripts to generate ICs
EXTRA_DIST = Cooling/CoolingBox/coolingBox.yml Cooling/CoolingBox/energy_plot.py Cooling/CoolingBox/makeIC.py Cooling/CoolingBox/run.sh \
             Cosmology/ConstantCosmoVolume/run.sh Cosmology/ConstantCosmoVolume/makeIC.py Cosmology/ConstantCosmoVolume/plotSolution.py Cosmology/ConstantCosmoVolume/constant_volume.yml \
             Cosmology/ZeldovichPancake_3D/makeIC.py Cosmology/ZeldovichPancake_3D/zeldovichPancake.yml Cosmology/ZeldovichPancake_3D/run.sh Cosmology/ZeldovichPancake_3D/plotSolution.py \
             EAGLE_low_z/EAGLE_6/eagle_6.yml EAGLE_low_z/EAGLE_6/getIC.sh EAGLE_low_z/EAGLE_6/README EAGLE_low_z/EAGLE_6/run.sh \
	     EAGLE_low_z/EAGLE_12/eagle_12.yml EAGLE_low_z/EAGLE_12/getIC.sh EAGLE_low_z/EAGLE_12/README EAGLE_low_z/EAGLE_12/run.sh \
	     EAGLE_low_z/EAGLE_25/eagle_25.yml EAGLE_low_z/EAGLE_25/getIC.sh EAGLE_low_z/EAGLE_25/README EAGLE_low_z/EAGLE_25/run.sh \
	     EAGLE_low_z/EAGLE_50/eagle_50.yml EAGLE_low_z/EAGLE_50/getIC.sh EAGLE_low_z/EAGLE_50/README EAGLE_low_z/EAGLE_50/run.sh \
	     EAGLE_low_z/EAGLE_100/eagle_100.yml EAGLE_low_z/EAGLE_100/getIC.sh EAGLE_low_z/EAGLE_100/README EAGLE_low_z/EAGLE_100/run.sh \
	     EAGLE_DMO_low_z/EAGLE_DMO_12/eagle_12.yml EAGLE_DMO_low_z/EAGLE_DMO_12/getIC.sh EAGLE_DMO_low_z/EAGLE_DMO_12/README EAGLE_DMO_low_z/EAGLE_DMO_12/run.sh \
	     EAGLE_DMO_low_z/EAGLE_DMO_25/eagle_25.yml EAGLE_DMO_low_z/EAGLE_DMO_25/getIC.sh EAGLE_DMO_low_z/EAGLE_DMO_25/README EAGLE_DMO_low_z/EAGLE_DMO_25/run.sh \
	     EAGLE_DMO_low_z/EAGLE_DMO_50/eagle_50.yml EAGLE_DMO_low_z/EAGLE_DMO_50/getIC.sh EAGLE_DMO_low_z/EAGLE_DMO_50/README EAGLE_DMO_low_z/EAGLE_DMO_50/run.sh \
	     EAGLE_DMO_low_z/EAGLE_DMO_100/eagle_100.yml EAGLE_DMO_low_z/EAGLE_DMO_100/getIC.sh EAGLE_DMO_low_z/EAGLE_DMO_100/README EAGLE_DMO_low_z/EAGLE_DMO_100/run.sh \
	     GravityTests/ExternalPointMass/externalPointMass.yml GravityTests/ExternalPointMass/makeIC.py GravityTests/ExternalPointMass/run.sh GravityTests/ExternalPointMass/energy_plot.py \
	     GravityTests/HydrostaticHalo/README GravityTests/HydrostaticHalo/hydrostatic.yml GravityTests/HydrostaticHalo/makeIC.py GravityTests/HydrostaticHalo/run.sh \
	     GravityTests/HydrostaticHalo/density_profile.py GravityTests/HydrostaticHalo/velocity_profile.py GravityTests/HydrostaticHalo/internal_energy_profile.py GravityTests/HydrostaticHalo/test_energy_conservation.py \
             GravityTests/IsothermalPotential/README GravityTests/IsothermalPotential/run.sh GravityTests/IsothermalPotential/energy_plot.py GravityTests/IsothermalPotential/isothermal.yml GravityTests/IsothermalPotential/makeIC.py \
	     HydroTests/GreshoVortex_2D/getGlass.sh HydroTests/GreshoVortex_2D/gresho.yml HydroTests/GreshoVortex_2D/makeIC.py HydroTests/GreshoVortex_2D/plotSolution.py HydroTests/GreshoVortex_2D/run.sh \
	     HydroTests/GreshoVortex_3D/getGlass.sh HydroTests/GreshoVortex_3D/gresho.yml HydroTests/GreshoVortex_3D/makeIC.py HydroTests/GreshoVortex_3D/plotSolution.py HydroTests/GreshoVortex_3D/run.sh \
	     HydroTests/EvrardCollapse_3D/evrard.yml HydroTests/EvrardCollapse_3D/makeIC.py HydroTests/EvrardCollapse_3D/plotSolution.py HydroTests/EvrardCollapse_3D/run.sh HydroTests/EvrardCollapse_3D/getReference.sh \
	     HydroTests/InteractingBlastWaves_1D/run.sh HydroTests/InteractingBlastWaves_1D/makeIC.py HydroTests/InteractingBlastWaves_1D/plotSolution.py HydroTests/InteractingBlastWaves_1D/interactingBlastWaves.yml HydroTests/InteractingBlastWaves_1D/getReference.sh \
	     HydroTests/KelvinHelmholtz_2D/kelvinHelmholtz.yml HydroTests/KelvinHelmholtz_2D/makeIC.py HydroTests/KelvinHelmholtz_2D/plotSolution.py HydroTests/KelvinHelmholtz_2D/run.sh \
             HydroTests/Noh_1D/makeIC.py HydroTests/Noh_1D/noh.yml HydroTests/Noh_1D/plotSolution.py HydroTests/Noh_1D/run.sh \
             HydroTests/Noh_2D/makeIC.py HydroTests/Noh_2D/noh.yml HydroTests/Noh_2D/plotSolution.py HydroTests/Noh_2D/run.sh HydroTests/Noh_2D/getGlass.sh \
             HydroTests/Noh_3D/makeIC.py HydroTests/Noh_3D/noh.yml HydroTests/Noh_3D/plotSolution.py HydroTests/Noh_3D/run.sh HydroTests/Noh_3D/getGlass.sh \
	     HydroTests/PerturbedBox_2D/makeIC.py HydroTests/PerturbedBox_2D/perturbedPlane.yml \
	     HydroTests/PerturbedBox_3D/makeIC.py HydroTests/PerturbedBox_3D/perturbedBox.yml HydroTests/PerturbedBox_3D/run.sh \
	     HydroTests/SedovBlast_1D/makeIC.py HydroTests/SedovBlast_1D/plotSolution.py HydroTests/SedovBlast_1D/run.sh HydroTests/SedovBlast_1D/sedov.yml \
	     HydroTests/SedovBlast_2D/getGlass.sh HydroTests/SedovBlast_2D/makeIC.py HydroTests/SedovBlast_2D/plotSolution.py HydroTests/SedovBlast_2D/run.sh HydroTests/SedovBlast_2D/sedov.yml \
	     HydroTests/SedovBlast_3D/getGlass.sh HydroTests/SedovBlast_3D/makeIC.py HydroTests/SedovBlast_3D/plotSolution.py HydroTests/SedovBlast_3D/run.sh HydroTests/SedovBlast_3D/sedov.yml \
             HydroTests/SineWavePotential_1D/makeIC.py HydroTests/SineWavePotential_1D/plotSolution.py HydroTests/SineWavePotential_1D/run.sh HydroTests/SineWavePotential_1D/sineWavePotential.yml \
             HydroTests/SineWavePotential_2D/makeIC.py HydroTests/SineWavePotential_2D/plotSolution.py HydroTests/SineWavePotential_2D/run.sh HydroTests/SineWavePotential_2D/sineWavePotential.yml \
             HydroTests/SineWavePotential_3D/makeIC.py HydroTests/SineWavePotential_3D/plotSolution.py HydroTests/SineWavePotential_3D/run.sh HydroTests/SineWavePotential_3D/sineWavePotential.yml \
	     HydroTests/SodShock_1D/makeIC.py HydroTests/SodShock_1D/plotSolution.py HydroTests/SodShock_1D/run.sh HydroTests/SodShock_1D/sodShock.yml \
	     HydroTests/SodShock_2D/getGlass.sh HydroTests/SodShock_2D/makeIC.py HydroTests/SodShock_2D/plotSolution.py HydroTests/SodShock_2D/run.sh HydroTests/SodShock_2D/sodShock.yml \
	     HydroTests/SodShock_3D/getGlass.sh HydroTests/SodShock_3D/makeIC.py HydroTests/SodShock_3D/plotSolution.py HydroTests/SodShock_3D/run.sh HydroTests/SodShock_3D/sodShock.yml \
	     HydroTests/SquareTest_2D/makeIC.py HydroTests/SquareTest_2D/plotSolution.py HydroTests/SquareTest_2D/run.sh HydroTests/SquareTest_2D/square.yml \
	     HydroTests/UniformBox_2D/makeIC.py HydroTests/UniformBox_2D/run.sh HydroTests/UniformBox_2D/uniformPlane.yml \
	     HydroTests/UniformBox_3D/makeICbig.py HydroTests/UniformBox_3D/makeIC.py HydroTests/UniformBox_3D/run.sh HydroTests/UniformBox_3D/uniformBox.yml \
             SmallCosmoVolume/SmallCosmoVolume_hydro/README SmallCosmoVolume/SmallCosmoVolume_hydro/getIC.sh SmallCosmoVolume/SmallCosmoVolume_hydro/run.sh SmallCosmoVolume/SmallCosmoVolume_hydro/small_cosmo_volume.yml SmallCosmoVolume/SmallCosmoVolume_hydro/plotTempEvolution.py \
             SmallCosmoVolume/SmallCosmoVolume_DM/README SmallCosmoVolume/SmallCosmoVolume_DM/getIC.sh SmallCosmoVolume/SmallCosmoVolume_DM/run.sh SmallCosmoVolume/SmallCosmoVolume_DM/small_cosmo_volume_dm.yml SmallCosmoVolume/SmallCosmoVolume_DM/stf_input_6dfof_dmonly_sub.cfg \
	     PMillennium/PMillennium-384/p-mill-384.yml \
             PMillennium/PMillennium-768/p-mill-768.yml \
             PMillennium/PMillennium-1536/p-mill-1536.yml \
	     SantaBarbara/SantaBarbara-256/README SantaBarbara/SantaBarbara-256/getIC.sh SantaBarbara/SantaBarbara-256/santa_barbara.yml SantaBarbara/SantaBarbara-256/run.sh \
	     SantaBarbara/SantaBarbara-128/README SantaBarbara/SantaBarbara-128/getIC.sh SantaBarbara/SantaBarbara-128/santa_barbara.yml SantaBarbara/SantaBarbara-128/run.sh

# Default parameter file
EXTRA_DIST += parameter_example.yml
