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

# Add the non-standard paths to the included library headers
AM_CFLAGS = $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS) $(GRACKLE_INCS)

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS)

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Additional dependencies for shared libraries.
EXTRA_LIBS = $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS) $(GSL_LIBS)

# MPI libraries.
MPI_LIBS = $(PARMETIS_LIBS) $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(PARMETIS_INCS) $(METIS_INCS)

# Build the libswiftsim library and a convenience library just for the gravity tasks
lib_LTLIBRARIES = libswiftsim.la
noinst_LTLIBRARIES = libgrav.la 
# Build a MPI-enabled version too?
if HAVEMPI
lib_LTLIBRARIES += libswiftsim_mpi.la
noinst_LTLIBRARIES += libgrav_mpi.la
endif

# List required headers
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h \
    engine.h swift.h serial_io.h timers.h debug.h scheduler.h proxy.h parallel_io.h \
    common_io.h single_io.h multipole.h map.h tools.h partition.h partition_fixed_costs.h \
    clocks.h parser.h physical_constants.h physical_constants_cgs.h potential.h version.h \
    hydro_properties.h riemann.h threadpool.h cooling_io.h cooling.h cooling_struct.h \
    statistics.h memswap.h cache.h runner_doiact_hydro_vec.h profiler.h entropy_floor.h \
    dump.h logger.h active.h timeline.h xmf.h gravity_properties.h gravity_derivatives.h \
    gravity_softened_derivatives.h vector_power.h collectgroup.h hydro_space.h sort_part.h \
    chemistry.h chemistry_io.h chemistry_struct.h cosmology.h restart.h space_getsid.h utilities.h \
    mesh_gravity.h cbrt.h exp10.h velociraptor_interface.h swift_velociraptor_part.h outputlist.h \
    logger_io.h tracers_io.h tracers.h tracers_struct.h star_formation_io.h fof.h fof_struct.h fof_io.h \
    star_formation_struct.h star_formation.h star_formation_iact.h \
    star_formation_logger.h star_formation_logger_struct.h \
    pressure_floor.h pressure_floor_struct.h pressure_floor_iact.h \
    velociraptor_struct.h velociraptor_io.h random.h memuse.h mpiuse.h memuse_rnodes.h \
    black_holes.h black_holes_io.h black_holes_properties.h black_holes_struct.h \
    feedback.h feedback_struct.h feedback_properties.h task_order.h

# source files for EAGLE cooling
EAGLE_COOLING_SOURCES =
if HAVEEAGLECOOLING
EAGLE_COOLING_SOURCES += cooling/EAGLE/cooling.c cooling/EAGLE/cooling_tables.c
endif

# source files for EAGLE feedback
EAGLE_FEEDBACK_SOURCES =
if HAVEEAGLEFEEDBACK
EAGLE_FEEDBACK_SOURCES += feedback/EAGLE/feedback.c
endif

# Common source files
AM_SOURCES = space.c runner_main.c runner_doiact_hydro.c \
    runner_doiact_stars.c runner_doiact_black_holes.c runner_ghost.c runner_recv.c \
    runner_sort.c runner_drift.c runner_black_holes.c runner_time_integration.c \
    runner_doiact_hydro_vec.c runner_others.c\
    queue.c task.c cell.c engine.c engine_maketasks.c \
    engine_marktasks.c engine_drift.c engine_unskip.c engine_collect_end_of_step.c \
    engine_redistribute.c engine_fof.c serial_io.c timers.c debug.c scheduler.c \
    proxy.c parallel_io.c units.c common_io.c single_io.c multipole.c version.c map.c \
    kernel_hydro.c tools.c part.c partition.c clocks.c parser.c \
    physical_constants.c potential.c hydro_properties.c \
    threadpool.c cooling.c star_formation.c \
    statistics.c profiler.c dump.c logger.c \
    part_type.c xmf.c gravity_properties.c gravity.c \
    collectgroup.c hydro_space.c equation_of_state.c \
    chemistry.c cosmology.c restart.c mesh_gravity.c velociraptor_interface.c \
    outputlist.c velociraptor_dummy.c logger_io.c memuse.c mpiuse.c memuse_rnodes.c fof.c \
    hashmap.c pressure_floor.c \
    $(EAGLE_COOLING_SOURCES) $(EAGLE_FEEDBACK_SOURCES)

# Include files for distribution, not installation.
nobase_noinst_HEADERS = align.h approx_math.h atomic.h barrier.h cycle.h error.h inline.h kernel_hydro.h kernel_gravity.h \
		 gravity_iact.h kernel_long_gravity.h vector.h cache.h \
	         runner_doiact_nosort.h runner_doiact_hydro.h runner_doiact_stars.h runner_doiact_black_holes.h runner_doiact_grav.h \
                 runner_doiact_functions_hydro.h runner_doiact_functions_stars.h runner_doiact_functions_black_holes.h \
		 units.h intrinsics.h minmax.h \
                 kick.h timestep.h drift.h adiabatic_index.h io_properties.h dimension.h part_type.h periodic.h memswap.h \
                 dump.h logger.h sign.h logger_io.h timestep_limiter.h hashmap.h \
		 gravity.h gravity_io.h gravity_cache.h \
		 gravity/Default/gravity.h gravity/Default/gravity_iact.h gravity/Default/gravity_io.h \
		 gravity/Default/gravity_debug.h gravity/Default/gravity_part.h  \
		 gravity/Potential/gravity.h gravity/Potential/gravity_iact.h gravity/Potential/gravity_io.h \
		 gravity/Potential/gravity_debug.h gravity/Potential/gravity_part.h  \
		 gravity/MultiSoftening/gravity.h gravity/MultiSoftening/gravity_iact.h gravity/MultiSoftening/gravity_io.h \
		 gravity/MultiSoftening/gravity_debug.h gravity/MultiSoftening/gravity_part.h  \
		 equation_of_state.h \
		 equation_of_state/ideal_gas/equation_of_state.h equation_of_state/isothermal/equation_of_state.h \
	 	 hydro.h hydro_io.h hydro_parameters.h \
		 hydro/Minimal/hydro.h hydro/Minimal/hydro_iact.h hydro/Minimal/hydro_io.h \
                 hydro/Minimal/hydro_debug.h hydro/Minimal/hydro_part.h \
		 hydro/Minimal/hydro_parameters.h \
		 hydro/Default/hydro.h hydro/Default/hydro_iact.h hydro/Default/hydro_io.h \
                 hydro/Default/hydro_debug.h hydro/Default/hydro_part.h \
		 hydro/Default/hydro_parameters.h \
		 hydro/Gadget2/hydro.h hydro/Gadget2/hydro_iact.h hydro/Gadget2/hydro_io.h \
                 hydro/Gadget2/hydro_debug.h hydro/Gadget2/hydro_part.h \
		 hydro/Gadget2/hydro_parameters.h \
		 hydro/PressureEntropy/hydro.h hydro/PressureEntropy/hydro_iact.h hydro/PressureEntropy/hydro_io.h \
                 hydro/PressureEntropy/hydro_debug.h hydro/PressureEntropy/hydro_part.h \
		 hydro/PressureEntropy/hydro_parameters.h \
		 hydro/PressureEnergy/hydro.h hydro/PressureEnergy/hydro_iact.h hydro/PressureEnergy/hydro_io.h \
                 hydro/PressureEnergy/hydro_debug.h hydro/PressureEnergy/hydro_part.h \
		 hydro/PressureEnergy/hydro_parameters.h \
		 hydro/PressureEnergyMorrisMonaghanAV/hydro.h hydro/PressureEnergyMorrisMonaghanAV/hydro_iact.h hydro/PressureEnergyMorrisMonaghanAV/hydro_io.h \
                 hydro/PressureEnergyMorrisMonaghanAV/hydro_debug.h hydro/PressureEnergyMorrisMonaghanAV/hydro_part.h \
		 hydro/PressureEnergyMorrisMonaghanAV/hydro_parameters.h \
		 hydro/AnarchyPU/hydro.h hydro/PressureEnergy/hydro_iact.h hydro/PressureEnergy/hydro_io.h \
                 hydro/AnarchyPU/hydro_debug.h hydro/PressureEnergy/hydro_part.h \
		 hydro/AnarchyPU/hydro_parameters.h \
		 hydro/AnarchyDU/hydro.h hydro/PressureEnergy/hydro_iact.h hydro/PressureEnergy/hydro_io.h \
                 hydro/AnarchyDU/hydro_debug.h hydro/PressureEnergy/hydro_part.h \
		 hydro/AnarchyDU/hydro_parameters.h \
		 hydro/GizmoMFV/hydro.h hydro/GizmoMFV/hydro_iact.h \
                 hydro/GizmoMFV/hydro_io.h hydro/GizmoMFV/hydro_debug.h \
                 hydro/GizmoMFV/hydro_part.h \
                 hydro/GizmoMFV/hydro_gradients_gizmo.h \
                 hydro/GizmoMFV/hydro_gradients.h \
                 hydro/GizmoMFV/hydro_gradients_sph.h \
                 hydro/GizmoMFV/hydro_slope_limiters_cell.h \
                 hydro/GizmoMFV/hydro_slope_limiters_face.h \
                 hydro/GizmoMFV/hydro_slope_limiters.h \
                 hydro/GizmoMFV/hydro_unphysical.h \
                 hydro/GizmoMFV/hydro_velocities.h \
		 hydro/GizmoMFV/hydro_parameters.h \
		 hydro/GizmoMFM/hydro.h hydro/GizmoMFM/hydro_iact.h \
                 hydro/GizmoMFM/hydro_io.h hydro/GizmoMFM/hydro_debug.h \
                 hydro/GizmoMFM/hydro_part.h \
                 hydro/GizmoMFM/hydro_gradients_gizmo.h \
                 hydro/GizmoMFM/hydro_gradients.h \
                 hydro/GizmoMFM/hydro_gradients_sph.h \
                 hydro/GizmoMFM/hydro_slope_limiters_cell.h \
                 hydro/GizmoMFM/hydro_slope_limiters_face.h \
                 hydro/GizmoMFM/hydro_slope_limiters.h \
                 hydro/GizmoMFM/hydro_unphysical.h \
		 hydro/GizmoMFM/hydro_parameters.h \
                 hydro/Shadowswift/hydro_debug.h \
                 hydro/Shadowswift/hydro_gradients.h hydro/Shadowswift/hydro.h \
                 hydro/Shadowswift/hydro_iact.h \
                 hydro/Shadowswift/hydro_io.h \
                 hydro/Shadowswift/hydro_part.h \
                 hydro/Shadowswift/hydro_slope_limiters_cell.h \
                 hydro/Shadowswift/hydro_slope_limiters_face.h \
                 hydro/Shadowswift/hydro_slope_limiters.h \
                 hydro/Shadowswift/voronoi1d_algorithm.h \
                 hydro/Shadowswift/voronoi1d_cell.h \
                 hydro/Shadowswift/voronoi2d_algorithm.h \
                 hydro/Shadowswift/voronoi2d_cell.h \
                 hydro/Shadowswift/voronoi3d_algorithm.h \
                 hydro/Shadowswift/voronoi3d_cell.h \
                 hydro/Shadowswift/voronoi_algorithm.h \
                 hydro/Shadowswift/voronoi_cell.h \
                 hydro/Shadowswift/hydro_parameters.h \
	         riemann/riemann_hllc.h riemann/riemann_trrs.h \
		 riemann/riemann_exact.h riemann/riemann_vacuum.h \
                 riemann/riemann_checks.h \
	 	 stars.h stars_io.h \
		 stars/Default/stars.h stars/Default/stars_iact.h stars/Default/stars_io.h \
		 stars/Default/stars_debug.h stars/Default/stars_part.h  \
		 stars/EAGLE/stars.h stars/EAGLE/stars_iact.h stars/EAGLE/stars_io.h \
		 stars/EAGLE/stars_debug.h stars/EAGLE/stars_part.h \
	         potential/none/potential.h potential/point_mass/potential.h \
                 potential/isothermal/potential.h potential/disc_patch/potential.h \
                 potential/sine_wave/potential.h \
		 star_formation/none/star_formation.h star_formation/none/star_formation_struct.h \
		 star_formation/none/star_formation_io.h star_formation/none/star_formation_iact.h \
		 star_formation/EAGLE/star_formation.h star_formation/EAGLE/star_formation_struct.h \
		 star_formation/EAGLE/star_formation_io.h star_formation/EAGLE/star_formation_iact.h \
		 star_formation/GEAR/star_formation.h star_formation/GEAR/star_formation_struct.h \
		 star_formation/GEAR/star_formation_io.h star_formation/GEAR/star_formation_iact.h \
                 star_formation/EAGLE/star_formation_logger.h star_formation/EAGLE/star_formation_logger_struct.h \
                 star_formation/GEAR/star_formation_logger.h star_formation/GEAR/star_formation_logger_struct.h \
                 star_formation/none/star_formation_logger.h star_formation/none/star_formation_logger_struct.h \
                 cooling/none/cooling.h cooling/none/cooling_struct.h \
                 cooling/none/cooling_io.h \
		 cooling/Compton/cooling.h cooling/Compton/cooling_struct.h \
                 cooling/Compton/cooling_io.h \
	         cooling/const_du/cooling.h cooling/const_du/cooling_struct.h \
                 cooling/const_du/cooling_io.h \
                 cooling/const_lambda/cooling.h cooling/const_lambda/cooling_struct.h \
                 cooling/const_lambda/cooling_io.h \
                 cooling/grackle/cooling.h cooling/grackle/cooling_struct.h \
                 cooling/grackle/cooling_io.h \
		 cooling/EAGLE/cooling.h cooling/EAGLE/cooling_struct.h cooling/EAGLE/cooling_tables.h \
                 cooling/EAGLE/cooling_io.h cooling/EAGLE/interpolate.h cooling/EAGLE/cooling_rates.h \
                 chemistry/none/chemistry.h \
		 chemistry/none/chemistry_io.h \
		 chemistry/none/chemistry_struct.h \
		 chemistry/none/chemistry_iact.h \
                 chemistry/GEAR/chemistry.h \
		 chemistry/GEAR/chemistry_io.h \
		 chemistry/GEAR/chemistry_struct.h \
		 chemistry/GEAR/chemistry_iact.h \
                 chemistry/EAGLE/chemistry.h \
		 chemistry/EAGLE/chemistry_io.h \
		 chemistry/EAGLE/chemistry_struct.h\
		 chemistry/EAGLE/chemistry_iact.h \
	         entropy_floor/none/entropy_floor.h \
                 entropy_floor/EAGLE/entropy_floor.h \
		 tracers/none/tracers.h tracers/none/tracers_struct.h \
                 tracers/none/tracers_io.h \
		 tracers/EAGLE/tracers.h tracers/EAGLE/tracers_struct.h \
                 tracers/EAGLE/tracers_io.h \
	         feedback/none/feedback.h feedback/none/feedback_struct.h feedback/none/feedback_iact.h \
                 feedback/none/feedback_properties.h \
	         feedback/EAGLE/feedback.h feedback/EAGLE/feedback_struct.h feedback/EAGLE/feedback_iact.h \
                 feedback/EAGLE/feedback_properties.h feedback/EAGLE/imf.h feedback/EAGLE/interpolate.h \
                 feedback/EAGLE/yield_tables.h \
                 black_holes/Default/black_holes.h black_holes/Default/black_holes_io.h \
		 black_holes/Default/black_holes_part.h black_holes/Default/black_holes_iact.h \
                 black_holes/Default/black_holes_properties.h \
                 black_holes/Default/black_holes_struct.h \
                 black_holes/EAGLE/black_holes.h black_holes/EAGLE/black_holes_io.h \
		 black_holes/EAGLE/black_holes_part.h black_holes/EAGLE/black_holes_iact.h \
                 black_holes/EAGLE/black_holes_properties.h \
                 black_holes/EAGLE/black_holes_struct.h \
		 pressure_floor/GEAR/pressure_floor.h pressure_floor/none/pressure_floor.h \
		 pressure_floor/GEAR/pressure_floor_iact.h pressure_floor/none/pressure_floor_iact.h \
		 pressure_floor/GEAR/pressure_floor_struct.h pressure_floor/none/pressure_floor_struct.h \
	   task_order/GEAR/task_order.h task_order/EAGLE/task_order.h task_order/none/task_order.h


# Sources and special flags for the gravity library
libgrav_la_SOURCES = runner_doiact_grav.c
libgrav_la_CFLAGS = $(AM_CFLAGS) $(GRAVITY_CFLAGS) 
libgrav_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS)

# Sources and special flags for the gravity MPI library
libgrav_mpi_la_SOURCES = runner_doiact_grav.c
libgrav_mpi_la_CFLAGS = $(AM_CFLAGS) $(GRAVITY_CFLAGS) -DWITH_MPI
libgrav_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS)

# Sources and flags for regular library
libswiftsim_la_SOURCES = $(AM_SOURCES)
libswiftsim_la_CFLAGS = $(AM_CFLAGS)
libswiftsim_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS) -version-info 0:0:0
libswiftsim_la_LIBADD = $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS) libgrav.la

# Sources and flags for MPI library
libswiftsim_mpi_la_SOURCES = $(AM_SOURCES)
libswiftsim_mpi_la_CFLAGS = $(AM_CFLAGS) $(MPI_FLAGS)
libswiftsim_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS) -version-info 0:0:0
libswiftsim_mpi_la_SHORTNAME = mpi
libswiftsim_mpi_la_LIBADD = $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS) $(MPI_LIBS) libgrav_mpi.la

# Versioning. If any sources change then update the version_string.h file with
# the current git revision and package version.
# May have a checkout without a version_string.h file and no git command (tar/zip
# download), allow that, but make sure we know it.
version_string.h: version_string.h.in Makefile $(AM_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	if test "X$(GIT_CMD)" != "X"; then \
	    GIT_REVISION=`$(GIT_CMD) describe --abbrev=8  --always --tags --dirty 2>/dev/null`; \
	    GIT_BRANCH=`$(GIT_CMD) branch 2>/dev/null | sed -n 's/^\* \(.*\)/\1/p' 2>/dev/null`; \
            GIT_DATE=`$(GIT_CMD) log -1 --format=%ci 2>/dev/null`; \
	    sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	        -e "s,@GIT_REVISION\@,$${GIT_REVISION}," \
	        -e "s|@GIT_BRANCH\@|$${GIT_BRANCH}|" \
	        -e "s|@GIT_DATE\@|$${GIT_DATE}|" \
	        -e "s|@SWIFT_CFLAGS\@|$(CFLAGS)|" $< > version_string.h; \
	else \
	    if test ! -f version_string.h; then \
	        sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	            -e "s,@GIT_REVISION\@,unknown," \
		    -e "s,@GIT_BRANCH\@,unknown," \
		    -e "s,@GIT_DATE\@,unknown," \
	            -e "s|@SWIFT_CFLAGS\@|$(CFLAGS)|" $< > version_string.h; \
	    fi; \
	fi

#  Make sure version_string.h is built first.
BUILT_SOURCES = version_string.h

#  And distribute the built files.
EXTRA_DIST = version_string.h version_string.h.in
