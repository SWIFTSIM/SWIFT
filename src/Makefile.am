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
AM_CFLAGS = $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS)

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS) -version-info 0:0:0

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Additional dependencies for shared libraries.
EXTRA_LIBS = $(HDF5_LIBS) $(FFTW_LIBS) $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIB) $(GSL_LIBS)

# MPI libraries.
MPI_LIBS = $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(METIS_INCS)

# Build the libswiftsim library
lib_LTLIBRARIES = libswiftsim.la
# Build a MPI-enabled version too?
if HAVEMPI
lib_LTLIBRARIES += libswiftsim_mpi.la
endif

# List required headers
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h \
    engine.h swift.h serial_io.h timers.h debug.h scheduler.h proxy.h parallel_io.h \
    common_io.h single_io.h multipole.h map.h tools.h partition.h clocks.h parser.h \
    physical_constants.h physical_constants_cgs.h potential.h version.h \
    hydro_properties.h riemann.h threadpool.h cooling_io.h cooling.h cooling_struct.h \
    sourceterms.h sourceterms_struct.h statistics.h memswap.h cache.h runner_doiact_vec.h profiler.h \
    dump.h logger.h active.h timeline.h xmf.h gravity_properties.h gravity_derivatives.h \
    gravity_softened_derivatives.h vector_power.h collectgroup.h hydro_space.h sort_part.h \
    chemistry.h chemistry_io.h chemistry_struct.h cosmology.h restart.h space_getsid.h utilities.h \
    mesh_gravity.h cbrt.h velociraptor_interface.h swift_velociraptor_part.h outputlist.h

# Common source files
AM_SOURCES = space.c runner.c queue.c task.c cell.c engine.c \
    serial_io.c timers.c debug.c scheduler.c proxy.c parallel_io.c \
    units.c common_io.c single_io.c multipole.c version.c map.c \
    kernel_hydro.c tools.c part.c partition.c clocks.c parser.c \
    physical_constants.c potential.c hydro_properties.c \
    threadpool.c cooling.c sourceterms.c \
    statistics.c runner_doiact_vec.c profiler.c dump.c logger.c \
    part_type.c xmf.c gravity_properties.c gravity.c \
    collectgroup.c hydro_space.c equation_of_state.c \
    chemistry.c cosmology.c restart.c mesh_gravity.c velociraptor_interface.c \
    outputlist.c

# Include files for distribution, not installation.
nobase_noinst_HEADERS = align.h approx_math.h atomic.h barrier.h cycle.h error.h inline.h kernel_hydro.h kernel_gravity.h \
		 gravity_iact.h kernel_long_gravity.h vector.h cache.h runner_doiact.h runner_doiact_vec.h runner_doiact_grav.h  \
                 runner_doiact_nosort.h units.h intrinsics.h minmax.h kick.h timestep.h drift.h adiabatic_index.h io_properties.h \
		 dimension.h part_type.h periodic.h memswap.h dump.h logger.h sign.h \
		 gravity.h gravity_io.h gravity_cache.h \
		 gravity/Default/gravity.h gravity/Default/gravity_iact.h gravity/Default/gravity_io.h \
		 gravity/Default/gravity_debug.h gravity/Default/gravity_part.h  \
		 gravity/Potential/gravity.h gravity/Potential/gravity_iact.h gravity/Potential/gravity_io.h \
		 gravity/Potential/gravity_debug.h gravity/Potential/gravity_part.h  \
		 sourceterms.h \
		 equation_of_state.h \
		 equation_of_state/ideal_gas/equation_of_state.h equation_of_state/isothermal/equation_of_state.h \
	 	 hydro.h hydro_io.h \
		 hydro/Minimal/hydro.h hydro/Minimal/hydro_iact.h hydro/Minimal/hydro_io.h \
                 hydro/Minimal/hydro_debug.h hydro/Minimal/hydro_part.h \
		 hydro/Default/hydro.h hydro/Default/hydro_iact.h hydro/Default/hydro_io.h \
                 hydro/Default/hydro_debug.h hydro/Default/hydro_part.h \
		 hydro/Gadget2/hydro.h hydro/Gadget2/hydro_iact.h hydro/Gadget2/hydro_io.h \
                 hydro/Gadget2/hydro_debug.h hydro/Gadget2/hydro_part.h \
		 hydro/PressureEntropy/hydro.h hydro/PressureEntropy/hydro_iact.h hydro/PressureEntropy/hydro_io.h \
                 hydro/PressureEntropy/hydro_debug.h hydro/PressureEntropy/hydro_part.h \
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
                 hydro/GizmoMFM/hydro_velocities.h \
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
	         riemann/riemann_hllc.h riemann/riemann_trrs.h \
		 riemann/riemann_exact.h riemann/riemann_vacuum.h \
                 riemann/riemann_checks.h \
	 	 stars.h stars_io.h \
		 stars/Default/star.h stars/Default/star_iact.h stars/Default/star_io.h \
		 stars/Default/star_debug.h stars/Default/star_part.h  \
	         potential/none/potential.h potential/point_mass/potential.h \
                 potential/isothermal/potential.h potential/disc_patch/potential.h \
                 potential/sine_wave/potential.h \
		 cooling/none/cooling.h cooling/none/cooling_struct.h \
                 cooling/none/cooling_io.h \
	         cooling/const_du/cooling.h cooling/const_du/cooling_struct.h \
                 cooling/const_du/cooling_io.h \
                 cooling/const_lambda/cooling.h cooling/const_lambda/cooling_struct.h \
                 cooling/const_lambda/cooling_io.h \
                 cooling/grackle/cooling.h cooling/grackle/cooling_struct.h \
                 cooling/grackle/cooling_io.h \
		 cooling/EAGLE/cooling.h cooling/EAGLE/cooling_struct.h \
                 cooling/EAGLE/cooling_io.h \
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
		 chemistry/EAGLE/chemistry_iact.h


# Sources and flags for regular library
libswiftsim_la_SOURCES = $(AM_SOURCES)
libswiftsim_la_CFLAGS = $(AM_CFLAGS)
libswiftsim_la_LDFLAGS = $(AM_LDFLAGS) $(EXTRA_LIBS)
libswiftsim_la_LIBADD = $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS)

# Sources and flags for MPI library
libswiftsim_mpi_la_SOURCES = $(AM_SOURCES)
libswiftsim_mpi_la_CFLAGS = $(AM_CFLAGS) $(MPI_FLAGS)
libswiftsim_mpi_la_LDFLAGS = $(AM_LDFLAGS) $(MPI_LIBS) $(EXTRA_LIBS)
libswiftsim_mpi_la_SHORTNAME = mpi
libswiftsim_mpi_la_LIBADD = $(GRACKLE_LIBS) $(VELOCIRAPTOR_LIBS)


# Versioning. If any sources change then update the version_string.h file with
# the current git revision and package version.
# May have a checkout without a version_string.h file and no git command (tar/zip
# download), allow that, but make sure we know it.
version_string.h: version_string.h.in $(AM_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	if test "X$(GIT_CMD)" != "X"; then \
	    GIT_REVISION=`$(GIT_CMD) describe --abbrev=8  --always --tags --dirty`; \
	    GIT_BRANCH=`$(GIT_CMD) branch | sed -n 's/^\* \(.*\)/\1/p'`; \
            GIT_DATE=`$(GIT_CMD) log -1 --format=%ci`; \
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
