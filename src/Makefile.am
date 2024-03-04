# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
AM_CFLAGS = $(HDF5_CPPFLAGS) $(GSL_INCS) $(FFTW_INCS) $(NUMA_INCS) $(GRACKLE_INCS)  $(SUNDIALS_INCS) $(CHEALPIX_CFLAGS)

# Assign a "safe" version number
AM_LDFLAGS = $(HDF5_LDFLAGS) $(FFTW_LIBS)

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Additional dependencies for shared libraries.
EXTRA_LIBS = $(GSL_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) $(NUMA_LIBS) $(PROFILER_LIBS) $(TCMALLOC_LIBS) $(JEMALLOC_LIBS) $(TBBMALLOC_LIBS) $(GRACKLE_LIBS)  $(SUNDIALS_LIBS) $(CHEALPIX_LIBS)

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
include_HEADERS = space.h runner.h queue.h task.h lock.h cell.h part.h const.h 
include_HEADERS += cell_hydro.h cell_stars.h cell_grav.h cell_sinks.h cell_black_holes.h cell_rt.h
include_HEADERS += engine.h swift.h serial_io.h timers.h debug.h scheduler.h proxy.h parallel_io.h 
include_HEADERS += common_io.h single_io.h distributed_io.h map.h tools.h  partition_fixed_costs.h 
include_HEADERS += partition.h clocks.h parser.h physical_constants.h physical_constants_cgs.h potential.h version.h 
include_HEADERS += hydro_properties.h riemann.h threadpool.h cooling_io.h cooling.h cooling_struct.h cooling_properties.h cooling_debug.h
include_HEADERS += statistics.h memswap.h cache.h runner_doiact_hydro_vec.h runner_doiact_undef.h profiler.h entropy_floor.h 
include_HEADERS += csds.h active.h timeline.h xmf.h gravity_properties.h gravity_derivatives.h 
include_HEADERS += gravity_softened_derivatives.h vector_power.h collectgroup.h hydro_space.h sort_part.h 
include_HEADERS += chemistry.h chemistry_additions.h chemistry_io.h chemistry_struct.h chemistry_debug.h
include_HEADERS += cosmology.h restart.h space_getsid.h utilities.h
include_HEADERS += cbrt.h exp10.h velociraptor_interface.h swift_velociraptor_part.h output_list.h 
include_HEADERS += csds_io.h
include_HEADERS += tracers_io.h tracers.h tracers_triggers.h tracers_struct.h tracers_debug.h
include_HEADERS += star_formation_io.h star_formation_debug.h extra_io.h
include_HEADERS += fof.h fof_struct.h fof_io.h fof_catalogue_io.h
include_HEADERS += multipole.h multipole_accept.h multipole_struct.h binomial.h integer_power.h sincos.h 
include_HEADERS += star_formation_struct.h star_formation.h star_formation_iact.h 
include_HEADERS += star_formation_logger.h star_formation_logger_struct.h 
include_HEADERS += pressure_floor.h pressure_floor_struct.h pressure_floor_iact.h pressure_floor_debug.h
include_HEADERS += velociraptor_struct.h velociraptor_io.h random.h memuse.h mpiuse.h memuse_rnodes.h 
include_HEADERS += black_holes.h black_holes_iact.h black_holes_io.h black_holes_properties.h black_holes_struct.h black_holes_debug.h
include_HEADERS += feedback.h feedback_new_stars.h feedback_struct.h feedback_properties.h feedback_debug.h feedback_iact.h
include_HEADERS += space_unique_id.h line_of_sight.h io_compression.h
include_HEADERS += rays.h rays_struct.h
include_HEADERS += sink.h sink_struct.h sink_io.h sink_properties.h sink_debug.h
include_HEADERS += particle_splitting.h particle_splitting_struct.h
include_HEADERS += chemistry_csds.h star_formation_csds.h
include_HEADERS += mesh_gravity.h mesh_gravity_mpi.h mesh_gravity_patch.h mesh_gravity_sort.h row_major_id.h
include_HEADERS += hdf5_object_to_blob.h ic_info.h particle_buffer.h exchange_structs.h
include_HEADERS += lightcone/lightcone.h lightcone/lightcone_particle_io.h lightcone/lightcone_replications.h
include_HEADERS += lightcone/lightcone_crossing.h lightcone/lightcone_array.h lightcone/lightcone_map.h
include_HEADERS += lightcone/lightcone_map_types.h lightcone/projected_kernel.h lightcone/lightcone_shell.h
include_HEADERS += lightcone/healpix_util.h lightcone/pixel_index.h
include_HEADERS += adaptive_softening.h adaptive_softening_iact.h adaptive_softening_struct.h
include_HEADERS += forcing.h
include_HEADERS += power_spectrum.h
include_HEADERS += ghost_stats.h

# source files for EAGLE extra I/O
EAGLE_EXTRA_IO_SOURCES=
if HAVEEAGLEEXTRAIO
EAGLE_EXTRA_IO_SOURCES += extra_io/EAGLE/extra_lightcone_map_types.c
endif

# source files for QLA (Ploeckinger+20) cooling
QLA_COOLING_SOURCES =
if HAVEQLACOOLING
QLA_COOLING_SOURCES += cooling/QLA/cooling.c cooling/QLA/cooling_tables.c
endif

# source files for QLA (Wiersma+09) cooling
QLA_EAGLE_COOLING_SOURCES =
if HAVEQLAEAGLECOOLING
QLA_COOLING_SOURCES += cooling/QLA_EAGLE/cooling.c cooling/QLA_EAGLE/cooling_tables.c
endif

# source files for EAGLE cooling
EAGLE_COOLING_SOURCES =
if HAVEEAGLECOOLING
EAGLE_COOLING_SOURCES += cooling/EAGLE/cooling.c cooling/EAGLE/cooling_tables.c
endif

# source files for PS2020 cooling
PS2020_COOLING_SOURCES =
if HAVEPS2020COOLING
PS2020_COOLING_SOURCES += cooling/PS2020/cooling.c cooling/PS2020/cooling_tables.c
endif

# source files for EAGLE feedback
EAGLE_FEEDBACK_SOURCES =
if HAVEEAGLETHERMALFEEDBACK
EAGLE_FEEDBACK_SOURCES += feedback/EAGLE_thermal/feedback.c
endif
if HAVEEAGLEKINETICFEEDBACK
EAGLE_FEEDBACK_SOURCES += feedback/EAGLE_kinetic/feedback.c
endif

# source files for GRACKLE cooling
GRACKLE_COOLING_SOURCES =
if HAVEGRACKLECOOLING
GRACKLE_COOLING_SOURCES += cooling/grackle/cooling.c
endif

# source files for GEAR feedback
GEAR_FEEDBACK_SOURCES =
if HAVEGEARFEEDBACK
GEAR_FEEDBACK_SOURCES += feedback/GEAR/stellar_evolution.c feedback/GEAR/feedback.c \
	feedback/GEAR/initial_mass_function.c feedback/GEAR/supernovae_ia.c feedback/GEAR/supernovae_ii.c
endif

# source files for AGORA feedback
AGORA_FEEDBACK_SOURCES =
if HAVEAGORAFEEDBACK
AGORA_FEEDBACK_SOURCES += feedback/AGORA/feedback.c
endif

# source files for GEAR RT
GEAR_RT_SOURCES =
if HAVEGEARRT
GEAR_RT_SOURCES += rt/GEAR/rt_interaction_cross_sections.c
endif

# source files for SPHM1RT cooling
SPHM1RT_RT_SOURCES =
if HAVESPHM1RTRT
SPHM1RT_RT_SOURCES += rt/SPHM1RT/rt_rate_equations.c
SPHM1RT_RT_SOURCES += rt/SPHM1RT/rt_cooling.c
endif

# Common source files
AM_SOURCES = space.c space_rebuild.c space_regrid.c space_unique_id.c 
AM_SOURCES += space_sort.c space_split.c space_extras.c space_first_init.c space_init.c 
AM_SOURCES += space_cell_index.c space_recycle.c 
AM_SOURCES += runner_main.c runner_doiact_hydro.c runner_doiact_limiter.c 
AM_SOURCES += runner_doiact_stars.c runner_doiact_black_holes.c runner_ghost.c
AM_SOURCES += runner_recv.c runner_pack.c
AM_SOURCES += runner_sort.c runner_drift.c runner_black_holes.c runner_time_integration.c 
AM_SOURCES += runner_doiact_hydro_vec.c runner_others.c
AM_SOURCES += runner_sinks.c
AM_SOURCES += cell.c cell_convert_part.c cell_drift.c cell_lock.c cell_pack.c cell_split.c 
AM_SOURCES += cell_unskip.c 
AM_SOURCES += engine.c engine_maketasks.c engine_split_particles.c engine_strays.c 
AM_SOURCES += engine_marktasks.c engine_drift.c engine_unskip.c engine_collect_end_of_step.c 
AM_SOURCES += engine_redistribute.c engine_fof.c engine_proxy.c engine_io.c engine_config.c 
AM_SOURCES += queue.c task.c timers.c debug.c scheduler.c proxy.c version.c 
AM_SOURCES += common_io.c common_io_copy.c common_io_cells.c common_io_fields.c 
AM_SOURCES += single_io.c serial_io.c distributed_io.c parallel_io.c 
AM_SOURCES += output_options.c line_of_sight.c restart.c parser.c xmf.c 
AM_SOURCES += kernel_hydro.c tools.c map.c part.c partition.c clocks.c  
AM_SOURCES += physical_constants.c units.c potential.c hydro_properties.c 
AM_SOURCES += threadpool.c cooling.c star_formation.c 
AM_SOURCES += hydro.c stars.c
AM_SOURCES += statistics.c profiler.c csds.c part_type.c 
AM_SOURCES += gravity_properties.c gravity.c multipole.c 
AM_SOURCES += collectgroup.c hydro_space.c equation_of_state.c io_compression.c 
AM_SOURCES += chemistry.c cosmology.c velociraptor_interface.c 
AM_SOURCES += output_list.c csds_io.c memuse.c mpiuse.c memuse_rnodes.c
AM_SOURCES += fof.c fof_catalogue_io.c
AM_SOURCES += hashmap.c
AM_SOURCES += mesh_gravity.c mesh_gravity_mpi.c mesh_gravity_patch.c mesh_gravity_sort.c
AM_SOURCES += runner_neutrino.c
AM_SOURCES += neutrino/Default/fermi_dirac.c neutrino/Default/neutrino.c neutrino/Default/neutrino_response.c 
AM_SOURCES += rt_parameters.c hdf5_object_to_blob.c ic_info.c exchange_structs.c particle_buffer.c
AM_SOURCES += lightcone/lightcone.c lightcone/lightcone_particle_io.c lightcone/lightcone_replications.c
AM_SOURCES += lightcone/healpix_util.c lightcone/lightcone_array.c lightcone/lightcone_map.c
AM_SOURCES += lightcone/lightcone_map_types.c lightcone/projected_kernel.c lightcone/lightcone_shell.c
AM_SOURCES += power_spectrum.c
AM_SOURCES += forcing.c
AM_SOURCES += ghost_stats.c
AM_SOURCES += $(EAGLE_EXTRA_IO_SOURCES)
AM_SOURCES += $(QLA_COOLING_SOURCES) $(QLA_EAGLE_COOLING_SOURCES) 
AM_SOURCES += $(EAGLE_COOLING_SOURCES) $(EAGLE_FEEDBACK_SOURCES) 
AM_SOURCES += $(GRACKLE_COOLING_SOURCES) $(GEAR_FEEDBACK_SOURCES) 
AM_SOURCES += $(AGORA_FEEDBACK_SOURCES)
AM_SOURCES += $(PS2020_COOLING_SOURCES)
AM_SOURCES += $(SPHM1RT_RT_SOURCES)
AM_SOURCES += $(GEAR_RT_SOURCES)

# Include files for distribution, not installation.
nobase_noinst_HEADERS = align.h approx_math.h atomic.h barrier.h cycle.h error.h inline.h kernel_hydro.h kernel_gravity.h 
nobase_noinst_HEADERS += gravity_iact.h kernel_long_gravity.h vector.h accumulate.h cache.h exp.h log.h
nobase_noinst_HEADERS += runner_doiact_nosort.h runner_doiact_hydro.h runner_doiact_stars.h runner_doiact_black_holes.h runner_doiact_grav.h 
nobase_noinst_HEADERS += runner_doiact_functions_hydro.h runner_doiact_functions_stars.h runner_doiact_functions_black_holes.h 
nobase_noinst_HEADERS += runner_doiact_functions_limiter.h runner_doiact_limiter.h units.h intrinsics.h minmax.h 
nobase_noinst_HEADERS += runner_doiact_sinks.h
nobase_noinst_HEADERS += kick.h timestep.h drift.h adiabatic_index.h io_properties.h dimension.h part_type.h periodic.h memswap.h 
nobase_noinst_HEADERS += timestep_limiter.h timestep_limiter_iact.h timestep_sync.h timestep_sync_part.h timestep_limiter_struct.h 
nobase_noinst_HEADERS += csds.h sign.h csds_io.h hashmap.h gravity.h gravity_io.h gravity_csds.h  gravity_cache.h output_options.h
nobase_noinst_HEADERS += gravity/Default/gravity.h gravity/Default/gravity_iact.h gravity/Default/gravity_io.h 
nobase_noinst_HEADERS += gravity/Default/gravity_debug.h gravity/Default/gravity_part.h  
nobase_noinst_HEADERS += gravity/MultiSoftening/gravity.h gravity/MultiSoftening/gravity_iact.h gravity/MultiSoftening/gravity_io.h 
nobase_noinst_HEADERS += gravity/MultiSoftening/gravity_debug.h gravity/MultiSoftening/gravity_part.h 
nobase_noinst_HEADERS += gravity/MultiSoftening/gravity_csds.h 
nobase_noinst_HEADERS += equation_of_state.h 
nobase_noinst_HEADERS += equation_of_state/ideal_gas/equation_of_state.h equation_of_state/isothermal/equation_of_state.h equation_of_state/barotropic/equation_of_state.h
nobase_noinst_HEADERS += signal_velocity.h
nobase_noinst_HEADERS += hydro.h hydro_io.h hydro_csds.h hydro_parameters.h 
nobase_noinst_HEADERS += hydro/None/hydro.h hydro/None/hydro_iact.h hydro/None/hydro_io.h 
nobase_noinst_HEADERS += hydro/None/hydro_debug.h hydro/None/hydro_part.h 
nobase_noinst_HEADERS += hydro/None/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/Minimal/hydro.h hydro/Minimal/hydro_iact.h hydro/Minimal/hydro_io.h 
nobase_noinst_HEADERS += hydro/Minimal/hydro_debug.h hydro/Minimal/hydro_part.h 
nobase_noinst_HEADERS += hydro/Minimal/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/Phantom/hydro.h hydro/Phantom/hydro_iact.h hydro/Phantom/hydro_io.h 
nobase_noinst_HEADERS += hydro/Phantom/hydro_debug.h hydro/Phantom/hydro_part.h 
nobase_noinst_HEADERS += hydro/Phantom/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/Gadget2/hydro.h hydro/Gadget2/hydro_iact.h hydro/Gadget2/hydro_io.h 
nobase_noinst_HEADERS += hydro/Gadget2/hydro_debug.h hydro/Gadget2/hydro_part.h 
nobase_noinst_HEADERS += hydro/Gadget2/hydro_parameters.h hydro/Gadget2/hydro_csds.h 
nobase_noinst_HEADERS += hydro/PressureEntropy/hydro.h hydro/PressureEntropy/hydro_iact.h hydro/PressureEntropy/hydro_io.h 
nobase_noinst_HEADERS += hydro/PressureEntropy/hydro_debug.h hydro/PressureEntropy/hydro_part.h 
nobase_noinst_HEADERS += hydro/PressureEntropy/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/PressureEnergy/hydro.h hydro/PressureEnergy/hydro_iact.h hydro/PressureEnergy/hydro_io.h 
nobase_noinst_HEADERS += hydro/PressureEnergy/hydro_debug.h hydro/PressureEnergy/hydro_part.h 
nobase_noinst_HEADERS += hydro/PressureEnergy/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/PressureEnergyMorrisMonaghanAV/hydro.h hydro/PressureEnergyMorrisMonaghanAV/hydro_iact.h hydro/PressureEnergyMorrisMonaghanAV/hydro_io.h 
nobase_noinst_HEADERS += hydro/PressureEnergyMorrisMonaghanAV/hydro_debug.h hydro/PressureEnergyMorrisMonaghanAV/hydro_part.h 
nobase_noinst_HEADERS += hydro/PressureEnergyMorrisMonaghanAV/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/AnarchyPU/hydro.h hydro/AnarchyPU/hydro_iact.h hydro/AnarchyPU/hydro_io.h 
nobase_noinst_HEADERS += hydro/AnarchyPU/hydro_debug.h hydro/AnarchyPU/hydro_part.h 
nobase_noinst_HEADERS += hydro/AnarchyPU/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/SPHENIX/hydro.h hydro/SPHENIX/hydro_iact.h hydro/SPHENIX/hydro_io.h 
nobase_noinst_HEADERS += hydro/SPHENIX/hydro_debug.h hydro/SPHENIX/hydro_part.h hydro/SPHENIX/hydro_csds.h
nobase_noinst_HEADERS += hydro/SPHENIX/hydro_parameters.h
nobase_noinst_HEADERS += hydro/Gasoline/hydro.h hydro/Gasoline/hydro_iact.h hydro/Gasoline/hydro_io.h 
nobase_noinst_HEADERS += hydro/Gasoline/hydro_debug.h hydro/Gasoline/hydro_part.h
nobase_noinst_HEADERS += hydro/Gasoline/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_parameters.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_io.h hydro/Gizmo/hydro_debug.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro.h hydro/Gizmo/hydro_iact.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_part.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_gradients.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_getters.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_setters.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_flux.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_slope_limiters.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_slope_limiters_face.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_slope_limiters_cell.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_unphysical.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_gradients_sph.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_gradients_gizmo.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_velocities.h 
nobase_noinst_HEADERS += hydro/Gizmo/hydro_lloyd.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFV/hydro_debug.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFV/hydro_part.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFV/hydro_velocities.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFV/hydro_flux.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFM/hydro_debug.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFM/hydro_part.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFM/hydro_flux.h 
nobase_noinst_HEADERS += hydro/Gizmo/MFM/hydro_velocities.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_debug.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_gradients.h hydro/Shadowswift/hydro.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_iact.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_io.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_part.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_slope_limiters_cell.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_slope_limiters_face.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_slope_limiters.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi1d_algorithm.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi1d_cell.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi2d_algorithm.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi2d_cell.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi3d_algorithm.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi3d_cell.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi_algorithm.h 
nobase_noinst_HEADERS += hydro/Shadowswift/voronoi_cell.h 
nobase_noinst_HEADERS += hydro/Shadowswift/hydro_parameters.h 
nobase_noinst_HEADERS += mhd.h mhd_struct.h mhd_io.h
nobase_noinst_HEADERS += mhd/None/mhd.h mhd/None/mhd_iact.h mhd/None/mhd_struct.h mhd/None/mhd_io.h mhd/None/mhd_debug.h mhd/None/mhd_parameters.h
nobase_noinst_HEADERS += riemann/riemann_hllc.h riemann/riemann_trrs.h 
nobase_noinst_HEADERS += riemann/riemann_exact.h riemann/riemann_vacuum.h 
nobase_noinst_HEADERS += riemann/riemann_checks.h 
nobase_noinst_HEADERS += rt.h  
nobase_noinst_HEADERS += rt_additions.h  
nobase_noinst_HEADERS += rt_io.h 
nobase_noinst_HEADERS += rt_parameters.h 
nobase_noinst_HEADERS += rt_properties.h 
nobase_noinst_HEADERS += rt_struct.h 
nobase_noinst_HEADERS += rt/none/rt.h  
nobase_noinst_HEADERS += rt/none/rt_additions.h 
nobase_noinst_HEADERS += rt/none/rt_iact.h 
nobase_noinst_HEADERS += rt/none/rt_io.h 
nobase_noinst_HEADERS += rt/none/rt_parameters.h
nobase_noinst_HEADERS += rt/none/rt_properties.h 
nobase_noinst_HEADERS += rt/none/rt_struct.h
nobase_noinst_HEADERS += rt/debug/rt.h  
nobase_noinst_HEADERS += rt/debug/rt_additions.h 
nobase_noinst_HEADERS += rt/debug/rt_debugging.h 
nobase_noinst_HEADERS += rt/debug/rt_gradients.h 
nobase_noinst_HEADERS += rt/debug/rt_iact.h 
nobase_noinst_HEADERS += rt/debug/rt_io.h 
nobase_noinst_HEADERS += rt/debug/rt_parameters.h
nobase_noinst_HEADERS += rt/debug/rt_properties.h 
nobase_noinst_HEADERS += rt/debug/rt_struct.h
nobase_noinst_HEADERS += rt/GEAR/rt_additions.h 
nobase_noinst_HEADERS += rt/GEAR/rt_blackbody.h 
nobase_noinst_HEADERS += rt/GEAR/rt_debugging.h 
nobase_noinst_HEADERS += rt/GEAR/rt_flux.h 
nobase_noinst_HEADERS += rt/GEAR/rt_getters.h 
nobase_noinst_HEADERS += rt/GEAR/rt_grackle_utils.h 
nobase_noinst_HEADERS += rt/GEAR/rt_gradients.h 
nobase_noinst_HEADERS += rt/GEAR/rt.h  
nobase_noinst_HEADERS += rt/GEAR/rt_iact.h 
nobase_noinst_HEADERS += rt/GEAR/rt_interaction_cross_sections.h 
nobase_noinst_HEADERS += rt/GEAR/rt_interaction_rates.h 
nobase_noinst_HEADERS += rt/GEAR/rt_io.h 
nobase_noinst_HEADERS += rt/GEAR/rt_ionization_equilibrium.h 
nobase_noinst_HEADERS += rt/GEAR/rt_parameters.h
nobase_noinst_HEADERS += rt/GEAR/rt_properties.h 
nobase_noinst_HEADERS += rt/GEAR/rt_riemann_GLF.h 
nobase_noinst_HEADERS += rt/GEAR/rt_riemann_HLL_eigenvalues.h 
nobase_noinst_HEADERS += rt/GEAR/rt_riemann_HLL.h 
nobase_noinst_HEADERS += rt/GEAR/rt_slope_limiters_cell.h
nobase_noinst_HEADERS += rt/GEAR/rt_slope_limiters_face.h 
nobase_noinst_HEADERS += rt/GEAR/rt_species.h 
nobase_noinst_HEADERS += rt/GEAR/rt_stellar_emission_model.h 
nobase_noinst_HEADERS += rt/GEAR/rt_struct.h
nobase_noinst_HEADERS += rt/GEAR/rt_thermochemistry.h
nobase_noinst_HEADERS += rt/GEAR/rt_thermochemistry_utils.h
nobase_noinst_HEADERS += rt/GEAR/rt_unphysical.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt.h  
nobase_noinst_HEADERS += rt/SPHM1RT/rt_getters.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_setters.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_iact.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_io.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_parameters.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt_properties.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_struct.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt_gradients.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_additions.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_cooling.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt_cooling_rates.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt_unphysical.h
nobase_noinst_HEADERS += rt/SPHM1RT/rt_species_and_elements.h 
nobase_noinst_HEADERS += rt/SPHM1RT/rt_stellar_emission_rate.h 
nobase_noinst_HEADERS += stars.h stars_io.h stars_csds.h 
nobase_noinst_HEADERS += stars/None/stars.h stars/None/stars_iact.h stars/None/stars_io.h 
nobase_noinst_HEADERS += stars/None/stars_debug.h stars/None/stars_part.h
nobase_noinst_HEADERS += stars/Basic/stars.h stars/Basic/stars_iact.h stars/Basic/stars_io.h 
nobase_noinst_HEADERS += stars/Basic/stars_debug.h stars/Basic/stars_part.h stars/Basic/stars_csds.h  
nobase_noinst_HEADERS += stars/EAGLE/stars.h stars/EAGLE/stars_iact.h stars/EAGLE/stars_io.h 
nobase_noinst_HEADERS += stars/EAGLE/stars_debug.h stars/EAGLE/stars_part.h 
nobase_noinst_HEADERS += stars/GEAR/stars.h stars/GEAR/stars_iact.h stars/GEAR/stars_io.h 
nobase_noinst_HEADERS += stars/GEAR/stars_debug.h stars/GEAR/stars_csds.h stars/GEAR/stars_part.h
nobase_noinst_HEADERS += forcing/none/forcing.h forcing/roberts_flow/forcing.h forcing/roberts_flow_acceleration/forcing.h 
nobase_noinst_HEADERS += potential/none/potential.h potential/point_mass/potential.h 
nobase_noinst_HEADERS += potential/isothermal/potential.h potential/disc_patch/potential.h 
nobase_noinst_HEADERS += potential/sine_wave/potential.h potential/constant/potential.h 
nobase_noinst_HEADERS += potential/hernquist/potential.h potential/nfw/potential.h 
nobase_noinst_HEADERS += potential/nfw_mn/potential.h potential/point_mass_softened/potential.h 
nobase_noinst_HEADERS += star_formation/none/star_formation.h star_formation/none/star_formation_struct.h 
nobase_noinst_HEADERS += star_formation/none/star_formation_io.h star_formation/none/star_formation_iact.h
nobase_noinst_HEADERS += star_formation/none/star_formation_csds.h star_formation/none/star_formation_debug.h
nobase_noinst_HEADERS += star_formation/QLA/star_formation.h star_formation/QLA/star_formation_struct.h 
nobase_noinst_HEADERS += star_formation/QLA/star_formation_io.h star_formation/QLA/star_formation_iact.h 
nobase_noinst_HEADERS += star_formation/QLA/star_formation_debug.h
nobase_noinst_HEADERS += star_formation/EAGLE/star_formation.h star_formation/EAGLE/star_formation_struct.h 
nobase_noinst_HEADERS += star_formation/EAGLE/star_formation_io.h star_formation/EAGLE/star_formation_iact.h 
nobase_noinst_HEADERS += star_formation/EAGLE/star_formation_debug.h
nobase_noinst_HEADERS += star_formation/GEAR/star_formation.h star_formation/GEAR/star_formation_struct.h 
nobase_noinst_HEADERS += star_formation/GEAR/star_formation_io.h star_formation/GEAR/star_formation_iact.h 
nobase_noinst_HEADERS += star_formation/GEAR/star_formation_csds.h star_formation/GEAR/star_formation_debug.h
nobase_noinst_HEADERS += star_formation/EAGLE/star_formation_logger.h star_formation/EAGLE/star_formation_logger_struct.h 
nobase_noinst_HEADERS += star_formation/GEAR/star_formation_logger.h star_formation/GEAR/star_formation_logger_struct.h 
nobase_noinst_HEADERS += star_formation/none/star_formation_logger.h star_formation/none/star_formation_logger_struct.h 
nobase_noinst_HEADERS += cooling/none/cooling.h cooling/none/cooling_struct.h 
nobase_noinst_HEADERS += cooling/none/cooling_io.h cooling/none/cooling_properties.h  cooling/none/cooling_debug.h
nobase_noinst_HEADERS += cooling/const_du/cooling.h cooling/const_du/cooling_struct.h 
nobase_noinst_HEADERS += cooling/const_du/cooling_io.h cooling/const_du/cooling_properties.h cooling/const_du/cooling_debug.h
nobase_noinst_HEADERS += cooling/const_lambda/cooling.h cooling/const_lambda/cooling_struct.h 
nobase_noinst_HEADERS += cooling/const_lambda/cooling_io.h cooling/const_lambda/cooling_properties.h cooling/const_lambda/cooling_debug.h
nobase_noinst_HEADERS += cooling/grackle/cooling.h cooling/grackle/cooling_struct.h 
nobase_noinst_HEADERS += cooling/grackle/cooling_io.h cooling/grackle/cooling_properties.h cooling/grackle/cooling_debug.h
nobase_noinst_HEADERS += cooling/EAGLE/cooling.h cooling/EAGLE/cooling_struct.h cooling/EAGLE/cooling_tables.h 
nobase_noinst_HEADERS += cooling/EAGLE/cooling_io.h cooling/EAGLE/interpolate.h cooling/EAGLE/cooling_rates.h 
nobase_noinst_HEADERS += cooling/EAGLE/cooling_properties.h cooling/EAGLE/cooling_debug.h
nobase_noinst_HEADERS += cooling/QLA_EAGLE/cooling.h cooling/QLA_EAGLE/cooling_struct.h cooling/QLA_EAGLE/cooling_tables.h 
nobase_noinst_HEADERS += cooling/QLA_EAGLE/cooling_io.h cooling/QLA_EAGLE/interpolate.h cooling/QLA_EAGLE/cooling_rates.h 
nobase_noinst_HEADERS += cooling/QLA_EAGLE/cooling_properties.h cooling/QLA_EAGLE/cooling_debug.h
nobase_noinst_HEADERS += cooling/QLA/cooling.h cooling/QLA/cooling_struct.h cooling/QLA/cooling_tables.h 
nobase_noinst_HEADERS += cooling/QLA/cooling_io.h cooling/QLA/interpolate.h cooling/QLA/cooling_rates.h 
nobase_noinst_HEADERS += cooling/QLA/cooling_properties.h cooling/QLA/cooling_debug.h
nobase_noinst_HEADERS += cooling/PS2020/cooling.h cooling/PS2020/cooling_struct.h cooling/PS2020/cooling_subgrid.h 
nobase_noinst_HEADERS += cooling/PS2020/cooling_io.h cooling/PS2020/interpolate.h cooling/PS2020/cooling_rates.h 
nobase_noinst_HEADERS += cooling/PS2020/cooling_tables.h cooling/PS2020/cooling_subgrid.h
nobase_noinst_HEADERS += cooling/PS2020/cooling_properties.h cooling/PS2020/cooling_debug.h
nobase_noinst_HEADERS += chemistry/none/chemistry.h 
nobase_noinst_HEADERS += chemistry/none/chemistry_additions.h
nobase_noinst_HEADERS += chemistry/none/chemistry_io.h
nobase_noinst_HEADERS += chemistry/none/chemistry_csds.h
nobase_noinst_HEADERS += chemistry/none/chemistry_struct.h 
nobase_noinst_HEADERS += chemistry/none/chemistry_iact.h 
nobase_noinst_HEADERS += chemistry/none/chemistry_debug.h 
nobase_noinst_HEADERS += chemistry/GEAR/chemistry.h 
nobase_noinst_HEADERS += chemistry/GEAR/chemistry_io.h
nobase_noinst_HEADERS += chemistry/GEAR/chemistry_csds.h
nobase_noinst_HEADERS += chemistry/GEAR/chemistry_struct.h 
nobase_noinst_HEADERS += chemistry/GEAR/chemistry_iact.h 
nobase_noinst_HEADERS += chemistry/GEAR/chemistry_debug.h 
nobase_noinst_HEADERS += chemistry/GEAR_DIFFUSION/chemistry.h
nobase_noinst_HEADERS += chemistry/GEAR_DIFFUSION/chemistry_io.h
nobase_noinst_HEADERS += chemistry/GEAR_DIFFUSION/chemistry_struct.h
nobase_noinst_HEADERS += chemistry/GEAR_DIFFUSION/chemistry_iact.h
nobase_noinst_HEADERS += chemistry/GEAR_DIFFUSION/chemistry_debug.h
nobase_noinst_HEADERS += chemistry/AGORA/chemistry.h 
nobase_noinst_HEADERS += chemistry/AGORA/chemistry_io.h
nobase_noinst_HEADERS += chemistry/AGORA/chemistry_csds.h
nobase_noinst_HEADERS += chemistry/AGORA/chemistry_struct.h 
nobase_noinst_HEADERS += chemistry/AGORA/chemistry_iact.h 
nobase_noinst_HEADERS += chemistry/AGORA/chemistry_debug.h 
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry.h 
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry_additions.h
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry_io.h
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry_struct.h
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry_iact.h 
nobase_noinst_HEADERS += chemistry/EAGLE/chemistry_debug.h 
nobase_noinst_HEADERS += chemistry/QLA/chemistry.h 
nobase_noinst_HEADERS += chemistry/QLA/chemistry_io.h
nobase_noinst_HEADERS += chemistry/QLA/chemistry_struct.h
nobase_noinst_HEADERS += chemistry/QLA/chemistry_iact.h 
nobase_noinst_HEADERS += chemistry/QLA/chemistry_debug.h 
nobase_noinst_HEADERS += entropy_floor/none/entropy_floor.h 
nobase_noinst_HEADERS += entropy_floor/EAGLE/entropy_floor.h 
nobase_noinst_HEADERS += entropy_floor/QLA/entropy_floor.h 
nobase_noinst_HEADERS += tracers/none/tracers.h tracers/none/tracers_struct.h 
nobase_noinst_HEADERS += tracers/none/tracers_io.h tracers/none/tracers_debug.h
nobase_noinst_HEADERS += tracers/EAGLE/tracers.h tracers/EAGLE/tracers_struct.h 
nobase_noinst_HEADERS += tracers/EAGLE/tracers_io.h tracers/EAGLE/tracers_debug.h
nobase_noinst_HEADERS += extra_io/EAGLE/extra_io.h extra_io/EAGLE/extra.h
nobase_noinst_HEADERS += feedback/none/feedback.h feedback/none/feedback_struct.h feedback/none/feedback_iact.h 
nobase_noinst_HEADERS += feedback/none/feedback_properties.h feedback/none/feedback.h
nobase_noinst_HEADERS += feedback/none/feedback_debug.h
nobase_noinst_HEADERS += feedback/AGORA/feedback.h feedback/AGORA/feedback_struct.h feedback/AGORA/feedback_iact.h 
nobase_noinst_HEADERS += feedback/AGORA/feedback_properties.h feedback/AGORA/feedback.h
nobase_noinst_HEADERS += feedback/AGORA/feedback_debug.h
nobase_noinst_HEADERS += feedback/EAGLE_kinetic/feedback.h feedback/EAGLE_kinetic/feedback_struct.h 
nobase_noinst_HEADERS += feedback/EAGLE_kinetic/feedback_properties.h feedback/EAGLE_kinetic/feedback_iact.h 
nobase_noinst_HEADERS += feedback/EAGLE_kinetic/feedback_debug.h
nobase_noinst_HEADERS += feedback/EAGLE_thermal/feedback.h feedback/EAGLE_thermal/feedback_struct.h 
nobase_noinst_HEADERS += feedback/EAGLE_thermal/feedback_properties.h feedback/EAGLE_thermal/feedback_iact.h 
nobase_noinst_HEADERS += feedback/EAGLE_thermal/feedback_debug.h
nobase_noinst_HEADERS += feedback/EAGLE/yield_tables.h feedback/EAGLE/imf.h feedback/EAGLE/interpolate.h
nobase_noinst_HEADERS += feedback/EAGLE/enrichment.h
nobase_noinst_HEADERS += feedback/GEAR/stellar_evolution_struct.h feedback/GEAR/stellar_evolution.h 
nobase_noinst_HEADERS += feedback/GEAR/feedback.h feedback/GEAR/feedback_iact.h 
nobase_noinst_HEADERS += feedback/GEAR/feedback_properties.h feedback/GEAR/feedback_struct.h 
nobase_noinst_HEADERS += feedback/GEAR/initial_mass_function.h feedback/GEAR/supernovae_ia.h feedback/GEAR/supernovae_ii.h 
nobase_noinst_HEADERS += feedback/GEAR/lifetime.h feedback/GEAR/hdf5_functions.h feedback/GEAR/interpolation.h 
nobase_noinst_HEADERS += feedback/GEAR/feedback_debug.h
nobase_noinst_HEADERS += black_holes/Default/black_holes.h black_holes/Default/black_holes_io.h 
nobase_noinst_HEADERS += black_holes/Default/black_holes_part.h black_holes/Default/black_holes_iact.h 
nobase_noinst_HEADERS += black_holes/Default/black_holes_properties.h 
nobase_noinst_HEADERS += black_holes/Default/black_holes_struct.h
nobase_noinst_HEADERS += black_holes/Default/black_holes_debug.h
nobase_noinst_HEADERS += black_holes/EAGLE/black_holes.h black_holes/EAGLE/black_holes_io.h 
nobase_noinst_HEADERS += black_holes/EAGLE/black_holes_part.h black_holes/EAGLE/black_holes_iact.h 
nobase_noinst_HEADERS += black_holes/EAGLE/black_holes_properties.h black_holes/EAGLE/black_holes_parameters.h 
nobase_noinst_HEADERS += black_holes/EAGLE/black_holes_struct.h black_holes/EAGLE/black_holes_debug.h
nobase_noinst_HEADERS += black_holes/SPIN_JET/black_holes.h black_holes/SPIN_JET/black_holes_io.h 
nobase_noinst_HEADERS += black_holes/SPIN_JET/black_holes_part.h black_holes/SPIN_JET/black_holes_iact.h 
nobase_noinst_HEADERS += black_holes/SPIN_JET/black_holes_properties.h black_holes/SPIN_JET/black_holes_parameters.h
nobase_noinst_HEADERS += black_holes/SPIN_JET/black_holes_spin.h black_holes/SPIN_JET/black_holes_struct.h 
nobase_noinst_HEADERS += black_holes/SPIN_JET/black_holes_debug.h 
nobase_noinst_HEADERS += pressure_floor/GEAR/pressure_floor.h pressure_floor/none/pressure_floor.h 
nobase_noinst_HEADERS += pressure_floor/GEAR/pressure_floor_iact.h pressure_floor/none/pressure_floor_iact.h 
nobase_noinst_HEADERS += pressure_floor/GEAR/pressure_floor_struct.h pressure_floor/none/pressure_floor_struct.h 
nobase_noinst_HEADERS += pressure_floor/GEAR/pressure_floor_debug.h pressure_floor/none/pressure_floor_debug.h
nobase_noinst_HEADERS += sink/Default/sink.h sink/Default/sink_io.h sink/Default/sink_part.h sink/Default/sink_properties.h 
nobase_noinst_HEADERS += sink/Default/sink_iact.h sink/Default/sink_struct.h sink/Default/sink_debug.h
nobase_noinst_HEADERS += sink/GEAR/sink.h sink/GEAR/sink_io.h sink/GEAR/sink_part.h sink/GEAR/sink_properties.h 
nobase_noinst_HEADERS += sink/GEAR/sink_iact.h sink/GEAR/sink_struct.h sink/GEAR/sink_debug.h
nobase_noinst_HEADERS += neutrino.h neutrino_properties.h neutrino_io.h
nobase_noinst_HEADERS += neutrino/Default/neutrino.h neutrino/Default/relativity.h neutrino/Default/fermi_dirac.h
nobase_noinst_HEADERS += neutrino/Default/neutrino_properties.h neutrino/Default/neutrino_io.h
nobase_noinst_HEADERS += neutrino/Default/neutrino_response.h

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
