/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2025 Katy Proctor (katy.proctor@fysik.su.se).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Some standard headers. */
#include <config.h>

/* Includes. */
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  // const char *base_name = "testSelectSIDMOutput";
  size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nspart = 0, Nbpart = 0,
         Nsink = 0, Nnupart = 0, Nsipart = 0;
  int flag_entropy_ICs = -1;
  int periodic = 1;
  double dim[3];
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;
  struct sink *sinks = NULL;
  struct sipart *siparts = NULL;

  struct gravity_props gravity_properties;
  struct hydro_props hydro_properties;
  struct stars_props stars_properties;
  struct sink_props sink_properties;
  struct neutrino_props neutrino_properties;
  struct neutrino_response neutrino_response;
  struct sidm_props sidm_properties;
  struct feedback_props feedback_properties;
  struct rt_props rt_properties;
  struct entropy_floor_properties entropy_floor;
  struct pressure_floor_props pressure_floor_props;
  struct black_holes_props black_holes_properties;
  struct fof_props fof_properties;
  struct lightcone_array_props lightcone_array_properties;
  struct ic_info ics_metadata;

  struct external_potential potential;
  struct forcing_terms forcing_terms;
  struct extra_io_properties extra_io_props;
  struct star_formation starform;
  struct pm_mesh mesh;
  struct power_spectrum_data pow_data;
  struct chemistry_global_data chemistry;
  struct cooling_function_data cooling_func;
  struct los_props los_properties;
  struct output_options output_options;
  struct repartition reparttype;

  strcpy(ics_metadata.group_name, "NoSUCH");

  if (argc < 2) error("Missing filename argument!");

  /* parse parameters */
  message("Reading parameters.");
  struct swift_params param_file;
  parser_read_file(argv[1], &param_file);

  /* Default unit system */
  message("Initialization of the unit system.");
  struct unit_system us;
  units_init_cgs(&us);

  /* Default physical constants */
  message("Initialization of the physical constants.");
  struct phys_const prog_const;
  phys_const_init(&us, &param_file, &prog_const);

  /* Read data */
  message("Reading initial conditions.");
  read_ic_single("input_SIDM.hdf5", &us, dim, &parts, &gparts, &sinks, &sparts,
                 &bparts, &siparts, &Ngas, &Ngpart, &Ngpart_background,
                 &Nnupart, &Nsink, &Nspart, &Nbpart, &Nsipart,
                 &flag_entropy_ICs,
                 /*with_hydro=*/1,
                 /*with_gravity=*/0,
                 /*with_sink=*/0,
                 /*with_stars=*/0,
                 /*with_black_holes=*/0,
                 /*with_sidm=*/1,
                 /*with_cosmology=*/0,
                 /*cleanup_h=*/0,
                 /*cleanup_sqrt_a=*/0,
                 /*h=*/1., /*a=*/1., /*n_threads=*/1, /*dry_run=*/0,
                 /*remap_ids=*/0, &ics_metadata);

  /* initialization of cosmology */
  message("Initialization of the cosmology.");
  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  /* pseudo initialization of the hydro properties */
  message("Initialization of the gas.");
  hydro_props_init(&hydro_properties, &prog_const, &us, &param_file);

  /* pseudo initialization of SIDM */
  message("Initialization of the SIDM.");
  sidm_props_init(&sidm_properties, &prog_const, &us, &param_file,
                  &hydro_properties, &cosmo);

  /* pseudo initialization of the space */
  message("Initialization of the space.");
  struct space s;
  space_init(&s, &param_file, &cosmo, dim, &hydro_properties, parts, gparts,
             sinks, sparts, bparts, siparts, Ngas, Ngpart, Nsink, Nspart,
             Nbpart, Nnupart, Nsipart, periodic,
             /*replicate=*/1,
             /*remap_ids=*/0,
             /*generate_gas_in_ics=*/0,
             /*with_hydro=*/1,
             /*with_self_gravity=*/0,
             /*with_star_formation=*/0,
             /*with_sinks=*/0,
             /*with_DM_particles=*/0,
             /*with_DM_background_particles=*/0,
             /*with_neutrinos=*/0,
             /*talking=*/0, /*dry_run=*/0, /*nr_nodes=*/1);

  message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
          s.dim[2]);
  message("space %s periodic.", s.periodic ? "is" : "isn't");
  message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
          s.cdim[1], s.cdim[2]);
  message("%zi siparts in %i cells.", s.nr_siparts, s.tot_cells);
  message("%zi parts in %i cells.", s.nr_parts, s.tot_cells);

  /* Verify that each particle is in its proper cell. */
  int icount = 0;
  space_map_cells_pre(&s, 0, &map_cellcheck, &icount);
  message("map_cellcheck picked up %i parts.", icount);

  /* Verify the maximal depth of cells. */
  int data[2] = {s.maxdepth, 0};
  space_map_cells_pre(&s, 0, &map_maxdepth, data);
  message("nr of cells at depth %i is %i.", data[0], data[1]);

  /* Construct the engine policy */
  int engine_policies = ENGINE_POLICY | engine_policy_steal;
  engine_policies |= engine_policy_sidm;
  engine_policies |= engine_policy_hydro;

  /* pseudo initialization of the engine */
  message("Initialization of the engine.");
  struct engine e;

  e.physical_constants = &prog_const;

  /* Initialize the engine with the space and policies. */
  engine_init(&e, &s, &param_file, &output_options, Ngas,
              /*Ngparts=*/0, /*Nsinks=*/0, /*Nstars=*/0,
              /*Nblackholes=*/0, /*Nbackground_gparts=*/0, /*Nnuparts=*/0,
              Nsipart, engine_policies,
              /*talking=*/0, &us, &prog_const, &cosmo, &hydro_properties,
              &entropy_floor, &gravity_properties, &stars_properties,
              &black_holes_properties, &sink_properties, &neutrino_properties,
              &neutrino_response, &sidm_properties, &feedback_properties,
              &pressure_floor_props, &rt_properties, &mesh, &pow_data,
              &potential, &forcing_terms, &cooling_func, &starform, &chemistry,
              &extra_io_props, &fof_properties, &los_properties,
              &lightcone_array_properties, &ics_metadata);

  engine_config(/*restart=*/0, /*fof=*/0, &e, &param_file,
                /*nr_nodes=*/1, /*myrank=*/0, /*nr_threads=*/1,
                /*nr_pool_threads=*/1, /*with_aff=*/0, /*talking=*/0, NULL,
                NULL, &reparttype);

  message(
      "Running on %zu gas particles, %zu sink particles, %zu stars "
      "particles %zu black hole particles, %zu neutrino particles, "
      "%zu CDM particles and %zu SIDM particles",
      Ngas, Nsink, Nspart, Nbpart, Nnupart, Ngpart, Nsipart);

  /* Clean-up */
  message("Cleaning memory.");
  cosmology_clean(&cosmo);
  engine_clean(&e, /*fof=*/0, /*restart=*/0);

  return 0;
}
