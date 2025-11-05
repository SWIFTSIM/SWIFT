/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

void select_output_engine_init(struct engine *e, struct space *s,
                               struct cosmology *cosmo,
                               struct swift_params *params,
                               struct output_options *output,
                               struct cooling_function_data *cooling,
                               struct hydro_props *hydro_properties) {
  /* set structures */
  e->s = s;
  e->cooling_func = cooling;
  e->parameter_file = params;
  e->output_options = output;
  e->cosmology = cosmo;
  e->policy = engine_policy_hydro;
  e->hydro_properties = hydro_properties;

  /* initialization of threadpool */
  threadpool_init(&e->threadpool, 1);

  /* set parameters */
  e->verbose = 1;
  e->time = 0;
  e->snapshot_output_count = 0;
  e->snapshot_compression = 0;
};

void select_output_space_init(struct space *s, double *dim, int periodic,
                              size_t Ngas, size_t Nspart, size_t Ngpart,
                              struct part *parts, struct spart *sparts,
                              struct gpart *gparts) {
  s->periodic = periodic;
  for (int i = 0; i < 3; i++) {
    s->dim[i] = dim[i];
  }

  /* init space particles */
  s->nr_parts = Ngas;
  s->nr_sparts = Nspart;
  s->nr_gparts = Ngpart;

  s->parts = parts;
  s->gparts = gparts;
  s->sparts = sparts;

  /* Allocate the extra parts array for the gas particles. */
  if (posix_memalign((void **)&s->xparts, xpart_align,
                     Ngas * sizeof(struct xpart)) != 0)
    error("Failed to allocate xparts.");
  bzero(s->xparts, Ngas * sizeof(struct xpart));
};

void select_output_space_clean(struct space *s) { free(s->xparts); };

void select_output_engine_clean(struct engine *e) {
  threadpool_clean(&e->threadpool);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  // const char *base_name = "testSelectOutput";
  size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nspart = 0, Nbpart = 0,
         Nsink = 0, Nnupart = 0;
  int flag_entropy_ICs = -1;
  int periodic = 1;
  double dim[3];
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;
  struct sink *sinks = NULL;
  struct ic_info ics_metadata;
  strcpy(ics_metadata.group_name, "NoSUCH");

  /* parse parameters */
  message("Reading parameters.");
  struct swift_params param_file;
  const char *input_file = "selectOutputParameters.yml";
  parser_read_file(input_file, &param_file);

  struct output_options output_options;
  output_options_init(&param_file, 0, &output_options);

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
  read_ic_single("input.hdf5", &us, dim, &parts, &gparts, &sinks, &sparts,
                 &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart, &Nsink,
                 &Nspart, &Nbpart, &flag_entropy_ICs,
                 /*with_hydro=*/1,
                 /*with_gravity=*/0,
                 /*with_sink=*/0,
                 /*with_stars=*/0,
                 /*with_black_holes=*/0,
                 /*with_cosmology=*/0,
                 /*cleanup_h=*/0,
                 /*cleanup_sqrt_a=*/0,
                 /*h=*/1., /*a=*/1., /*n_threads=*/1, /*dry_run=*/0,
                 /*remap_ids=*/0, &ics_metadata);

  /* pseudo initialization of the space */
  message("Initialization of the space.");
  struct space s;
  select_output_space_init(&s, dim, periodic, Ngas, Nspart, Ngpart, parts,
                           sparts, gparts);

  /* initialization of cosmology */
  message("Initialization of the cosmology.");
  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  /* pseudo initialization of cooling */
  message("Initialization of the cooling.");
  struct cooling_function_data cooling;

  /* pseudo initialization of hydro */
  message("Initialization of the hydro.");
  struct hydro_props hydro_properties;
  hydro_props_init(&hydro_properties, &prog_const, &us, &param_file);

  /* pseudo initialization of the engine */
  message("Initialization of the engine.");
  struct engine e;
  e.physical_constants = &prog_const;
  sprintf(e.snapshot_base_name, "testSelectOutput");
  sprintf(e.run_name, "Select Output Test");
  select_output_engine_init(&e, &s, &cosmo, &param_file, &output_options,
                            &cooling, &hydro_properties);

  /* check output selection */
  message("Checking output parameters.");
  io_prepare_output_fields(&output_options, /*with_cosmology=*/0,
                           /*with_fof=*/0,
                           /*with_structure_finding=*/0,
                           /*verbose=*/1);

  /* write output file */
  message("Writing output.");
  write_output_single(&e, &us, &us, /*fof=*/0);

  /* Clean-up */
  message("Cleaning memory.");
  select_output_engine_clean(&e);
  select_output_space_clean(&s);
  cosmology_clean(&cosmo);
  free(parts);
  free(gparts);

  return 0;
}
