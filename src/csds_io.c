/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/* Config parameters. */
#include "../config.h"

#if defined(WITH_CSDS)

/* Some standard headers. */
#include "common_io.h"

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

/* This object's header. */
#include "csds_io.h"

/* Local includes. */
#include "version.h"

/**
 * @brief Write the parameters into a yaml file.
 *
 * @params log The #csds.
 * @params e The #engine.
 */
void csds_write_description(struct csds_writer* log, struct engine* e) {

  /* Only the master writes the description */
  if (e->nodeID != 0) {
    return;
  }

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%.100s.yml", e->csds->base_name);

  /* Open file */
  FILE* f = NULL;
  f = fopen(fileName, "w");

  if (f == NULL) {
    error("Failed to open file %s", fileName);
  }

  /* Write the header section */
  fprintf(f, "Header:\n");
  fprintf(f, "  BoxSize: [%g, %g, %g]\n", e->s->dim[0], e->s->dim[1],
          e->s->dim[2]);
  fprintf(f, "  Dimension: %i\n", (int)hydro_dimension);
  fprintf(f, "  RunName: %s\n", e->run_name);
  fprintf(f, "  Periodic: %i\n", e->s->periodic);
  /* Here I am using the local number of particles as
   * it is used to estimate the number of particles per logfile. */
  fprintf(f, "  NumberParts: %li\n", e->s->nr_parts);
  fprintf(f, "  NumberSParts: %li\n", e->s->nr_sparts);
  fprintf(f, "  NumberGParts: %li\n", e->s->nr_gparts);
  fprintf(f, "\n");

  /* Write the cosmology */
  fprintf(f, "Cosmology:\n");
  fprintf(f, "  Omega_cdm: %g\n", e->cosmology->Omega_cdm);
  fprintf(f, "  Omega_lambda: %g\n", e->cosmology->Omega_lambda);
  fprintf(f, "  Omega_b: %g\n", e->cosmology->Omega_b);
  fprintf(f, "  Omega_r: %g\n", e->cosmology->Omega_r);
  fprintf(f, "  Omega_k: %g\n", e->cosmology->Omega_k);
  fprintf(f, "  Omega_nu_0: %g\n", e->cosmology->Omega_nu_0);
  fprintf(f, "  w_0: %g\n", e->cosmology->w_0);
  fprintf(f, "  w_a: %g\n", e->cosmology->w_a);
  fprintf(f, "  Hubble0: %g\n", e->cosmology->H0);
  fprintf(f, "\n");

  /* Write unit system */
  const struct unit_system* us = e->internal_units;
  fprintf(f, "InternalUnitSystem:\n");
  fprintf(f, "  UnitMass_in_cgs: %g\n", units_get_base_unit(us, UNIT_MASS));
  fprintf(f, "  UnitLength_in_cgs: %g\n", units_get_base_unit(us, UNIT_LENGTH));
  fprintf(f, "  UnitTime_in_cgs: %g\n", units_get_base_unit(us, UNIT_TIME));
  fprintf(f, "  UnitCurrent_in_cgs: %g\n",
          units_get_base_unit(us, UNIT_CURRENT));
  fprintf(f, "  UnitTemp_in_cgs: %g\n",
          units_get_base_unit(us, UNIT_TEMPERATURE));
  fprintf(f, "\n");

  /* Write the code section */
  fprintf(f, "Code:\n");
  fprintf(f, "  Code: SWIFT\n");
  fprintf(f, "  CodeVersion: %s\n", package_version());
  fprintf(f, "  CompilerName: %s\n", compiler_name());
  fprintf(f, "  CompilerVersion: %s\n", compiler_version());
  fprintf(f, "  GitBranch: %s\n", git_branch());
  fprintf(f, "  GitRevision: %s\n", git_revision());
  fprintf(f, "  GitDate: %s\n", git_date());
  fprintf(f, "  ConfigurationOptions: %s\n", configuration_options());
  fprintf(f, "  RandomSeed: %i\n", SWIFT_RANDOM_SEED_XOR);
  fprintf(f, "\n");

  /* Write the policy section */
  fprintf(f, "Policy:\n");
  for (int i = 1; i < engine_maxpolicy; i++) {
    const int value = e->policy & (1 << i) ? 1 : 0;
    fprintf(f, "  %s: %i\n", engine_policy_names[i + 1], value);
  }
  fprintf(f, "\n");

  /* Close file */
  fclose(f);
}

#endif /* WITH_CSDS */
