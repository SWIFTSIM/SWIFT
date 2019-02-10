/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 James Willis (james.s.willis@durham.ac.uk).
 *               2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include "../config.h"

/* Includes. */
#include "swift.h"

#define Ntest 3
#define tol 1e-12
#define filename "output_list_params.yml"

/* Expected values from file */
const double time_values[Ntest] = {
    0.,
    10.,
    12.,
};

/* Expected values from file */
const double a_values[Ntest] = {
    0.01,
    0.1,
    0.5,
};

void test_no_cosmo(struct engine *e, char *name, int with_assert) {
  message("Testing output time for %s without cosmology", name);

  struct output_list *list;
  double delta_time = 0;
  double output_time = 0;

  /* Test Time */
  e->time_begin = 0;
  e->time_end = 14;
  e->time_base = (e->time_end - e->time_begin) / max_nr_timesteps;
  e->ti_current = 0;
  e->policy = !engine_policy_cosmology;

  /* initialize output_list */
  output_list_init(&list, e, name, &delta_time, &output_time);
  output_list_print(list);

  for (int i = 0; i < Ntest; i++) {
    /* Test last value */
    if (with_assert) {
      assert(fabs(output_time - time_values[i]) < tol);
    }

    /* Set current time */
    e->ti_current = (output_time - e->time_begin) / e->time_base;
    e->ti_current += 1;

    /* Read next value */
    integertime_t ti_next;
    output_list_read_next_time(list, e, name, &ti_next);

    output_time = (double)(ti_next * e->time_base) + e->time_begin;
  }

  output_list_clean(&list);
};

void test_cosmo(struct engine *e, char *name, int with_assert) {
  message("Testing output time for %s with cosmology", name);

  struct output_list *list;
  double delta_time = 0;
  double output_time = 0;

  /* Test Time */
  e->time_base = log(e->time_end / e->cosmology->a_begin) / max_nr_timesteps;
  e->ti_current = 0;
  e->policy = engine_policy_cosmology;

  /* initialize output_list */
  output_list_init(&list, e, name, &delta_time, &output_time);
  output_list_print(list);

  for (int i = 0; i < Ntest; i++) {
    /* Test last value */
    if (with_assert) {
      assert(fabs(output_time - a_values[i]) < tol);
    }

    /* Set current time */
    e->ti_current = log(output_time / e->cosmology->a_begin) / e->time_base;
    e->ti_current += 1;

    /* Read next value */
    integertime_t ti_next;
    output_list_read_next_time(list, e, name, &ti_next);

    output_time = (double)exp(ti_next * e->time_base) * e->cosmology->a_begin;
  }

  output_list_clean(&list);
};

int main(int argc, char *argv[]) {
  /* Create a structure to read file into. */
  struct swift_params params;

  /* Read the parameter file. */
  parser_read_file(filename, &params);

  /* initialization of unit system */
  struct unit_system us;
  units_init_from_params(&us, &params, "Units");

  /* initialization of phys_const */
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);

  /* initialization of cosmo */
  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);

  /* Pseudo initialization of engine */
  struct engine e;
  e.cosmology = &cosmo;
  e.parameter_file = &params;
  e.physical_constants = &phys_const;
  e.internal_units = &us;

  int with_assert = 1;
  int without_assert = 0;
  /* Test without cosmo */
  test_no_cosmo(&e, "Time", with_assert);

  /* Test with cosmo */
  test_cosmo(&e, "Redshift", with_assert);
  test_cosmo(&e, "ScaleFactor", with_assert);
  test_cosmo(&e, "Time", without_assert);

  cosmology_clean(&cosmo);

  /* Write message and leave */
  message("Test done");
  return 0;
}
