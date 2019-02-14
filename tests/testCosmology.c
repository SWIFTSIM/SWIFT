/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausamman (loic.hausammann@epfl.ch)
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
#include "../config.h"

/* Includes. */
#include "swift.h"

#define N_CHECK 20
#define TOLERANCE 1e-7

void test_params_init(struct swift_params *params) {
  parser_init("", params);
  parser_set_param(params, "Cosmology:Omega_m:0.3075");
  parser_set_param(params, "Cosmology:Omega_lambda:0.6910");
  parser_set_param(params, "Cosmology:Omega_b:0.0486");
  parser_set_param(params, "Cosmology:h:0.6774");
  parser_set_param(params, "Cosmology:a_begin:0.1");
  parser_set_param(params, "Cosmology:a_end:1.0");
}

int main(int argc, char *argv[]) {

  message("Initialization...");

  /* pseudo initialization of params */
  struct swift_params params;
  test_params_init(&params);

  /* initialization of unit system */
  struct unit_system us;
  units_init_cgs(&us);

  /* initialization of phys_const */
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);

  /* initialization of cosmo */
  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);

  message("Start checking computation...");

  for (int i = 0; i < N_CHECK; i++) {
    double a = 0.1 + 0.9 * i / (N_CHECK - 1.);
    /* Compute a(t(a)) and check if same results */
    double tmp = cosmology_get_time_since_big_bang(&cosmo, a);
    tmp = cosmology_get_scale_factor(&cosmo, tmp);

    /* check accuracy */
    tmp = (tmp - a) / a;
    message("Accuracy of %g at a=%g", tmp, a);
    assert(fabs(tmp) < TOLERANCE);
  }

  message("Everything seems fine with cosmology.");

  cosmology_clean(&cosmo);
  return 0;
}
