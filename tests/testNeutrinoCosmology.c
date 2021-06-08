/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Willem Elbers (whe@willemelbers.com)
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
#include <fenv.h>
#include <math.h>

/* Includes. */
#include "../config.h"
#include "swift.h"

#define N_CHECK 20
#define TOLERANCE 1e-6
#define EXTERNAL_TOLERANCE 1e-5  // comparing against CLASS

void test_params_init(struct swift_params *params, int testnr) {
  switch (testnr) {
    /* One doubly degenerate light neutrino (total 0.10 eV) and one massless */
    case 0:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_cdm:0.263069938");
      parser_set_param(params, "Cosmology:Omega_lambda:6.853646e-01");
      parser_set_param(params, "Cosmology:Omega_b:0.049198926683");
      parser_set_param(params, "Cosmology:h:0.6737");
      parser_set_param(params, "Cosmology:a_begin:1e-2");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:T_nu_0:1.95176");
      parser_set_param(params, "Cosmology:N_nu:1");
      parser_set_param(params, "Cosmology:N_ur:1.0196");
      parser_set_param(params, "Cosmology:M_nu_eV:0.0486");
      parser_set_param(params, "Cosmology:deg_nu:2");
      break;
    /* Three very massive neutrinos */
    case 1:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_cdm:0.241989232");
      parser_set_param(params, "Cosmology:Omega_lambda:0.685365");
      parser_set_param(params, "Cosmology:Omega_b:0.049199");
      parser_set_param(params, "Cosmology:h:0.6737");
      parser_set_param(params, "Cosmology:a_begin:1e-2");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:M_nu_eV:0.2,0.3,0.4");
      break;
    /* One light massive neutrino (0.05 eV) + two massless */
    case 2:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_cdm:0.261887181");
      parser_set_param(params, "Cosmology:Omega_lambda:0.685365");
      parser_set_param(params, "Cosmology:Omega_b:0.049199");
      parser_set_param(params, "Cosmology:h:0.6737");
      parser_set_param(params, "Cosmology:a_begin:1e-2");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:T_nu_0:1.95176");
      parser_set_param(params, "Cosmology:N_ur:2");
      parser_set_param(params, "Cosmology:N_nu:1");
      parser_set_param(params, "Cosmology:M_nu_eV:0.05");
      break;
    /* Three very massive neutrinos (same as 1), but starting earlier */
    case 3:
      parser_init("", params);
      parser_set_param(params, "Cosmology:Omega_cdm:0.241989232");
      parser_set_param(params, "Cosmology:Omega_lambda:0.685365");
      parser_set_param(params, "Cosmology:Omega_b:0.049199");
      parser_set_param(params, "Cosmology:h:0.6737");
      parser_set_param(params, "Cosmology:a_begin:1e-3");
      parser_set_param(params, "Cosmology:a_end:1.0");
      parser_set_param(params, "Cosmology:N_nu:3");
      parser_set_param(params, "Cosmology:M_nu_eV:0.2,0.3,0.4");
      break;
  }

  parser_set_param(params, "InternalUnitSystem:UnitMass_in_cgs:1.988435e+43");
  parser_set_param(params,
                   "InternalUnitSystem:UnitLength_in_cgs:3.0856775815e+24");
  parser_set_param(params,
                   "InternalUnitSystem:UnitVelocity_in_cgs:9.7846194238e+07");
  parser_set_param(params, "InternalUnitSystem:UnitCurrent_in_cgs:1");
  parser_set_param(params, "InternalUnitSystem:UnitTemp_in_cgs:1");
}

int main(int argc, char *argv[]) {

  /* Load a CLASS data table of background quantities for cosmology 0 */
  const int rows = 51;
  const int cols = 3;
  double CLASS_table[rows * cols];
  FILE *stream = fopen("testNeutrinoCosmology.dat", "r");
  if (stream == NULL) error("Could not open reference solution file!");
  char line[1024];
  int row = 0;
  while (fgets(line, 1024, stream)) {
    if (line[0] == '#') continue;
    char *tmp = strdup(line);
    char *ptr = NULL;
    CLASS_table[0 + row * 3] = strtod(tmp, &ptr);
    CLASS_table[1 + row * 3] = strtod(ptr + 1, &ptr);
    CLASS_table[2 + row * 3] = strtod(ptr + 1, &ptr);
    row++;
    free(tmp);
  }
  fclose(stream);

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Test a number of different cosmologies */
  int N_cosmos = 4;
  double times1[N_CHECK];  // compare t(a) for cosmology 1
  double times2[N_CHECK];  // with t(a) for cosmology 3
  for (int testnr = 0; testnr < N_cosmos; testnr++) {
    message("Initialization of cosmology %i", testnr);

    /* pseudo initialization of params */
    struct swift_params params;
    test_params_init(&params, testnr);

    /* initialization of unit system */
    struct unit_system us;
    units_init_from_params(&us, &params, "InternalUnitSystem");

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
      double time = cosmology_get_time_since_big_bang(&cosmo, a);

      /* Store the value for cosmologies 1 and 3 to compare later */
      if (testnr == 1) {
        times1[i] = time;
      } else if (testnr == 3) {
        times2[i] = time;
      }

      double my_a = cosmology_get_scale_factor(&cosmo, time);

      /* check accuracy */
      double rel_err = (my_a - a) / a;
      message("Accuracy of %g at a=%g", rel_err, a);
      assert(fabs(rel_err) < TOLERANCE);
    }

    if (testnr == 0) {
      /* For cosmology 0, compare against independently computed array */
      message("Relative to CLASS values (a, t, Omega_nu(a))");
      const double delta_a = (cosmo.log_a_end - cosmo.log_a_begin) / 10000;
      for (int j = 0; j < 50; j++) {
        double a1 = exp(cosmo.log_a_begin + delta_a * (j * 200 + 1));
        double time1 = cosmology_get_time_since_big_bang(&cosmo, a1);
        // double time1 = cosmo.time_interp_table[j*200] +
        // cosmo.time_interp_table_offset;
        double Onu1 = cosmology_get_neutrino_density(&cosmo, a1);

        double a2 = CLASS_table[0 + 3 * j];
        double time2 = CLASS_table[1 + 3 * j];
        double Onu2 = CLASS_table[2 + 3 * j];

        /* Class defines years as 365.25 days, we use 365 days (in this test).
         */
        time2 *= 365.25 / 365.0;

        assert(fabs(a1 - a2) / (a1 + a2) < EXTERNAL_TOLERANCE);
        assert(fabs(time1 - time2) / (time1 + time2) < EXTERNAL_TOLERANCE);
        assert(fabs(Onu1 - Onu2) / (Onu1 + Onu2) < EXTERNAL_TOLERANCE);

        message("At a = %e: (%e, %e, %e)", a1, a1 / a2, time1 / time2,
                Onu1 / Onu2);
      }
    }

    message("Everything seems fine with this cosmology.");

    cosmology_clean(&cosmo);
  }

  /* Compare cosmologies 0 and 3 */
  for (int j = 0; j < N_CHECK; j++) {
    /* check accuracy */
    double err = (times1[j] - times2[j]) / times1[j];
    message("Agreement of %g at step %i", err, j);
    assert(fabs(err) < TOLERANCE);
  }

  /* For the massive neutrino cosmology 3, check that it reproduces
   *  the relativistic (early) and non-relativistic (late) limits */

  message("Initialization...");

  /* pseudo initialization of params */
  struct swift_params params;
  test_params_init(&params, 3);

  /* initialization of unit system */
  struct unit_system us;
  units_init_cgs(&us);

  /* initialization of phys_const */
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);

  /* initialization of cosmo */
  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);

  message("Start checking the fermion integration...");

  /* Relativistic limit */
  double Omega_nu_early = cosmology_get_neutrino_density(&cosmo, 1e-10);
  double Omega_nu_r =
      cosmo.Omega_g * cosmo.N_eff * (7. / 8.) * pow(4. / 11., 4. / 3.);
  double err = (Omega_nu_early - Omega_nu_r) / Omega_nu_early;
  message("Accuracy of %g of the relativistic limit", err);
  assert(fabs(err) < TOLERANCE);

  /* Zeta function and other constants */
  const double zeta3 = 1.202056903159594;
  const double kb = phys_const.const_boltzmann_k;
  const double hbar = phys_const.const_planck_h / (2 * M_PI);
  const double cvel = phys_const.const_speed_light_c;
  const double eV = phys_const.const_electron_volt;

  /* Non-relativistic limit (limit not reached, so lower tolerance is fine)*/
  double Omega_nu_late = cosmology_get_neutrino_density(&cosmo, 1.0);
  double Omega_nu_nr = 0.;
  for (int i = 0; i < cosmo.N_nu; i++) {
    Omega_nu_nr += 6 * zeta3 / (11 * M_PI * M_PI) * pow(kb * cosmo.T_CMB_0, 3) /
                   pow(cvel * hbar, 3) * cosmo.M_nu_eV[i] * eV /
                   (cosmo.critical_density_0 * cvel * cvel);
  }
  double err2 = (Omega_nu_late - Omega_nu_nr) / Omega_nu_late;
  message("Accuracy of %g of the non-relativistic limit", err2);
  assert(fabs(err2) < 1e-5);

  message("Neutrino cosmology test successful.");

  cosmology_clean(&cosmo);
  return 0;
}
