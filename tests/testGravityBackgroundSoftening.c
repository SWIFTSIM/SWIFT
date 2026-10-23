/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Softening factor of the background (PartType2) particles. */
int main(int argc, char *argv[]) {

  struct swift_params params;
  parser_init("", &params);
  parser_set_param(&params, "Cosmology:Omega_cdm:0.2589");
  parser_set_param(&params, "Cosmology:Omega_lambda:0.6910");
  parser_set_param(&params, "Cosmology:Omega_b:0.0486");
  parser_set_param(&params, "Cosmology:h:0.6774");
  parser_set_param(&params, "Cosmology:a_begin:0.1");
  parser_set_param(&params, "Cosmology:a_end:1.0");
  parser_set_param(&params, "Gravity:eta:0.025");
  parser_set_param(&params, "Gravity:MAC:geometric");
  parser_set_param(&params, "Gravity:theta_cr:0.7");
  parser_set_param(&params, "Gravity:mesh_side_length:32");
  parser_set_param(&params, "Gravity:comoving_DM_softening:0.002");
  parser_set_param(&params, "Gravity:max_physical_DM_softening:0.002");
  parser_set_param(&params, "Gravity:comoving_baryon_softening:0.002");
  parser_set_param(&params, "Gravity:max_physical_baryon_softening:0.002");
  parser_set_param(&params, "Gravity:softening_ratio_background:0.2");
  parser_set_param(&params, "Gravity:max_adaptive_softening:1.0");
  parser_set_param(&params, "Gravity:min_adaptive_softening:0.0");

  struct unit_system us;
  units_init_cgs(&us);
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);

  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);
  struct cosmology no_cosmo;
  cosmology_init_no_cosmo(&no_cosmo);

  const double ratio = 0.2;
  const double dim[3] = {1.5e25, 2.0e25, 2.5e25};
  const int cdim[3] = {8, 8, 8};
  const double volume = dim[0] * dim[1] * dim[2];
  const double rho_mean_cosmo =
      (cosmo.Omega_cdm + cosmo.Omega_b) * cosmo.critical_density_0;

  struct gravity_props p_cosmo;
  bzero(&p_cosmo, sizeof(struct gravity_props));
  gravity_props_init(&p_cosmo, &params, &phys_const, &cosmo,
                     /*with_cosmology=*/1, /*with_external_potential=*/0,
                     /*has_baryons=*/1, /*has_DM=*/1, /*has_neutrinos=*/0,
                     /*is_zoom_simulation=*/1, /*total_matter_mass=*/0.,
                     /*dry_run=*/0, /*periodic=*/1, dim, cdim);
  const float fac_cosmo_expected = kernel_gravity_softening_plummer_equivalent *
                                   ratio * cbrt(1. / rho_mean_cosmo);
  if (fabsf(p_cosmo.epsilon_background_fac - fac_cosmo_expected) >
      1e-6f * fac_cosmo_expected)
    error("Cosmological factor %.9e, expected %.9e",
          p_cosmo.epsilon_background_fac, fac_cosmo_expected);

  /* Non-cosmological path with the same mean density as above */
  const double total_mass = rho_mean_cosmo * volume;
  struct gravity_props p_no_cosmo;
  bzero(&p_no_cosmo, sizeof(struct gravity_props));
  gravity_props_init(&p_no_cosmo, &params, &phys_const, &no_cosmo,
                     /*with_cosmology=*/0, /*with_external_potential=*/0,
                     /*has_baryons=*/1, /*has_DM=*/1, /*has_neutrinos=*/0,
                     /*is_zoom_simulation=*/1, total_mass, /*dry_run=*/0,
                     /*periodic=*/1, dim, cdim);
  const float fac_no_cosmo_expected =
      kernel_gravity_softening_plummer_equivalent * ratio *
      cbrt(volume / total_mass);
  if (!(p_no_cosmo.epsilon_background_fac > 0.f))
    error("Non-cosmological factor is not positive: %e",
          p_no_cosmo.epsilon_background_fac);
  if (fabsf(p_no_cosmo.epsilon_background_fac - fac_no_cosmo_expected) >
      1e-6f * fac_no_cosmo_expected)
    error("Non-cosmological factor %.9e, expected %.9e",
          p_no_cosmo.epsilon_background_fac, fac_no_cosmo_expected);
  if (fabsf(p_no_cosmo.epsilon_background_fac - fac_cosmo_expected) >
      1e-6f * fac_cosmo_expected)
    error("Same mean density gives %.9e without and %.9e with cosmology",
          p_no_cosmo.epsilon_background_fac, fac_cosmo_expected);

  /* Softening of a particle of mass m: ratio times its mean spacing */
  const double m = 1e-3 * total_mass;
  const double eps = p_no_cosmo.epsilon_background_fac * cbrt(m);
  const double eps_expected = kernel_gravity_softening_plummer_equivalent *
                              ratio * cbrt(m * volume / total_mass);
  if (fabs(eps - eps_expected) > 1e-6 * eps_expected)
    error("Softening %.9e, expected %.9e", eps, eps_expected);

  message("Background softening factor: cosmo %.9e, no cosmo %.9e. OK.",
          p_cosmo.epsilon_background_fac, p_no_cosmo.epsilon_background_fac);

  return 0;
}
