/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

/* Local headers */
#include "neutrino.h"

/* Compute the dimensionless neutrino momentum (units of kb*T).
 *
 * @param v The internal 3-velocity
 * @param m_eV The neutrino mass in electron-volts
 * @param fac Conversion factor = 1. / (speed_of_light * T_nu_eV)
 */
INLINE static double neutrino_momentum(const float v[3], const double m_eV,
                                       const double fac) {

  float v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  float vmag = sqrtf(v2);
  double p = vmag * fac * m_eV;
  return p;
}

void gather_neutrino_data(const struct space *s,
                          struct neutrino_data *nu_data) {
  nu_data->use_delta_f = 1;
  nu_data->m_eV_array = s->e->cosmology->M_nu_eV;
  nu_data->N_nu = s->e->cosmology->N_nu;
  nu_data->fac = 1.0 / (s->e->physical_constants->const_speed_light_c *
                        s->e->cosmology->T_nu_0_eV);
  nu_data->neutrino_seed = s->e->neutrino_properties->neutrino_seed;
}

/**
 * @brief Compute delta-f weight of a neutrino particle
 *
 * @param gp The #gpart.
 * @param nu_data Properties of the neutrino model
 * @param weight The resulting weight (output)
 */
void gpart_neutrino_weight(const struct gpart *gp,
                           const struct neutrino_data *nu_data,
                           double *weight) {

  /* Unpack neutrino model properties */
  const double *m_eV_array = nu_data->m_eV_array;
  const int N_nu = nu_data->N_nu;
  const double fac = nu_data->fac;
  const long long neutrino_seed = nu_data->neutrino_seed;

  /* Use a particle id dependent seed */
  const long long seed = gp->id_or_neg_offset + neutrino_seed;

  /* Compute the initial dimensionless momentum from the seed */
  const double pi = neutrino_seed_to_fermi_dirac(seed);

  /* The neutrino mass and degeneracy (we cycle based on the seed) */
  const double m_eV = neutrino_seed_to_mass(N_nu, m_eV_array, seed);

  /* Compute the current dimensionless momentum */
  double p = neutrino_momentum(gp->v_full, m_eV, fac);

  /* Compute the initial and current background phase-space density */
  double fi = fermi_dirac_density(pi);
  double f = fermi_dirac_density(p);
  *weight = 1.0 - f / fi;
}