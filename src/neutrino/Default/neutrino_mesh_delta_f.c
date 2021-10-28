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

/**
 * @brief Gather neutrino constants
 *
 * @param s The #space for this run.
 * @param ncs Struct with neutrino constants
 */
void gather_neutrino_consts(const struct space *s,
                            struct neutrino_consts *ncs) {
  ncs->use_mesh_delta_f = s->e->neutrino_properties->use_delta_f_mesh_only;
  ncs->m_eV_array = s->e->cosmology->M_nu_eV;
  ncs->N_nu = s->e->cosmology->N_nu;
  ncs->fac = 1.0 / (s->e->physical_constants->const_speed_light_c *
                    s->e->cosmology->T_nu_0_eV);
  ncs->neutrino_seed = s->e->neutrino_properties->neutrino_seed;
}

/**
 * @brief Compute delta-f weight of a neutrino particle
 *
 * @param gp The #gpart.
 * @param nu_consts Properties of the neutrino model
 * @param weight The resulting weight (output)
 */
void gpart_neutrino_weight(const struct gpart *gp,
                           const struct neutrino_consts *nu_consts,
                           double *weight) {

  /* Unpack neutrino model properties */
  const double *m_eV_array = nu_consts->m_eV_array;
  const int N_nu = nu_consts->N_nu;
  const double fac = nu_consts->fac;
  const long long neutrino_seed = nu_consts->neutrino_seed;

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