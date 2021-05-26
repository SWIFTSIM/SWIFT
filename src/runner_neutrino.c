/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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

/* This object's header. */
#include "runner.h"

/* Phase space density functions needed */
#include "neutrino.h"
#include "neutrino_properties.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"

/* Compute the dimensionless neutrino momentum (units of kb*T).
 *
 * @param v The internal 3-velocity
 * @param m_eV The neutrino mass in electron-volts
 * @param fac Conversion factor = 1. / (speed_of_light * T_nu_eV)
 */
INLINE static double neutrino_momentum(float v[3], double m_eV, double fac) {

  float v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  float vmag = sqrtf(v2);
  double p = vmag * fac * m_eV;
  return p;
}

/**
 * @brief Weight the active neutrino particles in a cell using the delta-f
 * method.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_neutrino_weighting(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->grav.count == 0) return;
  if (!cell_is_starting_gravity(c, e)) return;
  if (!with_cosmology)
    error("Phase space weighting without cosmology not implemented.");

  /* Retrieve physical and cosmological constants */
  const double c_vel = e->physical_constants->const_speed_light_c;
  const double *m_eV_array = e->cosmology->M_nu_eV;
  const int N_nu = e->cosmology->N_nu;
  const double T_eV = e->cosmology->T_nu_0_eV;
  const double fac = 1.0 / (c_vel * T_eV);
  const double inv_mass_factor = 1. / e->neutrino_mass_conversion_factor;
  const long long neutrino_seed = e->neutrino_properties->neutrino_seed;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_neutrino_weighting(r, c->progeny[k], 0);
  } else {
    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {
      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* Only act on neutrinos that needed to be kicked */
      if (!(gp->type == swift_type_neutrino && gpart_is_starting(gp, e)))
        continue;

      /* Use a particle id dependent seed */
      const long long seed = gp->id_or_neg_offset + neutrino_seed;

      /* Compute the initial dimensionless momentum from the seed */
      const double pi = neutrino_seed_to_fermi_dirac(seed);

      /* The neutrino mass (we cycle based on the neutrino seed) */
      const double m_eV = neutrino_seed_to_mass(N_nu, m_eV_array, seed);
      const double mass = m_eV * inv_mass_factor;

      /* Compute the current dimensionless momentum */
      double p = neutrino_momentum(gp->v_full, m_eV, fac);

      /* Compute the initial and current background phase-space density */
      double fi = fermi_dirac_density(pi);
      double f = fermi_dirac_density(p);
      double weight = 1.0 - f / fi;

      /* Set the statistically weighted mass */
      gp->mass = mass * weight;

      /* Prevent degeneracies */
      if (gp->mass == 0.) {
        gp->mass = FLT_MIN;
      }
    }
  }

  if (timer) TIMER_TOC(timer_weight);
}
