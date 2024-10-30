/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Willem Elbers (whe@willemelbers.com)
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
#include <config.h>

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"

/**
 * @brief Weight the active dcdm particles in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_dcdm_weighting(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->grav.count == 0) return;
  if (!cell_is_starting_gravity(c, e)) return;
  if (!with_cosmology)
    error("Decaying dark matter without cosmology not implemented.");

  /* Is it a hydro run? */
  const int is_hydro = (e->total_nr_parts > 0) ? 1 : 0;
  
  /* Compute the decaying + non-decaying (+ baryon if not hydro) matter density */
  const double Omega_cdm = e->cosmology->Omega_cdm;
  const double Omega_dcdm = cosmology_get_dcdm_density(e->cosmology, e->cosmology->a);
  const double Omega_b = e->cosmology->Omega_b;
  
  /* Also compute the matter density at the starting redshift */
  const double Omega_dcdm_begin = cosmology_get_dcdm_density(e->cosmology, e->cosmology->a_begin);
  
  double Omega_dm;
  double Omega_dm_begin;
  
  if (is_hydro) {
      Omega_dm = Omega_cdm + Omega_dcdm;
      Omega_dm_begin = Omega_cdm + Omega_dcdm_begin;
  } else {
      Omega_dm = Omega_cdm + Omega_dcdm + Omega_b;
      Omega_dm_begin = Omega_cdm + Omega_dcdm_begin + Omega_b;
  }
  
  /* The scaling factor for the particle masses */
  const double mass_scaling_factor = Omega_dm / Omega_dm_begin;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_dcdm_weighting(r, c->progeny[k], 0);
  } else {
    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {
      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* Only act on dcdm particles that needed to be kicked */
      if (!(gp->type == swift_type_dark_matter && gpart_is_starting(gp, e)))
        continue;

      /* Set the statistically weighted mass */
      gp->mass = gp->mass_ini * mass_scaling_factor;
    }
  }

  if (timer) TIMER_TOC(timer_weight);
}
