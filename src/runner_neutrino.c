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
#include <config.h>

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
  struct neutrino_model nu_model;
  gather_neutrino_consts(e->s, &nu_model);

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

      /* Compute the mass and delta-f weight */
      double mass, weight;
      gpart_neutrino_mass_weight(gp, &nu_model, &mass, &weight);

      /* Set the statistically weighted mass */
      gp->mass = mass * weight;

      /* Prevent degeneracies */
      if (gp->mass == 0.) {
        gp->mass = FLT_MIN;
      }
    }
  }

  if (timer) TIMER_TOC(timer_neutrino_weighting);
}
