/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifndef SWIFT_PART_INIT_H
#define SWIFT_PART_INIT_H

/* Config parameters. */
#include <config.h>

/* Local headers */
#include "adaptive_softening.h"
#include "black_holes.h"
#include "chemistry.h"
#include "engine.h"
#include "feedback.h"
#include "gravity.h"
#include "mhd.h"
#include "rt.h"
#include "sink.h"
#include "star_formation.h"
#include "stars.h"
#include "threadpool.h"
#include "tracers.h"

/**
 * @brief Call all the initialisation functions of the start of a step for a gas
 * particle.
 * @param p The #part.
 * @param xp The extended particle data.
 * @param e The #engine.
 */
static INLINE void part_init(struct part *p, struct xpart *xp,
                             const struct engine *e) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);

  hydro_init_part(p, &e->s->hs);
  adaptive_softening_init_part(p);
  mhd_init_part(p);
  black_holes_init_potential(&p->black_holes_data);
  chemistry_init_part(p, e->chemistry);
  star_formation_init_part(p, e->star_formation);
  tracers_after_init(p, xp, e->internal_units, e->physical_constants,
                     with_cosmology, e->cosmology, e->hydro_properties,
                     e->cooling_func, e->time);
  sink_init_part(p, e->sink_properties);
  rt_init_part(p);
}

/**
 * @brief Call all the initialisation functions of the start of a step for a
 * gravity particle.
 * @param gp The #gpart.
 * @param e The #engine.
 */
static INLINE void gpart_init(struct gpart *gp, const struct engine *e) {

  gravity_init_gpart(gp);
}

/**
 * @brief Call all the initialisation functions of the start of a step for a
 * star particle.
 * @param sp The #spart.
 * @param e The #engine.
 */
static INLINE void spart_init(struct spart *sp, const struct engine *e) {

  stars_init_spart(sp);
  feedback_init_spart(sp);
  rt_init_spart(sp);
}

/**
 * @brief Call all the initialisation functions of the start of a step for a
 * black hole particle.
 * @param bp The #bpart.
 * @param e The #engine.
 */
static INLINE void bpart_init(struct bpart *bp, const struct engine *e) {

  black_holes_init_bpart(bp);
}

/**
 * @brief Call all the initialisation functions of the start of a step for a
 * sink particle.
 * @param si The #sink.
 * @param e The #engine.
 */
static INLINE void sink_init(struct sink *si, const struct engine *e) {
  sink_init_sink(si);
}

#endif /* SWIFT_PART_INIT_H */
