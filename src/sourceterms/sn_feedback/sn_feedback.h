/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_SN_FEEDBACK_H
#define SWIFT_SN_FEEDBACK_H
#include <float.h>
/* Config parameters. */
#include "../config.h"

#include "engine.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "runner.h"
#include "timestep.h"

/**
 * @file src/sourceterms/sn_feedback.h
 *
 * @brief Routines related to sourceterms (supernova feedback): determine if
 * feedback occurs in this cell
 *
 * @param cell_min: corner of cell to test
 * @param cell_width: width of cell to test
 * @param sourceterms: properties of source terms to test
 * @param dimen: dimensionality of the problem
 *
 * This routine tests whether a source term should be applied to this cell
 * return: 1 if yes, return: 0 if no
 */
__attribute__((always_inline)) INLINE static int supernova_feedback_test_cell(
    const double cell_min[], const double cell_width[],
    struct sourceterms* sourceterms, const int dimen) {
  if (sourceterms->supernova.status == supernova_is_done) return 0;

  const double location[3] = {sourceterms->supernova.x,
                              sourceterms->supernova.y,
                              sourceterms->supernova.z};
  for (int i = 0; i < dimen; i++) {
    if (cell_min[i] > location[i]) return 0;
    if ((cell_min[i] + cell_width[i]) <= location[i]) return 0;
  };
  return 1;
};

/**
 * @file src/sourceterms/sn_feedback.h
 *
 * @brief Routines related to source terms (supernova feedback): perform
 * feedback in this cell
 * @param r: the runner
 * @param sourceterms the structure describing the source terms properties
 * @param c the cell to apply feedback to
 *
 * This routine heats an individual particle (p), increasing its thermal energy
 * per unit mass
 *      by supernova energy / particle mass.
 */
__attribute__((always_inline)) INLINE static void supernova_feedback_apply(
    struct runner* restrict r, struct sourceterms* restrict sourceterms,
    struct cell* restrict c) {

  const int count = c->count;
  struct part* restrict parts = c->parts;
  struct xpart* restrict xparts = c->xparts;
  const double timeBase = r->e->timeBase;
  const int ti_current = r->e->ti_current;

  /* inject SN energy into the particle with highest id in this cell if it is
   * active */
  int imax = 0;
  struct part* restrict p_sn = NULL;
  struct xpart* restrict xp_sn = NULL;

  for (int i = 0; i < count; i++) {

    /* Get a direct pointer on the part. */
    struct part* restrict p = &parts[i];
    if (p->id > imax) {
      imax = p->id;
      p_sn = p;
      xp_sn = &xparts[i];
    }
  }

  /* Is this part within the time step? */
  if (p_sn->ti_begin == ti_current) {

    /* Does this time step straddle the feedback injection time? */
    const float t_begin = p_sn->ti_begin * timeBase;
    const float t_end = p_sn->ti_end * timeBase;
    if (t_begin <= sourceterms->supernova.time &&
        t_end > sourceterms->supernova.time) {

      /* store old time step */
      const int dti_old = p_sn->ti_end - p_sn->ti_begin;

      /* add supernova feedback */
      const float u_old = hydro_get_internal_energy(p_sn, 0);
      const float ent_old = hydro_get_entropy(p_sn, 0.0);
      const float u_new =
          u_old + sourceterms->supernova.energy / hydro_get_mass(p_sn);
      hydro_set_internal_energy(p_sn, u_new);
      const float u_set = hydro_get_internal_energy(p_sn, 0.0);
      const float ent_set = hydro_get_entropy(p_sn, 0.0);
      message(
          " applied super nova, time = %e, location= %e %e %e velocity= %e %e "
          "%e",
          ti_current * timeBase, p_sn->x[0], p_sn->x[1], p_sn->x[2], p_sn->v[0],
          p_sn->v[1], p_sn->v[2]);
      message(
          " injected SN energy in particle = %lld, increased energy from %e to "
          "%e and is notw %e, entropy from %e to %e",
          p_sn->id, u_old, u_new, u_set, ent_old, ent_set);

      /* label supernova as done */
      sourceterms->supernova.status = supernova_is_done;

      /* update timestep if new time step shorter than old time step */
      const int dti = get_part_timestep(p_sn, xp_sn, r->e);
      if (dti < dti_old) {
        p_sn->ti_end = p_sn->ti_begin + dti;
        message(" changed timestep from %d to %d", dti_old, dti);

        /* apply simple time-step limiter on all particles in same cell:
         */
        int i_limit = 0;
        for (int i = 0; i < count; i++) {
          struct part* restrict p = &parts[i];
          const int dti_old = p->ti_end - p->ti_begin;
          if (dti_old > 2 * dti) {
            i_limit++;
            const int dti_new = 2 * dti;
            p->ti_end = p->ti_begin + dti_new;
            message(" old step = %d new step = %d", dti_old, dti_new);
          } else
            message(" old step = %d", dti_old);
        }
        message(" count= %d limited timestep of %d particles ", count, i_limit);
      } /* end of limiter */
      error("end");
    }
  }
};

/**
 * @file src/sourceterms/sn_feedback.h
 *
 * @brief Routine to initialise supernova feedback
 * @param parameterfile: the parse parmeter file
 * @param us: the unit system in use
 * @param sourceterms the structure describing the source terms properties
 *
 * This routine heats an individual particle (p), increasing its thermal energy
 * per unit mass
 *      by supernova energy / particle mass.
 */

__attribute__((always_inline)) INLINE static void supernova_init(
    struct swift_params* parameter_file, struct unit_system* us,
    struct sourceterms* source) {
  source->supernova.time = parser_get_param_double(parameter_file, "SN:time");
  source->supernova.energy =
      parser_get_param_double(parameter_file, "SN:energy");
  source->supernova.x = parser_get_param_double(parameter_file, "SN:x");
  source->supernova.y = parser_get_param_double(parameter_file, "SN:y");
  source->supernova.z = parser_get_param_double(parameter_file, "SN:z");
  source->supernova.status = supernova_is_not_done;
}
__attribute__((always_inline)) INLINE static void supernova_print(
    struct sourceterms* source) {
  message(
      " Single SNe of energy= %e will explode at time= %e at location "
      "(%e,%e,%e)",
      source->supernova.energy, source->supernova.time, source->supernova.x,
      source->supernova.y, source->supernova.z);
}
#endif /* SWIFT_SN_FEEDBACK_H */
