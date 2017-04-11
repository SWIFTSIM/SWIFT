/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_DOIACT_GRAV_H
#define SWIFT_RUNNER_DOIACT_GRAV_H

/* Includes. */
#include "cell.h"
#include "gravity.h"
#include "part.h"

/**
 * @brief Recursively propagate the multipoles down the tree by applying the
 * L2L and L2P kernels.
 *
 * @param r The #runner.
 * @param c The #cell we are working on.
 * @param timer Are we timing this ?
 */
void runner_do_grav_down(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int periodic = e->s->periodic;
  struct gpart *gparts = c->gparts;
  const int gcount = c->gcount;

  TIMER_TIC;

  if (c->split) { /* Node case */

    /* Add the field-tensor to all the 8 progenitors */
    for (int k = 0; k < 8; ++k) {
      struct cell *cp = c->progeny[k];
      struct grav_tensor temp;

      /* Do we have a progenitor with any active g-particles ? */
      if (cp != NULL && cell_is_active(cp, e)) {

        /* Shift the field tensor */
        gravity_L2L(&temp, &c->multipole->pot, cp->multipole->CoM,
                    c->multipole->CoM, 0 * periodic);

        /* Add it to this level's tensor */
        gravity_field_tensors_add(&cp->multipole->pot, &temp);

        /* Recurse */
        runner_do_grav_down(r, cp, 0);
      }
    }

  } else { /* Leaf case */

    /* Apply accelerations to the particles */
    for (int i = 0; i < gcount; ++i) {

      /* Get a handle on the gpart */
      struct gpart *gp = &gparts[i];

      /* Update if active */
      if (gpart_is_active(gp, e))
        gravity_L2P(&c->multipole->pot, c->multipole->CoM, gp);
    }
  }

  if (timer) TIMER_TOC(timer_dograv_down);
}

/**
 * @brief Computes the interaction of the field tensor in a cell with the
 * multipole of another cell.
 *
 * @param r The #runner.
 * @param ci The #cell with field tensor to interact.
 * @param cj The #cell with the multipole.
 */
void runner_dopair_grav_mm(const struct runner *r, struct cell *restrict ci,
                           struct cell *restrict cj) {

  const struct engine *e = r->e;
  const int periodic = e->s->periodic;
  const struct multipole *multi_j = &cj->multipole->m_pole;
  // const float a_smooth = e->gravity_properties->a_smooth;
  // const float rlr_inv = 1. / (a_smooth * ci->super->width[0]);

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_j->M_000 == 0.f) error("Multipole does not seem to have been set.");
#endif

  /* Anything to do here? */
  if (!cell_is_active(ci, e)) return;

  /* Do we need to drift the multipole ? */
  if (cj->ti_old_multipole != e->ti_current) cell_drift_multipole(cj, e);

  /* Let's interact at this level */
  gravity_M2L(&ci->multipole->pot, multi_j, ci->multipole->CoM,
              cj->multipole->CoM, periodic * 0);

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Computes the interaction of all the particles in a cell with the
 * multipole of another cell.
 *
 * @param r The #runner.
 * @param ci The #cell with particles to interct.
 * @param cj The #cell with the multipole.
 */
void runner_dopair_grav_pm(const struct runner *r,
                           const struct cell *restrict ci,
                           const struct cell *restrict cj) {

  error("Function should not be called");
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 *
 * @todo Use a local cache for the particles.
 */
void runner_dopair_grav_pp(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;
  const float a_smooth = e->gravity_properties->a_smooth;
  const float rlr_inv = 1. / (a_smooth * ci->super->width[0]);

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->width[0] != cj->width[0])
    error("Non matching cell sizes !! h_i=%f h_j=%f", ci->width[0],
          cj->width[0]);
#endif

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

#if ICHECK > 0
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &gparts_i[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count= %d",
              gp->id_or_neg_offset, cj->loc[0], cj->loc[1], cj->loc[2],
              cj->width[0], cj->gcount);
  }

  for (int pid = 0; pid < gcount_j; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &gparts_j[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count=%d",
              gp->id_or_neg_offset, ci->loc[0], ci->loc[1], ci->loc[2],
              ci->width[0], ci->gcount);
  }
#endif

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gpi = &gparts_i[pid];

    if (!gpart_is_active(gpi, e)) continue;

    /* Loop over every particle in the other cell. */
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Get a hold of the jth part in cj. */
      const struct gpart *restrict gpj = &gparts_j[pjd];

      /* Compute the pairwise distance. */
      const float dx[3] = {gpi->x[0] - gpj->x[0],   // x
                           gpi->x[1] - gpj->x[1],   // y
                           gpi->x[2] - gpj->x[2]};  // z
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Interact ! */
      runner_iact_grav_pp_nonsym(rlr_inv, r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
      gpi->num_interacted++;
#endif
    }
  }

  /* Loop over all particles in cj... */
  for (int pjd = 0; pjd < gcount_j; pjd++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gpj = &gparts_j[pjd];

    if (!gpart_is_active(gpj, e)) continue;

    /* Loop over every particle in the other cell. */
    for (int pid = 0; pid < gcount_i; pid++) {

      /* Get a hold of the ith part in ci. */
      const struct gpart *restrict gpi = &gparts_i[pid];

      /* Compute the pairwise distance. */
      const float dx[3] = {gpj->x[0] - gpi->x[0],   // x
                           gpj->x[1] - gpi->x[1],   // y
                           gpj->x[2] - gpi->x[2]};  // z
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Interact ! */
      runner_iact_grav_pp_nonsym(rlr_inv, r2, dx, gpj, gpi);

#ifdef SWIFT_DEBUG_CHECKS
      gpj->num_interacted++;
#endif
    }
  }

  TIMER_TOC(timer_dopair_grav_pp);
}

/**
 * @brief Computes the interaction of all the particles in a cell directly
 *
 * @param r The #runner.
 * @param c The #cell.
 *
 * @todo Use a local cache for the particles.
 */
void runner_doself_grav_pp(struct runner *r, struct cell *c) {

  const struct engine *e = r->e;
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;
  const float a_smooth = e->gravity_properties->a_smooth;
  const float rlr_inv = 1. / (a_smooth * c->super->width[0]);

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->gcount == 0) error("Doing self gravity on an empty cell !");
#endif

  /* Anything to do here? */
  if (!cell_is_active(c, e)) return;

#if ICHECK > 0
  for (int pid = 0; pid < gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &gparts[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count= %d",
              gp->id_or_neg_offset, c->loc[0], c->loc[1], c->loc[2],
              c->width[0], c->gcount);
  }
#endif

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gpi = &gparts[pid];

    /* Loop over every particle in the other cell. */
    for (int pjd = pid + 1; pjd < gcount; pjd++) {

      /* Get a hold of the jth part in ci. */
      struct gpart *restrict gpj = &gparts[pjd];

      /* Compute the pairwise distance. */
      float dx[3] = {gpi->x[0] - gpj->x[0],   // x
                     gpi->x[1] - gpj->x[1],   // y
                     gpi->x[2] - gpj->x[2]};  // z
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Interact ! */
      if (gpart_is_active(gpi, e) && gpart_is_active(gpj, e)) {

        runner_iact_grav_pp(rlr_inv, r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
        gpi->num_interacted++;
        gpj->num_interacted++;
#endif

      } else {

        if (gpart_is_active(gpi, e)) {

          runner_iact_grav_pp_nonsym(rlr_inv, r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
          gpi->num_interacted++;
#endif

        } else if (gpart_is_active(gpj, e)) {

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          runner_iact_grav_pp_nonsym(rlr_inv, r2, dx, gpj, gpi);

#ifdef SWIFT_DEBUG_CHECKS
          gpj->num_interacted++;
#endif
        }
      }
    }
  }

  TIMER_TOC(timer_doself_grav_pp);
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param gettimer Are we timing this ?
 *
 * @todo Use a local cache for the particles.
 */
void runner_dopair_grav(struct runner *r, struct cell *ci, struct cell *cj,
                        int gettimer) {

#ifdef SWIFT_DEBUG_CHECKS

  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;

  /* Early abort? */
  if (gcount_i == 0 || gcount_j == 0)
    error("Doing pair gravity on an empty cell !");

  /* Bad stuff will happen if cell sizes are different */
  if (ci->width[0] != cj->width[0])
    error("Non matching cell sizes !! h_i=%f h_j=%f", ci->width[0],
          cj->width[0]);

  /* Sanity check */
  if (ci == cj)
    error(
        "The impossible has happened: pair interaction between a cell and "
        "itself.");

  /* Are the cells direct neighbours? */
  if (!cell_are_neighbours(ci, cj))
    error(
        "Non-neighbouring cells ! ci->x=[%f %f %f] ci->width=%f cj->loc=[%f %f "
        "%f] "
        "cj->width=%f",
        ci->loc[0], ci->loc[1], ci->loc[2], ci->width[0], cj->loc[0],
        cj->loc[1], cj->loc[2], cj->width[0]);

#endif

#if ICHECK > 0
  for (int pid = 0; pid < ci->gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &ci->gparts[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count= %d",
              gp->id_or_neg_offset, cj->loc[0], cj->loc[1], cj->loc[2],
              cj->width[0], cj->gcount);
  }

  for (int pid = 0; pid < cj->gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &cj->gparts[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count= %d",
              gp->id_or_neg_offset, ci->loc[0], ci->loc[1], ci->loc[2],
              ci->width[0], ci->gcount);
  }
#endif

  TIMER_TIC;

  /* Are both cells split ? */
  if (ci->split && cj->split) {

    for (int j = 0; j < 8; j++) {
      if (ci->progeny[j] != NULL) {

        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL) {

            if (cell_are_neighbours(ci->progeny[j], cj->progeny[k])) {

              /* Recurse */
              runner_dopair_grav(r, ci->progeny[j], cj->progeny[k], 0);

            } else {

              /* Ok, here we can go for multipole-multipole interactions */
              runner_dopair_grav_mm(r, ci->progeny[j], cj->progeny[k]);
              runner_dopair_grav_mm(r, cj->progeny[k], ci->progeny[j]);
            }
          }
        }
      }
    }
  } else { /* Not split */

    /* Compute the interactions at this level directly. */
    runner_dopair_grav_pp(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
}

/**
 * @brief Computes the interaction of all the particles in a cell
 *
 * @param r The #runner.
 * @param c The first #cell.
 * @param gettimer Are we timing this ?
 *
 * @todo Use a local cache for the particles.
 */
void runner_doself_grav(struct runner *r, struct cell *c, int gettimer) {

#ifdef SWIFT_DEBUG_CHECKS

  /* Early abort? */
  if (c->gcount == 0) error("Doing self gravity on an empty cell !");
#endif

  TIMER_TIC;

  /* If the cell is split, interact each progeny with itself, and with
     each of its siblings. */
  if (c->split) {

    for (int j = 0; j < 8; j++) {
      if (c->progeny[j] != NULL) {

        runner_doself_grav(r, c->progeny[j], 0);

        for (int k = j + 1; k < 8; k++) {
          if (c->progeny[k] != NULL) {

            runner_dopair_grav(r, c->progeny[j], c->progeny[k], 0);
          }
        }
      }
    }
  }

  /* If the cell is not split, then just go for it... */
  else {

    runner_doself_grav_pp(r, c);
  }

  if (gettimer) TIMER_TOC(timer_dosub_self_grav);
}

void runner_dosub_grav(struct runner *r, struct cell *ci, struct cell *cj,
                       int timer) {

  /* Is this a single cell? */
  if (cj == NULL) {

    runner_doself_grav(r, ci, 1);

  } else {

#ifdef SWIFT_DEBUG_CHECKS
    if (!cell_are_neighbours(ci, cj))
      error("Non-neighbouring cells in pair task !");
#endif

    runner_dopair_grav(r, ci, cj, 1);
  }
}

void runner_do_grav_long_range(struct runner *r, struct cell *ci, int timer) {

#if ICHECK > 0
  for (int pid = 0; pid < ci->gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &ci->gparts[pid];

    if (gp->id_or_neg_offset == ICHECK)
      message("id=%lld loc=[ %f %f %f ] size= %f count= %d",
              gp->id_or_neg_offset, ci->loc[0], ci->loc[1], ci->loc[2],
              ci->width[0], ci->gcount);
  }
#endif

  TIMER_TIC;

  /* Recover the list of top-level cells */
  const struct engine *e = r->e;
  struct cell *cells = e->s->cells_top;
  const int nr_cells = e->s->nr_cells;
  /* const double max_d = */
  /*     const_gravity_a_smooth * const_gravity_r_cut * ci->width[0]; */
  /* const double max_d2 = max_d * max_d; */
  // const double pos_i[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};

  /* Anything to do here? */
  if (!cell_is_active(ci, e)) return;

  /* Drift our own multipole if need be */
  if (ci->ti_old_multipole != e->ti_current) cell_drift_multipole(ci, e);

  /* Loop over all the cells and go for a p-m interaction if far enough but not
   * too far */
  for (int i = 0; i < nr_cells; ++i) {

    struct cell *cj = &cells[i];

    if (ci == cj) continue;
    if (cj->gcount == 0) continue;

    /* const double dx[3] = {cj->loc[0] - pos_i[0],   // x */
    /*                       cj->loc[1] - pos_i[1],   // y */
    /*                       cj->loc[2] - pos_i[2]};  // z */
    /* const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; */
    /* if (r2 > max_d2) continue; */

    if (!cell_are_neighbours(ci, cj)) runner_dopair_grav_mm(r, ci, cj);
  }

  if (timer) TIMER_TOC(timer_dograv_long_range);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
