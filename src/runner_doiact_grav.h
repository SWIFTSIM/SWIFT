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
#include "inline.h"
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

  /* Some constants */
  const struct engine *e = r->e;

  /* Cell properties */
  struct gpart *gparts = c->gparts;
  const int gcount = c->gcount;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->ti_old_multipole != e->ti_current) error("c->multipole not drifted.");
#endif

  if (c->split) { /* Node case */

    /* Add the field-tensor to all the 8 progenitors */
    for (int k = 0; k < 8; ++k) {
      struct cell *cp = c->progeny[k];

      /* Do we have a progenitor with any active g-particles ? */
      if (cp != NULL && cell_is_active(cp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        if (cp->ti_old_multipole != e->ti_current)
          error("cp->multipole not drifted.");
#endif
        struct grav_tensor shifted_tensor;

        /* Shift the field tensor */
        gravity_L2L(&shifted_tensor, &c->multipole->pot, cp->multipole->CoM,
                    c->multipole->CoM);

        /* Add it to this level's tensor */
        gravity_field_tensors_add(&cp->multipole->pot, &shifted_tensor);

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
      if (gpart_is_active(gp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (gp->ti_drift != e->ti_current)
          error("gpart not drifted to current time");
#endif

        /* Apply the kernel */
        gravity_L2P(&c->multipole->pot, c->multipole->CoM, gp);
      }
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

  /* Some constants */
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const struct gravity_props *props = e->gravity_properties;
  // const float a_smooth = e->gravity_properties->a_smooth;
  // const float rlr_inv = 1. / (a_smooth * ci->super->width[0]);

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e)) return;

  /* Short-cut to the multipole */
  const struct multipole *multi_j = &cj->multipole->m_pole;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_j->M_000 == 0.f) error("Multipole does not seem to have been set.");

  if (ci->ti_old_multipole != e->ti_current)
    error("ci->multipole not drifted.");
#endif

  /* Do we need to drift the multipole ? */
  if (cj->ti_old_multipole != e->ti_current) cell_drift_multipole(cj, e);

  /* Let's interact at this level */
  gravity_M2L(&ci->multipole->pot, multi_j, ci->multipole->CoM,
              cj->multipole->CoM, props, periodic, dim);

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell using the full Newtonian potential
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param shift The distance vector (periodically wrapped) between the cell
 * centres.
 */
void runner_dopair_grav_pp_full(struct runner *r, struct cell *ci,
                                struct cell *cj, double shift[3]) {

  /* Some constants */
  const struct engine *const e = r->e;
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;
  struct gravity_cache *const cj_cache = &r->cj_gravity_cache;

  /* Cell properties */
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;
  const int ci_active = cell_is_active(ci, e);
  const int cj_active = cell_is_active(cj, e);
  const double loc_i[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const double loc_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};
  const double loc_mean[3] = {0.5 * (loc_i[0] + loc_j[0]),
                              0.5 * (loc_i[1] + loc_j[1]),
                              0.5 * (loc_i[2] + loc_j[2])};

  /* Anything to do here ?*/
  if (!ci_active && !cj_active) return;

  /* Check that we fit in cache */
  if (gcount_i > ci_cache->count || gcount_j > cj_cache->count)
    error("Not enough space in the caches! gcount_i=%d gcount_j=%d", gcount_i,
          gcount_j);

  /* Computed the padded counts */
  const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;
  const int gcount_padded_j = gcount_j - (gcount_j % VEC_SIZE) + VEC_SIZE;

  /* Fill the caches */
  gravity_cache_populate(ci_cache, gparts_i, gcount_i, gcount_padded_i,
                         loc_mean);
  gravity_cache_populate(cj_cache, gparts_j, gcount_j, gcount_padded_j,
                         loc_mean);

  /* Ok... Here we go ! */

  if (ci_active) {

    /* Loop over all particles in ci... */
    for (int pid = 0; pid < gcount_i; pid++) {

      /* Skip inactive particles */
      if (!gpart_is_active(&gparts_i[pid], e)) continue;

      const float x_i = ci_cache->x[pid];
      const float y_i = ci_cache->y[pid];
      const float z_i = ci_cache->z[pid];

      /* Some powers of the softening length */
      const float h_i = ci_cache->epsilon[pid];
      const float h2_i = h_i * h_i;
      const float h_inv_i = 1.f / h_i;
      const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

      /* Local accumulators for the acceleration */
      float a_x = 0.f, a_y = 0.f, a_z = 0.f;

      /* Make the compiler understand we are in happy vectorization land */
      swift_align_information(cj_cache->x, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->y, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->z, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->m, SWIFT_CACHE_ALIGNMENT);
      swift_assume_size(gcount_padded_j, VEC_SIZE);

      /* Loop over every particle in the other cell. */
      for (int pjd = 0; pjd < gcount_padded_j; pjd++) {

        /* Get info about j */
        const float x_j = cj_cache->x[pjd];
        const float y_j = cj_cache->y[pjd];
        const float z_j = cj_cache->z[pjd];
        const float mass_j = cj_cache->m[pjd];

        /* Compute the pairwise (square) distance. */
        const float dx = x_i - x_j;
        const float dy = y_i - y_j;
        const float dz = z_i - z_j;
        const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
        if (r2 == 0.f) error("Interacting particles with 0 distance");

        /* Check that particles have been drifted to the current time */
        if (gparts_i[pid].ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (pjd < gcount_j && gparts_j[pjd].ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Get the inverse distance */
        const float r_inv = 1.f / sqrtf(r2);

        float f_ij, W_ij;

        if (r2 >= h2_i) {

          /* Get Newtonian gravity */
          f_ij = mass_j * r_inv * r_inv * r_inv;

        } else {

          const float r = r2 * r_inv;
          const float ui = r * h_inv_i;

          kernel_grav_eval(ui, &W_ij);

          /* Get softened gravity */
          f_ij = mass_j * h_inv3_i * W_ij;
        }

        /* Store it back */
        a_x -= f_ij * dx;
        a_y -= f_ij * dy;
        a_z -= f_ij * dz;

#ifdef SWIFT_DEBUG_CHECKS
        /* Update the interaction counter if it's not a padded gpart */
        if (pjd < gcount_j) gparts_i[pid].num_interacted++;
#endif
      }

      /* Store everything back in cache */
      ci_cache->a_x[pid] = a_x;
      ci_cache->a_y[pid] = a_y;
      ci_cache->a_z[pid] = a_z;
    }
  }

  /* Now do the opposite loop */
  if (cj_active) {

    /* Loop over all particles in ci... */
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Skip inactive particles */
      if (!gpart_is_active(&gparts_j[pjd], e)) continue;

      const float x_j = cj_cache->x[pjd];
      const float y_j = cj_cache->y[pjd];
      const float z_j = cj_cache->z[pjd];

      /* Some powers of the softening length */
      const float h_j = cj_cache->epsilon[pjd];
      const float h2_j = h_j * h_j;
      const float h_inv_j = 1.f / h_j;
      const float h_inv3_j = h_inv_j * h_inv_j * h_inv_j;

      /* Local accumulators for the acceleration */
      float a_x = 0.f, a_y = 0.f, a_z = 0.f;

      /* Make the compiler understand we are in happy vectorization land */
      swift_align_information(ci_cache->x, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->y, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->z, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->m, SWIFT_CACHE_ALIGNMENT);
      swift_assume_size(gcount_padded_i, VEC_SIZE);

      /* Loop over every particle in the other cell. */
      for (int pid = 0; pid < gcount_padded_i; pid++) {

        /* Get info about j */
        const float x_i = ci_cache->x[pid];
        const float y_i = ci_cache->y[pid];
        const float z_i = ci_cache->z[pid];
        const float mass_i = ci_cache->m[pid];

        /* Compute the pairwise (square) distance. */
        const float dx = x_j - x_i;
        const float dy = y_j - y_i;
        const float dz = z_j - z_i;
        const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
        if (r2 == 0.f) error("Interacting particles with 0 distance");

        /* Check that particles have been drifted to the current time */
        if (gparts_j[pjd].ti_drift != e->ti_current)
          error("gpj not drifted to current time");
        if (pid < gcount_i && gparts_i[pid].ti_drift != e->ti_current)
          error("gpi not drifted to current time");
#endif

        /* Get the inverse distance */
        const float r_inv = 1.f / sqrtf(r2);

        float f_ji, W_ji;

        if (r2 >= h2_j) {

          /* Get Newtonian gravity */
          f_ji = mass_i * r_inv * r_inv * r_inv;

        } else {

          const float r = r2 * r_inv;
          const float uj = r * h_inv_j;

          kernel_grav_eval(uj, &W_ji);

          /* Get softened gravity */
          f_ji = mass_i * h_inv3_j * W_ji;
        }

        /* Store it back */
        a_x -= f_ji * dx;
        a_y -= f_ji * dy;
        a_z -= f_ji * dz;

#ifdef SWIFT_DEBUG_CHECKS
        /* Update the interaction counter if it's not a padded gpart */
        if (pid < gcount_i) gparts_j[pjd].num_interacted++;
#endif
      }

      /* Store everything back in cache */
      cj_cache->a_x[pjd] = a_x;
      cj_cache->a_y[pjd] = a_y;
      cj_cache->a_z[pjd] = a_z;
    }
  }

  /* Write back to the particles */
  if (ci_active) gravity_cache_write_back(ci_cache, gparts_i, gcount_i);
  if (cj_active) gravity_cache_write_back(cj_cache, gparts_j, gcount_j);

#ifdef MATTHIEU_OLD_STUFF

  /* Some constants */
  const struct engine *const e = r->e;

  /* Cell properties */
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;

  /* MATTHIEU: Should we use local DP accumulators ? */

  /* Loop over all particles in ci... */
  if (cell_is_active(ci, e)) {
    for (int pid = 0; pid < gcount_i; pid++) {

      /* Get a hold of the ith part in ci. */
      struct gpart *restrict gpi = &gparts_i[pid];

      if (!gpart_is_active(gpi, e)) continue;

      /* Apply boundary condition */
      const double pix[3] = {gpi->x[0] - shift[0], gpi->x[1] - shift[1],
                             gpi->x[2] - shift[2]};

      /* Loop over every particle in the other cell. */
      for (int pjd = 0; pjd < gcount_j; pjd++) {

        /* Get a hold of the jth part in cj. */
        const struct gpart *restrict gpj = &gparts_j[pjd];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix[0] - gpj->x[0],   // x
                             pix[1] - gpj->x[1],   // y
                             pix[2] - gpj->x[2]};  // z
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (gpi->ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (gpj->ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Interact ! */
        runner_iact_grav_pp_nonsym(r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
        gpi->num_interacted++;
#endif
      }
    }
  }

  /* Loop over all particles in cj... */
  if (cell_is_active(cj, e)) {
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Get a hold of the ith part in ci. */
      struct gpart *restrict gpj = &gparts_j[pjd];

      if (!gpart_is_active(gpj, e)) continue;

      /* Apply boundary condition */
      const double pjx[3] = {gpj->x[0] + shift[0], gpj->x[1] + shift[1],
                             gpj->x[2] + shift[2]};

      /* Loop over every particle in the other cell. */
      for (int pid = 0; pid < gcount_i; pid++) {

        /* Get a hold of the ith part in ci. */
        const struct gpart *restrict gpi = &gparts_i[pid];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx[0] - gpi->x[0],   // x
                             pjx[1] - gpi->x[1],   // y
                             pjx[2] - gpi->x[2]};  // z
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (gpi->ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (gpj->ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Interact ! */
        runner_iact_grav_pp_nonsym(r2, dx, gpj, gpi);

#ifdef SWIFT_DEBUG_CHECKS
        gpj->num_interacted++;
#endif
      }
    }
  }
#endif
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell using the truncated Newtonian potential
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param shift The distance vector (periodically wrapped) between the cell
 * centres.
 */
void runner_dopair_grav_pp_truncated(struct runner *r, struct cell *ci,
                                     struct cell *cj, double shift[3]) {

  /* Some constants */
  const struct engine *const e = r->e;
  const struct space *s = e->s;
  const double cell_width = s->width[0];
  const double a_smooth = e->gravity_properties->a_smooth;
  const double rlr = cell_width * a_smooth;
  const float rlr_inv = 1. / rlr;

  /* Caches to play with */
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;
  struct gravity_cache *const cj_cache = &r->cj_gravity_cache;

  /* Cell properties */
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;
  const int ci_active = cell_is_active(ci, e);
  const int cj_active = cell_is_active(cj, e);
  const double loc_i[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const double loc_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};
  const double loc_mean[3] = {0.5 * (loc_i[0] + loc_j[0]),
                              0.5 * (loc_i[1] + loc_j[1]),
                              0.5 * (loc_i[2] + loc_j[2])};

  /* Anything to do here ?*/
  if (!ci_active && !cj_active) return;

  /* Check that we fit in cache */
  if (gcount_i > ci_cache->count || gcount_j > cj_cache->count)
    error("Not enough space in the caches! gcount_i=%d gcount_j=%d", gcount_i,
          gcount_j);

  /* Computed the padded counts */
  const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;
  const int gcount_padded_j = gcount_j - (gcount_j % VEC_SIZE) + VEC_SIZE;

  /* Fill the caches */
  gravity_cache_populate(ci_cache, gparts_i, gcount_i, gcount_padded_i,
                         loc_mean);
  gravity_cache_populate(cj_cache, gparts_j, gcount_j, gcount_padded_j,
                         loc_mean);

  /* Ok... Here we go ! */

  if (ci_active) {

    /* Loop over all particles in ci... */
    for (int pid = 0; pid < gcount_i; pid++) {

      /* Skip inactive particles */
      if (!gpart_is_active(&gparts_i[pid], e)) continue;

      const float x_i = ci_cache->x[pid];
      const float y_i = ci_cache->y[pid];
      const float z_i = ci_cache->z[pid];

      /* Some powers of the softening length */
      const float h_i = ci_cache->epsilon[pid];
      const float h2_i = h_i * h_i;
      const float h_inv_i = 1.f / h_i;
      const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

      /* Local accumulators for the acceleration */
      float a_x = 0.f, a_y = 0.f, a_z = 0.f;

      /* Make the compiler understand we are in happy vectorization land */
      swift_align_information(cj_cache->x, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->y, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->z, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(cj_cache->m, SWIFT_CACHE_ALIGNMENT);
      swift_assume_size(gcount_padded_j, VEC_SIZE);

      /* Loop over every particle in the other cell. */
      for (int pjd = 0; pjd < gcount_padded_j; pjd++) {

        /* Get info about j */
        const float x_j = cj_cache->x[pjd];
        const float y_j = cj_cache->y[pjd];
        const float z_j = cj_cache->z[pjd];
        const float mass_j = cj_cache->m[pjd];

        /* Compute the pairwise (square) distance. */
        const float dx = x_i - x_j;
        const float dy = y_i - y_j;
        const float dz = z_i - z_j;
        const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
        if (r2 == 0.f) error("Interacting particles with 0 distance");

        /* Check that particles have been drifted to the current time */
        if (gparts_i[pid].ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (pjd < gcount_j && gparts_j[pjd].ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Get the inverse distance */
        const float r_inv = 1.f / sqrtf(r2);
        const float r = r2 * r_inv;

        float f_ij, W_ij, corr_lr;

        if (r2 >= h2_i) {

          /* Get Newtonian gravity */
          f_ij = mass_j * r_inv * r_inv * r_inv;

        } else {

          const float ui = r * h_inv_i;

          kernel_grav_eval(ui, &W_ij);

          /* Get softened gravity */
          f_ij = mass_j * h_inv3_i * W_ij;
        }

        /* Get long-range correction */
        const float u_lr = r * rlr_inv;
        kernel_long_grav_eval(u_lr, &corr_lr);
        f_ij *= corr_lr;

        /* Store it back */
        a_x -= f_ij * dx;
        a_y -= f_ij * dy;
        a_z -= f_ij * dz;

#ifdef SWIFT_DEBUG_CHECKS
        /* Update the interaction counter if it's not a padded gpart */
        if (pjd < gcount_j) gparts_i[pid].num_interacted++;
#endif
      }

      /* Store everything back in cache */
      ci_cache->a_x[pid] = a_x;
      ci_cache->a_y[pid] = a_y;
      ci_cache->a_z[pid] = a_z;
    }
  }

  /* Now do the opposite loop */
  if (cj_active) {

    /* Loop over all particles in ci... */
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Skip inactive particles */
      if (!gpart_is_active(&gparts_j[pjd], e)) continue;

      const float x_j = cj_cache->x[pjd];
      const float y_j = cj_cache->y[pjd];
      const float z_j = cj_cache->z[pjd];

      /* Some powers of the softening length */
      const float h_j = cj_cache->epsilon[pjd];
      const float h2_j = h_j * h_j;
      const float h_inv_j = 1.f / h_j;
      const float h_inv3_j = h_inv_j * h_inv_j * h_inv_j;

      /* Local accumulators for the acceleration */
      float a_x = 0.f, a_y = 0.f, a_z = 0.f;

      /* Make the compiler understand we are in happy vectorization land */
      swift_align_information(ci_cache->x, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->y, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->z, SWIFT_CACHE_ALIGNMENT);
      swift_align_information(ci_cache->m, SWIFT_CACHE_ALIGNMENT);
      swift_assume_size(gcount_padded_i, VEC_SIZE);

      /* Loop over every particle in the other cell. */
      for (int pid = 0; pid < gcount_padded_i; pid++) {

        /* Get info about j */
        const float x_i = ci_cache->x[pid];
        const float y_i = ci_cache->y[pid];
        const float z_i = ci_cache->z[pid];
        const float mass_i = ci_cache->m[pid];

        /* Compute the pairwise (square) distance. */
        const float dx = x_j - x_i;
        const float dy = y_j - y_i;
        const float dz = z_j - z_i;
        const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
        if (r2 == 0.f) error("Interacting particles with 0 distance");

        /* Check that particles have been drifted to the current time */
        if (gparts_j[pjd].ti_drift != e->ti_current)
          error("gpj not drifted to current time");
        if (pid < gcount_i && gparts_i[pid].ti_drift != e->ti_current)
          error("gpi not drifted to current time");
#endif

        /* Get the inverse distance */
        const float r_inv = 1.f / sqrtf(r2);
        const float r = r2 * r_inv;

        float f_ji, W_ji, corr_lr;

        if (r2 >= h2_j) {

          /* Get Newtonian gravity */
          f_ji = mass_i * r_inv * r_inv * r_inv;

        } else {

          const float uj = r * h_inv_j;

          kernel_grav_eval(uj, &W_ji);

          /* Get softened gravity */
          f_ji = mass_i * h_inv3_j * W_ji;
        }

        /* Get long-range correction */
        const float u_lr = r * rlr_inv;
        kernel_long_grav_eval(u_lr, &corr_lr);
        f_ji *= corr_lr;

        /* Store it back */
        a_x -= f_ji * dx;
        a_y -= f_ji * dy;
        a_z -= f_ji * dz;

#ifdef SWIFT_DEBUG_CHECKS
        /* Update the interaction counter if it's not a padded gpart */
        if (pid < gcount_i) gparts_j[pjd].num_interacted++;
#endif
      }

      /* Store everything back in cache */
      cj_cache->a_x[pjd] = a_x;
      cj_cache->a_y[pjd] = a_y;
      cj_cache->a_z[pjd] = a_z;
    }
  }

  /* Write back to the particles */
  if (ci_active) gravity_cache_write_back(ci_cache, gparts_i, gcount_i);
  if (cj_active) gravity_cache_write_back(cj_cache, gparts_j, gcount_j);

#ifdef MATTHIEU_OLD_STUFF
  /* Some constants */
  const struct engine *const e = r->e;
  const struct space *s = e->s;
  const double cell_width = s->width[0];
  const double a_smooth = e->gravity_properties->a_smooth;
  const double rlr = cell_width * a_smooth;
  const float rlr_inv = 1. / rlr;

  /* Cell properties */
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;

  /* MATTHIEU: Should we use local DP accumulators ? */

  /* Loop over all particles in ci... */
  if (cell_is_active(ci, e)) {
    for (int pid = 0; pid < gcount_i; pid++) {

      /* Get a hold of the ith part in ci. */
      struct gpart *restrict gpi = &gparts_i[pid];

      if (!gpart_is_active(gpi, e)) continue;

      /* Apply boundary condition */
      const double pix[3] = {gpi->x[0] - shift[0], gpi->x[1] - shift[1],
                             gpi->x[2] - shift[2]};

      /* Loop over every particle in the other cell. */
      for (int pjd = 0; pjd < gcount_j; pjd++) {

        /* Get a hold of the jth part in cj. */
        const struct gpart *restrict gpj = &gparts_j[pjd];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix[0] - gpj->x[0],   // x
                             pix[1] - gpj->x[1],   // y
                             pix[2] - gpj->x[2]};  // z
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (gpi->ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (gpj->ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Interact ! */
        runner_iact_grav_pp_truncated_nonsym(r2, dx, gpi, gpj, rlr_inv);

#ifdef SWIFT_DEBUG_CHECKS
        gpi->num_interacted++;
#endif
      }
    }
  }

  /* Loop over all particles in cj... */
  if (cell_is_active(cj, e)) {
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Get a hold of the ith part in ci. */
      struct gpart *restrict gpj = &gparts_j[pjd];

      if (!gpart_is_active(gpj, e)) continue;

      /* Apply boundary condition */
      const double pjx[3] = {gpj->x[0] + shift[0], gpj->x[1] + shift[1],
                             gpj->x[2] + shift[2]};

      /* Loop over every particle in the other cell. */
      for (int pid = 0; pid < gcount_i; pid++) {

        /* Get a hold of the ith part in ci. */
        const struct gpart *restrict gpi = &gparts_i[pid];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx[0] - gpi->x[0],   // x
                             pjx[1] - gpi->x[1],   // y
                             pjx[2] - gpi->x[2]};  // z
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (gpi->ti_drift != e->ti_current)
          error("gpi not drifted to current time");
        if (gpj->ti_drift != e->ti_current)
          error("gpj not drifted to current time");
#endif

        /* Interact ! */
        runner_iact_grav_pp_truncated_nonsym(r2, dx, gpj, gpi, rlr_inv);

#ifdef SWIFT_DEBUG_CHECKS
        gpj->num_interacted++;
#endif
      }
    }
  }

#endif
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell (switching function between full and truncated).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 */
void runner_dopair_grav_pp(struct runner *r, struct cell *ci, struct cell *cj) {

  /* Some properties of the space */
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const int periodic = s->periodic;
  const double cell_width = s->width[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double a_smooth = e->gravity_properties->a_smooth;
  const double r_cut_min = e->gravity_properties->r_cut_min;
  const double min_trunc = cell_width * r_cut_min * a_smooth;
  double shift[3] = {0.0, 0.0, 0.0};

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  /* Let's start by drifting things */
  if (!cell_are_gpart_drifted(ci, e)) cell_drift_gpart(ci, e);
  if (!cell_are_gpart_drifted(cj, e)) cell_drift_gpart(cj, e);

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {
    runner_dopair_grav_pp_full(r, ci, cj, shift);
  } else {

    /* Get the relative distance between the pairs, wrapping. */
    shift[0] = nearest(cj->loc[0] - ci->loc[0], dim[0]);
    shift[1] = nearest(cj->loc[1] - ci->loc[1], dim[1]);
    shift[2] = nearest(cj->loc[2] - ci->loc[2], dim[2]);
    const double r2 =
        shift[0] * shift[0] + shift[1] * shift[1] + shift[2] * shift[2];

    /* Get the maximal distance between any two particles */
    const double max_r = sqrt(r2) + ci->multipole->r_max + cj->multipole->r_max;

    /* Do we need to use the truncated interactions ? */
    if (max_r > min_trunc)
      runner_dopair_grav_pp_truncated(r, ci, cj, shift);
    else
      runner_dopair_grav_pp_full(r, ci, cj, shift);
  }

  TIMER_TOC(timer_dopair_grav_pp);
}

/**
 * @brief Computes the interaction of all the particles in a cell using the
 * full Newtonian potential.
 *
 * @param r The #runner.
 * @param c The #cell.
 *
 * @todo Use a local cache for the particles.
 */
void runner_doself_grav_pp_full(struct runner *r, struct cell *c) {

  /* Some constants */
  const struct engine *const e = r->e;
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;

  /* Cell properties */
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;
  const int c_active = cell_is_active(c, e);
  const double loc[3] = {c->loc[0] + 0.5 * c->width[0],
                         c->loc[1] + 0.5 * c->width[1],
                         c->loc[2] + 0.5 * c->width[2]};

  /* Anything to do here ?*/
  if (!c_active) return;

  /* Check that we fit in cache */
  if (gcount > ci_cache->count)
    error("Not enough space in the cache! gcount=%d", gcount);

  /* Computed the padded counts */
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

  gravity_cache_populate(ci_cache, gparts, gcount, gcount_padded, loc);

  /* Ok... Here we go ! */

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Skip inactive particles */
    if (!gpart_is_active(&gparts[pid], e)) continue;

    const float x_i = ci_cache->x[pid];
    const float y_i = ci_cache->y[pid];
    const float z_i = ci_cache->z[pid];

    /* Some powers of the softening length */
    const float h_i = ci_cache->epsilon[pid];
    const float h2_i = h_i * h_i;
    const float h_inv_i = 1.f / h_i;
    const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

    /* Local accumulators for the acceleration */
    float a_x = 0.f, a_y = 0.f, a_z = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->m, SWIFT_CACHE_ALIGNMENT);
    swift_assume_size(gcount_padded, VEC_SIZE);

    /* Loop over every other particle in the cell. */
    for (int pjd = 0; pjd < gcount_padded; pjd++) {

      /* No self interaction */
      if (pid == pjd) continue;

      /* Get info about j */
      const float x_j = ci_cache->x[pjd];
      const float y_j = ci_cache->y[pjd];
      const float z_j = ci_cache->z[pjd];
      const float mass_j = ci_cache->m[pjd];

      /* Compute the pairwise (square) distance. */
      const float dx = x_i - x_j;
      const float dy = y_i - y_j;
      const float dz = z_i - z_j;
      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f) error("Interacting particles with 0 distance");

      /* Check that particles have been drifted to the current time */
      if (gparts[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount && gparts[pjd].ti_drift != e->ti_current)
        error("gpj not drifted to current time");
#endif

      /* Get the inverse distance */
      const float r_inv = 1.f / sqrtf(r2);

      float f_ij, W_ij;

      if (r2 >= h2_i) {

        /* Get Newtonian gravity */
        f_ij = mass_j * r_inv * r_inv * r_inv;

      } else {

        const float r = r2 * r_inv;
        const float ui = r * h_inv_i;

        kernel_grav_eval(ui, &W_ij);

        /* Get softened gravity */
        f_ij = mass_j * h_inv3_i * W_ij;
      }

      /* Store it back */
      a_x -= f_ij * dx;
      a_y -= f_ij * dy;
      a_z -= f_ij * dz;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount) gparts[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] = a_x;
    ci_cache->a_y[pid] = a_y;
    ci_cache->a_z[pid] = a_z;
  }

  /* Write back to the particles */
  gravity_cache_write_back(ci_cache, gparts, gcount);

#ifdef MATTHIEU_OLD_STUFF

  /* Some constants */
  const struct engine *const e = r->e;

  /* Cell properties */
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;

  /* MATTHIEU: Should we use local DP accumulators ? */

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

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (gpi->ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (gpj->ti_drift != e->ti_current)
        error("gpj not drifted to current time");
#endif

      /* Interact ! */
      if (gpart_is_active(gpi, e) && gpart_is_active(gpj, e)) {

        runner_iact_grav_pp(r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
        gpi->num_interacted++;
        gpj->num_interacted++;
#endif

      } else {

        if (gpart_is_active(gpi, e)) {

          runner_iact_grav_pp_nonsym(r2, dx, gpi, gpj);

#ifdef SWIFT_DEBUG_CHECKS
          gpi->num_interacted++;
#endif

        } else if (gpart_is_active(gpj, e)) {

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          runner_iact_grav_pp_nonsym(r2, dx, gpj, gpi);

#ifdef SWIFT_DEBUG_CHECKS
          gpj->num_interacted++;
#endif
        }
      }
    }
  }

#endif
}

/**
 * @brief Computes the interaction of all the particles in a cell using the
 * truncated Newtonian potential.
 *
 * @param r The #runner.
 * @param c The #cell.
 *
 * @todo Use a local cache for the particles.
 */
void runner_doself_grav_pp_truncated(struct runner *r, struct cell *c) {

  /* Some constants */
  const struct engine *const e = r->e;
  const struct space *s = e->s;
  const double cell_width = s->width[0];
  const double a_smooth = e->gravity_properties->a_smooth;
  const double rlr = cell_width * a_smooth;
  const float rlr_inv = 1. / rlr;

  /* Caches to play with */
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;

  /* Cell properties */
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;
  const int c_active = cell_is_active(c, e);
  const double loc[3] = {c->loc[0] + 0.5 * c->width[0],
                         c->loc[1] + 0.5 * c->width[1],
                         c->loc[2] + 0.5 * c->width[2]};

  /* Anything to do here ?*/
  if (!c_active) return;

  /* Check that we fit in cache */
  if (gcount > ci_cache->count)
    error("Not enough space in the caches! gcount=%d", gcount);

  /* Computed the padded counts */
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

  gravity_cache_populate(ci_cache, gparts, gcount, gcount_padded, loc);

  /* Ok... Here we go ! */

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Skip inactive particles */
    if (!gpart_is_active(&gparts[pid], e)) continue;

    const float x_i = ci_cache->x[pid];
    const float y_i = ci_cache->y[pid];
    const float z_i = ci_cache->z[pid];

    /* Some powers of the softening length */
    const float h_i = ci_cache->epsilon[pid];
    const float h2_i = h_i * h_i;
    const float h_inv_i = 1.f / h_i;
    const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

    /* Local accumulators for the acceleration */
    float a_x = 0.f, a_y = 0.f, a_z = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(ci_cache->m, SWIFT_CACHE_ALIGNMENT);
    swift_assume_size(gcount_padded, VEC_SIZE);

    /* Loop over every other particle in the cell. */
    for (int pjd = 0; pjd < gcount_padded; pjd++) {

      /* No self interaction */
      if (pid == pjd) continue;

      /* Get info about j */
      const float x_j = ci_cache->x[pjd];
      const float y_j = ci_cache->y[pjd];
      const float z_j = ci_cache->z[pjd];
      const float mass_j = ci_cache->m[pjd];

      /* Compute the pairwise (square) distance. */
      const float dx = x_i - x_j;
      const float dy = y_i - y_j;
      const float dz = z_i - z_j;
      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f) error("Interacting particles with 0 distance");

      /* Check that particles have been drifted to the current time */
      if (gparts[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount && gparts[pjd].ti_drift != e->ti_current)
        error("gpj not drifted to current time");
#endif

      /* Get the inverse distance */
      const float r_inv = 1.f / sqrtf(r2);
      const float r = r2 * r_inv;

      float f_ij, W_ij, corr_lr;

      if (r2 >= h2_i) {

        /* Get Newtonian gravity */
        f_ij = mass_j * r_inv * r_inv * r_inv;

      } else {

        const float ui = r * h_inv_i;

        kernel_grav_eval(ui, &W_ij);

        /* Get softened gravity */
        f_ij = mass_j * h_inv3_i * W_ij;
      }

      /* Get long-range correction */
      const float u_lr = r * rlr_inv;
      kernel_long_grav_eval(u_lr, &corr_lr);
      f_ij *= corr_lr;

      /* Store it back */
      a_x -= f_ij * dx;
      a_y -= f_ij * dy;
      a_z -= f_ij * dz;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount) gparts[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] = a_x;
    ci_cache->a_y[pid] = a_y;
    ci_cache->a_z[pid] = a_z;
  }

  /* Write back to the particles */
  gravity_cache_write_back(ci_cache, gparts, gcount);

#ifdef MATTHIEU_OLD_STUFF
  /* Some constants */
  const struct engine *const e = r->e;
  const struct space *s = e->s;
  const double cell_width = s->width[0];
  const double a_smooth = e->gravity_properties->a_smooth;
  const double rlr = cell_width * a_smooth;
  const float rlr_inv = 1. / rlr;

  /* Cell properties */
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;

  /* MATTHIEU: Should we use local DP accumulators ? */

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

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (gpi->ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (gpj->ti_drift != e->ti_current)
        error("gpj not drifted to current time");
#endif

      /* Interact ! */
      if (gpart_is_active(gpi, e) && gpart_is_active(gpj, e)) {

        runner_iact_grav_pp_truncated(r2, dx, gpi, gpj, rlr_inv);

#ifdef SWIFT_DEBUG_CHECKS
        gpi->num_interacted++;
        gpj->num_interacted++;
#endif

      } else {

        if (gpart_is_active(gpi, e)) {

          runner_iact_grav_pp_truncated_nonsym(r2, dx, gpi, gpj, rlr_inv);

#ifdef SWIFT_DEBUG_CHECKS
          gpi->num_interacted++;
#endif

        } else if (gpart_is_active(gpj, e)) {

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          runner_iact_grav_pp_truncated_nonsym(r2, dx, gpj, gpi, rlr_inv);

#ifdef SWIFT_DEBUG_CHECKS
          gpj->num_interacted++;
#endif
        }
      }
    }
  }
#endif
}

/**
 * @brief Computes the interaction of all the particles in a cell directly
 * (Switching function between truncated and full)
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void runner_doself_grav_pp(struct runner *r, struct cell *c) {

  /* Some properties of the space */
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const int periodic = s->periodic;
  const double cell_width = s->width[0];
  const double a_smooth = e->gravity_properties->a_smooth;
  const double r_cut_min = e->gravity_properties->r_cut_min;
  const double min_trunc = cell_width * r_cut_min * a_smooth;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->gcount == 0) error("Doing self gravity on an empty cell !");
#endif

  /* Anything to do here? */
  if (!cell_is_active(c, e)) return;

  /* Do we need to start by drifting things ? */
  if (!cell_are_gpart_drifted(c, e)) cell_drift_gpart(c, e);

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {
    runner_doself_grav_pp_full(r, c);
  } else {

    /* Get the maximal distance between any two particles */
    const double max_r = 2 * c->multipole->r_max;

    /* Do we need to use the truncated interactions ? */
    if (max_r > min_trunc)
      runner_doself_grav_pp_truncated(r, c);
    else
      runner_doself_grav_pp_full(r, c);
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

  /* Some constants */
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const int periodic = s->periodic;
  const double cell_width = s->width[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const struct gravity_props *props = e->gravity_properties;
  const double theta_crit_inv = props->theta_crit_inv;
  const double max_distance = props->a_smooth * props->r_cut_max * cell_width;
  const double max_distance2 = max_distance * max_distance;

#ifdef SWIFT_DEBUG_CHECKS

  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;

  /* Early abort? */
  if (gcount_i == 0 || gcount_j == 0)
    error("Doing pair gravity on an empty cell !");

  /* Sanity check */
  if (ci == cj) error("Pair interaction between a cell and itself.");

  if (cell_is_active(ci, e) && ci->ti_old_multipole != e->ti_current)
    error("ci->multipole not drifted.");
  if (cell_is_active(cj, e) && cj->ti_old_multipole != e->ti_current)
    error("cj->multipole not drifted.");
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  /* Recover the multipole information */
  struct gravity_tensors *const multi_i = ci->multipole;
  struct gravity_tensors *const multi_j = cj->multipole;

  /* Get the distance between the CoMs */
  double dx = multi_i->CoM[0] - multi_j->CoM[0];
  double dy = multi_i->CoM[1] - multi_j->CoM[1];
  double dz = multi_i->CoM[2] - multi_j->CoM[2];

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  /* Are we beyond the distance where the truncated forces are 0? */
  if (periodic && r2 > max_distance2) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Need to account for the interactions we missed */
    if (cell_is_active(ci, e))
      multi_i->pot.num_interacted += multi_j->m_pole.num_gpart;
    if (cell_is_active(cj, e))
      multi_j->pot.num_interacted += multi_i->m_pole.num_gpart;
#endif
    return;
  }

  /* OK, we actually need to compute this pair. Let's find the cheapest
   * option... */

  /* Can we use M-M interactions ? */
  if (gravity_multipole_accept(multi_i, multi_j, theta_crit_inv, r2)) {

    /* MATTHIEU: make a symmetric M-M interaction function ! */
    runner_dopair_grav_mm(r, ci, cj);
    runner_dopair_grav_mm(r, cj, ci);
  }
  /* We have two leaves. Go P-P. */
  else if (!ci->split && !cj->split) {
    runner_dopair_grav_pp(r, ci, cj);
  }
  /* Alright, we'll have to split and recurse. */
  else {

    const double ri_max = multi_i->r_max;
    const double rj_max = multi_j->r_max;

    /* Split the larger of the two cells and start over again */
    if (ri_max > rj_max) {

      /* Can we actually split that interaction ? */
      if (ci->split) {

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL)
            runner_dopair_grav(r, ci->progeny[k], cj, 0);
        }

      } else if (cj->split) {
        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL)
            runner_dopair_grav(r, ci, cj->progeny[k], 0);
        }

      } else {
        error("Fundamental error in the logic");
      }
    } else {

      /* Can we actually split that interaction ? */
      if (cj->split) {

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL)
            runner_dopair_grav(r, ci, cj->progeny[k], 0);
        }

      } else if (ci->split) {
        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL)
            runner_dopair_grav(r, ci->progeny[k], cj, 0);
        }

      } else {
        error("Fundamental error in the logic");
      }
    }
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

  /* Some constants */
  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  /* Early abort? */
  if (c->gcount == 0) error("Doing self gravity on an empty cell !");
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(c, e)) return;

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

    runner_dopair_grav(r, ci, cj, 1);
  }
}

/**
 * @brief Performs all M-M interactions between a given top-level cell and all
 * the other top-levels that are far enough.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param timer Are we timing this ?
 */
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

  /* Some constants */
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const struct gravity_props *props = e->gravity_properties;
  const int periodic = s->periodic;
  const double cell_width = s->width[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double theta_crit_inv = props->theta_crit_inv;
  const double max_distance = props->a_smooth * props->r_cut_max * cell_width;
  const double max_distance2 = max_distance * max_distance;

  TIMER_TIC;

  /* Recover the list of top-level cells */
  struct cell *cells = e->s->cells_top;
  const int nr_cells = e->s->nr_cells;

  /* Anything to do here? */
  if (!cell_is_active(ci, e)) return;

  /* Check multipole has been drifted */
  if (ci->ti_old_multipole != e->ti_current)
    error("Interacting un-drifted multipole");

  /* Recover the local multipole */
  struct gravity_tensors *const multi_i = ci->multipole;
  const double CoM_i[3] = {multi_i->CoM[0], multi_i->CoM[1], multi_i->CoM[2]};
  const double CoM_rebuild_i[3] = {multi_i->CoM_rebuild[0],
                                   multi_i->CoM_rebuild[1],
                                   multi_i->CoM_rebuild[2]};

  /* Loop over all the top-level cells and go for a M-M interaction if
   * well-separated */
  for (int i = 0; i < nr_cells; ++i) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[i];
    const struct gravity_tensors *const multi_j = cj->multipole;

    /* Avoid stupid cases */
    if (ci == cj || cj->gcount == 0) continue;

    /* Get the distance between the CoMs */
    double dx = CoM_i[0] - multi_j->CoM[0];
    double dy = CoM_i[1] - multi_j->CoM[1];
    double dz = CoM_i[2] - multi_j->CoM[2];

    /* Apply BC */
    if (periodic) {
      dx = nearest(dx, dim[0]);
      dy = nearest(dy, dim[1]);
      dz = nearest(dz, dim[2]);
    }
    const double r2 = dx * dx + dy * dy + dz * dz;

    /* Are we beyond the distance where the truncated forces are 0 ?*/
    if (periodic && r2 > max_distance2) {

#ifdef SWIFT_DEBUG_CHECKS
      /* Need to account for the interactions we missed */
      multi_i->pot.num_interacted += multi_j->m_pole.num_gpart;
#endif
      continue;
    }

    /* Check the multipole acceptance criterion */
    if (gravity_multipole_accept(multi_i, multi_j, theta_crit_inv, r2)) {

      /* Go for a (non-symmetric) M-M calculation */
      runner_dopair_grav_mm(r, ci, cj);

    } else {

      /* Let's check whether we need to still operate on this pair */

      /* Get the distance between the CoMs at the last rebuild*/
      double dx = CoM_rebuild_i[0] - multi_j->CoM_rebuild[0];
      double dy = CoM_rebuild_i[1] - multi_j->CoM_rebuild[1];
      double dz = CoM_rebuild_i[2] - multi_j->CoM_rebuild[2];

      /* Apply BC */
      if (periodic) {
        dx = nearest(dx, dim[0]);
        dy = nearest(dy, dim[1]);
        dz = nearest(dz, dim[2]);
      }
      const double r2_rebuild = dx * dx + dy * dy + dz * dz;

      /* Is the criterion violated now but was OK at the last rebuild ? */
      if (gravity_multipole_accept_rebuild(multi_i, multi_j, theta_crit_inv,
                                           r2_rebuild)) {

        /* Alright, we have to take charge of that pair in a different way. */
        // MATTHIEU: We should actually open the tree-node here and recurse.
        runner_dopair_grav_mm(r, ci, cj);
      }
    }
  } /* Loop over top-level cells */

  if (timer) TIMER_TOC(timer_dograv_long_range);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
