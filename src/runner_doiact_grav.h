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
  if (c->multipole->pot.ti_init != e->ti_current)
    error("c->field tensor not initialised");
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
        if (cp->multipole->pot.ti_init != e->ti_current)
          error("cp->field tensor not initialised");
#endif
        struct grav_tensor shifted_tensor;

        /* If the tensor received any contribution, push it down */
        if (c->multipole->pot.interacted) {

          /* Shift the field tensor */
          gravity_L2L(&shifted_tensor, &c->multipole->pot, cp->multipole->CoM,
                      c->multipole->CoM);

          /* Add it to this level's tensor */
          gravity_field_tensors_add(&cp->multipole->pot, &shifted_tensor);
        }

        /* Recurse */
        runner_do_grav_down(r, cp, 0);
      }
    }

  } else { /* Leaf case */

    /* We can abort early if no interactions via multipole happened */
    if (!c->multipole->pot.interacted) return;

    if (!cell_are_gpart_drifted(c, e)) error("Un-drifted gparts");

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
        if (c->multipole->pot.ti_init != e->ti_current)
          error("c->field tensor not initialised");
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

  if (ci->multipole->pot.ti_init != e->ti_current)
    error("ci->grav tensor not initialised.");
#endif

  /* Do we need to drift the multipole ? */
  if (cj->ti_old_multipole != e->ti_current) cell_drift_multipole(cj, e);

  /* Let's interact at this level */
  gravity_M2L(&ci->multipole->pot, multi_j, ci->multipole->CoM,
              cj->multipole->CoM, props, periodic, dim);

  TIMER_TOC(timer_dopair_grav_mm);
}

static INLINE void runner_dopair_grav_pp_full(const struct engine *e,
                                              struct gravity_cache *ci_cache,
                                              struct gravity_cache *cj_cache,
                                              int gcount_i, int gcount_j,
                                              int gcount_padded_j,
                                              struct gpart *restrict gparts_i,
                                              struct gpart *restrict gparts_j) {

  TIMER_TIC;

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

    /* Skip particle that can use the multipole */
    if (ci_cache->use_mpole[pid]) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (!gpart_is_active(&gparts_i[pid], e))
      error("Active particle went through the cache");
#endif

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

      /* Interact! */
      float f_ij;
      runner_iact_grav_pp_full(r2, h2_i, h_inv_i, h_inv3_i, mass_j, &f_ij);

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

  TIMER_TOC(timer_dopair_grav_pp);
}

static INLINE void runner_dopair_grav_pp_truncated(
    const struct engine *e, const float rlr_inv, struct gravity_cache *ci_cache,
    struct gravity_cache *cj_cache, int gcount_i, int gcount_j,
    int gcount_padded_j, struct gpart *restrict gparts_i,
    struct gpart *restrict gparts_j) {

  TIMER_TIC;

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

    /* Skip particle that can use the multipole */
    if (ci_cache->use_mpole[pid]) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (!gpart_is_active(&gparts_i[pid], e))
      error("Active particle went through the cache");
#endif

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

      /* Interact! */
      float f_ij;
      runner_iact_grav_pp_truncated(r2, h2_i, h_inv_i, h_inv3_i, mass_j,
                                    rlr_inv, &f_ij);

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

  TIMER_TOC(timer_dopair_grav_pp);
}

static INLINE void runner_dopair_grav_pm(
    const struct engine *restrict e, struct gravity_cache *ci_cache,
    int gcount_i, int gcount_padded_i, struct gpart *restrict gparts_i,
    const float CoM_j[3], const struct multipole *restrict multi_j,
    struct cell *restrict cj) {

  TIMER_TIC;

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, epsilon, ci_cache->epsilon,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_x, ci_cache->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, ci_cache->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, ci_cache->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, active, ci_cache->active,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, use_mpole, ci_cache->use_mpole,
                            SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded_i, VEC_SIZE);

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_padded_i; pid++) {

    /* Skip inactive particles */
    if (!active[pid]) continue;

    /* Skip particle that cannot use the multipole */
    if (!use_mpole[pid]) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pid < gcount_i && !gpart_is_active(&gparts_i[pid], e))
      error("Active particle went through the cache");
#endif

    const float x_i = x[pid];
    const float y_i = y[pid];
    const float z_i = z[pid];

    /* Some powers of the softening length */
    const float h_i = epsilon[pid];
    const float h_inv_i = 1.f / h_i;

    /* Distance to the Multipole */
    const float dx = x_i - CoM_j[0];
    const float dy = y_i - CoM_j[1];
    const float dz = z_i - CoM_j[2];
    const float r2 = dx * dx + dy * dy + dz * dz;

    /* Interact! */
    float f_x, f_y, f_z;
    runner_iact_grav_pm(dx, dy, dz, r2, h_i, h_inv_i, multi_j, &f_x, &f_y,
                        &f_z);

    /* Store it back */
    a_x[pid] = f_x;
    a_y[pid] = f_y;
    a_z[pid] = f_z;

#ifdef SWIFT_DEBUG_CHECKS
    /* Update the interaction counter */
    if (pid < gcount_i)
      gparts_i[pid].num_interacted += cj->multipole->m_pole.num_gpart;
#endif
  }

  TIMER_TOC(timer_dopair_grav_pm);
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

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  /* Check that we are not doing something stupid */
  if (ci->split || cj->split) error("Running P-P on splitable cells");

  /* Let's start by drifting things */
  if (!cell_are_gpart_drifted(ci, e)) error("Un-drifted gparts");
  if (!cell_are_gpart_drifted(cj, e)) error("Un-drifted gparts");

  /* Recover some useful constants */
  struct space *s = e->s;
  const int periodic = s->periodic;
  const double cell_width = s->width[0];
  const float theta_crit2 = e->gravity_properties->theta_crit2;
  const double a_smooth = e->gravity_properties->a_smooth;
  const double r_cut_min = e->gravity_properties->r_cut_min;
  const double rlr = cell_width * a_smooth;
  const double min_trunc = rlr * r_cut_min;
  const float rlr_inv = 1. / rlr;

  /* Caches to play with */
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;
  struct gravity_cache *const cj_cache = &r->cj_gravity_cache;

  /* Get the distance vector between the pairs, wrapping. */
  double cell_shift[3];
  space_getsid(s, &ci, &cj, cell_shift);

  /* Record activity status */
  const int ci_active = cell_is_active(ci, e);
  const int cj_active = cell_is_active(cj, e);

  /* Do we need to drift the multipoles ? */
  if (cj_active && ci->ti_old_multipole != e->ti_current)
    cell_drift_multipole(ci, e);
  if (ci_active && cj->ti_old_multipole != e->ti_current)
    cell_drift_multipole(cj, e);

  /* Centre of the cell pair */
  const double loc[3] = {ci->loc[0],   // + 0. * ci->width[0],
                         ci->loc[1],   // + 0. * ci->width[1],
                         ci->loc[2]};  // + 0. * ci->width[2]};

  /* Shift to apply to the particles in each cell */
  const double shift_i[3] = {loc[0] + cell_shift[0], loc[1] + cell_shift[1],
                             loc[2] + cell_shift[2]};
  const double shift_j[3] = {loc[0], loc[1], loc[2]};

  /* Recover the multipole info and shift the CoM locations */
  const float rmax_i = ci->multipole->r_max;
  const float rmax_j = cj->multipole->r_max;
  const float rmax2_i = rmax_i * rmax_i;
  const float rmax2_j = rmax_j * rmax_j;
  const struct multipole *multi_i = &ci->multipole->m_pole;
  const struct multipole *multi_j = &cj->multipole->m_pole;
  const float CoM_i[3] = {ci->multipole->CoM[0] - shift_i[0],
                          ci->multipole->CoM[1] - shift_i[1],
                          ci->multipole->CoM[2] - shift_i[2]};
  const float CoM_j[3] = {cj->multipole->CoM[0] - shift_j[0],
                          cj->multipole->CoM[1] - shift_j[1],
                          cj->multipole->CoM[2] - shift_j[2]};

  /* Start by constructing particle caches */

  /* Computed the padded counts */
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;
  const int gcount_padded_j = gcount_j - (gcount_j % VEC_SIZE) + VEC_SIZE;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we fit in cache */
  if (gcount_i > ci_cache->count || gcount_j > cj_cache->count)
    error("Not enough space in the caches! gcount_i=%d gcount_j=%d", gcount_i,
          gcount_j);
#endif

  /* Fill the caches */
  gravity_cache_populate(e->max_active_bin, ci_cache, ci->gparts, gcount_i,
                         gcount_padded_i, shift_i, CoM_j, rmax2_j, theta_crit2,
                         ci);
  gravity_cache_populate(e->max_active_bin, cj_cache, cj->gparts, gcount_j,
                         gcount_padded_j, shift_j, CoM_i, rmax2_i, theta_crit2,
                         cj);

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {

    /* Not periodic -> Can always use Newtonian potential */

    /* Let's updated the active cell(s) only */
    if (ci_active) {

      /* First the P2P */
      runner_dopair_grav_pp_full(e, ci_cache, cj_cache, gcount_i, gcount_j,
                                 gcount_padded_j, ci->gparts, cj->gparts);

      /* Then the M2P */
      runner_dopair_grav_pm(e, ci_cache, gcount_i, gcount_padded_i, ci->gparts,
                            CoM_j, multi_j, cj);
    }
    if (cj_active) {

      /* First the P2P */
      runner_dopair_grav_pp_full(e, cj_cache, ci_cache, gcount_j, gcount_i,
                                 gcount_padded_i, cj->gparts, ci->gparts);
      /* Then the M2P */
      runner_dopair_grav_pm(e, cj_cache, gcount_j, gcount_padded_j, cj->gparts,
                            CoM_i, multi_i, ci);
    }

  } else { /* Periodic BC */

    /* Get the relative distance between the CoMs */
    const double dx[3] = {CoM_j[0] - CoM_i[0], CoM_j[1] - CoM_i[1],
                          CoM_j[2] - CoM_i[2]};
    const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    /* Get the maximal distance between any two particles */
    const double max_r = sqrt(r2) + rmax_i + rmax_j;

    /* Do we need to use the truncated interactions ? */
    if (max_r > min_trunc) {

      /* Periodic but far-away cells must use the truncated potential */

      /* Let's updated the active cell(s) only */
      if (ci_active) {

        /* First the (truncated) P2P */
        runner_dopair_grav_pp_truncated(e, rlr_inv, ci_cache, cj_cache,
                                        gcount_i, gcount_j, gcount_padded_j,
                                        ci->gparts, cj->gparts);

        /* Then the M2P */
        runner_dopair_grav_pm(e, ci_cache, gcount_i, gcount_padded_i,
                              ci->gparts, CoM_j, multi_j, cj);
      }
      if (cj_active) {

        /* First the (truncated) P2P */
        runner_dopair_grav_pp_truncated(e, rlr_inv, cj_cache, ci_cache,
                                        gcount_j, gcount_i, gcount_padded_i,
                                        cj->gparts, ci->gparts);

        /* Then the M2P */
        runner_dopair_grav_pm(e, cj_cache, gcount_j, gcount_padded_j,
                              cj->gparts, CoM_i, multi_i, ci);
      }

    } else {

      /* Periodic but close-by cells can use the full Newtonian potential */

      /* Let's updated the active cell(s) only */
      if (ci_active) {

        /* First the (Newtonian) P2P */
        runner_dopair_grav_pp_full(e, ci_cache, cj_cache, gcount_i, gcount_j,
                                   gcount_padded_j, ci->gparts, cj->gparts);

        /* Then the M2P */
        runner_dopair_grav_pm(e, ci_cache, gcount_i, gcount_padded_i,
                              ci->gparts, CoM_j, multi_j, cj);
      }
      if (cj_active) {

        /* First the (Newtonian) P2P */
        runner_dopair_grav_pp_full(e, cj_cache, ci_cache, gcount_j, gcount_i,
                                   gcount_padded_i, cj->gparts, ci->gparts);

        /* Then the M2P */
        runner_dopair_grav_pm(e, cj_cache, gcount_j, gcount_padded_j,
                              cj->gparts, CoM_i, multi_i, ci);
      }
    }
  }

  /* Write back to the particles */
  if (ci_active) gravity_cache_write_back(ci_cache, ci->gparts, gcount_i);
  if (cj_active) gravity_cache_write_back(cj_cache, cj->gparts, gcount_j);

  TIMER_TOC(timer_dopair_grav_branch);
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

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we fit in cache */
  if (gcount > ci_cache->count)
    error("Not enough space in the cache! gcount=%d", gcount);
#endif

  /* Computed the padded counts */
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

  gravity_cache_populate_no_mpole(e->max_active_bin, ci_cache, gparts, gcount,
                                  gcount_padded, loc, c);

  /* Ok... Here we go ! */

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

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

      /* Interact! */
      float f_ij;
      runner_iact_grav_pp_full(r2, h2_i, h_inv_i, h_inv3_i, mass_j, &f_ij);

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

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we fit in cache */
  if (gcount > ci_cache->count)
    error("Not enough space in the caches! gcount=%d", gcount);
#endif

  /* Computed the padded counts */
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

  gravity_cache_populate_no_mpole(e->max_active_bin, ci_cache, gparts, gcount,
                                  gcount_padded, loc, c);

  /* Ok... Here we go ! */

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

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

      /* Interact! */
      float f_ij;
      runner_iact_grav_pp_truncated(r2, h2_i, h_inv_i, h_inv3_i, mass_j,
                                    rlr_inv, &f_ij);

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

  /* Check that we are not doing something stupid */
  if (c->split) error("Running P-P on a splitable cell");

  /* Do we need to start by drifting things ? */
  if (!cell_are_gpart_drifted(c, e)) error("Un-drifted gparts");

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {
    runner_doself_grav_pp_full(r, c);
  } else {

    /* Get the maximal distance between any two particles */
    const double max_r = 2. * c->multipole->r_max;

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
  const double theta_crit2 = props->theta_crit2;
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
  if (gravity_M2L_accept(multi_i->r_max, multi_j->r_max, theta_crit2, r2)) {

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
  const double theta_crit2 = props->theta_crit2;
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
    if (gravity_M2L_accept(multi_i->r_max, multi_j->r_max, theta_crit2, r2)) {

      /* Go for a (non-symmetric) M-M calculation */
      runner_dopair_grav_mm(r, ci, cj);

    } else {

      /* Let's check whether we need to still operate on this pair */

      /* Get the distance between the CoMs at the last rebuild*/
      double dx_rebuild = CoM_rebuild_i[0] - multi_j->CoM_rebuild[0];
      double dy_rebuild = CoM_rebuild_i[1] - multi_j->CoM_rebuild[1];
      double dz_rebuild = CoM_rebuild_i[2] - multi_j->CoM_rebuild[2];

      /* Apply BC */
      if (periodic) {
        dx_rebuild = nearest(dx_rebuild, dim[0]);
        dy_rebuild = nearest(dy_rebuild, dim[1]);
        dz_rebuild = nearest(dz_rebuild, dim[2]);
      }
      const double r2_rebuild = dx_rebuild * dx_rebuild +
                                dy_rebuild * dy_rebuild +
                                dz_rebuild * dz_rebuild;

      /* Is the criterion violated now but was OK at the last rebuild ? */
      if (gravity_M2L_accept(multi_i->r_max_rebuild, multi_j->r_max_rebuild,
                             theta_crit2, r2_rebuild)) {

        /* Alright, we have to take charge of that pair in a different way. */
        // MATTHIEU: We should actually open the tree-node here and recurse.
        runner_dopair_grav_mm(r, ci, cj);
      }
    }
  } /* Loop over top-level cells */

  if (timer) TIMER_TOC(timer_dograv_long_range);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
