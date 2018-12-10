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
#include "active.h"
#include "cell.h"
#include "gravity.h"
#include "gravity_cache.h"
#include "gravity_iact.h"
#include "inline.h"
#include "part.h"
#include "space_getsid.h"
#include "timers.h"

/**
 * @brief Recursively propagate the multipoles down the tree by applying the
 * L2L and L2P kernels.
 *
 * @param r The #runner.
 * @param c The #cell we are working on.
 * @param timer Are we timing this ?
 */
static INLINE void runner_do_grav_down(struct runner *r, struct cell *c,
                                       int timer) {

  /* Some constants */
  const struct engine *e = r->e;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.ti_old_multipole != e->ti_current)
    error("c->multipole not drifted.");
  if (c->grav.multipole->pot.ti_init != e->ti_current)
    error("c->field tensor not initialised");
#endif

  if (c->split) {

    /* Node case */

    /* Add the field-tensor to all the 8 progenitors */
    for (int k = 0; k < 8; ++k) {
      struct cell *cp = c->progeny[k];

      /* Do we have a progenitor with any active g-particles ? */
      if (cp != NULL && cell_is_active_gravity(cp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        if (cp->grav.ti_old_multipole != e->ti_current)
          error("cp->multipole not drifted.");
        if (cp->grav.multipole->pot.ti_init != e->ti_current)
          error("cp->field tensor not initialised");
#endif
        /* If the tensor received any contribution, push it down */
        if (c->grav.multipole->pot.interacted) {

          struct grav_tensor shifted_tensor;

          /* Shift the field tensor */
          gravity_L2L(&shifted_tensor, &c->grav.multipole->pot,
                      cp->grav.multipole->CoM, c->grav.multipole->CoM);

          /* Add it to this level's tensor */
          gravity_field_tensors_add(&cp->grav.multipole->pot, &shifted_tensor);
        }

        /* Recurse */
        runner_do_grav_down(r, cp, 0);
      }
    }

  } else {

    /* Leaf case */

    /* We can abort early if no interactions via multipole happened */
    if (!c->grav.multipole->pot.interacted) return;

    if (!cell_are_gpart_drifted(c, e)) error("Un-drifted gparts");

    /* Cell properties */
    struct gpart *gparts = c->grav.parts;
    const int gcount = c->grav.count;
    const struct grav_tensor *pot = &c->grav.multipole->pot;
    const double CoM[3] = {c->grav.multipole->CoM[0], c->grav.multipole->CoM[1],
                           c->grav.multipole->CoM[2]};

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
        if (c->grav.multipole->pot.ti_init != e->ti_current)
          error("c->field tensor not initialised");

        /* Check that we are not updated an inhibited particle */
        if (gpart_is_inhibited(gp, e)) error("Updating an inhibited particle!");

        /* Check that the particle was initialised */
        if (gp->initialised == 0)
          error("Adding forces to an un-initialised gpart.");
#endif
        /* Apply the kernel */
        gravity_L2P(pot, CoM, gp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_dograv_down);
}

/**
 * @brief Compute the non-truncated gravity interactions between all particles
 * of a cell and the particles of the other cell.
 *
 * The calculation is performed non-symmetrically using the pre-filled
 * #gravity_cache structures. The loop over the j cache should auto-vectorize.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param cj_cache #gravity_cache contaning the source particles.
 * @param gcount_i The number of particles in the cell i.
 * @param gcount_padded_j The number of particles in the cell j padded to the
 * vector length.
 * @param periodic Is the calculation using periodic BCs ?
 * @param dim The size of the simulation volume.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts_i The #gpart in cell i (for debugging checks only).
 * @param gparts_j The #gpart in cell j (for debugging checks only).
 * @param gcount_j The number of particles in the cell j (for debugging checks
 * only).
 */
static INLINE void runner_dopair_grav_pp_full(
    struct gravity_cache *restrict ci_cache,
    struct gravity_cache *restrict cj_cache, const int gcount_i,
    const int gcount_j, const int gcount_padded_j, const int periodic,
    const float dim[3], const struct engine *restrict e,
    struct gpart *restrict gparts_i, const struct gpart *restrict gparts_j) {

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

    /* Skip particle that can use the multipole */
    if (ci_cache->use_mpole[pid]) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (!gpart_is_active(&gparts_i[pid], e))
      error("Inactive particle went through the cache");
#endif

    const float x_i = ci_cache->x[pid];
    const float y_i = ci_cache->y[pid];
    const float z_i = ci_cache->z[pid];

    /* Some powers of the softening length */
    const float h_i = ci_cache->epsilon[pid];
    const float h2_i = h_i * h_i;
    const float h_inv_i = 1.f / h_i;
    const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

    /* Local accumulators for the acceleration and potential */
    float a_x = 0.f, a_y = 0.f, a_z = 0.f, pot = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(float, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->m, SWIFT_CACHE_ALIGNMENT);
    swift_assume_size(gcount_padded_j, VEC_SIZE);

    /* Loop over every particle in the other cell. */
    for (int pjd = 0; pjd < gcount_padded_j; pjd++) {

      /* Get info about j */
      const float x_j = cj_cache->x[pjd];
      const float y_j = cj_cache->y[pjd];
      const float z_j = cj_cache->z[pjd];
      const float mass_j = cj_cache->m[pjd];

      /* Compute the pairwise distance. */
      float dx = x_j - x_i;
      float dy = y_j - y_i;
      float dz = z_j - z_i;

      /* Correct for periodic BCs */
      if (periodic) {
        dx = nearestf(dx, dim[0]);
        dy = nearestf(dy, dim[1]);
        dz = nearestf(dz, dim[2]);
      }

      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f && h2_i == 0.)
        error("Interacting particles with 0 distance and 0 softening.");

      /* Check that particles have been drifted to the current time */
      if (gparts_i[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount_j && gparts_j[pjd].ti_drift != e->ti_current &&
          !gpart_is_inhibited(&gparts_j[pjd], e))
        error("gpj not drifted to current time");

      /* Check that we are not updated an inhibited particle */
      if (gpart_is_inhibited(&gparts_i[pid], e))
        error("Updating an inhibited particle!");

      /* Check that the particle we interact with was not inhibited */
      if (pjd < gcount_j && gpart_is_inhibited(&gparts_j[pjd], e) &&
          mass_j != 0.f)
        error("Inhibited particle used as gravity source.");

      /* Check that the particle was initialised */
      if (gparts_i[pid].initialised == 0)
        error("Adding forces to an un-initialised gpart.");
#endif

      /* Interact! */
      float f_ij, pot_ij;
      runner_iact_grav_pp_full(r2, h2_i, h_inv_i, h_inv3_i, mass_j, &f_ij,
                               &pot_ij);

      /* Store it back */
      a_x += f_ij * dx;
      a_y += f_ij * dy;
      a_z += f_ij * dz;
      pot += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount_j && !gpart_is_inhibited(&gparts_j[pjd], e))
        gparts_i[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] += a_x;
    ci_cache->a_y[pid] += a_y;
    ci_cache->a_z[pid] += a_z;
    ci_cache->pot[pid] += pot;
  }
}

/**
 * @brief Compute the truncated gravity interactions between all particles
 * of a cell and the particles of the other cell.
 *
 * The calculation is performed non-symmetrically using the pre-filled
 * #gravity_cache structures. The loop over the j cache should auto-vectorize.
 *
 * This function only makes sense in periodic BCs.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param cj_cache #gravity_cache contaning the source particles.
 * @param gcount_i The number of particles in the cell i.
 * @param gcount_padded_j The number of particles in the cell j padded to the
 * vector length.
 * @param dim The size of the simulation volume.
 * @param r_s_inv The inverse of the gravity-mesh smoothing-scale.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts_i The #gpart in cell i (for debugging checks only).
 * @param gparts_j The #gpart in cell j (for debugging checks only).
 * @param gcount_j The number of particles in the cell j (for debugging checks
 * only).
 */
static INLINE void runner_dopair_grav_pp_truncated(
    struct gravity_cache *restrict ci_cache,
    struct gravity_cache *restrict cj_cache, const int gcount_i,
    const int gcount_j, const int gcount_padded_j, const float dim[3],
    const float r_s_inv, const struct engine *restrict e,
    struct gpart *restrict gparts_i, const struct gpart *restrict gparts_j) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!e->s->periodic)
    error("Calling truncated PP function in non-periodic setup.");
#endif

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Skip inactive particles */
    if (!ci_cache->active[pid]) continue;

    /* Skip particle that can use the multipole */
    if (ci_cache->use_mpole[pid]) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (!gpart_is_active(&gparts_i[pid], e))
      error("Inactive particle went through the cache");
#endif

    const float x_i = ci_cache->x[pid];
    const float y_i = ci_cache->y[pid];
    const float z_i = ci_cache->z[pid];

    /* Some powers of the softening length */
    const float h_i = ci_cache->epsilon[pid];
    const float h2_i = h_i * h_i;
    const float h_inv_i = 1.f / h_i;
    const float h_inv3_i = h_inv_i * h_inv_i * h_inv_i;

    /* Local accumulators for the acceleration and potential */
    float a_x = 0.f, a_y = 0.f, a_z = 0.f, pot = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(float, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, cj_cache->m, SWIFT_CACHE_ALIGNMENT);
    swift_assume_size(gcount_padded_j, VEC_SIZE);

    /* Loop over every particle in the other cell. */
    for (int pjd = 0; pjd < gcount_padded_j; pjd++) {

      /* Get info about j */
      const float x_j = cj_cache->x[pjd];
      const float y_j = cj_cache->y[pjd];
      const float z_j = cj_cache->z[pjd];
      const float mass_j = cj_cache->m[pjd];

      /* Compute the pairwise distance. */
      float dx = x_j - x_i;
      float dy = y_j - y_i;
      float dz = z_j - z_i;

      /* Correct for periodic BCs */
      dx = nearestf(dx, dim[0]);
      dy = nearestf(dy, dim[1]);
      dz = nearestf(dz, dim[2]);

      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f && h2_i == 0.)
        error("Interacting particles with 0 distance and 0 softening.");

      /* Check that particles have been drifted to the current time */
      if (gparts_i[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount_j && gparts_j[pjd].ti_drift != e->ti_current &&
          !gpart_is_inhibited(&gparts_j[pjd], e))
        error("gpj not drifted to current time");

      /* Check that we are not updated an inhibited particle */
      if (gpart_is_inhibited(&gparts_i[pid], e))
        error("Updating an inhibited particle!");

      /* Check that the particle we interact with was not inhibited */
      if (pjd < gcount_j && gpart_is_inhibited(&gparts_j[pjd], e) &&
          mass_j != 0.f)
        error("Inhibited particle used as gravity source.");

      /* Check that the particle was initialised */
      if (gparts_i[pid].initialised == 0)
        error("Adding forces to an un-initialised gpart.");
#endif

      /* Interact! */
      float f_ij, pot_ij;
      runner_iact_grav_pp_truncated(r2, h2_i, h_inv_i, h_inv3_i, mass_j,
                                    r_s_inv, &f_ij, &pot_ij);

      /* Store it back */
      a_x += f_ij * dx;
      a_y += f_ij * dy;
      a_z += f_ij * dz;
      pot += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount_j && !gpart_is_inhibited(&gparts_j[pjd], e))
        gparts_i[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] += a_x;
    ci_cache->a_y[pid] += a_y;
    ci_cache->a_z[pid] += a_z;
    ci_cache->pot[pid] += pot;
  }
}

/**
 * @brief Compute the gravity interactions between all particles
 * of a cell and the multipole of the other cell.
 *
 * The calculation is performedusing the pre-filled
 * #gravity_cache structure. The loop over the i cache should auto-vectorize.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param gcount_padded_i The number of particles in the cell i padded to the
 * vector length.
 * @param CoM_j Position of the #multipole in #cell j.
 * @param multi_j The #multipole in #cell j.
 * @param periodic Is the calculation using periodic BCs ?
 * @param dim The size of the simulation volume.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts_i The #gpart in cell i (for debugging checks only).
 * @param gcount_i The number of particles in the cell i (for debugging checks
 * only).
 * @param cj The #cell j (for debugging checks only).
 */
static INLINE void runner_dopair_grav_pm_full(
    struct gravity_cache *ci_cache, const int gcount_padded_i,
    const float CoM_j[3], const struct multipole *restrict multi_j,
    const int periodic, const float dim[3], const struct engine *restrict e,
    struct gpart *restrict gparts_i, const int gcount_i,
    const struct cell *restrict cj) {

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, epsilon, ci_cache->epsilon,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_x, ci_cache->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, ci_cache->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, ci_cache->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pot, ci_cache->pot, SWIFT_CACHE_ALIGNMENT);
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

    /* Check that particles have been drifted to the current time */
    if (gparts_i[pid].ti_drift != e->ti_current)
      error("gpi not drifted to current time");

    /* Check that we are not updated an inhibited particle */
    if (gpart_is_inhibited(&gparts_i[pid], e))
      error("Updating an inhibited particle!");

    /* Check that the particle was initialised */
    if (gparts_i[pid].initialised == 0)
      error("Adding forces to an un-initialised gpart.");

    if (pid >= gcount_i) error("Adding forces to padded particle");
#endif

    const float x_i = x[pid];
    const float y_i = y[pid];
    const float z_i = z[pid];

    /* Some powers of the softening length */
    const float h_i = epsilon[pid];
    const float h_inv_i = 1.f / h_i;

    /* Distance to the Multipole */
    float dx = CoM_j[0] - x_i;
    float dy = CoM_j[1] - y_i;
    float dz = CoM_j[2] - z_i;

    /* Apply periodic BCs? */
    if (periodic) {
      dx = nearestf(dx, dim[0]);
      dy = nearestf(dy, dim[1]);
      dz = nearestf(dz, dim[2]);
    }

    const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
    const float r_max_j = cj->grav.multipole->r_max;
    const float r_max2 = r_max_j * r_max_j;
    const float theta_crit2 = e->gravity_properties->theta_crit2;

    /* Note: 1.1 to avoid FP rounding false-positives */
    if (!gravity_M2P_accept(r_max2, theta_crit2 * 1.1, r2))
      error(
          "use_mpole[i] set when M2P accept fails CoM=[%e %e %e] pos=[%e %e "
          "%e], rmax=%e",
          CoM_j[0], CoM_j[1], CoM_j[2], x_i, y_i, z_i, r_max_j);
#endif

    /* Interact! */
    float f_x, f_y, f_z, pot_ij;
    runner_iact_grav_pm_full(dx, dy, dz, r2, h_i, h_inv_i, multi_j, &f_x, &f_y,
                             &f_z, &pot_ij);

    /* Store it back */
    a_x[pid] += f_x;
    a_y[pid] += f_y;
    a_z[pid] += f_z;
    pot[pid] += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
    /* Update the interaction counter */
    if (pid < gcount_i)
      gparts_i[pid].num_interacted += cj->grav.multipole->m_pole.num_gpart;
#endif
  }
}

/**
 * @brief Compute the gravity interactions between all particles
 * of a cell and the multipole of the other cell.
 *
 * The calculation is performedusing the pre-filled
 * #gravity_cache structure. The loop over the i cache should auto-vectorize.
 *
 * This function only makes sense in periodic BCs.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param gcount_padded_i The number of particles in the cell i padded to the
 * vector length.
 * @param CoM_j Position of the #multipole in #cell j.
 * @param multi_j The #multipole in #cell j.
 * @param dim The size of the simulation volume.
 * @param r_s_inv The inverse of the gravity-mesh smoothing-scale.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts_i The #gpart in cell i (for debugging checks only).
 * @param gcount_i The number of particles in the cell i (for debugging checks
 * only).
 * @param cj The #cell j (for debugging checks only).
 */
static INLINE void runner_dopair_grav_pm_truncated(
    struct gravity_cache *ci_cache, const int gcount_padded_i,
    const float CoM_j[3], const struct multipole *restrict multi_j,
    const float dim[3], const float r_s_inv, const struct engine *restrict e,
    struct gpart *restrict gparts_i, const int gcount_i,
    const struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!e->s->periodic)
    error("Calling truncated PP function in non-periodic setup.");
#endif

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, epsilon, ci_cache->epsilon,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_x, ci_cache->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, ci_cache->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, ci_cache->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pot, ci_cache->pot, SWIFT_CACHE_ALIGNMENT);
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

    /* Check that particles have been drifted to the current time */
    if (gparts_i[pid].ti_drift != e->ti_current)
      error("gpi not drifted to current time");

    /* Check that we are not updated an inhibited particle */
    if (gpart_is_inhibited(&gparts_i[pid], e))
      error("Updating an inhibited particle!");

    /* Check that the particle was initialised */
    if (gparts_i[pid].initialised == 0)
      error("Adding forces to an un-initialised gpart.");

    if (pid >= gcount_i) error("Adding forces to padded particle");
#endif

    const float x_i = x[pid];
    const float y_i = y[pid];
    const float z_i = z[pid];

    /* Some powers of the softening length */
    const float h_i = epsilon[pid];
    const float h_inv_i = 1.f / h_i;

    /* Distance to the Multipole */
    float dx = CoM_j[0] - x_i;
    float dy = CoM_j[1] - y_i;
    float dz = CoM_j[2] - z_i;

    /* Apply periodic BCs */
    dx = nearestf(dx, dim[0]);
    dy = nearestf(dy, dim[1]);
    dz = nearestf(dz, dim[2]);

    const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
    const float r_max_j = cj->grav.multipole->r_max;
    const float r_max2 = r_max_j * r_max_j;
    const float theta_crit2 = e->gravity_properties->theta_crit2;

    /* 1.1 to avoid FP rounding false-positives */
    if (!gravity_M2P_accept(r_max2, theta_crit2 * 1.1, r2))
      error(
          "use_mpole[i] set when M2P accept fails CoM=[%e %e %e] pos=[%e %e "
          "%e], rmax=%e",
          CoM_j[0], CoM_j[1], CoM_j[2], x_i, y_i, z_i, r_max_j);
#endif

    /* Interact! */
    float f_x, f_y, f_z, pot_ij;
    runner_iact_grav_pm_truncated(dx, dy, dz, r2, h_i, h_inv_i, r_s_inv,
                                  multi_j, &f_x, &f_y, &f_z, &pot_ij);

    /* Store it back */
    a_x[pid] += f_x;
    a_y[pid] += f_y;
    a_z[pid] += f_z;
    pot[pid] += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
    /* Update the interaction counter */
    if (pid < gcount_i)
      gparts_i[pid].num_interacted += cj->grav.multipole->m_pole.num_gpart;
#endif
  }
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * This function switches between the full potential and the truncated one
 * depending on needs. It will also use the M2P (multipole) interaction
 * for the subset of particles in either cell for which the distance criterion
 * is valid.
 *
 * This function starts by constructing the require #gravity_cache for both
 * cells and then call the specialised functions doing the actual work on
 * the caches. It then write the data back to the particles.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param symmetric Are we updating both cells (1) or just ci (0) ?
 * @param allow_mpole Are we allowing the use of P2M interactions ?
 */
static INLINE void runner_dopair_grav_pp(struct runner *r, struct cell *ci,
                                         struct cell *cj, const int symmetric,
                                         const int allow_mpole) {

  /* Recover some useful constants */
  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float dim[3] = {(float)e->mesh->dim[0], (float)e->mesh->dim[1],
                        (float)e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;
  const double min_trunc = e->mesh->r_cut_min;

  TIMER_TIC;

  /* Record activity status */
  const int ci_active =
      cell_is_active_gravity(ci, e) && (ci->nodeID == e->nodeID);
  const int cj_active =
      cell_is_active_gravity(cj, e) && (cj->nodeID == e->nodeID);

  /* Anything to do here? */
  if (!ci_active && !cj_active) return;
  if (!ci_active && !symmetric) return;

  /* Check that we are not doing something stupid */
  if (ci->split || cj->split) error("Running P-P on splitable cells");

  /* Let's start by checking things are drifted */
  if (!cell_are_gpart_drifted(ci, e)) error("Un-drifted gparts");
  if (!cell_are_gpart_drifted(cj, e)) error("Un-drifted gparts");
  if (cj_active && ci->grav.ti_old_multipole != e->ti_current)
    error("Un-drifted multipole");
  if (ci_active && cj->grav.ti_old_multipole != e->ti_current)
    error("Un-drifted multipole");

  /* Caches to play with */
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;
  struct gravity_cache *const cj_cache = &r->cj_gravity_cache;

  /* Shift to apply to the particles in each cell */
  const double shift_i[3] = {0., 0., 0.};
  const double shift_j[3] = {0., 0., 0.};

  /* Recover the multipole info and shift the CoM locations */
  const float rmax_i = ci->grav.multipole->r_max;
  const float rmax_j = cj->grav.multipole->r_max;
  const float rmax2_i = rmax_i * rmax_i;
  const float rmax2_j = rmax_j * rmax_j;
  const struct multipole *multi_i = &ci->grav.multipole->m_pole;
  const struct multipole *multi_j = &cj->grav.multipole->m_pole;
  const float CoM_i[3] = {(float)(ci->grav.multipole->CoM[0] - shift_i[0]),
                          (float)(ci->grav.multipole->CoM[1] - shift_i[1]),
                          (float)(ci->grav.multipole->CoM[2] - shift_i[2])};
  const float CoM_j[3] = {(float)(cj->grav.multipole->CoM[0] - shift_j[0]),
                          (float)(cj->grav.multipole->CoM[1] - shift_j[1]),
                          (float)(cj->grav.multipole->CoM[2] - shift_j[2])};

  /* Start by constructing particle caches */

  /* Computed the padded counts */
  const int gcount_i = ci->grav.count;
  const int gcount_j = cj->grav.count;
  const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;
  const int gcount_padded_j = gcount_j - (gcount_j % VEC_SIZE) + VEC_SIZE;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we fit in cache */
  if (gcount_i > ci_cache->count || gcount_j > cj_cache->count)
    error("Not enough space in the caches! gcount_i=%d gcount_j=%d", gcount_i,
          gcount_j);
#endif

  /* Fill the caches */
  gravity_cache_populate(e->max_active_bin, allow_mpole, periodic, dim,
                         ci_cache, ci->grav.parts, gcount_i, gcount_padded_i,
                         shift_i, CoM_j, rmax2_j, ci, e->gravity_properties);
  gravity_cache_populate(e->max_active_bin, allow_mpole, periodic, dim,
                         cj_cache, cj->grav.parts, gcount_j, gcount_padded_j,
                         shift_j, CoM_i, rmax2_i, cj, e->gravity_properties);

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {

    /* Not periodic -> Can always use Newtonian potential */

    /* Let's updated the active cell(s) only */
    if (ci_active) {

      /* First the P2P */
      runner_dopair_grav_pp_full(ci_cache, cj_cache, gcount_i, gcount_j,
                                 gcount_padded_j, periodic, dim, e,
                                 ci->grav.parts, cj->grav.parts);

      /* Then the M2P */
      if (allow_mpole)
        runner_dopair_grav_pm_full(ci_cache, gcount_padded_i, CoM_j, multi_j,
                                   periodic, dim, e, ci->grav.parts, gcount_i,
                                   cj);
    }
    if (cj_active && symmetric) {

      /* First the P2P */
      runner_dopair_grav_pp_full(cj_cache, ci_cache, gcount_j, gcount_i,
                                 gcount_padded_i, periodic, dim, e,
                                 cj->grav.parts, ci->grav.parts);

      /* Then the M2P */
      if (allow_mpole)
        runner_dopair_grav_pm_full(cj_cache, gcount_padded_j, CoM_i, multi_i,
                                   periodic, dim, e, cj->grav.parts, gcount_j,
                                   ci);
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
        runner_dopair_grav_pp_truncated(ci_cache, cj_cache, gcount_i, gcount_j,
                                        gcount_padded_j, dim, r_s_inv, e,
                                        ci->grav.parts, cj->grav.parts);

        /* Then the M2P */
        if (allow_mpole)
          runner_dopair_grav_pm_truncated(ci_cache, gcount_padded_i, CoM_j,
                                          multi_j, dim, r_s_inv, e,
                                          ci->grav.parts, gcount_i, cj);
      }
      if (cj_active && symmetric) {

        /* First the (truncated) P2P */
        runner_dopair_grav_pp_truncated(cj_cache, ci_cache, gcount_j, gcount_i,
                                        gcount_padded_i, dim, r_s_inv, e,
                                        cj->grav.parts, ci->grav.parts);

        /* Then the M2P */
        if (allow_mpole)
          runner_dopair_grav_pm_truncated(cj_cache, gcount_padded_j, CoM_i,
                                          multi_i, dim, r_s_inv, e,
                                          cj->grav.parts, gcount_j, ci);
      }

    } else {

      /* Periodic but close-by cells can use the full Newtonian potential */

      /* Let's updated the active cell(s) only */
      if (ci_active) {

        /* First the (Newtonian) P2P */
        runner_dopair_grav_pp_full(ci_cache, cj_cache, gcount_i, gcount_j,
                                   gcount_padded_j, periodic, dim, e,
                                   ci->grav.parts, cj->grav.parts);

        /* Then the M2P */
        if (allow_mpole)
          runner_dopair_grav_pm_full(ci_cache, gcount_padded_i, CoM_j, multi_j,
                                     periodic, dim, e, ci->grav.parts, gcount_i,
                                     cj);
      }
      if (cj_active && symmetric) {

        /* First the (Newtonian) P2P */
        runner_dopair_grav_pp_full(cj_cache, ci_cache, gcount_j, gcount_i,
                                   gcount_padded_i, periodic, dim, e,
                                   cj->grav.parts, ci->grav.parts);

        /* Then the M2P */
        if (allow_mpole)
          runner_dopair_grav_pm_full(cj_cache, gcount_padded_j, CoM_i, multi_i,
                                     periodic, dim, e, cj->grav.parts, gcount_j,
                                     ci);
      }
    }
  }

  /* Write back to the particles */
  if (ci_active) gravity_cache_write_back(ci_cache, ci->grav.parts, gcount_i);
  if (cj_active && symmetric)
    gravity_cache_write_back(cj_cache, cj->grav.parts, gcount_j);

  TIMER_TOC(timer_dopair_grav_pp);
}

/**
 * @brief Compute the non-truncated gravity interactions between all particles
 * of a cell and the particles of the other cell.
 *
 * The calculation is performed non-symmetrically using the pre-filled
 * #gravity_cache structures. The loop over the j cache should auto-vectorize.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param gcount The number of particles in the cell.
 * @param gcount_padded The number of particles in the cell padded to the
 * vector length.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts The #gpart in the cell (for debugging checks only).
 */
static INLINE void runner_doself_grav_pp_full(
    struct gravity_cache *restrict ci_cache, const int gcount,
    const int gcount_padded, const struct engine *e, struct gpart *gparts) {

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
    float a_x = 0.f, a_y = 0.f, a_z = 0.f, pot = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(float, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
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
      /* Note: no need for periodic wrapping inside a cell */
      const float dx = x_j - x_i;
      const float dy = y_j - y_i;
      const float dz = z_j - z_i;
      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f && h2_i == 0.)
        error("Interacting particles with 0 distance and 0 softening.");

      /* Check that particles have been drifted to the current time */
      if (gparts[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount && gparts[pjd].ti_drift != e->ti_current &&
          !gpart_is_inhibited(&gparts[pjd], e))
        error("gpj not drifted to current time");

      /* Check that we are not updated an inhibited particle */
      if (gpart_is_inhibited(&gparts[pid], e))
        error("Updating an inhibited particle!");

      /* Check that the particle we interact with was not inhibited */
      if (pjd < gcount && gpart_is_inhibited(&gparts[pjd], e) && mass_j != 0.f)
        error("Inhibited particle used as gravity source.");

      /* Check that the particle was initialised */
      if (gparts[pid].initialised == 0)
        error("Adding forces to an un-initialised gpart.");
#endif

      /* Interact! */
      float f_ij, pot_ij;
      runner_iact_grav_pp_full(r2, h2_i, h_inv_i, h_inv3_i, mass_j, &f_ij,
                               &pot_ij);

      /* Store it back */
      a_x += f_ij * dx;
      a_y += f_ij * dy;
      a_z += f_ij * dz;
      pot += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount && !gpart_is_inhibited(&gparts[pjd], e))
        gparts[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] += a_x;
    ci_cache->a_y[pid] += a_y;
    ci_cache->a_z[pid] += a_z;
    ci_cache->pot[pid] += pot;
  }
}

/**
 * @brief Compute the truncated gravity interactions between all particles
 * of a cell and the particles of the other cell.
 *
 * The calculation is performed non-symmetrically using the pre-filled
 * #gravity_cache structures. The loop over the j cache should auto-vectorize.
 *
 * This function only makes sense in periodic BCs.
 *
 * @param ci_cache #gravity_cache contaning the particles to be updated.
 * @param gcount The number of particles in the cell.
 * @param gcount_padded The number of particles in the cell padded to the
 * vector length.
 * @param r_s_inv The inverse of the gravity-mesh smoothing-scale.
 *
 * @param e The #engine (for debugging checks only).
 * @param gparts The #gpart in the cell (for debugging checks only).
 */
static INLINE void runner_doself_grav_pp_truncated(
    struct gravity_cache *restrict ci_cache, const int gcount,
    const int gcount_padded, const float r_s_inv, const struct engine *e,
    struct gpart *gparts) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!e->s->periodic)
    error("Calling truncated PP function in non-periodic setup.");
#endif

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

    /* Local accumulators for the acceleration and potential */
    float a_x = 0.f, a_y = 0.f, a_z = 0.f, pot = 0.f;

    /* Make the compiler understand we are in happy vectorization land */
    swift_align_information(float, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_align_information(float, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
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
      /* Note: no need for periodic wrapping inside a cell */
      const float dx = x_j - x_i;
      const float dy = y_j - y_i;
      const float dz = z_j - z_i;

      const float r2 = dx * dx + dy * dy + dz * dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (r2 == 0.f && h2_i == 0.)
        error("Interacting particles with 0 distance and 0 softening.");

      /* Check that particles have been drifted to the current time */
      if (gparts[pid].ti_drift != e->ti_current)
        error("gpi not drifted to current time");
      if (pjd < gcount && gparts[pjd].ti_drift != e->ti_current &&
          !gpart_is_inhibited(&gparts[pjd], e))
        error("gpj not drifted to current time");

      /* Check that we are not updated an inhibited particle */
      if (gpart_is_inhibited(&gparts[pid], e))
        error("Updating an inhibited particle!");

      /* Check that the particle we interact with was not inhibited */
      if (pjd < gcount && gpart_is_inhibited(&gparts[pjd], e) && mass_j != 0.f)
        error("Inhibited particle used as gravity source.");

      /* Check that the particle was initialised */
      if (gparts[pid].initialised == 0)
        error("Adding forces to an un-initialised gpart.");
#endif

      /* Interact! */
      float f_ij, pot_ij;
      runner_iact_grav_pp_truncated(r2, h2_i, h_inv_i, h_inv3_i, mass_j,
                                    r_s_inv, &f_ij, &pot_ij);

      /* Store it back */
      a_x += f_ij * dx;
      a_y += f_ij * dy;
      a_z += f_ij * dz;
      pot += pot_ij;

#ifdef SWIFT_DEBUG_CHECKS
      /* Update the interaction counter if it's not a padded gpart */
      if (pjd < gcount && !gpart_is_inhibited(&gparts[pjd], e))
        gparts[pid].num_interacted++;
#endif
    }

    /* Store everything back in cache */
    ci_cache->a_x[pid] += a_x;
    ci_cache->a_y[pid] += a_y;
    ci_cache->a_z[pid] += a_z;
    ci_cache->pot[pid] += pot;
  }
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * other ones.
 *
 * This function switches between the full potential and the truncated one
 * depending on needs.
 *
 * This function starts by constructing the require #gravity_cache for the
 * cell and then call the specialised functions doing the actual work on
 * the cache. It then write the data back to the particles.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
static INLINE void runner_doself_grav_pp(struct runner *r, struct cell *c) {

  /* Recover some useful constants */
  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float r_s_inv = e->mesh->r_s_inv;
  const double min_trunc = e->mesh->r_cut_min;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.count == 0) error("Doing self gravity on an empty cell !");
#endif

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Check that we are not doing something stupid */
  if (c->split) error("Running P-P on a splitable cell");

  /* Do we need to start by drifting things ? */
  if (!cell_are_gpart_drifted(c, e)) error("Un-drifted gparts");

  /* Start by constructing a cache for the particles */
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;

  /* Shift to apply to the particles in the cell */
  const double loc[3] = {c->loc[0] + 0.5 * c->width[0],
                         c->loc[1] + 0.5 * c->width[1],
                         c->loc[2] + 0.5 * c->width[2]};

  /* Computed the padded counts */
  const int gcount = c->grav.count;
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we fit in cache */
  if (gcount > ci_cache->count)
    error("Not enough space in the cache! gcount=%d", gcount);
#endif

  /* Fill the cache */
  gravity_cache_populate_no_mpole(e->max_active_bin, ci_cache, c->grav.parts,
                                  gcount, gcount_padded, loc, c,
                                  e->gravity_properties);

  /* Can we use the Newtonian version or do we need the truncated one ? */
  if (!periodic) {

    /* Not periodic -> Can always use Newtonian potential */
    runner_doself_grav_pp_full(ci_cache, gcount, gcount_padded, e,
                               c->grav.parts);

  } else {

    /* Get the maximal distance between any two particles */
    const double max_r = 2. * c->grav.multipole->r_max;

    /* Do we need to use the truncated interactions ? */
    if (max_r > min_trunc) {

      /* Periodic but far-away cells must use the truncated potential */
      runner_doself_grav_pp_truncated(ci_cache, gcount, gcount_padded, r_s_inv,
                                      e, c->grav.parts);

    } else {

      /* Periodic but close-by cells can use the full Newtonian potential */
      runner_doself_grav_pp_full(ci_cache, gcount, gcount_padded, e,
                                 c->grav.parts);
    }
  }

  /* Write back to the particles */
  gravity_cache_write_back(ci_cache, c->grav.parts, gcount);

  TIMER_TOC(timer_doself_grav_pp);
}

/**
 * @brief Computes the interaction of the field tensor and multipole
 * of two cells symmetrically.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
static INLINE void runner_dopair_grav_mm_symmetric(struct runner *r,
                                                   struct cell *restrict ci,
                                                   struct cell *restrict cj) {

  /* Some constants */
  const struct engine *e = r->e;
  const struct gravity_props *props = e->gravity_properties;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;

  TIMER_TIC;

  /* Anything to do here? */
  if ((!cell_is_active_gravity_mm(ci, e) || ci->nodeID != engine_rank) ||
      (!cell_is_active_gravity_mm(cj, e) || cj->nodeID != engine_rank))
    error("Invalid state in symmetric M-M calculation!");

  /* Short-cut to the multipole */
  const struct multipole *multi_i = &ci->grav.multipole->m_pole;
  const struct multipole *multi_j = &cj->grav.multipole->m_pole;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_i->num_gpart == 0)
    error("Multipole i does not seem to have been set.");

  if (multi_j->num_gpart == 0)
    error("Multipole j does not seem to have been set.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("ci->grav tensor not initialised.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("cj->grav tensor not initialised.");

  if (ci->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole ci->grav.ti_old_multipole=%lld ci->nodeID=%d "
        "cj->nodeID=%d e->ti_current=%lld",
        ci->grav.ti_old_multipole, ci->nodeID, cj->nodeID, e->ti_current);

  if (cj->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole cj->grav.ti_old_multipole=%lld cj->nodeID=%d "
        "ci->nodeID=%d e->ti_current=%lld",
        cj->grav.ti_old_multipole, cj->nodeID, ci->nodeID, e->ti_current);
#endif

  /* Let's interact at this level */
  gravity_M2L_symmetric(&ci->grav.multipole->pot, &cj->grav.multipole->pot,
                        multi_i, multi_j, ci->grav.multipole->CoM,
                        cj->grav.multipole->CoM, props, periodic, dim, r_s_inv);

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Computes the interaction of the field tensor in a cell with the
 * multipole of another cell.
 *
 * @param r The #runner.
 * @param ci The #cell with field tensor to interact.
 * @param cj The #cell with the multipole.
 */
static INLINE void runner_dopair_grav_mm_nonsym(
    struct runner *r, struct cell *restrict ci,
    const struct cell *restrict cj) {

  /* Some constants */
  const struct engine *e = r->e;
  const struct gravity_props *props = e->gravity_properties;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity_mm(ci, e) || ci->nodeID != engine_rank) return;

  /* Short-cut to the multipole */
  const struct multipole *multi_j = &cj->grav.multipole->m_pole;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_j->num_gpart == 0)
    error("Multipole does not seem to have been set.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("ci->grav tensor not initialised.");

  if (cj->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole cj->grav.ti_old_multipole=%lld cj->nodeID=%d "
        "ci->nodeID=%d e->ti_current=%lld",
        cj->grav.ti_old_multipole, cj->nodeID, ci->nodeID, e->ti_current);
#endif

  /* Let's interact at this level */
  gravity_M2L_nonsym(&ci->grav.multipole->pot, multi_j, ci->grav.multipole->CoM,
                     cj->grav.multipole->CoM, props, periodic, dim, r_s_inv);

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Call the M-M calculation on two cells if active.
 *
 * @param r The #runner object.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
static INLINE void runner_dopair_grav_mm(struct runner *r,
                                         struct cell *restrict ci,
                                         struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* What do we need to do? */
  const int do_i =
      cell_is_active_gravity_mm(ci, e) && (ci->nodeID == e->nodeID);
  const int do_j =
      cell_is_active_gravity_mm(cj, e) && (cj->nodeID == e->nodeID);

  /* Do we need drifting first? */
  if (ci->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(ci, e);
  if (cj->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(cj, e);

  /* Interact! */
  if (do_i && do_j)
    runner_dopair_grav_mm_symmetric(r, ci, cj);
  else if (do_i)
    runner_dopair_grav_mm_nonsym(r, ci, cj);
  else if (do_j)
    runner_dopair_grav_mm_nonsym(r, cj, ci);
}

/**
 * @brief Computes all the M-M interactions between all the well-separated (at
 * rebuild) pairs of progenies of the two cells.
 *
 * @param r The #runner thread.
 * @param flags The task flag containing the list of well-separated pairs as a
 * bit-field.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
static INLINE void runner_dopair_grav_mm_progenies(struct runner *r,
                                                   const long long flags,
                                                   struct cell *restrict ci,
                                                   struct cell *restrict cj) {

  /* Loop over all pairs of progenies */
  for (int i = 0; i < 8; i++) {
    if (ci->progeny[i] != NULL) {
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] != NULL) {

          struct cell *cpi = ci->progeny[i];
          struct cell *cpj = cj->progeny[j];

          const int flag = i * 8 + j;

          /* Did we agree to use an M-M interaction here at the last rebuild? */
          if (flags & (1ULL << flag)) runner_dopair_grav_mm(r, cpi, cpj);
        }
      }
    }
  }
}

static INLINE void runner_dopair_recursive_grav_pm(struct runner *r,
                                                   struct cell *ci,
                                                   const struct cell *cj) {
  /* Some constants */
  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float dim[3] = {(float)e->mesh->dim[0], (float)e->mesh->dim[1],
                        (float)e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;

  /* Anything to do here? */
  if (!(cell_is_active_gravity(ci, e) && ci->nodeID == e->nodeID)) return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Early abort? */
  if (ci->grav.count == 0 || cj->grav.count == 0)
    error("Doing pair gravity on an empty cell !");

  /* Sanity check */
  if (ci == cj) error("Pair interaction between a cell and itself.");

  if (cj->grav.ti_old_multipole != e->ti_current)
    error("cj->grav.multipole not drifted.");
#endif

  /* Can we recurse further? */
  if (ci->split) {

    /* Loop over ci's children */
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        runner_dopair_recursive_grav_pm(r, ci->progeny[k], cj);
    }

    /* Ok, let's do the interaction here */
  } else {

    /* Start by constructing particle caches */

    /* Cache to play with */
    struct gravity_cache *const ci_cache = &r->ci_gravity_cache;

    /* Computed the padded counts */
    const int gcount_i = ci->grav.count;
    const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we fit in cache */
    if (gcount_i > ci_cache->count)
      error("Not enough space in the cache! gcount_i=%d", gcount_i);
#endif

    /* Recover the multipole info and the CoM locations */
    const struct multipole *multi_j = &cj->grav.multipole->m_pole;
    const float r_max = cj->grav.multipole->r_max;
    const float CoM_j[3] = {(float)(cj->grav.multipole->CoM[0]),
                            (float)(cj->grav.multipole->CoM[1]),
                            (float)(cj->grav.multipole->CoM[2])};

    /* Fill the cache */
    gravity_cache_populate_all_mpole(
        e->max_active_bin, periodic, dim, ci_cache, ci->grav.parts, gcount_i,
        gcount_padded_i, ci, CoM_j, r_max * r_max, e->gravity_properties);

    /* Can we use the Newtonian version or do we need the truncated one ? */
    if (!periodic) {

      runner_dopair_grav_pm_full(ci_cache, gcount_padded_i, CoM_j, multi_j,
                                 periodic, dim, e, ci->grav.parts, gcount_i,
                                 cj);

    } else {

      runner_dopair_grav_pm_truncated(ci_cache, gcount_padded_i, CoM_j, multi_j,
                                      dim, r_s_inv, e, ci->grav.parts, gcount_i,
                                      cj);
    }

    /* Write back to the particles */
    gravity_cache_write_back(ci_cache, ci->grav.parts, gcount_i);
  }
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * This function will try to recurse as far down the tree as possible and only
 * default to direct summation if there is no better option.
 *
 * If using periodic BCs, we will abort the recursion if th distance between the
 * cells is larger than the set threshold.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param gettimer Are we timing this ?
 */
static INLINE void runner_dopair_recursive_grav(struct runner *r,
                                                struct cell *ci,
                                                struct cell *cj, int gettimer) {

  /* Some constants */
  const struct engine *e = r->e;
  const int nodeID = e->nodeID;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const double theta_crit2 = e->gravity_properties->theta_crit2;
  const double max_distance = e->mesh->r_cut_max;

  /* Anything to do here? */
  if (!((cell_is_active_gravity(ci, e) && ci->nodeID == nodeID) ||
        (cell_is_active_gravity(cj, e) && cj->nodeID == nodeID)))
    return;

#ifdef SWIFT_DEBUG_CHECKS

  const int gcount_i = ci->grav.count;
  const int gcount_j = cj->grav.count;

  /* Early abort? */
  if (gcount_i == 0 || gcount_j == 0)
    error("Doing pair gravity on an empty cell !");

  /* Sanity check */
  if (ci == cj) error("Pair interaction between a cell and itself.");

  if (cell_is_active_gravity(ci, e) &&
      ci->grav.ti_old_multipole != e->ti_current)
    error("ci->grav.multipole not drifted.");
  if (cell_is_active_gravity(cj, e) &&
      cj->grav.ti_old_multipole != e->ti_current)
    error("cj->grav.multipole not drifted.");
#endif

  TIMER_TIC;

  /* Recover the multipole information */
  struct gravity_tensors *const multi_i = ci->grav.multipole;
  struct gravity_tensors *const multi_j = cj->grav.multipole;

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

  /* Minimal distance between any 2 particles in the two cells */
  const double r_lr_check = sqrt(r2) - (multi_i->r_max + multi_j->r_max);

  /* Are we beyond the distance where the truncated forces are 0? */
  if (periodic && r_lr_check > max_distance) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Need to account for the interactions we missed */
    if (cell_is_active_gravity(ci, e))
      multi_i->pot.num_interacted += multi_j->m_pole.num_gpart;
    if (cell_is_active_gravity(cj, e))
      multi_j->pot.num_interacted += multi_i->m_pole.num_gpart;
#endif
    return;
  }

  /* OK, we actually need to compute this pair. Let's find the cheapest
   * option... */

  /* Can we use M-M interactions ? */
  if (gravity_M2L_accept(multi_i->r_max, multi_j->r_max, theta_crit2, r2)) {

    /* Go M-M */
    runner_dopair_grav_mm(r, ci, cj);

  } else if (!ci->split && !cj->split) {

    /* We have two leaves. Go P-P. */
    runner_dopair_grav_pp(r, ci, cj, /*symmetric*/ 1, /*allow_mpoles*/ 1);

  } else {

    /* Alright, we'll have to split and recurse. */
    /* We know at least one of ci and cj is splittable */

    const double ri_max = multi_i->r_max;
    const double rj_max = multi_j->r_max;

    /* Split the larger of the two cells and start over again */
    if (ri_max > rj_max) {

      /* Can we actually split that interaction ? */
      if (ci->split) {

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL)
            runner_dopair_recursive_grav(r, ci->progeny[k], cj, 0);
        }

      } else {
        /* cj is split */

        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL)
            runner_dopair_recursive_grav(r, ci, cj->progeny[k], 0);
        }
      }
    } else {

      /* Can we actually split that interaction ? */
      if (cj->split) {

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL)
            runner_dopair_recursive_grav(r, ci, cj->progeny[k], 0);
        }

      } else {
        /* ci is split */

        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL)
            runner_dopair_recursive_grav(r, ci->progeny[k], cj, 0);
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
}

/**
 * @brief Computes the interaction of all the particles in a cell.
 *
 * This function will try to recurse as far down the tree as possible and only
 * default to direct summation if there is no better option.
 *
 * @param r The #runner.
 * @param c The first #cell.
 * @param gettimer Are we timing this ?
 */
static INLINE void runner_doself_recursive_grav(struct runner *r,
                                                struct cell *c, int gettimer) {

  /* Some constants */
  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  /* Early abort? */
  if (c->grav.count == 0) error("Doing self gravity on an empty cell !");
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* If the cell is split, interact each progeny with itself, and with
     each of its siblings. */
  if (c->split) {

    for (int j = 0; j < 8; j++) {
      if (c->progeny[j] != NULL) {

        runner_doself_recursive_grav(r, c->progeny[j], 0);

        for (int k = j + 1; k < 8; k++) {
          if (c->progeny[k] != NULL) {

            runner_dopair_recursive_grav(r, c->progeny[j], c->progeny[k], 0);
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
static INLINE void runner_do_grav_long_range(struct runner *r, struct cell *ci,
                                             int timer) {

  /* Some constants */
  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const double theta_crit2 = e->gravity_properties->theta_crit2;
  const double max_distance2 = e->mesh->r_cut_max * e->mesh->r_cut_max;

  TIMER_TIC;

  /* Recover the list of top-level cells */
  struct cell *cells = e->s->cells_top;
  int *cells_with_particles = e->s->cells_with_particles_top;
  const int nr_cells_with_particles = e->s->nr_cells_with_particles;

  /* Anything to do here? */
  if (!cell_is_active_gravity(ci, e)) return;

  if (ci->nodeID != engine_rank)
    error("Non-local cell in long-range gravity task!");

  /* Check multipole has been drifted */
  if (ci->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(ci, e);

  /* Get this cell's multipole information */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Find this cell's top-level (great-)parent */
  struct cell *top = ci;
  while (top->parent != NULL) top = top->parent;

  /* Recover the top-level multipole (for distance checks) */
  struct gravity_tensors *const multi_top = top->grav.multipole;
  const double CoM_rebuild_top[3] = {multi_top->CoM_rebuild[0],
                                     multi_top->CoM_rebuild[1],
                                     multi_top->CoM_rebuild[2]};

  /* Loop over all the top-level cells and go for a M-M interaction if
   * well-separated */
  for (int n = 0; n < nr_cells_with_particles; ++n) {

    /* Handle on the top-level cell and it's gravity business*/
    const struct cell *cj = &cells[cells_with_particles[n]];
    const struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Can we escape early in the periodic BC case? */
    if (periodic) {

      /* Minimal distance between any pair of particles */
      const double min_radius2 =
          cell_min_dist2_same_size(top, cj, periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0 ?*/
      if (min_radius2 > max_distance2) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Need to account for the interactions we missed */
        multi_i->pot.num_interacted += multi_j->m_pole.num_gpart;
#endif

        /* Record that this multipole received a contribution */
        multi_i->pot.interacted = 1;

        /* We are done here. */
        continue;
      }
    }

    /* Get the distance between the CoMs at the last rebuild*/
    double dx_r = CoM_rebuild_top[0] - multi_j->CoM_rebuild[0];
    double dy_r = CoM_rebuild_top[1] - multi_j->CoM_rebuild[1];
    double dz_r = CoM_rebuild_top[2] - multi_j->CoM_rebuild[2];

    /* Apply BC */
    if (periodic) {
      dx_r = nearest(dx_r, dim[0]);
      dy_r = nearest(dy_r, dim[1]);
      dz_r = nearest(dz_r, dim[2]);
    }
    const double r2_rebuild = dx_r * dx_r + dy_r * dy_r + dz_r * dz_r;

    /* Are we in charge of this cell pair? */
    if (gravity_M2L_accept(multi_top->r_max_rebuild, multi_j->r_max_rebuild,
                           theta_crit2, r2_rebuild)) {

      /* Call the PM interaction fucntion on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  }   /* Loop over top-level cells */

  if (timer) TIMER_TOC(timer_dograv_long_range);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
