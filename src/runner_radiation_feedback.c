/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "fof.h"
#include "gravity.h"
#include "hydro.h"
#include "runner_doiact_sinks.h"
#include "space.h"
#include "star_formation.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"

/**
 * @brief Top-level function for HII ionization feedback.
 *
 * This function recurses down the cell hierarchy until it reaches a level
 * where the cell size is comparable to the maximum HII radius of the stars
 * within it. At that point, it calls the 'branch' function to handle the
 * actual neighbor gathering.
 *
 * @param r The #runner thread.
 * @param c The #cell to process.
 * @param timer If true, records the timing of this operation.
 */
void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer) {
#ifdef IONIZATION_FEEDBACK_LOOP
  struct engine *e = r->e;

  /* TODO: Or h_max_old? */
  /* TODO: Which cell h_hii_max shall I look at? The current code changes the
     value as we do deeper in the cell hierarchy. Ideally, we want to use the
     smallest hmax for a cell, process the stars and parts at this level for
     this cell, and go higher for cells/stars with higher h_hii */
  const float r_hii_max = c->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;
  const int can_recurse = c->split && (interaction_limit < 0.5f * c->dmin);

  if (c->stars.count == 0 || c->hydro.count == 0 || r_hii_max <= 0.f ||
      !cell_is_active_stars(c, e))
    return;

  /* message("[%lld], hydro super = %lld", c->cellID, c->hydro.super->cellID);
   */

  TIMER_TIC;

  /*
   * Decision: Should we process this cell or go deeper?
   * We stop recursing if:
   * 1. The cell is a leaf (!c->split).
   * 2. The cell is small enough that we can't efficiently split the
   *    interaction radius further (similar to Mosaics logic).
   * 3. The R_HII_max covers a significant fraction of the cell.
   */
  /* Is the cell split and not smaller than the smoothing length? */
  if (can_recurse) {
    /* Keep recursing deeper into the super-cell hierarchy. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        runner_do_stars_hii_ionization_feedback(r, c->progeny[k], 0);
    }
  } else {
    /* We have reached the 'Working Level' */
    runner_do_stars_hii_ionization_feedback_branch(r, c);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Branching point to coordinate HII ionization for stars in a cell.
 *
 * This function manages the neighbor gathering for all active stars in the
 * cell. It allocates a temporary buffer, gathers neighbors from the cell itself
 * (self) and its neighbors (pair) by climbing the cell hierarchy, and
 * prepares the gathered particles for sorting and ionization.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the stars to process.
 */
void runner_do_stars_hii_ionization_feedback_branch(struct runner *r,
                                                    struct cell *c) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  /* const struct feedback_props *fb_props = e->feedback_props; */
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;

  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  /* OR h_max_old? */
  const float r_hii_max = c->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;

  if (c->stars.count == 0 || c->hydro.count == 0 || r_hii_max <= 0.f ||
      !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > c->dmin)
    error("Cell smaller than HII interaction length");
#endif

  /* message("[c = %lld, super = %lld] interaction_limit = %e, cell_dmin = %e",
   */
  /*         c->cellID, c->hydro.super->cellID, interaction_limit, c->dmin); */

  /* TODO: Handle the case where this value is too small. */
  /* TODO: Add multiple tries if the star has not exhausted its photons */
  const int max_ngbs = 1024;
  struct hii_neighbor *ngb_buffer =
      malloc(max_ngbs * sizeof(struct hii_neighbor));

  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[sid];

    /* Is this part within the timestep? */
    if (!spart_is_active(si, e) && !feedback_is_active(si, e) &&
        feedback_is_HII_ionization_active(si, e))
      return;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that particles have been drifted to the current time */
    if (si->ti_drift != e->ti_current)
      error("Particle si (%lld) not drifted to current time c = %lld, "
              "c->super = %lld", si->id, c->cellID, c->super->cellID);

#endif

    int count_found = 0;

    /***************************************************/
    /* First loop over particles in the current cell */
    runner_do_stars_hii_ionization_feedback_self(r, c, si, ngb_buffer, max_ngbs,
                                                 &count_found);

    /***************************************************/
    /* Now loop over particles in the neighboring cells */

    /* Climb up the cell hierarchy. */
    for (struct cell *finger = c; finger != NULL; finger = finger->parent) {

      /* These are defined at the super level... When we reach the progeny,
         these will be NULL... So we need to grab the finger first */
      for (struct link *l = finger->stars.density; l != NULL; l = l->next) {
        /* We have already handled the self case */
        if (l->t->type == task_type_self) continue;

        struct cell *c_in = l->t->cj;

        /* Now, find the correct level... */
        const int can_recurse_j =
            c_in->split && (interaction_limit < 0.5f * c_in->dmin);
        if (can_recurse_j) {
          /* Keep recursing deeper into the super-cell hierarchy. */
          for (int k = 0; k < 8; k++)
            if (c_in->progeny[k] != NULL)
              runner_do_stars_hii_ionization_feedback_pair(
                  r, c, c_in->progeny[k], si, ngb_buffer, max_ngbs,
                  &count_found);
        } else {
          /* We have reached the 'Working Level' */
          runner_do_stars_hii_ionization_feedback_pair(
              r, c, c_in, si, ngb_buffer, max_ngbs, &count_found);
        } /* Recurse */
      } /* Neighbour search */
    } /* Climb up in the cell hierarchy */

    /***************************************************/
    /* It's time to sort the gas particles */
    if (count_found > 0) {
      runner_sort_hii_neighbors(ngb_buffer, count_found);

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify that the neighbors are properly sorted by distance */
      for (int k = 0; k < count_found - 1; k++) {
        if (ngb_buffer[k].r2 > ngb_buffer[k + 1].r2) {
          error(
              "HII neighbor buffer not properly sorted! "
              "Index %d (r2=%e) is larger than index %d (r2=%e).",
              k, ngb_buffer[k].r2, k + 1, ngb_buffer[k + 1].r2);
        }
      }
#endif

      /* Now let's ionize the gas particles! */
      for (int k = 0; k < count_found; k++) {
        struct part *pj = ngb_buffer[k].p;
        struct xpart *xpj = ngb_buffer[k].xp;

        /* Do the ionization */
	/* Photoionization - HII region:
	   From step 1 we know N_dot_ion = L_ion / (h nu_ion) . We will use a
	   simple stromgren sphere approximation. For each particle:
	   a) Test if the particle is already ionized : T > 10^4 or particle was
	   flagged to be in an ionized region.
	   b) If it is not ionized, compute the ioninzing rate needed to fully
	   ionize:
	   \Delta N_dot_j = N(H)_j beta n_e_j
	   N(H)_j = X_H m_part / (mu m_proton) (the number of H atoms)
	   with beta = 3e-13 cm^3 / s is the recombination coefficient, n_e_j the
	   electron number density assuming full ionization, X_H is the hydrogen
	   mass fraction and mu the molecular weight.
	   c) If \Delta N_dot_j <= N_dot_ion:
	   tag the particle as being in a HII region
	   consume the photons: N_ion -= Delta N_dot_j

	   else:
	   determine randomly if the particle is ionized by computing the
	   proba p = N_dot_io / \Delta N_dot_j
	   If rand_number <= proba :
	   tag the particle as being in a HII region
	   consume the photons: N_ion -= Delta N_dot_j

	   d) For the particles tagged as ionized:
	   set the temperature (internal energy) to the
	   min(current temperature + heat added from the energy of the ionisation,
	   equilibrium HII region tem from collisional cooling)
	   Set the incident rad and FUV flux to the stromgren value --> to
	   compute the inonizing

	   Concretely,
	   u_new = min(u + delta U, U_collisional),
	   delta U = N_H * E_ion / m_gas,
	   E_ion = 13.6 eV = 2.18e-11 erg
	   Gamma = \Delta N_dot_j / N_H.

	*/
        /* 1. Tag gas as ionized (atomics)
           2. Flag to be synchronized (atomics)
           3. Consume photons from the star
           4. Compute M_HII and r_HII for the star
           5. Update the stars r_HII.
           6. Update the cell max r_HII after finising processing all star
           particles
        */
        /* feedback_do_HII_ionization(p, xp); */
        const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(
            phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

        if (Delta_dot_N_ion <= radiation_get_star_ionization_rate(si) &&
            atomic_cas(&xpj->feedback_data.radiation.is_ionized, 0, 1) == 0) {

          /* Flag the particle to be synchronized on the timeline. */
          atomic_or(&pj->limiter_data.to_be_synchronized, 1);

          /* Consume photons from the star */
          radiation_consume_ionizing_photons(si, Delta_dot_N_ion);

          /* Update the star's HII mass */
          /* si->feedback_data.mass_HII += hydro_get_mass(pj); */

          /* Update the star's r_HII to the distance of this particle */
          si->h_hii = sqrtf(ngb_buffer[k].r2) * kernel_gamma_inv;
        } else {
          /* Particle was already ionized by another star in this sub-cycle.
           * We skip it and move to the next nearest neighbor. */

	  /* TODO: Handle that */
	  /* /\* If we cannot fully ionize, compute a probability to determine if */
	  /*    we */
	  /*    fully ionize pj or not and draw the random number.  *\/ */
	  /* const float proba = dot_N_ion / Delta_dot_N_ion; */
	  /* const float random_number = */
          /*   random_unit_interval(sp->id, ti_begin, */
	  /* 			 random_number_HII_regions); */

	  /* /\* If we are lucky, do the ionization *\/ */
	  /* if (random_number <= proba) { */
	  /*   /\* Update the Stromgren sphere radius *\/ */
	  /*   sp->feedback_data.radiation.R_stromgren = stromgren[i].distance; */
	  /* } */

	  /* /\* Consume the photons in all cases *\/ */
	  /* radiation_consume_ionizing_photons(sp, Delta_dot_N_ion); */
          continue;
        }
      }
    } /* Loop over the sorted particles */

  } /* Loop over sparts */
  free(ngb_buffer);
}

/**
 * @brief Gather gas particles within the HII radius from the star's own cell.
 *
 * Performs a brute-force search over all hydro particles in cell @p c that
 * fall within the HII smoothing length of the star @p si. Candidates are
 * added to the @p buffer if they are not already inhibited or ionized.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_do_stars_hii_ionization_feedback_self(
    struct runner *r, struct cell *c, struct spart *si,
    struct hii_neighbor *buffer, int max_size, int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  const float hi = si->h_hii;
  const float hig2 = hi * hi * kernel_gamma2;
  const float six[3] = {(float)(si->x[0] - c->loc[0]),
                        (float)(si->x[1] - c->loc[1]),
                        (float)(si->x[2] - c->loc[2])};

  /* Loop over the parts in c. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = &xparts[pjd];

    /* const float hj = pj->h; */

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error(
              "Particle pj (%lld) not drifted to current time. c = %lld, "
              "c->super = %lld",
              pj->id, c->cellID, c->super->cellID);
#endif

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                          (float)(pj->x[1] - c->loc[1]),
                          (float)(pj->x[2] - c->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if (r2 < hig2 && *count_found < max_size) {
      /* Gather */
      buffer[*count_found].r2 = r2;
      buffer[*count_found].p = pj;
      buffer[*count_found].xp = xpj;
#ifdef SWIFT_DEBUG_CHECKS
      buffer[*count_found].c = c;
#endif
      (*count_found)++;
    }
    /* TODO: If we reach the max_size, we shoudl print a warning. Maybe add a
       flag to increase the stack size? Or to iterate again later by discarding
       these particles, since they will be tagged as ionized */
  } /* Loop in current cell */
}

/**
 * @brief Gather gas particles within the HII radius from a neighboring cell.
 *
 * Similar to the 'self' version, but handles the distance calculations between
 * two different cells (@p ci and @p cj), including periodic boundary
 * conditions/wrapping.
 *
 * @param r The #runner thread.
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_do_stars_hii_ionization_feedback_pair(
    struct runner *r, struct cell *ci, struct cell *cj, struct spart *si,
    struct hii_neighbor *buffer, int max_size, int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  /* OR h_max_old? */
  const float r_hii_max = ci->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > cj->dmin)
    error("Cell smaller than HII interaction length");

  /* Call the function that will do the work */
  if (interaction_limit > cj->dmin) {
    warning(
        "[c = %lld, cj = %lld, star = %lld] cj size is too small. We need to "
        "go up in the cell hierarchy!",
        ci->cellID, cj->cellID, si->id);
  }

  /* message( */
  /*     "[c = %lld, cj = %lld, cj->super = %lld, star = %lld] Neighbour search,
   * " */
  /*     "interaction_limit = %e, cell_dmin = %e", */
  /*     ci->cellID, cj->cellID, cj->hydro.super->cellID, si->id, */
  /*     interaction_limit, cj->dmin); */
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  const float hi = si->h_hii;
  const float hig2 = hi * hi * kernel_gamma2;
  const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                        (float)(si->x[1] - (cj->loc[1] + shift[1])),
                        (float)(si->x[2] - (cj->loc[2] + shift[2]))};

  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];
    struct xpart *restrict xpj = &xparts_j[pjd];

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error(
              "Particle pj (%lld) not drifted to current time. c = %lld, "
              "c->super = %lld",
              pj->id, cj->cellID, cj->super->cellID);
#endif

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                          (float)(pj->x[1] - cj->loc[1]),
                          (float)(pj->x[2] - cj->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if (r2 < hig2 && *count_found < max_size) {
      /* Gather */
      buffer[*count_found].r2 = r2;
      buffer[*count_found].p = pj;
      buffer[*count_found].xp = xpj;
#ifdef SWIFT_DEBUG_CHECKS
      buffer[*count_found].c = cj;
#endif
      (*count_found)++;
    }
    /* TODO: If we reach the max_size, we shoudl print a warning. Maybe add a
       flag to increase the stack size? Or to iterate again later by discarding
       these particles, since they will be tagged as ionized */
  }
}

/**
 * @brief Sort the gathered HII neighbors by distance (ascending).
 *
 * Uses an insertion sort to order the gas particle candidates stored in the
 * buffer based on their squared distance from the star. This is necessary
 * to ionize gas particles from the inside out.
 *
 * @param buffer The array of #hii_neighbor structures to sort.
 * @param N The number of elements in the buffer.
 */
void runner_sort_hii_neighbors(struct hii_neighbor *buffer, int N) {
  for (int i = 1; i < N; i++) {
    struct hii_neighbor key = buffer[i];
    int j = i - 1;
    while (j >= 0 && buffer[j].r2 > key.r2) {
      buffer[j + 1] = buffer[j];
      j--;
    }
    buffer[j + 1] = key;
  }
}
