/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_RADIATION_FEEBDACK_H
#define SWIFT_RUNNER_RADIATION_FEEBDACK_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "hydro.h"
#include "stars.h"

/* Avoid cyclic inclusions. */
struct runner;
struct cell;

#define search_radius_factor 1.2f
#define max_ngbs 128
#define max_retry_full_buffer 10

/**
 * @brief Temporary structure to store gas particles found within the HII
 * interaction radius of a star.
 *
 * This is used to gather and sort potential ionization candidates before
 * performing the feedback.
 */
struct hii_neighbor {

  /*! Squared distance between the star and the gas particle. */
  float r2;

  /*! Pointer to the gas particle data. */
  struct part *p;

  /*! Pointer to the gas particle extra data. */
  struct xpart *xp;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Pointer to the cell this particle belongs to */
  struct cell *c;
#endif
};

int runner_hii_check_cell_can_be_reached(const struct cell *ci,
                                         const struct cell *cj, const int sid,
                                         const int flipped,
                                         const double shift[3],
                                         const struct spart *si,
                                         const float search_radius);

void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer);
void runner_dosub_stars_hii_ionization_feedback(struct runner *r,
                                                struct cell *c,
                                                const float interaction_limit);
void runner_do_stars_hii_ionization_feedback_branch(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *ngb_buffer, int max_size,
    int *count_found);

void runner_doself_stars_hii_ionization_feedback(
    struct runner *r, struct cell *c, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found);
void runner_dopair_naive_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const double shift[3],
    struct spart *si, const float search_radius, struct hii_neighbor *buffer,
    int max_size, int *count_found);
void runner_dopair_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found);

/**
 * @brief Maintain a sorted buffer by inserting a new neighbor at the correct
 * position.
 *
 * If the buffer is full, it replaces the furthest element if the new one is
 * closer.
 */
__attribute__((always_inline)) INLINE static void runner_hii_buffer_insert(
    struct hii_neighbor *buffer, int max_size, int *count_found, float r2,
    struct part *p, struct xpart *xp, struct cell *c) {

  /* Case A: Buffer is not yet full */
  if (*count_found < max_size) {
    int i = *count_found - 1;
    /* Shift elements to make room (standard insertion) */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
    (*count_found)++;
  }
  /* Case B: Buffer is full, check if new particle is closer than the furthest
   */
  else if (r2 < buffer[max_size - 1].r2) {
    int i = max_size - 2;
    /* Shift elements to replace the furthest */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
  }
}

/**
 * @brief Verify that the gathered HII neighbor buffer is properly sorted.
 *
 * @param ngb_buffer The array of gathered #hii_neighbor structures.
 * @param count_found The number of valid neighbors stored in the buffer.
 */
__attribute__((always_inline)) INLINE static void
runner_do_stars_hii_ionization_feedback_check_sort(
    const struct hii_neighbor *ngb_buffer, const int count_found) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < count_found - 1; k++) {
    if (ngb_buffer[k].r2 > ngb_buffer[k + 1].r2) {
      error(
          "HII neighbor buffer not properly sorted! "
          "Index %d (r2=%e) is larger than index %d (r2=%e).",
          k, ngb_buffer[k].r2, k + 1, ngb_buffer[k + 1].r2);
    }
  }
#endif
}
#endif /* SWIFT_RUNNER_RADIATION_FEEBDACK_H */
