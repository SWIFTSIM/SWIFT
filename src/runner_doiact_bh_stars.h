/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Robert McGibbon (robjmcgibbon@gmail.com)
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
#ifndef SWIFT_RUNNER_DOIACT_BH_STARS_H
#define SWIFT_RUNNER_DOIACT_BH_STARS_H

/* Config parameters. */
#include <config.h>

struct runner;
struct cell;
struct bpart;

/* Loops computing the density of star particles around black holes.
 *
 * The black holes gather from the stars only (non-symmetric); the star
 * particles are never updated. The black holes are not sorted (they are
 * rare), but the pair interactions walk the star sort arrays where
 * available. */

void runner_doself_bh_stars_density(struct runner *r, struct cell *c,
                                    int timer);
void runner_dopair_bh_stars_density(struct runner *r, struct cell *ci,
                                    struct cell *cj, int timer);

void runner_doself_branch_bh_stars_density(struct runner *r, struct cell *c);
void runner_dopair_branch_bh_stars_density(struct runner *r, struct cell *ci,
                                           struct cell *cj);

void runner_doself_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts,
                                                  int *ind, const int bcount);
void runner_dopair_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts_i,
                                                  int *ind, const int bcount,
                                                  struct cell *cj);
void runner_dosub_subset_bh_stars_density(struct runner *r, struct cell *ci,
                                          struct bpart *bparts, int *ind,
                                          const int bcount, struct cell *cj,
                                          int gettimer);

void runner_dosub_self_bh_stars_density(struct runner *r, struct cell *c,
                                        int gettimer);
void runner_dosub_pair_bh_stars_density(struct runner *r, struct cell *ci,
                                        struct cell *cj, int gettimer);

#endif /* SWIFT_RUNNER_DOIACT_BH_STARS_H */
