/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

void runner_doself_branch_sinks_swallow(struct runner *r, struct cell *c);
void runner_dopair_branch_sinks_swallow(struct runner *r, struct cell *ci,
                                        struct cell *cj);
void runner_dosub_self_sinks_swallow(struct runner *r, struct cell *ci,
                                     int gettimer);
void runner_dosub_pair_sinks_swallow(struct runner *r, struct cell *ci,
                                     struct cell *cj, int gettimer);

void runner_do_sinks_gas_swallow_self(struct runner *r, struct cell *c,
                                      int timer);
void runner_do_sinks_gas_swallow_pair(struct runner *r, struct cell *ci,
                                      struct cell *cj, int timer);

void runner_do_sinks_sink_swallow_self(struct runner *r, struct cell *c,
                                       int timer);
void runner_do_sinks_sink_swallow_pair(struct runner *r, struct cell *ci,
                                       struct cell *cj, int timer);

void runner_do_prepare_part_sink_formation(struct runner *r, struct cell *c,
                                           struct part *restrict p,
                                           struct xpart *restrict xp);
