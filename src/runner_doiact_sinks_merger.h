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

#define PASTE(x, y) x##_##y

#define _DOSELF1_SINKS_MERGER(f) PASTE(runner_doself_sinks, f)
#define DOSELF1_SINKS_MERGER _DOSELF1_SINKS_MERGER(FUNCTION)

#define _DO_SYM_PAIR1_SINKS_MERGER(f) PASTE(runner_do_sym_pair_sinks, f)
#define DO_SYM_PAIR1_SINKS_MERGER _DO_SYM_PAIR1_SINKS_MERGER(FUNCTION)

#define _DOSUB_PAIR1_SINKS_MERGER(f) PASTE(runner_dosub_pair_sinks, f)
#define DOSUB_PAIR1_SINKS_MERGER _DOSUB_PAIR1_SINKS_MERGER(FUNCTION)

#define _DOSUB_SELF1_SINKS_MERGER(f) PASTE(runner_dosub_self_sinks, f)
#define DOSUB_SELF1_SINKS_MERGER _DOSUB_SELF1_SINKS_MERGER(FUNCTION)

#define _IACT_SINK_MERGER(f) PASTE(runner_iact_sym_sinks, f)
#define IACT_SINK_MERGER _IACT_SINK_MERGER(FUNCTION)

void DOSELF1_SINKS_MERGER(struct runner *r, struct cell *c);
void DO_SYM_PAIR1_SINKS_MERGER(struct runner *r, struct cell *ci,
                               struct cell *cj);

void DOSUB_SELF1_SINKS_MERGER(struct runner *r, struct cell *ci);
void DOSUB_PAIR1_SINKS_MERGER(struct runner *r, struct cell *ci,
                              struct cell *cj);
