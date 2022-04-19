/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
#include "../config.h"

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "sink.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"




/**
 * @brief Calculate the number density of #part around the #bpart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_doself_sinks_swallow(struct runner *r, struct cell *c, int timer) {


}

/**
 * @brief Calculate the number density of cj #part around the ci #bpart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void runner_do_nonsym_pair_sinks_naive_swallow(struct runner *r, struct cell *restrict ci,
                              struct cell *restrict cj) {

}

void runner_dopair_sinks_naive_swallow(struct runner *r, struct cell *restrict ci,
                      struct cell *restrict cj, int timer) {

}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts_i The #bpart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void runner_dopair_subset_sinks_naive_swallow(struct runner *r, struct cell *restrict ci,
                             struct bpart *restrict bparts_i, int *restrict ind,
                             const int bcount, struct cell *restrict cj,
                             const double *shift) {

}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 */
void runner_doself_subset_sinks_swallow(struct runner *r, struct cell *restrict ci,
                       struct bpart *restrict bparts, int *restrict ind,
                       const int bcount) {

}

/**
 * @brief Determine which version of DOSELF1_SUBSET_BH needs to be called
 * depending on the optimisation level.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 */
void runner_doself_subset_branch_sinks_swallow(struct runner *r, struct cell *restrict ci,
                              struct bpart *restrict bparts, int *restrict ind,
                              const int bcount) {

}

/**
 * @brief Determine which version of DOPAIR1_SUBSET_BH needs to be called
 * depending on the orientation of the cells or whether DOPAIR1_SUBSET_BH
 * needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts_i The #bpart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void runner_dopair_subset_branch_sinks_swallow(struct runner *r, struct cell *restrict ci,
                              struct bpart *restrict bparts_i,
                              int *restrict ind, int const bcount,
                              struct cell *restrict cj) {

}

void runner_dosub_subset_sinks_swallow(struct runner *r, struct cell *ci, struct bpart *bparts,
                     int *ind, const int bcount, struct cell *cj,
                     int gettimer) {


}

/**
 * @brief Determine which version of DOSELF1_BH needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void runner_doself_branch_sinks_swallow(struct runner *r, struct cell *c) {


}

/**
 * @brief Determine which version of DOPAIR1_BH needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_BH needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void runner_dopair_branch_sinks_swallow(struct runner *r, struct cell *ci, struct cell *cj) {


}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void runner_dosub_pair_sinks_swallow(struct runner *r, struct cell *ci, struct cell *cj,
                    int gettimer) {


}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_self_sinks_swallow(struct runner *r, struct cell *ci, int gettimer) {


}




/////////////////////////////////////////////////////////////////////////////









/**
 * @brief Process all the gas particles in a cell that have been flagged for
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the gas particle and the gas particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the gas particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow(struct runner *r, struct cell *c, int timer) {

}

/**
 * @brief Processing of gas particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_self(struct runner *r, struct cell *c, int timer) {


}

/**
 * @brief Processing of gas particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_pair(struct runner *r, struct cell *ci,
                                struct cell *cj, int timer) {

 
}

/**
 * @brief Process all the BH particles in a cell that have been flagged for
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the BH particle and the BH particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the BH particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_bh_swallow(struct runner *r, struct cell *c, int timer) {

}

/**
 * @brief Processing of bh particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_bh_swallow_self(struct runner *r, struct cell *c, int timer) {


}

/**
 * @brief Processing of bh particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_bh_swallow_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer) {


}
