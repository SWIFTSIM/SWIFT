/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Josh Borrow (josh@joshborrow.com)
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

#ifndef SWIFT_PARTICLE_SPLITTING_H
#define SWIFT_PARTICLE_SPLITTING_H

#include "inline.h"
#include "io_properties.h"
#include "particle_splitting_struct.h"

#include <stdint.h>

/**
 * @brief Initialise a particle_splitting_data struct
 *        at the start of a run, given an initial
 *        progenitor ID.
 *
 * @param splitting_data the uninitialised particle struct.
 * @param id the ID of the particle (used for progenitor_id)
 */
__attribute__((always_inline)) INLINE static void
particle_splitting_mark_part_as_not_split(
    struct particle_splitting_data* restrict splitting_data, long long id) {

  splitting_data->progenitor_id = id;
  splitting_data->split_tree = 0;
  splitting_data->split_count = 0;
}

/**
 * @brief Updates the binary trees of the particles that
 *        have been split. pi should retain its original ID,
 *        and hence gets the relevant part of the tree set
 *        to zero.
 *
 * @param sdi first particle_splitting_data* resulting from
 *           the splitting event.
 * @param sdj second particle_splitting_data* resulting from
 *           the splitting event.
 * @param id_i ID of the first particle.
 * @param id_j ID of the first particle.
 * @param extra_split_logger File pointer (opened) where to log
 *           the resetting of the split counters.
 */
__attribute__((always_inline)) INLINE static void
particle_splitting_update_binary_tree(
    struct particle_splitting_data* restrict sdi,
    struct particle_splitting_data* restrict sdj, const long long id_i,
    const long long id_j, FILE* extra_split_logger,
    swift_lock_type* file_lock) {

  /* Print warnings if we have split these particles more
   * than the number of times the tree can accommodate.
   * Warning is only printed once for each particle */
  if (sdi->split_count > 0 &&
      sdi->split_count % (8 * sizeof(sdi->split_tree)) == 0) {
    message(
        "Warning: Particle (%lld) with progenitor ID %lld with binary tree "
        "%lld has been split over the maximum %zu times, making its binary "
        "tree invalid.",
        id_i, sdi->progenitor_id, sdi->split_tree, sizeof(sdi->split_tree));

    /* Are we logging this? */
    if (extra_split_logger != NULL) {

      /* Log the old state before reseting. Use the lock to prevent multiple
       * threads from writing at the same time. */
      lock_lock(file_lock);
      fprintf(extra_split_logger, "  %12d %20lld %20lld %20d %20lld\n",
              engine_current_step, id_i, sdi->progenitor_id, sdi->split_count,
              sdi->split_tree);
      fflush(extra_split_logger);

      /* Release the lock and continue in parallel */
      if (lock_unlock(file_lock) != 0)
        error("Impossible to unlock particle splitting");
    }

    /* Reset both counters and trees */
    sdi->split_tree = 0LL;
    sdj->split_tree = 0LL;

    /* Set both particles as having particle i as their progenitor */
    sdi->progenitor_id = id_i;
    sdj->progenitor_id = id_i;
  }

  /* Update the binary tree */
  sdj->split_tree |= 1LL << sdj->split_count % (8 * sizeof(sdi->split_tree));

  /* Increase counters on both; sdi implicitly has a zero
   * in the relevant spot in its binary tree */
  sdj->split_count++;
  sdi->split_count++;
}

/**
 * @brief Specifies the particle-splitting related fields
 *        to write to a dataset.
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int particle_splitting_write_particles(const struct part* parts,
                                                     const struct xpart* xparts,
                                                     struct io_props* list,
                                                     const int with_cosmology) {

  list[0] = io_make_output_field(
      "ProgenitorParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      split_data.progenitor_id,
      "ID of the progenitor of this particle. If this particle is the result "
      "of one (or many) splitting events, this ID corresponds to the ID of the "
      "particle in the initial conditions that its lineage can be traced back "
      "to. If the particle was never split, this is the same as ParticleIDs.");

  list[1] = io_make_output_field(
      "SplitCounts", UINT8, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      split_data.split_count,
      "Number of times this particle has been split. Note that both particles "
      "that take part in the splitting have their counter incremented.");

  list[2] = io_make_output_field(
      "SplitTrees", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      split_data.split_tree,
      "Binary tree describing splitting events. Particles that keep the "
      "original ID have a value of zero in a splitting event, whereas"
      "particles given a new ID have a value of one.");

  return 3;
}

/**
 * @brief Specifies which star particle fields to write to a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int particle_splitting_write_sparticles(
    const struct spart* sparts, struct io_props* list) {

  list[0] = io_make_output_field(
      "ProgenitorParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      split_data.progenitor_id,
      "Progenitor ID of the gas particle that became this star. If this "
      "particle is the result of one (or many) splitting events, this ID "
      "corresponds to the ID of the particle in the initial conditions that "
      "its lineage can be traced back to. If the particle was never split, "
      "this is the same as ParticleIDs.");

  list[1] = io_make_output_field(
      "SplitCounts", UINT8, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      split_data.split_count,
      "Number of times the gas particle that turned into this star particle "
      "was split. Note that both particles that take part in the splitting "
      "have their counter incremented.");

  list[2] = io_make_output_field(
      "SplitTrees", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      split_data.split_tree,
      "Binary tree describing splitting events. Particles that keep the "
      "original ID have a value of zero in a splitting event, whereas"
      "particles given a new ID have a value of one.");

  return 3;
}

/**
 * @brief Specifies which black hole particle fields to write to a dataset
 *
 * @param bparts The black hole particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int particle_splitting_write_bparticles(
    const struct bpart* bparts, struct io_props* list) {
  list[0] = io_make_output_field(
      "ProgenitorParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      split_data.progenitor_id,
      "Progenitor ID of the gas particle that became the seed BH. If this "
      "particle is the result of one (or many) splitting events, this ID "
      "corresponds to the ID of the particle in the initial conditions that "
      "its lineage can be traced back to. If the particle was never split, "
      "this is the same as ParticleIDs.");

  list[1] = io_make_output_field(
      "SplitCounts", UINT8, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      split_data.split_count,
      "Number of times the gas particle that became this BH seed "
      "was split. Note that both particles that take part in the splitting "
      "have their counter incremented.");

  list[2] = io_make_output_field(
      "SplitTrees", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      split_data.split_tree,
      "Binary tree describing splitting events prior to BH seeding. Particles "
      "that keep the original ID have a value of zero in a splitting event, "
      "whereas particles given a new ID have a value of one.");

  return 3;
}

/**
 * @brief Write particle splitting data to the stdout for debugging purposes.
 *
 * @param p Particle data.
 * @param xp Extra particle data.
 */
__attribute__((always_inline)) INLINE static void
particle_splitting_debug_particle(const struct part* p,
                                  const struct xpart* xp) {

  if (xp != NULL) {
    warning("[PID%lld] particle_splitting_data:", p->id);
    warning(
        "[PID%lld] progenitor_id = %lli, split_tree = %lli, "
        "split_count = %hhu",
        p->id, xp->split_data.progenitor_id, xp->split_data.split_tree,
        xp->split_data.split_count);
  }
}

#endif /* SWIFT_PARTICLE_SPLITTING_H */
