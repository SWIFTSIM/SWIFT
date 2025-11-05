/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_BLACK_HOLES_STRUCT_DEFAULT_H
#define SWIFT_BLACK_HOLES_STRUCT_DEFAULT_H

/**
 * @brief Black holes-related fields carried by each *gas* particle.
 */
struct black_holes_part_data {};

/**
 * @brief Black holes-related fields carried by each *BH* particle.
 */
struct black_holes_bpart_data {};

/**
 * @brief Update a given #part's BH data field to mark the particle has
 * not yet been swallowed.
 *
 * @param p_data The #part's #black_holes_part_data structure.
 */
__attribute__((always_inline)) INLINE static void
black_holes_mark_part_as_not_swallowed(struct black_holes_part_data* p_data) {

  /* Nothing to do here: No swallowing in the default model */
}

/**
 * @brief Reset the particle-carried potential at the start of a time-step.
 *
 * Nothing to do here.
 *
 * @param p_data The #part's black hole data.
 */
__attribute__((always_inline)) INLINE static void black_holes_init_potential(
    struct black_holes_part_data* p_data) {}

/**
 * @brief Update a given #part's BH data field to mark the particle has
 * having been been swallowed.
 *
 * @param p_data The #part's #black_holes_part_data structure.
 */
__attribute__((always_inline)) INLINE static void
black_holes_mark_part_as_swallowed(struct black_holes_part_data* p_data) {

  /* Nothing to do here: No swallowing in the default model */
}

/**
 * @brief Return the ID of the BH that should swallow this #part.
 *
 * @param p_data The #part's #black_holes_part_data structure.
 */
__attribute__((always_inline)) INLINE static long long
black_holes_get_part_swallow_id(struct black_holes_part_data* p_data) {

  /* Return a non-existing ID */
  return -1;
}

/**
 * @brief Update a given #bpart's BH data field to mark the particle has
 * not yet been swallowed.
 *
 * @param p_data The #bpart's #black_holes_bpart_data structure.
 */
__attribute__((always_inline)) INLINE static void
black_holes_mark_bpart_as_not_swallowed(struct black_holes_bpart_data* p_data) {

  /* Nothing to do here: No merging in the default model */
}

/**
 * @brief Update a given #bpart's BH data field to mark the particle has
 * having been been swallowed.
 *
 * @param p_data The #bpart's #black_holes_bpart_data structure.
 */
__attribute__((always_inline)) INLINE static void
black_holes_mark_bpart_as_merged(struct black_holes_bpart_data* p_data) {

  /* Nothing to do here: No merging in the default model */
}

/**
 * @brief Return the ID of the BH that should swallow this #bpart.
 *
 * @param p_data The #bpart's #black_holes_bpart_data structure.
 */
__attribute__((always_inline)) INLINE static long long
black_holes_get_bpart_swallow_id(struct black_holes_bpart_data* p_data) {

  /* Return a non-existing ID */
  return -1;
}

#endif /* SWIFT_BLACK_HOLES_STRUCT_DEFAULT_H */
