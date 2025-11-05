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
#ifndef SWIFT_SINK_STRUCT_DEFAULT_H
#define SWIFT_SINK_STRUCT_DEFAULT_H

/**
 * @brief Sink-related fields carried by each *gas* particle.
 */
struct sink_part_data {};

/**
 * @brief Sink-related fields carried by each *sink* particle.
 */
struct sink_sink_data {};

/**
 * @brief Return the ID of the sink that should swallow this #part.
 *
 * @param s_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static long long sink_get_part_swallow_id(
    struct sink_part_data* s_data) {

  /* Return a non-existing ID */
  return -1;
}

/**
 * @brief Update a given #part's sink data field to mark the particle has
 * not yet been swallowed.
 *
 * @param s_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static void
sink_mark_part_as_not_swallowed(struct sink_part_data* s_data) {}

/**
 * @brief Update a given #part's sink data field to mark the particle has
 * having been been swallowed.
 *
 * @param p_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static void sink_mark_part_as_swallowed(
    struct sink_part_data* s_data) {}

/**
 * @brief Update a given #sink's sink data field to mark the particle has
 * not yet been swallowed.
 *
 * @param s_data The #sink's #sink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static void
sink_mark_sink_as_not_swallowed(struct sink_sink_data* s_data) {}

/**
 * @brief Update a given #sink's sink data field to mark the particle has
 * having been been swallowed.
 *
 * @param s_data The #sink's #bsink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static void sink_mark_sink_as_merged(
    struct sink_sink_data* s_data) {}

/**
 * @brief Return the ID of the sink that should swallow this #sink.
 *
 * @param s_data The #sink's #sink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static long long sink_get_sink_swallow_id(
    struct sink_sink_data* s_data) {

  /* Return a non-existing ID */
  return -1;
}

#endif /* SWIFT_SINK_STRUCT_DEFAULT_H */
