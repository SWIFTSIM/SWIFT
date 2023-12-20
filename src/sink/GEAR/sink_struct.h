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
#include <sys/types.h>
struct sink_part_data {

  /*! ID of the sink that will swallow this #part. */
  long long swallow_id;

  /*! Gravitational potential of the particle */
  uint8_t can_form_sink;

  /* Total kinetic energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_kin_neighbours;

  /* Total interal energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_int_neighbours;

  /* Total radiated energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_rad_neighbours;

  /* Total self potential energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_pot_self_neighbours;

  /* Total external potential energy of the neigbouring gas particles (i.e.
   * inside sink_cut_off_radius) */
  double E_pot_ext_neighbours;

  /* Total magnetic energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_mag_neighbours;

  /* Total rotational energy of the neigbouring gas particles (i.e. inside
   * sink_cut_off_radius) */
  double E_rot_neighbours;

  /* Potential of the particle copied from the #gpart */
  float potential;

  /* Mechanical energy between the part and the sink with swallow_id.
   * This is used to check that this part is, out of all sinks, the most bound to
     the sink with swallow_id. */
  double E_mec_bound;

  /* Does the future sink overalp an existing one ? */
  uint8_t is_overlapping_sink;
};

/**
 * @brief Sink-related fields carried by each *sink* particle.
 */
struct sink_sink_data {

  /*! ID of the sink that will swallow this #sink. */
  long long swallow_id;

  /*! Mass of the sink that will swallow this #sink. */
  float swallow_mass;
};

/**
 * @brief Return the ID of the sink that should swallow this #part.
 *
 * @param s_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static long long sink_get_part_swallow_id(
    struct sink_part_data* s_data) {

  return s_data->swallow_id;
}

/**
 * @brief Update a given #part's sink data field to mark the particle has
 * not yet been swallowed.
 *
 * @param s_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static void
sink_mark_part_as_not_swallowed(struct sink_part_data* s_data) {

  s_data->swallow_id = -1;
}

/**
 * @brief Update a given #part's sink data field to mark the particle has
 * having been been swallowed.
 *
 * @param p_data The #part's #sink_part_data structure.
 */
__attribute__((always_inline)) INLINE static void sink_mark_part_as_swallowed(
    struct sink_part_data* s_data) {

  s_data->swallow_id = -2;
}

/**
 * @brief Update a given #sink's sink data field to mark the particle has
 * not yet been swallowed.
 *
 * @param s_data The #sink's #sink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static void
sink_mark_sink_as_not_swallowed(struct sink_sink_data* s_data) {

  s_data->swallow_id = -1;
  s_data->swallow_mass = 0.f;
}

/**
 * @brief Update a given #sink's sink data field to mark the particle has
 * having been been swallowed.
 *
 * @param s_data The #sink's #bsink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static void sink_mark_sink_as_merged(
    struct sink_sink_data* s_data) {

  s_data->swallow_id = -2;
  s_data->swallow_mass = -1.f;
}

/**
 * @brief Return the ID of the sink that should swallow this #sink.
 *
 * @param s_data The #sink's #sink_sink_data structure.
 */
__attribute__((always_inline)) INLINE static long long sink_get_sink_swallow_id(
    struct sink_sink_data* s_data) {

  return s_data->swallow_id;
}

#endif /* SWIFT_SINK_STRUCT_DEFAULT_H */
