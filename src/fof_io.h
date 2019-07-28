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
#ifndef SWIFT_FOF_IO_H
#define SWIFT_FOF_IO_H

/* Config parameters. */
#include "../config.h"

#ifdef WITH_FOF

INLINE static void convert_part_group_id(const struct engine* e,
                                         const struct part* p,
                                         const struct xpart* xp,
                                         long long* ret) {
  ret[0] = p->gpart->fof_data.group_id;
}

INLINE static void convert_spart_group_id(const struct engine* e,
                                          const struct spart* sp,
                                          long long* ret) {
  ret[0] = sp->gpart->fof_data.group_id;
}

INLINE static void convert_bpart_group_id(const struct engine* e,
                                          const struct bpart* bp,
                                          long long* ret) {
  ret[0] = bp->gpart->fof_data.group_id;
}

#endif /* WITH_FOF */

/**
 * @brief Specifies which FoF-related particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 *
 * @return The number of fields to write.
 */
INLINE static int fof_write_parts(const struct part* parts,
                                  const struct xpart* xparts,
                                  struct io_props* list) {

#ifdef WITH_FOF

  list[0] = io_make_output_field_convert_part("GroupIDs", LONGLONG, 1,
                                              UNIT_CONV_NO_UNITS, parts, xparts,
                                              convert_part_group_id);
  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Specifies which FoF-related g-particle fields to write to a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to write.
 *
 * @return The number of fields to write.
 */
INLINE static int fof_write_gparts(const struct gpart* gparts,
                                   struct io_props* list) {

#ifdef WITH_FOF

  list[0] = io_make_output_field("GroupIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 gparts, fof_data.group_id);

  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Specifies which FoF-related s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 *
 * @return The number of fields to write.
 */
INLINE static int fof_write_sparts(const struct spart* sparts,
                                   struct io_props* list) {

#ifdef WITH_FOF

  list[0] = io_make_output_field_convert_spart("GroupIDs", LONGLONG, 1,
                                               UNIT_CONV_NO_UNITS, sparts,
                                               convert_spart_group_id);
  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Specifies which FoF-related b-particle fields to write to a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to write.
 *
 * @return The number of fields to write.
 */
INLINE static int fof_write_bparts(const struct bpart* bparts,
                                   struct io_props* list) {

#ifdef WITH_FOF

  list[0] = io_make_output_field_convert_bpart("GroupIDs", LONGLONG, 1,
                                               UNIT_CONV_NO_UNITS, bparts,
                                               convert_bpart_group_id);
  return 1;
#else
  return 0;
#endif
}

#endif /* SWIFT_FOF_IO_H */
