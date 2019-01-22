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
#ifndef SWIFT_VELOCIRAPTOR_IO_H
#define SWIFT_VELOCIRAPTOR_IO_H

/* Config parameters. */
#include "../config.h"

INLINE static void velociraptor_convert_part_groupID(const struct engine* e,
                                                     const struct part* p,
                                                     const struct xpart* xp,
                                                     long long* ret) {
  if (p->gpart == NULL)
    ret[0] = 0.f;
  else {
    const ptrdiff_t offset = p->gpart - e->s->gparts;
    *ret = (e->s->gpart_group_data + offset)->groupID;
  }
}

INLINE static void velociraptor_convert_spart_groupID(const struct engine* e,
                                                      const struct spart* sp,
                                                      long long* ret) {
  if (sp->gpart == NULL)
    ret[0] = 0.f;
  else {
    const ptrdiff_t offset = sp->gpart - e->s->gparts;
    *ret = (e->s->gpart_group_data + offset)->groupID;
  }
}

__attribute__((always_inline)) INLINE static int velociraptor_write_parts(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS, parts, xparts,
      velociraptor_convert_part_groupID);

  return 1;
}

__attribute__((always_inline)) INLINE static int velociraptor_write_gparts(
    const struct velociraptor_gpart_data* group_data, struct io_props* list) {

  list[0] = io_make_output_field("GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 group_data, groupID);

  return 1;
}

__attribute__((always_inline)) INLINE static int velociraptor_write_sparts(
    const struct spart* sparts, struct io_props* list) {

  list[0] = io_make_output_field_convert_spart(
      "GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS, sparts,
      velociraptor_convert_spart_groupID);

  return 1;
}

#endif /* SWIFT_VELOCIRAPTOR_IO_H */
