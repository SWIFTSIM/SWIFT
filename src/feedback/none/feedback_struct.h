/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_NONE_H
#define SWIFT_FEEDBACK_STRUCT_NONE_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each hydro particles
 */
struct feedback_part_data {};

/**
 * @brief Report-back payload for task_subtype_part_hii_tag. Empty stub:
 * this module has no HII radiation, this type exists only so core files
 * (scheduler.c, cell_pack.c) compile against the subtype unconditionally.
 */
struct hii_tag_report {};

/**
 * @brief Post-cooling state payload for task_subtype_part_hii_state. Empty
 * stub, same rationale as hii_tag_report above.
 */
struct hii_state_update {};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {};

/**
 * @brief Feedback fields carried by each star particles
 *
 * Nothing here since this is a no-feedback model.
 */
struct feedback_spart_data {};

#endif /* SWIFT_FEEDBACK_STRUCT_NONE_H */
