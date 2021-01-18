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

/* Include this object's header */
#include "logger_star_formation.h"

/* Local headers  */
#include "logger_tools.h"

const int star_formation_logger_field_size[star_formation_logger_field_count] =
    {
        member_size(struct spart, sf_data.birth_density) +
            member_size(struct spart, sf_data.birth_mass) +
            member_size(struct spart, sf_data.progenitor_id),
};

int star_formation_logger_local_to_global[star_formation_logger_field_count];
