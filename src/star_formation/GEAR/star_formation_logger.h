/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_STARFORMATION_LOGGER_H
#define SWIFT_GEAR_STARFORMATION_LOGGER_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cell.h"
#include "cosmology.h"
#include "hydro.h"
#include "part.h"
#include "star_formation_logger_struct.h"

/**
 * @brief Update the star foramtion history in the current cell after creating
 * the new star particle spart sp
 *
 * @param sp new created star particle
 * @param sf the star_formation_history struct of the current cell
 */
INLINE static void star_formation_update_SFH(
    struct spart *sp, struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history struct
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_SFH(struct star_formation_history *sf) {}

/**
 * @brief function to add the progeny SFH to the parent SFH.
 *
 * @param sf parent SFH struct
 * @param sfprogeny progeny SFH struct
 */
INLINE static void star_formation_add_progeny_SFH(
    struct star_formation_history *sf,
    const struct star_formation_history *sfprogeny) {}

/**
 * @brief Get the total star formation in this cell and add it to the star
 * formation history struct
 *
 * @param c the cell of which we want to know the star formation
 * @param sf the star formation structure to which we want to add the star
 * formation
 * @param cosmo the cosmology struct
 * @param with_cosmology if we run with cosmology
 */
INLINE static void star_formation_get_total_cell(
    struct cell *c, struct star_formation_history *sf) {}

/**
 * @brief Clear the total star formation in this cell
 *
 * @param c the cell of which we want to know the star formation
 */
INLINE static void star_formation_clear_total_cell(struct cell *c) {}

/**
 * @brief add the star formation to the parent cell
 *
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_add_to_parent_cell(
    struct cell *c, struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history structure
 *
 * @param The pointer to the star formation history structure
 * */
INLINE static void star_formation_init_SFH_engine(
    struct star_formation_history *sfh) {}

/**
 * @brief Write the final SFH to a file
 *
 * @param time the simulation time
 * @param a the scale factor
 * @param z the redshift
 * @param sf the star_formation_history struct
 */
INLINE static void star_formation_write_to_file(
    const double time, const double a, const double z,
    struct star_formation_history sf) {}

/**
 * @brief Initialize the SFH logger file
 *
 * @param none
 */
INLINE static void star_formation_init_file_writer(void) {}

#endif /* SWIFT_GEAR_STARFORMATION_LOGGER_H */
