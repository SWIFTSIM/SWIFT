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
#include "hydro.h"
#include "part.h"
#include "star_formation_logger_struct.h"

/**
 * @brief Update the stellar mass in the current cell after creating
 * the new star particle spart sp
 *
 * @param sp new created star particle
 * @param sf the star_formation_history struct of the current cell
 */
INLINE static void star_formation_logger_log_new_spart(
    struct spart *sp, struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history struct in the case the cell is
 * inactive
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_logger_log_inactive_cell(
    struct star_formation_history *sf) {}

/**
 * @brief add a star formation history struct to an other star formation history
 * struct
 *
 * @param sf_add the star formation struct which we want to add to the star
 * formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void star_formation_logger_add(
    struct star_formation_history *sf_update,
    const struct star_formation_history *sf_add) {}

/**
 * @brief Initialize the star formation history structure in the #engine
 *
 * @param sfh The pointer to the star formation history structure
 */
INLINE static void star_formation_logger_init(
    struct star_formation_history *sfh) {}

/**
 * @brief Write the final SFH to a file
 *
 * @param fp The file to write to.
 * @param time the simulation time (time since Big Bang) in internal units.
 * @param a the scale factor.
 * @param z the redshift.
 * @param sf the #star_formation_history struct.
 * @param step The time-step of the simulation.
 */
INLINE static void star_formation_logger_write_to_log_file(
    FILE *fp, const double time, const double a, const double z,
    const struct star_formation_history sf, const int step) {}

/**
 * @brief Initialize the SFH logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void star_formation_logger_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {}

/**
 * @brief Add the SFR tracer to the total active SFR of this cell
 *
 * @param p the #part
 * @param xp the #xpart
 * @param sf the SFH logger struct
 * @param dt_star The length of the time-step in physical internal units.
 */
INLINE static void star_formation_logger_log_active_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf, const double dt_star) {}

/**
 * @brief Add the SFR tracer to the total inactive SFR of this cell as long as
 * the SFR tracer is larger than 0
 *
 * @param p the #part
 * @param xp the #xpart
 * @param sf the SFH logger struct
 */
INLINE static void star_formation_logger_log_inactive_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf) {}

#endif /* SWIFT_GEAR_STARFORMATION_LOGGER_H */
