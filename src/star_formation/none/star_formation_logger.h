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


/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cell.h"
#include "cosmology.h"
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
INLINE static void star_formation_update_stellar_mass(
    struct spart *sp, struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history struct for the stellar mass only
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_stellar_mass(struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history struct if the cell is active
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_SFH_active(struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history struct in the case the cell is inactive
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_init_SFH_inactive(struct star_formation_history *sf) {}

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
 * formation history struct in the #engine
 *
 * @param c the cell of which we want to know the star formation
 * @param sf the star formation structure to which we want to add the star
 * formation
 */
INLINE static void star_formation_get_total_cell(
    struct cell *c, struct star_formation_history *sf) {}

/**
 * @brief add the star formation to the parent cell in the #engine
 *
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_add_to_parent_cell(
    struct cell *c, struct star_formation_history *sf) {}

/**
 * @brief Initialize the star formation history structure in the #engine
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
    struct star_formation_history sf, const char* baseName) {}

/**
 * @brief Initialize the SFH logger file
 *
 * @param none
 */
INLINE static void star_formation_init_file_writer(const char* baseName) {}

/**
 * @brief Add the SFR tracer to the total active SFR of this cell
 *
 * @param p the #part
 * @param xp the #xpart 
 * @param sf the SFH logger struct 
 */
INLINE static void star_formation_log_for_active_particles(
    const struct part* p, const struct xpart* xp, struct star_formation_history *sf, const double dt_star){}

/**
 * @brief Add the SFR tracer to the total inactive SFR of this cell as long as the SFR tracer is larger than 0
 *
 * @param p the #part
 * @param xp the #xpart 
 * @param sf the SFH logger struct 
 */
INLINE static void star_formation_log_for_inactive_particles(
    const struct part* p, const struct xpart* xp, struct star_formation_history *sf){}

/**
 * @brief Initialize the inactive SFR for the cells in the rebuilt 
 *
 * @param p the #part
 * @param xp the #xpart 
 * @param sf the SFH logger struct 
 */
INLINE static void star_formation_SFR_rebuilt(const struct part* p, const struct xpart* xp, struct star_formation_history *sf){}


/**
 * @brief do the recurse after the rebuilt to update the inactive SFR
 *
 * @param p the #part
 * @param xp the #xpart 
 * @param sf the SFH logger struct 
 */
INLINE static void star_formation_recurse_SFR_rebuilt(struct cell *c, const struct cell *cp){}

#endif /* SWIFT_GEAR_STARFORMATION_LOGGER_H */
