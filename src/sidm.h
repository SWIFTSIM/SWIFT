/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_SIDM_H
#define SWIFT_SIDM_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "inline.h"
#include "part.h"
#include "sidm_properties.h"
#include "sidm_iact.h"

/**
 * @brief Sets the SIDM properties of the g-particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void sidm_first_init_gpart(struct gpart* restrict gp,
                                                                        const struct sidm_props* sidm_props) {
    
    /*! Flag to indicate the particle has been scattered yes(1)/no(0) */
    gp->sidm_data.sidm_flag = 0.0f;
    
    /*! Particle search radius */
    gp->sidm_data.h_sidm = sidm_props->h_search_radius;
    
    /*! Number of DM-DM particle collisions */
    gp->sidm_data.num_sidm = 0.0f;
    
    /* Particle velocity */
    gp->sidm_data.si_v_full[0] = 0.0f;
    gp->sidm_data.si_v_full[1] = 0.0f;
    gp->sidm_data.si_v_full[2] = 0.0f;
    
    gp->sidm_data.test_flag = 0.0f;
    
    const double v2 = gp->v_full[0] * gp->v_full[0] + gp->v_full[1] * gp->v_full[1] + gp->v_full[2] * gp->v_full[2];
    if (v2 > 0.001) gp->sidm_data.test_flag = 1.0f; /* 1 indicates parts in cube */
}

#endif /* SWIFT_SIDM_H */
