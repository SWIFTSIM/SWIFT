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
#ifndef SWIFT_SIDM_PROPERTIES_H
#define SWIFT_SIDM_PROPERTIES_H

/**
 * @file sidm_properties.h
 * @brief Contains all the constants and parameters of the SIDM scheme
 */

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Forward declarations */
struct phys_const;
struct unit_system;
struct swift_params;
struct cosmology;

/**
 * @brief Properties of the Self-interacting dark matter model.
 */
struct sidm_props {
    
    /* ------------ Main operation modes ------------- */
    
    /*! Are we using velocity-dependent scattering? */
    /*int with_sigma_velocity_function;*/
    
    /*! Are we using constant scattering? */
    /*int with_sigma_constant;*/

    /* Scattering cross section (in physical units: cm^2/g) */
    double sigma_cgs;

    /* Scattering cross section (in internal units) */
    double sigma;
    
    float h_search_radius;
    
};

/**
 * @brief extra particle data for #gpart in the SIDM model.
 */
struct sidm_gpart_data {

    /*! Velocity changed due to DM-DM self-interactions. */
    float si_v_full[3];
    
    /*! flag indicating if particle in given time-step has been scattered*/
    float sidm_flag;
    
    float test_flag;
    
    /*! Particle search radius */
    float h_sidm;
    
    /*! Number of DM-DM collisions */
    float num_sidm;
};

void sidm_props_init(struct sidm_props* sidm_props,
                     const struct phys_const* phys_const,
                     const struct unit_system* us,
                     struct swift_params* params,
                     const struct cosmology* cosmo);

/* Dump/restore. */
/*void sidm_props_struct_dump(const struct sidm_props *p, FILE *stream);
void sidm_props_struct_restore(struct sidm_props *p, FILE *stream);
*/

#endif /* SWIFT_SIDM_PROPERTIES_H */
