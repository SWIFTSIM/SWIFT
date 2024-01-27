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
#ifndef SWIFT_DEFAULT_DARK_MATTER_PROPERTIES_H
#define SWIFT_DEFAULT_DARK_MATTER_PROPERTIES_H

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
struct gravity_props;

/**
 * @brief Properties of the Self-interacting dark matter model.
 */
struct sidm_props {
    
    /* ------------ Main operation modes ------------- */
    
    /*! Are we using velocity-dependent scattering? */
    int with_velocity_dependent_sigma;
    int with_momentum_transfer_sigma;
    
    /* If using velocity-dependent scattering, give values */
    /* for dark matter particle mass and mediator mass */
    /* (in physical units: MeV) */
    double mx;
    double mphi;
    double alphax;

    /*! Are we using constant scattering? */
    int with_constant_sigma;
    
    /*! Are we using angle dependency when particles collide? */
    int with_isotropic_scattering;
    int with_anisotropic_scattering;

    /* Scattering cross section (in physical units: cm^2/g) */
    double sigma_cgs;

    /* Scattering cross section (in internal units) */
    double sigma;
    
    /* Parameters for smoothing length calculation */
    /*! Maximal smoothing length (internal units) */
    float h_max;
    
    /*! Minimal smoothing length (internal units) */
    float h_min;
    
    /*! Smoothing length tolerance */
    float h_tolerance;
    
    /*! Minimal smoothing length expressed as ratio to softening length */
    float h_min_ratio;
    
    /*! Resolution parameter */
    float eta_neighbours;

    /*! Target neighbours */
    float target_neighbours;

    /*!  */
    float delta_neighbours;

    /*! Maximal number of iterations to converge h */
    int max_smoothing_iterations;
    
    /* ------ Neighbour number definition ------------ */
    
    /*! Are we using the mass-weighted definition of neighbour number? */
    int use_mass_weighted_num_ngb;
    
    /* ------ Time integration parameters ------------ */
    
    /*! Maximal change of h over one time-step */
    float log_max_h_change;

    /*! -- Minimum allowed time-step of DMpart in internal units */
    float time_step_min;

    /*! -- Maximum allowed time-step of DMpart in internal units */
    float time_step_max;

};

/**
 * @brief extra particle data for #gpart in the SIDM model.
 */
struct sidm_dmpart_data {

    /*! Velocity changed due to DM-DM collisions. */
    float v_full[3];
    
    /*! Flag indicating if particle in given time-step has been scattered */
    float sidm_flag;
    
    /* Cross section */
    double sigma;
    
    /*! Cumulative number of DM-DM collisions */
    float number_of_sidm_events;
    
    /* Maximum number of SIDM events per timestep */
    float max_sidm_events_per_timestep;
    
    /* Counter of SIDM events per timestep */
    float sidm_events_per_timestep;

    /*! Time step size of particle's drift */
    float dt_drift;
};

void sidm_props_print(struct sidm_props* sidm_props);

void sidm_props_init(struct sidm_props* sidm_props,
                     const struct phys_const* phys_const,
                     const struct unit_system* us,
                     struct swift_params* params,
                     const struct cosmology* cosmo);

void sidm_props_update(struct sidm_props* sidm_props,
                    const struct gravity_props *gp,
                    const struct cosmology *cosmo);


#if defined(HAVE_HDF5)
void sidm_props_print_snapshot(hid_t h_grpsph, const struct sidm_props *p);
#endif

/* Dump/restore. */
void sidm_props_struct_dump(const struct sidm_props *p, FILE *stream);
void sidm_props_struct_restore(const struct sidm_props *p, FILE *stream);

#endif /* SWIFT_DARK_MATTER_PROPERTIES_H */
