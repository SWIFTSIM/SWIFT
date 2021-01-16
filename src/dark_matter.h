/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_DARK_MATTER_H
#define SWIFT_DARK_MATTER_H

/* Config parameters. */
#include "../config.h"
#include "approx_math.h"
#include "./dark_matter_part.h"
#include "./dark_matter_iact.h"

/**
 * @brief Prepares a dm-particle for its interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void dark_matter_init_dmpart(struct dmpart* gp) {
    
    gp->density.wcount = 0.f;
    gp->density.wcount_dh = 0.f;
    gp->rho = 0.f;
    gp->density.rho_dh = 0.f;
    
    gp->sidm_probability = 0.f;
    gp->time_step_size = 0.f;
    gp->num_neighbours = 0.f;
}


/**
 * @brief Prepares a dm-particle for its SIDM interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_velocities(struct dmpart* gp, float dt_drift) {
    
    /* No SIDM flag at the beginning */
    gp->sidm_data.sidm_flag = 0.0f;
    
    /* Reset counter */
    gp->sidm_data.sidm_events_per_timestep = 0.0f;

    /* Drift time */
    gp->sidm_data.dt_drift = dt_drift;
    
    /* Make a copy of particle velocity */
    gp->sidm_data.v_full[0] = gp->v_full[0];
    gp->sidm_data.v_full[1] = gp->v_full[1];
    gp->sidm_data.v_full[2] = gp->v_full[2];

}


/**
 * @brief Initialises the dark matter particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param gp The particle to act upon.
 * @param sidm_properties Properties of the self-interacting dark matter model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_first_init_dmpart(
     struct dmpart* dmp, const struct sidm_props* sidm_props) {
    
    /*! Flags to indicate if the particle has been scattered yes(1)/no(0) */
    dmp->sidm_data.sidm_flag = 0.0f;
    
    /*! Number of DM-DM particle collisions */
    dmp->sidm_data.number_of_sidm_events = 0.0f;
    dmp->sidm_data.sidm_events_per_timestep = 0.0f;
    dmp->sidm_data.max_sidm_events_per_timestep = 0.0f;

    /*! Fixed search radius */
    dmp->h_sidm = dmp->h;
    
    /* Particle velocity */
    dmp->sidm_data.v_full[0] = dmp->v_full[0];
    dmp->sidm_data.v_full[1] = dmp->v_full[1];
    dmp->sidm_data.v_full[2] = dmp->v_full[2];

    dmp->sidm_data.dt_drift = 0.0f;

    dmp->time_bin = 0;

    dark_matter_init_dmpart(dmp);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param gp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void dark_matter_predict_extra(
               struct dmpart* restrict dmp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param gp The particle.
 */
__attribute__((always_inline)) INLINE static void dark_matter_reset_predicted_values(
                struct dmpart* restrict gp) {}


/**
 * @brief Finishes the calculation of density on dark matter particles
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_end_density(
    struct dmpart* gp, const struct cosmology* cosmo, const struct sidm_props* sidm_props,
    const double dt) {
        
    /* Some smoothing length multiples. */
    const float h = gp->h;
    const float h_inv = 1.0f / h;                       /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
    const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    
    /* Final operation on the density (add self-contribution). */
    gp->rho += gp->mass * dm_kernel_root;
    gp->density.rho_dh -= hydro_dimension * gp->mass * dm_kernel_root;
    gp->density.wcount += dm_kernel_root;
    gp->density.wcount_dh -= hydro_dimension * dm_kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    gp->rho *= h_inv_dim;
    gp->density.rho_dh *= h_inv_dim_plus_one;
    gp->density.wcount *= h_inv_dim;
    gp->density.wcount_dh *= h_inv_dim_plus_one;
    
    gp->time_step_size = dt;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_part_has_no_neighbours(
    struct dmpart* restrict gp, const struct cosmology* cosmo) {
    
    /* Some smoothing length multiples. */
    const float h = gp->h;
    const float h_inv = 1.0f / h;                 /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

    /* Re-set problematic values */
    gp->rho = gp->mass * dm_kernel_root * h_inv_dim;
    gp->density.wcount = dm_kernel_root * h_inv_dim;
    gp->density.rho_dh = 0.f;
    gp->density.wcount_dh = 0.f;
    gp->num_neighbours = 0.f;
}



/**
 * @brief Updates #dmparts velocities
 *
 * @param dmp #dmpart
 *
 */
__attribute__((always_inline)) INLINE static void sidm_kick_to_dmpart(struct dmpart *restrict dmp) {
    
    if (dmp->sidm_data.sidm_flag > 3) {
        
        /* Reverse recent half drift */
        dmp->x[0] -= dmp->v_full[0] * dmp->sidm_data.dt_drift / 2.f;
        dmp->x[1] -= dmp->v_full[1] * dmp->sidm_data.dt_drift / 2.f;
        dmp->x[2] -= dmp->v_full[2] * dmp->sidm_data.dt_drift / 2.f;

        /* Update velocity due to collision */
        dmp->v_full[0] = dmp->sidm_data.v_full[0];
        dmp->v_full[1] = dmp->sidm_data.v_full[1];
        dmp->v_full[2] = dmp->sidm_data.v_full[2];
        
        /* Now add half drift with new velocity */
        dmp->x[0] += dmp->v_full[0] * dmp->sidm_data.dt_drift / 2.f;
        dmp->x[1] += dmp->v_full[1] * dmp->sidm_data.dt_drift / 2.f;
        dmp->x[2] += dmp->v_full[2] * dmp->sidm_data.dt_drift / 2.f;

        /* Get its gravity friend */
        struct gpart *gp = dmp->gpart;
        
        /* Synchronize particle's velocities and positions */
        gp->v_full[0] = dmp->v_full[0];
        gp->v_full[1] = dmp->v_full[1];
        gp->v_full[2] = dmp->v_full[2];

        gp->x[0] = dmp->x[0];
        gp->x[1] = dmp->x[1];
        gp->x[2] = dmp->x[2];

        /* Remove flag */
        dmp->sidm_data.sidm_flag = 0.f;
        
        /* Update counters */
        if (dmp->sidm_data.max_sidm_events_per_timestep < dmp->sidm_data.sidm_events_per_timestep)
            dmp->sidm_data.max_sidm_events_per_timestep = dmp->sidm_data.sidm_events_per_timestep;
        
    }
}


/**
 * @brief Computes the dark matter time-step of a given particle
 *
 * This function returns the time-step of a particle given its DM-dynamical
 * state?
 *
 * @param p Pointer to the particle data
 * @param sidm_properties The SIDM parameters
 */
__attribute__((always_inline)) INLINE static float dark_matter_compute_timestep(
    const struct dmpart *restrict dmp, const struct sidm_props* sidm_props) {
    
    /* Constant to limit probability from being too large? */
    const float kappa = 1e-2;
    
    /* Timestep limiter (internal units) */
    const float dm_timestep = kappa / dmp->sidm_probability;
    
    return dm_timestep;
}

/**
 * @brief Kick the additional variables
 *
 * @param dmp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void dark_matter_kick_extra(struct dmpart* dmp, float dt) {
}


#endif

