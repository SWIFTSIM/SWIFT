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
#ifndef SWIFT_DEFAULT_DARK_MATTER_H
#define SWIFT_DEFAULT_DARK_MATTER_H

/* Config parameters. */
#include "../config.h"
#include "approx_math.h"
#include "dark_matter_part.h"
#include "dark_matter_iact.h"

/**
 * @brief Prepares a dm-particle for its interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void dark_matter_init_dmpart(struct dmpart* dmp) {
    
    dmp->density.wcount = 0.f;
    dmp->density.wcount_dh = 0.f;
    dmp->rho = 0.f;
    dmp->density.rho_dh = 0.f;
    
    dmp->sidm_probability = 0.f;
    dmp->time_step_size = 0.f;
    dmp->num_neighbours = 0.f;

    dmp->velocity_ngb[0] = 0.f;
    dmp->velocity_ngb[1] = 0.f;
    dmp->velocity_ngb[2] = 0.f;
    dmp->velocity_dispersion = 0.f;

}


/**
 * @brief Prepares a dm-particle for its SIDM interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_velocities(struct dmpart* dmp, float dt_drift) {
    
    /* No SIDM flag at the beginning */
    dmp->sidm_data.sidm_flag = 0.0f;
    
    /* Reset counter */
    dmp->sidm_data.sidm_events_per_timestep = 0.0f;

    /* Drift time */
    dmp->sidm_data.dt_drift = dt_drift;

    /* Make a copy of particle velocity */
    dmp->sidm_data.v_full[0] = dmp->v_full[0];
    dmp->sidm_data.v_full[1] = dmp->v_full[1];
    dmp->sidm_data.v_full[2] = dmp->v_full[2];

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

    /* Particle velocity */
    dmp->sidm_data.v_full[0] = dmp->v_full[0];
    dmp->sidm_data.v_full[1] = dmp->v_full[1];
    dmp->sidm_data.v_full[2] = dmp->v_full[2];

    dmp->sidm_data.dt_drift = 0.0f;

    dmp->time_bin = 0;

    dark_matter_init_dmpart(dmp);

    if (sidm_props->with_constant_sigma) {

      dmp->sidm_data.sigma = sidm_props->sigma;

    } else {

      dmp->sidm_data.sigma = 0.;

    }

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
                struct dmpart* restrict dmp) {}


/**
 * @brief Finishes the calculation of density on dark matter particles
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_end_density(
    struct dmpart* dmp, const struct cosmology* cosmo, const struct sidm_props* sidm_props,
    const double dt) {

    /* Some smoothing length multiples. */
    const float h = dmp->h;
    const float h_inv = 1.0f / h;                        /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv);        /* 1/h^d */
    const float h_inv_dim_plus_one = h_inv_dim * h_inv;  /* 1/h^(d+1) */
    
    /* Final operation on the density (add self-contribution). */
    dmp->rho += dmp->mass * dm_kernel_root;
    /*const float rho_inv = 1.f / dmp->rho;*/

    dmp->density.rho_dh -= hydro_dimension * dmp->mass * dm_kernel_root;
    dmp->density.wcount += dm_kernel_root;
    dmp->density.wcount_dh -= hydro_dimension * dm_kernel_root;

    /* Finish the calculation by inserting the missing h-factors */
    dmp->rho *= h_inv_dim;
    dmp->density.rho_dh *= h_inv_dim_plus_one;
    dmp->density.wcount *= h_inv_dim;
    dmp->density.wcount_dh *= h_inv_dim_plus_one;
    /*/*dmp->velocity_dispersion *= h_inv_dim * rho_inv;
    dmp->velocity_ngb[0] *= h_inv_dim * rho_inv;
    dmp->velocity_ngb[1] *= h_inv_dim * rho_inv;
    dmp->velocity_ngb[2] *= h_inv_dim * rho_inv;*/

    dmp->velocity_dispersion /= dmp->num_neighbours;
    dmp->velocity_ngb[0] /= dmp->num_neighbours;
    dmp->velocity_ngb[1] /= dmp->num_neighbours;
    dmp->velocity_ngb[2] /= dmp->num_neighbours;

    dmp->time_step_size = dt;

    /* Finish velocity dispersion calculation. Currently, the variable
    *  velocity_ngb holds <v^2> instead. */
    const double v2 = dmp->velocity_ngb[0] * dmp->velocity_ngb[0] +
                      dmp->velocity_ngb[1] * dmp->velocity_ngb[1] +
                      dmp->velocity_ngb[2] * dmp->velocity_ngb[2];

    dmp->velocity_dispersion -= v2;

    /* Low number of neighbours may cause this, avoid .. */
    if (dmp->velocity_dispersion < 0.f) dmp->velocity_dispersion *= -1.f;

    dmp->velocity_dispersion /= 3.f;
    dmp->velocity_dispersion = sqrt(dmp->velocity_dispersion);

    /* Let's initialize here sigma if we have velocity dependence */
    if (!sidm_props->with_constant_sigma) {
      dmp->sidm_data.sigma = 0.;
    }
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_part_has_no_neighbours(
    struct dmpart* restrict dmp, const struct cosmology* cosmo) {
    
    /* Some smoothing length multiples. */
    const float h = dmp->h;
    const float h_inv = 1.0f / h;                 /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

    /* Re-set problematic values */
    dmp->rho = dmp->mass * dm_kernel_root * h_inv_dim;
    dmp->density.wcount = dm_kernel_root * h_inv_dim;
    dmp->density.rho_dh = 0.f;
    dmp->density.wcount_dh = 0.f;
    dmp->num_neighbours = 0.f;

    dmp->velocity_ngb[0] = FLT_MAX;
    dmp->velocity_ngb[1] = FLT_MAX;
    dmp->velocity_ngb[2] = FLT_MAX;
    dmp->velocity_dispersion = FLT_MAX;
}



/**
 * @brief Updates #dmparts velocities
 *
 * @param dmp #dmpart
 *
 */
__attribute__((always_inline)) INLINE static void sidm_kick_to_dmpart(struct dmpart *restrict dmp) {

    if (dmp->sidm_data.sidm_flag > 0) {

#ifdef SWIFT_DEBUG_CHECKS
        message("SIDM kick for dmpart %lld h=%e", dmp->id_or_neg_offset, dmp->h);
#endif
        
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
    /* Remove flag */
    dmp->sidm_data.sidm_flag = 0.f;
    
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
    const struct dmpart *restrict dmp, const struct sidm_props* sidm_props,
    const struct cosmology* cosmo) {

    double dm_timestep;

    /* Only limit the time step if you are surrounded by lots of
     * neighbours, otherwise do not */
    if (dmp->num_neighbours > 5) {

      /* Constant to limit probability from being too large? */
      const float kappa = 1e-2;
      const float rho_phys = dmp->rho * cosmo->a3_inv;
      const float veldisp_phys = dmp->velocity_dispersion * cosmo->a_inv;

      /* Scattering cross section is already in physical units */
      const double probability = rho_phys * dmp->sidm_data.sigma * veldisp_phys;

      /* Time-step limiter (internal units) */
      dm_timestep = kappa / probability;

      /* The new timestep of the DMpart cannot be smaller than the minimum allowed
         * time-step */
      dm_timestep = max(dm_timestep, sidm_props->time_step_min);

    } else {

      dm_timestep = sidm_props->time_step_max;

    }

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

