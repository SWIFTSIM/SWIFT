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
#ifndef SWIFT_DARK_MATTER_IACT_H
#define SWIFT_DARK_MATTER_IACT_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "random.h"
#include "dark_matter.h"
#include "kernel_dark_matter.h"
#include "dark_matter_logger.h"
#include "timestep_sync_part.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_dark_matter_density(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H) {
    
    float wi, wj, wi_dx, wj_dx;
    
    const float r = sqrtf(r2);
    
    /* Get the masses. */
    const float mi = pi->mass;
    const float mj = pj->mass;
    
    /* Compute density of pi. */
    const float hi_inv = 1.f / hi;
    const float ui = r * hi_inv;
    
    dm_kernel_deval(ui, &wi, &wi_dx);
    
    pi->rho += mj * wi;
    pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

    pi->density.wcount += wi;
    pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

    /* Compute density of pj. */
    const float hj_inv = 1.f / hj;
    const float uj = r * hj_inv;
    dm_kernel_deval(uj, &wj, &wj_dx);
    
    pj->rho += mi * wj;
    pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
    pj->density.wcount += wj;
    pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);
    
    /* Velocities of interacting particles */
    const double dv[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    double vij = sqrt(v2);
    
    pj->avg_pair_v += vij;
    pi->avg_pair_v += vij;
    
    /* Increasing counters */
    ++pi->num_neighbours;
    ++pj->num_neighbours;

}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_density(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    const struct dmpart* pj, float a, float H) {
    
    float wi, wi_dx;
    
    /* Get the masses. */
    const float mj = pj->mass;
    
    /* Get r and r inverse. */
    const float r = sqrtf(r2);
    
    const float h_inv = 1.f / hi;
    const float ui = r * h_inv;
    dm_kernel_deval(ui, &wi, &wi_dx);
    
    pi->rho += mj * wi;
    pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
    
    pi->density.wcount += wi;
    pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
    
    /* Neighbour's (drifted) velocity in the frame of the dark matter particle
     * (we don't include a Hubble term since we are interested in the
     * velocity contribution at the location of the particle) */
    const double dv[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    double vij = sqrt(v2);    
    pi->avg_pair_v += vij;
    
    /* Contribution to the smoothed velocity */
    pi->velocity_ngb[0] += mj * dv[0] * wi;
    pi->velocity_ngb[1] += mj * dv[1] * wi;
    pi->velocity_ngb[2] += mj * dv[2] * wi;

    /* Contribution to the smoothed squared relative velocity (for dispersion)
     * We will convert this to actual dispersion later. */
    pi->velocity_dispersion += mj * wi * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

    
    /* Increasing counter */
    ++pi->num_neighbours;
}

/**
 * @brief Perform the 'kick' operation on both #gparts
 *
 * @param gpj #gpart
 * @param gpi #gpart
 * @param ti_current Current integer time (for random numbers).
 *
 */
__attribute__((always_inline)) INLINE static void sidm_do_kick(struct dmpart *restrict pj,
                                                               struct dmpart *restrict pi, const integertime_t ti_current,
                                                               struct sidm_history* sidm_history) {
    
    /* Center of Mass Velocity of interacting particles */
    const double VCM[3] = {(pi->sidm_data.v_full[0] + pj->sidm_data.v_full[0])/2.0, (pi->sidm_data.v_full[1] + pj->sidm_data.v_full[1])/2.0, (pi->sidm_data.v_full[2] + pj->sidm_data.v_full[2])/2.0};
    double dw[3] = {pi->sidm_data.v_full[0] - pj->sidm_data.v_full[0], pi->sidm_data.v_full[1] - pj->sidm_data.v_full[1], pi->sidm_data.v_full[2] - pj->sidm_data.v_full[2]};
    double dv2 = dw[0] * dw[0] + dw[1] * dw[1] + dw[2] * dw[2];
    double dv = sqrt(dv2) / 2.0;
    
    /* Direction of kick is randomly chosen */
    
    /* Draw a random number */
    const float rand_theta = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM_theta);
    
    /* Transform to random number in [0, pi] */
    const float theta = M_PI * rand_theta;
    
    /* Random number for other angle */
    const float rand_phi = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM_phi);
    
    /* Transform to random number in [-pi, pi] range */
    const float phi = 2.f * M_PI * rand_phi - M_PI;
    
    /* Randomly oriented unit vector */
    float e[3] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    
    double p_prev_i[3] = {pi->sidm_data.v_full[0], pi->sidm_data.v_full[1], pi->sidm_data.v_full[2]};
    double p_prev_j[3] = {pj->sidm_data.v_full[0], pj->sidm_data.v_full[1], pj->sidm_data.v_full[2]};
    
    double p_before[3], p_after[3];

    double energy_before, energy_after;
    double energy_prev_i = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];

    double energy_prev_j = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];

    pj->sidm_data.v_full[0] = VCM[0] + dv * e[0];
    pj->sidm_data.v_full[1] = VCM[1] + dv * e[1];
    pj->sidm_data.v_full[2] = VCM[2] + dv * e[2];
    
    pi->sidm_data.v_full[0] = VCM[0] - dv * e[0];
    pi->sidm_data.v_full[1] = VCM[1] - dv * e[1];
    pi->sidm_data.v_full[2] = VCM[2] - dv * e[2];

    /* Communicating this kick to logger */
    if (pi->sidm_data.sidm_flag > 0) {
        energy_before = 0.;
        p_before[0] = 0.;
        p_before[1] = 0.;
        p_before[2] = 0.;
        
        energy_after = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];
        
        p_after[0] = pi->sidm_data.v_full[0];
        p_after[1] = pi->sidm_data.v_full[1];
        p_after[2] = pi->sidm_data.v_full[2];
        
        energy_after -= energy_prev_i;
        p_after[0] -= p_prev_i[0];
        p_after[1] -= p_prev_i[1];
        p_after[2] -= p_prev_i[2];
        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after, p_before, p_after);

    } else {
        energy_before = energy_prev_i;
        p_before[0] = p_prev_i[0];
        p_before[1] = p_prev_i[1];
        p_before[2] = p_prev_i[2];
        
        energy_after = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];
        
        p_after[0] = pi->sidm_data.v_full[0];
        p_after[1] = pi->sidm_data.v_full[1];
        p_after[2] = pi->sidm_data.v_full[2];

        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after, p_before, p_after);
    }

    if (pj->sidm_data.sidm_flag > 0) {
        energy_before = 0.;
        p_before[0] = 0.;
        p_before[1] = 0.;
        p_before[2] = 0.;

        energy_after = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];
        
        p_after[0] = pj->sidm_data.v_full[0];
        p_after[1] = pj->sidm_data.v_full[1];
        p_after[2] = pj->sidm_data.v_full[2];
        
        energy_after -= energy_prev_j;
        p_after[0] -= p_prev_j[0];
        p_after[1] -= p_prev_j[1];
        p_after[2] -= p_prev_j[2];
        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after, p_before, p_after);
   
    } else {
        
        energy_before = energy_prev_j;
        p_before[0] = p_prev_j[0];
        p_before[1] = p_prev_j[1];
        p_before[2] = p_prev_j[2];
        
        energy_after = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];
        
        p_after[0] = pj->sidm_data.v_full[0];
        p_after[1] = pj->sidm_data.v_full[1];
        p_after[2] = pj->sidm_data.v_full[2];

        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after, p_before, p_after);
    }

    
    /*! change flag to indicate the particle has been scattered */
    pj->sidm_data.sidm_flag = 1;
    pi->sidm_data.sidm_flag = 1;
    
    /* Add counter of DM-DM collisions of individual particles */
    pj->sidm_data.num_sidm += 1.f;
    pi->sidm_data.num_sidm += 1.f;
    
}

/**
 * @brief Interaction between two dark matter particles during force loop
 * It calculates the probability of DM particles i & j of scattering within the next time step
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param ti_current Current integer time (for random numbers).
 */
__attribute__((always_inline)) INLINE static void runner_iact_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, const double dti, const double dtj,
    const integertime_t ti_current, const struct sidm_props* sidm_props, const struct unit_system* us,
    struct sidm_history* sidm_history) {
    
    /* Velocities of interacting particles */
    const double dv[3] = {pi->sidm_data.v_full[0] - pj->sidm_data.v_full[0], pi->sidm_data.v_full[1] - pj->sidm_data.v_full[1], pi->sidm_data.v_full[2] - pj->sidm_data.v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    double vij = sqrt(v2);
    
    /* Scattering cross section per unit mass (in internal units) */
    const double sigma = sidm_props->sigma;
    
    /* DM particle mass */
    const double mass_i = pi->mass;
    const double mass_j = pj->mass;

    /* DM-DM distance */
    float hi_3 = hi * hi * hi;
    float hj_3 = hj * hj * hj;

    float a_inv = 1.0f / a;
    float a_inv4 = a_inv * a_inv * a_inv * a_inv;
    
    /* Calculate scattering rate */
    float Rate_SIDM_i = sigma * mass_i * vij * a_inv4 / (4.0f * M_PI * hi_3 / 3.0f);
    float Rate_SIDM_j = sigma * mass_j * vij * a_inv4 / (4.0f * M_PI * hj_3 / 3.0f);

    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;
    float Probability_SIDM_j = Rate_SIDM_j * dtj;

    /* Draw a random number */
    const float randi = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM);
    const float randj = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM);

    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability_SIDM_i > randi || Probability_SIDM_j > randj) {
        
        /* Doing SIDM kick */
        sidm_do_kick(pi, pj, ti_current, sidm_history);
        
        /* Log the kick */
        dark_matter_log_num_events(sidm_history, 1);

    }
}

/**
 * @brief Interaction between two dark matter particles during force loop (non-symmetric).
 * It calculates the probability of DM particles i & j of scattering within the next time step
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not active).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, const double dti, const double dtj,
    const integertime_t ti_current, const struct sidm_props* sidm_props, const struct unit_system* us,
    struct sidm_history* sidm_history) {
    
    /* Velocities of interacting particles */
    const double dv[3] = {pi->sidm_data.v_full[0] - pj->sidm_data.v_full[0], pi->sidm_data.v_full[1] - pj->sidm_data.v_full[1], pi->sidm_data.v_full[2] - pj->sidm_data.v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    double vij = sqrt(v2);
    
    /* Scattering cross section per unit mass (in internal units) */
    const double sigma = sidm_props->sigma;
    
    /* DM particle mass */
    const double mass_i = pi->mass;
    
    /* DM-DM distance */
    float hi_3 = hi * hi * hi;
    
    float a_inv = 1.0f / a;
    float a_inv4 = a_inv * a_inv * a_inv * a_inv;
    
    /* Calculate scattering rate */
    float Rate_SIDM_i = sigma * mass_i * vij * a_inv4 / (4.0f * M_PI * hi_3 / 3.0f);
    
    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;
    
    /* Draw a random number */
    const float rand = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM);
    
    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability_SIDM_i > rand) {
        
        /* Part j is not within the timestep. Let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pj);
        
        /* Doing SIDM kick */
        sidm_do_kick(pi, pj, ti_current, sidm_history);
        
        /* Log the kick */
        dark_matter_log_num_events(sidm_history, 1);
    }
}


#endif

