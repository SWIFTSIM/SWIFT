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
#include "random.h"
#include "dark_matter.h"


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
    
    message("entering iact DM density");
    float wi, wj, wi_dx, wj_dx;
    
    const float r = sqrtf(r2);
    
    /* Get the masses. */
    const float mi = pi->mass;
    const float mj = pj->mass;
    
    /* Compute density of pi. */
    const float hi_inv = 1.f / hi;
    const float ui = r * hi_inv;
    
    kernel_deval(ui, &wi, &wi_dx);
    
    pi->rho += mj * wi;
    pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

    pi->density.wcount += wi;
    pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

    /* Compute density of pj. */
    const float hj_inv = 1.f / hj;
    const float uj = r * hj_inv;
    kernel_deval(uj, &wj, &wj_dx);
    
    pj->rho += mi * wj;
    pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
    pj->density.wcount += wj;
    pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);
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
    
    message("entering iact nonsym DM density");

    float wi, wi_dx;
    
    /* Get the masses. */
    const float mj = pj->mass;
    
    /* Get r and r inverse. */
    const float r = sqrtf(r2);
    
    const float h_inv = 1.f / hi;
    const float ui = r * h_inv;
    kernel_deval(ui, &wi, &wi_dx);
    
    pi->rho += mj * wi;
    pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
    
    pi->density.wcount += wi;
    pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

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
                                                               struct dmpart *restrict pi, const integertime_t ti_current) {
    
    /* Center of Mass Velocity of interacting particles */
    const double VCM[3] = {(pi->v_full[0] + pj->v_full[0])/2.0, (pi->v_full[1] + pj->v_full[1])/2.0, (pi->v_full[2] + pj->v_full[2])/2.0};
    
    double dw[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
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
    
    pj->sidm_data.si_v_full[0] = VCM[0] + dv * e[0];
    pj->sidm_data.si_v_full[1] = VCM[1] + dv * e[1];
    pj->sidm_data.si_v_full[2] = VCM[2] + dv * e[2];
    
    pi->sidm_data.si_v_full[0] = VCM[0] - dv * e[0];
    pi->sidm_data.si_v_full[1] = VCM[1] - dv * e[1];
    pi->sidm_data.si_v_full[2] = VCM[2] - dv * e[2];
    
    /*! change flag to indicate the particle has been scattered */
    pj->sidm_data.sidm_flag = 1;
    pi->sidm_data.sidm_flag = 1;
    
    /* Add counter of DM-DM collisions */
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
    const integertime_t ti_current, const struct sidm_props* sidm_props, const struct unit_system* us) {
    
    /* Velocities of interacting particles */
    const double dv[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
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

    /*printf("Rateij %f & Probij %f & sigma %f",Rate_SIDM_ij,Probability_SIDM_ij, sigma);*/
    
    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability_SIDM_i > randi) sidm_do_kick(pi, pj, ti_current);
    if (Probability_SIDM_j > randj) sidm_do_kick(pj, pi, ti_current);

}

/**
 * @brief Interaction between two dark matter particles during force loop (non-symmetric).
 * It calculates the probability of DM particles i & j of scattering within the next time step
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, const double dti, const double dtj,
    const integertime_t ti_current, const struct sidm_props* sidm_props, const struct unit_system* us) {
    
    /* Velocities of interacting particles */
    const double dv[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
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
        sidm_do_kick(pi, pj, ti_current);
    }
}


#endif

