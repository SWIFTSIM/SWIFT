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
 * @brief Integrate the kernels using the trapezoidal rule.
 */
INLINE static double integrate_kernels(float r2, float hi, float hj) {
    
    float h_max = hi;
    if (hj > h_max) h_max = hj;
    h_max *= dm_kernel_gamma;
    
    /* Bin spacing. Assumes uniform spacing. */
    const float r = sqrtf(r2);
    const int N_bins = 50;
    const float bin_size = h_max / N_bins;
    
    /* Array for the integrand */
    double integrand[N_bins];
    const int i_min = 0;
    const int i_max = N_bins - 1;

    float wi, wj, ui, uj, r_int;
    const float hi_inv = 1.f / hi;
    const float hj_inv = 1.f / hj;

    const float hi_inv3 = hi_inv * hi_inv * hi_inv;
    const float hj_inv3 = hj_inv * hj_inv * hj_inv;

    /* Calculate integral function */
    for (int i = i_min; i < i_max + 1; i++) {
        
        r_int = (i + 1) * bin_size;
        
        ui = r_int * hi_inv;
        dm_kernel_eval(ui, &wi);
        
        uj = (r_int + r) * hj_inv;
        dm_kernel_eval(uj, &wj);
        
        integrand[i] = wi * wj * r_int * r_int;
    }

    /* Integrate using trapezoidal rule */
    double result = 0.;
    
    for (int i = i_min; i < i_max + 1; i++) {
        result += integrand[i];
    }
    
    /* Adding missing h factors */
    result *= hi_inv3;
    result *= hj_inv3;
    
    /* Update end bins since contribution was overcounted when summing up all
     * entries */
    result -= 0.5 * (integrand[i_min] + integrand[i_max]);
    result *= bin_size * 4.f * M_PI;
    
    /* Done */
    return result;
}

/**
 * @brief Calculate the norm of double kernels integral.
 */
INLINE static double norm_for_kernels_integral(float hi, float hj) {
    
    float h_max = hi;
    if (hj > h_max) h_max = hj;
    h_max *= dm_kernel_gamma;
    
    /* Bin spacing. Assumes uniform spacing. */
    const int N_bins = 50;
    const float bin_size = h_max / N_bins;
    float r_int;
    
    /* Array for the integrand */
    double integrand[N_bins];
    const int i_min = 0;
    const int i_max = N_bins - 1;
    
    /* Calculate integral function */
    for (int i = i_min; i < i_max + 1; i++) {
        
        r_int = (i + 1) * bin_size;
        
        integrand[i] = integrate_kernels(r_int * r_int, hi, hj);
        integrand[i] *= r_int * r_int;
    }
    
    /* Integrate using trapezoidal rule */
    double result = 0.;
    
    for (int i = i_min; i < i_max + 1; i++) {
        result += integrand[i];
    }
    
    /* Update end bins since contribution was overcounted when summing up all
     * entries */
    result -= 0.5 * (integrand[i_min] + integrand[i_max]);
    result *= bin_size * 4.f * M_PI;
    
    /* Done */
    return result;
    
}

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
    
    /* Delta velocities :
     * (Note we don't include a Hubble term since we are interested in the
     * velocity contribution at the location of the particle) */
    const double dv[3] = {pj->v_full[0] - pi->v_full[0], pj->v_full[1] - pi->v_full[1], pj->v_full[2] - pi->v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    
    pj->avg_pair_v += sqrt(v2);
    pi->avg_pair_v += sqrt(v2);
    
    /* Contribution to the smoothed velocity */
    pi->velocity_ngb[0] += mj * dv[0] * wi;
    pi->velocity_ngb[1] += mj * dv[1] * wi;
    pi->velocity_ngb[2] += mj * dv[2] * wi;

    const double dvi[3] = {pi->v_full[0] - pj->v_full[0], pi->v_full[1] - pj->v_full[1], pi->v_full[2] - pj->v_full[2]};
    pj->velocity_ngb[0] += mi * dvi[0] * wj;
    pj->velocity_ngb[1] += mi * dvi[1] * wj;
    pj->velocity_ngb[2] += mi * dvi[2] * wj;

    /* Contribution to the smoothed squared relative velocity (for dispersion)
     * We will convert this to actual dispersion later. */
    pi->velocity_dispersion += mj * wi * v2;
    pj->velocity_dispersion += mi * wj * v2;
    
    /* Increasing counters */
    ++pi->num_neighbours;
    ++pj->num_neighbours;
    
    /*float gij = integrate_kernels(r2, hi, hj);
    float normed_gij = norm_for_kernels_integral(hi, hj);
    float gji = integrate_kernels(r2, hj, hi);
    float normed_gji = norm_for_kernels_integral(hj, hi);

    pi->sidm_probability += mj * sqrt(v2) * gij / normed_gij;
    pj->sidm_probability += mi * sqrt(v2) * gji / normed_gji;*/

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
    
    /* Neighbour's velocity in the frame of the dark matter particle
     * (we don't include a Hubble term since we are interested in the
     * velocity contribution at the location of the particle) */
    const double dv[3] = {pj->v_full[0] - pi->v_full[0], pj->v_full[1] - pi->v_full[1], pj->v_full[2] - pi->v_full[2]};
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    pi->avg_pair_v += sqrt(v2);
    
    /* Contribution to the smoothed velocity */
    pi->velocity_ngb[0] += mj * dv[0] * wi;
    pi->velocity_ngb[1] += mj * dv[1] * wi;
    pi->velocity_ngb[2] += mj * dv[2] * wi;

    /* Contribution to the smoothed squared relative velocity (for dispersion)
     * We will convert this to actual dispersion later. */
    pi->velocity_dispersion += mj * wi * v2;

    /* Increasing counter */
    ++pi->num_neighbours;
    
    /* Calculate probability */
    /*float gij = integrate_kernels(r2, hi, hj);
    float normed_gij = norm_for_kernels_integral(hi, hj);
    
    pi->sidm_probability += mj * sqrt(v2) * gij / normed_gij;*/
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
        
        energy_after = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];
        
        energy_after -= energy_prev_i;
        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after);

    } else {
        energy_before = energy_prev_i;

        energy_after = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];
        
        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after);
    }

    if (pj->sidm_data.sidm_flag > 0) {
        energy_before = 0.;

        energy_after = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];
        
        energy_after -= energy_prev_j;

        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after);
   
    } else {
        
        energy_before = energy_prev_j;
        
        energy_after = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];
        
        dark_matter_log_total_kinetic_energy(sidm_history, energy_before, energy_after);
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
    
    float gij = integrate_kernels(r2, hi, hj);
    float gji = integrate_kernels(r2, hj, hi);
    
    float normed_gij = norm_for_kernels_integral(hi, hj);
    float normed_gji = norm_for_kernels_integral(hj, hi);

    float Rate_SIDM_i = mass_j * sigma * vij * gij / normed_gij;
    float Rate_SIDM_j = mass_i * sigma * vij * gji / normed_gji;
    
    pi->sidm_probability += mass_j * sigma * vij * gij * dti / normed_gij;
    pj->sidm_probability += mass_i * sigma * vij * gji * dtj / normed_gji;

    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;
    float Probability_SIDM_j = Rate_SIDM_j * dtj;
    float Probability = Probability_SIDM_i + Probability_SIDM_j;
    
    /* Draw a random number */
    const float randi = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM);
    const float randj = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM);

    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability > randi || Probability > randj) {
        
        /* If part j is not within the timestep, let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pj);
        
        /* If part i is not within the timestep, let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pi);
        
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
    const double mass_j = pj->mass;

    /* Calculate scattering rate */
    float gij = integrate_kernels(r2, hi, hj);
    float normed_gij = norm_for_kernels_integral(hi, hj);
    float gji = integrate_kernels(r2, hj, hi);
    float normed_gji = norm_for_kernels_integral(hj, hi);

    float Rate_SIDM_i = mass_j * sigma * vij * gij / normed_gij;
    float Rate_SIDM_j = mass_i * sigma * vij * gji / normed_gji;
    
    pi->sidm_probability += mass_j * sigma * vij * gij * dti / normed_gij;

    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;
    float Probability_SIDM_j = Rate_SIDM_j * dtj;
    float Probability = Probability_SIDM_i  + Probability_SIDM_j;

    /* Draw a random number */
    const float randi = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM);
    const float randj = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM);

    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability > randi || Probability > randj) {

        /* If part j is not within the timestep, let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pj);
        
        /* If part i is not within the timestep, let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pi);
        
        /* Doing SIDM kick */
        sidm_do_kick(pi, pj, ti_current, sidm_history);
        
        /* Log the kick */
        dark_matter_log_num_events(sidm_history, 1);
    }
}


#endif

