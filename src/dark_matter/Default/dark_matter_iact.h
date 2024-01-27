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
#ifndef SWIFT_DEFAULT_DARK_MATTER_IACT_H
#define SWIFT_DEFAULT_DARK_MATTER_IACT_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "random.h"
#include "dark_matter.h"
#include "kernel_dark_matter.h"
#include "dark_matter_logger.h"
#include "timestep_sync_part.h"
#include "units.h"

/**
 * @brief Integrate the kernels using the trapezoidal rule.
 */
INLINE static double integrate_kernels(float r2, float hi, float hj) {
    
    float h_max = hi;
    if (hj > h_max) h_max = hj;
    h_max *= dm_kernel_gamma;
    
    /* Bin spacing. Assumes uniform spacing. */
    const float r = sqrtf(r2);
    const int N_bins = 20;
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
    const int N_bins = 20;
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
 * @brief Analytical form of Int[4 pi x^2 W1(x/hi) W1(x/hj+r/hj)]dx
 */
INLINE static double w1w1(float x, float hi, float hj, float r) {
    r /= hj;
    float x3 = x * x * x;
    float x2 = x * x;
    float hi3 = hi * hi * hi;
    float hj3 = hj * hj * hj;
    float f;

    f = 14. * hj3 * (1. - 6.*r*r + 6.*r*r*r) * (5.*hi3 - 18.*hi*x2 + 15.*x3);
    f += 45.*hj*hj*r*(-2. + 3.*r)*x*(7.*hi3 - 28.*hi*x2 + 24.*x3);
    f += 9.*hj*(-1. + 3.*r)*x2*(28.*hi3 - 120.*hi*x2 + 105.*x3);
    f += 105.*(2.*hi3*x3 - 9.*hi*x3*x2 + 8.*x3*x3);
    f *= 0.00119048*x3;
    f /= (hj3*hi3);
    f *= 4. * M_PI;
    return f;
}

/**
 * @brief Analytical form of Int[4 pi x^2 W1(x/hi) W2(x/hj+r/hj)]dx
 */
INLINE static double w1w2(float x, float hi, float hj, float r) {
    r /= hj;
    float x3 = x * x * x;
    float x2 = x * x;
    float hi3 = hi * hi * hi;
    float hj3 = hj * hj * hj;
    float f;

    f = 0.333333*hj3*hi3*(-1. + r)*(-1. + r)*(-1. + r)*x3 + 0.75*hj*hj*hi3*(-1. + r)*(-1. + r)*x2*x2;
    f -= 0.6*hj*hi*(-1.*hi*hi + 2.*hj*hj*(-1. + r)*(-1. + r))*(-1. + r)*x3*x2;
    f += 0.166667*(hi3 - 18.*hj*hj*hi*(-1. + r)*(-1. + r) + 6.*hj3*(-1. + r)*(-1. + r)*(-1. + r))*x3*x3;
    f += 2.57143*hj*(-1.*hi + hj*(-1. + r))*(-1. + r)*x3*x3*x;
    f -= 0.75*(hi - 3.*hj*(-1. + r))*x3*x3*x2;
    f += 0.666667*x3*x3*x3;
    f *= (-0.5);
    f /= (hj3*hi3);
    f *= 4. * M_PI;
    return f;
}

/**
 * @brief Analytical form of Int[4 pi x^2 W2(x/hi) W1(x/hj+r/hj)]dx
 */
INLINE static double w2w1(float x, float hi, float hj, float r) {
    r /= hj;
    float x3 = x * x * x;
    float x2 = x * x;
    float hi3 = hi * hi * hi;
    float hj3 = hj * hj * hj;
    float f;
    
    f = 5.*x3 * (84.*hi3 - 216.*hi*hi*x + 189.*hi*x2 - 56.*x3);
    f += 9.*hj*(-1. + 3.*r)*x2*(56.*hi3 - 140.*hi*hi*x + 120.*hi*x2 - 35.*x3);
    f += 18.*hj*hj*r*(-2. + 3.*r)*x*(35.*hi3 - 84.*hi*hi*x + 70.*hi*x2 - 20.*x3);
    f += 7.*hj3*(1. - 6.*r*r + 6.*r*r*r)*(20.*hi3 - 45.*hi*hi*x + 36.*hi*x2 - 10.*x3);
    f *= 0.00119048*x3;
    
    f /= (hj3*hi3);
    f *= 4. * M_PI;
    return f;
}

/**
 * @brief Analytical form of Int[4 pi x^2 W2(x/hi) W2(x/hj+r/hj)]dx
 */
INLINE static double w2w2(float x, float hi, float hj, float r) {
    r /= hj;
    float x3 = x * x * x;
    float x2 = x * x;
    float hi3 = hi * hi * hi;
    float hj3 = hj * hj * hj;
    float f;
    
    f = (hj3*hi3*(-1. + r)*(-1. + r)*(-1. + r)*x3)/3.;
    f += (3.*hj*hj*hi*hi*(-1. + r)*(-1. + r)*(hj + hi - hj*r)*x2*x2)/4.;
    f += (3.*hj*hi*(hi*hi - 3.*hj*hi*(-1. + r) + hj*hj*(-1. + r)*(-1. + r))*(-1. + r)*x3*x2)/5.;
    f += ((hi3 - 9.*hj*hi*hi*(-1. + r) + 9.*hj*hj*hi*(-1. + r)*(-1. + r) - hj3*(-1. + r)*(-1. + r)*(-1. + r))*x3*x3)/6.;
    f -= (3.*(hi*hi - 3.*hj*hi*(-1. + r) + hj*hj*(-1. + r)*(-1. + r))*x3*x2*x2)/7.;
    f += (3.*(hj + hi - hj*r)*x3*x3*x2)/8.f;
    f -= x3*x3*x3/9.;
    f /= (hj3*hi3);
    f *= -1.;
    f *= 4. * M_PI;
    return f;
}

    

/**
 * @brief Integrate the kernels using analytical solution
 */
INLINE static double integrate_kernels_analytical(float r2, float hi, float hj) {
    
    hi *= dm_kernel_gamma;
    hj *= dm_kernel_gamma;
    float h_max = hi;
    if (hj > h_max) h_max = hj;
    float hi3 = hi * hi * hi;
    float hj3 = hj * hj * hj;
    const float r = sqrtf(r2);
    
    double integrand = 0.;
    float x1 = hi/2.f;
    float x3 = hi;
    float x2, x4;
    
    x2 = hj/2.f-r;
    if (x2<0.) x2=0.f;
    x4 = hj-r;
    if (x4<0.) x4=0.f;
    
    if (x1 <= x2 && x2 <= x3 && x3 <= x4){
        integrand = w1w1(x1,hi,hj,r);
        integrand += w2w1(x2,hi,hj,r)-w2w1(x1,hi,hj,r);
        integrand += w2w2(x3,hi,hj,r)-w2w2(x2,hi,hj,r);
    }
    
    if (x1 <= x2 && x2 <= x4 && x4 <= x3){
        integrand = w1w1(x1,hi,hj,r);
        integrand += w2w1(x2,hi,hj,r)-w2w1(x1,hi,hj,r);
        integrand += w2w2(x4,hi,hj,r)-w2w2(x2,hi,hj,r);
    }
    
    if (x2 <= x1 && x1 <= x3 && x3 <= x4){
        integrand = w1w1(x2,hi,hj,r);
        integrand += w1w2(x1,hi,hj,r)-w1w2(x2,hi,hj,r);
        integrand += w2w2(x3,hi,hj,r)-w2w2(x1,hi,hj,r);
    }

    if (x2 <= x4 && x4 <= x1 && x1 <= x3){
        integrand = w1w1(x2,hi,hj,r);
        integrand += w1w2(x4,hi,hj,r)-w1w2(x2,hi,hj,r);
    }

    if (x1 <= x3 && x3 <= x2 && x2 <= x4){
        integrand = w1w1(x1,hi,hj,r);
        integrand += w2w1(x3,hi,hj,r)-w2w1(x1,hi,hj,r);
    }
    
    integrand *= dm_kernel_constant * dm_kernel_constant;
    integrand /= hi3;
    integrand /= hj3;

    return integrand;
}

/**
 * @brief Calculate the norm of double kernels integral.
 */
INLINE static double norm_for_kernels_analytical_integral(float hi, float hj) {
    
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
        
        integrand[i] = integrate_kernels_analytical(r_int * r_int, hi, hj);
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

    /* Neighbour's (drifted) velocity   */
    /* (we don't include a Hubble term since we are interested in the */
    /* velocity contribution at the location of the DM particle)      */
    const float dvi[3] = {pj->v_full[0] - pi->v_full[0],
                          pj->v_full[1] - pi->v_full[1],
                          pj->v_full[2] - pi->v_full[2]};

    pi->velocity_ngb[0] += dvi[0];
    pi->velocity_ngb[1] += dvi[1];
    pi->velocity_ngb[2] += dvi[2];

    pi->velocity_dispersion += (dvi[0] * dvi[0] + dvi[1] * dvi[1] + dvi[2] * dvi[2]);

    /* Compute density of pj. */
    const float hj_inv = 1.f / hj;
    const float uj = r * hj_inv;
    dm_kernel_deval(uj, &wj, &wj_dx);
    
    pj->rho += mi * wj;
    pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
    pj->density.wcount += wj;
    pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

    const float dvj[3] = {pi->v_full[0] - pj->v_full[0],
                          pi->v_full[1] - pj->v_full[1],
                          pi->v_full[2] - pj->v_full[2]};

    pj->velocity_ngb[0] += dvj[0];
    pj->velocity_ngb[1] += dvj[1];
    pj->velocity_ngb[2] += dvj[2];

    pj->velocity_dispersion += (dvj[0] * dvj[0] + dvj[1] * dvj[1] + dvj[2] * dvj[2]);

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

    /* Neighbour's (drifted) velocity   */
    /* (we don't include a Hubble term since we are interested in the */
    /* velocity contribution at the location of the DM particle)      */
    const float dvi[3] = {pj->v_full[0] - pi->v_full[0],
                          pj->v_full[1] - pi->v_full[1],
                          pj->v_full[2] - pi->v_full[2]};

    pi->velocity_ngb[0] += dvi[0];
    pi->velocity_ngb[1] += dvi[1];
    pi->velocity_ngb[2] += dvi[2];

    pi->velocity_dispersion += (dvi[0] * dvi[0] + dvi[1] * dvi[1] + dvi[2] * dvi[2]);

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
__attribute__((always_inline)) INLINE static void sidm_do_kick(struct dmpart *restrict pi,
                                                               struct dmpart *restrict pj, const integertime_t ti_current,
                                                               /*struct sidm_history* sidm_history,*/
                                                               const struct cosmology* cosmo,
                                                               const struct sidm_props* sidm_props,
                                                               const struct unit_system* us) {

    /* Center of Mass Velocity of interacting particles */
    double VCM[3];
    VCM[0] = cosmo->a_inv * (pi->sidm_data.v_full[0] + pj->sidm_data.v_full[0])/2.0;
    VCM[1] = cosmo->a_inv * (pi->sidm_data.v_full[1] + pj->sidm_data.v_full[1])/2.0;
    VCM[2] = cosmo->a_inv * (pi->sidm_data.v_full[2] + pj->sidm_data.v_full[2])/2.0;

    double dw[3];
    dw[0] = pi->sidm_data.v_full[0] - pj->sidm_data.v_full[0];
    dw[1] = pi->sidm_data.v_full[1] - pj->sidm_data.v_full[1];
    dw[2] = pi->sidm_data.v_full[2] - pj->sidm_data.v_full[2];

    double dv2 = dw[0] * dw[0] + dw[1] * dw[1] + dw[2] * dw[2];
    double dv = cosmo->a_inv * sqrt(dv2) / 2.0;
    float e[3] = {0.f,0.f,0.f};
    
    /* Direction of kick is randomly chosen or not, it depends on scattering model */
    if (sidm_props->with_isotropic_scattering) {
        
        /* Draw a random number between (0,1] */
        const float u = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM_theta);
        
        /* Calculate theta from prob. distribution */
        const float theta = acos(1.f - 2.f*u);
        
        /* Random number for other angle */
        const float rand_phi = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM_phi);
        
        /* Transform to random number in [0, 2 pi] range */
        const float phi = 2.f * M_PI * rand_phi;
        
        /* Randomly oriented unit vector */
        e[0] = sin(theta) * cos(phi);
        e[1] = sin(theta) * sin(phi);
        e[2] = cos(theta);

    } else if (sidm_props->with_anisotropic_scattering) {
        
        /* Draw a random number between (0,1] */
        const float u = random_unit_interval(pi->id_or_neg_offset, ti_current, random_number_SIDM_theta);
        
        double v = cosmo->a_inv * sqrt(dv2);
        double mx = sidm_props->mx / 10.;
        double mphi = sidm_props->mphi / 10.;
        double w = 300. * 1e5 * mphi / mx; /* physical units cm/s */
        w /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
        w *= units_cgs_conversion_factor(us, UNIT_CONV_TIME); /* physical but internal units now */
        double w2 = w * w;
        const float a = v * v / w2;
        const float u_a = u - (1.f + a ) / a;
        float cos_theta;
        
        if (a <= 0.f || u_a == 0.f){
            /* If we happen to be in this bad regime, let's go with isotropic scattering */
            cos_theta = 1.f - 2.f * u;

        } else {
            /* Otherwise let's calculate theta from prob. distribution */
            const float dx = u * a / (1.f + a) - 1.f;
            cos_theta = ( 1.f / dx + 1.f) * (2.f / a) + 1.f;

        }

        const float sin_theta = sqrt(1.f - cos_theta * cos_theta);

        /* Calculate theta from prob. distribution */
        /* const float theta = acos(x); */
        
        /* Random number for other angle */
        const float rand_phi = random_unit_interval(pj->id_or_neg_offset, ti_current, random_number_SIDM_phi);
        
        /* Transform to random number in [0, 2 pi] range */
        const float phi = 2.f * M_PI * rand_phi;
        
        /* Not so randomly oriented unit vector */
        /* e[0] = sin(theta) * cos(phi);
        e[1] = sin(theta) * sin(phi);
        e[2] = cos(theta);*/
        e[0] = sin_theta * cos(phi);
        e[1] = sin_theta * sin(phi);
        e[2] = cos_theta;
    }
    
    /*double energy_before, energy_after;
    double energy_prev_i = pi->sidm_data.v_full[0] * pi->sidm_data.v_full[0] + pi->sidm_data.v_full[1] * pi->sidm_data.v_full[1] + pi->sidm_data.v_full[2] * pi->sidm_data.v_full[2];

    double energy_prev_j = pj->sidm_data.v_full[0] * pj->sidm_data.v_full[0] + pj->sidm_data.v_full[1] * pj->sidm_data.v_full[1] + pj->sidm_data.v_full[2] * pj->sidm_data.v_full[2];
    */

    /* I'm doing the kicks here by updating the particle velocities.
     * Note that v_full = a^2 * dx/dt, with x the comoving coordinate.
     * At this point sidm->v_full is in physical units */
    pi->sidm_data.v_full[0] = VCM[0] - dv * e[0];
    pi->sidm_data.v_full[1] = VCM[1] - dv * e[1];
    pi->sidm_data.v_full[2] = VCM[2] - dv * e[2];

    pj->sidm_data.v_full[0] = VCM[0] + dv * e[0];
    pj->sidm_data.v_full[1] = VCM[1] + dv * e[1];
    pj->sidm_data.v_full[2] = VCM[2] + dv * e[2];
    
    /* Therefore, the code velocity kick needs to go back to comoving.. */
    pi->sidm_data.v_full[0] = cosmo->a * pi->sidm_data.v_full[0];
    pi->sidm_data.v_full[1] = cosmo->a * pi->sidm_data.v_full[1];
    pi->sidm_data.v_full[2] = cosmo->a * pi->sidm_data.v_full[2];
    
    pj->sidm_data.v_full[0] = cosmo->a * pj->sidm_data.v_full[0];
    pj->sidm_data.v_full[1] = cosmo->a * pj->sidm_data.v_full[1];
    pj->sidm_data.v_full[2] = cosmo->a * pj->sidm_data.v_full[2];

    /* Communicating this kick to logger */
    /*if (pi->sidm_data.sidm_flag > 0) {
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
    }*/

    
    /*! change flag to indicate the particle has been scattered */
    pj->sidm_data.sidm_flag = 1;
    pi->sidm_data.sidm_flag = 1;
    
    /* Add counter of DM-DM collisions of individual particles */
    pj->sidm_data.number_of_sidm_events += 1.f;
    pi->sidm_data.number_of_sidm_events += 1.f;
    
    pi->sidm_data.sidm_events_per_timestep += 1.f;
    pj->sidm_data.sidm_events_per_timestep += 1.f;
}

/**
 * @brief Calculate velocity-dependent cross section.
 * following the Yukawa Potential.
 */
INLINE static double momentum_transfer_sigma_model(float v, const struct sidm_props* sidm_props, const struct unit_system* us) {
    
    float alpha_x = sidm_props->alphax / 0.01;
    float alpha_x2 = alpha_x * alpha_x;
    float mx = sidm_props->mx / 10.;
    float mphi = sidm_props->mphi / 10.;
    float mphi_inv2 = (1./ mphi) * (1./ mphi);
    float mphi_inv4 = mphi_inv2 * mphi_inv2;
    float sigma0 = 274.868 * alpha_x2 * mx * mphi_inv4; /* physical units cm2/g */
    float w = 300. * 1e5 * mphi / mx; /* physical units cm/s */
    w /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    w *= units_cgs_conversion_factor(us, UNIT_CONV_TIME); /* physical but internal units now */
    float w2 = w * w;
    float w4 = w2 * w2;
    float v4 = v * v * v * v;
    double dln_sigma = 2. * log(1. + v * v * 0.5 /w2) - log(1. + v * v /w2);
    
    /* Scattering cross section in physical units, with extra factor of 2 since
     * it has proven to work better (Robertson +2018) */
    double sigma = sigma0 * 4. * w4 / v4;
    sigma *= dln_sigma;
    
    /* Scattering cross section in internal units */
    sigma *= units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    
    return sigma;
}

/**
 * @brief Calculate velocity-dependent cross section.
 * following the Yukawa Potential.
 */
INLINE static double velocity_dependent_sigma_model(float v, const struct sidm_props* sidm_props, const struct unit_system* us) {
    
    float alpha_x = sidm_props->alphax / 0.01;
    float alpha_x2 = alpha_x * alpha_x;
    float mx = sidm_props->mx / 10.f;
    float mphi = sidm_props->mphi / 10.f;
    float mphi_inv2 = (1.0f/ mphi) * (1.0f/ mphi);
    float mphi_inv4 = mphi_inv2 * mphi_inv2;
    float sigma0 = 274.868 * alpha_x2 * mx * mphi_inv4; /* physical units cm2/g */
    float w = 300.0f * 1e5 * mphi / mx; /* physical units cm/s */
    w /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    w *= units_cgs_conversion_factor(us, UNIT_CONV_TIME); /* physical but internal units now */
    float w2 = w * w;
    
    /* Scattering cross section in physical units */
    double sigma = sigma0 / (1.0 + v * v / w2);
    
    /* Scattering cross section in internal units */
    sigma *= units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    
    return sigma;
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, float time_base, integertime_t t_current,
    const struct cosmology* cosmo, const int with_cosmology, const struct sidm_props* sidm_props,
    const struct unit_system* us) {
/*    struct sidm_history* sidm_history, const struct cosmology* cosmo) {*/

    /* Velocities of interacting particles */

    double dv[3];
    dv[0] = pi->v_full[0] - pj->v_full[0];
    dv[1] = pi->v_full[1] - pj->v_full[1];
    dv[2] = pi->v_full[2] - pj->v_full[2];
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    const double vij = sqrt(v2) * cosmo->a_inv;

    /* Manage time interval of particle i */
    double dti;
    if (with_cosmology) {
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin =
                get_integer_time_begin(t_current - 1, pi->time_bin);
        dti = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
        dti = get_timestep(pi->time_bin, time_base);
    }

    double sigma = 0.;

    /* Scattering cross section per unit mass (in internal units) */

    if (sidm_props->with_constant_sigma) {

      sigma = sidm_props->sigma;

    } else if (sidm_props->with_velocity_dependent_sigma) {

      sigma = velocity_dependent_sigma_model(vij, sidm_props, us);

      if (pi->num_neighbours > 0) {
        pi->sidm_data.sigma += sigma / pi->num_neighbours;
      }
    } else if (sidm_props->with_momentum_transfer_sigma) {

      sigma = momentum_transfer_sigma_model(vij, sidm_props, us);

      if (pi->num_neighbours > 0){
        pi->sidm_data.sigma += sigma / pi->num_neighbours;
      }
    }

    /* DM particle mass */
    const double mass_j = pj->mass;

    /*float gij = integrate_kernels(r2, hi, hj);
    float normed_gij = norm_for_kernels_integral(hi, hj);*/

    float gij = integrate_kernels_analytical(r2, hi, hj) * cosmo->a3_inv;
    float normed_gij = norm_for_kernels_analytical_integral(hi, hj);

    float Rate_SIDM_i = mass_j * sigma * vij * gij / normed_gij;
    
    pi->sidm_probability += mass_j * sigma * vij * gij / normed_gij;

    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;

    /* Draw a random number */
    const float rand = random_unit_interval(pi->id_or_neg_offset, t_current, random_number_SIDM);

    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability_SIDM_i > rand) {

        /* Part j is not within the timestep, let's wake it up for the SIDM kick */
        timestep_sync_dmpart(pj);

        /* Doing SIDM kick */
        sidm_do_kick(pi, pj, t_current, cosmo, sidm_props, us);

        /* Log the kick */
        /*    dark_matter_log_num_events(sidm_history, 1); */
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
__attribute__((always_inline)) INLINE static void runner_iact_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, float time_base, integertime_t t_current,
    const struct cosmology* cosmo, const int with_cosmology, const struct sidm_props* sidm_props,
    const struct unit_system* us) {
/*    struct sidm_history* sidm_history, const struct cosmology* cosmo) {*/

    /* Velocities of interacting particles */
    double dv[3];
    dv[0] = pi->v_full[0] - pj->v_full[0];
    dv[1] = pi->v_full[1] - pj->v_full[1];
    dv[2] = pi->v_full[2] - pj->v_full[2];
    const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    const double vij = sqrt(v2) * cosmo->a_inv;

    /* Manage time interval of particles i & j */
    double dti;
    if (with_cosmology) {
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin =
                get_integer_time_begin(t_current - 1, pi->time_bin);
        dti = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
        dti = get_timestep(pi->time_bin, time_base);
    }

    /* Same for particle j */
    double dtj;
    if (with_cosmology) {
        const integertime_t tj_step = get_integer_timestep(pj->time_bin);
        const integertime_t tj_begin =
                get_integer_time_begin(t_current - 1, pj->time_bin);
        dtj = cosmology_get_delta_time(cosmo, tj_begin, tj_begin + tj_step);
    } else {
        dtj = get_timestep(pj->time_bin, time_base);
    }

    double sigma = 0.;

    /* Scattering cross section per unit mass (in internal units) */
    if (sidm_props->with_constant_sigma) {

      sigma = sidm_props->sigma;

    } else if (sidm_props->with_velocity_dependent_sigma) {

      sigma = velocity_dependent_sigma_model(vij, sidm_props, us);

      if (pi->num_neighbours > 0){
        pi->sidm_data.sigma += sigma / pi->num_neighbours;
      }

      if (pj->num_neighbours > 0) {
        pj->sidm_data.sigma += sigma / pj->num_neighbours;
      }

    } else if (sidm_props->with_momentum_transfer_sigma) {

      sigma = momentum_transfer_sigma_model(vij, sidm_props, us);

      if (pi->num_neighbours > 0) {
        pi->sidm_data.sigma += sigma / pi->num_neighbours;
      }

      if (pj->num_neighbours > 0){
        pj->sidm_data.sigma += sigma / pj->num_neighbours;
      }

    }

    /* DM particle mass */
    const double mass_i = pi->mass;
    const double mass_j = pj->mass;

    /* Calculate scattering rate */
    /*float gij = integrate_kernels(r2, hi, hj);
    float normed_gij = norm_for_kernels_integral(hi, hj);
    float gji = integrate_kernels(r2, hj, hi);
    float normed_gji = norm_for_kernels_integral(hj, hi);*/
    float gij = integrate_kernels_analytical(r2, hi, hj) * cosmo->a3_inv;
    float normed_gij = norm_for_kernels_analytical_integral(hi, hj);
    float gji = integrate_kernels_analytical(r2, hj, hi) * cosmo->a3_inv;
    float normed_gji = norm_for_kernels_analytical_integral(hj, hi);

    float Rate_SIDM_i = mass_j * sigma * vij * gij / normed_gij;
    float Rate_SIDM_j = mass_i * sigma * vij * gji / normed_gji;
    
    pi->sidm_probability += mass_j * sigma * vij * gij / normed_gij;
    pj->sidm_probability += mass_i * sigma * vij * gji / normed_gji;

    /* Calculate SIDM probability */
    float Probability_SIDM_i = Rate_SIDM_i * dti;
    float Probability_SIDM_j = Rate_SIDM_j * dtj;
    float Probability = 0.5 * (Probability_SIDM_i + Probability_SIDM_j);

    /* Draw a random number */
    const float rand = random_unit_interval(pi->id_or_neg_offset, t_current, random_number_SIDM);

    /* Are we lucky? If so we have DM-DM interactions */
    if (Probability > rand) {

        /* Doing SIDM kick */
        sidm_do_kick(pi, pj, t_current, cosmo, sidm_props, us);

        /* Log the kick */
        /*dark_matter_log_num_events(sidm_history, 1);*/
    }
}


#endif

