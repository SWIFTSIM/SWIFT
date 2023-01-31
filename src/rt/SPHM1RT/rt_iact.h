/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IACT_SPHM1RT_H
#define SWIFT_RT_IACT_SPHM1RT_H

#include "rt_gradients.h"

/**
 * @file src/rt/SPHM1RT/rt_iact.h
 * @brief Main header file for no radiative transfer scheme particle
 * interactions.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Preparation step for injection to gather necessary data.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle (not updated).
 * @param cosmo The cosmological model.
 * @param rt_props Properties of the RT scheme.
 */

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_injection_prep(const float r2, const float *dx,
                                     const float hi, const float hj,
                                     struct spart *si, const struct part *pj,
                                     const struct cosmology *cosmo,
                                     const struct rt_props *rt_props) {

  /* If the star doesn't have any neighbours, we
   * have nothing to do here. */
  if (si->density.wcount == 0.f) return;

  /* Compute the weight of the neighbouring particle */
  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);
  /* Get the gas density. */
  const float rhoj = hydro_get_comoving_density(pj);
  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* This is actually the inverse of the enrichment weight */
  /* we abuse the variable here */
  if (r2 != 0.f) {
#if defined(HYDRO_DIMENSION_3D)
    si->rt_data.injection_weight += mj / rhoj / r2;
#elif defined(HYDRO_DIMENSION_2D)
    si->rt_data.injection_weight += mj / rhoj / r;
#elif defined(HYDRO_DIMENSION_1D)
    si->rt_data.injection_weight += mj / rhoj;
#endif
  }
  /* get the radiation energy within injection radius */
  /* we need it only when we need to redistribute the radiation energy */
  if (rt_props->reinject) {
    float urad[RT_NGROUPS];
    rt_get_physical_urad_multifrequency(pj, cosmo, urad);
    for (int g = 0; g < RT_NGROUPS; g++) {
      si->rt_data.emission_reinject[g] += mj * urad[g];
    }
  }
}

/**
 * @brief Injection step interaction between star and hydro particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si Star particle.
 * @param pj Hydro particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param rt_props Properties of the RT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float *dx, const float hi, const float hj,
    struct spart *restrict si, struct part *restrict pj, float a, float H,
    const struct rt_props *rt_props) {

  /* If the star doesn't have any neighbours, we
   * have nothing to do here. */
  if (si->density.wcount == 0.f) return;
  if (si->rt_data.injection_weight == 0.f) return;

  /* the direction of the radiation injected */
  const float r = sqrtf(r2);
  const float minus_r_inv = -1.f / r;
  const float n_unit[3] = {dx[0] * minus_r_inv, dx[1] * minus_r_inv,
                           dx[2] * minus_r_inv};

  /* Get particle mass */
  const float mj = hydro_get_mass(pj);
  const float mj_inv = 1.f / mj;
  /* Get the gas density. */
  /* TK comment: we need to be careful in the cosmological case here: */
  const float rhoj = hydro_get_comoving_density(pj);
  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* collect the enrichment weights from the neighborhood */
  float tot_weight_inv;
  tot_weight_inv = 1.f / si->rt_data.injection_weight;

  float injection_weight = 0.f;
  /* the enrichment weight of individual gas particle */
  if (r2 != 0.f) {
#if defined(HYDRO_DIMENSION_3D)
    injection_weight = mj / rhoj / r2;
#elif defined(HYDRO_DIMENSION_2D)
    injection_weight = mj / rhoj / r;
#elif defined(HYDRO_DIMENSION_1D)
    injection_weight = mj / rhoj;
#endif
  }

  float new_urad, new_frad;
  float urad[RT_NGROUPS];
  float frad[RT_NGROUPS][3];
  const float cred = rt_get_comoving_cred(pj, a);
  for (int g = 0; g < RT_NGROUPS; g++) {
    /* Inject energy. */
    if (rt_props->reinject) {
      new_urad = (si->rt_data.emission_this_step[g] +
                  si->rt_data.emission_reinject[g]) *
                 injection_weight * tot_weight_inv * mj_inv;
    } else {
      new_urad = si->rt_data.emission_this_step[g] * injection_weight *
                     tot_weight_inv * mj_inv +
                 pj->rt_data.conserved[g].urad;
    }

    urad[g] = new_urad;

    /* Inject flux. */
    /* We assume the path from the star to the gas is optically thin */
    new_frad = new_urad * cred;
    frad[g][0] = new_frad * n_unit[0];
    frad[g][1] = new_frad * n_unit[1];
    frad[g][2] = new_frad * n_unit[2];
  }

  rt_set_comoving_urad_multifrequency(pj, urad);
  rt_set_comoving_frad_multifrequency(pj, frad);
}

/**
 * @brief do radiation gradient computation
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param mode: mode=1 symmetric, mode=0 non symmetric
 *
 *
 */
__attribute__((always_inline)) INLINE static void
radiation_gradient_loop_function(float r2, const float *dx, float hi, float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, float a, float H,
                                 int mode) {

  struct rt_part_data *rpi = &pi->rt_data;
  struct rt_part_data *rpj = &pj->rt_data;

  float wi, wj, wi_dx, wj_dx;
  /* Get r */
  const float r = sqrtf(r2);

  /* part j */
  /* Get the kernel for hj */
  const float hj_inv = 1.0f / hj;

  /* Compute the kernel function for pj */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* part i */
  /* Get the kernel for hi */
  const float hi_inv = 1.0f / hi;

  /* Compute the kernel function for pi */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get mass */
  const float mj = hydro_get_mass(pj);
  const float mi = hydro_get_mass(pi);
  const float rhoj = hydro_get_comoving_density(pj);
  const float rhoi = hydro_get_comoving_density(pi);
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wj_dr = hjd_inv * wj_dx;
  const float wi_dr = hid_inv * wi_dx;

  float uradmfi[RT_NGROUPS];
  float uradmfj[RT_NGROUPS];

  const float credi = rt_get_comoving_cred(pi, a);
  const float credj = rt_get_comoving_cred(pj, a);

  /* use urad * c instead */
  float uradci;
  float uradcj;

  float fradmfi[RT_NGROUPS][3];
  float fradmfj[RT_NGROUPS][3];

  float fradi[3];
  float fradj[3];

  /*******************************/
  /* Computer gradient of radiation field times c */
  /*******************************/
  float gradi[3], gradj[3];
  int diffmode = 2;
  float divfi, divfj;

  /* gas density should not be zero */
  if ((rhoi == 0.f) || (rhoi == 0.f)) return;

  rt_get_comoving_urad_multifrequency(pi, uradmfi);
  rt_get_comoving_urad_multifrequency(pj, uradmfj);

  rt_get_comoving_frad_multifrequency(pi, fradmfi);
  rt_get_comoving_frad_multifrequency(pj, fradmfj);

  for (int g = 0; g < RT_NGROUPS; g++) {

    uradci = uradmfi[g] * credi;
    uradcj = uradmfj[g] * credj;
    radiation_gradient_SPH(uradci, uradcj, mi, mj, rpi->force.f, rpj->force.f,
                           rhoi, rhoj, wi_dr, wj_dr, dx, r, diffmode, gradi,
                           gradj);

    rpi->diffusion[g].graduradc[0] += gradi[0];
    rpi->diffusion[g].graduradc[1] += gradi[1];
    rpi->diffusion[g].graduradc[2] += gradi[2];
    if (mode == 1) {
      rpj->diffusion[g].graduradc[0] += gradj[0];
      rpj->diffusion[g].graduradc[1] += gradj[1];
      rpj->diffusion[g].graduradc[2] += gradj[2];
    }

    /*******************************/
    /* Now we need to compute the div of f terms */
    /*******************************/
    divfi = 0.0f;
    divfj = 0.0f;
    fradi[0] = fradmfi[g][0];
    fradi[1] = fradmfi[g][1];
    fradi[2] = fradmfi[g][2];
    fradj[0] = fradmfj[g][0];
    fradj[1] = fradmfj[g][1];
    fradj[2] = fradmfj[g][2];
    radiation_divergence_SPH(fradi, fradj, mi, mj, rpi->force.f, rpj->force.f,
                             rhoi, rhoj, wi_dr, wj_dr, dx, r, diffmode, &divfi,
                             &divfj);
    rpi->viscosity[g].divf += divfi;
    if (mode == 1) {
      rpj->viscosity[g].divf += divfj;
    }
  }
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {
  radiation_gradient_loop_function(r2, dx, hi, hj, pi, pj, a, H, 1);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_gradient(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, float a, float H) {
  radiation_gradient_loop_function(r2, dx, hi, hj, pi, pj, a, H, 0);
}

/**
 * @brief do radiation force computation
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param mode: mode=1 symmetric, mode=0 non symmetric
 *
 */
__attribute__((always_inline)) INLINE static void radiation_force_loop_function(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, int mode) {

  struct rt_part_data *rpi = &pi->rt_data;
  struct rt_part_data *rpj = &pj->rt_data;

  float wi, wj, wi_dx, wj_dx;
  /* Get r */
  const float r = sqrtf(r2);

  /* part j */
  /* Get the kernel for hj */
  const float hj_inv = 1.0f / hj;

  /* Compute the kernel function for pj */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* part i */
  /* Get the kernel for hi */
  const float hi_inv = 1.0f / hi;

  /* Compute the kernel function for pi */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get mass */
  const float mj = hydro_get_mass(pj);
  const float mi = hydro_get_mass(pi);
  const float rhoj = hydro_get_comoving_density(pj);
  const float rhoi = hydro_get_comoving_density(pi);
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wj_dr = hjd_inv * wj_dx;
  const float wi_dr = hid_inv * wi_dx;

  float uradmfi[RT_NGROUPS];
  float uradmfj[RT_NGROUPS];

  const float credi = rt_get_comoving_cred(pi, a);
  const float credj = rt_get_comoving_cred(pj, a);

  float fradmfi[RT_NGROUPS][3];
  float fradmfj[RT_NGROUPS][3];

  rt_get_comoving_urad_multifrequency(pi, uradmfi);
  rt_get_comoving_urad_multifrequency(pj, uradmfj);

  rt_get_comoving_frad_multifrequency(pi, fradmfi);
  rt_get_comoving_frad_multifrequency(pj, fradmfj);

  /*******************************/
  /* CALCULATIONS OF TWO MOMENT EQUATIONS */
  /*******************************/

  /* Calculate the radiation energy term from the divergence of f */
  int diffmode;

  float diss_durad_term_i, diss_durad_term_j;
  float fradmagi, fradmagj;
  double fradmagidouble, fradmagjdouble;

  float hid_inv_temp, wi_dr_temp, hjd_inv_temp, wj_dr_temp;
  float drhou_low, diss_durad_term;

  float ddfi, ddfj;
  float diss_dfrad_term_i[3], diss_dfrad_term_j[3];
  float durad_dt_i, durad_dt_j;
  float dfrad_dt_i[3], dfrad_dt_j[3];
  float fraduniti[3], fradunitj[3];
  float vsig_diss_i, vsig_diss_j;

  float J_tensori[3][3];
  float J_tensorj[3][3];
  float ddi, ddj;
  float foxi, foxj;
  float sqi, sqj;
  float flimi, flimj;
  float F_tensori[3][3];
  float F_tensorj[3][3];
  float funiti[3], funitj[3];

  int diffmodeaniso;
  float diff_dfrad_term_i[3], diff_dfrad_term_j[3];
  const float r_inv = 1.f / r;

  float uradi, uradj;
  float uradi0, uradj0;
  float cred0 = fmaxf(credi, credj);
  float rhomean2;
  float fradi[3], fradj[3];
  double fradidouble[3], fradjdouble[3];
  float divfipar, divfjpar; /* divfipar is inside the loop, and divfi is summed
                               over particle */
  float divfi, divfj;
  float alpha_diss_i, alpha_f_diss_i, alpha_diss_j, alpha_f_diss_j;

  /* gas density should not be zero */
  if ((rhoi == 0.f) || (rhoi == 0.f)) return;

  for (int g = 0; g < RT_NGROUPS; g++) {
    uradi = uradmfi[g];
    uradj = uradmfj[g];
    uradi0 = uradi * credi / cred0;
    uradj0 = uradj * credj / cred0;
    fradi[0] = fradmfi[g][0];
    fradi[1] = fradmfi[g][1];
    fradi[2] = fradmfi[g][2];
    fradj[0] = fradmfj[g][0];
    fradj[1] = fradmfj[g][1];
    fradj[2] = fradmfj[g][2];

    alpha_diss_i = rpi->diffusion[g].alpha;
    alpha_f_diss_i = rpi->viscosity[g].alpha;
    alpha_diss_j = rpj->diffusion[g].alpha;
    alpha_f_diss_j = rpj->viscosity[g].alpha;

    divfi = rpi->viscosity[g].divf;
    divfj = rpj->viscosity[g].divf;

    /* do nothing if there is no radiation */
    if ((uradi == 0.f) && (uradj == 0.f)) continue;

#if defined(HYDRO_DIMENSION_1D)
    fradi[1] = 0.0f;
    fradi[2] = 0.0f;
    fradj[1] = 0.0f;
    fradj[2] = 0.0f;
#endif
#if defined(HYDRO_DIMENSION_2D)
    fradi[2] = 0.0f;
    fradj[2] = 0.0f;
#endif

    if ((fradi[0] == 0.f) && (fradi[1] == 0.f) && (fradi[2] == 0.f)) {
      fradmagi = 0.f;
    } else {
      fradidouble[0] = (double)(fradi[0]);
      fradidouble[1] = (double)(fradi[1]);
      fradidouble[2] = (double)(fradi[2]);
      fradmagidouble = sqrt(fradidouble[0] * fradidouble[0] +
                            fradidouble[1] * fradidouble[1] +
                            fradidouble[2] * fradidouble[2]);
      if (isinf(fradmagidouble) || isnan(fradmagidouble))
        error("Got inf/nan in fradmagi | %.6e| %.6e %.6e %.6e", fradmagidouble,
              fradi[0], fradi[1], fradi[2]);
      fradmagi = (float)(fradmagidouble);
    }

    if ((fradj[0] == 0.f) && (fradj[1] == 0.f) && (fradj[2] == 0.f)) {
      fradmagj = 0.f;
    } else {
      fradjdouble[0] = (double)(fradj[0]);
      fradjdouble[1] = (double)(fradj[1]);
      fradjdouble[2] = (double)(fradj[2]);
      fradmagjdouble = sqrt(fradjdouble[0] * fradjdouble[0] +
                            fradjdouble[1] * fradjdouble[1] +
                            fradjdouble[2] * fradjdouble[2]);
      if (isinf(fradmagjdouble) || isnan(fradmagjdouble))
        error("Got inf/nan in fradmagj | %.6e| %.6e %.6e %.6e", fradmagjdouble,
              fradj[0], fradj[1], fradj[2]);
      fradmagj = (float)(fradmagjdouble);
    }

    /*******************************/
    /* CALCULATIONS OF TWO MOMENT EQUATIONS */
    /*******************************/

    /* Calculate the radiation energy term from the divergence of f */
    diffmode = 2;
    divfipar = 0.0f;
    divfjpar = 0.0f;
    radiation_divergence_SPH(fradi, fradj, mi, mj, rpi->force.f, rpj->force.f,
                             rhoi, rhoj, wi_dr, wj_dr, dx, r, diffmode,
                             &divfipar, &divfjpar);

    /* Calculate the radiation flux term */
    /* funit is for Eddington tensor */
    if (fradmagi != 0.f) {
      funiti[0] = fradi[0] / fradmagi;
      funiti[1] = fradi[1] / fradmagi;
      funiti[2] = fradi[2] / fradmagi;
      funitj[0] = funiti[0];
      funitj[1] = funiti[1];
      funitj[2] = funiti[2];
    } else if (fradmagj != 0.f) {
      funitj[0] = fradj[0] / fradmagj;
      funitj[1] = fradj[1] / fradmagj;
      funitj[2] = fradj[2] / fradmagj;
      funiti[0] = funitj[0];
      funiti[1] = funitj[1];
      funiti[2] = funitj[2];
    } else {
      /* Nothing we can do */
      return;
    }

    /* Eddington factor (or optical thickness estimator?) */
    if (credi * uradi == 0.f) {
      foxi = expf(-rpi->params.chi[g] * rhoi * hi);
    } else {
      foxi = fmaxf(expf(-rpi->params.chi[g] * rhoi * hi),
                   fradmagi / (credi * uradi));
    }
    if (credj * uradj == 0.f) {
      foxj = expf(-rpj->params.chi[g] * rhoj * hj);
    } else {
      foxj = fmaxf(expf(-rpj->params.chi[g] * rhoj * hj),
                   fradmagj / (credj * uradj));
    }

    foxi = fminf(foxi, 1.0f);
    foxj = fminf(foxj, 1.0f);

    /* M1 closure with (modified) Eddington factor */
    /* protect against negative in the square root */
    sqi = 4.f - 3.f * foxi * foxi;
    sqj = 4.f - 3.f * foxj * foxj;
    flimi = fminf(1.f, (3.f + 4.f * foxi * foxi) / (5.f + 2.f * sqrtf(sqi)));
    flimj = fminf(1.f, (3.f + 4.f * foxj * foxj) / (5.f + 2.f * sqrtf(sqj)));

    /* compute the Eddington tensor (without radiation energy density yet) */

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        F_tensori[k][j] = 0.5f * (3.0f * flimi - 1.0f) * funiti[k] * funiti[j];
        F_tensorj[k][j] = 0.5f * (3.0f * flimj - 1.0f) * funitj[k] * funitj[j];
      }
    }

    F_tensori[0][0] += 0.5f * (1.0f - flimi);
    F_tensori[1][1] += 0.5f * (1.0f - flimi);
    F_tensori[2][2] += 0.5f * (1.0f - flimi);
    F_tensorj[0][0] += 0.5f * (1.0f - flimj);
    F_tensorj[1][1] += 0.5f * (1.0f - flimj);
    F_tensorj[2][2] += 0.5f * (1.0f - flimj);

#if defined(HYDRO_DIMENSION_1D)
    /* the eddington tensor is different in 1D */
    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        F_tensori[k][j] = 0.0f;
        F_tensorj[k][j] = 0.0f;
      }
    }
    F_tensori[0][0] += 1.0f;
    F_tensorj[0][0] += 1.0f;
#endif

    /* compute the contribution from the Eddington tensor to df/dt */
    diffmodeaniso = 2;
    diff_dfrad_term_i[0] = 0.0f;
    diff_dfrad_term_i[1] = 0.0f;
    diff_dfrad_term_i[2] = 0.0f;
    diff_dfrad_term_j[0] = 0.0f;
    diff_dfrad_term_j[1] = 0.0f;
    diff_dfrad_term_j[2] = 0.0f;
    radiation_gradient_aniso_SPH(uradi0, uradj0, mi, mj, rpi->force.f,
                                 rpj->force.f, rhoi, rhoj, wi_dr, wj_dr,
                                 F_tensori, F_tensorj, dx, r, diffmodeaniso,
                                 diff_dfrad_term_i, diff_dfrad_term_j);
    diff_dfrad_term_i[0] *= -cred0 * cred0 / mj;
    diff_dfrad_term_i[1] *= -cred0 * cred0 / mj;
    diff_dfrad_term_i[2] *= -cred0 * cred0 / mj;
    diff_dfrad_term_j[0] *= -cred0 * cred0 / mi;
    diff_dfrad_term_j[1] *= -cred0 * cred0 / mi;
    diff_dfrad_term_j[2] *= -cred0 * cred0 / mi;

    /*******************************/
    /* HERE COME THE CALCULATIONS OF ARTIFICIAL DISSIPATION */
    /*******************************/

    /* fradunit is for anisotropic artificial viscosity */

    if (fradmagi != 0.f) {
      fraduniti[0] = fradi[0] / fradmagi;
      fraduniti[1] = fradi[1] / fradmagi;
      fraduniti[2] = fradi[2] / fradmagi;
    } else {
      fraduniti[0] = 0.0f;
      fraduniti[1] = 0.0f;
      fraduniti[2] = 0.0f;
    }

    if (fradmagj != 0.f) {
      fradunitj[0] = fradj[0] / fradmagj;
      fradunitj[1] = fradj[1] / fradmagj;
      fradunitj[2] = fradj[2] / fradmagj;
    } else {
      fradunitj[0] = 0.0f;
      fradunitj[1] = 0.0f;
      fradunitj[2] = 0.0f;
    }

    vsig_diss_i = credi;
    vsig_diss_j = credj;

    /* compute the artificial diffusion tensor */
    ddi = alpha_diss_i * vsig_diss_i * hi;
    ddj = alpha_diss_j * vsig_diss_j * hj;

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        J_tensori[k][j] = fraduniti[k] * fraduniti[j];
        J_tensorj[k][j] = fradunitj[k] * fradunitj[j];
      }
    }

    /* compute the anisotropic diffusion */
    hid_inv_temp = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
    wi_dr_temp = hid_inv_temp * wi_dx;
    hjd_inv_temp = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
    wj_dr_temp = hjd_inv_temp * wj_dx;
    drhou_low = rhoi * uradi0 - rhoj * uradj0;
    if ((uradi == 0.f) && (uradj == 0.f)) {
      diss_durad_term = 0.0f;
    } else {
      rhomean2 = fminf(rhoi, rhoj) * fminf(rhoi, rhoj);
      diss_durad_term = (drhou_low / rhomean2) * (wi_dr_temp + wj_dr_temp);
      /* TK test: the interpolation is broken: need to fix later. */
      diss_durad_term *= (ddi + ddj);
      diss_durad_term *= 0.5f * r_inv;
    }
    diss_durad_term_i = mj * diss_durad_term;
    diss_durad_term_j = -mi * diss_durad_term;

    /* compute the anisotropic diffusion for radiation flux */
    ddfi = alpha_f_diss_i * vsig_diss_i * hi;
    ddfj = alpha_f_diss_j * vsig_diss_j * hj;

    diffmodeaniso = 2;
    diss_dfrad_term_i[0] = 0.0f;
    diss_dfrad_term_i[1] = 0.0f;
    diss_dfrad_term_i[2] = 0.0f;
    diss_dfrad_term_j[0] = 0.0f;
    diss_dfrad_term_j[1] = 0.0f;
    diss_dfrad_term_j[2] = 0.0f;
    radiation_gradient_aniso_SPH(
        ddfi * divfi, ddfj * divfj, mi, mj, rpi->force.f, rpj->force.f, rhoi,
        rhoj, wi_dr, wj_dr, J_tensori, J_tensorj, dx, r, diffmodeaniso,
        diss_dfrad_term_i, diss_dfrad_term_j);

    /* Assemble the radiation energy equation term */
    durad_dt_i = -divfipar / mj + diss_durad_term_i / mj;
    /* Assemble the radiation flux equation term */
    dfrad_dt_i[0] = diff_dfrad_term_i[0] + diss_dfrad_term_i[0] / mj;
    dfrad_dt_i[1] = diff_dfrad_term_i[1] + diss_dfrad_term_i[1] / mj;
    dfrad_dt_i[2] = diff_dfrad_term_i[2] + diss_dfrad_term_i[2] / mj;

    if (mode == 1) {
      durad_dt_j = -divfjpar / mi + diss_durad_term_j / mi;
      dfrad_dt_j[0] = diff_dfrad_term_j[0] + diss_dfrad_term_j[0] / mi;
      dfrad_dt_j[1] = diff_dfrad_term_j[1] + diss_dfrad_term_j[1] / mi;
      dfrad_dt_j[2] = diff_dfrad_term_j[2] + diss_dfrad_term_j[2] / mi;
    }

    rpi->dconserved_dt[g].urad += mj * durad_dt_i * cred0 / credi;
    rpi->dconserved_dt[g].frad[0] += mj * dfrad_dt_i[0];
    rpi->dconserved_dt[g].frad[1] += mj * dfrad_dt_i[1];
    rpi->dconserved_dt[g].frad[2] += mj * dfrad_dt_i[2];
    if (mode == 1) {
      rpj->dconserved_dt[g].urad += mi * durad_dt_j * cred0 / credj;
      rpj->dconserved_dt[g].frad[0] += mi * dfrad_dt_j[0];
      rpj->dconserved_dt[g].frad[1] += mi * dfrad_dt_j[1];
      rpj->dconserved_dt[g].frad[2] += mi * dfrad_dt_j[2];
    }
  }
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_transport(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  radiation_force_loop_function(r2, dx, hi, hj, pi, pj, a, H, 1);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_transport(float r2, const float *dx, float hi, float hj,
                                struct part *restrict pi,
                                struct part *restrict pj, float a, float H) {

  radiation_force_loop_function(r2, dx, hi, hj, pi, pj, a, H, 0);
}

#endif /* SWIFT_RT_IACT_SPHM1RT_H */
