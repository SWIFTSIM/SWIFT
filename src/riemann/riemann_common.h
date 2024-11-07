//
// Created by yuyttenh on 28/02/23.
//

#ifndef SWIFT_RIEMANN_COMMON_H
#define SWIFT_RIEMANN_COMMON_H

/**
 * @brief Functions (4.6) and (4.7) in Toro.
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 * @param a The left or right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_fb(float p,
                                                              const float* W,
                                                              float a) {

  float fval;
  if (p > W[4]) {
    const float A = hydro_two_over_gamma_plus_one / W[0];
    const float B = hydro_gamma_minus_one_over_gamma_plus_one * W[4];
    fval = (p - W[4]) * sqrtf(A / (p + B));
  } else {
    fval = hydro_two_over_gamma_minus_one * a *
           (pow_gamma_minus_one_over_two_gamma(p / W[4]) - 1.0f);
  }
  return fval;
}

/**
 * @brief Compute the maximal mach number of the shock waves in the Riemann
 * problem (if any).
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param P_star The middle state pressure.
 * @returns The maximal mach number.
 */
__attribute__((always_inline)) INLINE static float riemann_get_max_mach_number(
    const float* restrict WL, const float* restrict WR, float P_star) {
  float mach = 0.f;

  /* left shock? */
  if (P_star > WL[4]) {
    float p_frac = P_star / WL[4];
    mach = fmaxf(mach, sqrtf(hydro_gamma_plus_one_over_two_gamma * p_frac +
                             hydro_gamma_minus_one_over_two_gamma));
  }
  /* right shock? */
  if (P_star > WR[4]) {
    float p_frac = P_star / WR[4];
    mach = fmaxf(mach, sqrtf(hydro_gamma_plus_one_over_two_gamma * p_frac +
                             hydro_gamma_minus_one_over_two_gamma));
  }

  return mach;
}

/**
 * @brief Sample the solution for the given value of pstar, assuming the
 * relations of the exact Riemann solver, in the x / t = 0 direction (i.e. in
 * the rest frame of the interface).
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param n_unit The unit vector of the interface between the left and right
 * states.
 * @param vL, vR The left and right velocities along the interface normal.
 * @param aL, aR The left and right sound speeds.
 * @param P_star The middle state pressure.
 * @param Whalf (return) The solution at x/t = 0.
 */
__attribute__((always_inline)) INLINE static void riemann_sample(
    const float* restrict WL, const float* restrict WR,
    const float* restrict n_unit, float vL, float vR, float aL, float aR,
    float P_star, float* restrict Whalf) {
  /* calculate the velocity in the intermediate state */
  float v_star =
      0.5f * (vL + vR) + 0.5f * (riemann_fb(P_star, WR, aR) - riemann_fb(P_star, WL, aL));

  /* variables used for sampling the solution */
  float vhalf;
  float pdpR, SR;
  float SHR, STR;
  float pdpL, SL;
  float SHL, STL;

  /* sample the solution */
  /* This corresponds to the flow chart in Fig. 4.14 in Toro */
  if (v_star < 0.0f) {
    /* advect velocity components */
    Whalf[1] = WR[1];
    Whalf[2] = WR[2];
    Whalf[3] = WR[3];
    pdpR = P_star / WR[4];
    if (P_star > WR[4]) {
      /* shockwave */
      SR = vR + aR * sqrtf(hydro_gamma_plus_one_over_two_gamma * pdpR +
                           hydro_gamma_minus_one_over_two_gamma);
      if (SR > 0.0f) {
        Whalf[0] = WR[0] * (pdpR + hydro_gamma_minus_one_over_gamma_plus_one) /
                   (hydro_gamma_minus_one_over_gamma_plus_one * pdpR + 1.0f);
        vhalf = v_star - vR;
        Whalf[4] = P_star;
      } else {
        Whalf[0] = WR[0];
        vhalf = 0.0f;
        Whalf[4] = WR[4];
      }
    } else {
      /* rarefaction wave */
      SHR = vR + aR;
      if (SHR > 0.0f) {
        STR = v_star + aR * pow_gamma_minus_one_over_two_gamma(pdpR);
        if (STR <= 0.0f) {
          Whalf[0] =
              WR[0] * pow_two_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one -
                          hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
          vhalf = hydro_two_over_gamma_plus_one *
                      (-aR + hydro_gamma_minus_one_over_two * vR) -
                  vR;
          Whalf[4] =
              WR[4] * pow_two_gamma_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one -
                          hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
        } else {
          Whalf[0] = WR[0] * pow_one_over_gamma(pdpR);
          vhalf = v_star - vR;
          Whalf[4] = P_star;
        }
      } else {
        Whalf[0] = WR[0];
        vhalf = 0.0f;
        Whalf[4] = WR[4];
      }
    }
  } else {
    Whalf[1] = WL[1];
    Whalf[2] = WL[2];
    Whalf[3] = WL[3];
    pdpL = P_star / WL[4];
    if (P_star > WL[4]) {
      /* shockwave */
      SL = vL - aL * sqrtf(hydro_gamma_plus_one_over_two_gamma * pdpL +
                           hydro_gamma_minus_one_over_two_gamma);
      if (SL < 0.0f) {
        Whalf[0] = WL[0] * (pdpL + hydro_gamma_minus_one_over_gamma_plus_one) /
                   (hydro_gamma_minus_one_over_gamma_plus_one * pdpL + 1.0f);
        vhalf = v_star - vL;
        Whalf[4] = P_star;
      } else {
        Whalf[0] = WL[0];
        vhalf = 0.0f;
        Whalf[4] = WL[4];
      }
    } else {
      /* rarefaction wave */
      SHL = vL - aL;
      if (SHL < 0.0f) {
        STL = v_star - aL * pow_gamma_minus_one_over_two_gamma(pdpL);
        if (STL > 0.0f) {
          Whalf[0] =
              WL[0] * pow_two_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one +
                          hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
          vhalf = hydro_two_over_gamma_plus_one *
                      (aL + hydro_gamma_minus_one_over_two * vL) -
                  vL;
          Whalf[4] =
              WL[4] * pow_two_gamma_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one +
                          hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
        } else {
          Whalf[0] = WL[0] * pow_one_over_gamma(pdpL);
          vhalf = v_star - vL;
          Whalf[4] = P_star;
        }
      } else {
        Whalf[0] = WL[0];
        vhalf = 0.0f;
        Whalf[4] = WL[4];
      }
    }
  }

  /* add the velocity solution along the interface normal to the velocities */
  Whalf[1] += vhalf * n_unit[0];
  Whalf[2] += vhalf * n_unit[1];
  Whalf[3] += vhalf * n_unit[2];
}

__attribute__((always_inline)) INLINE static void riemann_flux_from_half_state(
    const float* Whalf, const float* vij, const float* n_unit, float* totflux) {

  float vtot[3];
  vtot[0] = Whalf[1] + vij[0];
  vtot[1] = Whalf[2] + vij[1];
  vtot[2] = Whalf[3] + vij[2];

  float flux[5][3];
  flux[0][0] = Whalf[0] * Whalf[1];
  flux[0][1] = Whalf[0] * Whalf[2];
  flux[0][2] = Whalf[0] * Whalf[3];

  flux[1][0] = Whalf[0] * vtot[0] * Whalf[1] + Whalf[4];
  flux[1][1] = Whalf[0] * vtot[0] * Whalf[2];
  flux[1][2] = Whalf[0] * vtot[0] * Whalf[3];
  flux[2][0] = Whalf[0] * vtot[1] * Whalf[1];
  flux[2][1] = Whalf[0] * vtot[1] * Whalf[2] + Whalf[4];
  flux[2][2] = Whalf[0] * vtot[1] * Whalf[3];
  flux[3][0] = Whalf[0] * vtot[2] * Whalf[1];
  flux[3][1] = Whalf[0] * vtot[2] * Whalf[2];
  flux[3][2] = Whalf[0] * vtot[2] * Whalf[3] + Whalf[4];

  /* eqn. (15) */
  /* F_P = \rho e ( \vec{v} - \vec{v_{ij}} ) + P \vec{v} */
  /* \rho e = P / (\gamma-1) + 1/2 \rho \vec{v}^2 */
  float rhoe = Whalf[4] / hydro_gamma_minus_one +
               0.5f * Whalf[0] *
                   (vtot[0] * vtot[0] + vtot[1] * vtot[1] + vtot[2] * vtot[2]);
  flux[4][0] = rhoe * Whalf[1] + Whalf[4] * vtot[0];
  flux[4][1] = rhoe * Whalf[2] + Whalf[4] * vtot[1];
  flux[4][2] = rhoe * Whalf[3] + Whalf[4] * vtot[2];

  totflux[0] =
      flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2];
  totflux[1] =
      flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2];
  totflux[2] =
      flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2];
  totflux[3] =
      flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2];
  totflux[4] =
      flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2];
}

#endif  // SWIFT_RIEMANN_COMMON_H
