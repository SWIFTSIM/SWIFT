/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

/******************************************************************************
 * The Riemann solver in this file was written by Bert Vandenbroucke as part of
 * the moving mesh code Shadowfax and adapted for use with SWIFT. It consists
 * of an exact Riemann solver as described in
 *  Toro, Eleuterio F., Riemann Solvers and Numerical Methods for Fluid
 *  Dynamics, Springer (2009, 3rd edition)
 *
 ******************************************************************************/
#ifndef SWIFT_RIEMANN_EXACT_H
#define SWIFT_RIEMANN_EXACT_H

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "error.h"
#include "minmax.h"
#include "riemann_checks.h"
#include "riemann_vacuum.h"

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
 * @brief Function (4.5) in Toro
 *
 * @param p The current guess for the pressure
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_f(
    float p, const float* WL, const float* WR, float vL, float vR, float aL,
    float aR) {

  return riemann_fb(p, WL, aL) + riemann_fb(p, WR, aR) + (vR - vL);
}

/**
 * @brief Function (4.37) in Toro
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 * @param a The left or right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_fprimeb(
    float p, const float* W, float a) {

  float fval;
  if (p > W[4]) {
    const float A = hydro_two_over_gamma_plus_one / W[0];
    const float B = hydro_gamma_minus_one_over_gamma_plus_one * W[4];
    fval = (1.0f - 0.5f * (p - W[4]) / (B + p)) * sqrtf(A / (p + B));
  } else {
    fval = 1.0f / W[0] / a * pow_minus_gamma_plus_one_over_two_gamma(p / W[4]);
  }
  return fval;
}

/**
 * @brief The derivative of riemann_f w.r.t. p
 *
 * @param p The current guess for the pressure
 * @param WL The left state vector
 * @param WR The right state vector
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_fprime(
    float p, const float* WL, const float* WR, float aL, float aR) {

  return riemann_fprimeb(p, WL, aL) + riemann_fprimeb(p, WR, aR);
}

/**
 * @brief Bottom function of (4.48) in Toro
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 */
__attribute__((always_inline)) INLINE static float riemann_gb(float p,
                                                              const float* W) {

  const float A = hydro_two_over_gamma_plus_one / W[0];
  const float B = hydro_gamma_minus_one_over_gamma_plus_one * W[4];
  return sqrtf(A / (p + B));
}

/**
 * @brief Get a good first guess for the pressure in the iterative scheme
 *
 * This function is based on (4.47) and (4.48) in Toro and on the
 * FORTRAN code provided in Toro p.156-157
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_guess_p(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR) {

  float pguess, pmin, pmax, qmax;
  float ppv;

  pmin = min(WL[4], WR[4]);
  pmax = max(WL[4], WR[4]);
  qmax = pmax / pmin;
  ppv =
      0.5f * (WL[4] + WR[4]) - 0.125f * (vR - vL) * (WL[0] + WR[0]) * (aL + aR);
  ppv = max(1.e-8f, ppv);
  if (qmax <= 2.0f && pmin <= ppv && ppv <= pmax) {
    pguess = ppv;
  } else {
    if (ppv < pmin) {
      /* two rarefactions */
      pguess = pow_two_gamma_over_gamma_minus_one(
          (aL + aR - hydro_gamma_minus_one_over_two * (vR - vL)) /
          (aL / pow_gamma_minus_one_over_two_gamma(WL[4]) +
           aR / pow_gamma_minus_one_over_two_gamma(WR[4])));
    } else {
      /* two shocks */
      pguess = (riemann_gb(ppv, WL) * WL[4] + riemann_gb(ppv, WR) * WR[4] - vR +
                vL) /
               (riemann_gb(ppv, WL) + riemann_gb(ppv, WR));
    }
  }
  /* Toro: "Not that approximate solutions may predict, incorrectly, a negative
     value for pressure (...).
     Thus in order to avoid negative guess values we introduce the small
     positive constant _tolerance" */
  pguess = max(1.e-8f, pguess);
  return pguess;
}

/**
 * @brief Find the zeropoint of riemann_f(p) using Brent's method
 *
 * @param lower_limit Lower limit for the method (riemann_f(lower_limit) < 0)
 * @param upper_limit Upper limit for the method (riemann_f(upper_limit) > 0)
 * @param lowf Function value riemann_f(lower_limit)
 * @param upf  Function value riemann_f(upper_limit)
 * @param error_tol Tolerance used to decide if the solution is converged
 * @param WL Left state vector
 * @param WR Right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
__attribute__((always_inline)) INLINE static float riemann_solve_brent(
    float lower_limit, float upper_limit, float lowf, float upf,
    float error_tol, const float* WL, const float* WR, float vL, float vR,
    float aL, float aR) {

  float a, b, c, d, e, s;
  float fa, fb, fc, fs;
  float tmp, tmp2;
  int mflag;

  a = lower_limit;
  b = upper_limit;
  c = 0.0f;     // previous value of b: b_{n-1}
  d = FLT_MAX;  // value of b from two iterations ago: b_{n-2}
  e = FLT_MAX;  // previous value of a: a_{n-1}

  fa = lowf;
  fb = upf;

  fc = 0.0f;
  s = 0.0f;
  fs = 0.0f;

  /* if f(a) f(b) >= 0 then error-exit */
  if (fa * fb >= 0.0f) {
    error(
        "Brent's method called with equal sign function values!\n"
        "f(%g) = %g, f(%g) = %g\n",
        a, fa, b, fb);
    /* return NaN */
    return 0.0f / 0.0f;
  }

  /* if |f(a)| < |f(b)| then swap (a,b) */
  if (fabs(fa) < fabs(fb)) {
    tmp = a;
    a = b;
    b = tmp;
    tmp = fa;
    fa = fb;
    fb = tmp;
  }

  c = a;
  fc = fa;
  mflag = 1;

  /* Loop until convergence, i.e. until an exact zero point is found, or the
   * interval is sufficiently small, or the interval is unchanged since the
   * previous iteration */
  int counter = 0;
  while ((fb != 0.0f) && (fabs(a - b) > error_tol * 0.5f * (a + b)) &&
         (a != e || b != c)) {
    counter++;
    if (counter > 1000) error("Brent's method did not converge!\n");

    if ((fa != fc) && (fb != fc)) /* Inverse quadratic interpolation */
      s = a * fb * fc / (fa - fb) / (fa - fc) +
          b * fa * fc / (fb - fa) / (fb - fc) +
          c * fa * fb / (fc - fa) / (fc - fb);
    else
      /* Secant Rule */
      s = b - fb * (b - a) / (fb - fa);

    tmp2 = 0.25f * (3.0f * a + b);
    if (!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b))) ||
        (mflag && (fabs(s - b) >= (0.5f * fabs(b - c)))) ||
        (!mflag && (fabs(s - b) >= (0.5f * fabs(c - d)))) ||
        (mflag && (fabs(b - c) < error_tol)) ||
        (!mflag && (fabs(c - d) < error_tol))) {
      s = 0.5f * (a + b);
      mflag = 1;
    } else {
      mflag = 0;
    }
    fs = riemann_f(s, WL, WR, vL, vR, aL, aR);
    d = c;
    c = b;
    e = a;
    fc = fb;
    if (fa * fs < 0.) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    /* if |f(a)| < |f(b)| then swap (a,b) */
    if (fabs(fa) < fabs(fb)) {
      tmp = a;
      a = b;
      b = tmp;
      tmp = fa;
      fa = fb;
      fb = tmp;
    }
  }
  return b;
}

/* Solve the Riemann problem between the states WL and WR and store the result
 * in Whalf
 * The Riemann problem is solved in the x-direction; the velocities in the y-
 * and z-direction
 * are simply advected.
 */
/**
 * @brief Solve the Riemann problem between the given left and right state and
 * along the given interface normal
 *
 * Based on chapter 4 in Toro
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param Whalf Empty state vector in which the result will be stored
 * @param n_unit Normal vector of the interface
 */
__attribute__((always_inline)) INLINE static void riemann_solver_solve(
    const float* WL, const float* WR, float* Whalf, const float* n_unit) {

  /* velocity of the left and right state in a frame aligned with n_unit */
  float vL, vR, vhalf;
  /* sound speeds */
  float aL, aR;
  /* variables used for finding pstar */
  float p, pguess, fp, fpguess;
  /* variables used for sampling the solution */
  float u;
  float pdpR, SR;
  float SHR, STR;
  float pdpL, SL;
  float SHL, STL;

  /* calculate velocities in interface frame */
  vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* calculate sound speeds */
  aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* check vacuum (generation) condition */
  if (riemann_is_vacuum(WL, WR, vL, vR, aL, aR)) {
    riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
    return;
  }

  /* values are ok: let's find pstar (riemann_f(pstar) = 0)! */
  /* We normally use a Newton-Raphson iteration to find the zeropoint
     of riemann_f(p), but if pstar is close to 0, we risk negative p values.
     Since riemann_f(p) is undefined for negative pressures, we don't
     want this to happen.
     We therefore use Brent's method if riemann_f(0) is larger than some
     value. -5 makes the iteration fail safe while almost never invoking
     the expensive Brent solver. */
  p = 0.;
  /* obtain a first guess for p */
  pguess = riemann_guess_p(WL, WR, vL, vR, aL, aR);
  fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
  fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
  /* ok, pstar is close to 0, better use Brent's method... */
  /* we use Newton-Raphson until we find a suitable interval */
  if (fp * fpguess >= 0.0f) {
    /* Newton-Raphson until convergence or until suitable interval is found
       to use Brent's method */
    unsigned int counter = 0;
    while (fabs(p - pguess) > 1.e-6f * 0.5f * (p + pguess) && fpguess < 0.0f) {
      p = pguess;
      pguess = pguess - fpguess / riemann_fprime(pguess, WL, WR, aL, aR);
      fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
      counter++;
      if (counter > 1000) {
        error("Stuck in Newton-Raphson!\n");
      }
    }
  }
  /* As soon as there is a suitable interval: use Brent's method */
  if (1.e6 * fabs(p - pguess) > 0.5f * (p + pguess) && fpguess > 0.0f) {
    p = 0.0f;
    fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
    /* use Brent's method to find the zeropoint */
    p = riemann_solve_brent(p, pguess, fp, fpguess, 1.e-6, WL, WR, vL, vR, aL,
                            aR);
  } else {
    p = pguess;
  }

  /* calculate the velocity in the intermediate state */
  u = 0.5f * (vL + vR) + 0.5f * (riemann_fb(p, WR, aR) - riemann_fb(p, WL, aL));

  /* sample the solution */
  /* This corresponds to the flow chart in Fig. 4.14 in Toro */
  if (u < 0.0f) {
    /* advect velocity components */
    Whalf[1] = WR[1];
    Whalf[2] = WR[2];
    Whalf[3] = WR[3];
    pdpR = p / WR[4];
    if (p > WR[4]) {
      /* shockwave */
      SR = vR + aR * sqrtf(hydro_gamma_plus_one_over_two_gamma * pdpR +
                           hydro_gamma_minus_one_over_two_gamma);
      if (SR > 0.0f) {
        Whalf[0] = WR[0] * (pdpR + hydro_gamma_minus_one_over_gamma_plus_one) /
                   (hydro_gamma_minus_one_over_gamma_plus_one * pdpR + 1.0f);
        vhalf = u - vR;
        Whalf[4] = p;
      } else {
        Whalf[0] = WR[0];
        vhalf = 0.0f;
        Whalf[4] = WR[4];
      }
    } else {
      /* rarefaction wave */
      SHR = vR + aR;
      if (SHR > 0.0f) {
        STR = u + aR * pow_gamma_minus_one_over_two_gamma(pdpR);
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
          vhalf = u - vR;
          Whalf[4] = p;
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
    pdpL = p / WL[4];
    if (p > WL[4]) {
      /* shockwave */
      SL = vL - aL * sqrtf(hydro_gamma_plus_one_over_two_gamma * pdpL +
                           hydro_gamma_minus_one_over_two_gamma);
      if (SL < 0.0f) {
        Whalf[0] = WL[0] * (pdpL + hydro_gamma_minus_one_over_gamma_plus_one) /
                   (hydro_gamma_minus_one_over_gamma_plus_one * pdpL + 1.0f);
        vhalf = u - vL;
        Whalf[4] = p;
      } else {
        Whalf[0] = WL[0];
        vhalf = 0.0f;
        Whalf[4] = WL[4];
      }
    } else {
      /* rarefaction wave */
      SHL = vL - aL;
      if (SHL < 0.0f) {
        STL = u - aL * pow_gamma_minus_one_over_two_gamma(pdpL);
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
          vhalf = u - vL;
          Whalf[4] = p;
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

/**
 * @brief Solve the Riemann problem between the given left and right state and
 * return the velocity and pressure in the middle state
 *
 * Based on chapter 4 in Toro
 *
 * @param WL Left state.
 * @param vL Left state velocity.
 * @param WR Right state.
 * @param vR Right state velocity.
 * @param vM Middle state velocity.
 * @param PM Middle state pressure.
 */
__attribute__((always_inline)) INLINE static void
riemann_solver_solve_middle_state(const float* WL, const float vL,
                                  const float* WR, const float vR, float* vM,
                                  float* PM) {

  /* sound speeds */
  float aL, aR;
  /* variables used for finding pstar */
  float p, pguess, fp, fpguess;

  /* calculate sound speeds */
  aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* vacuum */
  /* check vacuum (generation) condition */
  if (riemann_is_vacuum(WL, WR, vL, vR, aL, aR)) {
    *vM = 0.0f;
    *PM = 0.0f;
    return;
  }

  /* values are ok: let's find pstar (riemann_f(pstar) = 0)! */
  /* We normally use a Newton-Raphson iteration to find the zeropoint
     of riemann_f(p), but if pstar is close to 0, we risk negative p values.
     Since riemann_f(p) is undefined for negative pressures, we don't
     want this to happen.
     We therefore use Brent's method if riemann_f(0) is larger than some
     value. -5 makes the iteration fail safe while almost never invoking
     the expensive Brent solver. */
  p = 0.0f;
  /* obtain a first guess for p */
  pguess = riemann_guess_p(WL, WR, vL, vR, aL, aR);
  fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
  fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
  /* ok, pstar is close to 0, better use Brent's method... */
  /* we use Newton-Raphson until we find a suitable interval */
  if (fp * fpguess >= 0.0f) {
    /* Newton-Raphson until convergence or until suitable interval is found
       to use Brent's method */
    unsigned int counter = 0;
    while (fabs(p - pguess) > 1.e-6f * 0.5f * (p + pguess) && fpguess < 0.0f) {
      p = pguess;
      pguess = pguess - fpguess / riemann_fprime(pguess, WL, WR, aL, aR);
      fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
      counter++;
      if (counter > 1000) {
        error("Stuck in Newton-Raphson!\n");
      }
    }
  }
  /* As soon as there is a suitable interval: use Brent's method */
  if (1.e6 * fabs(p - pguess) > 0.5f * (p + pguess) && fpguess > 0.0f) {
    p = 0.0f;
    fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
    /* use Brent's method to find the zeropoint */
    p = riemann_solve_brent(p, pguess, fp, fpguess, 1.e-6, WL, WR, vL, vR, aL,
                            aR);
  } else {
    p = pguess;
  }

  *PM = p;
  /* calculate the velocity in the intermediate state */
  *vM =
      0.5f * (vL + vR) + 0.5f * (riemann_fb(p, WR, aR) - riemann_fb(p, WL, aL));
}

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    const float* Wi, const float* Wj, const float* n_unit, const float* vij,
    float* totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(Wi, Wj, n_unit, vij);
#endif

  float Whalf[5];
  float flux[5][3];
  float vtot[3];
  float rhoe;

  riemann_solver_solve(Wi, Wj, Whalf, n_unit);

  flux[0][0] = Whalf[0] * Whalf[1];
  flux[0][1] = Whalf[0] * Whalf[2];
  flux[0][2] = Whalf[0] * Whalf[3];

  vtot[0] = Whalf[1] + vij[0];
  vtot[1] = Whalf[2] + vij[1];
  vtot[2] = Whalf[3] + vij[2];
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
  rhoe = Whalf[4] / hydro_gamma_minus_one +
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

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(Wi, Wj, n_unit, vij, totflux);
#endif
}

__attribute__((always_inline)) INLINE static void
riemann_solve_for_middle_state_flux(const float* Wi, const float* Wj,
                                    const float* n_unit, const float* vij,
                                    float* totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(Wi, Wj, n_unit, vij);
#endif

  /* vacuum? */
  if (Wi[0] == 0.0f || Wj[0] == 0.0f) {
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  const float vL = Wi[1] * n_unit[0] + Wi[2] * n_unit[1] + Wi[3] * n_unit[2];
  const float vR = Wj[1] * n_unit[0] + Wj[2] * n_unit[1] + Wj[3] * n_unit[2];

  float vM, PM;
  riemann_solver_solve_middle_state(Wi, vL, Wj, vR, &vM, &PM);

  const float vface =
      vij[0] * n_unit[0] + vij[1] * n_unit[1] + vij[2] * n_unit[2];

  totflux[0] = 0.0f;
  totflux[1] = PM * n_unit[0];
  totflux[2] = PM * n_unit[1];
  totflux[3] = PM * n_unit[2];
  totflux[4] = (vM + vface) * PM;

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(Wi, Wj, n_unit, vij, totflux);
#endif
}

#endif /* SWIFT_RIEMANN_EXACT_H */
