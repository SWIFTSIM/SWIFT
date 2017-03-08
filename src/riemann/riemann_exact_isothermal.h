/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

#ifndef SWIFT_RIEMANN_EXACT_ISOTHERMAL_H
#define SWIFT_RIEMANN_EXACT_ISOTHERMAL_H

#include <float.h>
#include "adiabatic_index.h"
#include "minmax.h"
#include "riemann_vacuum.h"

#define const_isothermal_soundspeed \
  sqrtf(hydro_gamma_minus_one* const_isothermal_internal_energy)

/**
 * @brief Relative difference between the middle state velocity and the left or
 * right state velocity used in the middle state density iteration.
 *
 * @param rho Current estimate of the middle state density.
 * @param W Left or right state vector.
 * @return Density dependent part of the middle state velocity.
 */
__attribute__((always_inline)) INLINE static float riemann_fb(float rho,
                                                              float* W) {
  if (rho < W[0]) {
    return const_isothermal_soundspeed * logf(rho / W[0]);
  } else {
    return const_isothermal_soundspeed *
           (sqrtf(rho / W[0]) - sqrtf(W[0] / rho));
  }
}

/**
 * @brief Derivative w.r.t. rho of the function riemann_fb.
 *
 * @param rho Current estimate of the middle state density.
 * @param W Left or right state vector.
 * @return Derivative of riemann_fb.
 */
__attribute__((always_inline)) INLINE static float riemann_fprimeb(float rho,
                                                                   float* W) {
  if (rho < W[0]) {
    return const_isothermal_soundspeed * W[0] / rho;
  } else {
    return 0.5 * const_isothermal_soundspeed *
           (sqrtf(rho / W[0]) + sqrtf(W[0] / rho)) / rho;
  }
}

/**
 * @brief Difference between the left and right middle state velocity estimates.
 *
 * Since the middle state velocity takes on a constant value, we want to get
 * this difference as close to zero as possible.
 *
 * @param rho Current estimate of the middle state density.
 * @param WL Left state vector.
 * @param WR Right state vector.
 * @param vL Left state velocity along the interface normal.
 * @param vR Right state velocity along the interface normal.
 * @return Difference between the left and right middle state velocity
 * estimates.
 */
__attribute__((always_inline)) INLINE static float riemann_f(
    float rho, float* WL, float* WR, float vL, float vR) {
  return riemann_fb(rho, WR) + riemann_fb(rho, WL) + vR - vL;
}

/**
 * @brief Derivative of riemann_f w.r.t. rho.
 *
 * @param rho Current estimate of the middle state density.
 * @param WL Left state vector.
 * @param WR Right state vector.
 * @return Derivative of riemann_f.
 */
__attribute__((always_inline)) INLINE static float riemann_fprime(float rho,
                                                                  float* WL,
                                                                  float* WR) {
  return riemann_fprimeb(rho, WL) + riemann_fprimeb(rho, WR);
}

/**
 * @brief Get a good first guess for the middle state density.
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 */
__attribute__((always_inline)) INLINE static float riemann_guess_rho(float* WL,
                                                                     float* WR,
                                                                     float vL,
                                                                     float vR) {

  /* Currently three possibilities and not really an algorithm to decide which
     one to choose: */
  /* just the average */
  //  return 0.5f * (WL[0] + WR[0]);

  /* two rarefaction approximation */
  return sqrtf(WL[0] * WR[0] * expf((vL - vR) / const_isothermal_soundspeed));

  /* linearized primitive variable approximation */
  return 0.25f * (WL[0] + WR[0]) * (vL - vR) / const_isothermal_soundspeed +
         0.5f * (WL[0] + WR[0]);
}

/**
 * @brief Find the zeropoint of riemann_f(rho) using Brent's method.
 *
 * @param lower_limit Lower limit for the method (riemann_f(lower_limit) < 0)
 * @param upper_limit Upper limit for the method (riemann_f(upper_limit) > 0)
 * @param lowf Value of riemann_f(lower_limit).
 * @param upf  Value of riemann_f(upper_limit).
 * @param error_tol Tolerance used to decide if the solution is converged
 * @param WL Left state vector
 * @param WR Right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 */
__attribute__((always_inline)) INLINE static float riemann_solve_brent(
    float lower_limit, float upper_limit, float lowf, float upf,
    float error_tol, float* WL, float* WR, float vL, float vR) {

  float a, b, c, d, s;
  float fa, fb, fc, fs;
  float tmp, tmp2;
  int mflag;
  int i;

  a = lower_limit;
  b = upper_limit;
  c = 0.0f;
  d = FLT_MAX;

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
  i = 0;

  while (!(fb == 0.0f) && (fabs(a - b) > error_tol * 0.5f * (a + b))) {
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
    fs = riemann_f(s, WL, WR, vL, vR);
    d = c;
    c = b;
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
    i++;
  }
  return b;
}

/**
 * @brief Solve the Riemann problem between the given left and right state and
 * along the given interface normal
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param Whalf Empty state vector in which the result will be stored
 * @param n_unit Normal vector of the interface
 */
__attribute__((always_inline)) INLINE static void riemann_solver_solve(
    float* WL, float* WR, float* Whalf, float* n_unit) {

  /* velocity of the left and right state in a frame aligned with n_unit */
  float vL, vR, vhalf;
  /* variables used for finding rhostar */
  float rho, rhoguess, frho, frhoguess;
  /* variables used for sampling the solution */
  float u, S, SH, ST;

  int errorFlag = 0;

  /* sanity checks */
  if (WL[0] != WL[0]) {
    printf("NaN WL!\n");
    errorFlag = 1;
  }
  if (WR[0] != WR[0]) {
    printf("NaN WR!\n");
    errorFlag = 1;
  }
  if (WL[0] < 0.0f) {
    printf("Negative WL!\n");
    errorFlag = 1;
  }
  if (WR[0] < 0.0f) {
    printf("Negative WR!\n");
    errorFlag = 1;
  }
  if (errorFlag) {
    printf("WL: %g %g %g %g %g\n", WL[0], WL[1], WL[2], WL[3], WL[4]);
    printf("WR: %g %g %g %g %g\n", WR[0], WR[1], WR[2], WR[3], WR[4]);
    error("Riemman solver input error!\n");
  }

  /* calculate velocities in interface frame */
  vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  if (WL[0] == 0. || WR[0] == 0.) {
    error(
        "One of the states is vacuum, the isothermal solver cannot solve "
        "this!");
  }

  rho = 0.;
  /* obtain a first guess for p */
  rhoguess = riemann_guess_rho(WL, WR, vL, vR);

#ifdef SWIFT_DEBUG_CHECKS
  if (rhoguess <= 0.) {
    error("Zero or negative initial density guess.");
  }
#endif

  /* we know that the value of riemann_f for rho=0 is negative (it is -inf). */
  frho = -1.;
  frhoguess = riemann_f(rhoguess, WL, WR, vL, vR);
  /* ok, rhostar is close to 0, better use Brent's method... */
  /* we use Newton-Raphson until we find a suitable interval */
  if (frho * frhoguess >= 0.0f) {
    /* Newton-Raphson until convergence or until suitable interval is found
       to use Brent's method */
    unsigned int counter = 0;
    while (fabs(rho - rhoguess) > 5.e-7f * (rho + rhoguess) &&
           frhoguess < 0.0f) {
      rho = rhoguess;
      rhoguess = rhoguess - frhoguess / riemann_fprime(rhoguess, WL, WR);
      frhoguess = riemann_f(rhoguess, WL, WR, vL, vR);
      counter++;
      if (counter > 1000) {
        error(
            "Stuck in Newton-Raphson (rho: %g, rhoguess: %g, frhoguess: %g, "
            "fprime: %g, rho-rhoguess: %g, WL: %g %g %g, WR: %g %g %g)!\n",
            rho, rhoguess, frhoguess, riemann_fprime(rhoguess, WL, WR),
            (rho - rhoguess), WL[0], vL, WL[4], WR[0], vR, WR[4]);
      }
    }
  }
  /* As soon as there is a suitable interval: use Brent's method */
  if (1.e6 * fabs(rho - rhoguess) > 0.5f * (rho + rhoguess) &&
      frhoguess > 0.0f) {
    rho = 0.0f;
    frho = -1.;
    /* use Brent's method to find the zeropoint */
    rho = riemann_solve_brent(rho, rhoguess, frho, frhoguess, 1.e-6, WL, WR, vL,
                              vR);
  } else {
    rho = rhoguess;
  }

  /* calculate the middle state velocity */
  u = 0.5f * (vL - riemann_fb(rho, WL) + vR + riemann_fb(rho, WR));

  /* sample the solution */
  if (u > 0.0f) {
    /* left state */
    Whalf[1] = WL[1];
    Whalf[2] = WL[2];
    Whalf[3] = WL[3];
    if (WL[0] < rho) {
      /* left shock wave */
      S = vL - const_isothermal_soundspeed * sqrtf(rho / WL[0]);
      if (S >= 0.) {
        /* to the left of the shock */
        Whalf[0] = WL[0];
        vhalf = 0.0f;
      } else {
        /* to the right of the shock */
        Whalf[0] = rho;
        vhalf = u - vL;
      }
    } else {
      /* left rarefaction wave */
      SH = vL - const_isothermal_soundspeed;
      ST = u - const_isothermal_soundspeed;
      if (SH > 0.) {
        /* to the left of the rarefaction */
        Whalf[0] = WL[0];
        vhalf = 0.0f;
      } else if (ST > 0.0f) {
        /* inside the rarefaction */
        Whalf[0] = WL[0] * expf(vL / const_isothermal_soundspeed - 1.0f);
        vhalf = const_isothermal_soundspeed - vL;
      } else {
        /* to the right of the rarefaction */
        Whalf[0] = rho;
        vhalf = u - vL;
      }
    }
  } else {
    /* right state */
    Whalf[1] = WR[1];
    Whalf[2] = WR[2];
    Whalf[3] = WR[3];
    if (WR[0] < rho) {
      /* right shock wave */
      S = vR + const_isothermal_soundspeed * sqrtf(rho / WR[0]);
      if (S > 0.0f) {
        /* to the left of the shock wave: middle state */
        Whalf[0] = rho;
        vhalf = u - vR;
      } else {
        /* to the right of the shock wave: right state */
        Whalf[0] = WR[0];
        vhalf = 0.0f;
      }
    } else {
      /* right rarefaction wave */
      SH = vR + const_isothermal_soundspeed;
      ST = u + const_isothermal_soundspeed;
      if (ST > 0.0f) {
        /* to the left of rarefaction: middle state */
        Whalf[0] = rho;
        vhalf = u - vR;
      } else if (SH > 0.0f) {
        /* inside rarefaction */
        Whalf[0] = WR[0] * expf(-vR / const_isothermal_soundspeed - 1.0f);
        vhalf = -const_isothermal_soundspeed - vR;
      } else {
        /* to the right of rarefaction: right state */
        Whalf[0] = WR[0];
        vhalf = 0.0f;
      }
    }
  }

  /* add the velocity solution along the interface normal to the velocities */
  Whalf[1] += vhalf * n_unit[0];
  Whalf[2] += vhalf * n_unit[1];
  Whalf[3] += vhalf * n_unit[2];

  /* the pressure is completely irrelevant in this case */
  Whalf[4] =
      Whalf[0] * const_isothermal_soundspeed * const_isothermal_soundspeed;
}

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    float* Wi, float* Wj, float* n_unit, float* vij, float* totflux) {

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
}

#endif /* SWIFT_RIEMANN_EXACT_ISOTHERMAL_H */
