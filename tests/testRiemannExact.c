/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#include "error.h"
#include "riemann/riemann_exact.h"

/**
 * @brief Check that a and b are consistent (up to some error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(float a, float b, const char* s) {
  if (fabsf(a - b) / fabsf(a + b) > 1.e-5f && fabsf(a - b) > 1.e-5f) {
    error("Values are inconsistent: %g %g (%s)!", a, b, s);
  } else {
    message("Values are consistent: %g %g (%s).", a, b, s);
  }
}

struct riemann_statevector {
  /*! @brief Density */
  float rho;

  /*! @brief Fluid velocity */
  float v;

  /*! @brief Pressure */
  float P;
};

/**
 * @brief Check that the solution to the Riemann problem with given left and
 * right state is consistent with the given expected solution
 *
 * @param WL Left state
 * @param WR Right state
 * @param Whalf Expected solution
 * @param s String used to identify this check in messages
 */
void check_riemann_solution(struct riemann_statevector* WL,
                            struct riemann_statevector* WR,
                            struct riemann_statevector* Whalf, const char* s) {
  float WLarr[5], WRarr[5], Whalfarr[5], n_unit[3];

  n_unit[0] = 1.0f;
  n_unit[1] = 0.0f;
  n_unit[2] = 0.0f;

  WLarr[0] = WL->rho;
  WLarr[1] = WL->v;
  WLarr[2] = 0.0f;
  WLarr[3] = 0.0f;
  WLarr[4] = WL->P;

  WRarr[0] = WR->rho;
  WRarr[1] = WR->v;
  WRarr[2] = 0.0f;
  WRarr[3] = 0.0f;
  WRarr[4] = WR->P;

  riemann_solver_solve(WLarr, WRarr, Whalfarr, n_unit);

  message("Checking %s...", s);
  check_value(Whalfarr[0], Whalf->rho, "rho");
  check_value(Whalfarr[1], Whalf->v, "v");
  check_value(Whalfarr[4], Whalf->P, "P");
}

/**
 * @brief Check the exact Riemann solver on the Toro test problems
 */
void check_riemann_exact() {
  struct riemann_statevector WL, WR, Whalf;

  /* Test 1 */
  WL.rho = 1.0f;
  WL.v = 0.0f;
  WL.P = 1.0f;
  WR.rho = 0.125f;
  WR.v = 0.0f;
  WR.P = 0.1f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 0.47969f;
  Whalf.v = 0.841194f;
  Whalf.P = 0.293945f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.411437f;
  Whalf.v = 0.953205f;
  Whalf.P = 0.306011f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.534767f;
  Whalf.v = 0.760062f;
  Whalf.P = 0.285975f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 1");

  /* Test 2 */
  WL.rho = 1.0f;
  WL.v = -2.0f;
  WL.P = 0.4f;
  WR.rho = 1.0f;
  WR.v = 2.0f;
  WR.P = 0.4f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 0.00617903f;
  Whalf.v = 0.0f;
  Whalf.P = 8.32249e-5f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.0257933f;
  Whalf.v = 0.0f;
  Whalf.P = 0.00304838f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.0f;
  Whalf.v = 0.0f;
  Whalf.P = 0.0f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 2");

  /* Test 3 */
  WL.rho = 1.0f;
  WL.v = 0.0f;
  WL.P = 1000.0f;
  WR.rho = 1.0f;
  WR.v = 0.0f;
  WR.P = 0.01f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 0.615719f;
  Whalf.v = 18.2812f;
  Whalf.P = 445.626f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.563517f;
  Whalf.v = 19.9735f;
  Whalf.P = 465.453f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.656768f;
  Whalf.v = 16.9572f;
  Whalf.P = 431.345f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 3");

  /* Test 4 */
  WL.rho = 1.0f;
  WL.v = 0.0f;
  WL.P = 0.01f;
  WR.rho = 1.0f;
  WR.v = 0.0f;
  WR.P = 100.0f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 0.61577f;
  Whalf.v = -5.78022f;
  Whalf.P = 44.5687f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.563567f;
  Whalf.v = -6.31525f;
  Whalf.P = 46.5508f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.656819f;
  Whalf.v = -5.36146f;
  Whalf.P = 43.1412f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 4");

  /* Test 5 */
  WL.rho = 5.99924f;
  WL.v = 19.5975f;
  WL.P = 460.894f;
  WR.rho = 5.99242f;
  WR.v = -6.19633f;
  WR.P = 46.0950f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 12.743f;
  Whalf.v = 8.56045f;
  Whalf.P = 1841.82f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 5.99924f;
  Whalf.v = 19.5975f;
  Whalf.P = 460.894f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 11.5089f;
  Whalf.v = 8.42099f;
  Whalf.P = 2026.27f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 5");
}

/**
 * @brief Check adiabatic index constants and the various Riemann solvers
 */
int main() {

  /* check the exact Riemann solver */
  check_riemann_exact();

  return 0;
}
