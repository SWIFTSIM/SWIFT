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
#include "../config.h"

/* Local headers. */
#include <string.h>

/* Local includes */
#include "error.h"
#include "riemann/riemann_trrs.h"
#include "tools.h"

int opposite(float a, float b) {
  if ((a - b)) {
    return fabs((a + b) / (a - b)) < 1.e-4;
  } else {
    return a == 0.0f;
  }
}

int equal(float a, float b) {
  if ((a + b)) {
    return fabs((a - b) / (a + b)) < 1.e-4;
  } else {
    return a == 0.0f;
  }
}

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
 * @brief Check the TRRS Riemann solver on the Toro test problems
 */
void check_riemann_trrs(void) {
  struct riemann_statevector WL, WR, Whalf;

  /* Test 1 */
  WL.rho = 1.0f;
  WL.v = 0.0f;
  WL.P = 1.0f;
  WR.rho = 0.125f;
  WR.v = 0.0f;
  WR.P = 0.1f;
#if defined(HYDRO_GAMMA_5_3)
  Whalf.rho = 0.481167f;
  Whalf.v = 0.838085f;
  Whalf.P = 0.295456f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.41586f;
  Whalf.v = 0.942546f;
  Whalf.P = 0.310406f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.53478f;
  Whalf.v = 0.760037f;
  Whalf.P = 0.285989f;
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
  Whalf.rho = 0.013932f;
  Whalf.v = 0.0f;
  Whalf.P = 7.76405e-5f;
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
  Whalf.rho = 0.919498f;
  Whalf.v = 3.37884f;
  Whalf.P = 869.464f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.941258f;
  Whalf.v = 2.19945f;
  Whalf.P = 922.454f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.902032f;
  Whalf.v = 4.49417f;
  Whalf.P = 813.662f;
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
  Whalf.rho = 0.857525f;
  Whalf.v = -1.93434f;
  Whalf.P = 77.4007f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 0.880649f;
  Whalf.v = -1.45215f;
  Whalf.P = 84.4119f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 0.843058f;
  Whalf.v = -2.31417f;
  Whalf.P = 71.0747f;
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
  Whalf.rho = 5.99924f;
  Whalf.v = 19.5975f;
  Whalf.P = 460.894f;
#elif defined(HYDRO_GAMMA_4_3)
  Whalf.rho = 5.99924f;
  Whalf.v = 19.5975f;
  Whalf.P = 460.894f;
#elif defined(HYDRO_GAMMA_2_1)
  Whalf.rho = 5.99924f;
  Whalf.v = 19.5975f;
  Whalf.P = 460.894f;
#else
#error "Unsupported adiabatic index!"
#endif
  check_riemann_solution(&WL, &WR, &Whalf, "Test 5");
}

/**
 * @brief Check the symmetry of the TRRS Riemann solver
 */
void check_riemann_symmetry(void) {
  float WL[5], WR[5], Whalf1[5], Whalf2[5], n_unit1[3], n_unit2[3], n_norm,
      vij[3], totflux1[5], totflux2[5];

  WL[0] = random_uniform(0.1f, 1.0f);
  WL[1] = random_uniform(-10.0f, 10.0f);
  WL[2] = random_uniform(-10.0f, 10.0f);
  WL[3] = random_uniform(-10.0f, 10.0f);
  WL[4] = random_uniform(0.1f, 1.0f);
  WR[0] = random_uniform(0.1f, 1.0f);
  WR[1] = random_uniform(-10.0f, 10.0f);
  WR[2] = random_uniform(-10.0f, 10.0f);
  WR[3] = random_uniform(-10.0f, 10.0f);
  WR[4] = random_uniform(0.1f, 1.0f);

  n_unit1[0] = random_uniform(-1.0f, 1.0f);
  n_unit1[1] = random_uniform(-1.0f, 1.0f);
  n_unit1[2] = random_uniform(-1.0f, 1.0f);

  n_norm = sqrtf(n_unit1[0] * n_unit1[0] + n_unit1[1] * n_unit1[1] +
                 n_unit1[2] * n_unit1[2]);
  n_unit1[0] /= n_norm;
  n_unit1[1] /= n_norm;
  n_unit1[2] /= n_norm;

  n_unit2[0] = -n_unit1[0];
  n_unit2[1] = -n_unit1[1];
  n_unit2[2] = -n_unit1[2];

  riemann_solver_solve(WL, WR, Whalf1, n_unit1);
  riemann_solver_solve(WR, WL, Whalf2, n_unit2);

  if (!equal(Whalf1[0], Whalf2[0]) || !equal(Whalf1[1], Whalf2[1]) ||
      !equal(Whalf1[2], Whalf2[2]) || !equal(Whalf1[3], Whalf2[3]) ||
      !equal(Whalf1[4], Whalf2[4])) {
    message(
        "Solver asymmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == "
        "[%.3e,%.3e,%.3e,%.3e,%.3e]\n",
        Whalf1[0], Whalf1[1], Whalf1[2], Whalf1[3], Whalf1[4], Whalf2[0],
        Whalf2[1], Whalf2[2], Whalf2[3], Whalf2[4]);
    error("Asymmetry in solution!");
  } else {
    /* message( */
    /*     "Solver symmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == " */
    /*     "[%.3e,%.3e,%.3e,%.3e,%.3e]\n", */
    /*     Whalf1[0], Whalf1[1], Whalf1[2], Whalf1[3], Whalf1[4], Whalf2[0], */
    /*     Whalf2[1], Whalf2[2], Whalf2[3], Whalf2[4]); */
  }

  vij[0] = random_uniform(-10.0f, 10.0f);
  vij[1] = random_uniform(-10.0f, 10.0f);
  vij[2] = random_uniform(-10.0f, 10.0f);

  riemann_solve_for_flux(WL, WR, n_unit1, vij, totflux1);
  riemann_solve_for_flux(WR, WL, n_unit2, vij, totflux2);

  if (!opposite(totflux1[0], totflux2[0]) ||
      !opposite(totflux1[1], totflux2[1]) ||
      !opposite(totflux1[2], totflux2[2]) ||
      !opposite(totflux1[3], totflux2[3]) ||
      !opposite(totflux1[4], totflux2[4])) {
    message(
        "Solver asymmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == "
        "[%.3e,%.3e,%.3e,%.3e,%.3e]\n",
        totflux1[0], totflux1[1], totflux1[2], totflux1[3], totflux1[4],
        totflux2[0], totflux2[1], totflux2[2], totflux2[3], totflux2[4]);
    error("Asymmetry in solution!");
  } else {
    /* message( */
    /*     "Solver symmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == " */
    /*     "[%.3e,%.3e,%.3e,%.3e,%.3e]\n", */
    /*     totflux1[0], totflux1[1], totflux1[2], totflux1[3], totflux1[4], */
    /*     totflux2[0], totflux2[1], totflux2[2], totflux2[3], totflux2[4]); */
  }
}

/**
 * @brief Check the TRRS Riemann solver
 */
int main(int argc, char* argv[]) {

  /* check the TRRS Riemann solver */
  check_riemann_trrs();

  /* symmetry test */
  int i;
  for (i = 0; i < 100; i++) check_riemann_symmetry();

  return 0;
}
