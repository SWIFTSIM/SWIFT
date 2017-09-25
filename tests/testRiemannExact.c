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

/* Force use of exact Riemann solver */
#undef RIEMANN_SOLVER_TRRS
#undef RIEMANN_SOLVER_HLLC
#undef RIEMANN_SOLVER_EXACT
#define RIEMANN_SOLVER_EXACT 1

/* Some standard headers. */
#include <string.h>

/* Local headers. */
#include "riemann/riemann_exact.h"
#include "swift.h"

const float max_abs_error = 1e-3f;
const float max_rel_error = 1e-2f;
const float min_threshold = 1e-2f;

/**
 * @brief Checks whether two numbers are opposite of each others.
 */
int are_symmetric(float a, float b) {

  /* Check that the signs are different */
  if ((a * b) > 0.f) {
    message("Identical signs a=%.8e b=%.8e", a, b);
    return 0;
  }

  const float abs_a = fabsf(a);
  const float abs_b = fabsf(b);

  const float abs_error = fabsf(abs_a - abs_b);

  /* Avoid FPEs... */
  if (fabsf(abs_a + abs_b) == 0.f) {
    return 1;
  }

  /* Avoid things close to 0 */
  if ((abs_a < min_threshold) || (abs_b < min_threshold)) {
    return 1;
  }

  const float rel_error = 0.5f * abs_error / fabsf(abs_a + abs_b);

  /* Check that we do not breach the relative error limit */
  if (rel_error > max_rel_error) {
    message("Relative error too large a=%.8e b=%.8e rel=%.8e", a, b, rel_error);
    return 0;
  }

  /* All good */
  return 1;
}

int equal(float a, float b) {

  const float abs_error = fabsf(a - b);

  /* Avoid FPEs... */
  if (fabsf(a + b) == 0.f) {
    return 1;
  }

  /* Avoid things close to 0 */
  if ((fabsf(a) < min_threshold) || (fabsf(b) < min_threshold)) {
    return 1;
  }

  const float rel_error = 0.5f * abs_error / fabsf(a + b);

  /* Check that we do not breach the relative error limit */
  if (rel_error > max_rel_error) {
    message("Relative error too large a=%.8e b=%.8e rel=%.8e", a, b, rel_error);
    return 0;
  }

  /* All good */
  return 1;
}

/**
 * @brief Check that a and b are consistent (up to some error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(float a, float b, const char* s) {
  if (fabsf(a + b) != 0.f && fabsf(a - b) / fabsf(a + b) > max_rel_error &&
      fabsf(a - b) > max_abs_error) {
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
void check_riemann_exact(void) {
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
    message("Asymmetry in solution!\n");
    /* This asymmetry is to be expected, since we do an iteration. Are the
       results at least consistent? */
    check_value(Whalf1[0], Whalf2[0], "Rho solution");
    check_value(Whalf1[1], Whalf2[1], "V[0] solution");
    check_value(Whalf1[2], Whalf2[2], "V[1] solution");
    check_value(Whalf1[3], Whalf2[3], "V[2] solution");
    check_value(Whalf1[4], Whalf2[4], "Pressure solution");
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

  if (!are_symmetric(totflux1[0], totflux2[0]) ||
      !are_symmetric(totflux1[1], totflux2[1]) ||
      !are_symmetric(totflux1[2], totflux2[2]) ||
      !are_symmetric(totflux1[3], totflux2[3]) ||
      !are_symmetric(totflux1[4], totflux2[4])) {
    message("WL=[%.8e, %.8e, %.8e, %.8e, %.8e]", WL[0], WL[1], WL[2], WL[3],
            WL[4]);
    message("WR=[%.8e, %.8e, %.8e, %.8e, %.8e]", WR[0], WR[1], WR[2], WR[3],
            WR[4]);
    message("n_unit1=[%.8e, %.8e, %.8e]", n_unit1[0], n_unit1[1], n_unit1[2]);
    message("vij=[%.8e, %.8e, %.8e]\n", vij[0], vij[1], vij[2]);
    message(
        "Flux solver asymmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == "
        "[%.3e,%.3e,%.3e,%.3e,%.3e]\n",
        totflux1[0], totflux1[1], totflux1[2], totflux1[3], totflux1[4],
        totflux2[0], totflux2[1], totflux2[2], totflux2[3], totflux2[4]);
    /* This asymmetry is to be expected, since we do an iteration. Are the
       results at least consistent? */
    check_value(totflux1[0], totflux2[0], "Mass flux");
    check_value(totflux1[1], totflux2[1], "Momentum[0] flux");
    check_value(totflux1[2], totflux2[2], "Momentum[1] flux");
    check_value(totflux1[3], totflux2[3], "Momentum[2] flux");
    check_value(totflux1[4], totflux2[4], "Energy flux");
    error("Asymmetry in flux solution!");
  } else {
    /* message( */
    /*     "Flux solver symmetric: [%.3e,%.3e,%.3e,%.3e,%.3e] == " */
    /*     "[%.3e,%.3e,%.3e,%.3e,%.3e]\n", */
    /*     totflux1[0], totflux1[1], totflux1[2], totflux1[3], totflux1[4], */
    /*     totflux2[0], totflux2[1], totflux2[2], totflux2[3], totflux2[4]); */
  }
}

/**
 * @brief Check the exact Riemann solver
 */
int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* check the exact Riemann solver */
  check_riemann_exact();

  /* symmetry test */
  for (int i = 0; i < 100000; ++i) {
    check_riemann_symmetry();
  }

  return 0;
}
