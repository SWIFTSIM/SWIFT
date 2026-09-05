/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Compare against the LEBRON formula by hand, in cgs (units_init_cgs
 * makes internal == cgs so the conversion factors in radiation_pressure.c
 * collapse to 1). Formerly a python-side, gas-kinematics-based example check
 * (RadiationPressureMomentumCheck/radiation_pressure_momentum_check.py)
 * covered this, but that measurement is contaminated by the surrounding
 * SPH's own hydrodynamical response to the kick and cannot isolate the
 * formula in isolation; this unit test replaces it as the Tier-1 formula
 * regression check.
 */
static double expected_p_rad(float rho_gas, float h_star, float grad_rho_norm,
                             float Z_star, double L_bol, double Delta_t,
                             double c) {
  const double Z_sun = 0.02;
  const float h_gas = h_star * kernel_gamma;
  const double sobolev_length =
      grad_rho_norm > 0.0f
          ? fmin((double)rho_gas / (double)grad_rho_norm, (double)h_gas)
          : (double)h_gas;
  const double Sigma_gas = ((double)h_gas + sobolev_length) * (double)rho_gas;
  const double tau_IR = 10.0 * ((double)Z_star / Z_sun) * Sigma_gas;
  const double tau_NUV = 1800.0 * ((double)Z_star / Z_sun) * Sigma_gas;
  return Delta_t * L_bol / c * (1.0 - exp(-tau_NUV)) * (1.0 + tau_IR);
}

static void check_case(const char *name, float rho_gas, float h_star,
                       float grad_rho_norm, float Z_star, double L_bol,
                       double Delta_t, double c) {
  struct spart sp;
  bzero(&sp, sizeof(struct spart));
  sp.h = h_star;
  sp.feedback_data.enrichment_weight = rho_gas;
  /* Point the gradient along x only; the formula uses its norm. */
  sp.feedback_data.grad_rho_star[0] = grad_rho_norm;
  sp.feedback_data.grad_rho_star[1] = 0.0f;
  sp.feedback_data.grad_rho_star[2] = 0.0f;
  sp.feedback_data.Z_star = Z_star;
  sp.feedback_data.radiation.L_bol = L_bol;

  struct unit_system us;
  units_init_cgs(&us);

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_speed_light_c = c;

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a2_inv = 1.0;
  cosmo.a = 1.0;

  const double expected =
      expected_p_rad(rho_gas, h_star, grad_rho_norm, Z_star, L_bol, Delta_t, c);
  const float actual = radiation_get_star_physical_radiation_pressure(
      &sp, (float)Delta_t, &phys_const, &us, &cosmo);

  if (expected == 0.0) {
    if (actual != 0.0f)
      error("%s: expected exactly 0, got %.8e.", name, (double)actual);
    return;
  }

  const double rel_err = fabs((double)actual - expected) / fabs(expected);
  if (rel_err > 1e-4)
    error(
        "%s: radiation_get_star_physical_radiation_pressure mismatch -- "
        "got %.8e, expected %.8e (rel_err=%.3e).",
        name, (double)actual, expected, rel_err);
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  const double c_cgs = 2.99792458e10; /* cm/s */
  const double L_bol = 1.0e38;        /* erg/s */
  /* Delta_t * L_bol is computed as a single float32 product inside
   * radiation_get_star_physical_radiation_pressure before dividing by c;
   * keep it well under FLT_MAX (~3.4e38) here even though the final ratio
   * would be representable at a much larger Delta_t. */
  const double Delta_t = 1.0;   /* s */
  const float rho_gas = 1.0f;   /* g/cm^3 */
  const float h_star = 1.0e-4f; /* cm: keeps tau_NUV unsaturated */

  /* Case 1: solar metallicity, no density gradient -- exercises the
   * unsaturated (1-exp(-tau_NUV)) regime this whole factor was added for
   * (the regression this test exists to catch: dropping back to the old
   * "always 1" behaviour would only show up here, not in a saturated case).
   * A zero gradient makes the Sobolev length cap directly at h_gas (not 0),
   * but tau_NUV still lands well short of saturation at this h_star. */
  check_case("unsaturated tau_NUV, no gradient", rho_gas, h_star, 0.0f,
             /*Z_star=*/0.02f, L_bol, Delta_t, c_cgs);

  /* Case 2: same but at a higher column (larger h_star) so tau_NUV >> 1 --
   * the saturated regime, where (1-exp(-tau_NUV)) ~ 1 and the formula
   * degenerates to the pre-fix (1+tau_IR) form. A smaller L_bol here keeps
   * the L_bol*(1+tau_IR) product -ffast-math is free to reassociate before
   * the /c division comfortably under FLT_MAX. */
  check_case("saturated tau_NUV, no gradient", rho_gas, /*h_star=*/1.0e-2f,
             0.0f, /*Z_star=*/0.02f, /*L_bol=*/1.0e35, Delta_t, c_cgs);

  /* Case 3: Z_star == 0 -- must give exactly zero radiation pressure
   * (both opacities scale with Z_star, so tau_IR = tau_NUV = 0). */
  check_case("zero metallicity", rho_gas, h_star, 0.0f, /*Z_star=*/0.0f, L_bol,
             Delta_t, c_cgs);

  /* Case 4: a well-resolved density gradient -- the raw ratio rho_gas/
   * grad_rho_norm is well below h_gas, so the cap is inactive and the
   * Sobolev length reduces to the plain rho_gas/norm_grad_rho asymptote. */
  {
    const float h_gas = h_star * kernel_gamma;
    const float grad_rho_scale = rho_gas / h_gas;
    check_case("resolved density gradient, cap inactive", rho_gas, h_star,
               /*grad_rho_norm=*/10.0f * grad_rho_scale, /*Z_star=*/0.02f,
               L_bol, Delta_t, c_cgs);
  }

  /* Case 5: a small (sub-resolution) density gradient -- the raw ratio
   * rho_gas/grad_rho_norm exceeds h_gas, so the Sobolev length caps at
   * h_gas rather than being divided by the noisy gradient directly. */
  {
    const float h_gas = h_star * kernel_gamma;
    const float grad_rho_scale = rho_gas / h_gas;
    check_case("small density gradient, capped at h_gas", rho_gas, h_star,
               /*grad_rho_norm=*/0.5f * grad_rho_scale, /*Z_star=*/0.02f, L_bol,
               Delta_t, c_cgs);
  }

  /* Case 6: L_bol and tau_IR (~38.7, since a zero gradient here caps the
   * Sobolev length at h_gas rather than dropping to 0) both large enough
   * that L_bol*(1+tau_IR), an intermediate -ffast-math is free to form
   * before the /c division, would overflow float32 (~3.4e38) even though
   * the true result (~1.3e29) does not; regression test for that
   * reassociation-overflow class of bug. */
  check_case("large L_bol and tau_IR: no spurious float32 overflow",
             /*rho_gas=*/1.0f, /*h_star=*/1.0f, 0.0f, /*Z_star=*/0.02f,
             /*L_bol=*/1.0e38, Delta_t, c_cgs);

  /* Case 7: a star with no gas neighbours this step (enrichment_weight and
   * grad_rho_star both exactly 0) -- must give exactly zero radiation
   * pressure. The Sobolev length itself is well-defined (h_gas, from the
   * norm_grad_rho == 0 branch, not a 0/0), but the final rho_gas
   * multiplication zeroes the whole column density out regardless. */
  check_case("no gas neighbours (rho_gas == 0)", /*rho_gas=*/0.0f, h_star, 0.0f,
             /*Z_star=*/0.02f, L_bol, Delta_t, c_cgs);

  /* Case 8: sweep norm_grad_rho across the point where the raw ratio
   * rho_gas/norm_grad_rho crosses h_gas (where the fminf cap engages) and
   * assert p_rad changes smoothly there, with no large jump between
   * adjacent sample points. min() of two continuous functions is itself
   * continuous, so this is a standing regression guard against a future
   * hard cutoff being reintroduced, not a live bug. */
  {
    const float h_gas = h_star * kernel_gamma;
    const float grad_rho_scale = rho_gas / h_gas;
    const float factors[] = {0.1f,  0.5f, 0.9f, 0.99f, 1.0f,
                             1.01f, 1.1f, 2.0f, 10.0f};
    const int n_factors = sizeof(factors) / sizeof(factors[0]);
    double prev_p_rad = -1.0;
    for (int i = 0; i < n_factors; ++i) {
      struct spart sp;
      bzero(&sp, sizeof(struct spart));
      sp.h = h_star;
      sp.feedback_data.enrichment_weight = rho_gas;
      sp.feedback_data.grad_rho_star[0] = factors[i] * grad_rho_scale;
      sp.feedback_data.grad_rho_star[1] = 0.0f;
      sp.feedback_data.grad_rho_star[2] = 0.0f;
      sp.feedback_data.Z_star = 0.02f;
      sp.feedback_data.radiation.L_bol = L_bol;

      struct unit_system us;
      units_init_cgs(&us);
      struct phys_const phys_const;
      bzero(&phys_const, sizeof(struct phys_const));
      phys_const.const_speed_light_c = c_cgs;
      struct cosmology cosmo;
      bzero(&cosmo, sizeof(struct cosmology));
      cosmo.a2_inv = 1.0;
      cosmo.a = 1.0;

      const double p_rad =
          (double)radiation_get_star_physical_radiation_pressure(
              &sp, (float)Delta_t, &phys_const, &us, &cosmo);

      if (i > 0) {
        /* Measured max step across this sweep is ~1.28x (factor=2.0 ->
         * 10.0); 1.5 leaves comfortable margin above that while still far
         * below the ~4x jump the pre-fix hard cutoff produced. */
        const double ratio =
            p_rad > prev_p_rad ? p_rad / prev_p_rad : prev_p_rad / p_rad;
        if (ratio > 1.5)
          error(
              "continuity sweep: p_rad jumped by a factor %.3f between "
              "factor=%.3f and factor=%.3f around the Sobolev-length cap.",
              ratio, factors[i - 1], factors[i]);
      }
      prev_p_rad = p_rad;
    }
  }

  return 0;
}
