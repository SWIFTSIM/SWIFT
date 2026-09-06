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
#include "feedback/GEAR/radiation_iact.h"
#include "swift.h"

/* Electron volt in erg (CODATA), hardcoded here for the same reason
 * testRadiationPressureFormula.c hardcodes the speed of light:
 * physical_constants_cgs.h's own definitions are meant for a single
 * translation unit (phys_const_init) and cannot be included a second time
 * without a duplicate-symbol link error. */
static const double test_electron_volt_cgs = 1.602176634e-12;

/* A non-trivial internal unit system (1 Msun, 1 pc, 1 Myr): every
 * conversion-factor call below is then actually exercised, instead of
 * collapsing to 1 the way it would under units_init_cgs(). */
static void make_test_units(struct unit_system *us) {
  units_init(us, /*U_M_in_cgs=*/1.98892e33, /*U_L_in_cgs=*/3.08567758e18,
             /*U_t_in_cgs=*/3.15576e13, /*U_C_in_cgs=*/1.0,
             /*U_T_in_cgs=*/1.0);
}

static void assert_close(const char *name, double actual, double expected,
                         double rel_tol) {
  if (expected == 0.0) {
    if (actual != 0.0) error("%s: expected exactly 0, got %.8e.", name, actual);
    return;
  }
  const double rel_err = fabs(actual - expected) / fabs(expected);
  if (rel_err > rel_tol)
    error("%s: got %.8e, expected %.8e (rel_err=%.3e, tol=%.3e).", name, actual,
          expected, rel_err, rel_tol);
}

static void make_default_cooling(struct cooling_function_data *cooling) {
  bzero(cooling, sizeof(struct cooling_function_data));
  cooling->chemistry_data.local_dust_to_gas_ratio =
      RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO;
}

/* ---------------------------------------------------------------------
 * Receiver-side dust extinction (radiation_get_part_LW_FUV_extinction_
 * factors): band-specific exp(-kappa_eff*Sigma_gas) formula, exercised
 * directly against a hand-computation, not just "returns something in
 * (0, 1]".
 * ------------------------------------------------------------------- */

static double expected_extinction(const struct unit_system *us, float Z,
                                  double sigma_d_band_cgs, double Sigma_gas_p,
                                  double local_dust_to_gas_ratio) {
  const double D_relative =
      fmax((double)Z, 0.0) / RADIATION_GRACKLE_SOLAR_METAL_FRACTION *
      (local_dust_to_gas_ratio / RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);
  const double kappa_eff_cgs =
      sigma_d_band_cgs * D_relative /
      ((double)RADIATION_MU_H * (double)RADIATION_HYDROGEN_MASS_CGS);
  const double kappa_eff = kappa_eff_cgs *
                           units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
                           units_cgs_conversion_factor(us, UNIT_CONV_AREA);
  return exp(-kappa_eff * Sigma_gas_p);
}

static void check_extinction(const char *name, const struct unit_system *us,
                             float h, float rho, float Z) {
  struct part p;
  bzero(&p, sizeof(struct part));
  p.h = h;
  p.rho = rho;

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a = 1.0;
  cosmo.a2_inv = 1.0;

  struct cooling_function_data cooling;
  make_default_cooling(&cooling);

  const float Sigma_gas_c =
      radiation_get_comoving_gas_column_density_at_part(&p);
  const float Sigma_gas_p = Sigma_gas_c * (float)cosmo.a2_inv;

  const double expected_FUV =
      expected_extinction(us, Z, RADIATION_SIGMA_D_FUV_CGS, Sigma_gas_p,
                          RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);
  const double expected_LW =
      expected_extinction(us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p,
                          RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);

  float actual_FUV, actual_LW;
  radiation_get_part_LW_FUV_extinction_factors(us, &cosmo, &p, Z, &cooling,
                                               &actual_FUV, &actual_LW);

  char buf[128];
  snprintf(buf, sizeof(buf), "%s: FUV extinction", name);
  assert_close(buf, (double)actual_FUV, expected_FUV, 1e-4);
  snprintf(buf, sizeof(buf), "%s: LW extinction", name);
  assert_close(buf, (double)actual_LW, expected_LW, 1e-4);

  if (actual_FUV <= 0.0f || actual_FUV > 1.0f)
    error("%s: FUV extinction factor %.6e out of (0, 1].", name,
          (double)actual_FUV);
  if (actual_LW <= 0.0f || actual_LW > 1.0f)
    error("%s: LW extinction factor %.6e out of (0, 1].", name,
          (double)actual_LW);

  /* Zero metallicity must give exactly no extinction (D(Z)=0). */
  if (Z == 0.0f) {
    if (actual_FUV != 1.0f || actual_LW != 1.0f)
      error(
          "%s: zero metallicity did not give exactly 1.0 extinction "
          "(FUV=%.8e, LW=%.8e).",
          name, (double)actual_FUV, (double)actual_LW);
  }
}

/* ---------------------------------------------------------------------
 * Kernel-weighted injection (radiation_iact_nonsym_feedback_apply):
 * superposition of two simultaneously-illuminating stars in the same
 * step, and the LW_FUV_last_touch_ti stamp resetting the field to
 * instantaneous strength (not an ever-growing dose) on the next step any
 * star touches it.
 * ------------------------------------------------------------------- */

static void check_injection(const struct unit_system *us) {
  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a = 1.0;
  cosmo.a2_inv = 1.0;

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));

  const float hi = 1.0f;
  const float r = 0.3f;
  const float dx[3] = {0.3f, 0.0f, 0.0f};
  const float r2 = r * r;
  const float rho_star = 3.0f; /* si->feedback_data.enrichment_weight */
  const float mj = 2.0f;       /* gas particle mass */
  const float rho_gas = 5.0f;  /* gas particle density */
  const float Z_gas = 0.01f;
  const double time_base = 1.0;
  const timebin_t time_bin = 1; /* get_integer_timestep(1) == 4 */

  /* Hand-compute the SPH kernel weight the same way the production code
   * does (real kernel machinery, not the formula under test): */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv);
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;
  const double weight = (double)mj * (double)wi * (1.0 / (double)rho_star);
  const double Delta_t = get_timestep(time_bin, time_base);

  struct part pj;
  bzero(&pj, sizeof(struct part));
  pj.h = hi;
  pj.mass = mj;
  pj.rho = rho_gas;
  pj.chemistry_data
      .smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1] = Z_gas;
  pj.feedback_data.u_FUV = 0.f;
  pj.feedback_data.u_LW = 0.f;
  pj.feedback_data.LW_FUV_last_touch_ti = -1; /* never touched yet */

  struct cooling_function_data cooling;
  make_default_cooling(&cooling);

  float extinction_FUV, extinction_LW;
  radiation_get_part_LW_FUV_extinction_factors(us, &cosmo, &pj, Z_gas, &cooling,
                                               &extinction_FUV, &extinction_LW);

  struct xpart xpj;
  bzero(&xpj, sizeof(struct xpart));

  struct spart si;
  bzero(&si, sizeof(struct spart));
  si.time_bin = time_bin;
  si.feedback_data.enrichment_weight = rho_star;
  si.feedback_data.radiation.L_FUV = 1.0e5;
  si.feedback_data.radiation.L_LW = 5.0e4;

  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &si, &pj, &xpj, &cosmo, /*hydro_props=*/NULL,
      /*fb_props=*/NULL, &phys_const, us, &cooling, /*ti_current=*/0, time_base,
      /*with_cosmology=*/0);

  const double expected_u_FUV_1 = Delta_t * weight *
                                  si.feedback_data.radiation.L_FUV *
                                  extinction_FUV / (double)mj;
  const double expected_u_LW_1 = Delta_t * weight *
                                 si.feedback_data.radiation.L_LW *
                                 extinction_LW / (double)mj;

  /* First touch this step (LW_FUV_last_touch_ti went from -1 to 0): result
   * must equal this star's own deposit exactly. */
  assert_close("injection: first star, u_FUV", (double)pj.feedback_data.u_FUV,
               expected_u_FUV_1, 1e-4);
  assert_close("injection: first star, u_LW", (double)pj.feedback_data.u_LW,
               expected_u_LW_1, 1e-4);

  if (!pj.feedback_data.is_illuminated_LW_FUV)
    error("injection: is_illuminated_LW_FUV not set after first touch.");
  if (!pj.limiter_data.to_be_synchronized)
    error("injection: timestep_sync_part not triggered on first touch.");

  /* A second, independently-illuminating star, same step (ti_current=0
   * again): must ADD to the same particle, not overwrite it -- multiple
   * simultaneously-illuminating sources must superpose within one step. */
  pj.limiter_data.to_be_synchronized = 0;

  struct spart si2;
  bzero(&si2, sizeof(struct spart));
  si2.time_bin = time_bin;
  si2.feedback_data.enrichment_weight = rho_star;
  si2.feedback_data.radiation.L_FUV = 2.0e5;
  si2.feedback_data.radiation.L_LW = 1.0e5;

  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &si2, &pj, &xpj, &cosmo, /*hydro_props=*/NULL,
      /*fb_props=*/NULL, &phys_const, us, &cooling, /*ti_current=*/0, time_base,
      /*with_cosmology=*/0);

  const double expected_u_FUV_2 = Delta_t * weight *
                                  si2.feedback_data.radiation.L_FUV *
                                  extinction_FUV / (double)mj;
  const double expected_u_LW_2 = Delta_t * weight *
                                 si2.feedback_data.radiation.L_LW *
                                 extinction_LW / (double)mj;

  assert_close("injection: two stars, u_FUV", (double)pj.feedback_data.u_FUV,
               expected_u_FUV_1 + expected_u_FUV_2, 1e-4);
  assert_close("injection: two stars, u_LW", (double)pj.feedback_data.u_LW,
               expected_u_LW_1 + expected_u_LW_2, 1e-4);

  /* Already illuminated: no repeated sync on the second touch. */
  if (pj.limiter_data.to_be_synchronized)
    error(
        "injection: timestep_sync_part re-triggered on an already-"
        "illuminated particle.");

  /* A third star, a NEW step (ti_current=1): must reset to just this
   * star's own deposit, not add onto the previous step's two-star total --
   * an instantaneous field strength, not an ever-growing dose. */
  struct spart si3;
  bzero(&si3, sizeof(struct spart));
  si3.time_bin = time_bin;
  si3.feedback_data.enrichment_weight = rho_star;
  si3.feedback_data.radiation.L_FUV = 4.0e5;
  si3.feedback_data.radiation.L_LW = 3.0e4;

  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &si3, &pj, &xpj, &cosmo, /*hydro_props=*/NULL,
      /*fb_props=*/NULL, &phys_const, us, &cooling, /*ti_current=*/1, time_base,
      /*with_cosmology=*/0);

  const double expected_u_FUV_3 = Delta_t * weight *
                                  si3.feedback_data.radiation.L_FUV *
                                  extinction_FUV / (double)mj;
  const double expected_u_LW_3 = Delta_t * weight *
                                 si3.feedback_data.radiation.L_LW *
                                 extinction_LW / (double)mj;

  assert_close("injection: new step resets, u_FUV",
               (double)pj.feedback_data.u_FUV, expected_u_FUV_3, 1e-4);
  assert_close("injection: new step resets, u_LW",
               (double)pj.feedback_data.u_LW, expected_u_LW_3, 1e-4);

  message(
      "injection OK: extinction_FUV=%.6f extinction_LW=%.6f "
      "u_FUV=%.6e u_LW=%.6e (after new-step reset, third star only)",
      (double)extinction_FUV, (double)extinction_LW,
      (double)pj.feedback_data.u_FUV, (double)pj.feedback_data.u_LW);
}

/* ---------------------------------------------------------------------
 * Grackle coupling: isrf_habing (photoelectric heating field) and the H2
 * Lyman-Werner dissociation rate, both a formula-identity check against a
 * hand computation, a zero-field check, and a linearity-in-flux check
 * (doubling u_FUV/u_LW must double the resulting rate -- both formulas
 * are a plain photon-flux * cross-section/normalization product), plus
 * the cooling-side gate wrappers.
 * ------------------------------------------------------------------- */

static void check_grackle_coupling(const struct unit_system *us) {
  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a3_inv = 1.0;

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  const double c_cgs = 2.99792458e10;
  phys_const.const_speed_light_c =
      c_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
  phys_const.const_electron_volt =
      test_electron_volt_cgs /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* A typical diffuse-ISM point: n_H ~ 1 cm^-3, and u_FUV = u_LW chosen so
   * their sum gives G0 ~ 1 in Habing units (c*rho*u_sum ~
   * RADIATION_HABING_FLUX_CGS). */
  const double rho_cgs = 1.6726219e-24;
  const double u_band_cgs = 1.6e10;

  struct part p;
  bzero(&p, sizeof(struct part));
  p.rho = (float)(rho_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  p.feedback_data.u_FUV =
      (float)(u_band_cgs /
              units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  p.feedback_data.u_LW = p.feedback_data.u_FUV;

  /* --- isrf_habing: formula-identity check --- */
  const double u_sum_cgs = 2.0 * u_band_cgs;
  const double flux_cgs = c_cgs * rho_cgs * u_sum_cgs;
  const double expected_G0 = flux_cgs / RADIATION_HABING_FLUX_CGS;

  const double actual_G0 =
      radiation_get_part_isrf_habing(&phys_const, us, &cosmo, &p);
  assert_close("Grackle coupling: isrf_habing", actual_G0, expected_G0, 1e-3);

  /* Order of magnitude: this input was chosen to land near G0 ~ 1; a unit
   * bug anywhere in the internal -> cgs -> Habing chain would show up as
   * many orders of magnitude off, not a small mismatch. */
  if (actual_G0 < 0.1 || actual_G0 > 10.0)
    error(
        "Grackle coupling: isrf_habing = %.6e, expected O(1) for a "
        "typical n_H ~ 1 cm^-3 point.",
        actual_G0);

  /* Zero field -> zero G0, exactly. */
  struct part p_zero;
  bzero(&p_zero, sizeof(struct part));
  p_zero.rho = p.rho;
  const double G0_zero =
      radiation_get_part_isrf_habing(&phys_const, us, &cosmo, &p_zero);
  if (G0_zero != 0.0)
    error(
        "Grackle coupling: isrf_habing at u_FUV=u_LW=0 gave %.6e, "
        "expected exactly 0.",
        G0_zero);

  /* Linearity: doubling the input flux (both bands) must double G0 --
   * G0 = c*rho*(u_FUV+u_LW)/const is a plain linear map of its input. */
  struct part p_double;
  p_double = p;
  p_double.feedback_data.u_FUV *= 2.0f;
  p_double.feedback_data.u_LW *= 2.0f;
  const double G0_double =
      radiation_get_part_isrf_habing(&phys_const, us, &cosmo, &p_double);
  assert_close("Grackle coupling: isrf_habing linearity in u_FUV/u_LW",
               G0_double, 2.0 * actual_G0, 1e-6);

  /* The gate wrapper must pass the same value through unchanged. */
  struct cooling_function_data cooling;
  bzero(&cooling, sizeof(struct cooling_function_data));
  cooling.with_LW_FUV = 1;
  const double gated_G0 =
      cooling_get_isrf_habing_subgrid(&phys_const, us, &cosmo, &cooling, &p);
  assert_close("Grackle coupling: cooling_get_isrf_habing_subgrid gate on",
               gated_G0, actual_G0, 1e-8);
  cooling.with_LW_FUV = 0;
  const double ungated_G0 =
      cooling_get_isrf_habing_subgrid(&phys_const, us, &cosmo, &cooling, &p);
  if (ungated_G0 != 0.0)
    error(
        "Grackle coupling: cooling_get_isrf_habing_subgrid did not gate "
        "off at with_LW_FUV=0 (got %.6e).",
        ungated_G0);

  /* --- LW dissociation rate: formula-identity check --- */
  const double E_LW_photon_cgs =
      RADIATION_LW_PHOTON_ENERGY_EV * test_electron_volt_cgs;
  const double flux_LW_cgs = c_cgs * rho_cgs * u_band_cgs;
  const double expected_k_diss_cgs =
      RADIATION_SIGMA_H2_LW_CGS * (flux_LW_cgs / E_LW_photon_cgs);

  const double actual_k_diss_internal =
      radiation_get_part_LW_dissociation_rate_internal(&phys_const, us, &cosmo,
                                                       &p);
  const double actual_k_diss_cgs =
      actual_k_diss_internal *
      units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);
  assert_close("Grackle coupling: LW dissociation rate", actual_k_diss_cgs,
               expected_k_diss_cgs, 1e-3);

  /* Zero LW field -> zero k_diss, exactly (independent of u_FUV: k_diss
   * only ever reads u_LW). */
  struct part p_zero_LW = p;
  p_zero_LW.feedback_data.u_LW = 0.0f;
  const double k_diss_zero_LW =
      radiation_get_part_LW_dissociation_rate_internal(&phys_const, us, &cosmo,
                                                       &p_zero_LW);
  if (k_diss_zero_LW != 0.0)
    error(
        "Grackle coupling: LW dissociation rate at u_LW=0 gave %.6e, "
        "expected exactly 0.",
        k_diss_zero_LW);

  /* Linearity: doubling u_LW alone must double k_diss. */
  struct part p_double_LW = p;
  p_double_LW.feedback_data.u_LW *= 2.0f;
  const double k_diss_double = radiation_get_part_LW_dissociation_rate_internal(
      &phys_const, us, &cosmo, &p_double_LW);
  assert_close("Grackle coupling: LW dissociation rate linearity in u_LW",
               k_diss_double, 2.0 * actual_k_diss_internal, 1e-6);

  /* Order of magnitude against Draine & Bertoldi's k_LW ~ 1e-10*chi s^-1
   * at a comparable G0 ~ 1 (chi and G0 use slightly different
   * normalizations of the same local FUV/LW field, so this is an
   * order-of-magnitude check, not an exact identity; see
   * theory/GEAR/Radiation/verify_sigma_h2_lw_sternberg2014.py for the
   * fuller cross-check against Sternberg et al. 2014). */
  if (actual_k_diss_cgs < 1e-11 || actual_k_diss_cgs > 1e-9)
    error(
        "Grackle coupling: LW dissociation rate = %.6e s^-1, expected "
        "O(1e-10) at G0 ~ 1.",
        actual_k_diss_cgs);

#if COOLING_GRACKLE_MODE > 1
  const double gated_k_diss = cooling_get_LW_dissociation_rate_subgrid(
      &phys_const, us, &cosmo, &cooling /* with_LW_FUV == 0 here */, &p);
  if (gated_k_diss != 0.0)
    error(
        "Grackle coupling: cooling_get_LW_dissociation_rate_subgrid did "
        "not gate off at with_LW_FUV=0 (got %.6e).",
        gated_k_diss);
  cooling.with_LW_FUV = 1;
  const double gated_k_diss_on = cooling_get_LW_dissociation_rate_subgrid(
      &phys_const, us, &cosmo, &cooling, &p);
  assert_close(
      "Grackle coupling: cooling_get_LW_dissociation_rate_subgrid gate on",
      gated_k_diss_on, actual_k_diss_internal, 1e-8);
#else
  /* H2 untracked at this Grackle mode: the wrapper must return exactly 0
   * regardless of with_LW_FUV. */
  cooling.with_LW_FUV = 1;
  const double gated_k_diss_no_h2 = cooling_get_LW_dissociation_rate_subgrid(
      &phys_const, us, &cosmo, &cooling, &p);
  if (gated_k_diss_no_h2 != 0.0)
    error(
        "Grackle coupling: cooling_get_LW_dissociation_rate_subgrid "
        "returned nonzero (%.6e) at COOLING_GRACKLE_MODE <= 1.",
        gated_k_diss_no_h2);
#endif

  message(
      "Grackle coupling OK: isrf_habing=%.6e (Habing units), "
      "k_diss=%.6e s^-1 (cgs, Draine & Bertoldi comparison point ~1e-10)",
      actual_G0, actual_k_diss_cgs);
}

/* ---------------------------------------------------------------------
 * local_dust_to_gas_ratio scaling: both the injection-side extinction and
 * the propagation-side linear absorption rate must read the resolved
 * chemistry_data.local_dust_to_gas_ratio and scale D(Z) linearly by
 * (that value / RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO), not
 * silently assume Grackle's own compiled default.
 * ------------------------------------------------------------------- */

static void check_local_dust_to_gas_ratio_scaling(
    const struct unit_system *us) {
  struct part p;
  bzero(&p, sizeof(struct part));
  p.h = 1.0f;
  p.rho = 5.0f;
  const float Z = 0.02f;

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a = 1.0;
  cosmo.a2_inv = 1.0;

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(&p) *
      (float)cosmo.a2_inv;

  const double ratios[3] = {RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
                            2.0 * RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
                            0.5 * RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO};
  double kappa_eff_FUV[3];

  for (int i = 0; i < 3; i++) {
    struct cooling_function_data cooling;
    bzero(&cooling, sizeof(struct cooling_function_data));
    cooling.chemistry_data.local_dust_to_gas_ratio = ratios[i];

    /* Injection-side: extinction factor must match the analytic formula
     * scaled by this ratio, not the default. */
    const double expected_FUV = expected_extinction(
        us, Z, RADIATION_SIGMA_D_FUV_CGS, Sigma_gas_p, ratios[i]);
    const double expected_LW = expected_extinction(
        us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p, ratios[i]);
    float actual_FUV, actual_LW;
    radiation_get_part_LW_FUV_extinction_factors(us, &cosmo, &p, Z, &cooling,
                                                 &actual_FUV, &actual_LW);
    char buf[128];
    snprintf(buf, sizeof(buf), "ratio scaling: FUV extinction, ratio=%.6g",
             ratios[i]);
    assert_close(buf, (double)actual_FUV, expected_FUV, 1e-4);
    snprintf(buf, sizeof(buf), "ratio scaling: LW extinction, ratio=%.6g",
             ratios[i]);
    assert_close(buf, (double)actual_LW, expected_LW, 1e-4);

    /* Propagation-side: kappa_eff = -ln(extinction)/Sigma_gas_p must scale
     * linearly in the ratio (recovered independently of the extinction
     * formula, from the public linear-absorption-rate entry point). */
    const float rho_phys = p.rho;
    const float kappa = radiation_get_part_linear_absorption_rate(
        us, Z, rho_phys, RADIATION_SIGMA_D_FUV_CGS, (float)ratios[i]);
    kappa_eff_FUV[i] = (double)kappa / (double)rho_phys;
  }

  /* Doubling/halving the ratio relative to the default must double/halve
   * kappa_eff exactly: this is the linear-scaling claim the fix makes. */
  assert_close("ratio scaling: kappa_eff doubles at 2x ratio", kappa_eff_FUV[1],
               2.0 * kappa_eff_FUV[0], 1e-6);
  assert_close("ratio scaling: kappa_eff halves at 0.5x ratio",
               kappa_eff_FUV[2], 0.5 * kappa_eff_FUV[0], 1e-6);

  message(
      "local_dust_to_gas_ratio scaling OK: kappa_eff(default)=%.6e, "
      "kappa_eff(2x)=%.6e, kappa_eff(0.5x)=%.6e",
      kappa_eff_FUV[0], kappa_eff_FUV[1], kappa_eff_FUV[2]);
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  struct unit_system us;
  make_test_units(&us);

  check_extinction("solar metallicity, moderate column", &us, /*h=*/1.0f,
                   /*rho=*/5.0f, /*Z=*/0.02f);
  check_extinction("low metallicity, high column", &us, /*h=*/2.0f,
                   /*rho=*/50.0f, /*Z=*/0.002f);
  check_extinction("zero metallicity: no extinction", &us, /*h=*/1.0f,
                   /*rho=*/5.0f, /*Z=*/0.0f);

  check_injection(&us);

  check_grackle_coupling(&us);

  check_local_dust_to_gas_ratio_scaling(&us);

  return 0;
}
