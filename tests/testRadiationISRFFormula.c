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
#include "feedback/GEAR/radiation_isrf.h"
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
 * Receiver-side dust extinction (radiation_get_part_ISRF_extinction_
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

  const float Sigma_gas_c = radiation_get_comoving_gas_column_density_at_part(
      &p, 2.0f * (p.h * kernel_gamma));
  const float Sigma_gas_p = Sigma_gas_c * (float)cosmo.a2_inv;

  const double expected_PE =
      expected_extinction(us, Z, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p,
                          RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);
  const double expected_LW =
      expected_extinction(us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p,
                          RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);

  float actual[ISRF_BAND_COUNT];
  radiation_get_part_ISRF_extinction_factors(
      us, &cosmo, &p, Z, &cooling, 2.0f * (p.h * kernel_gamma), actual);
  const float actual_PE = actual[ISRF_BAND_PE];
  const float actual_LW = actual[ISRF_BAND_LW];

  char buf[128];
  snprintf(buf, sizeof(buf), "%s: PE extinction", name);
  assert_close(buf, (double)actual_PE, expected_PE, 1e-4);
  snprintf(buf, sizeof(buf), "%s: LW extinction", name);
  assert_close(buf, (double)actual_LW, expected_LW, 1e-4);

  if (actual_PE <= 0.0f || actual_PE > 1.0f)
    error("%s: PE extinction factor %.6e out of (0, 1].", name,
          (double)actual_PE);
  if (actual_LW <= 0.0f || actual_LW > 1.0f)
    error("%s: LW extinction factor %.6e out of (0, 1].", name,
          (double)actual_LW);

  /* Zero metallicity must give exactly no extinction (D(Z)=0). */
  if (Z == 0.0f) {
    if (actual_PE != 1.0f || actual_LW != 1.0f)
      error(
          "%s: zero metallicity did not give exactly 1.0 extinction "
          "(PE=%.8e, LW=%.8e).",
          name, (double)actual_PE, (double)actual_LW);
  }

  /* Halving the path (R = 1 vs. the R = 2 above)
   * must halve the column, and therefore halve the log-extinction. */
  const float Sigma_gas_c_half =
      radiation_get_comoving_gas_column_density_at_part(
          &p, 1.0f * (p.h * kernel_gamma));
  snprintf(buf, sizeof(buf), "%s: path=1.0 halves the column", name);
  assert_close(buf, (double)Sigma_gas_c_half, 0.5 * (double)Sigma_gas_c, 1e-6);

  float actual_half[ISRF_BAND_COUNT];
  radiation_get_part_ISRF_extinction_factors(
      us, &cosmo, &p, Z, &cooling, 1.0f * (p.h * kernel_gamma), actual_half);
  snprintf(buf, sizeof(buf), "%s: path=1.0 halves the PE log-extinction", name);
  assert_close(buf, log((double)actual_half[ISRF_BAND_PE]),
               0.5 * log((double)actual_PE), 1e-4);
  snprintf(buf, sizeof(buf), "%s: path=1.0 halves the LW log-extinction", name);
  assert_close(buf, log((double)actual_half[ISRF_BAND_LW]),
               0.5 * log((double)actual_LW), 1e-4);
}

/* ---------------------------------------------------------------------
 * Kernel-weighted injection (radiation_iact_nonsym_feedback_apply):
 * superposition of two simultaneously-illuminating stars in the same
 * step, and the ISRF_last_touch_ti stamp resetting the field to
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
  pj.feedback_data.isrf_band[ISRF_BAND_PE].u = 0.f;
  pj.feedback_data.isrf_band[ISRF_BAND_LW].u = 0.f;
  pj.feedback_data.ISRF_last_touch_ti = -1; /* never touched yet */

  struct cooling_function_data cooling;
  make_default_cooling(&cooling);

  float extinction[ISRF_BAND_COUNT];
  radiation_get_part_ISRF_extinction_factors(us, &cosmo, &pj, Z_gas, &cooling,
                                             2.0f * (pj.h * kernel_gamma),
                                             extinction);
  const float extinction_PE = extinction[ISRF_BAND_PE];
  const float extinction_LW = extinction[ISRF_BAND_LW];

  struct xpart xpj;
  bzero(&xpj, sizeof(struct xpart));

  /* ISRF_propagation off: this test exercises the instantaneous-field
   * injection path (the formulas below have no rescale/phi factor). A
   * zero-initialized struct, not NULL: radiation_iact_nonsym_feedback_apply
   * reads fb_props->ISRF_propagation unconditionally. */
  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  fb_props.ISRF_extinction_path_in_kernel_radii = 2.0f;

  struct spart si;
  bzero(&si, sizeof(struct spart));
  si.time_bin = time_bin;
  si.feedback_data.enrichment_weight = rho_star;
  si.feedback_data.radiation.L_band[ISRF_BAND_PE] = 1.0e5;
  si.feedback_data.radiation.L_band[ISRF_BAND_LW] = 5.0e4;

  radiation_iact_nonsym_feedback_apply(r2, dx, hi, /*hj=*/hi, &si, &pj, &xpj,
                                       &cosmo, /*hydro_props=*/NULL,
                                       /*fb_props=*/&fb_props, &phys_const, us,
                                       &cooling, /*ti_current=*/0, time_base,
                                       /*with_cosmology=*/0);

  const double expected_u_PE_1 =
      Delta_t * weight * si.feedback_data.radiation.L_band[ISRF_BAND_PE] *
      extinction_PE / (double)mj;
  const double expected_u_LW_1 =
      Delta_t * weight * si.feedback_data.radiation.L_band[ISRF_BAND_LW] *
      extinction_LW / (double)mj;

  /* First touch this step (ISRF_last_touch_ti went from -1 to 0): result
   * must equal this star's own deposit exactly. */
  assert_close("injection: first star, u_PE",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_PE].u,
               expected_u_PE_1, 1e-4);
  assert_close("injection: first star, u_LW",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_LW].u,
               expected_u_LW_1, 1e-4);

  if (!pj.feedback_data.is_illuminated_ISRF)
    error("injection: is_illuminated_ISRF not set after first touch.");
  if (!pj.limiter_data.to_be_synchronized)
    error("injection: timestep_sync_part not triggered on first touch.");

  /* A second, independently-illuminating star, same step (ti_current=0
   * again): must ADD to the same particle, not overwrite it. Multiple
   * simultaneously-illuminating sources must superpose within one step. */
  pj.limiter_data.to_be_synchronized = 0;

  struct spart si2;
  bzero(&si2, sizeof(struct spart));
  si2.time_bin = time_bin;
  si2.feedback_data.enrichment_weight = rho_star;
  si2.feedback_data.radiation.L_band[ISRF_BAND_PE] = 2.0e5;
  si2.feedback_data.radiation.L_band[ISRF_BAND_LW] = 1.0e5;

  radiation_iact_nonsym_feedback_apply(r2, dx, hi, /*hj=*/hi, &si2, &pj, &xpj,
                                       &cosmo, /*hydro_props=*/NULL,
                                       /*fb_props=*/&fb_props, &phys_const, us,
                                       &cooling, /*ti_current=*/0, time_base,
                                       /*with_cosmology=*/0);

  const double expected_u_PE_2 =
      Delta_t * weight * si2.feedback_data.radiation.L_band[ISRF_BAND_PE] *
      extinction_PE / (double)mj;
  const double expected_u_LW_2 =
      Delta_t * weight * si2.feedback_data.radiation.L_band[ISRF_BAND_LW] *
      extinction_LW / (double)mj;

  assert_close("injection: two stars, u_PE",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_PE].u,
               expected_u_PE_1 + expected_u_PE_2, 1e-4);
  assert_close("injection: two stars, u_LW",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_LW].u,
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
  si3.feedback_data.radiation.L_band[ISRF_BAND_PE] = 4.0e5;
  si3.feedback_data.radiation.L_band[ISRF_BAND_LW] = 3.0e4;

  radiation_iact_nonsym_feedback_apply(r2, dx, hi, /*hj=*/hi, &si3, &pj, &xpj,
                                       &cosmo, /*hydro_props=*/NULL,
                                       /*fb_props=*/&fb_props, &phys_const, us,
                                       &cooling, /*ti_current=*/1, time_base,
                                       /*with_cosmology=*/0);

  const double expected_u_PE_3 =
      Delta_t * weight * si3.feedback_data.radiation.L_band[ISRF_BAND_PE] *
      extinction_PE / (double)mj;
  const double expected_u_LW_3 =
      Delta_t * weight * si3.feedback_data.radiation.L_band[ISRF_BAND_LW] *
      extinction_LW / (double)mj;

  assert_close("injection: new step resets, u_PE",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_PE].u,
               expected_u_PE_3, 1e-4);
  assert_close("injection: new step resets, u_LW",
               (double)pj.feedback_data.isrf_band[ISRF_BAND_LW].u,
               expected_u_LW_3, 1e-4);

  message(
      "injection OK: extinction_PE=%.6f extinction_LW=%.6f "
      "u_PE=%.6e u_LW=%.6e (after new-step reset, third star only)",
      (double)extinction_PE, (double)extinction_LW,
      (double)pj.feedback_data.isrf_band[ISRF_BAND_PE].u,
      (double)pj.feedback_data.isrf_band[ISRF_BAND_LW].u);
}

/* ---------------------------------------------------------------------
 * Dose reservoir: exact accumulation across several stars on different
 * time bins, and the `f = delta/(Delta - k*delta)` drawdown schedule that
 * turns the deposited dose into a constant source rate over the
 * receiving particle's own sub-steps. Mirrors the equivalent cadence
 * check in theory/GEAR/Radiation/verify_isrf_injection_cadence.py. Zero
 * metallicity throughout: extinction is then exactly 1.0 (check_extinction
 * above), isolating the reservoir bookkeeping from the extinction formula.
 * ------------------------------------------------------------------- */

static void check_dose_reservoir(const struct unit_system *us) {
  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a = 1.0;
  cosmo.a2_inv = 1.0;
  cosmo.a3_inv = 1.0;

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));

  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  fb_props.ISRF_propagation = 1;
  fb_props.ISRF_extinction_path_in_kernel_radii = 2.0f;

  const float hi = 1.0f;
  const float r = 0.3f;
  const float dx[3] = {0.3f, 0.0f, 0.0f};
  const float r2 = r * r;
  const float rho_star = 3.0f;
  const float mj = 2.0f;
  const float rho_gas = 5.0f;
  const double time_base = 1.0;
  const timebin_t bin_A = 3; /* get_integer_timestep(3) == 16 */
  const timebin_t bin_B = 1; /* get_integer_timestep(1) == 4 */

  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv);
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;
  const double weight = (double)mj * (double)wi * (1.0 / (double)rho_star);

  struct cooling_function_data cooling;
  make_default_cooling(&cooling);

  struct xpart xpj;
  bzero(&xpj, sizeof(struct xpart));

  /* Two stars, different time bins, same touch step: must both add their
   * own dose, neither resetting nor overwriting the other's. */
  struct part pj;
  bzero(&pj, sizeof(struct part));
  pj.h = hi;
  pj.mass = mj;
  pj.rho = rho_gas;
  pj.feedback_data.ISRF_reservoir_end_ti = -1;

  struct spart siA;
  bzero(&siA, sizeof(struct spart));
  siA.time_bin = bin_A;
  siA.feedback_data.enrichment_weight = rho_star;
  siA.feedback_data.radiation.L_band[ISRF_BAND_PE] = 1.0e5;
  siA.feedback_data.radiation.L_band[ISRF_BAND_LW] = 5.0e4;

  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &siA, &pj, &xpj, &cosmo, /*hydro_props=*/NULL,
      &fb_props, &phys_const, us, &cooling, /*ti_current=*/0, time_base,
      /*with_cosmology=*/0);

  struct spart siB;
  bzero(&siB, sizeof(struct spart));
  siB.time_bin = bin_B;
  siB.feedback_data.enrichment_weight = rho_star;
  siB.feedback_data.radiation.L_band[ISRF_BAND_PE] = 2.0e5;
  siB.feedback_data.radiation.L_band[ISRF_BAND_LW] = 1.0e5;

  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &siB, &pj, &xpj, &cosmo, /*hydro_props=*/NULL,
      &fb_props, &phys_const, us, &cooling, /*ti_current=*/0, time_base,
      /*with_cosmology=*/0);

  const double Delta_A = get_timestep(bin_A, time_base);
  const double Delta_B = get_timestep(bin_B, time_base);
  const double dose_A_PE = Delta_A * weight *
                           siA.feedback_data.radiation.L_band[ISRF_BAND_PE] /
                           (double)mj;
  const double dose_A_LW = Delta_A * weight *
                           siA.feedback_data.radiation.L_band[ISRF_BAND_LW] /
                           (double)mj;
  const double dose_B_PE = Delta_B * weight *
                           siB.feedback_data.radiation.L_band[ISRF_BAND_PE] /
                           (double)mj;
  const double dose_B_LW = Delta_B * weight *
                           siB.feedback_data.radiation.L_band[ISRF_BAND_LW] /
                           (double)mj;

  assert_close(
      "dose reservoir: two stars, different bins, PE",
      (double)pj.feedback_data.isrf_band[ISRF_BAND_PE].u_dose_reservoir,
      dose_A_PE + dose_B_PE, 1e-4);
  assert_close(
      "dose reservoir: two stars, different bins, LW",
      (double)pj.feedback_data.isrf_band[ISRF_BAND_LW].u_dose_reservoir,
      dose_A_LW + dose_B_LW, 1e-4);

  const integertime_t ti_step_A = get_integer_timestep(bin_A);
  const integertime_t ti_step_B = get_integer_timestep(bin_B);
  const integertime_t expected_horizon =
      ti_step_A > ti_step_B ? ti_step_A : ti_step_B;
  if (pj.feedback_data.ISRF_reservoir_end_ti != expected_horizon)
    error(
        "dose reservoir: horizon=%lld, expected max(ti_step_A, ti_step_B)"
        "=%lld.",
        (long long)pj.feedback_data.ISRF_reservoir_end_ti,
        (long long)expected_horizon);

  /* Drawdown schedule: one coarse star's dose, drained by the gas's own
   * finer time bin over N = Delta_A/Delta_B sub-steps, must give a constant
   * rate S = D/Delta_A at every sub-step and drain the reservoir to exactly
   * 0 at the last one. */
  struct part pk;
  bzero(&pk, sizeof(struct part));
  pk.h = hi;
  pk.mass = mj;
  pk.rho = rho_gas;
  pk.feedback_data.ISRF_reservoir_end_ti = -1;
  pk.time_bin = bin_B;

  const integertime_t T = 16; /* a boundary of the coarse star's own bin */
  radiation_iact_nonsym_feedback_apply(
      r2, dx, hi, /*hj=*/hi, &siA, &pk, &xpj, &cosmo, /*hydro_props=*/NULL,
      &fb_props, &phys_const, us, &cooling, /*ti_current=*/T, time_base,
      /*with_cosmology=*/0);

  const double D0_PE =
      (double)pk.feedback_data.isrf_band[ISRF_BAND_PE].u_dose_reservoir;
  const double D0_LW =
      (double)pk.feedback_data.isrf_band[ISRF_BAND_LW].u_dose_reservoir;
  const double dt_gas = get_timestep(bin_B, time_base);
  const double S_PE = D0_PE / Delta_A;
  const double S_LW = D0_LW / Delta_A;
  const int N = (int)lround(Delta_A / dt_gas);

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.policy = 0; /* no cosmology: the plain get_timestep branch */
  e.time_base = time_base;
  e.max_active_bin = num_time_bins; /* every bin active */
  e.internal_units = us;
  e.physical_constants = &phys_const;
  e.cosmology = &cosmo;
  e.cooling_func = &cooling;
  e.feedback_props = &fb_props;

  for (int k = 1; k <= N; k++) {
    e.ti_current = T + (integertime_t)k * get_integer_timestep(bin_B);
    radiation_snapshot_part_propagation(&pk, &e);

    assert_close("dose reservoir: constant drain rate, PE",
                 (double)pk.feedback_data.isrf_band[ISRF_BAND_PE].u_source_rate,
                 S_PE, 1e-5);
    assert_close("dose reservoir: constant drain rate, LW",
                 (double)pk.feedback_data.isrf_band[ISRF_BAND_LW].u_source_rate,
                 S_LW, 1e-5);
  }

  if (fabs((double)pk.feedback_data.isrf_band[ISRF_BAND_PE].u_dose_reservoir) >
          1e-6 * D0_PE ||
      fabs((double)pk.feedback_data.isrf_band[ISRF_BAND_LW].u_dose_reservoir) >
          1e-6 * D0_LW)
    error(
        "dose reservoir: not fully drained at the horizon (PE=%.6e, "
        "LW=%.6e).",
        (double)pk.feedback_data.isrf_band[ISRF_BAND_PE].u_dose_reservoir,
        (double)pk.feedback_data.isrf_band[ISRF_BAND_LW].u_dose_reservoir);

  message(
      "dose reservoir OK: two-star accumulation exact, horizon="
      "max(bins), constant drain rate over N=%d sub-steps, drained to 0.",
      N);
}

/* ---------------------------------------------------------------------
 * Grackle coupling: isrf_habing (photoelectric heating field) and the H2
 * Lyman-Werner dissociation rate, both a formula-identity check against a
 * hand computation, a zero-field check, and a linearity-in-flux check
 * (doubling u_PE/u_LW must double the resulting rate, since both formulas
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

  /* A typical diffuse-ISM point: n_H ~ 1 cm^-3, and u_PE = u_LW chosen so
   * their sum gives G0 ~ 1 in Habing units (c*rho*u_sum ~
   * RADIATION_HABING_FLUX_CGS). */
  const double rho_cgs = 1.6726219e-24;
  const double u_band_cgs = 1.6e10;

  struct part p;
  bzero(&p, sizeof(struct part));
  p.rho = (float)(rho_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  p.feedback_data.isrf_band[ISRF_BAND_PE].u =
      (float)(u_band_cgs /
              units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  p.feedback_data.isrf_band[ISRF_BAND_LW].u =
      p.feedback_data.isrf_band[ISRF_BAND_PE].u;

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
        "Grackle coupling: isrf_habing at u_PE=u_LW=0 gave %.6e, "
        "expected exactly 0.",
        G0_zero);

  /* Linearity: doubling the input flux (both bands) must double G0 --
   * G0 = c*rho*(u_PE+u_LW)/const is a plain linear map of its input. */
  struct part p_double;
  p_double = p;
  p_double.feedback_data.isrf_band[ISRF_BAND_PE].u *= 2.0f;
  p_double.feedback_data.isrf_band[ISRF_BAND_LW].u *= 2.0f;
  const double G0_double =
      radiation_get_part_isrf_habing(&phys_const, us, &cosmo, &p_double);
  assert_close("Grackle coupling: isrf_habing linearity in u_PE/u_LW",
               G0_double, 2.0 * actual_G0, 1e-6);

  /* A band that has undershot below zero must contribute nothing, not be
   * subtracted from the other band: a negative LW energy is propagation
   * undershoot, whose physical contribution is zero illumination. So a
   * mixed-sign particle must give exactly what the positive band alone
   * gives. Clamping only the summed result would instead deliver a
   * silently suppressed field, since the sum itself stays positive and no
   * clamp event fires. */
  struct part p_mix = p;
  p_mix.feedback_data.isrf_band[ISRF_BAND_PE].u *= 4.0f;
  p_mix.feedback_data.isrf_band[ISRF_BAND_LW].u *= -1.0f;
  struct part p_positive_band_only = p;
  p_positive_band_only.feedback_data.isrf_band[ISRF_BAND_PE].u *= 4.0f;
  p_positive_band_only.feedback_data.isrf_band[ISRF_BAND_LW].u = 0.0f;

  const double G0_mix =
      radiation_get_part_isrf_habing(&phys_const, us, &cosmo, &p_mix);
  const double G0_positive_band_only = radiation_get_part_isrf_habing(
      &phys_const, us, &cosmo, &p_positive_band_only);

  /* Bracket rather than a bare lower bound: a non-finite result compares
   * false against every bound and would pass a one-sided gate silently. */
  if (!(G0_mix >= 0.0 && G0_mix < 1e30))
    error(
        "Grackle coupling: isrf_habing with a negative LW band gave a "
        "non-finite or out-of-range value (%.8e).",
        G0_mix);

  assert_close("Grackle coupling: isrf_habing clamps each band before summing",
               G0_mix, G0_positive_band_only, 1e-12);

  /* The same, one level down: k_diss reads the LW band alone, so a
   * negative LW energy must give exactly zero. */
  const double k_diss_negative_LW =
      radiation_get_part_LW_dissociation_rate_internal(&phys_const, us, &cosmo,
                                                       &p_mix);
  if (k_diss_negative_LW != 0.0)
    error(
        "Grackle coupling: LW dissociation rate at a negative u_LW gave "
        "%.8e, expected exactly 0.",
        k_diss_negative_LW);

  /* The gate wrapper must pass the same value through unchanged. */
  struct cooling_function_data cooling;
  bzero(&cooling, sizeof(struct cooling_function_data));
  cooling.with_ISRF = 1;
  const double gated_G0 =
      cooling_get_isrf_habing_subgrid(&phys_const, us, &cosmo, &cooling, &p);
  assert_close("Grackle coupling: cooling_get_isrf_habing_subgrid gate on",
               gated_G0, actual_G0, 1e-8);
  cooling.with_ISRF = 0;
  const double ungated_G0 =
      cooling_get_isrf_habing_subgrid(&phys_const, us, &cosmo, &cooling, &p);
  if (ungated_G0 != 0.0)
    error(
        "Grackle coupling: cooling_get_isrf_habing_subgrid did not gate "
        "off at with_ISRF=0 (got %.6e).",
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

  /* Zero LW field -> zero k_diss, exactly (independent of u_PE: k_diss
   * only ever reads u_LW). */
  struct part p_zero_LW = p;
  p_zero_LW.feedback_data.isrf_band[ISRF_BAND_LW].u = 0.0f;
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
  p_double_LW.feedback_data.isrf_band[ISRF_BAND_LW].u *= 2.0f;
  const double k_diss_double = radiation_get_part_LW_dissociation_rate_internal(
      &phys_const, us, &cosmo, &p_double_LW);
  assert_close("Grackle coupling: LW dissociation rate linearity in u_LW",
               k_diss_double, 2.0 * actual_k_diss_internal, 1e-6);

  /* Order of magnitude against Draine & Bertoldi's k_LW ~ 1e-10*chi s^-1
   * at a comparable G0 ~ 1 (chi and G0 use slightly different
   * normalizations of the same local PE/LW field, so this is an
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
      &phys_const, us, &cosmo, &cooling /* with_ISRF == 0 here */, &p);
  if (gated_k_diss != 0.0)
    error(
        "Grackle coupling: cooling_get_LW_dissociation_rate_subgrid did "
        "not gate off at with_ISRF=0 (got %.6e).",
        gated_k_diss);
  cooling.with_ISRF = 1;
  const double gated_k_diss_on = cooling_get_LW_dissociation_rate_subgrid(
      &phys_const, us, &cosmo, &cooling, &p);
  assert_close(
      "Grackle coupling: cooling_get_LW_dissociation_rate_subgrid gate on",
      gated_k_diss_on, actual_k_diss_internal, 1e-8);
#else
  /* H2 untracked at this Grackle mode: the wrapper must return exactly 0
   * regardless of with_ISRF. */
  cooling.with_ISRF = 1;
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

  const float Sigma_gas_p = radiation_get_comoving_gas_column_density_at_part(
                                &p, 2.0f * (p.h * kernel_gamma)) *
                            (float)cosmo.a2_inv;

  const double ratios[3] = {RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
                            2.0 * RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
                            0.5 * RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO};
  double kappa_eff_PE[3];

  for (int i = 0; i < 3; i++) {
    struct cooling_function_data cooling;
    bzero(&cooling, sizeof(struct cooling_function_data));
    cooling.chemistry_data.local_dust_to_gas_ratio = ratios[i];

    /* Injection-side: extinction factor must match the analytic formula
     * scaled by this ratio, not the default. */
    const double expected_PE = expected_extinction(
        us, Z, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p, ratios[i]);
    const double expected_LW = expected_extinction(
        us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p, ratios[i]);
    float actual[ISRF_BAND_COUNT];
    radiation_get_part_ISRF_extinction_factors(
        us, &cosmo, &p, Z, &cooling, 2.0f * (p.h * kernel_gamma), actual);
    const float actual_PE = actual[ISRF_BAND_PE];
    const float actual_LW = actual[ISRF_BAND_LW];
    char buf[128];
    snprintf(buf, sizeof(buf), "ratio scaling: PE extinction, ratio=%.6g",
             ratios[i]);
    assert_close(buf, (double)actual_PE, expected_PE, 1e-4);
    snprintf(buf, sizeof(buf), "ratio scaling: LW extinction, ratio=%.6g",
             ratios[i]);
    assert_close(buf, (double)actual_LW, expected_LW, 1e-4);

    /* Propagation-side: kappa_eff = -ln(extinction)/Sigma_gas_p must scale
     * linearly in the ratio (recovered independently of the extinction
     * formula, from the public linear-absorption-rate entry point). */
    const float rho_phys = p.rho;
    const float kappa = radiation_get_part_linear_absorption_rate(
        us, Z, rho_phys, RADIATION_SIGMA_D_PE_CGS, (float)ratios[i]);
    kappa_eff_PE[i] = (double)kappa / (double)rho_phys;
  }

  /* Doubling/halving the ratio relative to the default must double/halve
   * kappa_eff exactly: this is the linear-scaling claim the fix makes. */
  assert_close("ratio scaling: kappa_eff doubles at 2x ratio", kappa_eff_PE[1],
               2.0 * kappa_eff_PE[0], 1e-6);
  assert_close("ratio scaling: kappa_eff halves at 0.5x ratio", kappa_eff_PE[2],
               0.5 * kappa_eff_PE[0], 1e-6);

  message(
      "local_dust_to_gas_ratio scaling OK: kappa_eff(default)=%.6e, "
      "kappa_eff(2x)=%.6e, kappa_eff(0.5x)=%.6e",
      kappa_eff_PE[0], kappa_eff_PE[1], kappa_eff_PE[2]);
}

/* ---------------------------------------------------------------------
 * Energy-ledger accumulation (#feedback_isrf_band_data.cumulative_injected/
 * cumulative_absorbed, SWIFT_DEBUG_CHECKS only): #radiation_end_force_
 * propagation's own I/A split, driven for two steps and checked against a
 * hand-integration of the exact-relaxation ODE, not just a re-statement of
 * the code's own expressions.
 * ------------------------------------------------------------------- */

#ifdef SWIFT_DEBUG_CHECKS
static void check_energy_ledger_accumulation(void) {

  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  struct engine e;
  bzero(&cosmo, sizeof(cosmo));
  bzero(&fp, sizeof(fp));
  bzero(&pc, sizeof(pc));
  bzero(&e, sizeof(e));
  cosmo.a = 1.0;
  cosmo.H = 0.0; /* Non-cosmological: isolates the source/transport split. */
  fp.ISRF_propagation = 1;
  pc.const_speed_light_c = 1.0e4;
  e.cosmology = &cosmo;
  e.feedback_props = &fp;
  e.physical_constants = &pc;

  struct part p;
  bzero(&p, sizeof(p));
  struct feedback_part_data *fd = &p.feedback_data;
  const float dt = 0.37f;
  const float c_hyp = 3.2f;
  const float kappa = 0.8f;
  const float u_source_rate = 1.1f;
  const float div_specific_flux = -0.6f;
  const float dissipation_u = 0.4f;
  const double u_prev0 = 2.5;
  fd->dt_prev = dt;
  fd->c_hyp = c_hyp;
  fd->isrf_band[ISRF_BAND_PE].kappa = kappa;
  fd->isrf_band[ISRF_BAND_PE].u_source_rate = u_source_rate;
  fd->isrf_band[ISRF_BAND_PE].div_specific_flux = div_specific_flux;
  fd->isrf_band[ISRF_BAND_PE].dissipation_u = dissipation_u;
  fd->isrf_band[ISRF_BAND_PE].u_prev = (float)u_prev0;

  /* Hand-integration, in double, of the same exact-relaxation ODE
   * (radiation_isrf.c's #radiation_end_force_propagation): `a =
   * (c_hyp*kappa + H)*dt`, `decay = exp(-a)`, `phi = (1-decay)/a` (`a` is
   * not near 0 here, so the direct ratio is well-conditioned and does not
   * need the Taylor branch #radiation_relaxation_phi_factor takes there),
   * `rescale = c_hyp/c`. Two identical-input steps, `u_prev` snapshotted
   * from the previous step's output between them, as the real per-step
   * driver (#radiation_snapshot_part_propagation) does. */
  const double a = ((double)c_hyp * kappa + cosmo.H) * dt;
  const double decay = exp(-a);
  const double phi = (1.0 - decay) / a;
  const double rescale = (double)c_hyp / pc.const_speed_light_c;
  const double I_step = dt * rescale * u_source_rate;
  const double residual_step = dt * (phi * dissipation_u - div_specific_flux);

  double u_prev = u_prev0;
  double Inj_expected = 0.0, Abs_expected = 0.0;
  for (int step = 0; step < 2; step++) {
    const double A_step =
        (u_prev + dt * phi * dissipation_u) * (1.0 - decay) +
        (rescale * u_source_rate - div_specific_flux) * dt * (1.0 - phi);
    const double u_new =
        u_prev + I_step - A_step + residual_step; /* Identity, not restated
                                                       from the code. */
    Inj_expected += I_step;
    Abs_expected += A_step;

    radiation_end_force_propagation(&p, &e);
    assert_close("cumulative_injected (running)",
                 (double)fd->isrf_band[ISRF_BAND_PE].cumulative_injected,
                 Inj_expected, 1e-5);
    assert_close("cumulative_absorbed (running)",
                 (double)fd->isrf_band[ISRF_BAND_PE].cumulative_absorbed,
                 Abs_expected, 1e-5);
    assert_close("u after end_force_propagation",
                 (double)fd->isrf_band[ISRF_BAND_PE].u, u_new, 1e-5);

    u_prev = u_new;
    fd->isrf_band[ISRF_BAND_PE].u_prev = fd->isrf_band[ISRF_BAND_PE].u;
  }

  /* The ledger identity itself, `E + Abs - Inj`, hand-derived (not the
   * code's own arithmetic): must equal the initial field plus the two
   * steps' transport/dissipation residual. */
  const double E = fd->isrf_band[ISRF_BAND_PE].u;
  const double ledger_lhs = E + Abs_expected - Inj_expected;
  const double ledger_rhs = u_prev0 + 2.0 * residual_step;
  assert_close("ledger identity E+Abs-Inj", ledger_lhs, ledger_rhs, 1e-5);

  message("energy ledger OK: Inj=%.6e, Abs=%.6e, E=%.6e, E+Abs-Inj=%.6e",
          Inj_expected, Abs_expected, E, ledger_lhs);
}
#endif

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

  check_dose_reservoir(&us);

  check_grackle_coupling(&us);

  check_local_dust_to_gas_ratio_scaling(&us);

#ifdef SWIFT_DEBUG_CHECKS
  check_energy_ledger_accumulation();
#endif

  return 0;
}
