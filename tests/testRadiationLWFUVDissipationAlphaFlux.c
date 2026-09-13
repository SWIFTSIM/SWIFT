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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/* The Stage-3 anisotropic flux-dissipation switch
 * (#radiation_update_dissipation_alpha_flux_band, radiation_isrf.c) is
 * file-static, so it can only be reached through the real dispatch,
 * #radiation_end_gradient_propagation. That function only touches it inside
 * an `#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX` block, and the
 * macro also gates extra `feedback_part_data` fields (radiation_isrf.c's own
 * doxygen), so this test needs the WHOLE project -- library and this test
 * binary alike -- built with the macro defined, or `struct part`'s layout
 * disagrees between this test's compilation unit and the linked
 * libswiftsim.a. It is off by default (a commented-out `#define` in
 * radiation_dissipation_stages.h), so a default build's run of this test is
 * a no-op by construction: see this project's dev log for how it was
 * actually exercised. */
#if defined(FEEDBACK_GEAR) && defined(RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX)

/**
 * @brief Assemble the minimal engine #radiation_end_gradient_propagation
 * reads.
 *
 * @param e (return) The engine.
 * @param cosmo The cosmology it points at.
 * @param fp The feedback properties it points at.
 * @param pc The physical constants it points at.
 */
static void make_engine(struct engine *e, struct cosmology *cosmo,
                        struct feedback_props *fp, struct phys_const *pc) {

  bzero(cosmo, sizeof(struct cosmology));
  bzero(fp, sizeof(struct feedback_props));
  bzero(pc, sizeof(struct phys_const));
  bzero(e, sizeof(struct engine));

  /* Non-cosmological: H = 0 removes the redshift term from the relaxation
   * this ghost also performs, so only the switch under test can move
   * dissipation_alpha_flux_FUV/LW away from their initial value. */
  cosmo->a = 1.0;
  cosmo->H = 0.0;

  fp->LW_FUV_propagation = 1;
  /* Read by the co-located Stage-1 trigger/floor update, which this ghost
   * always runs alongside the Stage-3 switch; sane non-zero values keep
   * those divisions defined without affecting the field under test. */
  fp->LW_FUV_dissipation_alpha_max = 1.f;
  fp->LW_FUV_dissipation_negativity_threshold = 0.1f;
  fp->LW_FUV_dissipation_alpha_floor = 0.f;
  fp->LW_FUV_dissipation_floor_h_over_lambda = 0.5f;

  pc->const_speed_light_c = 1.e4;

  e->cosmology = cosmo;
  e->feedback_props = fp;
  e->physical_constants = pc;
}

/**
 * @brief Fail unless two values agree to a relative tolerance.
 *
 * @param name Name of the quantity, for the failure message.
 * @param expected The analytic value.
 * @param value The value the code produced.
 * @param tol Relative tolerance.
 */
static void check_close(const char *name, double expected, double value,
                        double tol) {

  const double scale = fmax(fabs(expected), 1e-30);
  if (fabs(value - expected) / scale > tol)
    error("%s: expected %.9e, got %.9e", name, expected, value);
}

int main(int argc, char *argv[]) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc);

  struct part p;
  bzero(&p, sizeof(struct part));
  p.h = 1.f;

  struct feedback_part_data *fd = &p.feedback_data;
  fd->dt_prev = 0.1f;
  fd->c_hyp = 2.0f;
  fd->kappa_FUV = 0.f;
  fd->kappa_LW = 0.f;
  /* Deliberately far from 1: this is the value the pre-fix bug used as
   * `u_V = rho_prev*u` instead of the specific `u`. If the division ever
   * regresses to reading `rho_prev` again, the checks below (computed
   * against the specific-`u` formula) catch it, since rho_prev != 1 makes
   * the two formulas disagree by a factor of 4. */
  fd->rho_prev = 4.0f;
  fd->ngb_mean_abs_u_V_FUV = 1.f;
  fd->ngb_mean_abs_u_V_LW = 1.f;

  /* This step's already-finalized div(F) (normally the density loop's own
   * accumulator) and last step's snapshot of it, chosen so `d(div F)/dt` is
   * small and negative (compression): the switch's `if (div_F < 0.f)`
   * branch is the one under test. */
  /* div_specific_flux_*_prev = 0 (rather than some other previous-step
   * value close to this step's) so `d(div F)/dt` is a single division with
   * no subtractive cancellation, keeping the float result close enough to
   * the double reference below for a tight tolerance. */
  fd->u_FUV = 0.5f;
  fd->div_specific_flux_FUV = -6e-4f;
  fd->div_specific_flux_FUV_prev = 0.f;
  fd->dissipation_alpha_flux_FUV = 0.f;

  fd->u_LW = 0.25f;
  fd->div_specific_flux_LW = -2e-4f;
  fd->div_specific_flux_LW_prev = 0.f;
  fd->dissipation_alpha_flux_LW = 0.f;

  radiation_end_gradient_propagation(&p, &e);

  /* Hand-computed reference, in double, from the corrected formula
   * (radiation_isrf.c's own doxygen): shock_estimate =
   * -AMPLITUDE*h_phys^2*d(div F)/dt / (u*c_hyp^2), clamped to [0, 1]; since
   * alpha_prev = 0 <= alpha_aim here, the switch returns alpha_aim directly
   * (no decay branch). */
  const double amplitude = RADIATION_LW_FUV_DISSIPATION_FLUX_SWITCH_AMPLITUDE;
  const double h_phys = (double)cosmo.a * (double)p.h;
  const double c_hyp = (double)fd->c_hyp;
  const double dt = (double)fd->dt_prev;

  const double div_F_rate_FUV = (double)-6e-4 / dt;
  const double shock_estimate_FUV = -amplitude * h_phys * h_phys *
                                    div_F_rate_FUV / ((double)0.5 * c_hyp * c_hyp);
  const double expected_alpha_FUV = fmax(fmin(shock_estimate_FUV, 1.0), 0.0);

  const double div_F_rate_LW = (double)-2e-4 / dt;
  const double shock_estimate_LW = -amplitude * h_phys * h_phys *
                                   div_F_rate_LW / ((double)0.25 * c_hyp * c_hyp);
  const double expected_alpha_LW = fmax(fmin(shock_estimate_LW, 1.0), 0.0);

  message("expected alpha_flux_FUV = %.9e (shock_estimate = %.9e)",
          expected_alpha_FUV, shock_estimate_FUV);
  message("expected alpha_flux_LW = %.9e (shock_estimate = %.9e)",
          expected_alpha_LW, shock_estimate_LW);
  message("code alpha_flux_FUV = %.9e, alpha_flux_LW = %.9e",
          (double)fd->dissipation_alpha_flux_FUV,
          (double)fd->dissipation_alpha_flux_LW);

  check_close("dissipation_alpha_flux_FUV", expected_alpha_FUV,
              (double)fd->dissipation_alpha_flux_FUV, 1e-5);
  check_close("dissipation_alpha_flux_LW", expected_alpha_LW,
              (double)fd->dissipation_alpha_flux_LW, 1e-5);

  /* Neither expected value is clamped to the 0 or 1 boundary: both checks
   * above are exercising the actual division, not a saturated end of it. */
  if (expected_alpha_FUV <= 0.0 || expected_alpha_FUV >= 1.0)
    error("expected_alpha_FUV = %.9e is clamped: the check is vacuous",
          expected_alpha_FUV);
  if (expected_alpha_LW <= 0.0 || expected_alpha_LW >= 1.0)
    error("expected_alpha_LW = %.9e is clamped: the check is vacuous",
          expected_alpha_LW);

  /* Discrimination check: what the pre-fix `u_V = rho_prev*u` formula would
   * have produced for the same inputs. rho_prev = 4 makes this differ from
   * the correct value by exactly that factor, so a silent regression back
   * to the volumetric denominator would land here instead of on the
   * checks above. */
  const double rho_prev = (double)fd->rho_prev;
  const double buggy_alpha_FUV =
      fmax(fmin(shock_estimate_FUV / rho_prev, 1.0), 0.0);
  if (fabs((double)fd->dissipation_alpha_flux_FUV - buggy_alpha_FUV) <
      0.1 * expected_alpha_FUV)
    error(
        "dissipation_alpha_flux_FUV = %.9e matches the pre-fix u_V=rho*u "
        "formula's %.9e rather than the specific-u formula's %.9e",
        (double)fd->dissipation_alpha_flux_FUV, buggy_alpha_FUV,
        expected_alpha_FUV);

  message(
      "Stage-3 flux-dissipation switch matches the specific-u formula, not "
      "the pre-fix volumetric one.");

  return 0;
}

#else

int main(int argc, char *argv[]) {
  message(
      "Skipping: needs FEEDBACK_GEAR and Stage 3 anisotropic-flux "
      "dissipation (RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX, off by "
      "default) enabled.");
  return 0;
}

#endif /* FEEDBACK_GEAR && RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX */
