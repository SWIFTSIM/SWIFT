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

/* The cosmological expansion term of the LW/PE hyperbolic propagation, AND
 * the band-edge transfer it now carries (band-edge transfer derivation,
 * section 2): each moment's relaxation depth is `a(m) = (c_hyp*kappa(m) +
 * lambda(m)*H_dilated)*dt`, `lambda(m) = 1` recovering the pre-transfer grey
 * result; the LW energy moment's own lower-edge loss is additionally
 * transferred into the PE energy moment's `u` (variant A, band-edge
 * transfer derivation section 2.9), which #ISRF_MOMENT_LW_PHOTON does not
 * receive (no PE-side photon-number moment exists to receive it) despite
 * sharing ISRF_OPERATOR_LW's `kappa` with #ISRF_MOMENT_LW.
 *
 * REFERENCE: solved independently here from the two decoupled per-moment
 * ODEs `du(m)/dt = -r(m)*u(m)` (LW, LW_PHOTON: r(m) = c_hyp*kappa(m) +
 * lambda(m)*H_dilated) and `du_PE/dt = -r_PE*u_PE + (lambda_LW - 1) *
 * H_dilated * u_LW(t)`, via the standard linear-ODE integrating factor,
 * NOT by re-deriving the C code's own recurrence: this is "variant A"'s own
 * CLOSED FORM (the transfer added undecayed to PE's own already-relaxed
 * `u`, matching #radiation_end_force_propagation's own doxygen), which for
 * a single step is
 *
 *   u_PE(t) = u_PE(0)*e^{-r_PE t}
 *             + (lambda_LW-1)*H_dilated*u_LW(0)*(1-e^{-r_LW t})/r_LW
 *
 * (the second term is `f_edge * absorbed_LW` worked out in closed form: see
 * radiation_isrf.c's own doxygen for `f_edge`/`absorbed_LW`), with no
 * degenerate limit needed (unlike the TRUE continuum solution, which has a
 * removable singularity at `r_PE = r_LW`: see the comment below for why a
 * degenerate leg is still exercised). Driven through the real ghost
 * functions, with transport, sources and the dissipation accumulators all
 * held at zero, so only this term can move the state.
 */
#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_isrf.h"

static const float u_PE_0 = 3.0f;
static const float u_LW_0 = 0.8f;
static const float u_LW_PHOTON_0 = 1.2f;
static const float F_PE_0[3] = {0.4f, -0.2f, 0.1f};
static const float F_LW_0[3] = {-0.15f, 0.05f, 0.3f};
static const float F_LW_PHOTON_0[3] = {0.2f, 0.2f, -0.1f};

/**
 * @brief Assemble the minimal engine the two ghost functions read.
 *
 * @param e (return) The engine.
 * @param cosmo The cosmology it points at.
 * @param fp The feedback properties it points at.
 * @param pc The physical constants it points at.
 * @param a The scale factor to report.
 * @param H The Hubble rate to report, in internal units of inverse physical
 * time.
 * @param lambda_pe #feedback_props.band_edge_weight_pe to set.
 * @param lambda_lw #feedback_props.band_edge_weight_lw to set.
 * @param lambda_n_lw #feedback_props.band_edge_photon_weight_lw to set.
 */
static void make_engine(struct engine *e, struct cosmology *cosmo,
                        struct feedback_props *fp, struct phys_const *pc,
                        double a, double H, double lambda_pe, double lambda_lw,
                        double lambda_n_lw) {

  bzero(cosmo, sizeof(struct cosmology));
  bzero(fp, sizeof(struct feedback_props));
  bzero(pc, sizeof(struct phys_const));
  bzero(e, sizeof(struct engine));

  cosmo->a = a;
  cosmo->H = H;

  fp->ISRF_propagation = 1;
  fp->ISRF_extinction_path_in_kernel_radii = 2.0f;
  /* Only the dissipation-coefficient update reads these, and it cannot touch
   * `u` or `F`; they are set to sane non-zero values purely so that update's
   * own divisions stay defined. */
  fp->ISRF_dissipation_alpha_max = 1.f;
  fp->ISRF_dissipation_negativity_threshold = 0.1f;
  fp->ISRF_dissipation_alpha_floor = 0.f;
  fp->ISRF_dissipation_floor_h_over_lambda = 0.5f;

  fp->band_edge_weight_pe = lambda_pe;
  fp->band_edge_weight_lw = lambda_lw;
  fp->band_edge_photon_weight_lw = lambda_n_lw;

  pc->const_speed_light_c = 1.e4;

  e->cosmology = cosmo;
  e->feedback_props = fp;
  e->physical_constants = pc;
}

/**
 * @brief Set a particle to the reference radiation state, with no transport,
 * no source and no dissipation.
 *
 * @param p (return) The particle.
 * @param kappa_pe The PE operator's dust absorption rate.
 * @param kappa_lw The LW operator's dust absorption rate (shared by
 * #ISRF_MOMENT_LW and #ISRF_MOMENT_LW_PHOTON).
 * @param dt The step the ghosts should integrate over.
 * @param c_hyp The hyperbolic propagation speed.
 */
static void set_part(struct part *p, float kappa_pe, float kappa_lw, float dt,
                     float c_hyp) {

  bzero(p, sizeof(struct part));
  p->h = 1.f;

  struct feedback_part_data *fd = &p->feedback_data;
  fd->dt_prev = dt;
  fd->c_hyp = c_hyp;
  fd->isrf_operator[ISRF_OPERATOR_PE].kappa = kappa_pe;
  fd->isrf_operator[ISRF_OPERATOR_LW].kappa = kappa_lw;
  fd->rho_prev = 1.f;
  fd->isrf_moment[ISRF_MOMENT_PE].u_prev = u_PE_0;
  fd->isrf_moment[ISRF_MOMENT_LW].u_prev = u_LW_0;
  fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u_prev = u_LW_PHOTON_0;
  fd->isrf_moment[ISRF_MOMENT_PE].u = u_PE_0;
  fd->isrf_moment[ISRF_MOMENT_LW].u = u_LW_0;
  fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u = u_LW_PHOTON_0;
  fd->isrf_operator[ISRF_OPERATOR_PE].ngb_mean_abs_u_V = 1.f;
  fd->isrf_operator[ISRF_OPERATOR_LW].ngb_mean_abs_u_V = 1.f;
  for (int k = 0; k < 3; k++) {
    fd->isrf_moment[ISRF_MOMENT_PE].specific_flux[k] = F_PE_0[k];
    fd->isrf_moment[ISRF_MOMENT_LW].specific_flux[k] = F_LW_0[k];
    fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[k] = F_LW_PHOTON_0[k];
  }
}

/**
 * @brief Fail unless two values agree to a relative tolerance, with an
 * absolute floor on the scale so a near-zero expected value does not demand
 * an unreasonably tight ABSOLUTE match.
 *
 * @param name Name of the quantity, for the failure message.
 * @param expected The analytic value.
 * @param value The value the code produced.
 * @param tol Relative tolerance.
 * @param floor Smallest scale the tolerance is applied to.
 */
static void check_close(const char *name, double expected, double value,
                        double tol, double floor) {

  const double scale = fmax(fabs(expected), floor);
  if (fabs(value - expected) / scale > tol)
    error("%s: expected %.9e, got %.9e (rel diff %.3e)", name, expected, value,
          fabs(value - expected) / scale);
}

/**
 * @brief One (a, H, kappa_pe, kappa_lw, lambda_pe, lambda_lw, lambda_n_lw)
 * case: run both ghosts once and compare `u` and `F` for all three moments
 * against this file's own closed-form reference (see the file header).
 *
 * @param label Name of this leg, for the log message.
 * @param a The scale factor.
 * @param H The Hubble rate.
 * @param kappa_pe The PE operator's dust absorption rate.
 * @param kappa_lw The LW operator's dust absorption rate.
 * @param lambda_pe #feedback_props.band_edge_weight_pe.
 * @param lambda_lw #feedback_props.band_edge_weight_lw.
 * @param lambda_n_lw #feedback_props.band_edge_photon_weight_lw.
 * @param dt The step length.
 * @param c_hyp The hyperbolic propagation speed.
 */
static void run_case(const char *label, double a, double H, float kappa_pe,
                     float kappa_lw, double lambda_pe, double lambda_lw,
                     double lambda_n_lw, float dt, float c_hyp) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, a, H, lambda_pe, lambda_lw, lambda_n_lw);

  struct part p;
  set_part(&p, kappa_pe, kappa_lw, dt, c_hyp);

  radiation_end_gradient_propagation(&p, &e);
  radiation_end_force_propagation(&p, &e);

  /* Independent double-precision reference: H_dilated is the SAME quantity
   * the code forms (c_hyp/c)*H, computed here from the same inputs but a
   * fresh expression, not copied from radiation_isrf.c. */
  const double H_dilated = ((double)c_hyp / (double)pc.const_speed_light_c) * H;
  const double r_pe = (double)c_hyp * (double)kappa_pe + lambda_pe * H_dilated;
  const double r_lw = (double)c_hyp * (double)kappa_lw + lambda_lw * H_dilated;
  const double r_lw_photon =
      (double)c_hyp * (double)kappa_lw + lambda_n_lw * H_dilated;
  const double t = (double)dt;

  const double decay_pe = exp(-r_pe * t);
  const double decay_lw = exp(-r_lw * t);
  const double decay_lw_photon = exp(-r_lw_photon * t);

  /* variant A's own closed form for the transfer (file header): 0 when
   * r_lw == 0 exactly (a_LW <= 0 in the code's own guard), matching
   * radiation_isrf.c's structural guard, not an approximation of it. */
  const double transfer = (r_lw > 0.) ? (lambda_lw - 1.) * H_dilated * u_LW_0 *
                                            (1. - decay_lw) / r_lw
                                      : 0.;

  const double u_pe_expected = u_PE_0 * decay_pe + transfer;
  const double u_lw_expected = u_LW_0 * decay_lw;
  const double u_lw_photon_expected = u_LW_PHOTON_0 * decay_lw_photon;

  /* For CONTEXT ONLY (not a pass/fail assertion): the gap between variant
   * A's closed form above and the TRUE continuum solution (band-edge
   * transfer derivation, section 2.9: exact away from the removable
   * r_pe == r_lw singularity, l'Hopital's rule there). Logged so a reader
   * can see the size of the approximation variant A makes, without the
   * test depending on an approximate bound to catch a genuine bug. */
  const double denom = r_pe - r_lw;
  const double continuum_transfer =
      (r_lw <= 0.)
          ? 0.
          : (fabs(denom) > 1e-6 * fmax(r_pe, r_lw)
                 ? (lambda_lw - 1.) * H_dilated * u_LW_0 *
                       (decay_lw - decay_pe) / denom
                 : (lambda_lw - 1.) * H_dilated * u_LW_0 * t * decay_pe);

  check_close("u_PE", u_pe_expected,
              (double)p.feedback_data.isrf_moment[ISRF_MOMENT_PE].u, 2e-6,
              1e-8);
  check_close("u_LW", u_lw_expected,
              (double)p.feedback_data.isrf_moment[ISRF_MOMENT_LW].u, 2e-6,
              1e-8);
  check_close("u_LW_PHOTON", u_lw_photon_expected,
              (double)p.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u,
              2e-6, 1e-8);
  for (int k = 0; k < 3; k++) {
    check_close("specific_flux_PE", decay_pe * F_PE_0[k],
                p.feedback_data.isrf_moment[ISRF_MOMENT_PE].specific_flux[k],
                2e-6, 1e-8);
    check_close("specific_flux_LW", decay_lw * F_LW_0[k],
                p.feedback_data.isrf_moment[ISRF_MOMENT_LW].specific_flux[k],
                2e-6, 1e-8);
    check_close(
        "specific_flux_LW_PHOTON", decay_lw_photon * F_LW_PHOTON_0[k],
        p.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[k],
        2e-6, 1e-8);
  }

  message(
      "%s: r_PE=%.6g r_LW=%.6g u_PE %.6e (variant-A/continuum transfer "
      "%.3e/%.3e)",
      label, r_pe, r_lw, (double)p.feedback_data.isrf_moment[ISRF_MOMENT_PE].u,
      transfer, continuum_transfer);
}

int main(int argc, char *argv[]) {

  /* Leg 1: ordinary case, non-degenerate r_PE != r_LW, moderate lambda. */
  run_case("moderate", /*a=*/1.0, /*H=*/0.3, /*kappa_pe=*/1.f, /*kappa_lw=*/3.f,
           /*lambda_pe=*/2.154, /*lambda_lw=*/6.508, /*lambda_n_lw=*/6.0,
           /*dt=*/0.05f, /*c_hyp=*/2.f);

  /* Leg 2: LARGE lambda spread between the two bands, so a copy-paste of
   * one band's lambda into the other's slot fails immediately. */
  run_case("large_lambda_spread", /*a=*/1.0, /*H=*/0.5, /*kappa_pe=*/2.f,
           /*kappa_lw=*/0.5f, /*lambda_pe=*/1.5, /*lambda_lw=*/120.0,
           /*lambda_n_lw=*/40.0, /*dt=*/0.02f, /*c_hyp=*/1.5f);

  /* Leg 3: the r_PE == r_LW DEGENERACY (band-edge transfer derivation,
   * section 7 Gate 2: unreachable in production since kappa_LW > kappa_PE
   * while lambda_LW > lambda_PE puts the two sides on opposite signs, so
   * set kappa directly to hit it here). c_hyp=2, c=1e4 => H_dilated =
   * (c_hyp/c)*H = 2e-4*H. Want c_hyp*(kappa_PE-kappa_LW) = (lambda_LW -
   * lambda_PE)*H_dilated: with lambda_PE=2, lambda_LW=6, H=1e4 (H_dilated =
   * 2), the right side is 4*2=8, so kappa_PE-kappa_LW = 8/c_hyp = 4;
   * kappa_LW=1, kappa_PE=5 gives r_PE = 2*5+2*2 = 14 = r_LW = 2*1+6*2, exact
   * to double precision (checked at the values above; asserted below too,
   * not just claimed here). This exercises run_case()'s own l'Hopital
   * branch, the SAME degeneracy the code's continuum reference (not the
   * code itself, which has no degenerate branch: variant A's closed form
   * used for the PASS/FAIL check above has none either) would need it for
   * if it were the assertion; kept here as the leg that would catch a
   * future switch to the continuum-exact treatment mishandling it. */
  run_case("degenerate_r_pe_eq_r_lw", /*a=*/0.5, /*H=*/1.e4, /*kappa_pe=*/5.f,
           /*kappa_lw=*/1.f, /*lambda_pe=*/2.0, /*lambda_lw=*/6.0,
           /*lambda_n_lw=*/5.5, /*dt=*/0.01f, /*c_hyp=*/2.f);

  /* Leg 4: ISRF_MOMENT_LW_PHOTON's own lambda_N(LW) far from lambda_E(LW),
   * sharing ISRF_OPERATOR_LW's kappa with ISRF_MOMENT_LW (band-edge
   * transfer derivation, section 1.5/2.8): a coefficient applied per
   * OPERATOR instead of per MOMENT collapses lambda_N(LW) to lambda_E(LW)
   * and fails u_LW_PHOTON's own check above at every leg, but this one
   * maximises the gap. */
  run_case("lw_photon_own_lambda", /*a=*/1.0, /*H=*/0.8, /*kappa_pe=*/0.5f,
           /*kappa_lw=*/0.5f, /*lambda_pe=*/2.0, /*lambda_lw=*/50.0,
           /*lambda_n_lw=*/3.0, /*dt=*/0.03f, /*c_hyp=*/3.f);

  /* Leg 5: stiff leg, absorption depth large alongside a redshift depth
   * that is comparatively small: this is the fix's own physical prediction
   * (the corrected Hubble term "does almost nothing" once dust absorption
   * dominates). Uses the pre-band-edge-transfer lambda = 1 grey values
   * (this file's own regression floor, kept from the pre-transfer version
   * of this test), so the transfer term is (lambda_LW-1) = 0 exactly and
   * this leg is a pure per-band exponential-decay check, at any a_PE. */
  run_case("stiff", /*a=*/0.5, /*H=*/10.0, /*kappa_pe=*/50.f, /*kappa_lw=*/50.f,
           /*lambda_pe=*/1.0, /*lambda_lw=*/1.0, /*lambda_n_lw=*/1.0,
           /*dt=*/0.1f, /*c_hyp=*/2.f);

  /* Leg 6: free-field regime (kappa = 0, no dust), grey lambda: the
   * ISRFCosmology example's `free_field` fixture, this file's own
   * regression floor from before the band-edge transfer landed. */
  run_case("free_field", /*a=*/1.0, /*H=*/1000.0, /*kappa_pe=*/0.f,
           /*kappa_lw=*/0.f, /*lambda_pe=*/1.0, /*lambda_lw=*/1.0,
           /*lambda_n_lw=*/1.0, /*dt=*/1.0f, /*c_hyp=*/1.f);

  /* Non-cosmological no-op: `cosmology_init_no_cosmo` leaves H exactly 0, so
   * the term (and the transfer, which the code's own a_LW <= 0.f guard --
   * see radiation_isrf.c -- makes the literal 0.f here since kappa_LW = 0
   * too) must not perturb the state at all, bit for bit, at ANY lambda. */
  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, /*a=*/1.0, /*H=*/0.0, /*lambda_pe=*/2.154,
              /*lambda_lw=*/6.508, /*lambda_n_lw=*/6.0);

  struct part p;
  set_part(&p, /*kappa_pe=*/0.f, /*kappa_lw=*/0.f, /*dt=*/0.5f, /*c_hyp=*/2.f);
  radiation_end_gradient_propagation(&p, &e);
  radiation_end_force_propagation(&p, &e);

  if (p.feedback_data.isrf_moment[ISRF_MOMENT_PE].u != u_PE_0 ||
      p.feedback_data.isrf_moment[ISRF_MOMENT_LW].u != u_LW_0 ||
      p.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u != u_LW_PHOTON_0)
    error("H = 0 is not a no-op on u: %.9e vs %.9e",
          p.feedback_data.isrf_moment[ISRF_MOMENT_PE].u, u_PE_0);
  for (int k = 0; k < 3; k++)
    if (p.feedback_data.isrf_moment[ISRF_MOMENT_PE].specific_flux[k] !=
            F_PE_0[k] ||
        p.feedback_data.isrf_moment[ISRF_MOMENT_LW].specific_flux[k] !=
            F_LW_0[k] ||
        p.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[k] !=
            F_LW_PHOTON_0[k])
      error("H = 0 is not a no-op on the specific flux");

  message("H = 0 leaves u and F untouched, bit for bit, at any lambda");

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
