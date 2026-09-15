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

/* The cosmological expansion term of the LW/FUV hyperbolic propagation:
 * `du/dt` and `dF/dt` each carry `-H` times the quantity itself, one power of
 * the Hubble rate, because `u` and `F` are mass-specific and therefore
 * already dilute with the physical gas density they are measured against.
 * Driven through the real ghost functions, with transport, sources and the
 * dissipation accumulators all held at zero, so only that term can move the
 * state. */
#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_isrf.h"

static const float u_FUV_0 = 3.0f;
static const float u_LW_0 = 0.8f;
static const float F_FUV_0[3] = {0.4f, -0.2f, 0.1f};
static const float F_LW_0[3] = {-0.15f, 0.05f, 0.3f};

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
 */
static void make_engine(struct engine *e, struct cosmology *cosmo,
                        struct feedback_props *fp, struct phys_const *pc,
                        double a, double H) {

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
 * @param kappa The dust absorption rate to give both bands.
 * @param dt The step the ghosts should integrate over.
 * @param c_hyp The hyperbolic propagation speed.
 */
static void set_part(struct part *p, float kappa, float dt, float c_hyp) {

  bzero(p, sizeof(struct part));
  p->h = 1.f;

  struct feedback_part_data *fd = &p->feedback_data;
  fd->dt_prev = dt;
  fd->c_hyp = c_hyp;
  fd->isrf_band[ISRF_BAND_FUV].kappa = kappa;
  fd->isrf_band[ISRF_BAND_LW].kappa = kappa;
  fd->rho_prev = 1.f;
  fd->isrf_band[ISRF_BAND_FUV].u_prev = u_FUV_0;
  fd->isrf_band[ISRF_BAND_LW].u_prev = u_LW_0;
  fd->isrf_band[ISRF_BAND_FUV].u = u_FUV_0;
  fd->isrf_band[ISRF_BAND_LW].u = u_LW_0;
  fd->isrf_band[ISRF_BAND_FUV].ngb_mean_abs_u_V = 1.f;
  fd->isrf_band[ISRF_BAND_LW].ngb_mean_abs_u_V = 1.f;
  for (int k = 0; k < 3; k++) {
    fd->isrf_band[ISRF_BAND_FUV].specific_flux[k] = F_FUV_0[k];
    fd->isrf_band[ISRF_BAND_LW].specific_flux[k] = F_LW_0[k];
  }
}

/**
 * @brief Fail unless two values agree to a relative tolerance.
 *
 * @param name Name of the quantity, for the failure message.
 * @param expected The analytic value.
 * @param value The value the code produced.
 * @param tol Relative tolerance.
 */
static void check_close(const char *name, float expected, float value,
                        float tol) {

  const float scale = fmaxf(fabsf(expected), 1e-30f);
  if (fabsf(value - expected) / scale > tol)
    error("%s: expected %.9e, got %.9e", name, expected, value);
}

/**
 * @brief One (a, H, kappa) case: run both ghosts once and compare `u` and `F`
 * against the analytic one-step decay at the total relaxation rate.
 *
 * @param a The scale factor.
 * @param H The Hubble rate.
 * @param kappa The dust absorption rate.
 * @param dt The step length.
 * @param c_hyp The hyperbolic propagation speed.
 */
static void run_case(double a, double H, float kappa, float dt, float c_hyp) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, a, H);

  struct part p;
  set_part(&p, kappa, dt, c_hyp);

  radiation_end_density_propagation(&p, &e);
  radiation_end_gradient_propagation(&p, &e);

  const float rate = c_hyp * kappa + (float)H;
  const float expected_decay = expf(-rate * dt);

  check_close("u_FUV", expected_decay * u_FUV_0,
              p.feedback_data.isrf_band[ISRF_BAND_FUV].u, 1e-5f);
  check_close("u_LW", expected_decay * u_LW_0,
              p.feedback_data.isrf_band[ISRF_BAND_LW].u, 1e-5f);
  for (int k = 0; k < 3; k++) {
    check_close("specific_flux_FUV", expected_decay * F_FUV_0[k],
                p.feedback_data.isrf_band[ISRF_BAND_FUV].specific_flux[k],
                1e-5f);
    check_close("specific_flux_LW", expected_decay * F_LW_0[k],
                p.feedback_data.isrf_band[ISRF_BAND_LW].specific_flux[k],
                1e-5f);
  }

  message("a = %g, H = %g, kappa = %g: decay %.8e as expected", a, H, kappa,
          expected_decay);
}

int main(int argc, char *argv[]) {

  /* Pure redshift, two epochs. `H` is NOT what vanishes at a = 1 (there it is
   * H_0), so both legs below carry a real term. */
  run_case(/*a=*/1.0, /*H=*/0.3, /*kappa=*/0.f, /*dt=*/0.5f, /*c_hyp=*/2.f);
  run_case(/*a=*/0.25, /*H=*/2.0, /*kappa=*/0.f, /*dt=*/0.1f, /*c_hyp=*/2.f);

  /* Stiff leg, the discriminating one: absorption depth 10 alongside a
   * redshift depth of 1. An explicit `-dt*phi*H*u` decrement instead of the
   * rate used here would land near `-0.1*u_0`, wrong sign and four orders
   * of magnitude out, rather than on `exp(-11)*u_0`. */
  run_case(/*a=*/0.5, /*H=*/10.0, /*kappa=*/50.f, /*dt=*/0.1f, /*c_hyp=*/2.f);

  /* Non-cosmological no-op: `cosmology_init_no_cosmo` leaves H exactly 0, so
   * the term must not perturb the state at all, bit for bit. */
  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, /*a=*/1.0, /*H=*/0.0);

  struct part p;
  set_part(&p, /*kappa=*/0.f, /*dt=*/0.5f, /*c_hyp=*/2.f);
  radiation_end_density_propagation(&p, &e);
  radiation_end_gradient_propagation(&p, &e);

  if (p.feedback_data.isrf_band[ISRF_BAND_FUV].u != u_FUV_0 ||
      p.feedback_data.isrf_band[ISRF_BAND_LW].u != u_LW_0)
    error("H = 0 is not a no-op on u: %.9e vs %.9e",
          p.feedback_data.isrf_band[ISRF_BAND_FUV].u, u_FUV_0);
  for (int k = 0; k < 3; k++)
    if (p.feedback_data.isrf_band[ISRF_BAND_FUV].specific_flux[k] !=
            F_FUV_0[k] ||
        p.feedback_data.isrf_band[ISRF_BAND_LW].specific_flux[k] != F_LW_0[k])
      error("H = 0 is not a no-op on the specific flux");

  message("H = 0 leaves u and F untouched, bit for bit");

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
