/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @file testMagmaPrepareForce.c
 * @brief Closed-form checks of the per-particle part of the MAGMA scheme.
 *
 * No neighbour loop: one particle whose gradient sums are prescribed (an
 * identity C-matrix, so that the accumulated gradient of v is the gradient
 * itself, a prescribed shock-presence sum, Laplacian of u and neighbour
 * maxima) goes through hydro_prepare_force(), and its outputs are compared
 * with the closed forms of the scheme:
 *
 * 1. the convergence fraction and the shear limiter (pure shear, pure
 *    compression, rigid rotation, mixed flows, the shock-presence boost at
 *    R = -1, 0 and +1 and inside the dead zone, the limiter floor), and
 *    random velocity gradients against a double-precision reference;
 * 2. the trigger and the decay of the viscosity coefficient, its floors
 *    (the force loop's and the conduction cap's), the seeding of the time
 *    derivative after a spell without a valid C-matrix;
 * 3. the conduction coefficient: growth, decay, floor and the cap by the
 *    particle's own and its neighbours' shock coefficient, which the
 *    viscosity floor does not weaken;
 * 4. the feedback hook: the reset, the conduction switched off through the
 *    cap on the particle and on its neighbours, the monotonic release as the
 *    coefficient decays;
 * 5. the kinetic-feedback raise of the signal speed and the time-step
 *    criteria (Courant and acceleration).
 *
 * Runs in the configured dimension and does nothing for the other schemes.
 */

#if !defined(MAGMA_SPH)

int main(int argc, char *argv[]) {
  message("Not a MAGMA build: nothing to test.");
  return 0;
}

#else

static const int dim = (int)hydro_dimension;
static const float eta = 1.5f;

/* Deterministic pseudo-random numbers in [0, 1] and [-1, 1] */
static double rand_unit(void) { return rand() / (double)RAND_MAX; }
static double rand_sym(void) { return 2.0 * rand_unit() - 1.0; }

/**
 * @brief Fills the gradient sums of a particle (the gradient union, which
 * hydro_prepare_force() overwrites with the force data: refill before every
 * call).
 *
 * @param V The velocity gradient, V[k][l] = d v_k / d x_l.
 * @param R The shock-presence factor (the sum in units of the density).
 * @param laplace_u The SPH Laplacian of u.
 * @param c_max_ngb The largest sound speed among the neighbours.
 * @param alpha_shock_max_ngb The largest shock coefficient among them.
 */
static void fill_gradients(struct part *p, const double V[3][3],
                           const double R, const float laplace_u,
                           const float c_max_ngb,
                           const float alpha_shock_max_ngb) {
  sym_matrix_identity(&p->gradient.c_matrix_inv);
  for (int l = 0; l < 3; ++l) {
    p->gradient.gradient_vx[l] = (float)V[0][l];
    p->gradient.gradient_vy[l] = (float)V[1][l];
    p->gradient.gradient_vz[l] = (float)V[2][l];
    p->gradient.gradient_u[l] = 0.f;
    p->gradient.gradient_P[l] = 0.f;
  }
  p->gradient.shock_indicator = (float)R * p->rho;
  p->gradient.c_max_ngb = c_max_ngb;
  p->gradient.alpha_shock_max_ngb = alpha_shock_max_ngb;
  p->laplace_u = laplace_u;
}

/**
 * @brief A particle of unit density, smoothing length and internal energy
 * at the resolution eta, with a valid previous gradient and its
 * coefficients at the floors of the given properties.
 */
static void init_particle(struct part *p, const struct hydro_props *hp) {
  bzero(p, sizeof(struct part));
  p->id = 1;
  p->time_bin = 1;
  p->rho = 1.f;
  p->h = 1.f;
  p->u = 1.f;
  p->mass = pow_dimension(p->h / eta) * p->rho;
  p->fallback_flags = magma_fallback_none;
  p->div_v_previous = 0.f;
  p->alpha_shock =
      fminf(hp->viscosity.alpha_min, hydro_props_default_viscosity_alpha_min);
  p->alpha_visc = hp->viscosity.alpha_min;
  p->alpha_diff = hp->diffusion.alpha_min;
}

/**
 * @brief The convergence fraction and the shear limiter in double precision
 * from the velocity gradient and the raw shock-presence factor.
 */
static void reference_limiter(const struct part *p, const double V[3][3],
                              const double R_raw, const double limiter_floor,
                              double *fraction, double *limiter) {
  double div = 0.;
  for (int k = 0; k < dim; ++k) div += V[k][k];
  const double iso = div / dim;
  double S2 = 0., rot2 = 0.;
  for (int k = 0; k < dim; ++k) {
    for (int l = 0; l < dim; ++l) {
      const double S = 0.5 * (V[k][l] + V[l][k]) - (k == l ? iso : 0.);
      S2 += S * S;
      if (l > k) {
        const double Rot = 0.5 * (V[k][l] - V[l][k]);
        rot2 += 2. * Rot * Rot;
      }
    }
  }
  const double div2 = div * div;
  const double gradient_floor = 1e-12 * (div2 + S2 + rot2);
  *fraction = (div2 + S2 > gradient_floor) ? div2 / (div2 + S2) : 1.;

  /* The dead zone of the shock-presence factor */
  const double m_over_rho_hd = p->mass * pow_dimension(1. / p->h) / p->rho;
  const double noise =
      kernel_root * m_over_rho_hd + sqrt(m_over_rho_hd / kernel_norm);
  const double dz = fmin(0.9, const_viscosity_boost_dead_zone * noise);
  const double Rc = fmin(1., fmax(-1., R_raw));
  const double R = Rc < -dz ? (Rc + dz) / (1. - dz)
                   : Rc > dz ? (Rc - dz) / (1. - dz)
                             : 0.;
  const double boost = 2. * pow(1. - R, 4);
  const double bdiv2 = boost * boost * div2;
  *limiter =
      fmax((bdiv2 + S2 > gradient_floor) ? bdiv2 / (bdiv2 + S2) : 1.,
           limiter_floor);
}

/* Half-width of the dead zone, as the code computes it */
static double dead_zone(const struct part *p) {
  const double m_over_rho_hd = p->mass * pow_dimension(1. / p->h) / p->rho;
  const double noise =
      kernel_root * m_over_rho_hd + sqrt(m_over_rho_hd / kernel_norm);
  return fmin(0.9, const_viscosity_boost_dead_zone * noise);
}

/* Absolute tolerance for values below 1, relative above (a small limiter or
 * coefficient comes from a float cancellation of order-one terms) */
static int check_close(const double value, const double expected,
                       const double tol, const char *what) {
  const double err = fabs(value - expected) / (fabs(expected) + 1.);
  if (err > tol) {
    message("%s: %.8e, expected %.8e (err %.2e)", what, value, expected, err);
    return 1;
  }
  return 0;
}

/**
 * @brief Checks 1: the convergence fraction and the limiter.
 */
static int check_limiter(const struct hydro_props *hp_default,
                         const struct cosmology *cosmo,
                         const struct pressure_floor_props *floor) {

  int failed = 0;
  const double floors[2] = {0., 0.3};
  const double D = 0.3, s = 0.8;
  char what[200];

  for (int f = 0; f < 2; ++f) {
    struct hydro_props hp = *hp_default;
    hp.viscosity.limiter_floor = floors[f];
    struct part p;
    struct xpart xp;
    bzero(&xp, sizeof(struct xpart));
    init_particle(&p, &hp);
    const double dz = dead_zone(&p);

    /* The cases: V, R, expected fraction, expected limiter */
    struct {
      const char *name;
      double V[3][3];
      double R;
      double fraction, limiter;
      int min_dim;
    } cases[] = {
        {"pure shear", {{0, 1, 0}, {0, 0, 0}, {0, 0, 0}}, 0., 0., floors[f], 2},
        {"pure compression", {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, 0., 1., 1.,
         1},
        {"rigid rotation", {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}, 0., 1., 1., 2},
        {"mixed, R = 0",
         {{-D, s / 2, 0}, {s / 2, -D, 0}, {0, 0, -D}},
         0.,
         dim * dim * D * D / (dim * dim * D * D + s * s / 2),
         fmax(4 * dim * dim * D * D / (4 * dim * dim * D * D + s * s / 2),
              floors[f]),
         2},
        {"mixed, R = -1",
         {{-D, s / 2, 0}, {s / 2, -D, 0}, {0, 0, -D}},
         -1.,
         dim * dim * D * D / (dim * dim * D * D + s * s / 2),
         fmax(1024 * dim * dim * D * D /
                  (1024 * dim * dim * D * D + s * s / 2),
              floors[f]),
         2},
        {"mixed, R = +1",
         {{-D, s / 2, 0}, {s / 2, -D, 0}, {0, 0, -D}},
         1.,
         dim * dim * D * D / (dim * dim * D * D + s * s / 2),
         floors[f],
         2},
        {"mixed, R in the dead zone",
         {{-D, s / 2, 0}, {s / 2, -D, 0}, {0, 0, -D}},
         0.5 * dz,
         dim * dim * D * D / (dim * dim * D * D + s * s / 2),
         fmax(4 * dim * dim * D * D / (4 * dim * dim * D * D + s * s / 2),
              floors[f]),
         2},
        {"compression, R = +1", {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, 1., 1.,
         1., 1},
        {"compression, R = -1", {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, -1., 1.,
         1., 1},
        {"uniform flow", {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, 0., 1., 1., 1},
    };
    const int n_cases = sizeof(cases) / sizeof(cases[0]);

    for (int c = 0; c < n_cases; ++c) {
      if (dim < cases[c].min_dim) continue;
      /* Only the components of the configured dimension */
      double V[3][3];
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          V[k][l] = (k < dim && l < dim) ? cases[c].V[k][l] : 0.;
      init_particle(&p, &hp);
      fill_gradients(&p, V, cases[c].R, 0.f, 0.f, 0.f);
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, 0.f, 0.f);
      sprintf(what, "floor %.1f, %s: convergence fraction", floors[f],
              cases[c].name);
      failed += check_close(p.force.convergence_fraction, cases[c].fraction,
                            1e-6, what);
      sprintf(what, "floor %.1f, %s: limiter", floors[f], cases[c].name);
      failed +=
          check_close(p.force.viscosity_limiter, cases[c].limiter, 1e-6, what);
      /* The divergence recorded for the neighbours */
      double div = 0.;
      for (int k = 0; k < dim; ++k) div += V[k][k];
      sprintf(what, "floor %.1f, %s: div v", floors[f], cases[c].name);
      failed += check_close(p.div_v, div, 1e-6, what);
    }

    /* Random gradients and shock factors against the reference */
    int random_failed = 0;
    for (int t = 0; t < 500; ++t) {
      double V[3][3];
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          V[k][l] = (k < dim && l < dim) ? rand_sym() : 0.;
      const double R = 1.2 * rand_sym();
      init_particle(&p, &hp);
      fill_gradients(&p, V, R, 0.f, 0.f, 0.f);
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, 0.f, 0.f);
      double fraction, limiter;
      reference_limiter(&p, V, R, floors[f], &fraction, &limiter);
      sprintf(what, "floor %.1f, random %d: convergence fraction", floors[f],
              t);
      random_failed +=
          check_close(p.force.convergence_fraction, fraction, 1e-5, what);
      sprintf(what, "floor %.1f, random %d: limiter", floors[f], t);
      random_failed +=
          check_close(p.force.viscosity_limiter, limiter, 1e-5, what);
    }
    failed += random_failed;
  }
  message("limiter: %d failed", failed);
  return failed;
}

/**
 * @brief Checks 2: the trigger and decay of the viscosity coefficient.
 */
static int check_viscosity(const struct hydro_props *hp_default,
                           const struct cosmology *cosmo,
                           const struct pressure_floor_props *floor) {

  int failed = 0;
  const float floors[3] = {hydro_props_default_viscosity_alpha_min, 0.3f, 1.f};
  const double D_values[4] = {0.001, 0.01, 0.1, 1.};
  const float dt = 0.1f;
  char what[200];

  for (int f = 0; f < 3; ++f) {
    struct hydro_props hp = *hp_default;
    hp.viscosity.alpha_min = floors[f];
    const double cap_floor =
        fmin(floors[f], hydro_props_default_viscosity_alpha_min);
    struct part p;
    struct xpart xp;
    bzero(&xp, sizeof(struct xpart));

    for (int d = 0; d < 4; ++d) {
      const double D = D_values[d];
      double V[3][3] = {{0}};
      for (int k = 0; k < dim; ++k) V[k][k] = -D;

      init_particle(&p, &hp);
      const double H = kernel_gamma * p.h;
      const double c = gas_soundspeed_from_internal_energy(p.rho, p.u);
      const double decay = exp(-dt * 2. * hp.viscosity.length * c / H);

      /* First call: the divergence drops from 0 to -d D in dt */
      fill_gradients(&p, V, 0., 0.f, 0.f, 0.f);
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, dt, dt);
      const double A = H * H * (dim * D / dt); /* limiter 1 */
      const double alpha_loc = const_viscosity_alpha * A / (A + c * c);
      double alpha_v = alpha_loc > cap_floor
                           ? alpha_loc
                           : alpha_loc + (cap_floor - alpha_loc) * decay;
      sprintf(what, "floor %.1f, D %.3f, trigger: alpha_shock", floors[f], D);
      failed +=
          check_close(p.alpha_shock, fmax(alpha_v, cap_floor), 1e-6, what);
      sprintf(what, "floor %.1f, D %.3f, trigger: alpha_visc", floors[f], D);
      failed += check_close(p.alpha_visc, fmax(alpha_v, floors[f]), 1e-6, what);
      sprintf(what, "floor %.1f, D %.3f, trigger: div_v_previous", floors[f],
              D);
      failed += check_close(p.div_v_previous, -dim * D, 1e-6, what);

      /* Second call with the same gradient: no source, decay only */
      const double alpha_prev = fmax(alpha_v, cap_floor);
      fill_gradients(&p, V, 0., 0.f, 0.f, 0.f);
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, dt, dt);
      alpha_v = alpha_prev * decay;
      sprintf(what, "floor %.1f, D %.3f, decay: alpha_shock", floors[f], D);
      failed +=
          check_close(p.alpha_shock, fmax(alpha_v, cap_floor), 1e-6, what);
      sprintf(what, "floor %.1f, D %.3f, decay: alpha_visc", floors[f], D);
      failed += check_close(p.alpha_visc, fmax(alpha_v, floors[f]), 1e-6, what);

      /* No derivative on the first valid step after an inversion failure */
      init_particle(&p, &hp);
      p.fallback_flags = magma_fallback_inversion;
      fill_gradients(&p, V, 0., 0.f, 0.f, 0.f);
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, dt, dt);
      sprintf(what, "floor %.1f, D %.3f, after a failed inversion: alpha_shock",
              floors[f], D);
      failed += check_close(p.alpha_shock, cap_floor, 1e-6, what);
      sprintf(what, "floor %.1f, D %.3f, after a failed inversion: flags",
              floors[f], D);
      failed += check_close(p.fallback_flags, magma_fallback_none, 1e-6, what);
    }

    /* A diverging flow (div v > 0) never triggers */
    double V[3][3] = {{0}};
    for (int k = 0; k < dim; ++k) V[k][k] = 1.;
    init_particle(&p, &hp);
    fill_gradients(&p, V, 0., 0.f, 0.f, 0.f);
    hydro_prepare_force(&p, &xp, cosmo, &hp, floor, dt, dt);
    sprintf(what, "floor %.1f, expansion: alpha_shock", floors[f]);
    failed += check_close(p.alpha_shock, cap_floor, 1e-6, what);

    /* Random sequences stay within the bounds */
    init_particle(&p, &hp);
    int out_of_bounds = 0;
    for (int t = 0; t < 300; ++t) {
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          V[k][l] = (k < dim && l < dim) ? 20. * rand_sym() : 0.;
      fill_gradients(&p, V, 1.2 * rand_sym(), (float)(5. * rand_sym()), 1.5f,
                     (float)(cap_floor + (1. - cap_floor) * rand_unit()));
      hydro_prepare_force(&p, &xp, cosmo, &hp, floor, dt * (float)rand_unit(),
                          dt);
      if (!(p.alpha_visc >= floors[f] && p.alpha_visc <= const_viscosity_alpha))
        ++out_of_bounds;
      if (!(p.alpha_shock >= cap_floor && p.alpha_shock <= const_viscosity_alpha))
        ++out_of_bounds;
      /* The cap takes the conduction coefficient below its floor */
      if (!(p.alpha_diff >= 0.f && p.alpha_diff <= 1.f)) ++out_of_bounds;
      if (!(p.force.viscosity_limiter >= hp.viscosity.limiter_floor &&
            p.force.viscosity_limiter <= 1.f))
        ++out_of_bounds;
      if (!(p.force.convergence_fraction >= 0.f &&
            p.force.convergence_fraction <= 1.f))
        ++out_of_bounds;
    }
    if (out_of_bounds) {
      message("floor %.1f: %d coefficients out of bounds in a random sequence",
              floors[f], out_of_bounds);
      ++failed;
    }
  }
  message("viscosity coefficient: %d failed", failed);
  return failed;
}

/**
 * @brief Checks 3: the conduction coefficient.
 */
static int check_conduction(const struct hydro_props *hp_default,
                            const struct cosmology *cosmo,
                            const struct pressure_floor_props *floor) {

  int failed = 0;
  const float dt = 0.05f;
  const float C = 1.7f, a0 = 0.2f;
  const double V0[3][3] = {{0}};
  char what[200];
  struct part p;
  struct xpart xp;
  bzero(&xp, sizeof(struct xpart));

  /* Floor 0: growth, decay and the lower bound */
  struct hydro_props hp0 = *hp_default;
  hp0.diffusion.alpha_min = 0.f;
  const double cap_floor = fmin(hp0.viscosity.alpha_min,
                                hydro_props_default_viscosity_alpha_min);
  const float L_values[3] = {2.f, -5.f, 0.f};
  for (int l = 0; l < 3; ++l) {
    init_particle(&p, &hp0);
    p.alpha_diff = a0;
    const double H = kernel_gamma * p.h;
    const double c = gas_soundspeed_from_internal_energy(p.rho, p.u);
    fill_gradients(&p, V0, 0., L_values[l], C, 0.f);
    hydro_prepare_force(&p, &xp, cosmo, &hp0, floor, dt, dt);
    const double source = hp0.diffusion.beta * H * L_values[l] / sqrt(p.u);
    const double ratio = dt * (c + C) / H;
    double expected = (a0 + dt * source) / (1. + ratio);
    expected = fmin(1., fmax(0., expected));
    sprintf(what, "floor 0, laplace u %.0f: alpha_diff", L_values[l]);
    failed += check_close(p.alpha_diff, expected, 1e-6, what);
  }

  /* Floor 1 (default): the floor whatever the source */
  for (int l = 0; l < 2; ++l) {
    init_particle(&p, hp_default);
    fill_gradients(&p, V0, 0., L_values[l], C, 0.f);
    hydro_prepare_force(&p, &xp, cosmo, hp_default, floor, dt, dt);
    sprintf(what, "floor 1, laplace u %.0f: alpha_diff", L_values[l]);
    failed += check_close(p.alpha_diff, 1., 1e-6, what);
  }

  /* The cap: a neighbour at alpha_max switches the conduction off, one at
   * the cap's floor leaves it alone, one half-way halves it; with the
   * viscosity floor at alpha_max as well (the cap does not weaken) */
  struct hydro_props hp1 = *hp_default;
  hp1.viscosity.alpha_min = const_viscosity_alpha;
  const struct hydro_props *props[3] = {&hp0, hp_default, &hp1};
  const char *prop_names[3] = {"conduction floor 0", "defaults",
                               "viscosity floor 1"};
  for (int h = 0; h < 3; ++h) {
    const struct hydro_props *hp = props[h];
    const float ngb_values[3] = {const_viscosity_alpha, (float)cap_floor,
                                 (float)(0.5 * (cap_floor + 1.))};
    for (int n = 0; n < 3; ++n) {
      init_particle(&p, hp);
      p.alpha_diff = 1.f;
      fill_gradients(&p, V0, 0., 2.f, C, ngb_values[n]);
      hydro_prepare_force(&p, &xp, cosmo, hp, floor, dt, dt);
      /* Without the cap the coefficient would be at least the floor, and
       * with a positive source and floor 1 exactly 1 */
      double uncapped;
      if (hp->diffusion.alpha_min == 1.f) {
        uncapped = 1.;
      } else {
        const double H = kernel_gamma * p.h;
        const double c = gas_soundspeed_from_internal_energy(p.rho, p.u);
        const double source = hp->diffusion.beta * H * 2. / sqrt(p.u);
        const double ratio = dt * (c + C) / H;
        uncapped = fmin(1., (1. + dt * source) / (1. + ratio));
      }
      const double cap = 1. - (ngb_values[n] - cap_floor) / (1. - cap_floor);
      sprintf(what, "%s, neighbour at %.2f: alpha_diff", prop_names[h],
              ngb_values[n]);
      failed += check_close(p.alpha_diff, fmin(uncapped, cap), 1e-6, what);
    }
  }
  message("conduction coefficient: %d failed", failed);
  return failed;
}

/**
 * @brief Checks 4: the feedback hook.
 */
static int check_feedback_hook(const struct hydro_props *hp,
                               const struct cosmology *cosmo,
                               const struct pressure_floor_props *floor) {

  int failed = 0;
  const double V0[3][3] = {{0}};
  struct part p, q;
  struct xpart xp;
  bzero(&xp, sizeof(struct xpart));
  const double cap_floor =
      fmin(hp->viscosity.alpha_min, hydro_props_default_viscosity_alpha_min);

  /* The reset */
  init_particle(&p, hp);
  hydro_diffusive_feedback_reset(&p);
  failed += check_close(p.alpha_visc, const_viscosity_alpha, 1e-7,
                        "hook: alpha_visc");
  failed += check_close(p.alpha_shock, const_viscosity_alpha, 1e-7,
                        "hook: alpha_shock");
  failed += check_close(p.alpha_diff, 0., 1e-7, "hook: alpha_diff");

  /* The next force preparation without a time step: the cap holds the
   * conduction at zero although the floor would restore it */
  fill_gradients(&p, V0, 0., 2.f, 1.7f, 0.f);
  hydro_prepare_force(&p, &xp, cosmo, hp, floor, 0.f, 0.f);
  failed += check_close(p.alpha_diff, 0., 1e-7, "hook, dt = 0: alpha_diff");
  failed += check_close(p.alpha_visc, const_viscosity_alpha, 1e-7,
                        "hook, dt = 0: alpha_visc");

  /* A neighbour of the heated particle */
  init_particle(&q, hp);
  fill_gradients(&q, V0, 0., 2.f, 1.7f, p.alpha_shock);
  hydro_prepare_force(&q, &xp, cosmo, hp, floor, 0.f, 0.f);
  failed += check_close(q.alpha_diff, 0., 1e-7, "hook, neighbour: alpha_diff");

  /* The release as the coefficient decays: monotonic, complete once the
   * coefficient has reached the cap's floor */
  const double H = kernel_gamma * p.h;
  const double c = gas_soundspeed_from_internal_energy(p.rho, p.u);
  const float dt = 0.1f * H / c;
  const double decay = exp(-dt * 2. * hp->viscosity.length * c / H);
  double previous = 0.;
  int non_monotonic = 0;
  fill_gradients(&p, V0, 0., 2.f, 1.7f, 0.f);
  hydro_prepare_force(&p, &xp, cosmo, hp, floor, dt, dt);
  const double first_cap = (1. - decay) / (1. - cap_floor);
  if (p.alpha_diff > first_cap + 1e-6) {
    message("hook, first decay: alpha_diff %.3e above the cap %.3e",
            p.alpha_diff, first_cap);
    ++failed;
  }
  previous = p.alpha_diff;
  for (int t = 0; t < 2000; ++t) {
    fill_gradients(&p, V0, 0., 2.f, 1.7f, 0.f);
    hydro_prepare_force(&p, &xp, cosmo, hp, floor, dt, dt);
    if (p.alpha_diff < previous - 1e-6) ++non_monotonic;
    previous = p.alpha_diff;
  }
  if (non_monotonic) {
    message("hook: the conduction coefficient decreased %d times during the "
            "release",
            non_monotonic);
    ++failed;
  }
  failed += check_close(p.alpha_diff, hp->diffusion.alpha_min, 1e-6,
                        "hook, released: alpha_diff");
  failed += check_close(p.alpha_shock, cap_floor, 1e-6,
                        "hook, released: alpha_shock");
  message("feedback hook: %d failed", failed);
  return failed;
}

/**
 * @brief Checks 5: the signal speed and the time step.
 */
static int check_timestep(const struct hydro_props *hp,
                          const struct cosmology *cosmo) {

  int failed = 0;
  struct part p;
  struct xpart xp;
  bzero(&xp, sizeof(struct xpart));
  init_particle(&p, hp);
  const double H = kernel_gamma * p.h;
  const float c = 2.f, mu = 0.3f, dv = 1.5f;

  /* The Courant criterion */
  p.force.c_sig = c;
  p.force.mu_tilde = mu;
  for (int k = 0; k < 3; ++k) p.a_hydro[k] = 0.f;
  const double v_sig =
      2. * (c + 0.6 * const_viscosity_alpha * (c + 2. * mu));
  failed += check_close(hydro_get_signal_velocity(&p), v_sig, 1e-6,
                        "signal velocity");
  const double dt_courant = hp->CFL_condition * 2. * H / v_sig;
  failed += check_close(hydro_compute_timestep(&p, &xp, hp, cosmo), dt_courant,
                        1e-5, "Courant time step");

  /* The kinetic-feedback raise: beta dv on the signal velocity */
  hydro_set_v_sig_based_on_velocity_kick(&p, cosmo, dv);
  failed += check_close(hydro_get_signal_velocity(&p),
                        v_sig + const_viscosity_beta * dv, 1e-5,
                        "signal velocity after a kick");
  failed += check_close(hydro_compute_timestep(&p, &xp, hp, cosmo),
                        hp->CFL_condition * 2. * H /
                            (v_sig + const_viscosity_beta * dv),
                        1e-5, "Courant time step after a kick");

  /* The acceleration criterion */
  p.force.c_sig = 1e-10f;
  p.force.mu_tilde = 0.f;
  p.a_hydro[0] = 0.3f;
  p.a_hydro[1] = (dim > 1) ? -0.4f : 0.f;
  p.a_hydro[2] = 0.f;
  const double norm_a = sqrt(p.a_hydro[0] * p.a_hydro[0] +
                             p.a_hydro[1] * p.a_hydro[1]);
  failed += check_close(hydro_compute_timestep(&p, &xp, hp, cosmo),
                        hp->CFL_condition * sqrt(2. * H / norm_a), 1e-5,
                        "acceleration time step");
  message("time step: %d failed", failed);
  return failed;
}

int main(int argc, char *argv[]) {

  /* Trap the usual floating-point errors */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.viscosity.eta_crit = 1.f / eta;
  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct pressure_floor_props floor;
  bzero(&floor, sizeof(struct pressure_floor_props));

  srand(42);
  int failed = 0;
  failed += check_limiter(&hp, &cosmo, &floor);
  failed += check_viscosity(&hp, &cosmo, &floor);
  failed += check_conduction(&hp, &cosmo, &floor);
  failed += check_feedback_hook(&hp, &cosmo, &floor);
  failed += check_timestep(&hp, &cosmo);

  if (failed > 0) {
    message("%d checks failed", failed);
    return 1;
  }
  message("All checks passed.");
  return 0;
}

#endif /* MAGMA_SPH */
