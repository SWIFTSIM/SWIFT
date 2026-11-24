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
 * @file testMagmaConservation.c
 * @brief Identities of the MAGMA force loop: conservation, symmetry and
 * invariance.
 *
 * Particles on a jittered lattice with random velocities of the order of the
 * sound speed, random internal energies and random (unequal) smoothing
 * lengths go through the full sequence of the scheme (density loop,
 * gradient loop and force preparation twice with a time step so that the
 * viscosity trigger, the shear limiter and the conduction cap take a range of
 * values, then the force loop). Without cosmology every pair term of the
 * force loop is antisymmetric in the momentum equation and its energy
 * partner cancels exactly between the two particles, in every formulation
 * and in the kernel-gradient fallbacks, so the checks are identities of the
 * discretisation rather than approximations:
 *
 * 1. the total momentum and the total energy (kinetic plus internal) are
 *    conserved to float round-off, in a periodic box (no boundary), with
 *    the particle-level fallback forced on a third of the particles and with
 *    the C-matrix regularisation forced on most of them;
 * 2. the symmetric and non-symmetric force functions agree with unequal
 *    smoothing lengths (testSymmetry uses equal ones);
 * 3. a constant velocity added to every particle leaves the accelerations,
 *    the heating rates and the time-step quantities unchanged (Galilean
 *    invariance: the reconstruction and the dissipation use differences
 *    only, the trigger the divergence only);
 * 4. a rotation of the configuration rotates the accelerations and leaves
 *    the scalars (heating, limiter, coefficients, fallbacks) unchanged,
 *    checked on the interior of a non-periodic block (a rotated periodic
 *    box is not periodic);
 * 5. the conduction term transports energy from the hotter to the colder
 *    particle of a pair and conserves it (single pairs with zero gradients
 *    and equal velocities, so that the pressure and viscous terms vanish).
 *
 * Runs in the configured dimension and does nothing for the other schemes.
 */

#if !defined(MAGMA_SPH)

int main(int argc, char *argv[]) {
  message("Not a MAGMA build: nothing to test.");
  return 0;
}

#else

/* Deterministic pseudo-random numbers in [0, 1] and [-1, 1] */
static double rand_unit(void) { return rand() / (double)RAND_MAX; }
static double rand_sym(void) { return 2.0 * rand_unit() - 1.0; }

static const double eta = 1.5;
static const double h_min_fac = 0.9, h_max_fac = 1.3;

/**
 * @brief Separation xi - xj, with the minimum image if periodic.
 */
static void separation(const struct part *pi, const struct part *pj,
                       const double L, const int periodic, float dxv[3]) {
  for (int k = 0; k < 3; ++k) {
    double d = pi->x[k] - pj->x[k];
    if (periodic) {
      if (d > 0.5 * L) d -= L;
      if (d < -0.5 * L) d += L;
    }
    dxv[k] = (float)d;
  }
}

/**
 * @brief Allocates and fills a particle set: jittered lattice, random
 * velocities, energies and smoothing lengths, equal masses.
 */
static void make_particles(struct part **parts_out, struct xpart **xparts_out,
                           float **v_out, const int n, const double L,
                           const double jitter) {

  const int N = (int)pow(n, hydro_dimension);
  const double dx = L / n;
  const double mass = pow(dx, hydro_dimension);

  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  if (posix_memalign((void **)&parts, part_align, N * sizeof(struct part)) != 0)
    error("Failed to allocate the particles");
  if (posix_memalign((void **)&xparts, xpart_align, N * sizeof(struct xpart)) !=
      0)
    error("Failed to allocate the extended particles");
  bzero(parts, N * sizeof(struct part));
  bzero(xparts, N * sizeof(struct xpart));
  float *v = (float *)malloc(3 * N * sizeof(float));

  for (int i = 0; i < N; ++i) {
    struct part *p = &parts[i];
    int idx = i;
    for (int k = 0; k < 3; ++k) p->x[k] = 0.;
    for (int k = 0; k < (int)hydro_dimension; ++k) {
      p->x[k] = (idx % n + 0.5) * dx + jitter * dx * rand_sym();
      idx /= n;
    }
    for (int k = 0; k < 3; ++k)
      v[3 * i + k] = (k < (int)hydro_dimension) ? 1.5 * rand_sym() : 0.f;
    p->u = 1. + 3. * rand_unit();
    p->mass = mass;
    p->h = eta * dx * (h_min_fac + (h_max_fac - h_min_fac) * rand_unit());
    p->id = i + 1;
  }

  *parts_out = parts;
  *xparts_out = xparts;
  *v_out = v;
}

/**
 * @brief Runs the scheme's sequence on a particle set.
 *
 * Density loop; gradient loop and force preparation with the velocities
 * divided by v_scale, then again with the velocities themselves (a time
 * derivative of the divergence for the viscosity trigger, and a valid
 * previous gradient for the shock-presence sums); optional forced fallbacks;
 * force loop (symmetric over unordered pairs, or non-symmetric over ordered
 * pairs, as the tasks call them); force finishing.
 *
 * @param parts The particles (positions, u, mass, h set).
 * @param xparts The extended particles.
 * @param N Number of particles.
 * @param L Box size.
 * @param periodic Minimum image or not.
 * @param v The velocities of the force loop.
 * @param v_scale The velocities of the first gradient pass are v / v_scale.
 * @param dt_alpha Time step of the viscosity trigger and decay.
 * @param nonsym_force Use the non-symmetric force function.
 * @param forced_fallbacks Force the inversion fallback on 30 % of the
 * particles after the force preparation.
 */
static void run_sequence(struct part *parts, struct xpart *xparts, const int N,
                         const double L, const int periodic, const float *v,
                         const float v_scale, const float dt_alpha,
                         const int nonsym_force, const int forced_fallbacks,
                         const struct hydro_props *hp,
                         const struct cosmology *cosmo,
                         const struct pressure_floor_props *floor) {

  for (int i = 0; i < N; ++i) {
    struct part *p = &parts[i];
    for (int k = 0; k < 3; ++k) p->v[k] = v[3 * i + k] / v_scale;
    hydro_first_init_part(p, &xparts[i]);
    p->time_bin = 1;
    hydro_init_part(p, NULL);
  }

  /* Density loop: every ordered pair within the kernel support of i */
  for (int i = 0; i < N; ++i) {
    struct part *pi = &parts[i];
    const float Hi = kernel_gamma * pi->h;
    const float hig2 = Hi * Hi;
    for (int j = 0; j < N; ++j) {
      if (i == j) continue;
      struct part *pj = &parts[j];
      float dxv[3];
      separation(pi, pj, L, periodic, dxv);
      const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
      if (r2 < hig2)
        runner_iact_nonsym_density(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
    }
  }
  for (int i = 0; i < N; ++i) {
    hydro_end_density(&parts[i], cosmo);
    hydro_prepare_gradient(&parts[i], &xparts[i], cosmo, hp, floor);
    hydro_reset_gradient(&parts[i]);
  }

  /* Two gradient passes */
  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 1) {
      for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) parts[i].v[k] = v[3 * i + k];
        hydro_reset_gradient(&parts[i]);
      }
    }
    for (int i = 0; i < N; ++i) {
      struct part *pi = &parts[i];
      const float Hi = kernel_gamma * pi->h;
      const float hig2 = Hi * Hi;
      for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        struct part *pj = &parts[j];
        float dxv[3];
        separation(pi, pj, L, periodic, dxv);
        const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
        if (r2 < hig2)
          runner_iact_nonsym_gradient(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
      }
    }
    for (int i = 0; i < N; ++i) {
      hydro_end_gradient(&parts[i]);
      hydro_prepare_force(&parts[i], &xparts[i], cosmo, hp, floor, dt_alpha,
                          dt_alpha);
    }
  }

  if (forced_fallbacks)
    for (int i = 0; i < N; ++i)
      if (i % 10 < 3) parts[i].fallback_flags |= magma_fallback_inversion;

  /* Force loop: pairs within the larger of the two supports */
  for (int i = 0; i < N; ++i) hydro_reset_acceleration(&parts[i]);
  for (int i = 0; i < N; ++i) {
    struct part *pi = &parts[i];
    const float Hi = kernel_gamma * pi->h;
    const float hig2 = Hi * Hi;
    for (int j = nonsym_force ? 0 : i + 1; j < N; ++j) {
      if (i == j) continue;
      struct part *pj = &parts[j];
      const float Hj = kernel_gamma * pj->h;
      const float hjg2 = Hj * Hj;
      float dxv[3];
      separation(pi, pj, L, periodic, dxv);
      const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
      if (r2 < hig2 || r2 < hjg2) {
        if (nonsym_force)
          runner_iact_nonsym_force(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
        else
          runner_iact_force(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
      }
    }
  }
  for (int i = 0; i < N; ++i) hydro_end_force(&parts[i], cosmo);
}

/**
 * @brief Prints how much of the scheme the configuration exercises.
 */
static void report_coverage(const struct part *parts, const int N,
                            const struct hydro_props *hp, const char *label) {
  int n_alpha_up = 0, n_limited = 0, n_regularised = 0, n_base = 0,
      n_capped = 0;
  long pair_fallbacks = 0;
  for (int i = 0; i < N; ++i) {
    const struct part *p = &parts[i];
    if (p->alpha_visc > hp->viscosity.alpha_min + 1e-6f) ++n_alpha_up;
    if (p->force.viscosity_limiter < 0.99f) ++n_limited;
    if (p->fallback_flags & magma_fallback_regularised) ++n_regularised;
    if (p->fallback_flags & magma_fallback_base_sph_mask) ++n_base;
    if (p->alpha_diff < 0.99f) ++n_capped;
    pair_fallbacks += p->force.n_pair_fallbacks;
  }
  message(
      "%s: %d particles; alpha above the floor %d, limiter below 1 %d, "
      "conduction capped %d, regularised %d, base SPH %d, pair fallbacks %ld",
      label, N, n_alpha_up, n_limited, n_capped, n_regularised, n_base,
      pair_fallbacks);
}

/**
 * @brief Checks the conservation of momentum and energy.
 *
 * @return The number of failed checks.
 */
static int check_conservation(const struct part *parts, const int N,
                              const double tol, const char *label) {

  double sum_ma[3] = {0., 0., 0.}, sum_m_abs_a[3] = {0., 0., 0.};
  double sum_E = 0., sum_abs_E = 0.;
  for (int i = 0; i < N; ++i) {
    const struct part *p = &parts[i];
    const double m = p->mass;
    double v_dot_a = 0., norm_v = 0., norm_a = 0.;
    for (int k = 0; k < 3; ++k) {
      sum_ma[k] += m * p->a_hydro[k];
      sum_m_abs_a[k] += m * fabs(p->a_hydro[k]);
      v_dot_a += p->v[k] * p->a_hydro[k];
      norm_v += p->v[k] * p->v[k];
      norm_a += p->a_hydro[k] * p->a_hydro[k];
    }
    sum_E += m * (v_dot_a + p->u_dt);
    sum_abs_E += m * (sqrt(norm_v * norm_a) + fabs(p->u_dt));
  }

  int failed = 0;
  double worst_p = 0.;
  for (int k = 0; k < (int)hydro_dimension; ++k)
    worst_p = fmax(worst_p, fabs(sum_ma[k]) / sum_m_abs_a[k]);
  const double err_E = fabs(sum_E) / sum_abs_E;
  message("%s: momentum error %.2e, energy error %.2e (tolerance %.0e)", label,
          worst_p, err_E, tol);
  if (worst_p > tol) {
    message("%s: momentum not conserved", label);
    ++failed;
  }
  if (err_E > tol) {
    message("%s: energy not conserved", label);
    ++failed;
  }
  return failed;
}

/**
 * @brief Compares two runs particle by particle.
 *
 * The accelerations of the second run must be Q times those of the first;
 * the scalars must agree. Vector quantities are compared relative to their
 * own norm plus the rms norm over the compared particles, scalars to their
 * magnitude plus one.
 *
 * @param Q The rotation (identity for a comparison without rotation).
 * @param mask Which particles to compare (NULL: all).
 * @param max_outliers Particles allowed beyond the tolerance (the scheme's
 * discrete switches can flip at round-off level on a borderline pair).
 * @return The number of failed checks.
 */
static int compare_runs(const struct part *a, const struct part *b,
                        const int N, const double Q[3][3], const int *mask,
                        const double tol, const int max_outliers,
                        const char *label) {

  /* Scales */
  double rms_a = 0., rms_udt = 0., rms_hdt = 0., rms_mu = 0.;
  int checked = 0;
  for (int i = 0; i < N; ++i) {
    if (mask && !mask[i]) continue;
    ++checked;
    for (int k = 0; k < 3; ++k) rms_a += a[i].a_hydro[k] * a[i].a_hydro[k];
    rms_udt += a[i].u_dt * a[i].u_dt;
    rms_hdt += a[i].force.h_dt * a[i].force.h_dt;
    rms_mu += a[i].force.mu_tilde * a[i].force.mu_tilde;
  }
  if (checked == 0) error("%s: nothing to compare", label);
  rms_a = sqrt(rms_a / checked);
  rms_udt = sqrt(rms_udt / checked);
  rms_hdt = sqrt(rms_hdt / checked);
  rms_mu = sqrt(rms_mu / checked);

  const char *names[] = {"a_hydro",
                         "u_dt",
                         "h_dt",
                         "mu_tilde",
                         "convergence_fraction",
                         "viscosity_limiter",
                         "alpha_visc",
                         "alpha_diff",
                         "div_v",
                         "c_sig"};
  const int n_quant = 10;
  double worst[10] = {0.};
  int outliers[10] = {0};
  int flag_mismatch = 0, fallback_mismatch = 0;

  for (int i = 0; i < N; ++i) {
    if (mask && !mask[i]) continue;
    const struct part *pa = &a[i];
    const struct part *pb = &b[i];

    /* Rotated acceleration */
    double Qa[3] = {0., 0., 0.}, norm_a = 0., diff_a = 0.;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) Qa[k] += Q[k][l] * pa->a_hydro[l];
      norm_a += pa->a_hydro[k] * pa->a_hydro[k];
    }
    for (int k = 0; k < 3; ++k)
      diff_a += (Qa[k] - pb->a_hydro[k]) * (Qa[k] - pb->a_hydro[k]);
    const double err[10] = {
        sqrt(diff_a) / (sqrt(norm_a) + rms_a),
        fabs(pa->u_dt - pb->u_dt) / (fabs(pa->u_dt) + rms_udt),
        fabs(pa->force.h_dt - pb->force.h_dt) /
            (fabs(pa->force.h_dt) + rms_hdt),
        fabs(pa->force.mu_tilde - pb->force.mu_tilde) /
            (fabs(pa->force.mu_tilde) + rms_mu),
        fabs(pa->force.convergence_fraction - pb->force.convergence_fraction),
        fabs(pa->force.viscosity_limiter - pb->force.viscosity_limiter),
        fabs(pa->alpha_visc - pb->alpha_visc),
        fabs(pa->alpha_diff - pb->alpha_diff),
        fabs(pa->div_v - pb->div_v) / (fabs(pa->div_v) + 1.),
        fabs(pa->force.c_sig - pb->force.c_sig) / (fabs(pa->force.c_sig) + 1.)};
    for (int q = 0; q < n_quant; ++q) {
      worst[q] = fmax(worst[q], err[q]);
      if (err[q] > tol) {
        if (outliers[q] < 3)
          message("%s: particle %lld: %s differs (err %.2e)", label, pa->id,
                  names[q], err[q]);
        ++outliers[q];
      }
    }
    if (pa->fallback_flags != pb->fallback_flags) ++flag_mismatch;
    if (pa->force.n_pair_fallbacks != pb->force.n_pair_fallbacks)
      ++fallback_mismatch;
  }

  int failed = 0;
  for (int q = 0; q < n_quant; ++q) {
    if (outliers[q] > max_outliers) {
      message("%s: %s differs on %d of %d particles (worst %.2e)", label,
              names[q], outliers[q], checked, worst[q]);
      ++failed;
    }
  }
  if (flag_mismatch > max_outliers) {
    message("%s: fallback flags differ on %d particles", label,
            flag_mismatch);
    ++failed;
  }
  if (fallback_mismatch > max_outliers) {
    message("%s: pair fallback counts differ on %d particles", label,
            fallback_mismatch);
    ++failed;
  }
  message(
      "%s: %d particles compared, worst errors: a %.1e, u_dt %.1e, h_dt %.1e, "
      "mu~ %.1e, xi_c %.1e, xi %.1e, alpha %.1e, alpha_u %.1e, div v %.1e, "
      "c_sig %.1e; flag mismatches %d, pair-fallback mismatches %d",
      label, checked, worst[0], worst[1], worst[2], worst[3], worst[4],
      worst[5], worst[6], worst[7], worst[8], worst[9], flag_mismatch,
      fallback_mismatch);
  return failed;
}

/**
 * @brief A random rotation of the configured dimension (a reflection in 1D).
 */
static void random_rotation(double Q[3][3]) {
  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l) Q[k][l] = (k == l) ? 1. : 0.;
#if defined(HYDRO_DIMENSION_1D)
  Q[0][0] = -1.;
#elif defined(HYDRO_DIMENSION_2D)
  const double theta = M_PI * rand_sym();
  Q[0][0] = cos(theta);
  Q[0][1] = -sin(theta);
  Q[1][0] = sin(theta);
  Q[1][1] = cos(theta);
#else
  /* From a random unit quaternion */
  double q[4];
  double norm = 0.;
  for (int k = 0; k < 4; ++k) {
    q[k] = rand_sym();
    norm += q[k] * q[k];
  }
  norm = sqrt(norm);
  for (int k = 0; k < 4; ++k) q[k] /= norm;
  const double w = q[0], x = q[1], y = q[2], z = q[3];
  Q[0][0] = 1. - 2. * (y * y + z * z);
  Q[0][1] = 2. * (x * y - z * w);
  Q[0][2] = 2. * (x * z + y * w);
  Q[1][0] = 2. * (x * y + z * w);
  Q[1][1] = 1. - 2. * (x * x + z * z);
  Q[1][2] = 2. * (y * z - x * w);
  Q[2][0] = 2. * (x * z - y * w);
  Q[2][1] = 2. * (y * z + x * w);
  Q[2][2] = 1. - 2. * (x * x + y * y);
#endif
}

/**
 * @brief Fills a particle for a single-pair interaction: no gradients, an
 * identity C-matrix, no limiting.
 */
static void make_pair_particle(struct part *p, const double x[3],
                               const float v[3], const float m,
                               const float rho, const float u, const float h,
                               const char flags) {
  bzero(p, sizeof(struct part));
  for (int k = 0; k < 3; ++k) {
    p->x[k] = x[k];
    p->v[k] = v[k];
  }
  p->mass = m;
  p->rho = rho;
  p->u = u;
  p->h = h;
  p->id = 1;
  p->time_bin = 1;
  sym_matrix_identity(&p->force.c_matrix);
  p->force.pressure = gas_pressure_from_internal_energy(rho, u);
  p->force.soundspeed = gas_soundspeed_from_internal_energy(rho, u);
  p->force.convergence_fraction = 1.f;
  p->force.viscosity_limiter = 1.f;
  p->force.eta_crit = 1.f / eta;
  p->alpha_visc = hydro_props_default_viscosity_alpha_min;
  p->alpha_shock = hydro_props_default_viscosity_alpha_min;
  p->alpha_diff = 1.f;
  p->fallback_flags = flags;
  hydro_reset_acceleration(p);
}

/**
 * @brief Second law of the conduction term on single pairs.
 *
 * Zero gradients and equal velocities: the pressure and viscous terms
 * vanish, the conduction acts on the raw jump of u with the pressure-residual
 * signal speed (the pressures differ by more than its threshold). The
 * hotter particle must lose energy, the colder gain it, and the pair must
 * conserve it, with the matrix gradient functions and with the kernel
 * gradients (fallback), in both orderings and in both force functions.
 *
 * @return The number of failed checks.
 */
static int check_second_law(const int n_pairs) {

  int failed = 0;
  for (int t = 0; t < n_pairs; ++t) {

    const char flags =
        (t % 2 == 0) ? magma_fallback_none : magma_fallback_inversion;
    const float mi = 0.5f + 1.5f * rand_unit();
    const float mj = 0.5f + 1.5f * rand_unit();
    const float rhoi = 0.5f + 1.5f * rand_unit();
    const float rhoj = 0.5f + 1.5f * rand_unit();
    const float hi = 0.8f + 0.4f * rand_unit();
    const float hj = 0.8f + 0.4f * rand_unit();
    float ui, uj;
    do {
      ui = 0.5f + 3.5f * rand_unit();
      uj = 0.5f + 3.5f * rand_unit();
    } while (fabsf(rhoi * ui - rhoj * uj) <
             0.2f * 0.5f * (rhoi * ui + rhoj * uj));

    /* A separation inside the larger support */
    const float H_max = kernel_gamma * fmaxf(hi, hj);
    float dxv[3] = {0.f, 0.f, 0.f};
    double norm = 0.;
    for (int k = 0; k < (int)hydro_dimension; ++k) {
      dxv[k] = rand_sym();
      norm += dxv[k] * dxv[k];
    }
    const float r = (0.3f + 0.6f * rand_unit()) * H_max;
    for (int k = 0; k < 3; ++k) dxv[k] *= r / sqrt(norm);
    const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
    const double xi[3] = {dxv[0], dxv[1], dxv[2]};
    const double xj[3] = {0., 0., 0.};
    const float v[3] = {0.7f, -0.4f, 0.2f};

    struct part pi, pj, pi2, pj2, pi3, pj3;
    make_pair_particle(&pi, xi, v, mi, rhoi, ui, hi, flags);
    make_pair_particle(&pj, xj, v, mj, rhoj, uj, hj, flags);
    pi2 = pi;
    pj2 = pj;
    pi3 = pi;
    pj3 = pj;

    /* Symmetric, then non-symmetric in both orderings */
    runner_iact_force(r2, dxv, hi, hj, &pi, &pj, 1.f, 0.f);
    runner_iact_nonsym_force(r2, dxv, hi, hj, &pi2, &pj2, 1.f, 0.f);
    const float dxv_rev[3] = {-dxv[0], -dxv[1], -dxv[2]};
    runner_iact_nonsym_force(r2, dxv_rev, hj, hi, &pj3, &pi3, 1.f, 0.f);

    const double E = mi * pi.u_dt + mj * pj.u_dt;
    const double E_abs = fabs(mi * pi.u_dt) + fabs(mj * pj.u_dt);
    const int sign_ok = (ui > uj) ? (pi.u_dt < 0.f && pj.u_dt > 0.f)
                                  : (pi.u_dt > 0.f && pj.u_dt < 0.f);
    const int conserved = fabs(E) <= 1e-5 * E_abs;
    const int nonsym_ok =
        fabs(pi2.u_dt - pi.u_dt) <= 1e-5 * fabs(pi.u_dt) &&
        fabs(pj3.u_dt - pj.u_dt) <= 1e-5 * fabs(pj.u_dt);
    /* The pressure force of the pair must be antisymmetric too */
    double p_sum = 0., p_abs = 0.;
    for (int k = 0; k < 3; ++k) {
      p_sum += fabs(mi * pi.a_hydro[k] + mj * pj.a_hydro[k]);
      p_abs += fabs(mi * pi.a_hydro[k]) + fabs(mj * pj.a_hydro[k]);
    }
    const int momentum_ok = p_sum <= 1e-5 * p_abs;

    if (!sign_ok || !conserved || !nonsym_ok || !momentum_ok) {
      message(
          "pair %d (%s): u_i %.3f u_j %.3f, u_dt_i %.3e u_dt_j %.3e (nonsym "
          "%.3e %.3e), energy %.2e, momentum %.2e: %s%s%s%s",
          t, flags ? "kernel gradients" : "gradient functions", ui, uj,
          pi.u_dt, pj.u_dt, pi2.u_dt, pj3.u_dt, E / E_abs, p_sum / p_abs,
          sign_ok ? "" : "wrong direction ", conserved ? "" : "not conserved ",
          nonsym_ok ? "" : "nonsym differs ",
          momentum_ok ? "" : "momentum not conserved");
      ++failed;
    }
  }
  message("second law: %d pairs, %d failed", n_pairs, failed);
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

  /* The box: the largest support below a quarter of the periodic box (16
   * per dimension with the cubic spline, 18 with the quintic) */
#if defined(HYDRO_DIMENSION_1D)
  const int n = 256, n_block = 256;
#elif defined(HYDRO_DIMENSION_2D)
  const int n = 48, n_block = 64;
#else
  const int n = (int)ceil(4. * kernel_gamma * h_max_fac * eta) + 1,
            n_block = 24;
#endif
  const double L = 1.0;
  const double jitter = 0.3;
  const double dx = L / n;
  const float dt_alpha = 0.3f * dx;
  const double tol_cons = 1e-5;
  const double tol_cmp = 1e-5;

  int failed = 0;

  /* The periodic configuration */
  srand(42);
  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  float *v = NULL;
  make_particles(&parts, &xparts, &v, n, L, jitter);
  const int N = (int)pow(n, hydro_dimension);
  const double L_min_support = 4. * kernel_gamma * h_max_fac * eta * dx;
  if (L_min_support > L)
    error("The box is too small for the periodic minimum image");

  struct part *ref = NULL;
  if (posix_memalign((void **)&ref, part_align, N * sizeof(struct part)) != 0)
    error("Failed to allocate the reference particles");
  const double I[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

  /* 1. Conservation */
  run_sequence(parts, xparts, N, L, /*periodic=*/1, v, 1.3f, dt_alpha,
               /*nonsym=*/0, /*forced=*/0, &hp, &cosmo, &floor);
  report_coverage(parts, N, &hp, "periodic box");
  failed += check_conservation(parts, N, tol_cons, "symmetric force loop");
  memcpy(ref, parts, N * sizeof(struct part));

  /* 2. Symmetric vs non-symmetric with unequal h */
  run_sequence(parts, xparts, N, L, 1, v, 1.3f, dt_alpha, /*nonsym=*/1, 0,
               &hp, &cosmo, &floor);
  failed += check_conservation(parts, N, tol_cons, "non-symmetric force loop");
  failed += compare_runs(ref, parts, N, I, NULL, tol_cmp, 0,
                         "symmetric vs non-symmetric");

  /* 3. Galilean invariance */
  float *v_boost = (float *)malloc(3 * N * sizeof(float));
  const float V[3] = {3.f, -2.f, 1.f};
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < 3; ++k)
      v_boost[3 * i + k] =
          v[3 * i + k] + ((k < (int)hydro_dimension) ? 1.5f * V[k] : 0.f);
  run_sequence(parts, xparts, N, L, 1, v_boost, 1.3f, dt_alpha, 0, 0, &hp,
               &cosmo, &floor);
  failed += compare_runs(ref, parts, N, I, NULL, tol_cmp, N / 1000,
                         "Galilean boost");

  /* 4. Conservation with forced fallbacks and with the regularisation */
  run_sequence(parts, xparts, N, L, 1, v, 1.3f, dt_alpha, 0, /*forced=*/1,
               &hp, &cosmo, &floor);
  report_coverage(parts, N, &hp, "forced fallbacks");
  failed += check_conservation(parts, N, tol_cons, "forced fallbacks");
  /* A condition-number limit just above 1: the jittered lattice is nearly
   * isotropic, so most of the inverses go through the regularised branch */
  struct hydro_props hp_reg = hp;
  hp_reg.viscosity.max_condition_number = 1.1f;
  run_sequence(parts, xparts, N, L, 1, v, 1.3f, dt_alpha, 0, 0, &hp_reg,
               &cosmo, &floor);
  report_coverage(parts, N, &hp_reg, "regularised C-matrices");
#if !defined(HYDRO_DIMENSION_1D)
  /* (in 1D the C-matrix is a scalar: nothing to regularise) */
  int n_regularised = 0;
  for (int i = 0; i < N; ++i)
    if (parts[i].fallback_flags & magma_fallback_regularised) ++n_regularised;
  if (n_regularised < N / 2) {
    message("regularised C-matrices: only %d of %d regularised", n_regularised,
            N);
    ++failed;
  }
#endif
  failed += check_conservation(parts, N, tol_cons, "regularised C-matrices");
  free(parts);
  free(xparts);
  free(v);
  free(v_boost);
  free(ref);

  /* 5. Rotational invariance on the interior of a non-periodic block */
  srand(43);
  make_particles(&parts, &xparts, &v, n_block, L, jitter);
  const int N_block = (int)pow(n_block, hydro_dimension);
  const double dx_block = L / n_block;
  const double margin = 2. * kernel_gamma * h_max_fac * eta * dx_block;
  int *mask = (int *)malloc(N_block * sizeof(int));
  int n_interior = 0;
  for (int i = 0; i < N_block; ++i) {
    mask[i] = 1;
    for (int k = 0; k < (int)hydro_dimension; ++k)
      if (parts[i].x[k] < margin || parts[i].x[k] > L - margin) mask[i] = 0;
    n_interior += mask[i];
  }
  if (n_interior == 0) error("No interior particle in the block");

  struct part *rotated = NULL;
  struct xpart *xrotated = NULL;
  if (posix_memalign((void **)&rotated, part_align,
                     N_block * sizeof(struct part)) != 0)
    error("Failed to allocate the rotated particles");
  if (posix_memalign((void **)&xrotated, xpart_align,
                     N_block * sizeof(struct xpart)) != 0)
    error("Failed to allocate the rotated extended particles");
  memcpy(rotated, parts, N_block * sizeof(struct part));
  bzero(xrotated, N_block * sizeof(struct xpart));
  float *v_rot = (float *)malloc(3 * N_block * sizeof(float));
  double Q[3][3];
  random_rotation(Q);
  for (int i = 0; i < N_block; ++i) {
    for (int k = 0; k < 3; ++k) {
      double xr = 0.5 * L, vr = 0.;
      for (int l = 0; l < 3; ++l) {
        xr += Q[k][l] * (parts[i].x[l] - 0.5 * L);
        vr += Q[k][l] * v[3 * i + l];
      }
      rotated[i].x[k] = xr;
      v_rot[3 * i + k] = (float)vr;
    }
  }

  run_sequence(parts, xparts, N_block, L, /*periodic=*/0, v, 1.3f, dt_alpha, 0,
               0, &hp, &cosmo, &floor);
  report_coverage(parts, N_block, &hp, "block");
  run_sequence(rotated, xrotated, N_block, L, 0, v_rot, 1.3f, dt_alpha, 0, 0,
               &hp, &cosmo, &floor);
  failed += compare_runs(parts, rotated, N_block, Q, mask, 1e-4,
                         n_interior / 1000, "rotation");
  free(parts);
  free(xparts);
  free(v);
  free(rotated);
  free(xrotated);
  free(v_rot);
  free(mask);

  /* 6. Second law of the conduction term */
  srand(44);
  failed += check_second_law(200);

  if (failed > 0) {
    message("%d checks failed", failed);
    return 1;
  }
  message("All checks passed.");
  return 0;
}

#endif /* MAGMA_SPH */
