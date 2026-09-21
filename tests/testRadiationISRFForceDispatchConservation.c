/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

/* The ISRF flux divergence and artificial dissipation through the real
 * type-2 force-loop dispatch, on two cells whose smoothing lengths differ by
 * a factor of two. Both operators are antisymmetric pair exchanges, so over
 * an all-active same-bin system sum_i m_i div(F)_i and sum_i m_i diss_i must
 * vanish, and each particle must receive every pair with r < max(H_i, H_j).
 * The cell pair is swept once without depth limits and once split across two
 * depth levels, as the sub-cell recursion does. */
#if defined(FEEDBACK_GEAR) && defined(SPHENIX_SPH)

#include "feedback/GEAR/radiation_propagation_iact.h"

#define NODE_ID 0
#define CELL_N 6

/* A lost mirrored side at a factor-two h step leaves an imbalance of order
 * 0.1 of sum m |term|, and a lost pair an error of order 1e-2 of a
 * particle's summed absolute pair terms. */
#define SUM_BAR 1e-5
#define PER_PART_BAR 1e-5

void runner_dopair2_branch_force(struct runner *r, struct cell *ci,
                                 struct cell *cj, int limit_h_min,
                                 int limit_h_max);
void runner_doself2_branch_force(struct runner *r, const struct cell *c,
                                 int limit_h_min, int limit_h_max);

/**
 * @brief Bit-exact float equality, via the raw bits rather than `==`: under
 * -ffast-math the compiler is free to fold a comparison, so read the bit
 * pattern explicitly instead of trusting the operator.
 *
 * @param x First value.
 * @param y Second value.
 * @return 1 if the two floats have identical bit patterns.
 */
static int bits_equal_f(float x, float y) {
  uint32_t bx, by;
  memcpy(&bx, &x, sizeof(bx));
  memcpy(&by, &y, sizeof(by));
  return bx == by;
}

/**
 * @brief Kernel-local #feedback_part_data.c_hyp: same-bin identity
 * (design-required to be bit-identical to the shipped per-particle
 * formula), the light-speed cap, and the debug pin, all through the real
 * #radiation_end_density_propagation dispatch.
 *
 * @param e The #engine (feedback_props, cosmology, physical_constants set).
 */
static void test_c_hyp_same_bin_cap_and_pin(struct engine *e) {

  /* Same-bin identity: every neighbour (there are none registered here)
   * shares this particle's own bin, so #feedback_part_data.max_ngb_time_bin
   * left at its own #part.time_bin reduces dt_max(i) to dt_i exactly, and
   * #radiation_end_density_propagation must give the shipped per-particle
   * `min(C_hyp*h/dt_i, c)` bit for bit. */
  struct part p;
  bzero(&p, sizeof(struct part));
  p.h = 0.37f;
  p.time_bin = 6;
  p.feedback_data.dt_prev = (float)get_timestep(p.time_bin, e->time_base);
  p.feedback_data.max_ngb_time_bin = p.time_bin;

  radiation_end_density_propagation(&p, e);

  const float h_phys = (float)e->cosmology->a * p.h;
  float expected =
      e->feedback_props->ISRF_c_hyp_margin * h_phys / p.feedback_data.dt_prev;
  expected = min(expected, (float)e->physical_constants->const_speed_light_c);
  if (e->feedback_props->ISRF_c_hyp_pin_for_debugging > 0.f)
    expected = e->feedback_props->ISRF_c_hyp_pin_for_debugging;

  if (!bits_equal_f(p.feedback_data.c_hyp, expected))
    error("same-bin identity: c_hyp = %.9g != shipped formula %.9g",
          (double)p.feedback_data.c_hyp, (double)expected);
  message(
      "same-bin identity: c_hyp = %.6g bit-identical to the shipped "
      "per-particle formula",
      (double)p.feedback_data.c_hyp);

  /* Light-speed cap: a huge h over a short step would give C_hyp*h/dt_max
   * far above c without the clamp. */
  struct part p_cap;
  bzero(&p_cap, sizeof(struct part));
  p_cap.h = 1e6f;
  p_cap.time_bin =
      1; /* get_integer_timestep(1) == 4, a short but nonzero step */
  p_cap.feedback_data.dt_prev =
      (float)get_timestep(p_cap.time_bin, e->time_base);
  p_cap.feedback_data.max_ngb_time_bin = p_cap.time_bin;

  radiation_end_density_propagation(&p_cap, e);

  const float c = (float)e->physical_constants->const_speed_light_c;
  if (!(p_cap.feedback_data.c_hyp <= c))
    error("light-speed cap: c_hyp = %.6g exceeds c = %.6g",
          (double)p_cap.feedback_data.c_hyp, (double)c);
  message("light-speed cap: c_hyp = %.6g capped at c = %.6g",
          (double)p_cap.feedback_data.c_hyp, (double)c);

  /* Debug pin: overrides the clamped value unconditionally. */
  const float saved_pin = e->feedback_props->ISRF_c_hyp_pin_for_debugging;
  e->feedback_props->ISRF_c_hyp_pin_for_debugging = 123.f;

  struct part p_pin;
  bzero(&p_pin, sizeof(struct part));
  p_pin.h = 1.f;
  p_pin.time_bin = 3;
  p_pin.feedback_data.dt_prev =
      (float)get_timestep(p_pin.time_bin, e->time_base);
  p_pin.feedback_data.max_ngb_time_bin = p_pin.time_bin;

  radiation_end_density_propagation(&p_pin, e);

  if (!bits_equal_f(p_pin.feedback_data.c_hyp, 123.f))
    error("debug pin: c_hyp = %.6g, expected the pinned value 123",
          (double)p_pin.feedback_data.c_hyp);
  message("debug pin: c_hyp = %.6g pinned as configured",
          (double)p_pin.feedback_data.c_hyp);

  e->feedback_props->ISRF_c_hyp_pin_for_debugging = saved_pin;
}

/**
 * @brief Receiver bound: for every registered neighbour j, `c_i*dt_j/h_i
 * <= C_hyp` must hold by construction, because #feedback_part_data.
 * max_ngb_time_bin is a maximum over the kernel (dt_max(i) >= dt_j for any
 * j that contributed to it). Sweeps time-bin differences 1, 2, 3 and
 * smoothing-length ratios 0.5, 1, 2 (h_j plays no role in #c_hyp's own
 * formula: only h_i and dt_max(i) do; the ratio is swept to document that
 * explicitly, not because the bound depends on it). Restricted to
 * neighbours that actually reached the kernel sum, as the design's own
 * "one coverage gap" note requires (a neighbour with `H_j > r >= H_i`
 * reads #c_hyp in the force loop without ever contributing to
 * #max_ngb_time_bin, and is not covered by this bound).
 *
 * @param e The #engine (feedback_props, cosmology, physical_constants set,
 * debug pin OFF).
 */
static void test_c_hyp_receiver_bound(struct engine *e) {

  if (e->feedback_props->ISRF_c_hyp_pin_for_debugging > 0.f)
    error("receiver bound: test setup requires the debug pin off");

  const float C_hyp = e->feedback_props->ISRF_c_hyp_margin;
  const float h_i = 0.42f;
  const timebin_t bin_i = 4;
  const int bin_diffs[3] = {1, 2, 3};
  const float h_ratios[3] = {0.5f, 1.f, 2.f};

  for (int bd = 0; bd < 3; bd++) {
    for (int hr = 0; hr < 3; hr++) {
      struct part p;
      bzero(&p, sizeof(struct part));
      p.h = h_i;
      p.time_bin = bin_i;
      p.feedback_data.dt_prev = (float)get_timestep(bin_i, e->time_base);
      p.feedback_data.max_ngb_time_bin = bin_i;

      const timebin_t bin_j = bin_i + bin_diffs[bd];
      /* h_j itself is unused below: it plays no role in the formula, only
       * in whether j is registered as a neighbour at all (simulated here
       * by directly folding bin_j into max_ngb_time_bin, as the density
       * loop's accumulation would for a genuine in-kernel neighbour). */
      (void)h_ratios[hr];
      p.feedback_data.max_ngb_time_bin =
          max(p.feedback_data.max_ngb_time_bin, bin_j);

      radiation_end_density_propagation(&p, e);

      const double dt_j = get_timestep(bin_j, e->time_base);
      const double bound = (double)p.feedback_data.c_hyp * dt_j / (double)h_i;
      if (!(bound <= (double)C_hyp * (1. + 1e-6)))
        error(
            "receiver bound: bin_diff=%d h_ratio=%.2f: c_i*dt_j/h_i = "
            "%.6e above C_hyp = %.6e",
            bin_diffs[bd], (double)h_ratios[hr], bound, (double)C_hyp);
    }
  }
  message(
      "receiver bound: c_i*dt_j/h_i <= C_hyp holds for bin differences "
      "1..3 and h ratios 0.5/1/2");
}

/**
 * @brief Density-loop #feedback_part_data.max_ngb_time_bin: the
 * per-h-iteration reset (#radiation_init_part_propagation) and the
 * running maximum through both the symmetric
 * (#runner_iact_isrf_propagation) and non-symmetric
 * (#runner_iact_nonsym_isrf_propagation) density-loop hooks.
 */
static void test_density_loop_max_ngb_time_bin(void) {

  struct part pi, pj, pk;
  bzero(&pi, sizeof(struct part));
  bzero(&pj, sizeof(struct part));
  bzero(&pk, sizeof(struct part));

  pi.h = 0.5f;
  pi.time_bin = 3;
  pj.h = 0.5f;
  pj.time_bin = 7;
  pk.h = 0.5f;
  pk.time_bin = 1;

  struct part *const triplet[3] = {&pi, &pj, &pk};
  for (int t = 0; t < 3; t++) {
    struct part *p = triplet[t];
    p->mass = 1.f;
    p->feedback_data.rho_prev = 1.f;
    for (int b = 0; b < ISRF_BAND_COUNT; b++)
      p->feedback_data.isrf_band[b].u_prev = 0.1f;
  }

  const float dx[3] = {0.01f, 0.f, 0.f};
  const float r2 = dx[0] * dx[0];

  /* Per-h-iteration reset: each particle's max starts at its own bin. */
  radiation_init_part_propagation(&pi);
  radiation_init_part_propagation(&pj);
  radiation_init_part_propagation(&pk);
  if (pi.feedback_data.max_ngb_time_bin != pi.time_bin ||
      pj.feedback_data.max_ngb_time_bin != pj.time_bin ||
      pk.feedback_data.max_ngb_time_bin != pk.time_bin)
    error(
        "density-loop max: reset did not seed max_ngb_time_bin from "
        "each particle's own time_bin");

  /* Symmetric hook: both sides pick up the other's bin. pi (3) vs pk (1):
   * pi's max stays 3 (pk is faster), pk's max rises to 3. */
  runner_iact_isrf_propagation(r2, dx, pi.h, pk.h, &pi, &pk, 1.f, 0.f, NULL);
  if (pi.feedback_data.max_ngb_time_bin != 3)
    error(
        "density-loop max: symmetric hook changed i's max from a slower "
        "neighbour (%d, expected 3)",
        (int)pi.feedback_data.max_ngb_time_bin);
  if (pk.feedback_data.max_ngb_time_bin != 3)
    error(
        "density-loop max: symmetric hook did not raise j's max to i's "
        "bin (%d, expected 3)",
        (int)pk.feedback_data.max_ngb_time_bin);

  /* Non-symmetric hook: only i's accumulator is touched. Two calls
   * accumulate (do not overwrite), against a slower and a faster
   * neighbour in turn; j's own field must never change. */
  radiation_init_part_propagation(&pi);
  const timebin_t pj_sentinel = pj.feedback_data.max_ngb_time_bin;
  runner_iact_nonsym_isrf_propagation(r2, dx, pi.h, pj.h, &pi, &pj, 1.f, 0.f,
                                      NULL);
  if (pi.feedback_data.max_ngb_time_bin != 7)
    error(
        "density-loop max: non-symmetric hook did not raise i's max to "
        "the slower neighbour's bin (%d, expected 7)",
        (int)pi.feedback_data.max_ngb_time_bin);
  if (pj.feedback_data.max_ngb_time_bin != pj_sentinel)
    error("density-loop max: non-symmetric hook wrote to j's own field");

  runner_iact_nonsym_isrf_propagation(r2, dx, pi.h, pk.h, &pi, &pk, 1.f, 0.f,
                                      NULL);
  if (pi.feedback_data.max_ngb_time_bin != 7)
    error(
        "density-loop max: non-symmetric hook overwrote the running "
        "maximum with a slower call's smaller bin (%d, expected 7)",
        (int)pi.feedback_data.max_ngb_time_bin);

  /* Redo path: a second h-iteration's reset must bring the maximum back
   * down to i's own bin, not leave the previous iteration's value. */
  radiation_init_part_propagation(&pi);
  if (pi.feedback_data.max_ngb_time_bin != pi.time_bin)
    error(
        "density-loop max: reset on a redo iteration left a stale "
        "maximum (%d, expected i's own bin %d)",
        (int)pi.feedback_data.max_ngb_time_bin, (int)pi.time_bin);

  message(
      "density-loop max: reset and both hooks' running maximum are "
      "correct");
}

/**
 * @brief Run the kernel-local #feedback_part_data.c_hyp unit tests: same-bin
 * identity, the light-speed cap, the debug pin, the receiver bound, and the
 * density-loop maximum/reset. Independent of the force-dispatch conservation
 * tests below, which exercise #feedback_part_data.c_hyp only as an opaque
 * per-particle input.
 */
static void test_kernel_local_c_hyp(void) {

  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  fb_props.ISRF_propagation = 1;
  /* radiation_end_density_propagation is a no-op for any other scheme
   * (see its own doxygen): these tests exercise it directly, so they must
   * select the kernel-local scheme explicitly rather than rely on the
   * struct's zero-init default (isrf_c_hyp_scheme_shipped). */
  fb_props.ISRF_c_hyp_scheme = isrf_c_hyp_scheme_kernel_local;
  fb_props.ISRF_c_hyp_margin = 0.5f;
  fb_props.ISRF_c_hyp_pin_for_debugging = 0.f;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_speed_light_c = 3e5; /* km/s-scale test units */

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.feedback_props = &fb_props;
  e.cosmology = &cosmo;
  e.physical_constants = &phys_const;
  e.time_base = 1e-3;
  e.ti_current = 8;
  e.policy = 0; /* no cosmology */

  test_c_hyp_same_bin_cap_and_pin(&e);
  test_c_hyp_receiver_bound(&e);
  test_density_loop_max_ngb_time_bin();
}

/**
 * @brief Build a cell of CELL_N^3 perturbed-lattice gas particles, all
 * active, with random ISRF state.
 *
 * @param offset The cell's lower corner.
 * @param h_spacing Base smoothing length in units of the particle spacing.
 * @param depth_h Depth tag given to every particle.
 * @param part_id (in/out) Running particle id.
 * @return The cell.
 */
static struct cell *make_cell(const double offset[3], double h_spacing,
                              char depth_h, long long *part_id) {

  const size_t count = CELL_N * CELL_N * CELL_N;
  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate the cell");
  bzero(c, sizeof(struct cell));
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("Couldn't allocate the particles");
  bzero(c->hydro.parts, count * sizeof(struct part));
  if (posix_memalign((void **)&c->hydro.xparts, part_align,
                     count * sizeof(struct xpart)) != 0)
    error("Couldn't allocate the extra particles");
  bzero(c->hydro.xparts, count * sizeof(struct xpart));

  float h_max = 0.f;
  struct part *p = c->hydro.parts;
  for (int x = 0; x < CELL_N; x++) {
    for (int y = 0; y < CELL_N; y++) {
      for (int z = 0; z < CELL_N; z++) {
        const int idx[3] = {x, y, z};
        for (int k = 0; k < 3; k++)
          p->x[k] =
              offset[k] + (idx[k] + 0.5 + random_uniform(-0.2, 0.2)) / CELL_N;
        p->h = h_spacing * random_uniform(1., 1.2) / CELL_N;
        h_max = max(h_max, p->h);
        p->id = ++(*part_id);
        p->depth_h = depth_h;
        p->mass = random_uniform(0.5, 2.);
        p->u = 1.f;
        p->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
        p->ti_drift = 8;
        p->ti_kick = 8;
#endif

        struct feedback_part_data *fd = &p->feedback_data;
        fd->rho_prev = random_uniform(0.5, 2.);
        fd->c_hyp = random_uniform(0.5, 1.5);
        for (int b = 0; b < ISRF_BAND_COUNT; b++) {
          struct feedback_isrf_band_data *band = &fd->isrf_band[b];
          band->u = random_uniform(-0.2, 1.);
          for (int k = 0; k < 3; k++)
            band->specific_flux[k] = random_uniform(-1., 1.);
          band->dissipation_alpha_trigger = random_uniform(0., 0.5);
          band->dissipation_alpha_floor = random_uniform(0., 0.5);
        }
        p++;
      }
    }
  }

  c->split = 0;
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max;
  c->hydro.count = count;
  for (int k = 0; k < 3; k++) {
    c->width[k] = 1.;
    c->loc[k] = offset[k];
  }
  c->dmin = 1.;
  c->hydro.super = c;
  c->hydro.ti_old_part = 8;
  c->hydro.ti_end_min = 8;
  c->nodeID = NODE_ID;
  return c;
}

/**
 * @brief Free a cell built by #make_cell.
 *
 * @param c The cell.
 */
static void clean_up(struct cell *c) {
  cell_free_hydro_sorts(c);
  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c);
}

/**
 * @brief Set the SPH force-loop inputs and zero both ISRF accumulators.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void prepare_force(struct cell *c, const struct engine *e) {

  for (int i = 0; i < c->hydro.count; i++) {
    struct part *p = &c->hydro.parts[i];
    p->rho = 1.f;
    p->density.rho_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
    p->force.pressure = hydro_get_comoving_pressure(p);
    p->viscosity.alpha = 0.8;
    p->viscosity.div_v = 0.f;
    p->viscosity.div_v_previous_step = 0.f;
    p->viscosity.v_sig = hydro_get_comoving_soundspeed(p);
    hydro_prepare_force(p, &c->hydro.xparts[i], e->cosmology,
                        e->hydro_properties, e->pressure_floor_props, 0., 0.);
    hydro_reset_acceleration(p);
    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      p->feedback_data.isrf_band[b].div_specific_flux = 0.f;
      p->feedback_data.isrf_band[b].dissipation_u = 0.f;
    }
  }
}

/**
 * @brief Set a cell's depth and the smoothing-length window its depth-limited
 * sweeps check against.
 *
 * @param c The cell.
 * @param depth The depth.
 * @param h_min_allowed Lower smoothing-length bound at this depth.
 * @param h_max_allowed Upper smoothing-length bound at this depth.
 */
static void set_level(struct cell *c, char depth, float h_min_allowed,
                      float h_max_allowed) {
  c->depth = depth;
  c->h_min_allowed = h_min_allowed;
  c->h_max_allowed = h_max_allowed;
}

/**
 * @brief Finiteness from the bit pattern: the build's -ffast-math lets the
 * compiler fold isfinite() to true.
 *
 * @param x The value.
 * @return 1 if the exponent field is not all ones.
 */
static int is_finite_bits(double x) {
  uint64_t bits;
  memcpy(&bits, &x, sizeof(bits));
  return ((bits >> 52) & 0x7FF) != 0x7FF;
}

/**
 * @brief Fail unless the mass-weighted sum of one accumulator vanishes
 * relative to its mass-weighted absolute sum.
 *
 * @param name The quantity, for messages.
 * @param sum sum_i m_i x_i.
 * @param abs_sum sum_i m_i |x_i|.
 * @param expect_zero 1 to require the ratio below #SUM_BAR, 0 to require it
 * above ten times #SUM_BAR (a sweep that drops one depth level on purpose).
 */
static void check_sum(const char *name, double sum, double abs_sum,
                      int expect_zero) {
  if (!is_finite_bits(sum) || !is_finite_bits(abs_sum) || !(abs_sum > 0.))
    error("%s: degenerate sums %e / %e", name, sum, abs_sum);
  const double ratio = fabs(sum) / abs_sum;
  if (expect_zero && !(ratio <= SUM_BAR))
    error("%s: sum m x / sum m |x| = %e above %e", name, ratio, SUM_BAR);
  if (!expect_zero && !(ratio > 10. * SUM_BAR))
    error("%s: a one-sided sweep gives %e, the metric cannot see a lost side",
          name, ratio);
  message("%s: sum m x / sum m |x| = %.3e", name, ratio);
}

/**
 * @brief Compare every particle's accumulators with a brute-force sum over
 * all pairs with r < max(H_i, H_j), and check both mass-weighted sums.
 *
 * @param cells The two cells.
 * @param label Case name, for messages.
 * @param expect_zero See #check_sum.
 */
static void check_cells(struct cell *cells[2], const char *label,
                        int expect_zero) {

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    double sum_div = 0., abs_div = 0., sum_diss = 0., abs_diss = 0.;
    double max_ref_div = 0., max_ref_diss = 0., max_err_div = 0.,
           max_err_diss = 0., min_cancel_div = 1.;

    double max_masked_err_div = 0., max_masked_err_diss = 0.;
    for (int ci = 0; ci < 2; ci++) {
      for (int i = 0; i < cells[ci]->hydro.count; i++) {
        const struct part *pi = &cells[ci]->hydro.parts[i];
        const struct feedback_isrf_band_data *bi =
            &pi->feedback_data.isrf_band[b];
        sum_div += pi->mass * (double)bi->div_specific_flux;
        abs_div += pi->mass * fabs((double)bi->div_specific_flux);
        sum_diss += pi->mass * (double)bi->dissipation_u;
        abs_diss += pi->mass * fabs((double)bi->dissipation_u);

        if (!expect_zero) continue;

        /* Reference in double from the same pair hook, one pair at a time.
         * The error is measured against the sum of the absolute pair terms:
         * the net is a cancelling sum of random-sign terms, so float32
         * rounding relative to the net alone is not a coverage signal. */
        double ref_div = 0., ref_diss = 0., abs_terms_div = 0.,
               abs_terms_diss = 0.;
        /* Negative control: a reference that masks every pair with
         * r >= H_i, as a dispatch that only ever visited i's own kernel
         * (missing the H_i <= r < H_j reach extension) would. Proves the
         * per-particle metric below is sensitive to that specific bug. */
        const float Hi = kernel_gamma * pi->h;
        double masked_ref_div = 0., masked_ref_diss = 0.;
        for (int cj = 0; cj < 2; cj++) {
          for (int j = 0; j < cells[cj]->hydro.count; j++) {
            const struct part *pj = &cells[cj]->hydro.parts[j];
            if (pj == pi) continue;
            const float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                                 (float)(pi->x[1] - pj->x[1]),
                                 (float)(pi->x[2] - pj->x[2])};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            const float H = kernel_gamma * max(pi->h, pj->h);
            if (!(r2 < H * H)) continue;
            struct part tmp = *pi;
            tmp.feedback_data.isrf_band[b].div_specific_flux = 0.f;
            tmp.feedback_data.isrf_band[b].dissipation_u = 0.f;
            runner_iact_nonsym_isrf_dissipation(r2, dx, pi->h, pj->h, &tmp, pj,
                                                1.f, 0.f);
            const double d_div =
                tmp.feedback_data.isrf_band[b].div_specific_flux;
            const double d_diss = tmp.feedback_data.isrf_band[b].dissipation_u;
            ref_div += d_div;
            ref_diss += d_diss;
            abs_terms_div += fabs(d_div);
            abs_terms_diss += fabs(d_diss);
            if (r2 < Hi * Hi) {
              masked_ref_div += d_div;
              masked_ref_diss += d_diss;
            }
          }
        }
        max_ref_div = max(max_ref_div, abs_terms_div);
        max_ref_diss = max(max_ref_diss, abs_terms_diss);
        const double scale_div = max(abs_terms_div, 1e-30);
        const double scale_diss = max(abs_terms_diss, 1e-30);
        const double err_div =
            fabs((double)bi->div_specific_flux - ref_div) / scale_div;
        const double err_diss =
            fabs((double)bi->dissipation_u - ref_diss) / scale_diss;
        max_err_div = max(max_err_div, err_div);
        max_err_diss = max(max_err_diss, err_diss);
        const double cancel_div = fabs(ref_div) / scale_div;
        min_cancel_div = min(min_cancel_div, cancel_div);
        const double masked_err_div =
            fabs((double)bi->div_specific_flux - masked_ref_div) / scale_div;
        const double masked_err_diss =
            fabs((double)bi->dissipation_u - masked_ref_diss) / scale_diss;
        max_masked_err_div = max(max_masked_err_div, masked_err_div);
        max_masked_err_diss = max(max_masked_err_diss, masked_err_diss);
      }
    }

    char name[64];
    sprintf(name, "%s band %d div(F)", label, b);
    check_sum(name, sum_div, abs_div, expect_zero);
    sprintf(name, "%s band %d dissipation", label, b);
    check_sum(name, sum_diss, abs_diss, expect_zero);

    if (expect_zero) {
      if (!(max_ref_div > 0.) || !(max_ref_diss > 0.))
        error("%s band %d: brute-force reference is zero", label, b);
      if (!(max_err_div <= PER_PART_BAR) || !(max_err_diss <= PER_PART_BAR))
        error("%s band %d: per-particle mismatch div %e diss %e above %e",
              label, b, max_err_div, max_err_diss, PER_PART_BAR);
      message(
          "%s band %d: per-particle max error / sum |pair terms| div %.2e "
          "diss %.2e (smallest |net| / sum |terms| for div %.2e)",
          label, b, max_err_div, max_err_diss, min_cancel_div);

      /* r >= H_i masked negative control: proves the per-particle metric
         above would have caught a dispatch that dropped the H_i <= r < H_j
         reach extension, by exceeding PER_PART_BAR many times over. */
      if (!(max_masked_err_div > 10. * PER_PART_BAR) ||
          !(max_masked_err_diss > 10. * PER_PART_BAR))
        error(
            "%s band %d: r>=H_i masked negative control margin div %.1fx "
            "diss %.1fx does not exceed 10x",
            label, b, max_masked_err_div / PER_PART_BAR,
            max_masked_err_diss / PER_PART_BAR);
      message(
          "%s band %d: r>=H_i masked negative control margin (x bar) div "
          "%.1f diss %.1f",
          label, b, max_masked_err_div / PER_PART_BAR,
          max_masked_err_diss / PER_PART_BAR);
    }
  }
}

/* ---------------------------------------------------------------------
 * The uniform reduced light-speed candidate
 * (ISRF_c_hyp_fixed_fraction_of_c): off leaves c_hyp bit-identical to the
 * shipped C_hyp*h/dt formula; on gives every particle exactly f*c,
 * independent of h and time bin; and the matching radiation timestep term
 * (radiation_isrf_part_timestep) returns C_hyp*h/(f*c) on an eligible
 * particle and is inert (FLT_MAX) whenever the fraction is off.
 * ------------------------------------------------------------------- */

/**
 * @brief Minimal fixture for radiation_snapshot_part_propagation /
 * radiation_isrf_part_timestep: a non-cosmological engine with just enough
 * state to reach the c_hyp closure (cooling_func/internal_units are needed
 * because the kappa computation ahead of it in
 * radiation_snapshot_part_propagation is unconditional whenever
 * ISRF_propagation is on), independent of the dispatch test's own engine
 * above.
 *
 * @param e (return) The engine.
 * @param cosmo (return) Non-cosmological (a=1) cosmology.
 * @param pc (return) Physical constants.
 * @param fp (return) Feedback properties.
 * @param cooling (return) Cooling function data.
 * @param us (return) Unit system.
 * @param fixed_fraction ISRF_c_hyp_fixed_fraction_of_c.
 * @param timestep_off_for_debugging
 * ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging.
 */
static void make_c_hyp_speed_test_engine(
    struct engine *e, struct cosmology *cosmo, struct phys_const *pc,
    struct feedback_props *fp, struct cooling_function_data *cooling,
    struct unit_system *us, float fixed_fraction,
    char timestep_off_for_debugging) {

  bzero(cosmo, sizeof(struct cosmology));
  cosmo->a = 1.0;
  cosmo->a2_inv = 1.0;
  cosmo->a3_inv = 1.0;

  bzero(pc, sizeof(struct phys_const));
  pc->const_speed_light_c = 1.e4;

  bzero(fp, sizeof(struct feedback_props));
  fp->ISRF_propagation = 1;
  fp->ISRF_extinction_path_in_kernel_radii = 2.0f;
  fp->ISRF_c_hyp_margin = 0.5f;
  /* The two schemes are alternatives (feedback_props_init() enforces this
   * at parse time; these tests call the dispatch directly, so they must
   * keep the same pairing by hand): a positive fraction only takes effect
   * under isrf_c_hyp_scheme_fixed_fraction, else radiation_snapshot_part_
   * propagation runs the shipped formula regardless of this value. */
  fp->ISRF_c_hyp_scheme = fixed_fraction > 0.f
                              ? isrf_c_hyp_scheme_fixed_fraction
                              : isrf_c_hyp_scheme_shipped;
  fp->ISRF_c_hyp_fixed_fraction_of_c = fixed_fraction;
  fp->ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging =
      timestep_off_for_debugging;

  bzero(cooling, sizeof(struct cooling_function_data));
  cooling->chemistry_data.local_dust_to_gas_ratio = 0.01;

  units_init(us, /*U_M_in_cgs=*/1.98892e33, /*U_L_in_cgs=*/3.08567758e18,
             /*U_t_in_cgs=*/3.15576e13, /*U_C_in_cgs=*/1.0,
             /*U_T_in_cgs=*/1.0);

  bzero(e, sizeof(struct engine));
  e->policy = 0; /* no cosmology: the plain get_timestep branch */
  e->ti_current = 0;
  e->time_base = 1.0;
  e->max_active_bin = num_time_bins; /* every bin active */
  e->internal_units = us;
  e->physical_constants = pc;
  e->cosmology = cosmo;
  e->cooling_func = cooling;
  e->feedback_props = fp;
}

/**
 * @brief Set a particle's h/time_bin/rho, leaving every ISRF band field at
 * its zero-init default.
 *
 * @param p (return) The particle.
 * @param h The smoothing length.
 * @param time_bin The time bin.
 */
static void set_c_hyp_test_part(struct part *p, float h, timebin_t time_bin) {
  bzero(p, sizeof(struct part));
  p->h = h;
  p->time_bin = time_bin;
  p->rho = 1.5f;
  p->mass = 1.f;
  p->feedback_data.ISRF_reservoir_end_ti = -1;
  p->feedback_data.ISRF_illumination_end_ti = -1;
}

/**
 * @brief Relative-tolerance check for a two-operator float formula
 * (multiply-then-divide) reproduced independently in this test: not a hard
 * `==`, since -ffast-math/-freciprocal-math may pick a different reciprocal
 * instruction across translation units for a bit-for-bit identical source
 * expression (see swift-knowledge.md's floating-point-hazards notes). A
 * single-multiply comparison (fraction*c) has no such ambiguity and is
 * checked with `==` instead.
 *
 * @param name Message label.
 * @param actual The value read back from the code under test.
 * @param expected This test's own independently-computed value.
 */
static void assert_close_c_hyp(const char *name, float actual, float expected) {
  const float rel_err = fabsf(actual - expected) / fabsf(expected);
  if (!(rel_err <= 1e-6f))
    error("%s: got %.8e, expected %.8e (rel_err=%.3e).", name, (double)actual,
          (double)expected, (double)rel_err);
}

/**
 * @brief Fraction off (default): c_hyp must be bit-identical to the shipped
 * min(C_hyp*h/dt, c) formula, for several (h, time_bin) combinations.
 */
static void test_c_hyp_fixed_fraction_off_matches_shipped_formula(void) {
  struct engine e;
  struct cosmology cosmo;
  struct phys_const pc;
  struct feedback_props fp;
  struct cooling_function_data cooling;
  struct unit_system us;
  make_c_hyp_speed_test_engine(&e, &cosmo, &pc, &fp, &cooling, &us,
                               /*fixed_fraction=*/0.f,
                               /*timestep_off_for_debugging=*/0);

  const float hs[3] = {0.5f, 1.0f, 2.3f};
  const timebin_t bins[3] = {1, 3, 6};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      struct part p;
      set_c_hyp_test_part(&p, hs[i], bins[j]);
      radiation_snapshot_part_propagation(&p, &e);

      const float dt_phys = (float)get_timestep(bins[j], e.time_base);
      const float h_phys = (float)cosmo.a * hs[i];
      float expected = fp.ISRF_c_hyp_margin * h_phys / dt_phys;
      expected = min(expected, (float)pc.const_speed_light_c);

      assert_close_c_hyp("fraction off: c_hyp vs shipped formula",
                         p.feedback_data.c_hyp, expected);
    }
  }
  message(
      "c_hyp fixed fraction off: matches the shipped formula "
      "at every (h, time_bin) tried.");
}

/**
 * @brief Fraction on: every particle must get exactly f*c, regardless of h
 * or time bin.
 */
static void test_c_hyp_fixed_fraction_on_gives_exact_fraction_of_c(void) {
  const float f = 0.02f;

  struct engine e;
  struct cosmology cosmo;
  struct phys_const pc;
  struct feedback_props fp;
  struct cooling_function_data cooling;
  struct unit_system us;
  make_c_hyp_speed_test_engine(&e, &cosmo, &pc, &fp, &cooling, &us,
                               /*fixed_fraction=*/f,
                               /*timestep_off_for_debugging=*/0);

  const float expected = f * (float)pc.const_speed_light_c;
  const float hs[3] = {0.5f, 1.0f, 2.3f};
  const timebin_t bins[3] = {1, 3, 6};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      struct part p;
      set_c_hyp_test_part(&p, hs[i], bins[j]);
      radiation_snapshot_part_propagation(&p, &e);

      if (p.feedback_data.c_hyp != expected)
        error(
            "fraction on, h=%g bin=%d: c_hyp=%.8e, expected f*c=%.8e "
            "exactly, independent of h/time_bin.",
            (double)hs[i], bins[j], (double)p.feedback_data.c_hyp,
            (double)expected);
    }
  }
  message(
      "c_hyp fixed fraction on: exactly f*c=%.6e at every (h, time_bin) "
      "tried.",
      (double)expected);
}

/**
 * @brief The radiation timestep term: C_hyp*h/(f*c) on an eligible particle
 * with the fraction on and the debug off-switch clear; FLT_MAX (no
 * constraint) whenever the fraction is off, whenever the debug off-switch
 * is set, and on a particle outside the eligible set.
 */
static void test_c_hyp_fixed_fraction_timestep_term(void) {
  const float f = 0.02f;
  const float h = 1.7f;

  /* Fraction off: FLT_MAX regardless of eligibility. */
  {
    struct engine e;
    struct cosmology cosmo;
    struct phys_const pc;
    struct feedback_props fp;
    struct cooling_function_data cooling;
    struct unit_system us;
    make_c_hyp_speed_test_engine(&e, &cosmo, &pc, &fp, &cooling, &us,
                                 /*fixed_fraction=*/0.f,
                                 /*timestep_off_for_debugging=*/0);
    struct part p;
    set_c_hyp_test_part(&p, h, /*time_bin=*/3);
    p.feedback_data.is_illuminated_ISRF = 1;
    const float dt_rad = radiation_isrf_part_timestep(&p, &e);
    if (dt_rad != FLT_MAX)
      error(
          "fraction off: radiation_isrf_part_timestep=%.8e, expected "
          "FLT_MAX (not applied).",
          (double)dt_rad);
  }

  /* Fraction on, debug off-switch clear, eligible particle (tagged
   * illuminated): must return the receiver-side CFL bound exactly. */
  {
    struct engine e;
    struct cosmology cosmo;
    struct phys_const pc;
    struct feedback_props fp;
    struct cooling_function_data cooling;
    struct unit_system us;
    make_c_hyp_speed_test_engine(&e, &cosmo, &pc, &fp, &cooling, &us,
                                 /*fixed_fraction=*/f,
                                 /*timestep_off_for_debugging=*/0);
    struct part p;
    set_c_hyp_test_part(&p, h, /*time_bin=*/3);
    p.feedback_data.is_illuminated_ISRF = 1;

    const float c_M = f * (float)pc.const_speed_light_c;
    const float expected = fp.ISRF_c_hyp_margin * (float)cosmo.a * h / c_M;
    const float dt_rad = radiation_isrf_part_timestep(&p, &e);
    assert_close_c_hyp("fraction on, eligible: dt_rad vs C_hyp*h/(f*c)", dt_rad,
                       expected);

    /* Same fixture, debug off-switch set: must go back to FLT_MAX. */
    fp.ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging = 1;
    const float dt_rad_off = radiation_isrf_part_timestep(&p, &e);
    if (dt_rad_off != FLT_MAX)
      error(
          "fraction on, timestep term off for debugging: "
          "radiation_isrf_part_timestep=%.8e, expected FLT_MAX.",
          (double)dt_rad_off);
    fp.ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging = 0;

    /* Same fixture, particle outside the eligible set (never illuminated,
     * no illuminated neighbour in kernel): must also be FLT_MAX. */
    struct part p_far;
    set_c_hyp_test_part(&p_far, h, /*time_bin=*/3);
    const float dt_rad_far = radiation_isrf_part_timestep(&p_far, &e);
    if (dt_rad_far != FLT_MAX)
      error(
          "fraction on, not near field: radiation_isrf_part_timestep="
          "%.8e, expected FLT_MAX.",
          (double)dt_rad_far);
  }
  message(
      "radiation timestep term: C_hyp*h/(f*c) on an eligible particle, "
      "FLT_MAX with the fraction off, with the debug off-switch set, and "
      "off the eligible set.");
}

/**
 * @brief Run @p scheme/@p fixed_fraction through
 * #feedback_props_check_c_hyp_scheme() in a forked child and assert it
 * aborts via error() (exit 1, or SIGABRT under SWIFT_DEVELOP_MODE), per
 * the #testRadiationRebuildCheck.c:365 pattern.
 *
 * @param label Message label for a failure report.
 * @param scheme Value to pass as ISRF_c_hyp_scheme.
 * @param fixed_fraction Value to pass as ISRF_c_hyp_fixed_fraction_of_c.
 */
static void assert_c_hyp_scheme_check_rejects(const char *label, int scheme,
                                              float fixed_fraction) {
  const pid_t pid = fork();
  if (pid == 0) {
    /* Child process: silence stderr (the error() message is expected
     * output, not a test failure), then trigger the check. */
    if (freopen("/dev/null", "w", stderr) == NULL) _exit(43);
    feedback_props_check_c_hyp_scheme(scheme, fixed_fraction);
    /* Reached only if the check did NOT reject -- signal failure with a
     * distinguishable exit code (real error() exits with status 1). */
    _exit(42);
  } else if (pid > 0) {
    int status;
    waitpid(pid, &status, 0);
    const int exited_with_error = WIFEXITED(status) && WEXITSTATUS(status) == 1;
    const int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
    if (!exited_with_error && !aborted)
      error(
          "%s: feedback_props_check_c_hyp_scheme(scheme=%d, "
          "fixed_fraction=%g) failed to reject (WIFEXITED=%d "
          "WEXITSTATUS=%d WIFSIGNALED=%d WTERMSIG=%d).",
          label, scheme, (double)fixed_fraction, WIFEXITED(status),
          WIFEXITED(status) ? WEXITSTATUS(status) : -1, WIFSIGNALED(status),
          WIFSIGNALED(status) ? WTERMSIG(status) : -1);
  } else {
    error("fork() failed in the c_hyp scheme/fraction mismatch test.");
  }
}

/**
 * @brief The two ISRF c_hyp speed schemes are alternatives, not layers
 * (requirement 4 of the comparison branch): every mismatched
 * (ISRF_c_hyp_scheme, ISRF_c_hyp_fixed_fraction_of_c) pairing must be
 * rejected by #feedback_props_check_c_hyp_scheme() at parse time, and
 * every matched pairing must return normally.
 */
static void test_c_hyp_scheme_mismatch_rejected(void) {

  /* Invalid: fixed-fraction scheme selected with no magnitude. */
  assert_c_hyp_scheme_check_rejects("scheme=fixed_fraction, fraction=0",
                                    isrf_c_hyp_scheme_fixed_fraction, 0.f);
  /* Invalid: a magnitude set while the kernel-local scheme is selected --
   * the two schemes stacked instead of chosen between. */
  assert_c_hyp_scheme_check_rejects("scheme=kernel_local, fraction=0.02",
                                    isrf_c_hyp_scheme_kernel_local, 0.02f);
  /* Invalid: a magnitude set while the shipped scheme is selected -- would
   * otherwise silently do nothing (see radiation_snapshot_part_propagation:
   * the fraction is only read when the scheme selects it). */
  assert_c_hyp_scheme_check_rejects("scheme=shipped, fraction=0.02",
                                    isrf_c_hyp_scheme_shipped, 0.02f);
  /* Invalid: a magnitude set while the combined kernel-local +
   * variable-c scheme is selected -- same stacking mistake as
   * scheme=kernel_local above. */
  assert_c_hyp_scheme_check_rejects(
      "scheme=kernel_local_plus_variable_c, fraction=0.02",
      isrf_c_hyp_scheme_kernel_local_plus_variable_c, 0.02f);

  /* Valid pairings must return normally (no fork needed: nothing to
   * observe but the absence of an abort). */
  feedback_props_check_c_hyp_scheme(isrf_c_hyp_scheme_shipped, 0.f);
  feedback_props_check_c_hyp_scheme(isrf_c_hyp_scheme_kernel_local, 0.f);
  feedback_props_check_c_hyp_scheme(isrf_c_hyp_scheme_fixed_fraction, 0.02f);
  feedback_props_check_c_hyp_scheme(isrf_c_hyp_scheme_consistent_variable_c,
                                    0.f);
  feedback_props_check_c_hyp_scheme(
      isrf_c_hyp_scheme_kernel_local_plus_variable_c, 0.f);

  message(
      "c_hyp scheme/fraction mismatch: every invalid pairing rejected at "
      "parse time, every valid pairing accepted.");
}

/**
 * @brief #radiation_end_density_propagation must be a no-op for the three
 * schemes it does not own (#isrf_c_hyp_scheme_shipped,
 * #isrf_c_hyp_scheme_fixed_fraction, #isrf_c_hyp_scheme_consistent_variable_c
 * -- the last keeps the shipped `dt_i` speed formula, only the operators
 * change): drift-time #radiation_snapshot_part_propagation already decided
 * #c_hyp for them, and this is the one guard standing between the default
 * scheme and silently losing bit-identity with the pre-comparison-branch
 * behaviour. Sets #max_ngb_time_bin to a bin the kernel-local formula would
 * definitely act on if reached, so a regression that dropped the scheme
 * gate (or one that accidentally widened it to include scheme 3) would flip
 * this test. The owning set (schemes 1 and 4) is covered separately by
 * #test_kernel_local_c_hyp and
 * #test_kernel_local_plus_variable_c_composition.
 */
static void test_c_hyp_end_density_no_clobber_for_other_schemes(void) {

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_speed_light_c = 3e5;

  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  fb_props.ISRF_propagation = 1;
  fb_props.ISRF_c_hyp_margin = 0.5f;

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.feedback_props = &fb_props;
  e.cosmology = &cosmo;
  e.physical_constants = &phys_const;
  e.time_base = 1e-3;
  e.ti_current = 8;
  e.policy = 0;

  const int schemes[3] = {isrf_c_hyp_scheme_shipped,
                          isrf_c_hyp_scheme_fixed_fraction,
                          isrf_c_hyp_scheme_consistent_variable_c};
  const float sentinel = 987.654f;

  for (int s = 0; s < 3; s++) {
    fb_props.ISRF_c_hyp_scheme = schemes[s];
    fb_props.ISRF_c_hyp_fixed_fraction_of_c =
        schemes[s] == isrf_c_hyp_scheme_fixed_fraction ? 0.02f : 0.f;

    struct part p;
    bzero(&p, sizeof(struct part));
    p.h = 0.37f;
    p.time_bin = 6;
    p.feedback_data.dt_prev = (float)get_timestep(p.time_bin, e.time_base);
    /* Different from time_bin: the kernel-local formula would certainly
     * move c_hyp if this scheme gate were ever bypassed. */
    p.feedback_data.max_ngb_time_bin = p.time_bin + 4;
    p.feedback_data.c_hyp = sentinel;

    radiation_end_density_propagation(&p, &e);

    if (!bits_equal_f(p.feedback_data.c_hyp, sentinel))
      error(
          "scheme=%d: radiation_end_density_propagation must not touch "
          "c_hyp (got %.9g, expected the untouched sentinel %.9g)",
          schemes[s], (double)p.feedback_data.c_hyp, (double)sentinel);
  }
  message(
      "c_hyp end-density no-clobber: radiation_end_density_propagation is "
      "a no-op for the shipped, fixed-fraction and consistent-variable-c "
      "schemes.");
}

/* ---------------------------------------------------------------------
 * isrf_c_hyp_scheme_consistent_variable_c (closing FABLE review,
 * .claude/dev/ISRF_HISTORY.md 2026-09-18 ~06:00, section C): every
 * operator at particle i becomes c_hyp_i/c times the true-speed equation,
 * with the stored flux becoming the reduced flux Ft = F_true/c_hyp. Three
 * decisive unit tests below.
 * ------------------------------------------------------------------- */

/**
 * @brief Write one band's pairwise flux input for the consistent-variable-c
 * bit-identity test: the true physical flux F_true when @p reduced is 0,
 * or its reduced-flux representation Ft = F_true/c_hyp when @p reduced is
 * 1 (this scheme's own stored state; see radiation_propagation_iact.h's
 * file header). @p c_hyp must be an exact power of two for the division to
 * be exact (see #test_consistent_variable_c_uniform_c_bit_identical's own
 * doxygen for why that matters).
 *
 * @param band (return) The band to write.
 * @param F_true The physical flux vector.
 * @param c_hyp The particle's own #feedback_part_data.c_hyp.
 * @param reduced 0 to store F_true directly, 1 to store F_true/c_hyp.
 */
static void set_consistent_c_test_band(struct feedback_isrf_band_data *band,
                                       const float F_true[3], float c_hyp,
                                       int reduced) {
  const float scale = reduced ? 1.f / c_hyp : 1.f;
  band->specific_flux[0] = F_true[0] * scale;
  band->specific_flux[1] = F_true[1] * scale;
  band->specific_flux[2] = F_true[2] * scale;
}

/**
 * @brief THE DECISIVE UNIT TEST for isrf_c_hyp_scheme_consistent_variable_c:
 * with a UNIFORM c_hyp, the new scheme's force-loop divergence and
 * dissipation accumulators must reduce BIT FOR BIT to the shipped scheme's,
 * through the real #runner_iact_isrf_dissipation dispatch. Checked on raw
 * float bits (#bits_equal_f), not a tolerance.
 *
 * c_hyp is chosen as an exact power of two (8.0): under this scheme,
 * #feedback_isrf_band_data.specific_flux stores Ft = F_true/c_hyp instead
 * of F_true, and IEEE-754 multiplication/division by an exact power of two
 * moves only the exponent field, leaving the significand untouched; a
 * later multiplication by the same power of two therefore reproduces the
 * pre-division rounding decision exactly (no new rounding is introduced by
 * the round trip, and scaling commutes exactly with rounding elsewhere in
 * the expression). An arbitrary c_hyp would only be numerically close, not
 * bit-identical, which is not what this test requires (see
 * swift-knowledge.md's -ffast-math notes: source-level operation order,
 * not just mathematical equivalence, decides bit-exactness). The
 * dissipation accumulator needs no such care: with c_i = c_j the shipped
 * `alpha_ij*min(c_i,c_j)` and this scheme's `alpha_ij*c_i`/`alpha_ij*c_j`
 * are literally the same expression evaluated on the same operands.
 */
static void test_consistent_variable_c_uniform_c_bit_identical(void) {

  const float c_hyp = 8.f; /* exact power of two */
  const float hi = 0.6f, hj = 0.9f;
  const float mi = 1.3f, mj = 0.7f;
  const float rho_i = 1.1f, rho_j = 0.6f;
  const float dx[3] = {0.31f, -0.12f, 0.05f};
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  const float F_true[ISRF_BAND_COUNT][2][3] = {
      {{0.4f, -0.2f, 0.1f}, {-0.15f, 0.25f, -0.05f}},
      {{0.05f, 0.6f, -0.3f}, {0.2f, -0.1f, 0.4f}}};

  struct part pi_A, pj_A, pi_B, pj_B;
  bzero(&pi_A, sizeof(struct part));
  bzero(&pj_A, sizeof(struct part));
  bzero(&pi_B, sizeof(struct part));
  bzero(&pj_B, sizeof(struct part));

  pi_A.h = pi_B.h = hi;
  pj_A.h = pj_B.h = hj;
  pi_A.mass = pi_B.mass = mi;
  pj_A.mass = pj_B.mass = mj;

  struct part *const parts_A[2] = {&pi_A, &pj_A};
  struct part *const parts_B[2] = {&pi_B, &pj_B};
  const float rhos[2] = {rho_i, rho_j};
  for (int s = 0; s < 2; s++) {
    parts_A[s]->feedback_data.rho_prev = rhos[s];
    parts_A[s]->feedback_data.c_hyp = c_hyp;
    parts_B[s]->feedback_data.rho_prev = rhos[s];
    parts_B[s]->feedback_data.c_hyp = c_hyp;
    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      const float u = 0.5f + 0.1f * s + 0.2f * b;
      const float alpha_trigger = 0.1f + 0.05f * s;
      parts_A[s]->feedback_data.isrf_band[b].u = u;
      parts_B[s]->feedback_data.isrf_band[b].u = u;
      parts_A[s]->feedback_data.isrf_band[b].dissipation_alpha_trigger =
          alpha_trigger;
      parts_B[s]->feedback_data.isrf_band[b].dissipation_alpha_trigger =
          alpha_trigger;
      parts_A[s]->feedback_data.isrf_band[b].dissipation_alpha_floor = 0.05f;
      parts_B[s]->feedback_data.isrf_band[b].dissipation_alpha_floor = 0.05f;
      set_consistent_c_test_band(&parts_A[s]->feedback_data.isrf_band[b],
                                 F_true[b][s], c_hyp, /*reduced=*/0);
      set_consistent_c_test_band(&parts_B[s]->feedback_data.isrf_band[b],
                                 F_true[b][s], c_hyp, /*reduced=*/1);
    }
  }

  isrf_c_hyp_consistent_variable_c = 0;
  runner_iact_isrf_dissipation(r2, dx, hi, hj, &pi_A, &pj_A, 1.f, 0.f);

  isrf_c_hyp_consistent_variable_c = 1;
  runner_iact_isrf_dissipation(r2, dx, hi, hj, &pi_B, &pj_B, 1.f, 0.f);
  isrf_c_hyp_consistent_variable_c = 0; /* restore the default */

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    const struct feedback_isrf_band_data *ai = &pi_A.feedback_data.isrf_band[b];
    const struct feedback_isrf_band_data *aj = &pj_A.feedback_data.isrf_band[b];
    const struct feedback_isrf_band_data *bi = &pi_B.feedback_data.isrf_band[b];
    const struct feedback_isrf_band_data *bj = &pj_B.feedback_data.isrf_band[b];

    if (!bits_equal_f(ai->div_specific_flux, bi->div_specific_flux))
      error(
          "uniform-c bit identity: band %d div_F_i shipped=%.9g "
          "consistent=%.9g",
          b, (double)ai->div_specific_flux, (double)bi->div_specific_flux);
    if (!bits_equal_f(aj->div_specific_flux, bj->div_specific_flux))
      error(
          "uniform-c bit identity: band %d div_F_j shipped=%.9g "
          "consistent=%.9g",
          b, (double)aj->div_specific_flux, (double)bj->div_specific_flux);
    if (!bits_equal_f(ai->dissipation_u, bi->dissipation_u))
      error(
          "uniform-c bit identity: band %d dissipation_u_i shipped=%.9g "
          "consistent=%.9g",
          b, (double)ai->dissipation_u, (double)bi->dissipation_u);
    if (!bits_equal_f(aj->dissipation_u, bj->dissipation_u))
      error(
          "uniform-c bit identity: band %d dissipation_u_j shipped=%.9g "
          "consistent=%.9g",
          b, (double)aj->dissipation_u, (double)bj->dissipation_u);
  }
  message(
      "uniform-c bit identity: consistent-variable-c operators reduce bit "
      "for bit to the shipped scheme's at uniform c_hyp");
}

/**
 * @brief Requirement 3 (closing review, section (C)): the receiver's
 * Courant number at a time-bin seam must not depend on the SENDER's bin at
 * all. c_hyp_i under isrf_c_hyp_scheme_consistent_variable_c is
 * #isrf_c_hyp_scheme_shipped's own formula (C_hyp*h_i/dt_i), which by
 * construction reads only particle i's own #part.time_bin/#part.h, never a
 * neighbour's, so `c_i*dt_i/h_j = C_hyp*(h_i/h_j)` holds for every sender
 * h_j and time bin, dt-free. Verified here through the real
 * #radiation_snapshot_part_propagation dispatch (not by re-deriving the
 * formula) for bin differences 1, 2, 3 crossed with h ratios 0.5, 1, 2 (9
 * combinations): the receiver's own c_i/dt_i must be UNCHANGED across
 * every combination (proving the "no dependence on the sender" half), and
 * the algebraic identity must hold in each (the seam-Courant half).
 *
 * Naming note: the task/review phrase this bound "C_hyp*h_j/h_i"; the
 * derivation actually checked here gives C_hyp*h_i/h_j (receiver over
 * sender), which is the dt-free form (c_i*dt_i/h_j =
 * (C_hyp*h_i/dt_i)*dt_i/h_j = C_hyp*h_i/h_j). The two differ only by which
 * particle's h sits in the numerator; this test is unambiguous about which
 * one (i) is the receiver, and reports the mismatch here rather than
 * silently relabeling the prose to match.
 */
static void test_consistent_variable_c_seam_courant_number(void) {

  struct engine e;
  struct cosmology cosmo;
  struct phys_const pc;
  struct feedback_props fp;
  struct cooling_function_data cooling;
  struct unit_system us;
  make_c_hyp_speed_test_engine(&e, &cosmo, &pc, &fp, &cooling, &us,
                               /*fixed_fraction=*/0.f,
                               /*timestep_off_for_debugging=*/0);
  /* Same c_hyp formula as shipped (see this test's own doxygen): only the
   * scheme selector itself needs overriding after the helper above. */
  fp.ISRF_c_hyp_scheme = isrf_c_hyp_scheme_consistent_variable_c;

  const float h_i = 0.42f;
  const timebin_t bin_i = 4;
  const int bin_diffs[3] = {1, 2, 3};
  const float h_ratios[3] = {0.5f, 1.f, 2.f};

  float c_i_ref = -1.f, dt_i_ref = -1.f;

  for (int bd = 0; bd < 3; bd++) {
    for (int hr = 0; hr < 3; hr++) {
      struct part p;
      set_c_hyp_test_part(&p, h_i, bin_i);
      radiation_snapshot_part_propagation(&p, &e);

      const float c_i = p.feedback_data.c_hyp;
      const float dt_i = p.feedback_data.dt_prev;
      if (c_i_ref < 0.f) {
        c_i_ref = c_i;
        dt_i_ref = dt_i;
      } else if (!bits_equal_f(c_i, c_i_ref) || !bits_equal_f(dt_i, dt_i_ref)) {
        error(
            "seam Courant number: receiver c_hyp/dt changed with the "
            "sender's bin difference=%d h_ratio=%.2f (c_i=%.9g dt_i=%.9g, "
            "expected %.9g/%.9g): the formula must not read the sender.",
            bin_diffs[bd], (double)h_ratios[hr], (double)c_i, (double)dt_i,
            (double)c_i_ref, (double)dt_i_ref);
      }

      const timebin_t bin_j = bin_i + bin_diffs[bd];
      (void)bin_j; /* only used to name the sender's bin in messages/errors */
      const float h_j = h_i * h_ratios[hr];

      const double lhs = (double)c_i * (double)dt_i / (double)h_j;
      const double rhs =
          (double)fp.ISRF_c_hyp_margin * (double)h_i / (double)h_j;
      if (fabs(lhs - rhs) > 1e-6 * fabs(rhs))
        error(
            "seam Courant number: bin_diff=%d h_ratio=%.2f: c_i*dt_i/h_j = "
            "%.9e != C_hyp*h_i/h_j = %.9e",
            bin_diffs[bd], (double)h_ratios[hr], lhs, rhs);
    }
  }
  message(
      "seam Courant number: c_i*dt_i/h_j = C_hyp*h_i/h_j for every bin "
      "difference (1..3) and h ratio (0.5/1/2), independent of the "
      "sender's own bin");
}

/**
 * @brief Requirement 4 (closing review, section (C)): the mass-weighted
 * exchange every other scheme conserves exactly, `m_i*X_i + m_j*X_j = 0`,
 * is no longer the invariant once each side carries its OWN c_hyp (see
 * radiation_propagation_iact.h's file header for the derivation);
 * `m_i*X_i/c_i + m_j*X_j/c_j = 0` is, for both the divergence and the
 * dissipation accumulators. Checked here at genuinely different c_i != c_j
 * (unlike the bit-identity test above, which deliberately keeps c_i = c_j)
 * across three trials, through the real #runner_iact_isrf_dissipation
 * dispatch.
 */
static void test_consistent_variable_c_conservation(void) {

  const float ci_values[3] = {4.f, 9.5f, 22.f};
  const float cj_values[3] = {12.f, 5.5f, 6.25f};
  const double SEAM_BAR = 1e-5;

  isrf_c_hyp_consistent_variable_c = 1;

  for (int t = 0; t < 3; t++) {
    struct part pi, pj;
    bzero(&pi, sizeof(struct part));
    bzero(&pj, sizeof(struct part));

    pi.h = 0.5f + 0.1f * t;
    pj.h = 0.8f - 0.05f * t;
    pi.mass = 1.1f + 0.3f * t;
    pj.mass = 0.6f + 0.2f * t;
    pi.feedback_data.rho_prev = 1.0f + 0.2f * t;
    pj.feedback_data.rho_prev = 0.7f + 0.1f * t;
    pi.feedback_data.c_hyp = ci_values[t];
    pj.feedback_data.c_hyp = cj_values[t];

    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      struct feedback_isrf_band_data *bi = &pi.feedback_data.isrf_band[b];
      struct feedback_isrf_band_data *bj = &pj.feedback_data.isrf_band[b];
      bi->specific_flux[0] = 0.3f + 0.05f * t + 0.1f * b;
      bi->specific_flux[1] = -0.2f + 0.02f * t;
      bi->specific_flux[2] = 0.15f;
      bj->specific_flux[0] = -0.1f + 0.03f * t;
      bj->specific_flux[1] = 0.25f - 0.01f * b;
      bj->specific_flux[2] = -0.05f;
      bi->u = 0.6f + 0.1f * t;
      bj->u = 0.4f + 0.05f * b;
      bi->dissipation_alpha_trigger = 0.2f;
      bj->dissipation_alpha_trigger = 0.15f;
      bi->dissipation_alpha_floor = 0.05f;
      bj->dissipation_alpha_floor = 0.05f;
    }

    const float dx[3] = {0.2f + 0.05f * t, -0.1f, 0.05f};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    runner_iact_isrf_dissipation(r2, dx, pi.h, pj.h, &pi, &pj, 1.f, 0.f);

    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      const struct feedback_isrf_band_data *bi = &pi.feedback_data.isrf_band[b];
      const struct feedback_isrf_band_data *bj = &pj.feedback_data.isrf_band[b];

      const double div_i = (double)pi.mass * (double)bi->div_specific_flux /
                           (double)pi.feedback_data.c_hyp;
      const double div_j = (double)pj.mass * (double)bj->div_specific_flux /
                           (double)pj.feedback_data.c_hyp;
      const double div_ratio =
          fabs(div_i + div_j) / max(fabs(div_i) + fabs(div_j), 1e-30);
      if (!(div_ratio <= SEAM_BAR))
        error(
            "conservation (divergence/c): trial %d band %d: "
            "m_i*div_i/c_i + m_j*div_j/c_j ratio = %.3e above %.3e",
            t, b, div_ratio, SEAM_BAR);

      const double diss_i = (double)pi.mass * (double)bi->dissipation_u /
                            (double)pi.feedback_data.c_hyp;
      const double diss_j = (double)pj.mass * (double)bj->dissipation_u /
                            (double)pj.feedback_data.c_hyp;
      const double diss_ratio =
          fabs(diss_i + diss_j) / max(fabs(diss_i) + fabs(diss_j), 1e-30);
      if (!(diss_ratio <= SEAM_BAR))
        error(
            "conservation (dissipation/c): trial %d band %d: "
            "m_i*diss_i/c_i + m_j*diss_j/c_j ratio = %.3e above %.3e",
            t, b, diss_ratio, SEAM_BAR);
    }
  }

  isrf_c_hyp_consistent_variable_c = 0; /* restore the default */
  message(
      "conservation under variable c: m_i*X_i/c_i + m_j*X_j/c_j = 0 for "
      "both the divergence and dissipation accumulators, at genuinely "
      "different c_i != c_j");
}

/* ---------------------------------------------------------------------
 * isrf_c_hyp_scheme_kernel_local_plus_variable_c: the two candidates above
 * fix different defects (kernel-local narrows the speed contrast between
 * neighbours; consistent-variable-c fixes the pairwise operators'
 * amplitude error) via disjoint gates (ISRF_c_hyp_scheme ==
 * isrf_c_hyp_scheme_kernel_local for the speed formula,
 * isrf_c_hyp_consistent_variable_c for the operators), so scheme 4 selects
 * both without a combined re-derivation. The decisive test below exercises
 * the ONE place they meet numerically: #radiation_cache_m1_closure_part's
 * `c_M`, pinned to 1 under the change of variable regardless of `c_hyp`'s
 * own value (see that function's own doxygen) -- so under scheme 4,
 * #radiation_end_density_propagation must still move `c_hyp` itself (the
 * kernel-local half, read by the pairwise operators' receiver-side
 * multiply), while leaving the just-rebuilt M1 closure numerically
 * unaffected by that move (the change-of-variable half). A regression that
 * made the two interact (e.g. `c_M` accidentally reading the NEW `c_hyp`)
 * would move `m1_closure_D` here and flip this test.
 * ------------------------------------------------------------------- */

/**
 * @brief THE DECISIVE UNIT TEST for
 * #isrf_c_hyp_scheme_kernel_local_plus_variable_c: through the real
 * #radiation_end_density_propagation dispatch, with
 * #isrf_c_hyp_consistent_variable_c set (as feedback_props_init() would set
 * it for this scheme), (1) `c_hyp` moves from a same-bin sentinel to the
 * kernel-local `dt_max(i)` formula's value -- proving the speed axis is
 * still engaged under scheme 4, exactly as it is under scheme 1 -- and (2)
 * the M1 closure tensor rebuilt in the same call is bit-identical before
 * and after, since #isrf_c_hyp_consistent_variable_c pins `c_M = 1`
 * independent of `c_hyp` -- proving the two axes do not interact through
 * this shared call site.
 */
static void test_kernel_local_plus_variable_c_composition(void) {

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_speed_light_c = 3e5;

  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  fb_props.ISRF_propagation = 1;
  fb_props.ISRF_c_hyp_scheme = isrf_c_hyp_scheme_kernel_local_plus_variable_c;
  fb_props.ISRF_c_hyp_margin = 0.5f;

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.feedback_props = &fb_props;
  e.cosmology = &cosmo;
  e.physical_constants = &phys_const;
  e.time_base = 1e-3;
  e.ti_current = 8;
  e.policy = 0;

  /* feedback_props_init() sets this global from ISRF_c_hyp_scheme; set it
   * directly here since this unit test bypasses the parser. */
  isrf_c_hyp_consistent_variable_c = 1;

  struct part p;
  bzero(&p, sizeof(struct part));
  p.h = 0.37f;
  p.time_bin = 6;
  p.feedback_data.dt_prev = (float)get_timestep(p.time_bin, e.time_base);
  /* Different from time_bin, exactly as the no-clobber test's sentinel is:
   * the kernel-local formula must move c_hyp away from the same-bin value
   * for this test to be decisive. */
  p.feedback_data.max_ngb_time_bin = p.time_bin + 4;
  const float sentinel_c_hyp = 42.f;
  p.feedback_data.c_hyp = sentinel_c_hyp;
  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    p.feedback_data.isrf_band[b].u = 0.6f + 0.1f * b;
    p.feedback_data.isrf_band[b].specific_flux[0] = 0.2f + 0.05f * b;
    p.feedback_data.isrf_band[b].specific_flux[1] = -0.1f;
    p.feedback_data.isrf_band[b].specific_flux[2] = 0.05f;
  }

  radiation_end_density_propagation(&p, &e);

  /* (1) Speed axis: c_hyp must equal the kernel-local dt_max(i) formula,
   * moved away from the same-bin sentinel -- same computation scheme 1's
   * own test_c_hyp_same_bin_cap_and_pin exercises directly (no cosmology
   * here, so dt_max(i) is the plain get_timestep at the neighbour-maximum
   * bin, matching radiation_end_density_propagation's own no-cosmology
   * branch). */
  float dt_max =
      (float)get_timestep(p.feedback_data.max_ngb_time_bin, e.time_base);
  dt_max = max(dt_max, FLT_MIN);
  const float h_phys = (float)e.cosmology->a * p.h;
  float expected_c_hyp = e.feedback_props->ISRF_c_hyp_margin * h_phys / dt_max;
  expected_c_hyp =
      min(expected_c_hyp, (float)e.physical_constants->const_speed_light_c);

  if (bits_equal_f(p.feedback_data.c_hyp, sentinel_c_hyp))
    error(
        "composition: c_hyp did not move off the same-bin sentinel under "
        "scheme 4 -- the kernel-local speed axis did not engage.");
  if (!bits_equal_f(p.feedback_data.c_hyp, expected_c_hyp))
    error("composition: c_hyp = %.9g != kernel-local dt_max(i) formula %.9g",
          (double)p.feedback_data.c_hyp, (double)expected_c_hyp);

  /* (2) Operator axis: the M1 closure rebuilt by the same call must match
   * one built directly with c_M = 1 (isrf_c_hyp_consistent_variable_c's own
   * formula), independent of the c_hyp value just computed above -- proving
   * the closure did not pick up the NEW c_hyp through this call.
   *
   * Compared to a tolerance, not on raw bits: the reference is the same
   * function re-run in a different inlining context, so on an FMA-capable
   * target the two contract differently and land a few ULP apart. The bar
   * sits ~30x above that round-off and ~1e5 below the deviation a closure
   * built on c_hyp instead of c_M = 1 produces (a factor of two). The
   * absolute floor covers the off-diagonal components, which sit near zero
   * where a purely relative bar carries no information. Finiteness is
   * tested first: a non-finite value compares false against every bar and
   * would slip through the tolerance gate silently. */
  const double closure_rel_bar = 1e-5;
  const double closure_abs_bar = 1e-8;
  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    float expected_D[3][3];
    radiation_get_m1_closure_tensor_band(
        p.feedback_data.isrf_band[b].u,
        p.feedback_data.isrf_band[b].specific_flux,
        /*c_M=*/1.f, expected_D);
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++) {
        const double got =
            (double)p.feedback_data.isrf_band[b].m1_closure_D[r][c];
        const double ref = (double)expected_D[r][c];
        if (!is_finite_bits(got) || !is_finite_bits(ref))
          error(
              "composition: band %d m1_closure_D[%d][%d] = %.9g against the "
              "c_M=1 closure %.9g: a non-finite closure component.",
              b, r, c, got, ref);
        const double dev = fabs(got - ref);
        const double bar = closure_rel_bar * fabs(ref) + closure_abs_bar;
        if (!(dev <= bar))
          error(
              "composition: band %d m1_closure_D[%d][%d] = %.9g != the "
              "c_M=1 closure %.9g (deviation %.3e above %.3e) -- the "
              "operator axis picked up the NEW c_hyp instead of staying "
              "pinned at c_M=1.",
              b, r, c, got, ref, dev, bar);
      }
  }

  isrf_c_hyp_consistent_variable_c = 0; /* restore the default */
  message(
      "composition: under scheme 4, c_hyp moves via the kernel-local "
      "dt_max(i) formula (speed axis engaged) while the M1 closure stays "
      "pinned at c_M=1, independent of that move (operator axis "
      "unaffected) -- the two do not interact through this shared call "
      "site.");
}

/**
 * @brief The default ISRF scheme (#isrf_c_hyp_scheme_shipped == 0) must be
 * unchanged by adding #isrf_c_hyp_scheme_kernel_local_plus_variable_c: the
 * enum ordinal stays 0, the parser's own default argument
 * (#feedback_props_init(), unchanged by this diff) still resolves to it,
 * and #isrf_c_hyp_consistent_variable_c must stay unset for it.
 */
static void test_default_scheme_unchanged(void) {

  if (isrf_c_hyp_scheme_shipped != 0)
    error(
        "default scheme regression: isrf_c_hyp_scheme_shipped = %d, "
        "expected 0 -- adding scheme 4 must not renumber the existing "
        "values.",
        (int)isrf_c_hyp_scheme_shipped);

  /* A zero-initialized feedback_props (the struct's own default before
   * parsing, and what every other test in this file relies on) selects the
   * shipped scheme and must not flag the change-of-variable operators. */
  struct feedback_props fb_props;
  bzero(&fb_props, sizeof(struct feedback_props));
  if (fb_props.ISRF_c_hyp_scheme != isrf_c_hyp_scheme_shipped)
    error(
        "default scheme regression: a zero-initialized feedback_props "
        "does not resolve to isrf_c_hyp_scheme_shipped.");

  isrf_c_hyp_consistent_variable_c =
      (fb_props.ISRF_c_hyp_scheme == isrf_c_hyp_scheme_consistent_variable_c ||
       fb_props.ISRF_c_hyp_scheme ==
           isrf_c_hyp_scheme_kernel_local_plus_variable_c);
  if (isrf_c_hyp_consistent_variable_c)
    error(
        "default scheme regression: the default scheme must not set "
        "isrf_c_hyp_consistent_variable_c.");
  isrf_c_hyp_consistent_variable_c = 0; /* restore the default */

  message(
      "default scheme unchanged: isrf_c_hyp_scheme_shipped is still 0 and "
      "a zero-initialized feedback_props still resolves to it, with the "
      "change-of-variable operators left off.");
}

/**
 * @brief Run one geometry: small-h cell at the origin, large-h cell at the
 * offset (or the reverse), swept without depth limits, then across two
 * depth levels, then with the large-h level dropped.
 *
 * @param r The runner.
 * @param offset Offset of the second cell.
 * @param small_first 1 to put the small-h cell at the origin.
 * @param part_id (in/out) Running particle id.
 */
static void run_geometry(struct runner *r, const double offset[3],
                         int small_first, long long *part_id) {

  const double origin[3] = {0., 0., 0.};
  const double h_small = 1.2348;
  const double h_large = 2. * h_small;
  /* Every small h is below this and every large h above it. */
  const float h_split = 1.5 * h_small / CELL_N;

  for (int mode = 0; mode < 3; mode++) {
    srand(1234 + 7 * mode + small_first);
    long long id = *part_id;
    struct cell *small = make_cell(small_first ? origin : offset, h_small,
                                   /*depth_h=*/1, &id);
    struct cell *large = make_cell(small_first ? offset : origin, h_large,
                                   /*depth_h=*/0, &id);
    struct cell *cells[2] = {small, large};

    for (int k = 0; k < 2; k++) {
      runner_do_hydro_sort(r, cells[k], 0x1FFF, 0, 0, 0, 0);
      prepare_force(cells[k], r->e);
    }

    if (mode == 0) {
      for (int k = 0; k < 2; k++) {
        cells[k]->depth = 0;
        for (int i = 0; i < cells[k]->hydro.count; i++)
          cells[k]->hydro.parts[i].depth_h = 0;
        runner_doself2_branch_force(r, cells[k], 0, 0);
      }
      runner_dopair2_branch_force(r, small, large, 0, 0);
      check_cells(cells, "no depth limit", 1);
    } else {
      /* Deeper level: particles whose h fits below the split. */
      for (int k = 0; k < 2; k++) set_level(cells[k], 1, 0.f, h_split);
      for (int k = 0; k < 2; k++)
        runner_doself2_branch_force(r, cells[k], 0, 1);
      runner_dopair2_branch_force(r, small, large, 0, 1);

      if (mode == 1) {
        /* Parent level: particles too large to recurse. */
        for (int k = 0; k < 2; k++)
          set_level(cells[k], 0, h_split, 1.f / kernel_gamma);
        for (int k = 0; k < 2; k++)
          runner_doself2_branch_force(r, cells[k], 1, 1);
        runner_dopair2_branch_force(r, small, large, 1, 1);
        check_cells(cells, "two depth levels", 1);
      } else {
        check_cells(cells, "parent level dropped", 0);
      }
    }

    clean_up(small);
    clean_up(large);
  }
  *part_id += 2 * CELL_N * CELL_N * CELL_N;
}

int main(int argc, char *argv[]) {

  clocks_set_cpufreq(0);

#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  test_c_hyp_fixed_fraction_off_matches_shipped_formula();
  test_c_hyp_fixed_fraction_on_gives_exact_fraction_of_c();
  test_c_hyp_fixed_fraction_timestep_term();
  test_kernel_local_c_hyp();
  test_c_hyp_scheme_mismatch_rejected();
  test_c_hyp_end_density_no_clobber_for_other_schemes();
  test_consistent_variable_c_uniform_c_bit_identical();
  test_consistent_variable_c_seam_courant_number();
  test_consistent_variable_c_conservation();
  test_kernel_local_plus_variable_c_composition();
  test_default_scheme_unchanged();

  struct space space;
  struct engine engine;
  struct cosmology cosmo;
  struct hydro_props hydro_props;
  struct pressure_floor_props pressure_floor;
  struct phys_const phys_const;
  bzero(&space, sizeof(struct space));
  bzero(&engine, sizeof(struct engine));
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  bzero(&phys_const, sizeof(struct phys_const));

  space.periodic = 0;
  for (int k = 0; k < 3; k++) space.dim[k] = 3.;

  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.nodeID = NODE_ID;
  engine.time_base = 1e-3;
  phys_const.const_vacuum_permeability = 1.0;
  engine.physical_constants = &phys_const;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;
  hydro_props_init_no_hydro(&hydro_props);
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;

  struct runner *runner = NULL;
  if (posix_memalign((void **)&runner, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct runner)) != 0)
    error("Couldn't allocate the runner");
  bzero(runner, sizeof(struct runner));
  runner->e = &engine;
#ifdef WITH_VECTORIZATION
  cache_init(&runner->ci_cache, 512);
  cache_init(&runner->cj_cache, 512);
#endif

  long long part_id = 0;
  const double offsets[3][3] = {{1., 0., 0.}, {1., 1., 0.}, {1., 1., 1.}};
  for (int g = 0; g < 3; g++)
    for (int small_first = 0; small_first < 2; small_first++)
      run_geometry(runner, offsets[g], small_first, &part_id);

#ifdef WITH_VECTORIZATION
  cache_clean(&runner->ci_cache);
  cache_clean(&runner->cj_cache);
#endif
  free(runner);
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR && SPHENIX_SPH */
