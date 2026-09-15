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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
