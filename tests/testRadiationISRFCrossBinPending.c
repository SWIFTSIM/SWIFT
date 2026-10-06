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

/* Energy conservation of the ISRF force-loop pair terms across time bins, at
 * a = 0 and uniform c_hyp, through the real drift, extra ghost, type-2 force
 * dispatch and end_force. Two cells of particles on mixed time bins are
 * stepped over two coarse steps. Per-owner integration of the pair rates
 * leaks L = c m_i m_j d [sum_k G(t_k) - R G(t_R)] per cross-bin pair; with
 * every pair amount booked once from the finer member's timeline, sum_i m_i
 * u_i is conserved to round-off at every synchronised time. */
#if defined(FEEDBACK_GEAR) && defined(SPHENIX_SPH)

#include "feedback/GEAR/radiation_propagation_iact.h"

/* The pending fields exist only with the cross-bin booking. */
#ifdef FEEDBACK_RESTART_ISRF_PART_LAYOUT
#define HAVE_CROSS_BIN_PENDING 1
#endif

#define NODE_ID 0
#define CELL_N 6
#define COUNT (CELL_N * CELL_N * CELL_N)

/** binary32 unit round-off. */
#define U_F32 (0.5 * (double)FLT_EPSILON)

/* Round-off allowance of one float accumulation of n terms: n - 1 additions
 * (Higham Thm 4.1) plus at most 10 roundings in forming each pair term. */
#define SUM_ERR_UNITS(n) ((double)(n) + 10.)

void runner_dopair2_branch_force(struct runner *r, struct cell *ci,
                                 struct cell *cj, int limit_h_min,
                                 int limit_h_max);
void runner_doself2_branch_force(struct runner *r, const struct cell *c,
                                 int limit_h_min, int limit_h_max);

/* The fixture, file scope: struct engine is large. */
static struct engine engine;
static struct space space;
static struct cosmology cosmo;
static struct phys_const phys_const;
static struct feedback_props fp;
static struct cooling_function_data cooling;
static struct unit_system us;
static struct hydro_props hydro_props;
static struct pressure_floor_props pressure_floor;
static struct runner *runner;

/**
 * @brief Finiteness from the bit pattern, robust to -ffinite-math-only.
 *
 * @param x The value.
 * @return 1 if the exponent field is not all ones.
 */
static int is_finite_bits(double x) {
  uint64_t bits;
  memcpy(&bits, &x, sizeof(bits));
  volatile uint64_t opaque = bits;
  return ((opaque >> 52) & 0x7FF) != 0x7FF;
}

/**
 * @brief Set up a non-cosmological engine with scheme 2 (uniform c_hyp = 1).
 *
 * @param nr_nodes Number of ranks the engine claims; above 1 forces the
 * legacy pair booking.
 * @param scheme The ISRF c_hyp scheme.
 * @param time_base The time base: the physical step of bin b is 2^(b+1) of it.
 */
static void make_engine(int nr_nodes, int scheme, double time_base) {

  bzero(&space, sizeof(struct space));
  space.periodic = 0;
  for (int k = 0; k < 3; k++) space.dim[k] = 3.;

  cosmology_init_no_cosmo(&cosmo);

  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_speed_light_c = 1.e4;
  phys_const.const_vacuum_permeability = 1.;
  phys_const.const_proton_mass = 1.;

  bzero(&fp, sizeof(struct feedback_props));
  fp.ISRF_propagation = 1;
  fp.ISRF_c_hyp_scheme = scheme;
  fp.ISRF_c_hyp_fixed_fraction_of_c =
      scheme == isrf_c_hyp_scheme_fixed_fraction ? 1.e-4f : 0.f;
  fp.ISRF_c_hyp_margin = 0.5f;
  fp.ISRF_dissipation_alpha_max = 0.5f;
  fp.ISRF_dissipation_negativity_threshold = 0.01f;
  fp.ISRF_dissipation_alpha_floor = 0.5f;
  fp.ISRF_dissipation_floor_h_over_lambda = 0.5f;
  fp.ISRF_dissipation_floor_relaxation_residual = 0.f;
  fp.band_edge_weight_pe = 1.;
  fp.band_edge_weight_lw = 1.;
  fp.band_edge_photon_weight_lw = 1.;

  bzero(&cooling, sizeof(struct cooling_function_data));
  cooling.chemistry_data.local_dust_to_gas_ratio = 0.01;

  units_init(&us, 1.98892e33, 3.08567758e18, 3.15576e13, 1., 1.);

  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  hydro_props_init_no_hydro(&hydro_props);

  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.nodeID = NODE_ID;
  engine.nr_nodes = nr_nodes;
  engine.policy = 0;
  engine.time_base = time_base;
  engine.ti_current = 0;
  engine.max_active_bin = num_time_bins;
  engine.physical_constants = &phys_const;
  engine.cosmology = &cosmo;
  engine.internal_units = &us;
  engine.cooling_func = &cooling;
  engine.feedback_props = &fp;
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  runner->e = &engine;
}

/** How the particle state of a cell is drawn. */
enum fill_kind {
  fill_random, /*!< Random u, F, masses and densities on mixed bins. */
  fill_hot,    /*!< Hot coarse particles: u = 1, F = 0, mass m. */
  fill_cold    /*!< Cold fine particles: u = 0, F = 0, mass 2m. */
};

/**
 * @brief Build a cell of COUNT perturbed-lattice gas particles.
 *
 * @param offset The cell's lower corner.
 * @param h_spacing Smoothing length in units of the particle spacing.
 * @param depth_h Depth tag given to every particle.
 * @param bins Time bins to draw from.
 * @param n_bins Number of entries in @p bins.
 * @param kind How the ISRF state is drawn.
 * @param part_id (in/out) Running particle id.
 * @return The cell.
 */
static struct cell *make_cell(const double offset[3], double h_spacing,
                              char depth_h, const timebin_t *bins, int n_bins,
                              enum fill_kind kind, long long *part_id) {

  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate the cell");
  bzero(c, sizeof(struct cell));
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     COUNT * sizeof(struct part)) != 0)
    error("Couldn't allocate the particles");
  bzero(c->hydro.parts, COUNT * sizeof(struct part));
  if (posix_memalign((void **)&c->hydro.xparts, part_align,
                     COUNT * sizeof(struct xpart)) != 0)
    error("Couldn't allocate the extra particles");
  bzero(c->hydro.xparts, COUNT * sizeof(struct xpart));

  float h_max = 0.f;
  struct part *p = c->hydro.parts;
  for (int x = 0; x < CELL_N; x++) {
    for (int y = 0; y < CELL_N; y++) {
      for (int z = 0; z < CELL_N; z++) {
        const int idx[3] = {x, y, z};
        /* On a 2^-20 grid, so the dispatch's frame-shifted dx and the
         * reference's (float)(x_i - x_j) are the same float. */
        for (int k = 0; k < 3; k++) {
          const double x_exact =
              offset[k] + (idx[k] + 0.5 + random_uniform(-0.2, 0.2)) / CELL_N;
          p->x[k] = round(x_exact * 1048576.) / 1048576.;
        }
        p->h = h_spacing * random_uniform(1., 1.2) / CELL_N;
        h_max = max(h_max, p->h);
        p->id = ++(*part_id);
        p->depth_h = depth_h;
        p->u = 1.f;
        p->rho = random_uniform(0.5, 2.);
        p->time_bin = bins[rand() % n_bins];

        struct feedback_part_data *fd = &p->feedback_data;
        fd->ISRF_reservoir_end_ti = -1;
        fd->ISRF_illumination_end_ti = -1;
        double u[ISRF_OPERATOR_COUNT];
        if (kind == fill_random) {
          p->mass = random_uniform(0.5, 2.) / COUNT;
          for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
            u[o] = random_uniform(0.2, 1.);
        } else {
          p->mass = (kind == fill_cold ? 2. : 1.) / COUNT;
          for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
            u[o] = kind == fill_hot ? 1. : 0.;
        }
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
          struct feedback_isrf_moment_data *mo = &fd->isrf_moment[m];
          const int o = radiation_isrf_moment_to_operator[m];
          mo->u = u[o];
          for (int k = 0; k < 3; k++)
            mo->specific_flux[k] =
                kind == fill_random ? 0.5 * u[o] * random_uniform(-1., 1.) : 0.;
        }
        for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
          fd->isrf_operator[o].dissipation_alpha_floor = 0.5f;
        p++;
      }
    }
  }

  c->split = 0;
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max;
  c->hydro.count = COUNT;
  for (int k = 0; k < 3; k++) {
    c->width[k] = 1.;
    c->loc[k] = offset[k];
  }
  c->dmin = 1.;
  c->hydro.super = c;
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
 * @brief Set the SPH force-loop inputs. Touches no ISRF field.
 *
 * @param c The cell.
 */
static void prepare_sph(struct cell *c) {
  for (int i = 0; i < c->hydro.count; i++) {
    struct part *p = &c->hydro.parts[i];
    p->density.rho_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
    p->force.pressure = hydro_get_comoving_pressure(p);
    p->viscosity.alpha = 0.8;
    p->viscosity.div_v = 0.f;
    p->viscosity.div_v_previous_step = 0.f;
    p->viscosity.v_sig = hydro_get_comoving_soundspeed(p);
    hydro_prepare_force(p, &c->hydro.xparts[i], &cosmo, &hydro_props,
                        &pressure_floor, 0., 0.);
    hydro_reset_acceleration(p);
  }
}

/**
 * @brief Set a cell's depth and the smoothing-length window of its sweeps.
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
 * @brief Absolute round-off of one #kernel_deval gradient, in units of
 * #U_F32 times the kernel's dW/dx normalisation (Horner conditioning).
 *
 * @param x The kernel argument r / H, in double.
 * @return The bound, 0 outside the kernel.
 */
static double kernel_deval_error_units(double x) {
  if (!(x < 1.)) return 0.;
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const c = &kernel_coeffs[ind * (kernel_degree + 1)];
  double abs_dp = 0., d2p = 0.;
  for (int j = 0; j < kernel_degree; j++) {
    const int pw = kernel_degree - j;
    abs_dp += pw * fabs((double)c[j]) * pow(x, pw - 1);
    if (pw >= 2) d2p += pw * (pw - 1) * (double)c[j] * pow(x, pw - 2);
  }
  return 2. * kernel_degree * abs_dp + 5.5 * x * fabs(d2p);
}

/**
 * @brief One side's kernel-gradient term `dW/dx * h^-(d+1)` of a pair and its
 * absolute round-off bound.
 *
 * @param r The pair distance, in double.
 * @param h The side's smoothing length.
 * @param w_dr (return) |w_dr| in double.
 * @return The absolute bound on the float32 w_dr.
 */
static double side_w_dr_error(double r, float h, double *w_dr) {
  const double x = r / (double)h * (double)kernel_gamma_inv;
  const double scale = (double)kernel_constant *
                       (double)kernel_gamma_inv_dim_plus_one *
                       pow((double)h, -(hydro_dimension + 1));
  double dw = 0.;
  if (x < 1.) {
    const int temp = (int)(x * kernel_ivals_f);
    const int ind = temp > kernel_ivals ? kernel_ivals : temp;
    const float *const c = &kernel_coeffs[ind * (kernel_degree + 1)];
    for (int j = 0; j < kernel_degree; j++)
      dw += (kernel_degree - j) * (double)c[j] * pow(x, kernel_degree - j - 1);
  }
  *w_dr = fabs(dw) * scale;
  return U_F32 * (scale * kernel_deval_error_units(x) + 10. * *w_dr);
}

/**
 * @brief Forward-error bound of the float32 pair terms written into pi. The
 * library and this file may evaluate the kernel gradient differently.
 *
 * Same bound as the dispatch conservation test (Horner conditioning of the
 * kernel gradient, then the products of the two pair terms).
 *
 * @param pi The receiving particle.
 * @param pj The neighbour.
 * @param r2 The pair's float32 squared distance.
 * @param dx The pair's float32 separation.
 * @param b The moment.
 * @return Bound on the sum of the divergence and dissipation term errors.
 */
static double pair_term_error(const struct part *pi, const struct part *pj,
                              float r2, const float dx[3], int b) {
  const struct feedback_part_data *fi = &pi->feedback_data;
  const struct feedback_part_data *fj = &pj->feedback_data;
  const double r = sqrt((double)r2);
  double wi_dr, wj_dr;
  const double ei = side_w_dr_error(r, pi->h, &wi_dr);
  const double ej = side_w_dr_error(r, pj->h, &wj_dr);
  const double rho_i = fi->rho_prev, rho_j = fj->rho_prev;
  const double ci_mj = (double)fi->c_hyp * (double)pj->mass;

  double dot_i = 0., dot_j = 0.;
  for (int k = 0; k < 3; k++) {
    dot_i += fabs((double)fi->isrf_moment[b].specific_flux[k] * dx[k]);
    dot_j += fabs((double)fj->isrf_moment[b].specific_flux[k] * dx[k]);
  }
  const double err_div = ci_mj / r *
                         (dot_i / rho_i * (13.5 * U_F32 * wi_dr + ei) +
                          dot_j / rho_j * (13.5 * U_F32 * wj_dr + ej));

  const int o = radiation_isrf_moment_to_operator[b];
  const double alpha =
      fmax(fmax(fi->isrf_operator[o].dissipation_alpha_trigger,
                fj->isrf_operator[o].dissipation_alpha_trigger),
           fmax(fi->isrf_operator[o].dissipation_alpha_floor,
                fj->isrf_operator[o].dissipation_alpha_floor));
  const double rui = fabs(rho_i * (double)fi->isrf_moment[b].u);
  const double ruj = fabs(rho_j * (double)fj->isrf_moment[b].u);
  const double d = fabs(rho_i * (double)fi->isrf_moment[b].u -
                        rho_j * (double)fj->isrf_moment[b].u);
  const double wbar = 0.5 * (wi_dr + wj_dr);
  const double err_diss = ci_mj * alpha / (rho_i * rho_j) *
                          (2. * U_F32 * (rui + ruj) * wbar +
                           8. * U_F32 * d * wbar + 0.5 * d * (ei + ej));
  return err_div + err_diss;
}

/** Per-particle reference sums of one step, in double. */
struct reference {
  double S[ISRF_MOMENT_COUNT];       /*!< Sum of |rate terms|. */
  int n;                             /*!< Number of rate terms. */
  double cross[ISRF_MOMENT_COUNT];   /*!< Cross-bin rate terms, diss - div. */
  double E_cross[ISRF_MOMENT_COUNT]; /*!< Error bound of #cross. */
  double Pt[ISRF_MOMENT_COUNT];      /*!< Expected transport pending. */
  double Pd[ISRF_MOMENT_COUNT];      /*!< Expected dissipation pending. */
  double SP[ISRF_MOMENT_COUNT];      /*!< Sum of |deposits|. */
  double E_P[ISRF_MOMENT_COUNT];     /*!< Error bound of the deposits. */
  int n_dep;                         /*!< Number of deposits. */
};

/** A two-cell system and its reference state. */
struct system {
  struct cell *cells[2];
  struct reference *ref[2];
  float h_split;
  int two_levels;
};

/**
 * @brief Brute-force pair terms of every active particle, from the legacy
 * symmetric iact on copies: the rate terms and the deposits owed to coarser
 * partners.
 *
 * @param sys The system.
 */
static void reference_pairs(struct system *sys) {

  for (int ci = 0; ci < 2; ci++) {
    for (int i = 0; i < COUNT; i++) {
      const struct part *pi = &sys->cells[ci]->hydro.parts[i];
      if (!part_is_active(pi, &engine)) continue;
      struct reference *ri = &sys->ref[ci][i];
      const double dt_i = pi->feedback_data.dt_prev;

      for (int cj = 0; cj < 2; cj++) {
        for (int j = 0; j < COUNT; j++) {
          const struct part *pj = &sys->cells[cj]->hydro.parts[j];
          if (pj == pi) continue;
          const float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                               (float)(pi->x[1] - pj->x[1]),
                               (float)(pi->x[2] - pj->x[2])};
          const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
          const float hig2 = pi->h * pi->h * kernel_gamma2;
          const float hjg2 = pj->h * pj->h * kernel_gamma2;
          if (!(r2 < hig2) && !(r2 < hjg2)) continue;

          struct part ti = *pi, tj = *pj;
          for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
            ti.feedback_data.isrf_moment[m].div_specific_flux = 0.f;
            ti.feedback_data.isrf_moment[m].dissipation_u = 0.f;
            tj.feedback_data.isrf_moment[m].div_specific_flux = 0.f;
            tj.feedback_data.isrf_moment[m].dissipation_u = 0.f;
          }
#ifdef HAVE_CROSS_BIN_PENDING
          ti.feedback_data.dt_active = 0.f;
          tj.feedback_data.dt_active = 0.f;
#endif
          runner_iact_isrf_dissipation(r2, dx, pi->h, pj->h, &ti, &tj, 1.f,
                                       0.f);

          struct reference *rj = &sys->ref[cj][j];
          const int coarser_j = pj->time_bin > pi->time_bin;
          ri->n++;
          if (coarser_j) rj->n_dep++;
          for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
            const struct feedback_isrf_moment_data *mi =
                &ti.feedback_data.isrf_moment[m];
            const struct feedback_isrf_moment_data *mj =
                &tj.feedback_data.isrf_moment[m];
            ri->S[m] += fabs(mi->div_specific_flux) + fabs(mi->dissipation_u);
            if (pj->time_bin != pi->time_bin) {
              ri->cross[m] +=
                  (double)mi->dissipation_u - (double)mi->div_specific_flux;
              ri->E_cross[m] += pair_term_error(pi, pj, r2, dx, m);
            }
            if (coarser_j) {
              const double c_j = pj->feedback_data.c_hyp;
              rj->Pt[m] += -(double)mj->div_specific_flux / c_j * dt_i;
              rj->Pd[m] += (double)mj->dissipation_u / c_j * dt_i;
              rj->SP[m] +=
                  (fabs(mj->div_specific_flux) + fabs(mj->dissipation_u)) /
                  c_j * dt_i;
              rj->E_P[m] += pair_term_error(pj, pi, r2, dx, m) / c_j * dt_i;
            }
          }
        }
      }
    }
  }
}

/**
 * @brief One pass of the real force dispatch over the two cells, with or
 * without the two-level depth split of the sub-cell recursion.
 *
 * @param sys The system.
 */
static void force_dispatch(struct system *sys) {
  struct cell **cells = sys->cells;
  if (!sys->two_levels) {
    for (int k = 0; k < 2; k++)
      runner_doself2_branch_force(runner, cells[k], 0, 0);
    runner_dopair2_branch_force(runner, cells[0], cells[1], 0, 0);
    return;
  }
  for (int k = 0; k < 2; k++) set_level(cells[k], 1, 0.f, sys->h_split);
  for (int k = 0; k < 2; k++)
    runner_doself2_branch_force(runner, cells[k], 0, 1);
  runner_dopair2_branch_force(runner, cells[0], cells[1], 0, 1);
  for (int k = 0; k < 2; k++)
    set_level(cells[k], 0, sys->h_split, 1.f / kernel_gamma);
  for (int k = 0; k < 2; k++)
    runner_doself2_branch_force(runner, cells[k], 1, 1);
  runner_dopair2_branch_force(runner, cells[0], cells[1], 1, 1);
}

/**
 * @brief Field energy sum_i m_i u_i of one moment, plus the pending amounts.
 *
 * @param sys The system.
 * @param m The moment.
 * @param with_pending Add c_hyp * m * pending.
 * @param abs_sum (return) sum_i m_i |u_i|.
 * @return The sum.
 */
static double ledger(const struct system *sys, int m, int with_pending,
                     double *abs_sum) {
  double sum = 0.;
  *abs_sum = 0.;
  for (int c = 0; c < 2; c++) {
    for (int i = 0; i < COUNT; i++) {
      const struct part *p = &sys->cells[c]->hydro.parts[i];
      const struct feedback_isrf_moment_data *mo =
          &p->feedback_data.isrf_moment[m];
      double x = (double)p->mass * mo->u;
#ifdef HAVE_CROSS_BIN_PENDING
      if (with_pending)
        x += (double)p->mass * (double)p->feedback_data.c_hyp *
             ((double)mo->pending_transport_u +
              (double)mo->pending_dissipation_u);
#endif
      if (!is_finite_bits(x)) error("particle %lld: non-finite ledger", p->id);
      sum += x;
      *abs_sum += fabs(x);
    }
  }
  return sum;
}

/** Outcome of one stepped case. */
struct outcome {
  double residual[ISRF_MOMENT_COUNT]; /*!< Final minus initial ledger. */
  double L_pred[ISRF_MOMENT_COUNT];   /*!< Closed-form cross-bin leak. */
  double bar[ISRF_MOMENT_COUNT];      /*!< Round-off bar of the residual. */
  double L_bar[ISRF_MOMENT_COUNT];    /*!< Extra bar of #L_pred, two TUs. */
  double E0[ISRF_MOMENT_COUNT];       /*!< Initial ledger. */
  double min_u_coarse;                /*!< Smallest PE u of the coarsest bin. */
  int routed;                         /*!< Cross-bin booking active. */
};

/**
 * @brief Step a two-cell system through the real drift, extra ghost, force
 * dispatch and end_force, and book the ledger.
 *
 * @param sys The system (cells built, sorted, SPH inputs set).
 * @param n_steps Number of fine steps (bin of the finest particle).
 * @param out (return) The outcome.
 */
static void step_system(struct system *sys, int n_steps, struct outcome *out) {

  timebin_t bin_min = num_time_bins, bin_max = 0;
  for (int c = 0; c < 2; c++)
    for (int i = 0; i < COUNT; i++) {
      bin_min = min(bin_min, sys->cells[c]->hydro.parts[i].time_bin);
      bin_max = max(bin_max, sys->cells[c]->hydro.parts[i].time_bin);
    }
  const integertime_t dti_fine = get_integer_timestep(bin_min);
  const integertime_t dti_coarse = get_integer_timestep(bin_max);

  double abs_sum;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    out->E0[m] = ledger(sys, m, 1, &abs_sum);
    out->L_pred[m] = 0.;
    out->L_bar[m] = 0.;
    out->bar[m] = 4. * DBL_EPSILON * abs_sum;
  }
  out->routed = 0;

  for (int c = 0; c < 2; c++) {
    sys->ref[c] = calloc(COUNT, sizeof(struct reference));
    if (sys->ref[c] == NULL) error("Couldn't allocate the reference");
  }

  for (int step = 1; step <= n_steps; step++) {
    const integertime_t T = step * dti_fine;
    engine.ti_current = T;
    engine.max_active_bin = get_max_active_bin(T);

    /* Drift: every particle, as every cell in a force pair is drifted. */
    for (int c = 0; c < 2; c++) {
      struct cell *cell = sys->cells[c];
      cell->hydro.ti_old_part = T;
      integertime_t ti_end_min = max_nr_timesteps;
      for (int i = 0; i < COUNT; i++) {
        struct part *p = &cell->hydro.parts[i];
#ifdef SWIFT_DEBUG_CHECKS
        p->ti_drift = T;
        p->ti_kick = T;
#endif
        radiation_snapshot_part_propagation(p, &engine);
        ti_end_min = min(ti_end_min, get_integer_time_end(T, p->time_bin));
      }
      cell->hydro.ti_end_min = ti_end_min;
    }

    /* Extra ghost on the active particles. */
    for (int c = 0; c < 2; c++)
      for (int i = 0; i < COUNT; i++) {
        struct part *p = &sys->cells[c]->hydro.parts[i];
        if (part_is_active(p, &engine))
          radiation_end_gradient_propagation(p, &engine);
#ifdef HAVE_CROSS_BIN_PENDING
        if (p->feedback_data.dt_active != 0.f) out->routed = 1;
#endif
      }

    reference_pairs(sys);
    force_dispatch(sys);

#ifdef HAVE_CROSS_BIN_PENDING
    /* Every particle's pending against its brute-force deposits. */
    for (int c = 0; c < 2; c++)
      for (int i = 0; i < COUNT; i++) {
        const struct part *p = &sys->cells[c]->hydro.parts[i];
        const struct reference *ri = &sys->ref[c][i];
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
          const struct feedback_isrf_moment_data *mo =
              &p->feedback_data.isrf_moment[m];
          const double eP = out->routed ? ri->Pt[m] : 0.;
          const double eD = out->routed ? ri->Pd[m] : 0.;
          const double bar =
              SUM_ERR_UNITS(ri->n_dep) * U_F32 * ri->SP[m] + 2. * ri->E_P[m];
          const double et = fabs(mo->pending_transport_u - eP);
          const double ed = fabs(mo->pending_dissipation_u - eD);
          if (!is_finite_bits(et) || !is_finite_bits(ed) || et > bar ||
              ed > bar)
            error(
                "particle %lld bin %d moment %d: pending (%e, %e) against "
                "brute force (%e, %e), bar %e",
                p->id, p->time_bin, m, mo->pending_transport_u,
                mo->pending_dissipation_u, eP, eD, bar);
        }
      }
#endif

    /* End force on the active particles; close their reference sums. */
    for (int c = 0; c < 2; c++)
      for (int i = 0; i < COUNT; i++) {
        struct part *p = &sys->cells[c]->hydro.parts[i];
        struct reference *ri = &sys->ref[c][i];
        if (!part_is_active(p, &engine)) continue;
        radiation_end_force_propagation(p, &engine);
        const double m_dt = (double)p->mass * p->feedback_data.dt_prev;
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
          out->L_pred[m] += m_dt * ri->cross[m];
          out->L_bar[m] += m_dt * 2. * ri->E_cross[m];
          out->bar[m] += m_dt * SUM_ERR_UNITS(ri->n) * U_F32 * ri->S[m] +
                         (double)p->mass * p->feedback_data.c_hyp *
                             SUM_ERR_UNITS(ri->n_dep) * U_F32 * ri->SP[m];
          ri->S[m] = ri->cross[m] = ri->E_cross[m] = 0.;
          ri->Pt[m] = ri->Pd[m] = ri->SP[m] = ri->E_P[m] = 0.;
#ifdef HAVE_CROSS_BIN_PENDING
          const struct feedback_isrf_moment_data *mo =
              &p->feedback_data.isrf_moment[m];
          if (mo->pending_transport_u != 0.f ||
              mo->pending_dissipation_u != 0.f)
            error("particle %lld: pending not consumed", p->id);
#endif
        }
        ri->n = ri->n_dep = 0;
      }

#ifdef HAVE_CROSS_BIN_PENDING
    /* With the pending amounts, the ledger holds after every step. */
    if (out->routed)
      for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
        double in_flight = 0.;
        for (int c = 0; c < 2; c++)
          for (int i = 0; i < COUNT; i++) {
            const struct part *p = &sys->cells[c]->hydro.parts[i];
            const struct reference *ri = &sys->ref[c][i];
            in_flight += (double)p->mass * p->feedback_data.c_hyp *
                         SUM_ERR_UNITS(ri->n_dep) * U_F32 * ri->SP[m];
          }
        const double E = ledger(sys, m, 1, &abs_sum);
        if (!(fabs(E - out->E0[m]) <= out->bar[m] + in_flight))
          error("step %d moment %d: ledger with pending off by %e, bar %e",
                step, m, E - out->E0[m], out->bar[m] + in_flight);
      }
#endif
  }

  if ((n_steps * dti_fine) % dti_coarse != 0)
    error("the case must end on a synchronised time");

  out->min_u_coarse = DBL_MAX;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    out->residual[m] = ledger(sys, m, 0, &abs_sum) - out->E0[m];
  for (int c = 0; c < 2; c++) {
    for (int i = 0; i < COUNT; i++) {
      const struct part *p = &sys->cells[c]->hydro.parts[i];
      if (p->time_bin == bin_max)
        out->min_u_coarse = min(out->min_u_coarse,
                                p->feedback_data.isrf_moment[ISRF_MOMENT_PE].u);
    }
    free(sys->ref[c]);
  }
}

/**
 * @brief Build two cells (small and large h), set them up for the dispatch.
 *
 * @param sys (return) The system.
 * @param offset Offset of the second cell.
 * @param small_first 1 to put the small-h cell at the origin.
 * @param two_levels 1 for the depth-split sweep.
 * @param bins Time bins to draw from.
 * @param n_bins Number of entries in @p bins.
 * @param kinds How each cell's state is drawn (small, large).
 * @param part_id (in/out) Running particle id.
 */
static void make_system(struct system *sys, const double offset[3],
                        int small_first, int two_levels, const timebin_t *bins,
                        int n_bins, const enum fill_kind kinds[2],
                        long long *part_id) {
  const double origin[3] = {0., 0., 0.};
  const double h_small = 1.2348;
  sys->h_split = 1.5 * h_small / CELL_N;
  sys->two_levels = two_levels;
  sys->cells[0] = make_cell(small_first ? origin : offset, h_small, 1, bins,
                            n_bins, kinds[0], part_id);
  sys->cells[1] = make_cell(small_first ? offset : origin, 2. * h_small, 0,
                            bins, n_bins, kinds[1], part_id);
  for (int k = 0; k < 2; k++) {
    struct cell *c = sys->cells[k];
    c->hydro.ti_old_part = 0;
    c->hydro.ti_end_min = 0;
    runner_do_hydro_sort(runner, c, 0x1FFF, 0, 0, 0, 0);
    prepare_sph(c);
    if (!two_levels) {
      c->depth = 0;
      for (int i = 0; i < COUNT; i++) c->hydro.parts[i].depth_h = 0;
    }
  }
}

/**
 * @brief Run one mixed-bin case and print its ledger.
 *
 * @param label Case name.
 * @param nr_nodes Engine rank count (above 1: legacy booking).
 * @param bins Time bins to draw from.
 * @param n_bins Number of entries in @p bins.
 * @param two_levels 1 for the depth-split sweep.
 * @param small_first 1 to put the small-h cell at the origin.
 * @param out (return) The outcome.
 */
static void run_case(const char *label, int nr_nodes, const timebin_t *bins,
                     int n_bins, int two_levels, int small_first,
                     struct outcome *out) {
  static long long part_id = 0;
  make_engine(nr_nodes, isrf_c_hyp_scheme_fixed_fraction, 1.e-2);
  srand(4321 + 17 * n_bins + 3 * two_levels + small_first);
  const double offset[3] = {1., 0., 0.};
  const enum fill_kind kinds[2] = {fill_random, fill_random};
  struct system sys;
  make_system(&sys, offset, small_first, two_levels, bins, n_bins, kinds,
              &part_id);

  timebin_t bin_max = 0;
  for (int b = 0; b < n_bins; b++) bin_max = max(bin_max, bins[b]);
  const int R = 1 << (bin_max - bins[0]);
  step_system(&sys, 2 * R, out);
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    if (!is_finite_bits(out->residual[m]) || !is_finite_bits(out->L_pred[m]) ||
        !is_finite_bits(out->bar[m]) || !(out->bar[m] > 0.))
      error("%s moment %d: non-finite or empty sums", label, m);
    message(
        "%s moment %d: residual %+.6e  closed form %+.6e  bar %.3e  (of E0: "
        "%.3e, %.3e, %.3e) routed %d",
        label, m, out->residual[m], out->L_pred[m], out->bar[m],
        out->residual[m] / out->E0[m], out->L_pred[m] / out->E0[m],
        out->bar[m] / out->E0[m], out->routed);
  }
  clean_up(sys.cells[0]);
  clean_up(sys.cells[1]);
}

/* Bin sets: R = 2, R = 4, and three bins. The first entry is the finest. */
static const timebin_t bins_r2[2] = {1, 2};
static const timebin_t bins_r4[2] = {1, 3};
static const timebin_t bins_r124[3] = {1, 2, 3};

/**
 * @brief Legacy booking (several ranks): its residual is the closed-form
 * cross-bin leak, which is well above the round-off bar.
 */
static void test_legacy_residual_is_closed_form(void) {
  struct outcome out;
  run_case("legacy R=4", 2, bins_r4, 2, 0, 1, &out);
  if (out.routed) error("legacy R=4: a particle left state 0 on two ranks");
  double resolved = 0.;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    if (!(fabs(out.residual[m] - out.L_pred[m]) <= out.bar[m] + out.L_bar[m]))
      error("legacy moment %d: residual %e is not the closed form %e (bar %e)",
            m, out.residual[m], out.L_pred[m], out.bar[m] + out.L_bar[m]);
    resolved = max(resolved, fabs(out.L_pred[m]) / out.bar[m]);
  }
  if (!(resolved > 10.)) error("legacy: no moment has a resolvable leak");
}

/**
 * @brief One rank, uniform c_hyp, a = 0: the cross-bin leak must vanish to
 * round-off at the synchronised end, in every geometry and depth split.
 */
static void test_routed_conserves(void) {
  const timebin_t *sets[3] = {bins_r2, bins_r4, bins_r124};
  const int n_bins[3] = {2, 2, 3};
  for (int s = 0; s < 3; s++)
    for (int two_levels = 0; two_levels < 2; two_levels++)
      for (int small_first = 0; small_first < 2; small_first++) {
        char label[64];
        sprintf(label, "routed bins %d-%d levels %d first %d", sets[s][0],
                sets[s][n_bins[s] - 1], two_levels + 1, small_first);
        struct outcome out;
        run_case(label, 1, sets[s], n_bins[s], two_levels, small_first, &out);
        double resolved = 0.;
        for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
          if (!(fabs(out.residual[m]) <= out.bar[m]))
            error("%s moment %d: residual %e above its round-off bar %e", label,
                  m, out.residual[m], out.bar[m]);
          resolved = max(resolved, fabs(out.L_pred[m]) / out.bar[m]);
        }
        if (!(resolved > 10.))
          error("%s: no moment has a resolvable leak to remove", label);
      }
}

/**
 * @brief Report only: a hot coarse cell next to a cold, twice heavier fine
 * cell, dissipation only. Smallest u of the coarse particles, both bookings.
 */
static void test_report_adverse_positivity(void) {
  for (int nr_nodes = 2; nr_nodes >= 1; nr_nodes--) {
    make_engine(nr_nodes, isrf_c_hyp_scheme_fixed_fraction, 1.e-2);
    srand(99);
    long long part_id = 1000000;
    const double offset[3] = {1., 0., 0.};
    const timebin_t coarse[1] = {4};
    const timebin_t fine[1] = {1};
    struct system sys;
    const double origin[3] = {0., 0., 0.};
    sys.h_split = 0.f;
    sys.two_levels = 0;
    sys.cells[0] = make_cell(origin, 2.4696, 0, coarse, 1, fill_hot, &part_id);
    sys.cells[1] = make_cell(offset, 2.4696, 0, fine, 1, fill_cold, &part_id);
    for (int k = 0; k < 2; k++) {
      sys.cells[k]->hydro.ti_old_part = 0;
      runner_do_hydro_sort(runner, sys.cells[k], 0x1FFF, 0, 0, 0, 0);
      prepare_sph(sys.cells[k]);
    }
    struct outcome out;
    step_system(&sys, 8, &out);
    message(
        "adverse geometry (report only), %s booking: min u of the hot coarse "
        "cell %+.4e, PE residual %+.3e of E0",
        nr_nodes > 1 ? "legacy" : "routed", out.min_u_coarse,
        out.residual[ISRF_MOMENT_PE] / out.E0[ISRF_MOMENT_PE]);
    clean_up(sys.cells[0]);
    clean_up(sys.cells[1]);
  }
}

#ifdef HAVE_CROSS_BIN_PENDING

/**
 * @brief Fork, run @p fn, and require that it stops with error() and prints
 * @p expected_text.
 *
 * @param label Message label.
 * @param fn The function to run in the child.
 * @param expected_text Substring of the expected error message.
 */
static void assert_aborts_with_error(const char *label, void (*fn)(void),
                                     const char *expected_text) {
  char path[] = "crossbinXXXXXX";
  const int fd = mkstemp(path);
  if (fd < 0) error("mkstemp() failed in %s.", label);
  close(fd);
  const pid_t pid = fork();
  if (pid == 0) {
    if (freopen(path, "w", stderr) == NULL) _exit(43);
    fn();
    _exit(42);
  } else if (pid > 0) {
    int status;
    waitpid(pid, &status, 0);
    const int failed = (WIFEXITED(status) && WEXITSTATUS(status) == 1) ||
                       (WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
    char text[4096];
    size_t n = 0;
    FILE *f = fopen(path, "r");
    if (f != NULL) {
      n = fread(text, 1, sizeof(text) - 1, f);
      fclose(f);
    }
    text[n] = '\0';
    unlink(path);
    if (!failed) error("%s: expected error(), got status %d", label, status);
    if (strstr(text, expected_text) == NULL)
      error("%s: wrong message: %s", label, text);
  } else {
    error("fork() failed in %s.", label);
  }
}

/**
 * @brief A cross-bin deposit from a dusty fine member (a > 0), which the a = 0
 * booking cannot weight.
 */
static void run_dusty_fine_deposit(void) {
  make_engine(1, isrf_c_hyp_scheme_fixed_fraction, 2.5e-3);
  srand(7);
  long long part_id = 2000000;
  const double offset[3] = {1., 0., 0.};
  const enum fill_kind kinds[2] = {fill_random, fill_random};
  struct system sys;
  make_system(&sys, offset, 1, 0, bins_r4, 2, kinds, &part_id);
  engine.ti_current = 4;
  engine.max_active_bin = get_max_active_bin(4);
  for (int c = 0; c < 2; c++) {
    sys.cells[c]->hydro.ti_old_part = 4;
    sys.cells[c]->hydro.ti_end_min = 4;
    for (int i = 0; i < COUNT; i++) {
      struct part *p = &sys.cells[c]->hydro.parts[i];
#ifdef SWIFT_DEBUG_CHECKS
      p->ti_drift = 4;
      p->ti_kick = 4;
#endif
      radiation_snapshot_part_propagation(p, &engine);
      for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
        p->feedback_data.isrf_operator[o].kappa = 1.f;
      if (part_is_active(p, &engine))
        radiation_end_gradient_propagation(p, &engine);
    }
  }
  force_dispatch(&sys);
}

/**
 * @brief The sentinel: -1 inactive and dt_prev active on one rank under
 * scheme 2; 0 on several ranks and under scheme 4. A deposit at a > 0 stops.
 */
static void test_sentinel_states(void) {
  const int schemes[2] = {isrf_c_hyp_scheme_fixed_fraction,
                          isrf_c_hyp_scheme_kernel_local_reduced_flux};
  for (int s = 0; s < 2; s++)
    for (int nr_nodes = 1; nr_nodes <= 2; nr_nodes++) {
      make_engine(nr_nodes, schemes[s], 2.5e-3);
      struct part p[2];
      for (int k = 0; k < 2; k++) {
        bzero(&p[k], sizeof(struct part));
        p[k].h = 0.3f;
        p[k].rho = 1.f;
        p[k].mass = 1.f;
        p[k].time_bin = k == 0 ? 1 : 3;
        p[k].feedback_data.ISRF_reservoir_end_ti = -1;
        p[k].feedback_data.c_hyp = 1.f;
      }
      engine.ti_current = 4;
      engine.max_active_bin = get_max_active_bin(4);
      for (int k = 0; k < 2; k++) {
        radiation_snapshot_part_propagation(&p[k], &engine);
        if (part_is_active(&p[k], &engine))
          radiation_end_gradient_propagation(&p[k], &engine);
      }
      const int routed = s == 0 && nr_nodes == 1;
      const float want_active = routed ? p[0].feedback_data.dt_prev : 0.f;
      const float want_inactive = routed ? -1.f : 0.f;
      if (p[0].feedback_data.dt_active != want_active ||
          p[1].feedback_data.dt_active != want_inactive ||
          (routed && !(want_active > 0.f)))
        error("scheme %d, %d ranks: sentinel (%e, %e), expected (%e, %e)",
              schemes[s], nr_nodes, p[0].feedback_data.dt_active,
              p[1].feedback_data.dt_active, want_active, want_inactive);
    }
  assert_aborts_with_error("dusty fine member deposit", run_dusty_fine_deposit,
                           "needs a = 0");
  message("sentinel states and the a > 0 stop: as expected");
}

#endif /* HAVE_CROSS_BIN_PENDING */

int main(int argc, char *argv[]) {

  clocks_set_cpufreq(0);

#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  if (posix_memalign((void **)&runner, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct runner)) != 0)
    error("Couldn't allocate the runner");
  bzero(runner, sizeof(struct runner));
#ifdef WITH_VECTORIZATION
  cache_init(&runner->ci_cache, 512);
  cache_init(&runner->cj_cache, 512);
#endif

#ifdef HAVE_CROSS_BIN_PENDING
  test_sentinel_states();
#endif
  test_legacy_residual_is_closed_form();
  test_report_adverse_positivity();
  test_routed_conserves();

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
