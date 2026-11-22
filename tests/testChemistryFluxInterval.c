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

/* Local headers. */
#include "swift.h"
#include "timestep_limiter.h"
#include "timestep_sync.h"

/* System headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file testChemistryFluxInterval.c
 * @brief Check that the FVPM chemistry integrates the flux of every pair over
 * the whole time line when time-steps are limited or synchronised.
 *
 * Three mutual neighbours start on steps of 32, 16 and 32 ticks (SWIFT issue
 * #832). At ti=20, particle 0 is woken by the limiter, synchronised, or both
 * particle 0 is woken and particle 1 synchronised. Each step runs the real
 * ghost-side preparation, the chemistry flux prep loop, the flux ghost, the
 * force loop and runner_do_end_hydro_force() on two cells, and kicks the
 * active particles.
 *
 * The state of the particles is restored before each step, so the flux rate
 * of each pair is constant. The metal mass moved in two elements with
 * different gradients gives 6 equations for the covered times of the 3
 * pairs, solved by least squares. Over [0, 32], each must be 32 ticks.
 */

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#define N_PARTS 3
#define TI_END 32
#define TI_EVENT 20

void runner_dosub_self2_force(struct runner *r, struct cell *c,
                              int recurse_below_h_max, const int gettimer);
void runner_dosub_pair2_force(struct runner *r, struct cell *ci,
                              struct cell *cj, int recurse_below_h_max,
                              const int gettimer);
void runner_dosub_self2_chemistry_flux_prep(struct runner *r, struct cell *c,
                                            int recurse_below_h_max,
                                            const int gettimer);
void runner_dosub_pair2_chemistry_flux_prep(struct runner *r, struct cell *ci,
                                            struct cell *cj,
                                            int recurse_below_h_max,
                                            const int gettimer);

enum schedule {
  schedule_no_event,
  schedule_limiter,
  schedule_sync,
  schedule_limiter_and_sync,
  schedule_count
};

static const char *schedule_names[schedule_count] = {
    "no event", "limiter on 0", "sync of 0", "limiter on 0, sync of 1"};

/* Time-bins of 32, 16 and 4 ticks. */
#define BIN_32 4
#define BIN_16 3
#define BIN_4 1

/* The pairs 0-1, 0-2 and 1-2, as (left, right) particles. */
static const int pair_i[3] = {0, 0, 1};
static const int pair_j[3] = {1, 2, 2};

/**
 * @brief The time-bin a particle takes at the end of a step at ti.
 */
static timebin_t next_bin(const int k, const integertime_t ti,
                          const timebin_t current) {
  if (ti >= TI_EVENT && current == BIN_4) return BIN_4;
  return (k == 1) ? BIN_16 : BIN_32;
}

/**
 * @brief Fill the cell properties derived from its particles.
 */
static void set_cell_properties(struct cell *c, const struct engine *e) {
  c->hydro.h_max = 0.f;
  c->hydro.h_max_active = 0.f;
  c->hydro.ti_end_min = max_nr_timesteps;
  for (int k = 0; k < c->hydro.count; k++) {
    const struct part *p = &c->hydro.parts[k];
    c->hydro.h_max = max(c->hydro.h_max, p->h);
    if (part_is_active(p, e))
      c->hydro.h_max_active = max(c->hydro.h_max_active, p->h);
    c->hydro.ti_end_min = min(c->hydro.ti_end_min,
                              get_integer_time_end(e->ti_current, p->time_bin));
  }
  c->hydro.ti_old_part = e->ti_current;
}

/**
 * @brief Kick a particle from ti_beg to ti_end.
 */
static void kick(struct part *p, struct xpart *xp, const struct engine *e,
                 const integertime_t ti_beg, const integertime_t ti_end) {
  const double tb = e->time_base;
  const struct cosmology *cosmo = e->cosmology;
  kick_part(p, xp, kick_get_hydro_kick_dt(ti_beg, ti_end, tb, 0, cosmo),
            kick_get_grav_kick_dt(ti_beg, ti_end, tb, 0, cosmo),
            /*dt_kick_mesh_grav=*/0.,
            kick_get_therm_kick_dt(ti_beg, ti_end, tb, 0, cosmo),
            kick_get_corr_kick_dt(ti_beg, ti_end, tb, 0, cosmo), cosmo,
            e->hydro_properties, e->entropy_floor, e->chemistry, ti_beg, ti_end,
            /*ti_start_mesh=*/-1, /*ti_end_mesh=*/-1);
}

/**
 * @brief Prepare an active particle for the force loop, as the ghost does.
 */
static void prepare_force(struct part *p, struct xpart *xp,
                          const struct engine *e) {
  const double dt_therm = get_timestep(p->time_bin, e->time_base);
  chemistry_prepare_force(p, xp, e->cosmology, dt_therm, dt_therm,
                          e->chemistry);
  timestep_limiter_prepare_force(p, xp);
}

/**
 * @brief Restore the physical state of a particle, keeping its time
 * variables.
 */
static void restore_state(struct part *p, const struct part *ref) {
  const timebin_t time_bin = p->time_bin;
  const struct timestep_limiter_data limiter_data = p->limiter_data;
  const struct chemistry_part_data chd = p->chemistry_data;
#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_kick = p->ti_kick, ti_drift = p->ti_drift;
#endif
  *p = *ref;
  p->time_bin = time_bin;
  p->limiter_data = limiter_data;
  p->chemistry_data.flux = chd.flux;
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_kick = ti_kick;
  p->ti_drift = ti_drift;
#endif
}

/**
 * @brief Run one schedule and measure the time covered by each pair.
 *
 * @param cover (return) Covered ticks of the pairs 0-1, 0-2 and 1-2.
 * @param residual (return) Relative residual of the fit.
 * @return 1 if the run kept the flux rates constant.
 */
static int run_schedule(const enum schedule sched, struct runner *r,
                        struct engine *e, double cover[3], double *residual) {

  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, 2 * sizeof(struct cell)) != 0)
    error("Could not allocate the cells.");
  bzero(c, 2 * sizeof(struct cell));
  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  if (posix_memalign((void **)&parts, part_align,
                     N_PARTS * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&xparts, xpart_align,
                     N_PARTS * sizeof(struct xpart)) != 0)
    error("Could not allocate the particles.");
  bzero(parts, N_PARTS * sizeof(struct part));
  bzero(xparts, N_PARTS * sizeof(struct xpart));

  /* Particle 0 in cell 0; particles 1 and 2 in cell 1. */
  const double pos[N_PARTS][3] = {
      {0.7, 0.5, 0.5}, {1.3, 0.5, 0.5}, {1.3, 0.8, 0.5}};
  for (int k = 0; k < N_PARTS; k++) {
    struct part *p = &parts[k];
    struct xpart *xp = &xparts[k];
    p->id = k + 1;
    for (int i = 0; i < 3; i++) p->x[i] = pos[k][i];
    p->h = 0.9 / kernel_gamma;
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
    p->conserved.mass = 1.f;
    p->conserved.energy = 1.f;
    hydro_first_init_part(p, xp);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) p->geometry.matrix_E[i][j] = (i == j);
    p->geometry.volume = 1.f;
    p->geometry.wcorr = 1.f;
    p->rho = 1.f;
    p->P = 1.f;
#else
    p->mass = 1.f;
    p->u = 1.f;
    p->rho = 1.f;
#endif
    p->geometry.condition_number = 1.f;
    hydro_prepare_force(p, xp, e->cosmology, e->hydro_properties,
                        e->pressure_floor_props, 0., 0.);
    hydro_reset_acceleration(p);

    struct chemistry_part_data *ch = &p->chemistry_data;
    chemistry_first_init_part(e->physical_constants, NULL, e->cosmology,
                              e->chemistry, p, xp);
    ch->limiter.maxr = 1.f;
    for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
      ch->limiter.rhoZ[m][0] = -1e10;
      ch->limiter.rhoZ[m][1] = 1e10;
      ch->limiter.Z[m][0] = -1e10;
      ch->limiter.Z[m][1] = 1e10;
      ch->metal_mass[m] = 1e-3;
    }
    p->time_bin = 0;
    p->limiter_data.wakeup = time_bin_not_awake;
    p->limiter_data.to_be_synchronized = 0;
#ifdef SWIFT_DEBUG_CHECKS
    p->ti_kick = 0;
    p->ti_drift = 0;
#endif
  }
  parts[0].chemistry_data.metal_mass[0] = 2e-3;
  parts[2].chemistry_data.metal_mass[1] = 2e-3;

  for (int k = 0; k < N_PARTS; k++) {
    struct chemistry_part_data *ch = &parts[k].chemistry_data;
    ch->gradients.rhoZ[0][0] = -1e-3 / 0.6;
    ch->gradients.Z[0][0] = -1e-3 / 0.6;
    ch->gradients.rhoZ[1][1] = 1e-3 / 0.3;
    ch->gradients.Z[1][1] = 1e-3 / 0.3;
  }

  /* The diffusion coefficient of an inactive particle is the one of its
     last ghost. */
  for (int k = 0; k < N_PARTS; k++)
    chemistry_prepare_force(&parts[k], &xparts[k], e->cosmology, 0., 0.,
                            e->chemistry);

  struct part ref[N_PARTS];
  memcpy(ref, parts, N_PARTS * sizeof(struct part));

  const int first_in_cell[2] = {0, 1};
  const int count_in_cell[2] = {1, 2};
  for (int n = 0; n < 2; n++) {
    struct cell *cn = &c[n];
    cn->loc[0] = n;
    cn->loc[1] = 0.;
    cn->loc[2] = 0.;
    for (int i = 0; i < 3; i++) cn->width[i] = 1.;
    cn->hydro.parts = parts + first_in_cell[n];
    cn->hydro.xparts = xparts + first_in_cell[n];
    cn->hydro.count = count_in_cell[n];
    cn->nodeID = e->nodeID;
    cn->top = cn;
    cn->hydro.super = cn;
    cn->dmin = 1.;
    cn->h_min_allowed = 0.5 * (1. / kernel_gamma);
    cn->h_max_allowed = 1. / kernel_gamma;
    for (int k = 0; k < cn->hydro.count; k++)
      cell_set_part_h_depth(&cn->hydro.parts[k], cn);
  }

  double delta[N_PARTS][2] = {{0.}};
  double rate[3][2] = {{0.}};
  int ok = 1;

  for (integertime_t ti = 0; ti <= TI_END; ti += 4) {

    e->ti_current = ti;
    e->max_active_bin = get_max_active_bin(ti);
    e->min_active_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
    for (int k = 0; k < N_PARTS; k++) parts[k].ti_drift = ti;
#endif

    /* The events at TI_EVENT happen in a step with no active particle. */
    if (ti == TI_EVENT) {
      for (int k = 0; k < N_PARTS; k++) {
        struct part *p = &parts[k];
        struct xpart *xp = &xparts[k];
        const int wake = (k == 0) && (sched == schedule_limiter ||
                                      sched == schedule_limiter_and_sync);
        const int sync = (k == 0 && sched == schedule_sync) ||
                         (k == 1 && sched == schedule_limiter_and_sync);
        if (sync) {
          timestep_process_sync_part(p, xp, e, e->cosmology);
          p->time_bin = min(BIN_4, e->max_active_bin);
        } else if (wake) {
          p->limiter_data.wakeup = 2 - BIN_4;
          timestep_limit_part(p, xp, e);
        }
        if (sync || wake) kick(p, xp, e, ti, ti + 2);
      }
      continue;
    }

    for (int k = 0; k < N_PARTS; k++) restore_state(&parts[k], &ref[k]);
    for (int n = 0; n < 2; n++) set_cell_properties(&c[n], e);
    if (!cell_is_active_hydro(&c[0], e) && !cell_is_active_hydro(&c[1], e))
      continue;

    for (int k = 0; k < N_PARTS; k++)
      if (part_is_active(&parts[k], e)) prepare_force(&parts[k], &xparts[k], e);

    /* The rate of each pair, from the state restored at the last step. */
    if (ti == TI_END) {
      for (int q = 0; q < 3; q++) {
        const struct part *pi = &parts[pair_i[q]];
        const struct part *pj = &parts[pair_j[q]];
        float dx[3];
        for (int i = 0; i < 3; i++) dx[i] = pi->x[i] - pj->x[i];
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4];
        float flux_dt = 0.f, vmax = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
        float delxbar_i = 0.f, delxbar_j = 0.f;
        const int status = chemistry_compute_pair_fluxes(
            r2, dx, pi->h, pj->h, pi, pj, e->time_base, ti, e->cosmology,
            /*with_cosmology=*/0, e->chemistry, totflux, &flux_dt, &vmax,
            &delxbar_i, &delxbar_j);
#else
        const int status = chemistry_compute_pair_fluxes(
            r2, dx, pi->h, pj->h, pi, pj, e->time_base, ti, e->cosmology,
            /*with_cosmology=*/0, e->chemistry, totflux, &flux_dt, &vmax);
#endif
        if (status != CHEMISTRY_PAIR_FLUX_FLUXES)
          error("Pair %d has no flux at the last step (status %d).", q, status);
        for (int m = 0; m < 2; m++) rate[q][m] = totflux[m][0];
      }
    }

    runner_dosub_self2_chemistry_flux_prep(r, &c[1], 0, 0);
    runner_dosub_pair2_chemistry_flux_prep(r, &c[0], &c[1], 0, 0);
    for (int n = 0; n < 2; n++) runner_do_chemistry_flux_ghost(r, &c[n], 0);
    for (int k = 0; k < N_PARTS; k++)
      for (int m = 0; m < 2; m++)
        if (parts[k].chemistry_data.fct_theta[m] != 1.0) ok = 0;

    runner_dosub_self2_force(r, &c[1], 0, 0);
    runner_dosub_pair2_force(r, &c[0], &c[1], 0, 0);
    for (int n = 0; n < 2; n++) runner_do_end_hydro_force(r, &c[n], 0);

    for (int k = 0; k < N_PARTS; k++)
      for (int m = 0; m < 2; m++)
        delta[k][m] += parts[k].chemistry_data.metal_mass[m] -
                       ref[k].chemistry_data.metal_mass[m];

    /* The time-step limiter is off: no neighbour wakes. */
    for (int k = 0; k < N_PARTS; k++) {
      struct part *p = &parts[k];
      p->limiter_data.wakeup = time_bin_not_awake;
      if (!part_is_active(p, e)) continue;
      const integertime_t dti_old = get_integer_timestep(p->time_bin);
      kick(p, &xparts[k], e, ti - dti_old / 2, ti);
      p->time_bin = next_bin(k, ti, p->time_bin);
      if (ti < TI_END)
        kick(p, &xparts[k], e, ti, ti + get_integer_timestep(p->time_bin) / 2);
    }
  }

  /* Least squares for delta[k][m] = tb * sum_q sign(k, q) rate[q][m] c_q,
     with sign -1 on the left particle of the pair and +1 on the right. */
  double M[3][3] = {{0.}}, b[3] = {0.};
  for (int m = 0; m < 2; m++) {
    for (int k = 0; k < N_PARTS; k++) {
      double row[3];
      for (int q = 0; q < 3; q++) {
        const double sign = (k == pair_i[q]) ? -1. : (k == pair_j[q]) ? 1. : 0.;
        row[q] = sign * rate[q][m] * e->time_base;
      }
      for (int q = 0; q < 3; q++) {
        b[q] += row[q] * delta[k][m];
        for (int s = 0; s < 3; s++) M[q][s] += row[q] * row[s];
      }
    }
  }
  const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
                     M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                     M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  for (int q = 0; q < 3; q++) {
    double Mq[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Mq[i][j] = (j == q) ? b[i] : M[i][j];
    cover[q] = (Mq[0][0] * (Mq[1][1] * Mq[2][2] - Mq[1][2] * Mq[2][1]) -
                Mq[0][1] * (Mq[1][0] * Mq[2][2] - Mq[1][2] * Mq[2][0]) +
                Mq[0][2] * (Mq[1][0] * Mq[2][1] - Mq[1][1] * Mq[2][0])) /
               det;
  }

  /* The fit must explain the moved mass: the flux is linear in time. */
  double res2 = 0., norm2 = 0.;
  for (int m = 0; m < 2; m++) {
    for (int k = 0; k < N_PARTS; k++) {
      double model = 0.;
      for (int q = 0; q < 3; q++) {
        const double sign = (k == pair_i[q]) ? -1. : (k == pair_j[q]) ? 1. : 0.;
        model += sign * rate[q][m] * e->time_base * cover[q];
      }
      res2 += (delta[k][m] - model) * (delta[k][m] - model);
      norm2 += delta[k][m] * delta[k][m];
    }
  }
  *residual = sqrt(res2 / norm2);

  for (int n = 0; n < 2; n++) cell_free_hydro_sorts(&c[n]);
  free(parts);
  free(xparts);
  free(c);
  return ok;
}

int main(int argc, char *argv[]) {

  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 0;
  space.dim[0] = 4.;
  space.dim[1] = 4.;
  space.dim[2] = 4.;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct hydro_props hydro_props;
  hydro_props_init_no_hydro(&hydro_props);
  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  struct entropy_floor_properties entropy_floor;
  bzero(&entropy_floor, sizeof(struct entropy_floor_properties));
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_vacuum_permeability = 1.;

  /* An inert rational limiter keeps the moved mass linear in the interval. */
  struct chemistry_global_data chem_data;
  bzero(&chem_data, sizeof(struct chemistry_global_data));
  chem_data.diffusion_coefficient = 1.f;
  chem_data.diffusion_mode = isotropic_constant;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  chem_data.tau = 1.;
  chem_data.relaxation_time_mode = constant_mode;
  chem_data.hyperbolic_limiter_scope = limiter_all_components;
#endif
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++)
    chem_data.initial_metallicities[m] = 0.f;
  chem_data.hll_riemann_solver_psi = 0.1f;
  chem_data.C_CFL_chemistry = 0.5f;
  chem_data.flux_limiter_noise_gate = 0.;
  chem_data.flux_limiter_safety = 1e30;
  chem_data.flux_limiter_sink_stability = 1e30;
  chem_data.flux_limiter_startup = 0.1;
  chem_data.use_positivity_limiter = 1;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time_base = 1. / (1 << 20);
  engine.nodeID = 0;
  engine.cosmology = &cosmo;
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  engine.entropy_floor = &entropy_floor;
  engine.physical_constants = &phys_const;
  engine.chemistry = &chem_data;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  int pass = 1;
  for (int s = 0; s < schedule_count; s++) {
    double cover[3], residual = 0.;
    const int ok =
        run_schedule((enum schedule)s, &runner, &engine, cover, &residual);
    int good = ok && residual <= 1e-9;
    for (int q = 0; q < 3; q++)
      good &= fabs(cover[q] - TI_END) <= 1e-5 * TI_END;
    message(
        "%-26s covered ticks: pair 0-1 %9.5f, pair 0-2 %9.5f, pair 1-2 "
        "%9.5f (bar %d), fit residual %.1e, theta == 1: %d -> %s",
        schedule_names[s], cover[0], cover[1], cover[2], TI_END, residual, ok,
        good ? "ok" : "WRONG");
    pass &= good;
  }

  if (!pass) error("testChemistryFluxInterval: FAILED");
  message("testChemistryFluxInterval: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
