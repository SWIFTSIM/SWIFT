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
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/* The Stage-1 LW/FUV artificial-dissipation term is GEAR-only physics and
 * lives in the hydro force loop. This test drives the real DOSELF2 force
 * dispatch (runner_doself2_branch_force), not the accumulate-band function
 * in isolation: the conservation property under test is a property of the
 * dispatch, and an isolated call cannot see it. */
#if defined(FEEDBACK_GEAR) && defined(SPHENIX_SPH)

#define TEST_NODE_ID 1

void runner_doself2_branch_force(struct runner *r, const struct cell *c,
                                 int limit_h_min, int limit_h_max);

/* A pair whose smoothing lengths differ enough that only one kernel
 * reaches the other particle. This is the geometry the relocation exists
 * for: the density-loop dispatch fires one side only there, leaving the
 * mirrored debit unwritten. */
static const float test_h_ratios[2] = {1.44f, 2.15f};

/* Particle state. Deliberately asymmetric in every input, so a mirrored
 * pair is the only way the mass-weighted sum can cancel. */
static const float test_mass_i = 1.0f;
static const float test_mass_j = 2.5f;
static const float test_rho_prev_i = 0.7f;
static const float test_rho_prev_j = 1.3f;
static const float test_u_FUV_i = 3.0f;
static const float test_u_FUV_j = -1.0f;
static const float test_u_LW_i = -0.4f;
static const float test_u_LW_j = 2.2f;
static const float test_c_hyp = 2.0f;
static const float test_alpha_pin = 0.3f;

/**
 * @brief Build a one-cell, two-particle setup with the requested h ratio.
 *
 * @param h_ratio Ratio h_i/h_j.
 * @param bulk_velocity Common velocity given to both particles.
 * @param cosmo The cosmology.
 * @param hydro_props The hydro properties.
 * @param pressure_floor The pressure-floor properties.
 * @return The cell, owning its own particle arrays.
 */
static struct cell *make_pair_cell(
    float h_ratio, float bulk_velocity, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  const double size = 40.;
  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate the cell");
  bzero(c, sizeof(struct cell));

  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     2 * sizeof(struct part)) != 0)
    error("Couldn't allocate the particles");
  bzero(c->hydro.parts, 2 * sizeof(struct part));
  if (posix_memalign((void **)&c->hydro.xparts, xpart_align,
                     2 * sizeof(struct xpart)) != 0)
    error("Couldn't allocate the extended particles");
  bzero(c->hydro.xparts, 2 * sizeof(struct xpart));

  const float hj = 1.0f;
  const float hi = h_ratio * hj;
  /* Half-way between the two kernel radii: inside i's reach, outside j's. */
  const float r = 0.5f * kernel_gamma * (hi + hj);

  struct part *pi = &c->hydro.parts[0];
  struct part *pj = &c->hydro.parts[1];

  pi->x[0] = 0.5 * size;
  pi->x[1] = 0.5 * size;
  pi->x[2] = 0.5 * size;
  pj->x[0] = pi->x[0] + (double)r;
  pj->x[1] = pi->x[1];
  pj->x[2] = pi->x[2];

  pi->h = hi;
  pj->h = hj;
  pi->mass = test_mass_i;
  pj->mass = test_mass_j;
  pi->u = 1.f;
  pj->u = 1.f;
  pi->id = 1;
  pj->id = 2;

  for (int k = 0; k < 2; k++) {
    struct part *p = &c->hydro.parts[k];
    p->time_bin = 1;
    p->depth_h = 0;
    p->v[0] = bulk_velocity;
    p->v[1] = bulk_velocity;
    p->v[2] = bulk_velocity;
    p->rho = 1.f;
    p->density.rho_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
    p->viscosity.alpha = 0.f;
    p->viscosity.div_v = 0.f;
    p->viscosity.div_v_previous_step = 0.f;
    p->viscosity.v_sig = hydro_get_comoving_soundspeed(p);
    p->force.pressure = hydro_get_comoving_pressure(p);
#ifdef SWIFT_DEBUG_CHECKS
    p->ti_drift = 8;
    p->ti_kick = 8;
#endif
    hydro_prepare_force(p, &c->hydro.xparts[k], cosmo, hydro_props,
                        pressure_floor, 0., 0.);
    hydro_reset_acceleration(p);
    mhd_init_part(p);
  }

  /* Radiation state. `dissipation_alpha_*` is pinned directly on the
   * particles rather than through the runtime parameter: the force-loop
   * interaction reads the field, never the parameter. */
  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  fdi->rho_prev = test_rho_prev_i;
  fdj->rho_prev = test_rho_prev_j;
  fdi->u_FUV = test_u_FUV_i;
  fdj->u_FUV = test_u_FUV_j;
  fdi->u_LW = test_u_LW_i;
  fdj->u_LW = test_u_LW_j;
  fdi->c_hyp = test_c_hyp;
  fdj->c_hyp = test_c_hyp;
  fdi->dissipation_alpha_FUV = test_alpha_pin;
  fdj->dissipation_alpha_FUV = test_alpha_pin;
  fdi->dissipation_alpha_LW = test_alpha_pin;
  fdj->dissipation_alpha_LW = test_alpha_pin;
  fdi->dissipation_u_FUV = 0.f;
  fdj->dissipation_u_FUV = 0.f;
  fdi->dissipation_u_LW = 0.f;
  fdj->dissipation_u_LW = 0.f;

  c->split = 0;
  c->depth = 0;
  c->hydro.count = 2;
  c->hydro.h_max = hi;
  c->hydro.h_max_active = hi;
  c->hydro.dx_max_part = 0.;
  c->hydro.dx_max_sort = 0.;
  c->width[0] = size;
  c->width[1] = size;
  c->width[2] = size;
  c->dmin = size;
  c->h_min_allowed = c->dmin * 0.5 * (1. / kernel_gamma);
  c->h_max_allowed = c->dmin * (1. / kernel_gamma);
  c->hydro.super = c;
  c->hydro.ti_old_part = 8;
  c->hydro.ti_end_min = 8;
  c->nodeID = TEST_NODE_ID;

  return c;
}

/**
 * @brief The expected pairwise exchange, computed in double from the same
 * inputs the C interaction reads.
 *
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param r Particle separation.
 * @param u_i Particle i's specific field (this band).
 * @param u_j Particle j's specific field (this band).
 * @return `Psi_ij`.
 */
static double expected_Psi(double hi, double hj, double r, double u_i,
                           double u_j) {

  float wi, wi_dx, wj, wj_dx;
  kernel_deval((float)(r / hi), &wi, &wi_dx);
  kernel_deval((float)(r / hj), &wj, &wj_dx);
  const double wi_dr = (double)wi_dx * pow_dimension_plus_one(1. / hi);
  const double wj_dr = (double)wj_dx * pow_dimension_plus_one(1. / hj);

  const double d_ij =
      (double)test_rho_prev_i * u_i - (double)test_rho_prev_j * u_j;
  const double Wbar = 0.5 * (wi_dr + wj_dr);
  const double v_sig = (double)test_alpha_pin * (double)test_c_hyp;
  return v_sig * d_ij * Wbar /
         ((double)test_rho_prev_i * (double)test_rho_prev_j);
}

/**
 * @brief Run one h ratio and check both assertions.
 *
 * @param h_ratio Ratio h_i/h_j.
 * @param bulk_velocity Common velocity given to both particles.
 * @param runner The runner.
 * @param cosmo The cosmology.
 * @param hydro_props The hydro properties.
 * @param pressure_floor The pressure-floor properties.
 * @param out (return) The two particles' `dissipation_u` values, FUV then
 * LW, for i then j.
 */
static void check_ratio(float h_ratio, float bulk_velocity,
                        struct runner *runner, const struct cosmology *cosmo,
                        const struct hydro_props *hydro_props,
                        const struct pressure_floor_props *pressure_floor,
                        float out[4]) {

  struct cell *c = make_pair_cell(h_ratio, bulk_velocity, cosmo, hydro_props,
                                  pressure_floor);

  const struct part *pi = &c->hydro.parts[0];
  const struct part *pj = &c->hydro.parts[1];
  const double dx = pj->x[0] - pi->x[0];
  const double r2 = dx * dx;
  const double hig2 = (double)pi->h * pi->h * kernel_gamma2;
  const double hjg2 = (double)pj->h * pj->h * kernel_gamma2;

  /* The whole point of the configuration: exactly one kernel reaches. A
   * mutually-reaching pair is conserved even under the old density-loop
   * dispatch, so a test run there would prove nothing. */
  if (!(r2 < hig2 && r2 >= hjg2))
    error(
        "h_ratio=%g: geometry is not asymmetric-reach (r2=%.8e, hig2=%.8e, "
        "hjg2=%.8e).",
        (double)h_ratio, r2, hig2, hjg2);

  runner_doself2_branch_force(runner, c, /*limit_h_min=*/0,
                              /*limit_h_max=*/0);

  const double mi = (double)test_mass_i;
  const double mj = (double)test_mass_j;
  const double r = sqrt(r2);

  const char *band_name[2] = {"FUV", "LW"};
  const double diss_i[2] = {(double)pi->feedback_data.dissipation_u_FUV,
                            (double)pi->feedback_data.dissipation_u_LW};
  const double diss_j[2] = {(double)pj->feedback_data.dissipation_u_FUV,
                            (double)pj->feedback_data.dissipation_u_LW};
  const double u_i[2] = {(double)test_u_FUV_i, (double)test_u_LW_i};
  const double u_j[2] = {(double)test_u_FUV_j, (double)test_u_LW_j};

  for (int b = 0; b < 2; b++) {

    const double credit = mi * diss_i[b];
    const double debit = mj * diss_j[b];
    const double sum = credit + debit;
    const double scale = fmax(fabs(credit), fabs(debit));

    const double Psi = expected_Psi((double)pi->h, (double)pj->h, r, u_i[b],
                                    u_j[b]);
    const double expected_credit = mi * mj * Psi;

    message(
        "h_ratio=%.3f band=%s: m_i*diss_i=%.10e m_j*diss_j=%.10e sum=%.10e "
        "|sum|/scale=%.3e expected m_i*diss_i=%.10e",
        (double)h_ratio, band_name[b], credit, debit, sum,
        (scale > 0. ? fabs(sum) / scale : 0.), expected_credit);

    /* (b) The term must actually have exchanged something. A guard that
     * silently skips the asymmetric-reach pair passes a sum-only check
     * trivially; this is what catches it. */
    if (scale == 0.)
      error("h_ratio=%g band=%s: no dissipation exchanged at all.",
            (double)h_ratio, band_name[b]);
    if (fabs(credit - expected_credit) > 1e-5 * fabs(expected_credit))
      error(
          "h_ratio=%g band=%s: exchanged magnitude %.10e differs from the "
          "expected %.10e (rel_err=%.3e).",
          (double)h_ratio, band_name[b], credit, expected_credit,
          fabs(credit - expected_credit) / fabs(expected_credit));

    /* (a) Exact conservation of `sum_i m_i*dissipation_u_i` to round-off.
     * Stated for a SAME-BIN ACTIVE-ACTIVE pair only: `doj` requires
     * PART_IS_ACTIVE, so an active-inactive pair stays one-sided under
     * this placement too (the separate, pre-existing cross-bin
     * asymmetry). */
    if (fabs(sum) > 1e-6 * scale)
      error(
          "h_ratio=%g band=%s: mass-weighted dissipation sum %.10e is not "
          "zero to round-off (scale=%.10e, ratio=%.3e).",
          (double)h_ratio, band_name[b], sum, scale, fabs(sum) / scale);
  }

  out[0] = pi->feedback_data.dissipation_u_FUV;
  out[1] = pj->feedback_data.dissipation_u_FUV;
  out[2] = pi->feedback_data.dissipation_u_LW;
  out[3] = pj->feedback_data.dissipation_u_LW;

  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c);
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
  struct gravity_props gravity_props;
  struct sink_props sink_props;
  struct chemistry_global_data chemistry;
  struct phys_const prog_const;
  struct unit_system us;
  struct runner *runner;

  bzero(&space, sizeof(struct space));
  bzero(&engine, sizeof(struct engine));
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  bzero(&gravity_props, sizeof(struct gravity_props));
  bzero(&sink_props, sizeof(struct sink_props));
  bzero(&chemistry, sizeof(struct chemistry_global_data));
  bzero(&prog_const, sizeof(struct phys_const));

  space.periodic = 0;
  space.dim[0] = 100.;
  space.dim[1] = 100.;
  space.dim[2] = 100.;

  units_init_cgs(&us);
  prog_const.const_vacuum_permeability = 1.0;
  gravity_props.G_Newton = 1.;

  engine.s = &space;
  engine.time = 0.1;
  engine.ti_current = 8;
  engine.time_base = 1e-10;
  engine.max_active_bin = num_time_bins;
  engine.nodeID = TEST_NODE_ID;
  engine.policy = engine_policy_hydro;
  engine.physical_constants = &prog_const;
  engine.internal_units = &us;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;
  hydro_props_init_no_hydro(&hydro_props);
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  engine.gravity_properties = &gravity_props;
  engine.sink_properties = &sink_props;
  engine.chemistry = &chemistry;

  if (posix_memalign((void **)&runner, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct runner)) != 0)
    error("Couldn't allocate the runner");
  bzero(runner, sizeof(struct runner));
  runner->e = &engine;

  float at_rest[2][4];
  for (int k = 0; k < 2; k++)
    check_ratio(test_h_ratios[k], /*bulk_velocity=*/0.f, runner, &cosmo,
                &hydro_props, &pressure_floor, at_rest[k]);

  /* Galilean spot check: the term reads positions, smoothing lengths,
   * densities, masses and `u` only, never a velocity, so a common bulk
   * velocity must leave it bitwise unchanged. */
  for (int k = 0; k < 2; k++) {
    float boosted[4];
    check_ratio(test_h_ratios[k], /*bulk_velocity=*/137.f, runner, &cosmo,
                &hydro_props, &pressure_floor, boosted);
    for (int m = 0; m < 4; m++) {
      if (boosted[m] != at_rest[k][m])
        error(
            "h_ratio=%g component %d: a common bulk velocity changed the "
            "dissipation from %.10e to %.10e.",
            (double)test_h_ratios[k], m, (double)at_rest[k][m],
            (double)boosted[m]);
    }
    message("h_ratio=%.3f: bitwise Galilean-invariant under a bulk boost.",
            (double)test_h_ratios[k]);
  }

  free(runner);
  message("All Stage-1 dissipation conservation checks passed.");
  return 0;
}

#else

int main(int argc, char *argv[]) {
  message(
      "Skipping: the Stage-1 LW/FUV dissipation term needs GEAR feedback "
      "and the SPHENIX hydro scheme.");
  return 0;
}

#endif /* FEEDBACK_GEAR && SPHENIX_SPH */
