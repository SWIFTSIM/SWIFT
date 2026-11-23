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
 * @file testChemistryFluxStamp.c
 * @brief Check that only the end of the force step moves the time stamp of
 * the last flux application of an FVPM particle.
 *
 * The kicks with a negative and then a positive time-step done by the
 * time-step limiter and by the synchronisation must leave the stamp
 * unchanged. chemistry_end_force() must set it to the current time.
 *
 * Without a limiter or a sync event, the flux interval of a pair computed from
 * the stamps must be bit for bit the interval of the ghost-based scheme it
 * replaces: the smaller of the two particles' step intervals, each rounded to
 * float and then multiplied by a^2. This is checked with and without
 * cosmology.
 */

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

/* A step of 32 ticks from 0, interrupted at ti = 20. */
#define BIN_32 4
#define BIN_4 1
#define TI_EVENT 20
#define TI_STAMP 8

/**
 * @brief Set up a particle on its step [0, 32], kicked to the middle.
 */
static void init_particle(struct part *p, struct xpart *xp,
                          const struct engine *e) {
  bzero(p, sizeof(struct part));
  bzero(xp, sizeof(struct xpart));
  p->id = 1;
  p->h = 1.f;
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
  p->conserved.mass = 1.f;
  p->conserved.energy = 1.f;
  hydro_first_init_part(p, xp);
  p->rho = 1.f;
  p->P = 1.f;
#else
  p->mass = 1.f;
  p->u = 1.f;
  p->rho = 1.f;
#endif
  p->geometry.volume = 1.f;
  chemistry_first_init_part(e->physical_constants, NULL, e->cosmology,
                            e->chemistry, p, xp);
  p->chemistry_data.flux.ti_flux_last = TI_STAMP;
  p->time_bin = BIN_32;
  p->limiter_data.wakeup = time_bin_not_awake;
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_kick = get_integer_timestep(BIN_32) / 2;
  p->ti_drift = TI_EVENT;
#endif
}

/**
 * @brief The flux interval of a particle in the ghost-based scheme.
 *
 * The ghost passed the step interval to chemistry_prepare_force() as a float,
 * which multiplied it by a^2.
 */
static float ghost_flux_dt(const integertime_t ti_begin,
                           const integertime_t ti_end, const double time_base,
                           const struct cosmology *cosmo,
                           const int with_cosmology) {
  const float dt_therm =
      with_cosmology ? cosmology_get_therm_kick_factor(cosmo, ti_begin, ti_end)
                     : (ti_end - ti_begin) * time_base;
  return dt_therm * cosmo->a * cosmo->a;
}

/**
 * @brief Count the pairs whose stamp-based flux interval is not bit for bit
 * the ghost-based one, over pairs of steps of several lengths ending at
 * several times.
 */
static int count_interval_mismatches(struct cosmology *cosmo,
                                     const struct phys_const *phys_const,
                                     const double time_base,
                                     const int with_cosmology,
                                     const struct engine *e) {

  /* Steps from 2^6 to 2^51 ticks, ending at times up to the end of the time
     line, so a covers the whole range of the cosmology. */
  const timebin_t bins[5] = {5, 17, 29, 41, 50};
  const integertime_t ends[5] = {1, 7, 23, 40, 64};
  const integertime_t longest = get_integer_timestep(bins[4]);

  struct part p[2];
  struct xpart xp[2];
  for (int k = 0; k < 2; k++) {
    init_particle(&p[k], &xp[k], e);
    p[k].x[0] = 0.5 * k;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) p[k].geometry.matrix_E[i][j] = (i == j);
    p[k].geometry.wcorr = 1.f;
    p[k].geometry.condition_number = 1.f;
    p[k].chemistry_data.kappa = 1.f;
  }
  float dx[3] = {p[0].x[0] - p[1].x[0], 0.f, 0.f};
  const float r2 = dx[0] * dx[0];

  int mismatches = 0;
  for (int n = 0; n < 5; n++) {
    const integertime_t ti = ends[n] * longest;
    if (with_cosmology) cosmology_update(cosmo, phys_const, ti);
    for (int s = 0; s < 5; s++) {
      for (int l = s; l < 5; l++) {
        const integertime_t begin_long = ti - get_integer_timestep(bins[l]);
        const integertime_t begin_short = ti - get_integer_timestep(bins[s]);
        p[0].chemistry_data.flux.ti_flux_last = begin_long;
        p[1].chemistry_data.flux.ti_flux_last = begin_short;

        double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4];
        float flux_dt = 0.f, vmax = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
        float delxbar_i = 0.f, delxbar_j = 0.f;
        chemistry_compute_pair_fluxes(r2, dx, p[0].h, p[1].h, &p[0], &p[1],
                                      time_base, ti, cosmo, with_cosmology,
                                      e->chemistry, totflux, &flux_dt, &vmax,
                                      &delxbar_i, &delxbar_j);
#else
        chemistry_compute_pair_fluxes(r2, dx, p[0].h, p[1].h, &p[0], &p[1],
                                      time_base, ti, cosmo, with_cosmology,
                                      e->chemistry, totflux, &flux_dt, &vmax);
#endif
        const float ghost = fminf(
            ghost_flux_dt(begin_long, ti, time_base, cosmo, with_cosmology),
            ghost_flux_dt(begin_short, ti, time_base, cosmo, with_cosmology));
        if (memcmp(&flux_dt, &ghost, sizeof(float)) != 0) {
          message("a=%.6f steps 2^%d, 2^%d ticks: interval %.9e, ghost %.9e",
                  cosmo->a, bins[l] + 1, bins[s] + 1, flux_dt, ghost);
          mismatches++;
        }
      }
    }
  }
  return mismatches;
}

int main(int argc, char *argv[]) {

  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct hydro_props hydro_props;
  hydro_props_init_no_hydro(&hydro_props);
  struct entropy_floor_properties entropy_floor;
  bzero(&entropy_floor, sizeof(struct entropy_floor_properties));
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  struct chemistry_global_data chem_data;
  bzero(&chem_data, sizeof(struct chemistry_global_data));
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  chem_data.tau = 1.;
  chem_data.relaxation_time_mode = constant_mode;
#endif

  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.time_base = 1e-3;
  e.ti_current = TI_EVENT;
  e.max_active_bin = get_max_active_bin(TI_EVENT);
  e.min_active_bin = 1;
  e.cosmology = &cosmo;
  e.hydro_properties = &hydro_props;
  e.entropy_floor = &entropy_floor;
  e.physical_constants = &phys_const;
  e.chemistry = &chem_data;

  struct part p;
  struct xpart xp;
  int pass = 1;

  /* The limiter kicks back by the old half step, then forward to now. */
  init_particle(&p, &xp, &e);
  p.limiter_data.wakeup = 2 - BIN_4;
  timestep_limit_part(&p, &xp, &e);
  const integertime_t after_limiter = p.chemistry_data.flux.ti_flux_last;
  pass &= (after_limiter == TI_STAMP);

  /* The sync kicks back by the old half step, then forward to now. */
  init_particle(&p, &xp, &e);
  timestep_process_sync_part(&p, &xp, &e, &cosmo);
  const integertime_t after_sync = p.chemistry_data.flux.ti_flux_last;
  pass &= (after_sync == TI_STAMP);

  /* The end of the force step moves the stamp to now. */
  chemistry_end_force(&p, &cosmo, /*with_cosmology=*/0, /*time=*/0.,
                      /*dt=*/0., e.ti_current, &chem_data);
  const integertime_t after_end_force = p.chemistry_data.flux.ti_flux_last;
  pass &= (after_end_force == TI_EVENT);

  message(
      "Stamp after the limiter kicks %lld, after the sync kicks %lld (both "
      "must be %d), after end_force %lld (must be %d)",
      after_limiter, after_sync, TI_STAMP, after_end_force, TI_EVENT);

  /* A time base that is not a power of two, so the products round. */
  const int no_cosmo_mismatches = count_interval_mismatches(
      &cosmo, &phys_const, 0.1 / max_nr_timesteps, /*with_cosmology=*/0, &e);

  struct unit_system us;
  units_init_cgs(&us);
  struct swift_params params;
  parser_init("", &params);
  parser_set_param(&params, "Cosmology:Omega_cdm:0.2589");
  parser_set_param(&params, "Cosmology:Omega_lambda:0.6910");
  parser_set_param(&params, "Cosmology:Omega_b:0.0486");
  parser_set_param(&params, "Cosmology:h:0.6774");
  parser_set_param(&params, "Cosmology:a_begin:0.0078125");
  parser_set_param(&params, "Cosmology:a_end:1.0");
  struct phys_const cosmo_const;
  phys_const_init(&us, &params, &cosmo_const);
  struct cosmology cosmo_on;
  cosmology_init(&params, &us, &cosmo_const, &cosmo_on);
  const int cosmo_mismatches = count_interval_mismatches(
      &cosmo_on, &cosmo_const, cosmo_on.time_base, /*with_cosmology=*/1, &e);
  cosmology_clean(&cosmo_on);

  message(
      "Pair intervals differing from the ghost-based ones: %d without "
      "cosmology, %d with cosmology (of 75 each, must be 0)",
      no_cosmo_mismatches, cosmo_mismatches);
  pass &= (no_cosmo_mismatches == 0 && cosmo_mismatches == 0);

  if (!pass) error("testChemistryFluxStamp: FAILED");
  message("testChemistryFluxStamp: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
