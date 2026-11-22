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

  if (!pass) error("testChemistryFluxStamp: FAILED");
  message("testChemistryFluxStamp: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
