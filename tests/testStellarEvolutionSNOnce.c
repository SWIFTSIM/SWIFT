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
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Build a two-point constant #interpolation_1d, so every lookup
 * returns `value` regardless of its argument.
 *
 * @param data Backing storage for the table (must outlive `interp`).
 * @param value The constant value every lookup must return.
 * @param interp (return) The interpolation table.
 */
static void make_constant_table(float data[2], float value,
                                struct interpolation_1d *interp) {
  data[0] = value;
  data[1] = value;
  interp->data = data;
  interp->xmin = 0.0f;
  interp->dx = 1.0f;
  interp->N = 2;
  interp->boundary_condition = boundary_condition_const;
}

/**
 * @brief A single_star's one core-collapse SN must inject exactly once, ever,
 * even when it is re-entered several times with a stale star_age_beg_step
 * that keeps reading as "not dead yet" (see stellar_evolution.c's is_dead
 * latch: star_age_beg_step is rebuilt from sp->time_bin every call, and the
 * observed HomogeneousBox run shows this rebuilt age can get pinned below
 * lifetime_myr for several calls in a row after the true death-crossing
 * step, while the step's own dt keeps doubling).
 */
int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  /* ------------------------------------------------------------------ */
  /* Build a minimal stellar_model: constant lifetime (independent of mass
     and metallicity), no radiation/winds, and flat SNII yield/energy
     tables so the real injection code path runs without needing any HDF5
     table. */
  struct stellar_model sm;
  bzero(&sm, sizeof(struct stellar_model));

  const float lifetime_myr = 10.0f;
  sm.lifetime.constant[2] = log10f(lifetime_myr);

  sm.rad.is_active = 0;
  sm.discrete_yields = 1;
  sm.discrete_star_minimal_gravity_mass = 1e-3f;
  sm.imf.mass_min = 0.1f;
  sm.imf.mass_max = 100.0f;

  float energy_table[2], mass_processed_table[2], mass_non_processed_table[2];
  make_constant_table(energy_table, 1.0f, &sm.snii.energy_per_progenitor_mass);
  make_constant_table(mass_processed_table, 0.5f,
                      &sm.snii.raw.ejected_mass_processed);
  make_constant_table(mass_non_processed_table, 0.0f,
                      &sm.snii.raw.ejected_mass_non_processed);

  float yield_tables[GEAR_CHEMISTRY_ELEMENT_COUNT][2];
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    make_constant_table(yield_tables[i], 0.0f, &sm.snii.raw.yields[i]);
  }

  /* ------------------------------------------------------------------ */
  struct spart sp;
  bzero(&sp, sizeof(struct spart));
  sp.id = 1;
  sp.star_type = single_star;
  sp.sf_data.birth_mass = 10.0;
  sp.mass = 10.0;

  struct unit_system us;
  units_init_cgs(&us);

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_solar_mass = 1.0;
  phys_const.const_year = 1.0 / 1e6; /* so 1 Myr == 1 internal time unit */

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(struct cosmology));
  cosmo.a = 1.0;
  cosmo.a3_inv = 1.0;

  const double conversion_to_myr = phys_const.const_year * 1e6;

  /* ------------------------------------------------------------------ */
  /* Step 1: the real death-crossing step. beg is just below lifetime_myr,
     end just above it: this must fire the star's one and only SN. */
  const double beg0 = (lifetime_myr - 0.01) * conversion_to_myr;
  double dt = 0.02 * conversion_to_myr;
  stellar_evolution_evolve_individual_star(&sp, &sm, /*with_cosmology=*/0,
                                           &cosmo, /*time=*/0.0, &us,
                                           &phys_const,
                                           /*with_stellar_wind_feedback=*/0,
                                           /*ti_begin=*/0, beg0, dt);

  if (sp.feedback_data.number_snii != 1)
    error("Step 1 should have fired the star's SN (number_snii=%g).",
          sp.feedback_data.number_snii);
  if (!sp.feedback_data.is_dead)
    error("Step 1 should have latched the star dead.");
  if (sp.tracers_data.snii_events.n_events != 1.0f)
    error("Step 1 should record exactly one SN event, got %g.",
          sp.tracers_data.snii_events.n_events);

  const double mass_after_first_explosion = sp.mass;
  if (mass_after_first_explosion >= 10.0)
    error("Step 1 should have ejected mass (mass=%g).",
          mass_after_first_explosion);

  /* Snapshot every field a repeat call could touch, to check afterwards
     that the repeat calls below were complete no-ops rather than assuming
     any one field in particular gets reset between calls. memcpy, not a
     struct assignment, so padding bytes are copied too and memcmp below
     is not comparing indeterminate padding. */
  struct feedback_spart_data snapshot;
  memcpy(&snapshot, &sp.feedback_data, sizeof(snapshot));

  /* ------------------------------------------------------------------ */
  /* Steps 2-5: replay the pattern measured in the HomogeneousBox log for a
     repeat-firing star, where star_age_beg_step stays PINNED at the same
     value as step 1 (not advancing) while dt keeps doubling. Without the
     one-shot latch, star_age_beg_step_myr stays below lifetime_myr on
     every one of these calls, so the age-only gate alone would fire again
     each time. */
  for (int i = 0; i < 4; i++) {
    dt *= 2.0;
    stellar_evolution_evolve_individual_star(&sp, &sm, /*with_cosmology=*/0,
                                             &cosmo, /*time=*/0.0, &us,
                                             &phys_const,
                                             /*with_stellar_wind_feedback=*/0,
                                             /*ti_begin=*/0, beg0, dt);
  }

  /* ------------------------------------------------------------------ */
  /* A single star has exactly one core-collapse SN, ever: the repeat calls
     above must have been complete no-ops. */
  if (sp.tracers_data.snii_events.n_events != 1.0f)
    error(
        "single_star fired %g SN events across the repeat calls; expected "
        "exactly 1 (the fix must latch the star dead at injection and gate "
        "on that latch, not on a fresh age recomputation).",
        sp.tracers_data.snii_events.n_events);

  if (sp.mass != mass_after_first_explosion)
    error(
        "Star mass changed on a repeat call (%g -> %g): SN ejecta was "
        "injected more than once.",
        mass_after_first_explosion, sp.mass);

  if (memcmp(&sp.feedback_data, &snapshot, sizeof(snapshot)) != 0)
    error(
        "feedback_data changed across the repeat (latched) calls: a "
        "dead star must not recompute or re-inject anything.");

  message("Single star fired its supernova exactly once, as expected.");
  return 0;
}
