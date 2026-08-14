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
#include <string.h>

/* Local headers. */
#include "swift.h"

/* IONIZATION_FEEDBACK_DEBUG_NO_COOLING compiles cooling_ionize_part_subgrid
   down to an unconditional "return 0" (cooling_gear_subgrid.h): under that
   flag the ionized-gas temperature floor this test targets never runs at
   all, by design (gas heated by feedback simply never cools back down, so
   the per-step re-floor this test exercises has nothing to do). Nothing to
   test in that configuration. */
#if defined(GEAR_COOLING) && !defined(IONIZATION_FEEDBACK_DEBUG_NO_COOLING)

/**
 * @brief Regression test: cooling_update_part_subgrid() must land the
 * ionized-gas temperature floor's energy exactly once.
 *
 * cooling_ionize_part_subgrid() returns the floor energy through an output
 * parameter without touching the particle's energy state directly, so the
 * du/dt machinery shared by every other cooling_new_energy() return path
 * lands it exactly once. This test drives the real
 * cooling_update_part_subgrid(), then replays cooling_cool_part's own du/dt
 * bookkeeping (cooling.c) and the real hydro_kick_extra() kick, and checks
 * the particle ends up AT u_floor, not overshot.
 */
static void test_ionized_temperature_floor_lands_once(void) {

  struct unit_system us;
  units_init_cgs(&us);

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(phys_const));
  phys_const.const_boltzmann_k = 1.380649e-16;      /* cgs */
  phys_const.const_proton_mass = 1.67262192369e-24; /* cgs */

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(cosmo));
  cosmo.a_factor_internal_energy = 1.0;
  /* Present-epoch (a=1) scale-factor state. a3_inv matters here beyond
     a_factor_internal_energy: hydro_get_physical_density() (SPHENIX
     hydro.h) returns cosmo->a3_inv * p->rho, and at COOLING_GRACKLE_MODE
     >= 1 cooling_get_mean_molecular_weight() computes every species number
     density from that physical density. Left at bzero's 0, the density is
     always 0 regardless of p.rho, so mu is a 0/0 that only stayed finite
     by -ffinite-math-only's instruction-selection luck -- not by design. */
  cosmo.a = 1.0;
  cosmo.a_inv = 1.0;
  cosmo.a3_inv = 1.0;

  struct hydro_props hydro_props;
  bzero(&hydro_props, sizeof(hydro_props));
  hydro_props.hydrogen_ionization_temperature = 1.0e4;
  hydro_props.mu_neutral = 1.22;
  hydro_props.mu_ionised = 0.59;
  hydro_props.minimal_internal_energy = 0.f;

  struct entropy_floor_properties floor_props;
  bzero(&floor_props, sizeof(floor_props));

  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(pressure_floor));

  struct cooling_function_data cooling;
  bzero(&cooling, sizeof(cooling));
  cooling.HydrogenFractionByMass = 0.75f;
  cooling.HII_couple_ionization_rate = 0;

  struct part p;
  struct xpart xp;
  bzero(&p, sizeof(p));
  bzero(&xp, sizeof(xp));

  p.mass = 1.0f;              /* Cancels out of the ionized-energy formula. */
  p.rho = 1.67262192369e-24f; /* ~1 particle/cm^3, cgs -- only read at
                                  COOLING_GRACKLE_MODE >= 1. */
#if COOLING_GRACKLE_MODE > 0
  /* At mode >= 1, cooling_get_mean_molecular_weight reads these species
     fractions instead of HydrogenFractionByMass; left at their zero-init
     they make mu a 0/0 that only stays finite by -ffinite-math-only's
     instruction-selection luck. Physical, X_H-consistent values instead. */
  xp.cooling_data.HI_frac = 0.75f;
  xp.cooling_data.HeI_frac = 0.25f;
#endif
  const float u_old = 1.0e8f; /* erg/g -- orders of magnitude below u_floor. */
  p.u = u_old;
  xp.u_full = u_old;
  p.u_dt = 1.0e5f; /* Nonzero adiabatic derivative, to exercise the
                       cancellation the fix relies on. */

  p.feedback_data.is_ionized = 1;
  p.feedback_data.star_id = 1;
  p.feedback_data.end_time = 1e30; /* Far future: tag must not lapse here. */

  const double time = 0.0;
  const double dt = 10.0;
  const double dt_therm = 10.0;

  /* Mirrors cooling_cool_part's pre-cooling adiabatic estimate
     (cooling.c:1132-1136). */
  const float u_ad_before =
      u_old + dt_therm * hydro_get_physical_internal_energy_dt(&p, &cosmo);

  float u_floor = 0.f;
  const int ionized_this_step = cooling_update_part_subgrid(
      &phys_const, &us, &cosmo, &hydro_props, &pressure_floor, &cooling, &p,
      &xp, dt, dt_therm, time, &u_floor);

  if (!ionized_this_step)
    error(
        "cooling_update_part_subgrid did not report the particle as "
        "ionized -- test setup error.");

  if (u_floor <= u_old)
    error(
        "Test setup error: u_floor (%e) must sit well above u_old (%e) for "
        "this test to detect a doubly-applied floor energy.",
        u_floor, u_old);

  /* The fix's actual contract: enforcing T_HII must not touch the
     particle's energy state directly -- only return u_floor. */
  if (xp.u_full != u_old)
    error(
        "cooling_update_part_subgrid wrote xp->u_full directly (got %e, "
        "expected it untouched at %e) -- enforcing T_HII must return "
        "u_floor instead of writing it, or the caller's own du/dt "
        "integration below double-applies it.",
        xp.u_full, u_old);

  /* Mirrors cooling_cool_part's post-cooling bookkeeping
     (cooling.c:1158-1169): derive du/dt from the returned floor energy and
     land it through the same machinery every other cooling path uses. */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(&p, &cosmo);
  const float cool_du_dt = (u_floor - u_ad_before) / dt_therm;
  const float du_dt = cool_du_dt + hydro_du_dt;
  hydro_set_physical_internal_energy_dt(&p, &cosmo, du_dt);

  /* Run the real kick, not a re-implementation of it. */
  hydro_kick_extra(&p, &xp, dt_therm, /*dt_grav=*/0.f, /*dt_grav_mesh=*/0.f,
                   /*dt_hydro=*/dt_therm, /*dt_kick_corr=*/0.f, &cosmo,
                   &hydro_props, &floor_props);

  const double rel_err = fabs((double)xp.u_full - (double)u_floor) / u_floor;
  message("u_old=%e u_floor=%e u_full_after_kick=%e rel_err=%e", (double)u_old,
          (double)u_floor, (double)xp.u_full, rel_err);

  if (rel_err > 1e-5)
    error(
        "Energy double-application regression: expected xp->u_full == "
        "u_floor (%e) after cooling+kick, got %e (rel_err=%e). The "
        "ionized-gas temperature floor must land its energy exactly once.",
        u_floor, xp.u_full, rel_err);
}

#endif /* GEAR_COOLING && !IONIZATION_FEEDBACK_DEBUG_NO_COOLING */

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

#if defined(GEAR_COOLING) && !defined(IONIZATION_FEEDBACK_DEBUG_NO_COOLING)
  test_ionized_temperature_floor_lands_once();
  message("testCoolingIonizedTemperatureFloor: PASS.");
#elif !defined(GEAR_COOLING)
  message(
      "testCoolingIonizedTemperatureFloor: skipped (not built with "
      "GEAR_COOLING).");
#else
  message(
      "testCoolingIonizedTemperatureFloor: skipped "
      "(IONIZATION_FEEDBACK_DEBUG_NO_COOLING disables the floor "
      "enforcement this test targets).");
#endif

  return 0;
}
