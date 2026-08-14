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

#ifdef GEAR_COOLING

/**
 * @brief Regression test for the ionized-gas temperature-pin energy
 * double-application bug (cooling_gear_subgrid.h).
 *
 * cooling_ionize_part_subgrid() used to write the pinned energy straight to
 * xp->u_full (via hydro_set_physical_internal_energy) AND zero p->u_dt, while
 * cooling_new_energy() then returned that same just-written value. Back in
 * cooling_cool_part(), that return value was turned into a du/dt relative to
 * the PRE-pin energy and hydro_kick_extra() applied it on top of the energy
 * already written -- landing the pinned energy twice
 * (u_end = 2*u_pin - u_ad_before instead of u_pin).
 *
 * The fix makes the pin path return u_pin without touching the particle's
 * energy state, so the existing du/dt machinery (shared with every other
 * cooling_new_energy() return path) lands it exactly once. This test drives
 * the real cooling_update_part_subgrid(), then replays cooling_cool_part's
 * own du/dt bookkeeping (cooling.c) and the real hydro_kick_extra() kick, and
 * checks the particle ends up AT u_pin, not overshot.
 */
static void test_ionized_pin_lands_once(void) {

  struct unit_system us;
  units_init_cgs(&us);

  struct phys_const phys_const;
  bzero(&phys_const, sizeof(phys_const));
  phys_const.const_boltzmann_k = 1.380649e-16;      /* cgs */
  phys_const.const_proton_mass = 1.67262192369e-24; /* cgs */

  struct cosmology cosmo;
  bzero(&cosmo, sizeof(cosmo));
  cosmo.a_factor_internal_energy = 1.0;

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
  const float u_old = 1.0e8f; /* erg/g -- orders of magnitude below u_pin. */
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

  float u_pin = 0.f;
  const int ionized_this_step = cooling_update_part_subgrid(
      &phys_const, &us, &cosmo, &hydro_props, &pressure_floor, &cooling, &p,
      &xp, dt, dt_therm, time, &u_pin);

  if (!ionized_this_step)
    error(
        "cooling_update_part_subgrid did not report the particle as "
        "ionized -- test setup error.");

  if (u_pin <= u_old)
    error(
        "Test setup error: u_pin (%e) must sit well above u_old (%e) for "
        "this test to discriminate the double-application bug.",
        u_pin, u_old);

  /* The fix's actual contract: the pin call must not have touched the
     particle's energy state directly -- only returned u_pin. */
  if (xp.u_full != u_old)
    error(
        "cooling_update_part_subgrid wrote xp->u_full directly (got %e, "
        "expected it untouched at %e) -- the pin must return u_pin instead "
        "of writing it, or the caller's own du/dt integration below "
        "double-applies it.",
        xp.u_full, u_old);

  /* Mirrors cooling_cool_part's post-cooling bookkeeping
     (cooling.c:1158-1169): derive du/dt from the returned pin energy and
     land it through the same machinery every other cooling path uses. */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(&p, &cosmo);
  const float cool_du_dt = (u_pin - u_ad_before) / dt_therm;
  const float du_dt = cool_du_dt + hydro_du_dt;
  hydro_set_physical_internal_energy_dt(&p, &cosmo, du_dt);

  /* Run the real kick, not a re-implementation of it. */
  hydro_kick_extra(&p, &xp, dt_therm, /*dt_grav=*/0.f, /*dt_grav_mesh=*/0.f,
                   /*dt_hydro=*/dt_therm, /*dt_kick_corr=*/0.f, &cosmo,
                   &hydro_props, &floor_props);

  const double rel_err = fabs((double)xp.u_full - (double)u_pin) / u_pin;
  message("u_old=%e u_pin=%e u_full_after_kick=%e rel_err=%e", (double)u_old,
          (double)u_pin, (double)xp.u_full, rel_err);

  if (rel_err > 1e-5)
    error(
        "Energy double-application regression: expected xp->u_full == "
        "u_pin (%e) after cooling+kick, got %e (rel_err=%e). The "
        "ionized-gas pin must land its energy exactly once.",
        u_pin, xp.u_full, rel_err);
}

#endif /* GEAR_COOLING */

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

#ifdef GEAR_COOLING
  test_ionized_pin_lands_once();
  message("testCoolingIonizedPin: PASS.");
#else
  message("testCoolingIonizedPin: skipped (not built with GEAR_COOLING).");
#endif

  return 0;
}
