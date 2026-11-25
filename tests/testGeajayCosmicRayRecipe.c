/*******************************************************************************
 * This file is part of SWIFT.
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

/* Depletion recipe of geajay_cooling.rst under fast math. Run by
 * testGeajayCosmicRayRecipe.sh. */

#define _GNU_SOURCE

/* Config parameters. */
#include <config.h>

#ifndef __FAST_MATH__
#error "This test must be built with -ffast-math."
#endif

/* Local headers. */
#include "swift.h"

#include <dlfcn.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geajay.h"
#include "geajay_particle.h"

#ifndef COOLING_GEAJAY_MODE
#error "This test needs --with-cooling=geajay_3 or geajay_4."
#endif

/* Library entry of the reduced cosmic-ray solve. With the variable FAIL set it
 * reports a failed solve, as the library does when the iterations run out. */
int geajay_particle_solve_cosmic_ray_reduced(
    const geajay_cosmic_ray_reduced_context *context,
    const geajay_species_layout *layout, void *state,
    const geajay_particle_inputs *inputs, double v, double dt, int n,
    geajay_cosmic_ray_output *output, geajay_particle_status *status) {

  static int (*real)(const geajay_cosmic_ray_reduced_context *,
                     const geajay_species_layout *, void *,
                     const geajay_particle_inputs *, double, double, int,
                     geajay_cosmic_ray_output *, geajay_particle_status *);
  if (!real)
    real = dlsym(RTLD_NEXT, "geajay_particle_solve_cosmic_ray_reduced");
  if (getenv("FAIL")) {
    status->reason = GEAJAY_PARTICLE_ITERATIONS_EXCEEDED;
    return GEAJAY_FAIL;
  }
  return real(context, layout, state, inputs, v, dt, n, output, status);
}

/* The recipe of geajay_cooling.rst, verbatim. */
static void deplete(struct xpart *xp, const double dt_s) {
  if (cooling_cosmic_ray_source_is_per_particle()) {
    const double loss = cooling_get_cosmic_ray_energy_loss_rate(xp);
    const double e = cooling_get_cosmic_ray_energy_density(xp);
    cooling_set_cosmic_ray_energy_density(xp, fmax(0., e - loss * dt_s));
  }
}

static struct swift_params params;
static struct unit_system us;
static struct phys_const pc;
static struct hydro_props hp;
static struct cooling_function_data cooling;

/* Usage: testGeajayCosmicRayRecipe <parameter file> <dt_s> [loss]. A third
 * argument sets the stored loss rate by hand (a finite control). */
int main(int argc, char **argv) {

  if (argc < 3) error("Usage: %s <parameter file> <dt_s> [loss]", argv[0]);
  const double dt_s = atof(argv[2]);

  parser_read_file(argv[1], &params);
  units_init_from_params(&us, &params, "InternalUnitSystem");
  phys_const_init(&us, &params, &pc);
  hydro_props_init(&hp, &pc, &us, &params);
  cooling_init(&params, &us, &pc, &hp, &cooling);

  struct part p;
  struct xpart xp;
  struct cosmology cosmo;
  struct pressure_floor_props pf;
  struct entropy_floor_properties ef;
  memset(&p, 0, sizeof(p));
  memset(&xp, 0, sizeof(xp));
  cosmology_init_no_cosmo(&cosmo);
  pressure_floor_init(&pf, &pc, &us, &hp, &params);
  entropy_floor_init(&ef, &pc, &us, &hp, &params);
  cooling_update(&pc, &cosmo, &pf, &cooling, NULL, 0.);

  p.id = 1;
  p.mass = 1.;
  p.h = 0.1f;
  p.rho = 3.3e5f;
  p.u = 1e3f;
  p.chemistry_data.smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT -
                                                1] = 0.01f;
  cooling_first_init_part(&pc, &us, &hp, NULL, &cooling, &p, &xp);
  cooling_set_cosmic_ray_energy_density(&xp, 1.0);

  const double dt_cool = 1e-6;
  cooling_cool_part(&pc, &us, &cosmo, &hp, &ef, &pf, &cooling, &p, &xp,
                    dt_cool, dt_cool, 0.);

  if (argc > 3) xp.cooling_data.cosmic_ray_energy_loss_rate = atof(argv[3]);

  const double e_before = cooling_get_cosmic_ray_energy_density(&xp);
  const double loss = cooling_get_cosmic_ray_energy_loss_rate(&xp);
  deplete(&xp, dt_s);
  const double e_after = cooling_get_cosmic_ray_energy_density(&xp);

  printf("FASTMATH_TEST %s\n", compilation_cflags());
  printf("RECIPE before=%.9g getter=%.9g after=%.9g dt=%.9g\n", e_before,
         loss, e_after, dt_s);
  return 0;
}
