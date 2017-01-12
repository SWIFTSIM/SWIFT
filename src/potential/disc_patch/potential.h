/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DISC_PATCH_H
#define SWIFT_DISC_PATCH_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Disc patch case
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948
 */
struct external_potential {

  /*! Surface density of the disc */
  double surface_density;

  /*! Disc scale-height */
  double scale_height;

  /*! Position of the disc along the z-axis */
  double z_disc;

  /*! Dynamical time of the system */
  double dynamical_time;

  /*! Time over which to grow the disk in units of the dynamical time */
  double growth_time;

  /*! Time-step condition pre-factor */
  double timestep_mult;
};

/**
 * @brief Computes the time-step from the acceleration due to a hydrostatic
 * disc.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948
 *
 * @param time The current time.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  /* initilize time step to disc dynamical time */
  const float dt_dyn = potential->dynamical_time;
  float dt = dt_dyn;

  /* absolute value of height above disc */
  const float dz = fabsf(g->x[2] - potential->z_disc);

  /* vertical acceleration */
  const float z_accel = 2.f * M_PI * phys_const->const_newton_G *
                        potential->surface_density *
                        tanhf(dz / potential->scale_height);

  /* demand that dt * velocity <  fraction of scale height of disc */
  float dt1 = FLT_MAX;
  if (g->v_full[2] != 0.f) {
    dt1 = potential->scale_height / fabsf(g->v_full[2]);
    if (dt1 < dt) dt = dt1;
  }

  /* demand that dt^2 * acceleration < fraction of scale height of disc */
  float dt2 = FLT_MAX;
  if (z_accel != 0.f) {
    dt2 = potential->scale_height / fabsf(z_accel);
    if (dt2 < dt * dt) dt = sqrtf(dt2);
  }

  /* demand that dt^3 * jerk < fraction of scale height of disc */
  float dt3 = FLT_MAX;
  if (g->v_full[2] != 0.f) {
    const float dz_accel_over_dt =
        2.f * M_PI * phys_const->const_newton_G * potential->surface_density /
        potential->scale_height / coshf(dz / potential->scale_height) /
        coshf(dz / potential->scale_height) * fabsf(g->v_full[2]);

    dt3 = potential->scale_height / fabsf(dz_accel_over_dt);
    if (dt3 < dt * dt * dt) dt = cbrtf(dt3);
  }

  return potential->timestep_mult * dt;
}

/**
 * @brief Computes the gravitational acceleration along z due to a hydrostatic
 * disc
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 17.
 *
 * @param time The current time in internal units.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dz = g->x[2] - potential->z_disc;
  const float t_dyn = potential->dynamical_time;

  float reduction_factor = 1.f;
  if (time < potential->growth_time * t_dyn)
    reduction_factor = time / (potential->growth_time * t_dyn);

  /* Accelerations. Note that they are multiplied by G later on */
  const float z_accel = reduction_factor * 2.f * M_PI *
                        potential->surface_density *
                        tanhf(fabsf(dz) / potential->scale_height);

  if (dz > 0) g->a_grav[2] -= z_accel;
  if (dz < 0) g->a_grav[2] += z_accel;
}

/**
 * @brief Computes the gravitational potential energy of a particle in the
 * disc patch potential.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 24.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* gp) {

  const float dz = gp->x[2] - potential->z_disc;
  const float t_dyn = potential->dynamical_time;

  float reduction_factor = 1.f;
  if (time < potential->growth_time * t_dyn)
    reduction_factor = time / (potential->growth_time * t_dyn);

  /* Accelerations. Note that they are multiplied by G later on */
  return reduction_factor * 2.f * M_PI * phys_const->const_newton_G *
         potential->surface_density * potential->scale_height *
         logf(coshf(dz / potential->scale_height));
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 22.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    const struct swift_params* parameter_file,
    const struct phys_const* phys_const, const struct UnitSystem* us,
    const struct space* s, struct external_potential* potential) {

  potential->surface_density = parser_get_param_double(
      parameter_file, "DiscPatchPotential:surface_density");
  potential->scale_height = parser_get_param_double(
      parameter_file, "DiscPatchPotential:scale_height");
  potential->z_disc =
      parser_get_param_double(parameter_file, "DiscPatchPotential:z_disc");
  potential->timestep_mult = parser_get_param_double(
      parameter_file, "DiscPatchPotential:timestep_mult");
  potential->growth_time = parser_get_opt_param_double(
      parameter_file, "DiscPatchPotential:growth_time", 0.);
  potential->dynamical_time =
      sqrt(potential->scale_height /
           (phys_const->const_newton_G * potential->surface_density));
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Disk-patch' with properties surface_density = %e "
      "disc height= %e scale height = %e timestep multiplier = %e.",
      potential->surface_density, potential->z_disc, potential->scale_height,
      potential->timestep_mult);

  if (potential->growth_time > 0.)
    message("Disc will grow for %f dynamical times.", potential->growth_time);
}

#endif /* SWIFT_DISC_PATCH_H */
