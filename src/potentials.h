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

#ifndef SWIFT_POTENTIALS_H
#define SWIFT_POTENTIALS_H

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
#include "units.h"

/* External Potential Properties */
struct external_potential {

#ifdef EXTERNAL_POTENTIAL_POINTMASS
  struct {
    double x, y, z;
    double mass;
    float timestep_mult;
  } point_mass;
#endif

#ifdef EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL
  struct {
    double x, y, z;
    double vrot;
    float timestep_mult;
  } isothermal_potential;
#endif

#ifdef EXTERNAL_POTENTIAL_DISK_PATCH
  struct {
    double surface_density;
    double scale_height;
    double z_disk;
    double dynamical_time;
    double timestep_mult;
  } disk_patch_potential;
#endif
};

/* ------------------------------------------------------------------------- */

/* external potential: disk_patch */
#ifdef EXTERNAL_POTENTIAL_DISK_PATCH

/**
 * @brief Computes the time-step due to the acceleration from a hydrostatic disk
 * (Creasy, Theuns & Bower, 2013)
 *
 * @param phys_cont The physical constants in internal units.
 * @param gp Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_disk_patch_timestep(const struct external_potential* potential,
                                     const struct phys_const* const phys_const,
                                     const struct gpart* const g) {

  /* initilize time step to disk dynamical time */
  const float dt_dyn = potential->disk_patch_potential.dynamical_time;
  float dt = dt_dyn;

  /* absolute value of height above disk */
  const float dz = fabs(g->x[2] - potential->disk_patch_potential.z_disk);

  /* vertical cceleration */
  const float z_accel = 2 * M_PI * phys_const->const_newton_G *
                        potential->disk_patch_potential.surface_density *
                        tanh(dz / potential->disk_patch_potential.scale_height);

  /* demand that dt * velocity <  fraction of scale height of disk */
  float dt1 = FLT_MAX;
  if (fabs(g->v_full[2]) > 0) {
    dt1 = potential->disk_patch_potential.scale_height / fabs(g->v_full[2]);
    if (dt1 < dt) dt = dt1;
  }

  /* demand that dt^2 * acceleration < fraction of scale height of disk */
  float dt2 = FLT_MAX;
  if (fabs(z_accel) > 0) {
    dt2 = potential->disk_patch_potential.scale_height / fabs(z_accel);
    if (dt2 < dt * dt) dt = sqrt(dt2);
  }

  /* demand that dt^3 jerk < fraction of scale height of disk */
  float dt3 = FLT_MAX;
  if (abs(g->v_full[2]) > 0) {
    const float dz_accel_over_dt =
        2 * M_PI * phys_const->const_newton_G *
        potential->disk_patch_potential.surface_density /
        potential->disk_patch_potential.scale_height /
        cosh(dz / potential->disk_patch_potential.scale_height) /
        cosh(dz / potential->disk_patch_potential.scale_height) *
        fabs(g->v_full[2]);

    dt3 = potential->disk_patch_potential.scale_height / fabs(dz_accel_over_dt);
    if (dt3 < dt * dt * dt) dt = pow(dt3, 1. / 3.);
  }

  return potential->disk_patch_potential.timestep_mult * dt;
}

/**
 * @brief Computes the gravitational acceleration from a hydrostatic disk
 * (Creasy, Theuns & Bower, 2013)
 *
 * @param phys_cont The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void
external_gravity_disk_patch_potential(
    const double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float G_newton = phys_const->const_newton_G;
  const float dz = g->x[2] - potential->disk_patch_potential.z_disk;
  g->a_grav[0] += 0;
  g->a_grav[1] += 0;

#ifdef EXTERNAL_POTENTIAL_DISK_PATCH_ICS
  float reduction_factor = 1.;
  if (time < 5 * t_dyn) reduction_factor = time / (5 * t_dyn);
#else
  const float reduction_factor = 1.;
#endif

  const float z_accel =
      reduction_factor * 2 * G_newton * M_PI *
      potential->disk_patch_potential.surface_density *
      tanh(fabs(dz) / potential->disk_patch_potential.scale_height);

  if (dz > 0)
    g->a_grav[2] -=
        z_accel /
        G_newton; /* returned acceleraton is multiplied by G later on */
  if (dz < 0) g->a_grav[2] += z_accel / G_newton;

#ifdef EXTERNAL_POTENTIAL_DISK_PATCH_ICS
  /* TT: add viscous drag */
  g->a_grav[0] -= g->v_full[0] / (2 * t_dyn) / G_newton;
  g->a_grav[1] -= g->v_full[1] / (2 * t_dyn) / G_newton;
  g->a_grav[2] -= g->v_full[2] / (2 * t_dyn) / G_newton;
#endif
}
#endif /* EXTERNAL_POTENTIAL_DISK_PATCH */

/* ------------------------------------------------------------------------- */

#ifdef EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL

/**
 * @brief Computes the time-step due to the acceleration from an isothermal
 * potential.
 *
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_isothermalpotential_timestep(
    const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* const g) {

  const float dx = g->x[0] - potential->isothermal_potential.x;
  const float dy = g->x[1] - potential->isothermal_potential.y;
  const float dz = g->x[2] - potential->isothermal_potential.z;

  const float rinv2 = 1.f / (dx * dx + dy * dy + dz * dz);
  const float drdv =
      dx * (g->v_full[0]) + dy * (g->v_full[1]) + dz * (g->v_full[2]);
  const double vrot = potential->isothermal_potential.vrot;

  const float dota_x =
      vrot * vrot * rinv2 * (g->v_full[0] - 2 * drdv * dx * rinv2);
  const float dota_y =
      vrot * vrot * rinv2 * (g->v_full[1] - 2 * drdv * dy * rinv2);
  const float dota_z =
      vrot * vrot * rinv2 * (g->v_full[2] - 2 * drdv * dz * rinv2);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return potential->isothermal_potential.timestep_mult * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration from an isothermal potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void
external_gravity_isothermalpotential(const struct external_potential* potential,
                                     const struct phys_const* const phys_const,
                                     struct gpart* g) {

  const float G_newton = phys_const->const_newton_G;
  const float dx = g->x[0] - potential->isothermal_potential.x;
  const float dy = g->x[1] - potential->isothermal_potential.y;
  const float dz = g->x[2] - potential->isothermal_potential.z;
  const float rinv2 = 1.f / (dx * dx + dy * dy + dz * dz);

  const double vrot = potential->isothermal_potential.vrot;
  const double term = -vrot * vrot * rinv2 / G_newton;

  g->a_grav[0] += term * dx;
  g->a_grav[1] += term * dy;
  g->a_grav[2] += term * dz;
}

#endif /* EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL */

/* ------------------------------------------------------------------------- */

/* Include exteral pointmass potential */
#ifdef EXTERNAL_POTENTIAL_POINTMASS

/**
 * @brief Computes the time-step due to the acceleration from a point mass
 *
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_pointmass_timestep(const struct external_potential* potential,
                                    const struct phys_const* const phys_const,
                                    const struct gpart* const g) {

  const float G_newton = phys_const->const_newton_G;
  const float dx = g->x[0] - potential->point_mass.x;
  const float dy = g->x[1] - potential->point_mass.y;
  const float dz = g->x[2] - potential->point_mass.z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float drdv = (g->x[0] - potential->point_mass.x) * (g->v_full[0]) +
                     (g->x[1] - potential->point_mass.y) * (g->v_full[1]) +
                     (g->x[2] - potential->point_mass.z) * (g->v_full[2]);
  const float dota_x = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[0] + 3.f * rinv * rinv * drdv * dx);
  const float dota_y = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[1] + 3.f * rinv * rinv * drdv * dy);
  const float dota_z = G_newton * potential->point_mass.mass * rinv * rinv *
                       rinv * (-g->v_full[2] + 3.f * rinv * rinv * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return potential->point_mass.timestep_mult * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration of a particle due to a
 * point mass
 *
 * Note that the accelerations are multiplied by Newton's G constant later
 * on.
 *
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_pointmass(
    const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float dx = g->x[0] - potential->point_mass.x;
  const float dy = g->x[1] - potential->point_mass.y;
  const float dz = g->x[2] - potential->point_mass.z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float rinv3 = rinv * rinv * rinv;

  g->a_grav[0] += -potential->point_mass.mass * dx * rinv3;
  g->a_grav[1] += -potential->point_mass.mass * dy * rinv3;
  g->a_grav[2] += -potential->point_mass.mass * dz * rinv3;
}

#endif /* EXTERNAL_POTENTIAL_POINTMASS */

/* ------------------------------------------------------------------------- */

/* Now, some generic functions, defined in the source file */
void potential_init(const struct swift_params* parameter_file,
                    const struct phys_const* phys_const,
                    const struct UnitSystem* us,
                    struct external_potential* potential);

void potential_print(const struct external_potential* potential);

#endif /* SWIFT_POTENTIALS_H */
