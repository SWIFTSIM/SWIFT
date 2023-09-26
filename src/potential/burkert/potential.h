//
// Created by yuyttenh on 12/08/22.
//

#ifndef SWIFT_POTENTIAL_BURKERT_H
#define SWIFT_POTENTIAL_BURKERT_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "gravity.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

struct external_potential {
  /*! Position of the centre of potential */
  double x[3];

  /*! The core radius of the halo */
  double r0;

  /*! Precomputed powers of the core radius */
  double r0_2;
  double r0_3;
  double one_over_r0;

  /*! The total mass inside the core radius of the halo */
  double core_mass;

  /*! The central density of the halo */
  double central_density;

  /*! The central potential */
  double central_potential;

  /*! Time-step condition pre_factor, this factor is used to multiply times the
   * orbital time, so in the case of 0.01 we take 1% of the orbital time as
   * the time integration steps */
  double timestep_mult;
};

/**
 * @brief Computes the enclosed mass due to the NFW potential.
 *
 * See Mori and Burkert (2000) eq. 3.
 *
 * @param potential The #external_potential used in the run.
 * @param radius The radius of the particle
 */
__attribute__((always_inline)) INLINE static double enclosed_mass_burkert(
    const struct external_potential* restrict potential,
    const double r_over_r0) {

  const double pre_factor = M_PI * potential->central_density * potential->r0_3;

  return pre_factor * (2. * log(1. + r_over_r0) +
                       log(1. + r_over_r0 * r_over_r0) - 2. * atan(r_over_r0));
}

/**
 * @brief Computes the Burkert potential for the given radius.
 *
 * See Mori and Burkert (2000) eq. 2. Note that we do not multiply with Newton's
 * G yet.
 *
 * @param potential The #external_potential used in the run.
 * @param r_over_r0 The rescaled distance to the potential.
 */
__attribute__((always_inline)) INLINE static double potential_burkert(
    const struct external_potential* restrict potential,
    const double r_over_r0) {

  const double r0_over_r = 1. / r_over_r0;
  const double pre_factor = M_PI * potential->central_density * potential->r0_2;

  return potential->central_potential +
         pre_factor * (2. * (1. + r0_over_r) * atan(r_over_r0) -
                       2. * (1. + r0_over_r) * log(1. + r_over_r0) +
                       (1. - r0_over_r) * log(1. + r_over_r0 * r_over_r0));
}

/**
 * @brief Computes the gravitational acceleration from an Burkert potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = M_encl(r) /r^3 * x
 * a_y = M_encl(r) /r^3 * y
 * a_z = M_encl(r) /r^3 * z
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  /* Determine the position relative to the centre of the potential */
  const double dx = g->x[0] - potential->x[0];
  const double dy = g->x[1] - potential->x[1];
  const double dz = g->x[2] - potential->x[2];

  /* Calculate the acceleration */
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double r = sqrt(r2);
  const double r_inv = r > 0. ? 1. / r : 0.f;
  const double r_over_r0 = r * potential->one_over_r0;
  const double M_encl = enclosed_mass_burkert(potential, r_over_r0);

  const double acc = -M_encl * r_inv * r_inv * r_inv;
  const float pot = (float)potential_burkert(potential, r_over_r0);

  g->a_grav[0] += (float)(acc * dx);
  g->a_grav[1] += (float)(acc * dy);
  g->a_grav[2] += (float)(acc * dz);
  gravity_add_comoving_potential(g, pot);
}

/** @brief Computes the time-step due to the acceleration from the NFW potential
 *        as a fraction (timestep_mult) of the circular orbital time of that
 *        particle.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  const double dx = g->x[0] - potential->x[0];
  const double dy = g->x[1] - potential->x[1];
  const double dz = g->x[2] - potential->x[2];
  const double r = sqrt(dx * dx + dy * dy + dz * dz);

  const double mr =
      enclosed_mass_burkert(potential, r * potential->one_over_r0);

  const double period =
      2. * M_PI * r * sqrt(r / (phys_const->const_newton_G * mr));

  /* Time-step as a fraction of the circular period */
  return (float)(potential->timestep_mult * period);
}

/**
 * @brief Computes the gravitational potential energy of a particle in a
 * Burkert potential.
 *
 * @param time The current time (unused here).
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  const double dx = g->x[0] - potential->x[0];
  const double dy = g->x[1] - potential->x[1];
  const double dz = g->x[2] - potential->x[2];
  const double r = sqrt(dx * dx + dy * dy + dz * dz);

  return (float)potential_burkert(potential, r * potential->one_over_r0);
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct space* s,
    struct external_potential* potential) {

  /* Read in the position of the centre of potential */
  parser_get_param_double_array(parameter_file, "BurkertPotential:position", 3,
                                potential->x);

  /* Is the position absolute or relative to the centre of the box? */
  const int useabspos =
      parser_get_param_int(parameter_file, "BurkertPotential:useabspos");

  if (!useabspos) {
    potential->x[0] += s->dim[0] / 2.;
    potential->x[1] += s->dim[1] / 2.;
    potential->x[2] += s->dim[2] / 2.;
  }

  /* Read the other parameters of the model */
  potential->timestep_mult =
      parser_get_param_double(parameter_file, "BurkertPotential:timestep_mult");
  potential->core_mass =
      parser_get_param_double(parameter_file, "BurkertPotential:M0");
  const double rescaled_mass =
      potential->core_mass / (1e9 * phys_const->const_solar_mass);

  /* Compute the core radius and its powers (Mori and Burkert (2000) eq. 4) */
  const double r0 = 3.09 * pow(rescaled_mass, 3. / 7.);  // kpc
  potential->r0 = r0;
  potential->r0_2 = r0 * r0;
  potential->r0_3 = r0 * r0 * r0;
  potential->one_over_r0 = 1. / r0;

  /* Compute the central density of the halo (Mori and Burkert (2000) eq. 5).
   * This is in g / cm^3 however, so we need a conversion factor. */
  const double units = us->UnitLength_in_cgs * us->UnitLength_in_cgs *
                       us->UnitLength_in_cgs / us->UnitMass_in_cgs;
  potential->central_density = 1.46e-24 * pow(rescaled_mass, -2. / 7.) * units;

  /* Compute the central potential (We do not multiply with G yet). */
  potential->central_potential =
      -M_PI_2 * potential->central_density * potential->r0_2;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'burkert' with properties: (x,y,z) = (%e, "
      "%e, %e), core radius = %e, core mass = %e, central density = %e, "
      "timestep multiplier = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->r0,
      potential->core_mass, potential->central_density,
      potential->timestep_mult);
}

#endif  // SWIFT_POTENTIAL_BURKERT_H
