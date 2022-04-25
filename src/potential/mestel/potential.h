#ifndef SWIFT_POTENTIAL_MESTEL_H
#define SWIFT_POTENTIAL_MESTEL_H

/* Config parameters. */
#include "../config.h"

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

/**
 * @brief External Potential Properties - Mestel Disk
 * [TODO]
 */
struct external_potential {
  /*! Position of the Mestel potential */
  double x[3];

  /*! Circular Velocity */
  double v0;

  /*! Scale Radius */
  double r0;

  /*! Time-step condition pre-factor */
  float timestep_mult;

  /*! Minimum time step based on the circular orbital time at the
  *  scale raidus times the timestep_mult */
  double mintime;
};

/**
 * @brief Computes the time-step due to the acceleration from a Mestel potential
 *
 * Given by a specified fraction of circular orbital time at radius of particle
 *
 * @param time The current time.
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  const float r = sqrtf(dx * dx + dy * dy + dz * dz);
  const float time_step = potential->timestep_mult * r/potential->v0;

  return max(time_step, potential->mintime);
}

/**
 * @brief Computes the gravitational acceleration of a particle due to a
 * Mestel potential
 *
 * Note that the accelerations are multiplied by Newton's G constant later
 * on.
 *
 * We pass in the time for simulations where the potential evolves with time.
 *
 * @param time The current time.
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  const float v02 = potential->v0*potential->v0;
  const float r = sqrtf(dx * dx + dy * dy + dz * dz);
  const float rinv2G = 1.f / ((dx * dx + dy * dy + dz * dz)*phys_const->const_newton_G);
  const float pot = v02 * logf(r/potential->r0) / phys_const->const_newton_G;

  g->a_grav[0] += -v02 * dx * rinv2G;
  g->a_grav[1] += -v02 * dy * rinv2G;
  g->a_grav[2] += -v02 * dz * rinv2G;
  gravity_add_comoving_potential(g, pot);
}

/**
 * @brief Computes the gravitational potential energy of a particle in a point
 * Mestel potential.
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

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  const float r = sqrtf(dx * dx + dy * dy + dz * dz);

  return potential->v0*potential->v0 * logf(r/potential->r0);
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space we run in.
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct space* s,
    struct external_potential* potential) {

  /* Read in the position of the centre of potential */
  parser_get_param_double_array(parameter_file, "MestelPotential:position",
                                3, potential->x);

  /* Read the other parameters of the model */
  potential->v0 =
      parser_get_param_double(parameter_file, "MestelPotential:v0");
  potential->timestep_mult = parser_get_opt_param_float(
      parameter_file, "MestelPotential:timestep_mult", FLT_MAX);
  potential->r0 =
      parser_get_param_float(parameter_file, "MestelPotential:r0");
  potential->mintime = potential->timestep_mult * potential->r0 / potential->v0;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Mestel Potential' with properties (x,y,z) = (%e, %e, "
      "%e), v0 = %e, r0 = %e timestep multiplier = %e.",
      potential->x[0], potential->x[1], potential->x[2], potential->v0,
      potential->r0, potential->timestep_mult);
}

#endif /* SWIFT_POTENTIAL_MESTEL_H */
