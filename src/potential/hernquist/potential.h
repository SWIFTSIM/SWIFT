/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_POTENTIAL_HERNQUIST_H
#define SWIFT_POTENTIAL_HERNQUIST_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Hernquist potential
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x[3];

  /*! Mass of the halo */
  double mass;

  /*! Scale length (often as a, to prevent confusion with the cosmological
   * scale-factor we use al) */
  double al;

  /*! Square of the softening length. Acceleration tends to zero within this
   * distance from the origin */
  double epsilon2;

  /* Minimum timestep of the potential given by the timestep multiple
   * times the orbital time at the softening length */
  double mintime;

  /*! Time-step condition pre-factor, is multiplied times the circular orbital
   * time to get the time steps */
  double timestep_mult;
};

/**
 * @brief Computes the time-step in a Hernquist potential based on a
 *        fraction of the circular orbital time
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

  const float G_newton = phys_const->const_newton_G;

  /* Calculate the relative potential with respect to the centre of the
   * potential */
  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  /* calculate the radius  */
  const float r = sqrtf(dx * dx + dy * dy + dz * dz + potential->epsilon2);
  const float sqrtgm_inv = 1.f / sqrtf(G_newton * potential->mass);

  /* Calculate the circular orbital period */
  const float period = 2.f * M_PI * sqrtf(r) * potential->al *
                       (1 + r / potential->al) * sqrtgm_inv;

  /* Time-step as a fraction of the cirecular orbital time */
  const float time_step = potential->timestep_mult * period;

  return max(time_step, potential->mintime);
}

/**
 * @brief Computes the gravitational acceleration from an Hernquist potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = - GM / (a+r)^2 * x/r
 * a_y = - GM / (a+r)^2 * y/r
 * a_z = - GM / (a+r)^2 * z/r
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  /* Determine the position relative to the centre of the potential */
  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  /* Calculate the acceleration */
  const float r = sqrtf(dx * dx + dy * dy + dz * dz + potential->epsilon2);
  const float r_plus_a_inv = 1.f / (r + potential->al);
  const float r_plus_a_inv2 = r_plus_a_inv * r_plus_a_inv;
  const float term = -potential->mass * r_plus_a_inv2 / r;

  g->a_grav[0] += term * dx;
  g->a_grav[1] += term * dy;
  g->a_grav[2] += term * dz;
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * Hernquist potential.
 *
 * phi = - GM/(r+a)
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
  const float r_plus_alinv = 1.f / (r + potential->al);
  return -phys_const->const_newton_G * potential->mass * r_plus_alinv;
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

  /* Define the default value */
  static const int idealized_disk_default = 0;
  static const double M200_default = 0.;
  static const double V200_default = 0.;
  static const double R200_default = 0.;

  /* Read in the position of the centre of potential */
  parser_get_param_double_array(parameter_file, "HernquistPotential:position",
                                3, potential->x);

  /* Is the position absolute or relative to the centre of the box? */
  const int useabspos =
      parser_get_param_int(parameter_file, "HernquistPotential:useabspos");

  if (!useabspos) {
    potential->x[0] += s->dim[0] / 2.;
    potential->x[1] += s->dim[1] / 2.;
    potential->x[2] += s->dim[2] / 2.;
  }

  /* check whether we use the more advanced idealized disk setting */
  const int usedisk = parser_get_opt_param_int(
      parameter_file, "HernquistPotential:idealizeddisk",
      idealized_disk_default);

  if (!usedisk) {
    /* Read the parameters of the model in the case of the simple
     * potential form \f$ \Phi = - \frac{GM}{r+a} \f$ */
    potential->mass =
        parser_get_param_double(parameter_file, "HernquistPotential:mass");
    potential->al = parser_get_param_double(parameter_file,
                                            "HernquistPotential:scalelength");
  } else {

    /* Read the parameters in the case of a idealized disk
     * There are 3 different possible input parameters M200, V200 and R200
     * First read in the mandatory parameters in this case */

    const float G_newton = phys_const->const_newton_G;
    const float kmoversoverMpc = phys_const->const_reduced_hubble;

    /* Initialize the variables */
    double M200 = parser_get_opt_param_double(
        parameter_file, "HernquistPotential:M200", M200_default);
    double V200 = parser_get_opt_param_double(
        parameter_file, "HernquistPotential:V200", V200_default);
    double R200 = parser_get_opt_param_double(
        parameter_file, "HernquistPotential:R200", R200_default);
    const double h =
        parser_get_param_double(parameter_file, "HernquistPotential:h");

    /* Hubble constant assumed for halo masses conversion */
    const double H0 = h * kmoversoverMpc;

    /* There are 3 legit runs possible with use disk,
     * with a known M200, V200 or R200 */
    if (M200 != 0.0) {
      /* Calculate V200 and R200 from M200 */
      V200 = cbrt(10. * M200 * G_newton * H0);
      R200 = V200 / (10 * H0);

    } else if (V200 != 0.0) {

      /* Calculate M200 and R200 from V200 */
      M200 = V200 * V200 * V200 / (10. * G_newton * H0);
      R200 = V200 / (10 * H0);
    } else if (R200 != 0.0) {

      /* Calculate M200 and V200 from R200 */
      V200 = 10. * H0 * R200;
      M200 = V200 * V200 * V200 / (10. * G_newton * H0);
    } else {
      error("Please specify one of the 3 variables M200, V200 or R200");
    }

    /* message("M200 = %g, R200 = %g, V200 = %g", M200, R200, V200); */
    /* message("H0 = %g", H0); */

    /* get the concentration from the parameter file */
    const double concentration = parser_get_param_double(
        parameter_file, "HernquistPotential:concentration");

    /* Calculate the Scale radius using the NFW definition */
    const double RS = R200 / concentration;

    /* Calculate the Hernquist equivalent scale length */
    potential->al = RS * sqrt(1. * (log(1. + concentration) -
                                    concentration / (1. + concentration)));

    /* Depending on the disk mass and and the bulge mass the halo
     * gets a different mass, because of this we read the fractions
     * from the parameter file and calculate the absolute mass*/
    const double diskfraction = parser_get_param_double(
        parameter_file, "HernquistPotential:diskfraction");
    const double bulgefraction = parser_get_param_double(
        parameter_file, "HernquistPotential:bulgefraction");
    /* Calculate the mass of the bulge and disk from the parameters  */
    const double Mdisk = M200 * diskfraction;
    const double Mbulge = M200 * bulgefraction;

    /* Store the mass of the DM halo */
    potential->mass = M200 - Mdisk - Mbulge;
  }

  /* Retrieve the timestep and softening of the potential */
  potential->timestep_mult = parser_get_param_float(
      parameter_file, "HernquistPotential:timestep_mult");
  const float epsilon =
      parser_get_param_double(parameter_file, "HernquistPotential:epsilon");
  potential->epsilon2 = epsilon * epsilon;

  /* Compute the minimal time-step. */
  /* This is the circular orbital time at the softened radius */
  const float sqrtgm = sqrtf(phys_const->const_newton_G * potential->mass);
  potential->mintime = 2.f * sqrtf(epsilon) * potential->al * M_PI *
                       (1. + epsilon / potential->al) / sqrtgm *
                       potential->timestep_mult;
}

/**
 * @brief prints the properties of the external potential to stdout.
 *
 * @param  potential the external potential properties.
 */
static inline void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "external potential is 'hernquist' with properties are (x,y,z) = (%e, "
      "%e, %e), mass = %e "
      "scale length = %e , minimum time = %e "
      "timestep multiplier = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->mass,
      potential->al, potential->mintime, potential->timestep_mult);
}

#endif /* SWIFT_POTENTIAL_HERNQUIST_H */
