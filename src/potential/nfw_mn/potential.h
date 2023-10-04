/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019  Alejandro Benitez-Llambay
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

#ifndef SWIFT_POTENTIAL_NFW_MN_H
#define SWIFT_POTENTIAL_NFW_MN_H

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

/**
 * @brief External Potential Properties - NFW Potential + Miyamoto-Nagai
 *
 * halo --> rho(r) = rho_0 / ( (r/R_s)*(1+r/R_s)^2 )
 * disk --> phi(R,z) = -G * Mdisk / (R^2 + (Rdisk +  (z^2+Zdisk^2)^1/2)^2)^(1/2)
 *
 * We however parameterise this in terms of c and virial_mass, Mdisk, Rdisk
 * and Zdisk
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x[3];

  /*! The scale radius of the NFW potential */
  double r_s;

  /*! The pre-factor \f$ 4 \pi G \rho_0 \r_s^3 \f$ */
  double pre_factor;

  /*! The critical density of the universe */
  double rho_c;

  /*! The concentration parameter */
  double c_200;

  /*! The virial mass */
  double M_200;

  /*! Disk Size */
  double Rdisk;

  /*! Disk height */
  double Zdisk;

  /*! Disk Mass */
  double Mdisk;

  /*! Time-step condition pre_factor, this factor is used to multiply times the
   * orbital time, so in the case of 0.01 we take 1% of the orbital time as
   * the time integration steps */
  double timestep_mult;

  /*! Minimum time step based on the orbital time at the softening times
   * the timestep_mult */
  double mintime;

  /*! Common log term \f$ \ln(1+c_{200}) - \frac{c_{200}}{1 + c_{200}} \f$ */
  double log_c200_term;

  /*! Softening length */
  double eps;
};

/**
 * @brief Computes the time-step due to the acceleration from the NFW + MN
 * potential as a fraction (timestep_mult) of the circular orbital time of that
 * particle.
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

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);

  const float mr = potential->M_200 *
                   (logf(1.f + r / potential->r_s) - r / (r + potential->r_s)) /
                   potential->log_c200_term;

  const float Vcirc_NFW = sqrtf((phys_const->const_newton_G * mr) / r);
  const float f1 =
      potential->Rdisk + sqrtf(potential->Zdisk * potential->Zdisk + dz * dz);
  const float Vcirc_MN = sqrtf(phys_const->const_newton_G * potential->Mdisk *
                               R2 / pow(R2 + f1 * f1, 3.0 / 2.0));
  const float Vcirc = sqrtf(Vcirc_NFW * Vcirc_NFW + Vcirc_MN * Vcirc_MN);

  const float period = 2 * M_PI * r / Vcirc;

  /* Time-step as a fraction of the circular period */
  const float time_step = potential->timestep_mult * period;

  return max(time_step, potential->mintime);
}

/**
 * @brief Computes the gravitational acceleration from an NFW Halo potential +
 * MN disk.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * x
 * a_y = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * y
 * a_z = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * z
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  /* First for the NFW part */
  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);
  const float r_inv = 1.f / r;
  const float term1 = potential->pre_factor;
  const float term2 =
      r / (r + potential->r_s) - logf(1.0f + r / potential->r_s);

  const float acc_nfw = term1 * term2 * r_inv * r_inv * r_inv;
  const float pot_nfw =
      -potential->pre_factor * logf(1.0f + r / potential->r_s) * r_inv;

  g->a_grav[0] += acc_nfw * dx;
  g->a_grav[1] += acc_nfw * dy;
  g->a_grav[2] += acc_nfw * dz;
  gravity_add_comoving_potential(g, pot_nfw);

  /* Now the the MN disk */
  const float f1 = sqrtf(potential->Zdisk * potential->Zdisk + dz * dz);
  const float f2 = potential->Rdisk + f1;
  const float f3 = powf(R2 + f2 * f2, -1.5f);
  const float mn_term = potential->Rdisk + sqrtf(potential->Zdisk + dz * dz);
  const float pot_mn = -potential->Mdisk / sqrtf(R2 + mn_term * mn_term);

  g->a_grav[0] -= potential->Mdisk * f3 * dx;
  g->a_grav[1] -= potential->Mdisk * f3 * dy;
  g->a_grav[2] -= potential->Mdisk * f3 * (f2 / f1) * dz;
  gravity_add_comoving_potential(g, pot_mn);
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * NFW potential + MN potential.
 *
 * phi = -4 * pi * G * rho_0 * r_s^3 * ln(1+r/r_s) - G * Mdisk / sqrt(R^2 +
 * (Rdisk + sqrt(z^2 + Zdisk^2))^2)
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

  /* First for the NFW profile */
  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);
  const float term1 = -potential->pre_factor / r;
  const float term2 = logf(1.0f + r / potential->r_s);

  /* Now for the MN disk */
  const float mn_term = potential->Rdisk + sqrtf(potential->Zdisk + dz * dz);
  const float mn_pot = -potential->Mdisk / sqrtf(R2 + mn_term * mn_term);

  return phys_const->const_newton_G * (term1 * term2 + mn_pot);
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
  parser_get_param_double_array(parameter_file, "NFW_MNPotential:position", 3,
                                potential->x);

  /* Is the position absolute or relative to the centre of the box? */
  const int useabspos =
      parser_get_param_int(parameter_file, "NFW_MNPotential:useabspos");

  if (!useabspos) {
    potential->x[0] += s->dim[0] / 2.;
    potential->x[1] += s->dim[1] / 2.;
    potential->x[2] += s->dim[2] / 2.;
  }

  /* Read the other parameters of the model */
  potential->timestep_mult =
      parser_get_param_double(parameter_file, "NFW_MNPotential:timestep_mult");
  potential->c_200 =
      parser_get_param_double(parameter_file, "NFW_MNPotential:concentration");
  potential->M_200 =
      parser_get_param_double(parameter_file, "NFW_MNPotential:M_200");
  potential->rho_c = parser_get_param_double(
      parameter_file, "NFW_MNPotential:critical_density");
  potential->Mdisk =
      parser_get_param_double(parameter_file, "NFW_MNPotential:Mdisk");
  potential->Rdisk =
      parser_get_param_double(parameter_file, "NFW_MNPotential:Rdisk");
  potential->Zdisk =
      parser_get_param_double(parameter_file, "NFW_MNPotential:Zdisk");

  potential->eps = 0.05;

  /* Compute R_200 */
  const double R_200 =
      cbrtf(3.0 * potential->M_200 / (4. * M_PI * 200.0 * potential->rho_c));

  /* NFW scale-radius */
  potential->r_s = R_200 / potential->c_200;
  const double r_s3 = potential->r_s * potential->r_s * potential->r_s;

  /* Log(c_200) term appearing in many expressions */
  potential->log_c200_term =
      log(1. + potential->c_200) - potential->c_200 / (1. + potential->c_200);

  const double rho_0 =
      potential->M_200 / (4.f * M_PI * r_s3 * potential->log_c200_term);

  /* Pre-factor for the accelerations (note G is multiplied in later on) */
  potential->pre_factor = 4.0f * M_PI * rho_0 * r_s3;

  /* Compute the orbital time at the softening radius */
  const double sqrtgm = sqrt(phys_const->const_newton_G * potential->M_200);
  const double epslnthing = log(1.f + potential->eps / potential->r_s) -
                            potential->eps / (potential->eps + potential->r_s);

  potential->mintime = 2. * M_PI * potential->eps * sqrtf(potential->eps) *
                       sqrtf(potential->log_c200_term / epslnthing) / sqrtgm *
                       potential->timestep_mult;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'NFW + MN disk' with properties are (x,y,z) = "
      "(%e, %e, %e), scale radius = %e timestep multiplier = %e, mintime = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->r_s,
      potential->timestep_mult, potential->mintime);
}

#endif /* SWIFT_POTENTIAL_NFW_MN_H */
