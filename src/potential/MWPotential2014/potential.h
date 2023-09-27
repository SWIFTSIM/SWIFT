/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023  Darwin Roduit
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

#ifndef SWIFT_POTENTIAL_MWPotential2014_H
#define SWIFT_POTENTIAL_MWPotential2014_H

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

#ifdef HAVE_LIBGSL
#include <gsl/gsl_sf_gamma.h>
#endif

/**
 * @brief External Potential Properties - MWPotential2014 composed by
 * NFW + Miyamoto-Nagai + Power Spherical cut-off potentials
 *
 * halo --> rho_NFW(r) = rho_0 / ( (r/R_s)*(1+r/R_s)^2 )
 * disk --> Phi_MN(R,z) = -G * Mdisk / (R^2 + (Rdisk +
 * (z^2+Zdisk^2)^1/2)^2)^(1/2) bulge --> rho_PSC(r) =
 * amplitude*(r_1/r)^alpha*exp(-(r/r_c)^2)
 *
 * We however parametrise this in terms of c and virial_mass, Mdisk, Rdisk
 * and Zdisk. Also, each potential is given a contribution amplitude such that
 * the resulting potential is:
 *      Phi_tot = f_1 * Phi_NFW + f_2 * Phi_MN + f_3 * Phi_PSC,
 * with f_1, f_2 and f_3 contained in the array f.
 *
 * This potential is inspired by the following article:
 * galpy: A Python Library for Galactic Dynamics, Jo Bovy (2015),
 * Astrophys. J. Supp., 216, 29 (arXiv/1412.3451).
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x[3];

  /*! The scale radius of the NFW potential */
  double r_s;

  /*! The pre-factor \f$ 4 \pi G \rho_0 \r_s^3 \f$ */
  double pre_factor;

  /*! Hubble parameter */
  double H;

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

  /*! Amplitude for the PSC potential */
  double amplitude;

  /*! Reference radius for amplitude */
  double r_1;

  /*! Inner power */
  double alpha;

  /*! Cut-off radius */
  double r_c;

  /*! Contribution of each potential : f[0]*NFW + f[1]*MN + f[2]*PSP */
  double f[3];

  /*! Prefactor \f$ 2 \pi amplitude r_1^\alpha r_c^(3-\alpha) \f$ */
  double prefactor_psc_1;

  /*! Prefactor \f$ 2 \pi amplitude r_1^\alpha r_c^(2-\alpha) \f$ */
  double prefactor_psc_2;

  /*! Gamma function evaluation \f$ \Gamma((3-\alpha)/2 \f$ */
  double gamma_psc;

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
 * @brief Computes the time-step due to the acceleration from the NFW + MN + PSC
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

#ifdef HAVE_LIBGSL

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);

  /* Vcirc for NFW */
  const float M_NFW = potential->pre_factor * (logf(1.0f + r / potential->r_s) -
                                               r / (r + potential->r_s));
  const float Vcirc_NFW = sqrtf((phys_const->const_newton_G * M_NFW) / r);

  /* Now for MN */
  const float R = sqrtf(R2);
  const float sqrt_term = sqrtf(dz * dz + potential->Zdisk * potential->Zdisk);
  const float MN_denominator =
      powf(R2 + powf(potential->Rdisk + sqrt_term, 2.0f), 1.5f);
  const float dPhi_dR_MN = potential->Mdisk * R / MN_denominator;
  const float dPhi_dz_MN = potential->Mdisk * dz *
                           (potential->Rdisk + sqrt_term) /
                           (sqrt_term * MN_denominator);
  const float Vcirc_MN = sqrtf(phys_const->const_newton_G * R * dPhi_dR_MN +
                               phys_const->const_newton_G * dz * dPhi_dz_MN);

  /* Now for PSC */
  const float r2 = r * r;
  const float M_psc =
      potential->prefactor_psc_1 *
      (potential->gamma_psc -
       gsl_sf_gamma_inc(1.5f - 0.5f * potential->alpha,
                        r2 / (potential->r_c * potential->r_c)));
  const float Vcirc_PSC = sqrtf(phys_const->const_newton_G * M_psc / r);

  /* Total circular velocity */
  const float Vcirc = sqrtf(potential->f[0] * Vcirc_NFW * Vcirc_NFW +
                            potential->f[1] * Vcirc_MN * Vcirc_MN +
                            potential->f[2] * Vcirc_PSC * Vcirc_PSC);

  const float period = 2.0f * M_PI * r / Vcirc;

  /* Time-step as a fraction of the circular period */
  const float time_step = potential->timestep_mult * period;

  return max(time_step, potential->mintime);

#else
  error("Code not compiled with GSL. Can't compute MWPotential2014.");
  return 0.0;
#endif
}

/**
 * @brief Computes the gravitational acceleration from an NFW Halo potential +
 * MN disk + PSC bulge.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * x
 *      - dphi_PSC/dr*x/r
 * a_y = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * y
 *      - dphi_PSC/dr*y/r
 * a_z = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * z
 *      - dphi_PSC/dr*z/r
 *
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

#ifdef HAVE_LIBGSL

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  /* First for the NFW part */
  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);

  const float r_inv = 1.0f / r;
  const float M_NFW = potential->pre_factor * (logf(1.0f + r / potential->r_s) -
                                               r / (r + potential->r_s));
  const float dpot_dr_NFW = M_NFW * r_inv * r_inv;
  const float pot_nfw =
      -potential->pre_factor * logf(1.0f + r / potential->r_s) * r_inv;
  g->a_grav[0] -= potential->f[0] * dpot_dr_NFW * dx * r_inv;
  g->a_grav[1] -= potential->f[0] * dpot_dr_NFW * dy * r_inv;
  g->a_grav[2] -= potential->f[0] * dpot_dr_NFW * dz * r_inv;
  gravity_add_comoving_potential(g, potential->f[0] * pot_nfw);

  /* Now the the MN disk */
  const float f1 = sqrtf(potential->Zdisk * potential->Zdisk + dz * dz);
  const float f2 = potential->Rdisk + f1;
  const float f3 = powf(R2 + f2 * f2, -1.5f);
  const float mn_term = potential->Rdisk + sqrtf(potential->Zdisk + dz * dz);
  const float pot_mn = -potential->Mdisk / sqrtf(R2 + mn_term * mn_term);

  g->a_grav[0] -= potential->f[1] * potential->Mdisk * f3 * dx;
  g->a_grav[1] -= potential->f[1] * potential->Mdisk * f3 * dy;
  g->a_grav[2] -= potential->f[1] * potential->Mdisk * f3 * (f2 / f1) * dz;
  gravity_add_comoving_potential(g, potential->f[1] * pot_mn);

  /* Now the the PSC bulge */
  const float r2 = r * r;
  const float M_psc =
      potential->prefactor_psc_1 *
      (potential->gamma_psc -
       gsl_sf_gamma_inc(1.5f - 0.5f * potential->alpha,
                        r2 / (potential->r_c * potential->r_c)));
  const float dpot_dr = M_psc / r2;
  const float pot_psc =
      -M_psc / r - potential->prefactor_psc_2 *
                       gsl_sf_gamma_inc(1.0f - 0.5f * potential->alpha,
                                        r2 / (potential->r_c * potential->r_c));

  g->a_grav[0] -= potential->f[2] * dpot_dr * dx * r_inv;
  g->a_grav[1] -= potential->f[2] * dpot_dr * dy * r_inv;
  g->a_grav[2] -= potential->f[2] * dpot_dr * dz * r_inv;
  gravity_add_comoving_potential(g, potential->f[2] * pot_psc);

#else
  error("Code not compiled with GSL. Can't compute MWPotential2014.");
#endif
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * NFW potential + MN potential.
 *
 * phi = f[0] * (-4 * pi * G * rho_0 * r_s^3 * ln(1+r/r_s)) - f[1] * (G * Mdisk
 * / sqrt(R^2 + (Rdisk + sqrt(z^2 + Zdisk^2))^2)) + f[2] * [- G / r * (2 * pi *
 * amplitude * r_1^alpha r_c^(3 - alpha) * gamma_inf((3 - alpha)/2, r^2 / r_c^2)
 * ) - 2 * pi * G * amplitude * r_1^alpha * r_c^(2 - alpha)
 * Gamma_sup((2-alpha)/2, r^2 / r_c^2 ) ]
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

#ifdef HAVE_LIBGSL

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  /* First for the NFW profile */
  const float R2 = dx * dx + dy * dy;
  const float r = sqrtf(R2 + dz * dz + potential->eps * potential->eps);
  const float r_inv = 1.0f / r;
  const float pot_nfw =
      -potential->pre_factor * logf(1.0f + r / potential->r_s) * r_inv;

  /* Now for the MN disk */
  const float sqrt_term = sqrtf(dz * dz + potential->Zdisk * potential->Zdisk);
  const float MN_denominator =
      sqrtf(R2 + powf(potential->Rdisk + sqrt_term, 2.0f));
  const float mn_pot = -potential->Mdisk / MN_denominator;

  /* Now for PSC bulge */
  const float r2 = r * r;
  const float M_psc =
      potential->prefactor_psc_1 *
      (potential->gamma_psc -
       gsl_sf_gamma_inc(1.5f - 0.5f * potential->alpha,
                        r2 / (potential->r_c * potential->r_c)));
  const float psc_pot =
      -M_psc / r - potential->prefactor_psc_2 *
                       gsl_sf_gamma_inc(1.0f - 0.5f * potential->alpha,
                                        r2 / (potential->r_c * potential->r_c));

  return phys_const->const_newton_G *
         (potential->f[0] * pot_nfw + potential->f[1] * mn_pot +
          potential->f[2] * psc_pot);

#else
  error("Code not compiled with GSL. Can't compute MWPotential2014.");
  return 0.0;
#endif
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

#ifdef HAVE_LIBGSL

  /* Read in the position of the centre of potential */
  parser_get_param_double_array(
      parameter_file, "MWPotential2014Potential:position", 3, potential->x);

  /* Is the position absolute or relative to the centre of the box? */
  const int useabspos = parser_get_param_int(
      parameter_file, "MWPotential2014Potential:useabspos");

  /* Define the default value in the above system of units*/
  const double c_200_default = 9.823403437774843;
  const double M_200_Msun_default = 147.41031542774076e10; /* M_sun  */
  const double H_default = 127.78254614201471e-2;          /* no unit  */
  const double Mdisk_Msun_default = 6.8e10;                /* M_sun  */
  const double Rdisk_kpc_default = 3.0;                    /* kpc  */
  const double Zdisk_kpc_default = 0.280;                  /* kpc  */
  const double amplitude_default = 1.0;                    /* no unit  */
  const double r_1_kpc_default = 1.0;                      /* kpc  */
  const double alpha_default = 1.8;                        /* no unit  */
  const double r_c_kpc_default = 1.9;                      /* kpc  */
  potential->f[0] = 0.4367419745056084;                    /* no unit  */
  potential->f[1] = 1.002641971008805;                     /* no unit */
  potential->f[2] = 0.022264787598364262;                  /* no unit */

  if (!useabspos) {
    potential->x[0] += s->dim[0] / 2.;
    potential->x[1] += s->dim[1] / 2.;
    potential->x[2] += s->dim[2] / 2.;
  }

  /* Read the other parameters of the model */
  potential->timestep_mult = parser_get_param_double(
      parameter_file, "MWPotential2014Potential:timestep_mult");

  /* Bug fix : Read the softening length from the params file */
  potential->eps = parser_get_param_double(parameter_file,
                                           "MWPotential2014Potential:epsilon");

  potential->c_200 = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:concentration", c_200_default);
  potential->M_200 = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:M_200_Msun",
      M_200_Msun_default);
  potential->H = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:H", H_default);
  potential->Mdisk = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:Mdisk_kpc", Mdisk_Msun_default);
  potential->Rdisk = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:Rdisk_kpc", Rdisk_kpc_default);
  potential->Zdisk = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:Zdisk_kpc", Zdisk_kpc_default);
  potential->amplitude = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:amplitude", amplitude_default);
  potential->r_1 = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:r_1_kpc", r_1_kpc_default);
  potential->alpha = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:alpha", alpha_default);
  potential->r_c = parser_get_opt_param_double(
      parameter_file, "MWPotential2014Potential:r_c_kpc", r_c_kpc_default);
  parser_get_opt_param_double_array(
      parameter_file, "MWPotential2014Potential:potential_factors", 3,
      potential->f);

  /* Convert to internal system of units by using the
   * physical constants defined in this system */
  potential->M_200 *= phys_const->const_solar_mass;
  potential->H *= phys_const->const_reduced_hubble;
  potential->Mdisk *= phys_const->const_solar_mass;
  potential->Rdisk *= 1000. * phys_const->const_parsec;
  potential->Zdisk *= 1000. * phys_const->const_parsec;
  potential->r_1 *= 1000. * phys_const->const_parsec;
  potential->r_c *= 1000. * phys_const->const_parsec;

  /* Compute rho_c */
  const double rho_c = 3.0 * potential->H * potential->H /
                       (8.0 * M_PI * phys_const->const_newton_G);

  /* Compute R_200 */
  const double R_200 =
      cbrtf(3.0 * potential->M_200 / (4. * M_PI * 200.0 * rho_c));

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

  /* Prefactor for the mass of the PSC profile */
  potential->prefactor_psc_1 = 2.0 * M_PI * potential->amplitude *
                               pow(potential->r_1, potential->alpha) *
                               pow(potential->r_c, 3.0 - potential->alpha);

  /* Gamma function value for the mass of the PSC profile */
  potential->gamma_psc = gsl_sf_gamma(1.5 - 0.5 * potential->alpha);

  /* Prefactor for the potential of the PSC profile */
  potential->prefactor_psc_2 = 2.0 * M_PI * potential->amplitude *
                               pow(potential->r_1, potential->alpha) *
                               pow(potential->r_c, 2.0 - potential->alpha);

  /* Compute the orbital time at the softening radius */
  const double sqrtgm = sqrt(phys_const->const_newton_G * potential->M_200);
  const double epslnthing = log(1.0 + potential->eps / potential->r_s) -
                            potential->eps / (potential->eps + potential->r_s);

  potential->mintime = 2. * M_PI * potential->eps * sqrtf(potential->eps) *
                       sqrtf(potential->log_c200_term / epslnthing) / sqrtgm *
                       potential->timestep_mult;
#else
  error("Code not compiled with GSL. Can't compute MWPotential2014.");
#endif
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'MWPotential2014' "
      "with properties (in internal units) are "
      "(x,y,z) = "
      "(%e, %e, %e), c_200 = %e, M_200 = %e, H = %e, M_disk = %e, R_disk = %e, "
      "z_disk = %e, amplitude = %e, r_1 = %e, alpha = %e, r_c = %e, timestep "
      "multiplier = %e mintime = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->c_200,
      potential->M_200, potential->H, potential->Mdisk, potential->Rdisk,
      potential->Zdisk, potential->amplitude, potential->r_1, potential->alpha,
      potential->r_c, potential->timestep_mult, potential->mintime);
}

#endif /* SWIFT_POTENTIAL_MWPotential2014_H */
