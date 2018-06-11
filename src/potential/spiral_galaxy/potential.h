/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_SPIRAL_GALAXY_H
#define SWIFT_SPIRAL_GALAXY_H

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
 * @brief External Potential Properties - Spiral galaxy case
 */
struct external_potential {

  /*! @brief Center of the potential. */
  double x[3];

  /*! @brief Number of spiral arms \f$N\f$. */
  int number_of_arms;

  /*! @brief Pitch angle of the spiral pattern \f$\alpha{}\f$. */
  double pitch_angle;

  /*! @brief Scale length of the arms \f$R_s\f$. */
  double scale_length;

  /*! @brief Mid-plane arm density \f$\rho{}_0\f$. */
  double mid_plane_density;

  /*! @brief Radius at which the arm density is maximal \f$r_0\f$. */
  double mid_plane_radius;

  /*! @brief Scale height of the arms \f$H\f$. */
  double scale_height;

  /*! @brief Start position angle \f$\phi{}_p(r_0)\f$. */
  double start_phi;

  /*! @brief Pattern rotation speed \f$\Omega{}\f$. */
  double pattern_speed;

  /*! @brief Central potential circular velocity squared \f$v_0^2\f$. */
  double circular_velocity2;

  /*! @brief Central potential characteristic radius squared \f$R_c^2\f$. */
  double characteristic_radius2;

  /*! @brief Central potential vertical scale factor squared \f$z_q^2\f$. */
  double vertical_scale_factor2;

  /*! @brief Time step limit for the potential. */
  double timestep_limit;
};

/**
 * @brief Computes the time-step from the acceleration due to a sine wave.
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

  return potential->timestep_limit;
}

/**
 * @brief Computes the gravitational acceleration along x given by the sine
 * wave.
 *
 * @param time The current time in internal units.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  static const double C[3] = {8. / (3. * M_PI), 0.5, 8. / (15. * M_PI)};

  const double newton_G_inv = 1. / phys_const->const_newton_G;

  /* get particle position */
  const double x[3] = {g->x[0] - potential->x[0], g->x[1] - potential->x[1],
                       g->x[2] - potential->x[2]};

  /* convert position to cylindrical coordinates */
  const double R2 = x[0] * x[0] + x[1] * x[1];
  const double Rinv = 1. / sqrt(R2);
  const double R = R2 * Rinv;
  const double phi = atan2(x[1], x[0]);
  const double z = x[2];
  const double z2 = z * z;

  /* compute central potential contribution */
  const double v_02 = potential->circular_velocity2;
  const double R_c2 = potential->characteristic_radius2;
  const double z_q2 = potential->vertical_scale_factor2;
  const double fac = v_02 / ((R_c2 + R2) * z_q2 + z2);
  const double dPhiCdr = fac * R * z_q2;
  const double dPhiCdz = fac * z;

  /* spiral arm contribution */
  const double N = potential->number_of_arms;
  const double alpha = potential->pitch_angle;
  const double R_s = potential->scale_length;
  const double rho_0 = potential->mid_plane_density;
  const double r_0 = potential->mid_plane_radius;
  const double H = potential->scale_height;
  const double phi_p = potential->start_phi;
  const double Omega = potential->pattern_speed;

  const double Gamma =
      N * (phi - phi_p + Omega * time - log(R / r_0) / tan(alpha));
  int in;
  double dPhidr = 0.;
  double dPhidphi = 0.;
  double dPhidz = 0.;
  for (in = 0; in < 3; ++in) {
    const double Cn = C[in];
    const double n = in + 1.;
    const double Kn = n * N * Rinv / sin(alpha);
    const double Kninv = 1. / Kn;
    const double dKndr = -Kn * Rinv;
    const double KnH = Kn * H;
    const double KnH2 = KnH * KnH;
    const double dKndrH = dKndr * H;
    const double betan = KnH + 0.4 * KnH2;
    const double betaninv = 1. / betan;
    const double dbetandr = (1. + 0.8 * KnH) * dKndrH;
    const double Dninv = (1. + 0.3 * KnH) / (1. + KnH + 0.3 * KnH2);
    const double dDndr =
        ((0.7 + 0.6 * KnH + 0.09 * KnH2) / (1. + 0.6 * KnH + 0.09 * KnH2)) *
        dKndrH;
    const double Knzbetan = Kn * z * betaninv;
    const double sechKnzbetan = 1. / cosh(Knzbetan);
    const double tanhKnzbetan = tanh(Knzbetan);
    const double nGamma = n * Gamma;
    const double tannGamma = tan(nGamma);

    const double Phin =
        (Cn * Dninv * Kninv) * cos(nGamma) * pow(sechKnzbetan, betan);
    const double dPhindr =
        Phin * (-1. / R_s - (Kninv + z * tanhKnzbetan) * dKndr +
                n * N * tannGamma / tan(alpha) * Rinv - dDndr * Dninv +
                (log(sechKnzbetan) + tanhKnzbetan * betaninv) * dbetandr);
    const double dPhindphi = -n * N * tannGamma * Phin;
    const double dPhindz = -Kn * tanhKnzbetan * Phin;

    dPhidr += dPhindr;
    dPhidphi += dPhindphi;
    dPhidz += dPhindz;
  }

  const double norm = -4. * M_PI * H * rho_0 * exp(-(R - r_0) / R_s);
  dPhidr *= norm;
  dPhidphi *= norm * Rinv;
  dPhidz *= norm;

  dPhidr += dPhiCdr * newton_G_inv;
  dPhidz += dPhiCdz * newton_G_inv;

  const double cosphi = x[0] * Rinv;
  const double sinphi = x[1] * Rinv;
  g->a_grav[0] = -dPhidr * cosphi + dPhidphi * sinphi;
  g->a_grav[1] = -dPhidr * sinphi - dPhidphi * cosphi;
  g->a_grav[2] = -dPhidz;
}

/**
 * @brief Computes the gravitational potential energy of a particle in the
 * sine wave.
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

  /* this potential does not really have a potential energy */
  return 0.;
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
    const struct swift_params* parameter_file,
    const struct phys_const* phys_const, const struct unit_system* us,
    const struct space* s, struct external_potential* potential) {

  /* unit conversion factors */
  const double lu = us->UnitLength_in_cgs;
  const double lu2 = lu * lu;
  const double kpc_to_length_unit = 3.086e21 / lu;
  const double degrees_to_radians = M_PI / 180.;
  const double g_cm3_to_density_unit = lu2 * lu / us->UnitMass_in_cgs;
  const double km_s_kpc_to_inverse_time_unit = 3.241e-17 * us->UnitTime_in_cgs;
  const double km_s_to_velocity_unit = 1.e5 * us->UnitTime_in_cgs / lu;

  potential->x[0] = 0.5 * s->dim[0];
  potential->x[1] = 0.5 * s->dim[1];
  potential->x[2] = 0.5 * s->dim[2];

  potential->number_of_arms = parser_get_opt_param_int(
      parameter_file, "SpiralGalaxyPotential:number_of_arms", 4);

  const double alpha = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:pitch_angle_in_degrees", 15.);
  /* convert to radians */
  potential->pitch_angle = alpha * degrees_to_radians;

  const double R_s = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:scale_length_in_kpc", 7.);
  /* convert to internal length unit */
  potential->scale_length = R_s * kpc_to_length_unit;

  const double rho_0 = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:mid_plane_density_in_g_cm3",
      2.13e-24);
  /* convert to internal density unit */
  potential->mid_plane_density = rho_0 * g_cm3_to_density_unit;

  const double r_0 = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:mid_plane_radius_in_kpc", 9.);
  /* convert to internal length unit */
  potential->mid_plane_radius = r_0 * kpc_to_length_unit;

  const double H = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:scale_height_in_kpc", 0.18);
  /* convert ot internal length unit */
  potential->scale_height = H * kpc_to_length_unit;

  const double phi_p = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:start_phi_in_degrees", 0.);
  /* convert to radians */
  potential->start_phi = phi_p * degrees_to_radians;

  const double Omega = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:pattern_speed_in_km_s_kpc", 20.);
  /* convert to inverse time unit */
  potential->pattern_speed = Omega * km_s_kpc_to_inverse_time_unit;

  double v_0 = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:circular_velocity_in_km_s", 220.);
  /* convert to internal velocity unit */
  v_0 *= km_s_to_velocity_unit;
  potential->circular_velocity2 = v_0 * v_0;

  double R_c = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:characteristic_radius_in_kpc",
      0.1);
  /* convert to internal length unit */
  R_c *= kpc_to_length_unit;
  potential->characteristic_radius2 = R_c * R_c;

  const double z_q = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:vertical_scale_factor", 0.7);
  potential->vertical_scale_factor2 = z_q * z_q;

  potential->timestep_limit = parser_get_opt_param_double(
      parameter_file, "SpiralGalaxyPotential:timestep_limit", 1.e-5);
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message("External potential is a spiral galaxy with following properties:");
  message("Number of arms (N): %i", potential->number_of_arms);
  message("Pitch angle (alpha, radians): %g", potential->pitch_angle);
  message("Arm scale length (R_s, U_L): %g", potential->scale_length);
  message("Mid-plane arm density (rho_0, U_M U_L^-3): %g",
          potential->mid_plane_density);
  message("Mid-plane peak radius (r_0, U_L): %g", potential->mid_plane_radius);
  message("Arm scale height (H, U_L): %g", potential->scale_height);
  message("Start position angle (phi_p, radians): %g", potential->start_phi);
  message("Pattern rotation speed (Omega, U_t^-1): %g",
          potential->pattern_speed);
  message("Central circular velocity (v_0, U_L U_t^-1): %g",
          sqrt(potential->circular_velocity2));
  message("Central characteristic radius (R_c, U_L): %g",
          sqrt(potential->characteristic_radius2));
  message("Central scale factor (z_q, dimensionless): %g",
          sqrt(potential->vertical_scale_factor2));
  message("Time step limit (U_t): %g", potential->timestep_limit);
}

#endif /* SWIFT_SINE_WAVE_H */
