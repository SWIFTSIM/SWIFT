/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_SINK_GETTERS_H
#define SWIFT_GEAR_SINK_GETTERS_H

#include "cosmology.h"
#include "sink_part.h"

/**
 * @file src/sink/GEAR/sink_getters.h
 * @brief Getter functions for GEAR sink scheme to avoid exposing
 * implementation details to the outer world. Keep the code clean and lean.
 */

/**
 * @brief Get the sink age in interal units.
 *
 * @param sink The #sink.
 * @param with_cosmology If we run with cosmology.
 * @param cosmo The co Birth scale-factor of the star.
 * @param with_cosmology If we run with cosmology.
 */

__attribute__((always_inline)) INLINE double sink_get_sink_age(
    const struct sink* restrict sink, const int with_cosmology,
    const struct cosmology* cosmo, const double time) {
  double sink_age;
  if (with_cosmology) {

    /* Deal with rounding issues */
    if (sink->birth_scale_factor >= cosmo->a) {
      sink_age = 0.;
    } else {
      sink_age = cosmology_get_delta_time_from_scale_factors(
          cosmo, sink->birth_scale_factor, cosmo->a);
    }
  } else {
    sink_age = time - sink->birth_time;
  }
  return sink_age;
}

/**
 * @brief Compute the rotational energy of the neighbouring gas particles.
 *
 * Note: This function must be used after having computed the rotational energy
 * per components, i.e. after sink_prepare_part_sink_formation().
 *
 * @param p The gas particle.
 *
 */
INLINE static double sink_compute_neighbour_rotation_energy_magnitude(
    const struct part* restrict p) {
  double E_rot_x = p->sink_data.E_rot_neighbours[0];
  double E_rot_y = p->sink_data.E_rot_neighbours[1];
  double E_rot_z = p->sink_data.E_rot_neighbours[2];
  double E_rot =
      sqrtf(E_rot_x * E_rot_x + E_rot_y * E_rot_y + E_rot_z * E_rot_z);
  return E_rot;
}

/**
 * @brief Retrieve the physical velocity divergence from the gas particle.
 *
 * @param p The gas particles.
 *
 */
INLINE static float sink_get_physical_div_v_from_part(
    const struct part* restrict p) {

  float div_v = 0.0;

  /* The implementation of div_v depends on the Hydro scheme. Furthermore, some
     add a Hubble flow term, some do not. We need to take care of this */
#ifdef SPHENIX_SPH
  /* SPHENIX is already including the Hubble flow. */
  div_v = hydro_get_div_v(p);
#elif GADGET2_SPH
  div_v = p->density.div_v;

  /* Add the missing term */
  div_v += hydro_dimension * cosmo->H;
#elif MINIMAL_SPH
  div_v = hydro_get_div_v(p);

  /* Add the missing term */
  div_v += hydro_dimension * cosmo->H;
#elif GASOLINE_SPH
  /* Copy the velocity divergence */
  div_v = (1. / 3.) * (p->viscosity.velocity_gradient[0][0] +
                       p->viscosity.velocity_gradient[1][1] +
                       p->viscosity.velocity_gradient[2][2]);
#elif HOPKINS_PU_SPH
  div_v = p->density.div_v;
#else
#error \
    "This scheme is not implemented. Note that Different scheme apply the Hubble flow in different places. Be careful about it."
#endif
  return div_v;
}

/**
 * @brief Compute the angular momentum-based criterion for sink-sink
 * interaction.
 *
 * This function calculates the angular momentum of a sink particle relative to
 * another particle (sink or gas) and evaluates the Keplerian angular momentum.
 *
 * @param dx Comoving vector separating the two particles (pi - pj).
 * @param dv_plus_H_flow Comoving relative velocity including the Hubble flow.
 * @param r Comoving distance between the two particles.
 * @param r_cut_i Comoving cut-off radius of particle i.
 * @param mass_i Mass of particle i.
 * @param cosmo The cosmological parameters and properties
 * @param grav_props The gravity scheme parameters and properties
 * @param L2_kepler (return) Keplerian angular momentum squared of particle i.
 * @param L2_j (return) Specific angular momentum squared relative to particle
 * j.
 */
__attribute__((always_inline)) INLINE static void
sink_compute_angular_momenta_criterion(
    const float dx[3], const float dv_plus_H_flow[3], const float r,
    const float r_cut_i, const float mass_i, const struct cosmology* cosmo,
    const struct gravity_props* grav_props, float* L2_kepler, float* L2_j) {

  /* Compute the physical relative velocity between the particles */
  const float dv_physical[3] = {dv_plus_H_flow[0] * cosmo->a_inv,
                                dv_plus_H_flow[1] * cosmo->a_inv,
                                dv_plus_H_flow[2] * cosmo->a_inv};

  /* Compute the physical distance between the particles */
  const float dx_physical[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a,
                                dx[2] * cosmo->a};
  const float r_physical = r * cosmo->a;

  /* Momentum check------------------------------------------------------- */
  /* Relative momentum of the gas */
  const float specific_angular_momentum[3] = {
      dx_physical[1] * dv_physical[2] - dx_physical[2] * dv_physical[1],
      dx_physical[2] * dv_physical[0] - dx_physical[0] * dv_physical[2],
      dx_physical[0] * dv_physical[1] - dx_physical[1] * dv_physical[0]};

  *L2_j = specific_angular_momentum[0] * specific_angular_momentum[0] +
          specific_angular_momentum[1] * specific_angular_momentum[1] +
          specific_angular_momentum[2] * specific_angular_momentum[2];

  /* Keplerian angular speed squared */
  const float omega_acc_2 =
      grav_props->G_Newton * mass_i / (r_physical * r_physical * r_physical);

  /*Keplerian angular momentum squared */
  *L2_kepler = (r_cut_i * r_cut_i * r_cut_i * r_cut_i) * omega_acc_2;
}

#endif /* SWIFT_GEAR_SINK_GETTERS_H */
