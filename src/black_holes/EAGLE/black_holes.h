/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_BLACK_HOLES_H
#define SWIFT_EAGLE_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "cosmology.h"
#include "dimension.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"

/* Standard includes */
#include <float.h>

/**
 * @brief Computes the gravity time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart* const bp) {

  return FLT_MAX;
}

/**
 * @brief Initialises the b-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_first_init_bpart(
    struct bpart* bp) {

  bp->time_bin = 0;
  bp->subgrid_mass = bp->mass;
  bp->energy_reservoir = 0.;
}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart* bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_density[i] = -1;
  bp->num_ngb_density = 0;
#endif

  bp->density.wcount = 0.f;
  bp->density.wcount_dh = 0.f;
  bp->rho_gas = 0.f;
  bp->sound_speed_gas = 0.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->ngb_mass = 0.f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param bp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void black_holes_predict_extra(
    struct bpart* restrict bp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param bp The particle.
 */
__attribute__((always_inline)) INLINE static void
black_holes_reset_predicted_values(struct bpart* restrict bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart* bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart* bp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  bp->sound_speed_gas *= h_inv_dim;
  bp->velocity_gas[0] *= h_inv_dim;
  bp->velocity_gas[1] *= h_inv_dim;
  bp->velocity_gas[2] *= h_inv_dim;

  const float rho_inv = 1.f / bp->rho_gas;

  /* Finish the calculation by undoing the mass smoothing */
  bp->sound_speed_gas *= rho_inv;
  bp->velocity_gas[0] *= rho_inv;
  bp->velocity_gas[1] *= rho_inv;
  bp->velocity_gas[2] *= rho_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_bpart_has_no_neighbours(struct bpart* restrict bp,
                                    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;
}

/**
 *
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const double dt) {

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd = props->f_Edd;
  const double epsilon_r = props->epsilon_r;
  const double epsilon_f = props->epsilon_f;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (all internal units)
   */
  const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
  const double gas_c_phys = bp->sound_speed_gas * cosmo->a_factor_sound_speed;
  const double gas_v_peculiar[3] = {bp->velocity_gas[0] * cosmo->a_inv,
                                    bp->velocity_gas[1] * cosmo->a_inv,
                                    bp->velocity_gas[2] * cosmo->a_inv};

  const double bh_v_peculiar[3] = {bp->v[0] * cosmo->a_inv,
                                   bp->v[1] * cosmo->a_inv,
                                   bp->v[2] * cosmo->a_inv};

  /* Difference in peculiar velocity between the gas and the BH */
  const double v_diff_peculiar[3] = {gas_v_peculiar[0] - bh_v_peculiar[0],
                                     gas_v_peculiar[1] - bh_v_peculiar[1],
                                     gas_v_peculiar[2] - bh_v_peculiar[2]};
  const double v_diff_norm2 = v_diff_peculiar[0] * v_diff_peculiar[0] +
                              v_diff_peculiar[1] * v_diff_peculiar[1] +
                              v_diff_peculiar[2] * v_diff_peculiar[2];

  /* We can now compute the Bondi accretion rate (internal units) */
  const double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double denominator2 = v_diff_norm2 + gas_c_phys2;
  const double denominator_inv = 1. / sqrt(denominator2);
  const double Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass *
                            gas_rho_phys * denominator_inv * denominator_inv *
                            denominator_inv;

  /* Compute the Eddington rate (internal units) */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (epsilon_r * c * sigma_Thomson);

  /* Limit the accretion rate to the Eddington fraction */
  const double accr_rate = max(Bondi_rate, f_Edd * Eddington_rate);

  /* Factor in the radiative efficiency */
  const double mass_rate = (1. - epsilon_r) * accr_rate;
  const double luminosity = epsilon_r * accr_rate * c * c;

  /* Integrate forward in time */
  bp->subgrid_mass += mass_rate * dt;
  bp->energy_reservoir += luminosity * epsilon_f * dt;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on black hole, therefore no need to use
 * it.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_reset_feedback(
    struct bpart* restrict bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

#endif /* SWIFT_EAGLE_BLACK_HOLES_H */
