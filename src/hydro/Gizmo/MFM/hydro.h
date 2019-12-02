/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016, 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_MFM_HYDRO_H
#define SWIFT_GIZMO_MFM_HYDRO_H

#include "../hydro_getters.h"
#include "../hydro_gradients.h"
#include "../hydro_setters.h"
#include "../hydro_unphysical.h"
#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "minmax.h"
#include "riemann.h"

#include <float.h>

/**
 * @brief Extra operations to be done during the drift
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float dt_drift, float dt_therm,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props) {

  const float h_inv = 1.0f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
  float h_corr;
  if (fabsf(w1) < 0.2f)
    h_corr = approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    h_corr = expf(w1);

  /* Limit the smoothing length correction (and make sure it is always
     positive). */
  if (h_corr < 2.0f && h_corr > 0.0f) {
    p->h *= h_corr;
  }

  /* drift the primitive variables based on the old fluxes */
  if (p->conserved.mass > 0.0f) {
    const float m_inv = 1.0f / p->conserved.mass;

    p->v[0] += p->flux.momentum[0] * dt_drift * m_inv;
    p->v[1] += p->flux.momentum[1] * dt_drift * m_inv;
    p->v[2] += p->flux.momentum[2] * dt_drift * m_inv;

#if !defined(EOS_ISOTHERMAL_GAS)
#ifdef GIZMO_TOTAL_ENERGY
    const float Etot = p->conserved.energy + p->flux.energy * dt_drift;
    const float v2 =
        (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
    const float u = (Etot * m_inv - 0.5f * v2);
#else
    const float u = (p->conserved.energy + p->flux.energy * dt_drift) * m_inv;
#endif
    p->P = hydro_gamma_minus_one * u * p->rho;
#endif
  }

    // MATTHIEU: Apply the entropy floor here.

#ifdef SWIFT_DEBUG_CHECKS
  if (p->h <= 0.0f) {
    error("Zero or negative smoothing length (%g)!", p->h);
  }
#endif

  gizmo_check_physical_quantities("density", "pressure", p->rho, p->v[0],
                                  p->v[1], p->v[2], p->P);
}

/**
 * @brief Set the particle acceleration after the flux loop
 *
 * We use the new conserved variables to calculate the new velocity of the
 * particle, and use that to derive the change of the velocity over the particle
 * time step.
 *
 * If the particle time step is zero, we set the accelerations to zero. This
 * should only happen at the start of the simulation.
 *
 * @param p Particle to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part* p, const struct cosmology* cosmo) {

  /* set the variables that are used to drift the primitive variables */

  // MATTHIEU: Bert is this correct? Do we need cosmology terms here?

  /* Add normalization to h_dt. */
  p->force.h_dt *= p->h * hydro_dimension_inv;
}

/**
 * @brief Extra operations done during the kick
 *
 * @param p Particle to act upon.
 * @param xp Extended particle data to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part* p, struct xpart* xp, float dt_therm, float dt_grav,
    float dt_hydro, float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props) {

  float a_grav[3];

  /* Update conserved variables (note: the mass does not change). */
  p->conserved.momentum[0] += p->flux.momentum[0] * dt_therm;
  p->conserved.momentum[1] += p->flux.momentum[1] * dt_therm;
  p->conserved.momentum[2] += p->flux.momentum[2] * dt_therm;
#if defined(EOS_ISOTHERMAL_GAS)
  /* We use the EoS equation in a sneaky way here just to get the constant u */
  p->conserved.energy =
      p->conserved.mass * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  p->conserved.energy += p->flux.energy * dt_therm;
#endif

#ifndef HYDRO_GAMMA_5_3
  const float Pcorr = (dt_hydro - dt_therm) * p->geometry.volume;
  p->conserved.momentum[0] -= Pcorr * p->gradients.P[0];
  p->conserved.momentum[1] -= Pcorr * p->gradients.P[1];
  p->conserved.momentum[2] -= Pcorr * p->gradients.P[2];
#endif

  /* Apply the minimal energy limit */
  const float min_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (p->conserved.energy < min_energy * p->conserved.mass) {
    p->conserved.energy = min_energy * p->conserved.mass;
    p->flux.energy = 0.0f;
  }

  // MATTHIEU: Apply the entropy floor here.

  gizmo_check_physical_quantities(
      "mass", "energy", p->conserved.mass, p->conserved.momentum[0],
      p->conserved.momentum[1], p->conserved.momentum[2], p->conserved.energy);

#ifdef SWIFT_DEBUG_CHECKS
  /* Note that this check will only have effect if no GIZMO_UNPHYSICAL option
     was selected. */
  if (p->conserved.energy < 0.0f) {
    error(
        "Negative energy after conserved variables update (energy: %g, "
        "denergy: %g)!",
        p->conserved.energy, p->flux.energy);
  }
#endif

  /* Add gravity. We only do this if we have gravity activated. */
  if (p->gpart) {
    /* Retrieve the current value of the gravitational acceleration from the
       gpart. We are only allowed to do this because this is the kick. We still
       need to check whether gpart exists though.*/
    a_grav[0] = p->gpart->a_grav[0];
    a_grav[1] = p->gpart->a_grav[1];
    a_grav[2] = p->gpart->a_grav[2];

    /* Kick the momentum for half a time step */
    /* Note that this also affects the particle movement, as the velocity for
       the particles is set after this. */
    p->conserved.momentum[0] += dt_grav * p->conserved.mass * a_grav[0];
    p->conserved.momentum[1] += dt_grav * p->conserved.mass * a_grav[1];
    p->conserved.momentum[2] += dt_grav * p->conserved.mass * a_grav[2];
  }

  /* Set the velocities: */
  /* We first set the particle velocity */
  if (p->conserved.mass > 0.0f && p->rho > 0.0f) {

    const float inverse_mass = 1.0f / p->conserved.mass;

    /* Normal case: set particle velocity to fluid velocity. */
    xp->v_full[0] = p->conserved.momentum[0] * inverse_mass;
    xp->v_full[1] = p->conserved.momentum[1] * inverse_mass;
    xp->v_full[2] = p->conserved.momentum[2] * inverse_mass;

  } else {
    /* Vacuum particles have no fluid velocity. */
    xp->v_full[0] = 0.0f;
    xp->v_full[1] = 0.0f;
    xp->v_full[2] = 0.0f;
  }

  if (p->gpart) {
    p->gpart->v_full[0] = xp->v_full[0];
    p->gpart->v_full[1] = xp->v_full[1];
    p->gpart->v_full[2] = xp->v_full[2];
  }

  /* reset wcorr */
  p->geometry.wcorr = 1.0f;
}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part* p, const struct xpart* xp) {}

#endif /* SWIFT_GIZMO_MFM_HYDRO_H */
