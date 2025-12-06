//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H

#include "hydro_getters.h"
#include "riemann.h"

/**
 * @brief Reset the fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_fluxes(
    struct part* restrict p) {

  p->flux.mass = 0.0f;
  p->flux.momentum[0] = 0.0f;
  p->flux.momentum[1] = 0.0f;
  p->flux.momentum[2] = 0.0f;
  p->flux.energy = 0.0f;
  p->flux.entropy = 0.0f;
  p->flux.dt = -1.f;

  p->gravity.mflux[0] = 0.0f;
  p->gravity.mflux[1] = 0.0f;
  p->gravity.mflux[2] = 0.0f;
}

/**
 * @brief Get the fluxes for the given particle.
 *
 * @param p Particle.
 * @param flux Fluxes for the particle (array of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_fluxes(
    const struct part* restrict p, float* flux) {

  flux[0] = p->flux.mass;
  flux[1] = p->flux.momentum[0];
  flux[2] = p->flux.momentum[1];
  flux[3] = p->flux.momentum[2];
  flux[4] = p->flux.energy;
  flux[5] = p->flux.entropy;
}

/**
 * @brief Compute the flux for the Riemann problem with the given left and right
 * state, and interface normal, surface area and velocity.
 *
 * @param WL Left state variables.
 * @param WR Right state variables.
 * @param n_unit Unit vector of the interface.
 * @param vLR Velocity of the interface.
 * @param Anorm Surface area of the interface.
 * @param fluxes Array to store the result in (of size 5 or more).
 * @return The maximal mach number of the signal speed of the shockwaves in the
 * solution of the Riemann problem (if any) or zero.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_flux(
    const float* WL, const float* WR, const float* n_unit, const float* vLR,
    const float Anorm, float* fluxes) {

  const float P_star = riemann_solve_for_flux(WL, WR, n_unit, vLR, fluxes);
  const float mach_number = riemann_get_max_mach_number(WL, WR, P_star);
  float entropy_flux;
  if (fluxes[0] > 0.f) {
    entropy_flux = fluxes[0] * WL[5];
  } else {
    entropy_flux = fluxes[0] * WR[5];
  }

  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
  fluxes[4] *= Anorm;
  fluxes[5] = Anorm * entropy_flux;

  return mach_number;
}

__attribute__((always_inline)) INLINE static void
hydro_part_positivity_limiter_fluxes(const struct part* pi,
                                     const struct part* pj, const float* n_unit,
                                     const float* vLR, const float Anorm,
                                     float epsilon_rho, float epsilon_P,
                                     float* fluxes) {

  epsilon_rho = fmin(fmin(pi->rho, pj->rho), epsilon_rho);
  epsilon_P = fmin(fmin(pi->P, pj->P), epsilon_P);

  double V_inv_i = 1. / pi->geometry.volume;
  double V_inv_j = 1. / pi->geometry.volume;
  double flux_fac_i =  pi->flux.dt * pi->geometry.nface;
  double flux_fac_j =  pj->flux.dt * pj->geometry.nface;
  //double flux_fac_i = (1. + 1e-5) * pi->geometry.area / Anorm * pi->flux.dt;
  //double flux_fac_j = (1. + 1e-5) * pj->geometry.area / Anorm * pj->flux.dt;
  double m_dagger_i = pi->conserved.mass - flux_fac_i * fluxes[0];
  double m_dagger_j = pj->conserved.mass + flux_fac_j * fluxes[0];
  double rho_i = m_dagger_i * V_inv_i;
  double rho_j = m_dagger_j * V_inv_j;

  /* Anything to do here? */
  if (rho_i > epsilon_rho && rho_j > epsilon_rho) return;

  float Wi[6], Wj[6];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);
  Wi[1] -= vLR[0];
  Wi[2] -= vLR[1];
  Wi[3] -= vLR[2];
  Wj[1] -= vLR[0];
  Wj[2] -= vLR[1];
  Wj[3] -= vLR[2];
  float fluxes_lo[6];
  hydro_compute_flux(Wi, Wj, n_unit, vLR, Anorm, fluxes_lo);
  double theta_rho = 1.;

  if (rho_i < epsilon_rho) {

    double rho_lo = (pi->conserved.mass - flux_fac_i * fluxes_lo[0]) * V_inv_i;
    theta_rho = fmin(theta_rho, (epsilon_rho - rho_lo) / (rho_i - rho_lo));
  }

  if (rho_j < epsilon_rho) {

    double rho_lo = (pj->conserved.mass + flux_fac_j * fluxes_lo[0]) * V_inv_j;
    theta_rho = fmin(theta_rho, (epsilon_rho - rho_lo) / (rho_j - rho_lo));
  }

  theta_rho = fmax(0., theta_rho);
  for (int k = 0; k < 6; ++k) {

    fluxes[k] = fluxes_lo[k] * (1. - theta_rho) + fluxes[k] * theta_rho;
  }

  float Qi[6], Qj[6];
  hydro_part_get_conserved_variables(pi, Qi);
  hydro_part_get_conserved_variables(pj, Qj);

  for (int k = 0; k < 6; ++k) {
    Qi[k] -= flux_fac_i * fluxes[k];
    Qj[k] += flux_fac_j * fluxes[k];
  }

  /* dE acts on entire cell, should not be multiplied by flux_fac */
  //Qi[4] -= pi->gravity.dE_prev;
  //Qj[4] -= pj->gravity.dE_prev;

  double mi_inv = 1. / Qi[0];
  double mj_inv = 1. / Qj[0];

  double ui =
      (Qi[4] - 0.5 * (Qi[1] * Qi[1] + Qi[2] * Qi[2] + Qi[3] * Qi[3]) * mi_inv) *
      mi_inv;

  double uj =
      (Qj[4] - 0.5 * (Qj[1] * Qj[1] + Qj[2] * Qj[2] + Qj[3] * Qj[3]) * mj_inv) *
      mj_inv;

  double Pi = gas_pressure_from_internal_energy(Qi[0] * V_inv_i, ui);
  double Pj = gas_pressure_from_internal_energy(Qj[0] * V_inv_j, uj);

  /* Anything to do here? */
  if (Pi > epsilon_P && Pj > epsilon_P) return;
  double theta_P = 1.;
  if (Pi < epsilon_P) {

    float Q_lo[6];
    hydro_part_get_conserved_variables(pi, Q_lo);

    for (int k = 0; k < 6; ++k) {
      Q_lo[k] -= flux_fac_i * fluxes_lo[k];
    }
    /* dE acts on entire cell, should not be multiplied by flux_fac */
    //Q_lo[4] -= pi->gravity.dE_prev;

    double m_inv = 1. / Q_lo[0];
    double u_lo =
        (Q_lo[4] -
         0.5 * (Q_lo[1] * Q_lo[1] + Q_lo[2] * Q_lo[2] + Q_lo[3] * Q_lo[3]) *
             m_inv) *
        m_inv;

    double P_lo = gas_pressure_from_internal_energy(Q_lo[0] * V_inv_i, u_lo);
    theta_P = fmin(theta_P, (epsilon_P - P_lo) / (Pi - P_lo));
  }
  if (Pj < epsilon_P) {

    float Q_lo[6];
    hydro_part_get_conserved_variables(pj, Q_lo);

    for (int k = 0; k < 6; ++k) {
      Q_lo[k] += flux_fac_j * fluxes_lo[k];
    }

    /* dE acts on entire cell, should not be multiplied by flux_fac */
    //Q_lo[4] -= pj->gravity.dE_prev;

    double m_inv = 1. / Q_lo[0];
    double u_lo =
        (Q_lo[4] -
         0.5 * (Q_lo[1] * Q_lo[1] + Q_lo[2] * Q_lo[2] + Q_lo[3] * Q_lo[3]) *
             m_inv) *
        m_inv;

    double P_lo = gas_pressure_from_internal_energy(Q_lo[0] * V_inv_j, u_lo);
    theta_P = fmin(theta_P, (epsilon_P - P_lo) / (Pj - P_lo));
  }

  theta_P = fmax(0., theta_P);

  for (int k = 0; k < 6; ++k) {

    fluxes[k] = fluxes_lo[k] * (1. - theta_P) + fluxes[k] * theta_P;
  }

  /* As an absolute last resort, just rescale fluxes if necessary */
  double theta = 1.;
  if (pi->conserved.mass <
      pi->geometry.area / Anorm * pi->flux.dt * fluxes[0]) {
    theta = fmin(theta, pi->conserved.mass / (flux_fac_i * fluxes[0]));
  }
  if (pj->conserved.mass <
      -pj->geometry.area / Anorm * pj->flux.dt * fluxes[0]) {
    theta = fmin(theta, pj->conserved.mass / (-flux_fac_j * fluxes[0]));
  }
  if (pi->conserved.energy <
      flux_fac_i * fluxes[4]) {
    theta = fmin(theta, (pi->conserved.energy) /
      (flux_fac_i * fluxes[4]));
  }
  if (pj->conserved.energy <
      -flux_fac_j * fluxes[4]) {
    theta = fmin(theta, (pj->conserved.energy) /
      (-flux_fac_j * fluxes[4]));
  }
  theta = fmax(0., theta);
  if (theta < 1.) {

    for (int k = 0; k < 6; ++k) {
      fluxes[k] = fluxes_lo[k];
    }
  }
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Time integrated fluxes across the interface.
 * @param dx Vector pointing from right particle to left particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_fluxes_left(
    struct part* restrict p, const float* restrict fluxes, const float dt) {

  p->flux.mass -= fluxes[0] * dt;
  p->flux.momentum[0] -= fluxes[1] * dt;
  p->flux.momentum[1] -= fluxes[2] * dt;
  p->flux.momentum[2] -= fluxes[3] * dt;
  p->flux.energy -= fluxes[4] * dt;
  p->flux.entropy -= fluxes[5] * dt;

  p->flux_count -= 1;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Time integrated fluxes across the interface.
 * @param dx Vector pointing from right particle to left particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_update_fluxes_right(struct part* restrict p,
                               const float* restrict fluxes, const float dt) {

  p->flux.mass += fluxes[0] * dt;
  p->flux.momentum[0] += fluxes[1] * dt;
  p->flux.momentum[1] += fluxes[2] * dt;
  p->flux.momentum[2] += fluxes[3] * dt;
  p->flux.energy += fluxes[4] * dt;
  p->flux.entropy += fluxes[5] * dt;

  p->flux_count += 1;
}

/**
 * @brief Get the drift term for the density based on the mass flux.
 *
 * @param mass_flux Current mass flux for the particle.
 * @param dt Drift time step (in co-moving units).
 * @param volume Volume of the particle.
 * @return Term that will drift the density.
 */
__attribute__((always_inline)) INLINE static float
hydro_flux_density_drift_term(const float mass_flux, const float dt,
                              const float volume) {

  if (volume > 0.0f) {
    return mass_flux * dt / volume;
  } else {
    return 0.0f;
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H
