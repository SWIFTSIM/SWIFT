/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H
#define SWIFT_CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H

#include "../chemistry_properties.h"
#include "../chemistry_riemann_checks.h"
#include "../chemistry_riemann_utils.h"
#include "../chemistry_struct.h"
#include "../chemistry_timesteps.h"
#include "../parabolic/chemistry_riemann_HLL.h"
#include "hydro.h"

/**
 * @brief HLL riemann solver.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param hyperFluxL The flux of the hyperbolic conservation law of the left
 * (in physical units).
 * @param hyperFluxR The flux of the hyperbolic conservation law of the right
 * (in physical units).
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param fluxes (return) The resulting flux at the interface (in physical
 * units).
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_solver_HLL(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const double hyperFluxL[4][3],
    const double hyperFluxR[4][3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double fluxes[4]) {

  /***************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* PVRS wavespeed approximation */
  const double uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  const double c_diff_L =
      chemistry_get_physical_hyperbolic_soundspeed(pi, chem_data, cosmo);
  const double c_diff_R =
      chemistry_get_physical_hyperbolic_soundspeed(pj, chem_data, cosmo);
  const double lambda_plus = max(uL + c_diff_L, uR + c_diff_R);  /* S_R */
  const double lambda_minus = min(uL - c_diff_L, uR - c_diff_R); /* S_L */

  /* Degenerate/near-degenerate wave speeds: 1/(lambda_plus - lambda_minus)
     below would divide by ~0. Guard exactly like the parabolic solver
     (parabolic/chemistry_riemann_HLL.h) and return immediately -- falling
     through into the general HLL branch here would overwrite these zeros
     with NaN (0/0 from lprod * one_over_dl). */
  if (fabs(lambda_plus - lambda_minus) <
      GEAR_FVPM_DIFF_WAVESPEED_ESTIMATE_DIFFERENCE_TOLERANCE) {
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
    return;
  }

  /***************************************************************************/
  /* Now project the fluxes */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 4 quantities */
  double fluxL[4];
  fluxL[0] = hyperFluxL[0][0] * n_unit[0] + hyperFluxL[0][1] * n_unit[1] +
             hyperFluxL[0][2] * n_unit[2];
  fluxL[1] = hyperFluxL[1][0] * n_unit[0] + hyperFluxL[1][1] * n_unit[1] +
             hyperFluxL[1][2] * n_unit[2];
  fluxL[2] = hyperFluxL[2][0] * n_unit[0] + hyperFluxL[2][1] * n_unit[1] +
             hyperFluxL[2][2] * n_unit[2];
  fluxL[3] = hyperFluxL[3][0] * n_unit[0] + hyperFluxL[3][1] * n_unit[1] +
             hyperFluxL[3][2] * n_unit[2];

  double fluxR[4];
  fluxR[0] = hyperFluxR[0][0] * n_unit[0] + hyperFluxR[0][1] * n_unit[1] +
             hyperFluxR[0][2] * n_unit[2];
  fluxR[1] = hyperFluxR[1][0] * n_unit[0] + hyperFluxR[1][1] * n_unit[1] +
             hyperFluxR[1][2] * n_unit[2];
  fluxR[2] = hyperFluxR[2][0] * n_unit[0] + hyperFluxR[2][1] * n_unit[1] +
             hyperFluxR[2][2] * n_unit[2];
  fluxR[3] = hyperFluxR[3][0] * n_unit[0] + hyperFluxR[3][1] * n_unit[1] +
             hyperFluxR[3][2] * n_unit[2];

  /***************************************************************************/
  /* Now solve the Riemann problem */

  if (lambda_minus > 0.0) {
    for (int i = 0; i < 4; i++) {
      fluxes[i] = fluxL[i];
    }
  } else if (lambda_minus <= 0.0 && lambda_plus >= 0.0) {

    /************************************************************************/
    /* HLL flux */
    const double one_over_dl = 1.f / (lambda_plus - lambda_minus);
    const double lprod = lambda_plus * lambda_minus;
    const double F_transport[4] = {
        (lambda_plus * fluxL[0] - lambda_minus * fluxR[0]) * one_over_dl,
        (lambda_plus * fluxL[1] - lambda_minus * fluxR[1]) * one_over_dl,
        (lambda_plus * fluxL[2] - lambda_minus * fluxR[2]) * one_over_dl,
        (lambda_plus * fluxL[3] - lambda_minus * fluxR[3]) * one_over_dl};
    double F_diss[4] = {lprod * (UR[0] - UL[0]) * one_over_dl,
                        lprod * (UR[1] - UL[1]) * one_over_dl,
                        lprod * (UR[2] - UL[2]) * one_over_dl,
                        lprod * (UR[3] - UL[3]) * one_over_dl};

    /* Hopkins-2017-style limiter on the raw HLL dissipation, based on how
       far the dynamical flux still is from its own Fickian relaxation
       target ("flux memory") -- see chemistry_riemann_utils.h. Always
       limits the mass-density dissipation (component 0); whether it also
       limits the flux-component dissipation (1-3) is a parameter-file
       choice, since limiting flux memory too aggressively could preserve
       it instead of letting it decay. */
    const double limiter_factor =
        chemistry_riemann_compute_hyperbolic_diffusivity_limiter(
            pi, pj, m, chem_data, cosmo);
    F_diss[0] *= limiter_factor;
    if (chem_data->hyperbolic_limiter_scope == limiter_all_components) {
      F_diss[1] *= limiter_factor;
      F_diss[2] *= limiter_factor;
      F_diss[3] *= limiter_factor;
    }

    /* Causal bound: rescale F_diss by zeta = min(1, c_hyp*dt_ij/length), so
       its transport rate can never exceed c_hyp -- see the theory document,
       "A Causal Bound on the Dissipation". The length scale is the smaller of
       the two particles' physical sizes, which sets the von Neumann /
       Lax-Wendroff stability floor for this term. zeta is floored at a small
       strictly positive value: it must not reach zero when c_hyp -> 0
       (turbulent_mode with vanishing local shear), a legitimate flow state
       that still needs some dissipation to stay stable. As with the
       flux-memory limiter above, component 0 is always rescaled while
       components 1-3 follow hyperbolic_limiter_scope, since unconditionally
       suppressing their dissipation would preserve rather than decay an
       already-stale flux memory. dt_ij uses the same order-independent
       mindt as chemistry.c's chemistry_gear_fvpm_compute_pair_fluxes() and
       the MUSCL predictor (chemistry_gradients.h): the smallest positive
       flux.dt among the two particles, whichever slot it ends up in after
       canonicalisation. */
    const double c_max = max(c_diff_L, c_diff_R);
    const float dt_i = pi->chemistry_data.flux.dt;
    const float dt_j = pj->chemistry_data.flux.dt;
    double dt_min;
    if (dt_i > 0.f && dt_j > 0.f) {
      dt_min = min(dt_i, dt_j);
    } else if (dt_i > 0.f) {
      dt_min = dt_i;
    } else if (dt_j > 0.f) {
      dt_min = dt_j;
    } else {
      dt_min = 0.f;
    }
    const double psize_i = chemistry_get_physical_particle_size(pi, cosmo);
    const double psize_j = chemistry_get_physical_particle_size(pj, cosmo);
    const double length_scale = min(psize_i, psize_j);
    const double causal_floor = 0.01;
    double causal_factor = 1.0;
    if (dt_min > 0. && length_scale > 0.) {
      const double zeta = min(1.0, c_max * dt_min / length_scale);
      causal_factor = max(causal_floor, zeta);
    }
    F_diss[0] *= causal_factor;
    if (chem_data->hyperbolic_limiter_scope == limiter_all_components) {
      F_diss[1] *= causal_factor;
      F_diss[2] *= causal_factor;
      F_diss[3] *= causal_factor;
    }

    /* The causal bound above is this solver's only cap on F_diss. It
       replaces the psi-scaled minmod TVD bound used by the parabolic
       solver, which would otherwise collide with it: both cap the same
       dissipation, and at a tight psi the minmod bound binds first and
       makes the causal bound inert. hll_riemann_solver_psi therefore acts
       on the parabolic solver only, including the parabolic branch of the
       Berthon blend below. */
    double fluxes_HLL[4];
    for (int i = 0; i < 4; i++) {
      fluxes_HLL[i] = F_transport[i] + F_diss[i];
    }

    /* The pure hyperbolic HLL flux is unstable when tau -> 0, i.e. in the
       parabolic diffusion regime. To make the solution stable, we use the
       parabolic diffusion solver and we blend the two solutions together.
       The blending uses the ratio of the physical relaxation time and the
       numerical one. */
    const double alpha = chemistry_riemann_compute_hyperbolic_blending_factor(
        dx, pi, pj, UL, UR, m, chem_data, cosmo);

    /* Compute the parabolic diffusion solution */
    double F_par_L[3], F_par_R[3];
    chemistry_get_physical_parabolic_flux(pi, m, F_par_L, chem_data, cosmo);
    chemistry_get_physical_parabolic_flux(pj, m, F_par_R, chem_data, cosmo);
    double metal_flux_parabolic = 0.0;
    chemistry_riemann_solver_hopkins2017_HLL(
        dx, pi, pj, UL[0], UR[0], WL, WR, F_par_L, F_par_R, Anorm, n_unit, m,
        chem_data, cosmo, &metal_flux_parabolic);

    /* Now blend the fluxes, similarly to Berthon et al (2007)
       (https://link.springer.com/10.1007/s10915-006-9108-6).
       In the hyperbolic regime, alpha ~ 1, so the parabolic flux is
       negligible. In parabolic regime, the hyperbolic flux is negligible and
       the parabolic flux dominates. */
    fluxes[0] = alpha * fluxes_HLL[0] + (1.0 - alpha) * metal_flux_parabolic;
    fluxes[1] = alpha * fluxes_HLL[1];
    fluxes[2] = alpha * fluxes_HLL[2];
    fluxes[3] = alpha * fluxes_HLL[3];

    /* fluxes[0] = fluxes_HLL[0]; */
    /* fluxes[1] = fluxes_HLL[1]; */
    /* fluxes[2] = fluxes_HLL[2]; */
    /* fluxes[3] = fluxes_HLL[3]; */

    /* message( */
    /*     "lambda_minus = %e, lambda_plus = %e, alpha = %e, flux_HLL = (%e %e
     * %e %e)" */
    /*     " fluxes_par = (%e), final_fluxes = (%e %e %e %e)", lambda_minus, */
    /*     lambda_plus, alpha, fluxes_HLL[0], fluxes_HLL[1], fluxes_HLL[2], */
    /*     fluxes_HLL[3], metal_flux_parabolic, fluxes[0], fluxes[1], fluxes[2],
     * fluxes[3]); */

  } else if (lambda_plus < 0.0) {
    for (int i = 0; i < 4; i++) {
      fluxes[i] = fluxR[i];
    }
  }
}

/**
 * @brief Solve the Riemann problem for the diffusion equations and return the
 * flux at the interface.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param hyperFluxL The flux of the hyperbolic conservation law of the left
 * (in physical units).
 * @param hyperFluxR The flux of the hyperbolic conservation law of the right
 * (in physical units).
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param fluxes (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solve_for_flux(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const double hyperFluxL[4][3],
    const double hyperFluxR[4][3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double fluxes[4]) {

  chemistry_riemann_check_input(WL, WR, UL, UR, n_unit);

  /* Handle pure vacuum */
  if ((!UL[0] && !UR[0]) || (!WL[0] && !WR[0])) {
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
    return;
  }

  /* No conversion to physical needed, everything is physical here */
  chemistry_riemann_solver_HLL(dx, pi, pj, UL, UR, WL, WR, hyperFluxL,
                               hyperFluxR, Anorm, n_unit, m, chem_data, cosmo,
                               fluxes);

  chemistry_riemann_check_output(WL, WR, UL, UR, n_unit, fluxes);
}

#endif /* SWIFT_CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H */
