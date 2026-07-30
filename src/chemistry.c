/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "chemistry.h"

/**
 * @brief Initialises the chemistry properties.
 *
 * Calls chemistry_init_backend for the chosen chemistry function.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
void chemistry_init(struct swift_params *parameter_file,
                    const struct unit_system *us,
                    const struct phys_const *phys_const,
                    struct chemistry_global_data *data) {

  chemistry_init_backend(parameter_file, us, phys_const, data);
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * Calls chemistry_print_backend for the chosen chemistry model.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
void chemistry_print(const struct chemistry_global_data *data) {
  chemistry_print_backend(data);
}

/**
 * @brief Write a chemistry struct to the given FILE as a stream of bytes.
 *
 * @param chemistry the struct
 * @param stream the file stream
 */
void chemistry_struct_dump(const struct chemistry_global_data *chemistry,
                           FILE *stream) {
  restart_write_blocks((void *)chemistry, sizeof(struct chemistry_global_data),
                       1, stream, "chemistry", "chemistry function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param chemistry the struct
 * @param stream the file stream
 */
void chemistry_struct_restore(const struct chemistry_global_data *chemistry,
                              FILE *stream) {
  restart_read_blocks((void *)chemistry, sizeof(struct chemistry_global_data),
                      1, stream, NULL, "chemistry function");
}

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#include "chemistry/GEAR_FVPM_DIFFUSION/chemistry_iact.h"

/**
 * @brief Computes the diffusion pair fluxes of all elements (side-effect
 * free).
 *
 * This is the computation part of runner_iact_chemistry_fluxes_common(),
 * shared by the FCT prep pass and the flux-applying force loop. It is
 * compiled exactly once and never inlined so that both passes provably run
 * identical code: the two callers live in different translation units, and
 * under fast-math/LTO two inlined copies would be free to round differently
 * and flip discontinuous limiter branches.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i (not modified).
 * @param pj Particle j (not modified).
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param totflux (return) Per-element fluxes across the interface.
 * @param mindt_out (return) Time step of the flux exchange.
 * @param vmax_out (return) Maximal signal velocity of the pair (hyperbolic
 * scheme; 0 otherwise).
 * @return a CHEMISTRY_PAIR_FLUX_* status describing which outputs are valid.
 */
__attribute__((noinline)) int chemistry_gear_fvpm_compute_pair_fluxes(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *restrict pi, const struct part *restrict pj,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo,
    double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4], float *mindt_out,
    float *vmax_out
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
    ,
    float *delxbar_i_out, float *delxbar_j_out
#endif
) {

  *mindt_out = 0.f;
  *vmax_out = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  *delxbar_i_out = 0.f;
  *delxbar_j_out = 0.f;
#endif

  /* If the masses are null, then there is nothing to diffuse. */
  if ((hydro_get_mass(pi) == 0.0 || hydro_get_mass(pj) == 0) ||
      (pi->chemistry_data.kappa == 0.0 && pj->chemistry_data.kappa == 0.0)) {
    return CHEMISTRY_PAIR_FLUX_NONE;
  }

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Initialize local variables */
  float Bi[3][3];
  float Bj[3][3];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  const float Vi = pi->geometry.volume;
  const float Vj = pj->geometry.volume;

#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  /* Calculate the maximal diffusion speed */
  const float ci =
      chemistry_get_physical_hyperbolic_soundspeed(pi, chem_data, cosmo);
  const float cj =
      chemistry_get_physical_hyperbolic_soundspeed(pj, chem_data, cosmo);
  float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];
  dvdr *= r_inv;
  dvdr *= cosmo->a_inv; /* Convert comoving velocities to physical units */
  *vmax_out = ci + cj - min(0.0, dvdr);
#endif

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute (square of) area */
  /* eqn. (7) */
  float Anorm2 = 0.0f;
  float A[3];
  if (fvpm_part_geometry_well_behaved(pi) &&
      fvpm_part_geometry_well_behaved(pj)) {
    /* in principle, we use Vi and Vj as weights for the left and right
     * contributions to the generalized surface vector.
     * However, if Vi and Vj are very different (because they have very
     * different smoothing lengths), then the expressions below are more
     * stable. */
    float Xi = Vi;
    float Xj = Vj;
#ifdef GIZMO_VOLUME_CORRECTION
    if (fabsf(Vi - Vj) / min(Vi, Vj) > 1.5f * hydro_dimension) {
      Xi = (Vi * hj + Vj * hi) / (hi + hj);
      Xj = Xi;
    }
#endif
    for (int k = 0; k < 3; k++) {
      /* we add a minus sign since dx is pi->x - pj->x */
      A[k] = -Xi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) *
                 wi * hi_inv_dim -
             Xj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) *
                 wj * hj_inv_dim;
      Anorm2 += A[k] * A[k];
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    const float hidp1 = pow_dimension_plus_one(hi_inv);
    const float hjdp1 = pow_dimension_plus_one(hj_inv);
    const float Anorm =
        -(hidp1 * Vi * Vi * wi_dx + hjdp1 * Vj * Vj * wj_dx) * r_inv;
    A[0] = -Anorm * dx[0];
    A[1] = -Anorm * dx[1];
    A[2] = -Anorm * dx[2];
    Anorm2 = Anorm * Anorm * r2;
  }

  /* if the interface has no area, nothing happens and we return */
  /* continuing results in dividing by zero and NaN's... */
  if (Anorm2 == 0.0f) {
    return CHEMISTRY_PAIR_FLUX_VMAX;
  }

  /* Compute the area */
  const float Anorm_inv = 1.0f / sqrtf(Anorm2);
  const float Anorm = Anorm2 * Anorm_inv;

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /* For stability reasons, we do require A and dx to have opposite
   * directions (basically meaning that the surface normal for the surface
   * always points from particle i to particle j, as it would in a real
   * moving-mesh code). If not, our scheme is no longer upwind and hence can
   * become unstable. */
  const float dA_dot_dx = A[0] * dx[0] + A[1] * dx[1] + A[2] * dx[2];
  /* In GIZMO, Phil Hopkins reverts to an SPH integration scheme if this
   * happens. We curently just ignore this case and display a message. */
  const float rdim = pow_dimension(r);
  if (dA_dot_dx > 1.e-6f * rdim) {
    error("Ill conditioned gradient matrix (%g %g %g %g %g)!", dA_dot_dx, Anorm,
          Vi, Vj, r);
  }
#endif

#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  /* Lanson & Vila (2008), eq. (58): Delta x_i = 1 / (2 * sum_l w_l ||A_il||),
   * where w_l is the neighbour's volume and ||A_il|| is the one-sided
   * renormalized-gradient face-area weight. The geometric equivalent of
   * each summand w_l ||A_il|| is the (symmetrized) interface area Anorm
   * divided by the neighbour's volume, Anorm / V_neighbour; the factor of
   * 1/2 in front of the sum in eq. (58) is applied once, in
   * chemistry_timesteps.h, to the completed sum delxbar. */
  *delxbar_i_out = Anorm / Vj;
  *delxbar_j_out = Anorm / Vi;
#endif

  /* Compute the normal vector of the interface */
  const float n_unit[3] = {A[0] * Anorm_inv, A[1] * Anorm_inv,
                           A[2] * Anorm_inv};

  /* Compute interface position (relative to pi, since we don't need
   * the actual position) eqn. (8) */
  const float xfac = -hi / (hi + hj);
  const float xij_i[3] = {xfac * dx[0], xfac * dx[1], xfac * dx[2]};

  /* Get the time step for the flux exchange. This is always the smallest time
   * step among the two particles. */
  const float mindt =
      (chj->flux.dt > 0.f) ? fminf(chi->flux.dt, chj->flux.dt) : chi->flux.dt;
  *mindt_out = mindt;

  /* Nothing to do */
  if (mindt == 0.f) {
    return CHEMISTRY_PAIR_FLUX_GEOMETRY;
  }

  /*****************************************/
  /* Predict the velocity at the interface to compute fluxes */
  /* Get the hydro W_L and W_R */
  float vi[3] = {pi->v[0], pi->v[1], pi->v[2]};
  float vj[3] = {pj->v[0], pj->v[1], pj->v[2]};

  /* Compute interface velocity */
  const float vij[3] = {vi[0] + (vi[0] - vj[0]) * xfac,
                        vi[1] + (vi[1] - vj[1]) * xfac,
                        vi[2] + (vi[2] - vj[2]) * xfac};

  /* Get the primitive variable of Euler eq */
  float Wi[5] = {hydro_get_comoving_density(pi), vi[0], vi[1], vi[2],
                 hydro_get_comoving_pressure(pi)};
  float Wj[5] = {hydro_get_comoving_density(pj), vj[0], vj[1], vj[2],
                 hydro_get_comoving_pressure(pj)};

  chemistry_gradients_predict_hydro(pi, pj, dx, r, xij_i, Wi, Wj);

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  /* Note: This is necessary to properly follow the fluid motion. */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* Convert to physical units */
  Wi[0] *= cosmo->a3_inv;
  Wi[1] /= cosmo->a;
  Wi[2] /= cosmo->a;
  Wi[3] /= cosmo->a;
  Wi[4] *= cosmo->a_factor_pressure;

  Wj[0] *= cosmo->a3_inv;
  Wj[1] /= cosmo->a;
  Wj[2] /= cosmo->a;
  Wj[3] /= cosmo->a;
  Wj[4] *= cosmo->a_factor_pressure;

  /*****************************************/
  /* Now solve the Riemann problem for each metal specie */
  /* Helper variable */
  const float a2 = cosmo->a * cosmo->a;
  const float Anorm_p = a2 * Anorm;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {

    /* Predict the diffusion state at the interface to compute fluxes */
    double Ui[4], Uj[4];
    chemistry_gradients_predict(pi, pj, m, dx, r, xij_i, cosmo, chem_data, Ui,
                                Uj);
    /* Note: The returned values are in physical units. No conversion needed */

    /* Solve the 1D Riemann problem at the interface A_ij _physical units_ */
    totflux[m][0] = 0.0;
    totflux[m][1] = 0.0;
    totflux[m][2] = 0.0;
    totflux[m][3] = 0.0;
    chemistry_compute_flux(dx, pi, pj, m, Ui, Uj, Wi, Wj, n_unit, Anorm_p,
                           chem_data, cosmo, totflux[m]);

    /* Flux limiter */
    /* First check that we won't have negative masses. If so, we have a check
       that will ensure masses are not negative and if so, it we set them to be
       positive. Then, we have metal mass creation. If this correction happen
       a lot, we will create a lot of metal mass. */
    /* Note: the limiter only reads interaction_mode for a debug message, so
       passing 0 here keeps the computation mode-independent. */
    chemistry_limit_metal_mass_flux(pi, pj, m, totflux[m], mindt,
                                    /*interaction_mode=*/0, chem_data);
  }

  return CHEMISTRY_PAIR_FLUX_FLUXES;
}

#endif /* GEAR FVPM diffusion chemistry */
