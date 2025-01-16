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
 ******************************************************************************/
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"

/**
 * @brief Compute the particle parabolic timestep proportional to h^2.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_parabolic_timestep(
    const struct part *restrict p,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

#if defined(GEAR_MF_HYPERBOLIC_DIFFUSION)
  const float CFL_condition = chem_data->C_CFL_chemistry;
  const float delta_x = cosmo->a * kernel_gamma * p->h;

  /* CFL condition */
  const float dt_cfl =
      CFL_condition * delta_x /
      chemistry_get_physical_hyperbolic_soundspeed(p, chem_data, cosmo);
  return dt_cfl;
#else
  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute the diffusion matrix K */
  double K[3][3];
  chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...), and q = (Z_1, Z_2,
     ...). Hence, the term norm(U)/norm(q) in eq (15) is abs(rho). */
  const float norm_U_over_norm_q =
      cosmo->a3_inv * fabs(chemistry_get_comoving_density(p));

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h * cosmo->a;
  float norm_q = 0.0;
  float norm_nabla_q = 0.0;
  float expression = 0.0;

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0) {
    return FLT_MAX;
  }

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_q += chemistry_get_metal_mass_fraction(p, i) *
              chemistry_get_metal_mass_fraction(p, i);

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q = Grad Z */
      norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt. Add the missing a^{-1} to have physical gradients */
  norm_q = sqrtf(norm_q);  // Physical
  norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv;

  /* If ||q|| is 0 (metal density = 0), then there is no metal to diffuse.
     Hence, no timestep limit is needed. */
  if (norm_q == 0) {
    return FLT_MAX;
  }

  if (chem_data->diffusion_mode == isotropic_constant) {
    /* Isotropic constant diffusion has the simple expression: */
    expression = delta_x;
  } else {
    /* Compute the expression in the square bracket in eq (15). Notice that I
       rewrote it to avoid division by 0 when norm_nabla_q = 0. */
    expression = norm_q * delta_x / (norm_nabla_q * delta_x + norm_q);
  }

  return expression * expression / norm_matrix_K * norm_U_over_norm_q;
#endif
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * This is equation (10) in Alexiades, Amiez and Gremaud (1996).
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param timestep_explicit Usual explicit parabolic timestep.
 */
__attribute__((always_inline)) INLINE static float chemistry_get_supertimestep(
    const struct part *restrict p, const struct chemistry_global_data *cd,
    float timestep_explicit) {
  const float N = cd->N_substeps;
  const float nu = cd->nu;
  const float nu_plus_term = pow(1 + sqrtf(nu), 2.0 * N);
  const float nu_minus_term = pow(1 - sqrtf(nu), 2.0 * N);
  const float left_term =
      (nu_plus_term - nu_minus_term) / (nu_plus_term + nu_minus_term);
  return timestep_explicit * N / (2.0 * sqrt(nu)) * left_term;
}

/**
 * @brief Compute the particle supertimestep with using a CFL-like condition,
 * proportioanl to h.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_advective_supertimestep(
    const struct part *restrict p, const struct chemistry_global_data *cd,
    const struct cosmology *cosmo) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute diffusion matrix K */
  double K[3][3];
  chemistry_get_physical_matrix_K(p, cd, cosmo, K);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h * cosmo->a;

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...). */
  float norm_U = 0.0;
  float norm_nabla_q = 0.0;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_U += chemistry_get_physical_metal_density(p, i, cosmo) *
              chemistry_get_physical_metal_density(p, i, cosmo);

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q */
      norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt and convert to physical units */
  norm_U = sqrtf(norm_U);
  norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv;

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0 || norm_nabla_q == 0.0 || norm_U == 0.0) {
    /* The firt two cases are derived from the equation. The last case comes
    from a physical argument: if there is no metal, there is no need to limit
    the timestep */
    return FLT_MAX;
  }

  return cd->C_CFL_chemistry * delta_x * norm_U /
         (norm_matrix_K * norm_nabla_q);
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param current_substep_number Current substep number.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_subtimestep(const struct part *restrict p,
                              const struct chemistry_global_data *cd,
                              int current_substep_number) {
  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* TODO: Check that the division is not an integer division */
  const float cos_argument =
      M_PI * (2.0 * current_substep_number - 1.0) / (2.0 * cd->N_substeps);
  const float expression = (1 + cd->nu) - (1 - cd->nu) * cos(cos_argument);
  return chd->timesteps.explicit_timestep / expression;
}

/**
 * @brief Compute a valid integer time-step form a given time-step
 *
 * TODO: This is a copy/paste from make_integer_timestep. This is bad. Improve
 * it.
 *
 * We consider the minimal time-bin of any neighbours and prevent particles
 * to differ from it by a fixed constant `time_bin_neighbour_max_delta_bin`.
 *
 * If min_ngb_bin is set to `num_time_bins`, then no limit from the neighbours
 * is imposed.
 *
 * @param new_dt The time-step to convert.
 * @param old_bin The old time bin.
 * @param min_ngb_bin Minimal time-bin of any neighbour of this particle.
 * @param ti_current The current time on the integer time-line.
 * @param time_base_inv The inverse of the system's minimal time-step.
 */
__attribute__((always_inline, const)) INLINE static integertime_t
chemistry_make_integer_timestep(const float new_dt, const timebin_t old_bin,
                                const timebin_t min_ngb_bin,
                                const integertime_t ti_current,
                                const double time_base_inv) {

  /* Convert to integer time */
  integertime_t new_dti = (integertime_t)(new_dt * time_base_inv);

  /* Are we allowed to use this bin given the neighbours? */
  timebin_t new_bin = get_time_bin(new_dti);
  new_bin = min(new_bin, min_ngb_bin + time_bin_neighbour_max_delta_bin);
  new_dti = get_integer_timestep(new_bin);

  /* Current time-step */
  const integertime_t current_dti = get_integer_timestep(old_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, old_bin);

  /* Limit timestep increase */
  if (old_bin > 0) new_dti = min(new_dti, 2 * current_dti);

  /* Put this timestep on the time line */
  integertime_t dti_timeline = max_nr_timesteps;
  while (new_dti < dti_timeline) dti_timeline /= ((integertime_t)2);
  new_dti = dti_timeline;

  /* Make sure we are allowed to increase the timestep size */
  if (new_dti > current_dti) {
    if ((max_nr_timesteps - ti_end) % new_dti > 0) new_dti = current_dti;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (new_dti == 0) error("Computed an integer time-step of size 0");
#endif

  return new_dti;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H  */
