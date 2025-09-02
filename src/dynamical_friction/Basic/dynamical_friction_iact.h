/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_BASIC_DYNAMICAL_FRICTION_IACT_H
#define SWIFT_BASIC_DYNAMICAL_FRICTION_IACT_H


/**
 * @brief Density interaction between star and DM gravity particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param si First sparticle.
 * @param gpj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_dm_density(const float r2, const float dx[3],
                                 const float hi,
                                 struct spart *restrict si,
                                 const struct gpart *restrict gpj, const int with_cosmology, 
                                 const struct cosmology *cosmo) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->df_data.density_dm.wcount += wi;
  si->df_data.density_dm.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float mj = gpj->mass;

  /* Contribution to the total neighbour mass */
  si->df_data.dm_ngb_mass += mj;

  /* Contribution to the density at the GC */
  si->df_data.dm_rho += mj * wi;
  
  /* Neighbour's velocity in the frame of the GC */
  const float dv[3] = {gpj->v_full[0] - si->v[0], gpj->v_full[1] - si->v[1],
                      gpj->v_full[2] - si->v[2]};

  /* Contribution to the smoothed velocity (surrounding medium rel. to GC) */
  si->df_data.dm_v_medium[0] += mj * dv[0] * wi;
  si->df_data.dm_v_medium[1] += mj * dv[1] * wi;
  si->df_data.dm_v_medium[2] += mj * dv[2] * wi;

  float norm_v2;
  if (with_cosmology){
    const float a = cosmo->a;
    const float H = cosmo->H;
    const float a2H = a * a * H;

    /* Calculate the velocity with the Hubble flow */
    const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                    a2H * dx[2] + dv[2]};

    norm_v2 = v_plus_H_flow[0] * v_plus_H_flow[0] +
              v_plus_H_flow[1] * v_plus_H_flow[1] +
              v_plus_H_flow[2] * v_plus_H_flow[2];

  } else {
    norm_v2 = dv[0] * dv[0] +
              dv[1] * dv[1] +
              dv[2] * dv[2];
  }

  si->df_data.dm_sigma2 += norm_v2 * wi * mj;

}

/**
 * @brief Density interaction between two star particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param sj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_stars_density(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct spart *restrict si,
                                 const struct spart *restrict sj, const int with_cosmology, 
                                 const struct cosmology *cosmo) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->df_data.density_stars.wcount += wi;
  si->df_data.density_stars.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

}

#endif /* SWIFT_BASIC_DYNAMICAL_FRICTION_IACT_H */
