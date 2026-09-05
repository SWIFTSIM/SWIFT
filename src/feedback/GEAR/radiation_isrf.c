/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
/**
 * @file src/feedback/GEAR/radiation_isrf.c
 * @brief Receiver-side LW/FUV dust extinction and Yukawa propagation
 * physics for GEAR.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "chemistry.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "radiation.h"
#include "radiation_isrf.h"

#include <float.h>
#include <math.h>

/**
 * @brief First-init of a #part's LW/FUV radiation-field state. Shared
 * across GEAR feedback variants: independent of the injection mechanism.
 *
 * @param p The #part to initialise.
 */
void radiation_first_init_part(struct part *restrict p) {
  p->feedback_data.u_FUV = 0.f;
  p->feedback_data.u_LW = 0.f;
  p->feedback_data.u_FUV_prev = 0.f;
  p->feedback_data.u_LW_prev = 0.f;
  p->feedback_data.LW_FUV_last_touch_ti = -1;
  radiation_init_part_propagation(p, NULL);
}

/**
 * @brief Snapshot #u_FUV/#u_LW once per step, before the density loop's
 * h-iterations begin (see #feedback_part_data.u_FUV_prev).
 *
 * @param p The #part to reset.
 */
void radiation_snapshot_part_propagation(struct part *p) {
  p->feedback_data.u_FUV_prev = p->feedback_data.u_FUV;
  p->feedback_data.u_LW_prev = p->feedback_data.u_LW;
}

/**
 * @brief Zero the Yukawa propagation's per-h-iteration mixing
 * accumulators and cache this iteration's per-band absorption rate.
 * Mirrors chemistry_init_part's own per-iteration reset (called from the
 * same sites: part_init.h and the ghost h-iteration redo path), so it is
 * safe to call once or several times per step. Pure scratch space (no
 * restart I/O). A NULL #engine (first-init only) leaves kappa at 0.
 *
 * @param p The #part to reset.
 * @param e The #engine, or NULL.
 */
void radiation_init_part_propagation(struct part *p, const struct engine *e) {
  p->feedback_data.isrf_prop_sum_w_FUV = 0.f;
  p->feedback_data.isrf_prop_sum_wu_FUV = 0.f;
  p->feedback_data.isrf_prop_sum_w_LW = 0.f;
  p->feedback_data.isrf_prop_sum_wu_LW = 0.f;

  if (e == NULL || !e->feedback_props->LW_FUV_propagation) {
    p->feedback_data.kappa_FUV = 0.f;
    p->feedback_data.kappa_LW = 0.f;
    return;
  }

  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  p->feedback_data.kappa_FUV = radiation_get_part_linear_absorption_rate(
      e->internal_units, Z, rho_phys, RADIATION_SIGMA_D_FUV_CGS);
  p->feedback_data.kappa_LW = radiation_get_part_linear_absorption_rate(
      e->internal_units, Z, rho_phys, RADIATION_SIGMA_D_LW_CGS);
}

/**
 * @brief Apply one step of the Yukawa propagation update from the
 * accumulators radiation_propagation_iact.h filled during the density
 * loop, which runs before star feedback: this always stamps
 * LW_FUV_last_touch_ti, so injection adds on top instead of resetting.
 * Idempotent: always recomputed from the stable #u_FUV_prev snapshot and
 * this h-iteration's accumulators, so repeated calls across h-iterations
 * converge to the same answer regardless of how many there are. No
 * mixing term (just decay) with no gas neighbours; no-op when
 * propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_density_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->LW_FUV_propagation) return;

  const float w_min = e->feedback_props->LW_FUV_yukawa_w_min;
  const float h = p->h;
  const struct feedback_part_data *fd = &p->feedback_data;

  const float alpha_FUV =
      radiation_get_isrf_propagation_alpha(h, fd->kappa_FUV, w_min);
  const float alpha_LW =
      radiation_get_isrf_propagation_alpha(h, fd->kappa_LW, w_min);
  const float lambda2_FUV = 1.0f / max(fd->kappa_FUV * fd->kappa_FUV, FLT_MIN);
  const float lambda2_LW = 1.0f / max(fd->kappa_LW * fd->kappa_LW, FLT_MIN);
  const float decay_FUV = expf(-alpha_FUV * h * h / lambda2_FUV);
  const float decay_LW = expf(-alpha_LW * h * h / lambda2_LW);

  const float mixed_FUV =
      fd->isrf_prop_sum_w_FUV > 0.0f
          ? fd->isrf_prop_sum_wu_FUV / fd->isrf_prop_sum_w_FUV
          : fd->u_FUV_prev;
  const float mixed_LW = fd->isrf_prop_sum_w_LW > 0.0f
                             ? fd->isrf_prop_sum_wu_LW / fd->isrf_prop_sum_w_LW
                             : fd->u_LW_prev;

  p->feedback_data.u_FUV =
      (1.0f - alpha_FUV) * fd->u_FUV_prev * decay_FUV + alpha_FUV * mixed_FUV;
  p->feedback_data.u_LW =
      (1.0f - alpha_LW) * fd->u_LW_prev * decay_LW + alpha_LW * mixed_LW;
  p->feedback_data.LW_FUV_last_touch_ti = e->ti_current;
}

/**
 * Comoving gas column density at a gas particle's own location: the
 * receiver-side analogue of the star-side Sobolev column
 * (#radiation_get_comoving_gas_column_density_at_star), used for LW/FUV
 * extinction. Simplified to the kernel-radius fallback (no resolved
 * density gradient on the gas side).
 *
 * @param p The #part.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_part(const struct part *p) {
  const float h_gas = p->h * kernel_gamma;
  return 2.0f * h_gas * p->rho;
}

/**
 * Dust-to-gas ratio relative to the Milky Way, matching Grackle's own
 * default convention (dust2gas = fgr * Z/z_solar): D(Z)/D(Z_sun) =
 * Z/z_solar, independent of fgr itself.
 *
 * @param Z Gas metal mass fraction.
 * @return Dust-to-gas ratio relative to the Milky Way.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_to_gas_ratio_relative_to_MW(float Z) {
  return max(Z, 0.f) / RADIATION_GRACKLE_SOLAR_METAL_FRACTION;
}

/**
 * Band-specific dust mass opacity (area/mass, internal units): kappa_eff
 * = sigma_d_band * D(Z) / (mu_H * m_H).
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @return Dust mass opacity, internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_mass_opacity(const struct unit_system *us, float Z,
                                float sigma_d_band_cgs) {

  const float D_relative = radiation_get_dust_to_gas_ratio_relative_to_MW(Z);
  const float kappa_eff_cgs = sigma_d_band_cgs * D_relative /
                              (RADIATION_MU_H * RADIATION_HYDROGEN_MASS_CGS);
  return kappa_eff_cgs * units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
         units_cgs_conversion_factor(us, UNIT_CONV_AREA);
}

/**
 * Band-specific dust extinction factor: exp(-kappa_eff * Sigma_gas_p).
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us, float Z,
                                     float sigma_d_band_cgs,
                                     float Sigma_gas_p) {

  const float kappa_eff =
      radiation_get_dust_mass_opacity(us, Z, sigma_d_band_cgs);
  return expf(-kappa_eff * Sigma_gas_p);
}

/**
 * Local linear dust absorption rate (1/length): kappa_eff(Z) * rho, the
 * quantity the Yukawa screening length and the propagation interface
 * coupling's harmonic mean are built from. Purely local: no column
 * density involved, unlike the injection extinction above.
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param rho_p Physical gas density, internal units.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @return Local linear dust absorption rate, internal units (1/length).
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_linear_absorption_rate(const struct unit_system *us, float Z,
                                          float rho_p, float sigma_d_band_cgs) {
  return radiation_get_dust_mass_opacity(us, Z, sigma_d_band_cgs) * rho_p;
}

/**
 * Per-particle mixing fraction ceiling for the Yukawa propagation update.
 * Evaluated non-circularly: x = exp(-h^2/lambda^2) uses the natural
 * unit-alpha diffusive scale, not the alpha being solved for; alpha =
 * min(1, (1+x)/(w_min+x)) then bounds the worst-case response of the
 * normalized mixing operator to a checkerboard perturbation. The decay
 * factor actually applied uses this alpha: exp(-alpha*h^2/lambda^2).
 *
 * @param h Comoving smoothing length.
 * @param kappa_i Local linear absorption rate, 1/length.
 * @param w_min See #radiation_compute_yukawa_w_min.
 * @return Mixing fraction alpha, in (0, 1].
 */
__attribute__((always_inline)) INLINE float
radiation_get_isrf_propagation_alpha(float h, float kappa_i, float w_min) {
  const float lambda2 = 1.0f / max(kappa_i * kappa_i, FLT_MIN);
  const float x = expf(-h * h / lambda2);
  return min(1.0f, (1.0f + x) / (w_min + x));
}

/**
 * @brief Receiver-side LW/FUV dust extinction factors for a gas particle.
 *
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param extinction_FUV (return) FUV-band extinction factor.
 * @param extinction_LW (return) Lyman-Werner-band extinction factor.
 */
__attribute__((always_inline)) INLINE void
radiation_get_part_LW_FUV_extinction_factors(const struct unit_system *us,
                                             const struct cosmology *cosmo,
                                             const struct part *p, float Z,
                                             float *extinction_FUV,
                                             float *extinction_LW) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p) * cosmo->a2_inv;

  *extinction_FUV = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_FUV_CGS, Sigma_gas_p);
  *extinction_LW = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p);
}

/* Maximum simple-cubic lattice neighbours considered by
   #radiation_compute_yukawa_w_min: generous for any kernel this codebase
   ships (largest support radius in use, Wendland C6, needs ~250). */
#define RADIATION_YUKAWA_LATTICE_MAX_NEIGHBOURS 2000

/**
 * Fourier symbol of the normalized kernel mixing operator w_ij =
 * W(r_ij,h)/sum_k W(r_ik,h) at wavevector k, on the lattice given by
 * offsets/weights (see #radiation_compute_yukawa_w_min).
 */
static float radiation_yukawa_what(const float offsets[][3],
                                   const float weights[], int n, float kx,
                                   float ky, float kz) {
  float s = 0.f;
  for (int i = 0; i < n; i++)
    s += weights[i] *
         cosf(kx * offsets[i][0] + ky * offsets[i][1] + kz * offsets[i][2]);
  return s;
}

/**
 * @brief Magnitude of the most negative Fourier response of the
 * kernel-normalized mixing operator, for this build's compiled-in kernel
 * and eta_neighbours: bounds how aggressively
 * #radiation_get_isrf_propagation_alpha can relax without the explicit
 * update oscillating. Measured once at start-up, on an idealized
 * simple-cubic lattice, by a coarse sweep of the operator's Fourier
 * symbol over the first Brillouin zone.
 *
 * @param hydro_props The runtime hydrodynamics scheme properties.
 * @return W_min > 0 (the operator's response ranges within [-W_min, 1]).
 */
float radiation_compute_yukawa_w_min(const struct hydro_props *hydro_props) {

  const float eta = hydro_props->eta_neighbours;
  const float h = eta;
  const float H = kernel_gamma * h;

  static float offsets[RADIATION_YUKAWA_LATTICE_MAX_NEIGHBOURS][3];
  static float weights[RADIATION_YUKAWA_LATTICE_MAX_NEIGHBOURS];
  int n = 0;
  float w_sum = 0.f;

  const int nmax = (int)ceilf(H) + 1;
  for (int ix = -nmax; ix <= nmax; ix++) {
    for (int iy = -nmax; iy <= nmax; iy++) {
      for (int iz = -nmax; iz <= nmax; iz++) {
        if (ix == 0 && iy == 0 && iz == 0) continue;
        const float r = sqrtf((float)(ix * ix + iy * iy + iz * iz));
        if (r >= H) continue;
        if (n >= RADIATION_YUKAWA_LATTICE_MAX_NEIGHBOURS)
          error(
              "Yukawa W_min lattice sum exceeded its fixed neighbour "
              "budget; raise RADIATION_YUKAWA_LATTICE_MAX_NEIGHBOURS.");
        float w;
        kernel_eval(r / h, &w);
        offsets[n][0] = (float)ix;
        offsets[n][1] = (float)iy;
        offsets[n][2] = (float)iz;
        weights[n] = w;
        w_sum += w;
        n++;
      }
    }
  }

  if (w_sum <= 0.f)
    error("Yukawa W_min: empty kernel support for eta_neighbours=%g", eta);
  for (int i = 0; i < n; i++) weights[i] /= w_sum;

  /* Coarse sweep of the Brillouin zone [-pi,pi]^3 (lattice spacing 1):
     the worst response of an isotropic, positive, decreasing kernel's
     mixing operator is a checkerboard-type high-frequency mode, so a
     modest grid already resolves it without needing gradient search. */
  const int n_grid = 12;
  float worst = 1.f;
  for (int a = 0; a < n_grid; a++) {
    const float kx = -(float)M_PI + 2.f * (float)M_PI * a / (n_grid - 1);
    for (int b = 0; b < n_grid; b++) {
      const float ky = -(float)M_PI + 2.f * (float)M_PI * b / (n_grid - 1);
      for (int c = 0; c < n_grid; c++) {
        const float kz = -(float)M_PI + 2.f * (float)M_PI * c / (n_grid - 1);
        const float what =
            radiation_yukawa_what(offsets, weights, n, kx, ky, kz);
        if (what < worst) worst = what;
      }
    }
  }

  return -worst;
}
