

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
 * @brief Receiver-side LW/PE dust extinction and per-step hooks of the packet
 * transport of the ISRF for GEAR.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "radiation.h"
#include "radiation_isrf.h"
#include "radiation_packet.h"
#include "radiation_propagation_iact.h"
#include "timeline.h"

#include <float.h>
#include <math.h>
#include <string.h>

/**
 * @brief First-init of a #part's LW/PE radiation-field state.
 *
 * Does not zero #feedback_isrf_moment_data.u, which an IC may supply through
 * the optional "PESpecificEnergy"/"LWSpecificEnergy" fields.
 *
 * @param p The #part to initialise.
 */
void radiation_first_init_part(struct part *restrict p) {
  struct feedback_part_data *fd = &p->feedback_data;
  /* The LW and LW photon IC fields are optional and independent. A zero one
   * is set equal to the other, in either direction, so the two stay
   * nonzero together. */
  if (fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u == 0.f &&
      fd->isrf_moment[ISRF_MOMENT_LW].u != 0.f) {
    fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u =
        fd->isrf_moment[ISRF_MOMENT_LW].u;
  } else if (fd->isrf_moment[ISRF_MOMENT_LW].u == 0.f &&
             fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u != 0.f) {
    fd->isrf_moment[ISRF_MOMENT_LW].u =
        fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u;
  }
  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
    fd->isrf_operator[o].kappa = 0.f;
  fd->ISRF_last_touch_ti = -1;
  /* MPI foreign particles are not bzero'd: keep them from reading as
     illuminated. */
  fd->is_illuminated_ISRF = 0;
  fd->ISRF_illumination_end_ti = -1;
  /* Placeholder density: a 0 seed would divide by zero before the first
   * snapshot. */
  fd->rho_prev = 1.0f;
  fd->dt_prev = 0.f;
  fd->ISRF_packet_transport = 0;
  memset(&fd->isrf_packet, 0, sizeof(struct feedback_isrf_packet_data));
}

/**
 * @brief Copy the transport choice onto the particle, so the neighbour loops
 * read it without the engine.
 *
 * @param p The #part.
 * @param e The #engine.
 */
static void radiation_set_part_transport(struct part *p,
                                         const struct engine *e) {
  struct feedback_part_data *fd = &p->feedback_data;
  fd->ISRF_packet_transport = e->feedback_props->ISRF_propagation;
  fd->isrf_packet.cone_power = e->feedback_props->ISRF_packet_cone_power;
  memcpy(fd->isrf_packet.octant_frame,
         e->feedback_props->ISRF_packet_octant_frame,
         sizeof(fd->isrf_packet.octant_frame));
}

/**
 * @brief Once per step, before the density loop: cache the absorption rates,
 * the comoving density and `dt_prev`, and merge the packet.
 *
 * Runs at drift time, the last point where the metal fraction and `p->rho`
 * still hold the previous step's converged values.
 *
 * @param p The #part to reset.
 * @param e The #engine.
 */
void radiation_snapshot_part_propagation(struct part *p,
                                         const struct engine *e) {
  /* Cached even without ISRF_propagation: it must never be 0. */
  const float rho_comoving = hydro_get_comoving_density(p);
  p->feedback_data.rho_prev = rho_comoving > 0.f ? rho_comoving : 1.0f;

  radiation_set_part_transport(p, e);
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
    p->feedback_data.isrf_packet.slot[s].weight_sum = 0.f;
    p->feedback_data.isrf_packet.slot[s].forward_advance = 0.f;
  }

  if (!e->feedback_props->ISRF_propagation) {
    p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa = 0.f;
    p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa = 0.f;
    return;
  }

  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  /* Resolved value, not the raw `-1` sentinel. */
  const float local_dust_to_gas_ratio =
      (float)e->cooling_func->chemistry_data.local_dust_to_gas_ratio;
  p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, Z, rho_phys,
          RADIATION_SIGMA_D_PE_CGS, local_dust_to_gas_ratio);
  p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, Z, rho_phys,
          RADIATION_SIGMA_D_LW_CGS, local_dust_to_gas_ratio);

  /* Physical timestep of this particle. Can be exactly 0 before its first
   * step: cached unfloored so readers can detect that case. */
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_step = get_integer_timestep(p->time_bin);
  /* Start of the step ending now. Using `ti_current` instead would
   * over-drain one sub-step early. */
  const integertime_t ti_begin =
      get_integer_time_begin(e->ti_current, p->time_bin);
  float dt_phys;
  if (with_cosmology) {
    dt_phys = (float)cosmology_get_delta_time(e->cosmology, ti_begin,
                                              ti_begin + ti_step);
  } else {
    dt_phys = (float)get_timestep(p->time_bin, e->time_base);
  }
  p->feedback_data.dt_prev = dt_phys;

  /* The drift runs for inactive particles too: they keep their inbox. */
  if (part_is_active(p, e)) {
    if (p->feedback_data.isrf_packet.step_end_ti == 0)
      radiation_packet_inject_test_beams(
          p, e->feedback_props->ISRF_packet_test_beam, e->ti_current);
    radiation_packet_merge(p, e->ti_current, ti_step, e->time_base,
                           with_cosmology ? e->cosmology : NULL);
    radiation_packet_reset_field_mean(p);
  }
}

/**
 * @brief Radiation timestep bound `f_dt a h_i / v` of the packet transport,
 * with `f_dt` the time step margin
 * (#feedback_props.ISRF_packet_light_speed_time_step_margin) and `v` the limit
 * speed (#feedback_props.ISRF_packet_light_speed_fraction_of_c times c).
 *
 * A packet makes one hop per own step, so its light is too slow unless the
 * step is short. Only the `fixed` light speed mode asks for that, for every
 * gas particle. The cap and the unlimited mode give FLT_MAX, and so does
 * ISRF_propagation off. Above `f_dt` of about 0.8 the bulk of the light lags
 * behind `v` while its front keeps the speed.
 *
 * @param p The #part to consider.
 * @param e The #engine.
 * @return The radiation timestep bound, or FLT_MAX if none applies.
 */
float radiation_isrf_part_timestep(const struct part *restrict p,
                                   const struct engine *e) {
  if (!e->feedback_props->ISRF_propagation) return FLT_MAX;
  if (e->feedback_props->ISRF_packet_light_speed !=
      isrf_packet_light_speed_fixed)
    return FLT_MAX;
  const float h_hop = (float)e->cosmology->a * p->h;
  if (h_hop <= 0.f) return FLT_MAX;
  const float c_limit =
      e->feedback_props->ISRF_packet_light_speed_fraction_of_c *
      (float)e->physical_constants->const_speed_light_c;
  return e->feedback_props->ISRF_packet_light_speed_time_step_margin * h_hop /
         c_limit;
}

/**
 * @brief Close the packet of the particle after its force loop and write the
 * kernel mean of its field into the band moments Grackle reads.
 *
 * No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_force_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  radiation_packet_end_force(p);
  radiation_packet_couple_field(p);
}

/**
 * @brief Close the gradient pass of the packet transport, once per step in the
 * extra ghost. No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_gradient_propagation(struct part *p,
                                        const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  /* Also set here: a cell may skip the drift before its first force loop. */
  radiation_set_part_transport(p, e);
  radiation_packet_end_gradient(
      p, (float)e->cosmology->a,
      (float)e->physical_constants->const_speed_light_c, fd->dt_prev,
      e->feedback_props->ISRF_packet_read_mean_time);
}

/**
 * @brief Comoving path length of the receiver-side LW/PE dust extinction
 * column, chosen by #isrf_extinction_path_mechanism.
 *
 * The arguments are the union over all mechanisms; each case uses a subset.
 * The call is inside the star-gas pair loop. The temperature_capped_jeans
 * mechanism is costly there (one temperature and a square root per pair).
 *
 * @param fb_props Properties of the feedback scheme.
 * @param p The receiving #part.
 * @param xp The receiving #xpart (tracked species, for the temperature).
 * @param r Comoving separation of the illuminating star and the receiver.
 * @param cosmo The current cosmological model.
 * @param phys_const The physical constants.
 * @param hydro_props The hydro scheme properties.
 * @param us Unit system.
 * @param cooling The cooling function properties.
 * @return Comoving extinction path length.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_extinction_path(
    const struct feedback_props *fb_props, const struct part *p,
    const struct xpart *xp, const float r, const struct cosmology *cosmo,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling) {

  /* One support radius: the longest admissible path, and the fallback for a
   * degenerate gas state. */
  const float h_gas = p->h * kernel_gamma;

  switch ((enum isrf_extinction_path_mechanism)
              fb_props->ISRF_extinction_path_mechanism) {

    case isrf_extinction_path_constant_kernel_path:
      return fb_props->ISRF_extinction_path_in_kernel_radii * h_gas;

    case isrf_extinction_path_pair_separation:
      return r;

    case isrf_extinction_path_temperature_capped_jeans: {
      const float T = cooling_get_temperature(phys_const, hydro_props, us,
                                              cosmo, cooling, p, xp);
      const float rho_phys = hydro_get_physical_density(p, cosmo);

      /* Return the cap: a clamped division is unsafe under -ffast-math. */
      if (T <= 0.f || rho_phys <= 0.f) return h_gas;

      /* c_s^2 scales with T, so the cap multiplies it by T_cap/T. */
      const float cs = hydro_get_physical_soundspeed(p, cosmo);
      const float T_ratio =
          min(1.f, fb_props->ISRF_extinction_jeans_temperature_cap_K / T);
      const float cs2_capped = cs * cs * T_ratio;

      /* lambda_J = sqrt(pi c_s^2 / (G rho)), physical, then comoving. */
      const float lambda_J_phys = sqrtf(
          (float)(M_PI * cs2_capped / (phys_const->const_newton_G * rho_phys)));
      return min(lambda_J_phys * cosmo->a_inv, h_gas);
    }
  }

  error("Unknown GEARFeedback:ISRF_extinction_path mechanism %d.",
        (int)fb_props->ISRF_extinction_path_mechanism);
  return 0.f;
}

/**
 * @brief Comoving gas column density at a gas particle: its own density times
 * the extinction path.
 *
 * @param p The #part.
 * @param extinction_path Comoving path length, from
 * #radiation_get_comoving_extinction_path.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_part(const struct part *p,
                                                  const float extinction_path) {
  return extinction_path * p->rho;
}

/**
 * @brief Dust-to-gas ratio relative to the Milky Way, following Grackle:
 * `(fgr/fgr_default) * Z/z_solar`.
 *
 * See #RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO.
 *
 * @param Z Gas metal mass fraction.
 * @param local_dust_to_gas_ratio Resolved chemistry_data value, not the raw
 * `-1` sentinel.
 * @return Dust-to-gas ratio relative to the Milky Way.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_to_gas_ratio_relative_to_MW(float Z,
                                               float local_dust_to_gas_ratio) {
  const float D_relative_Z =
      max(Z, 0.f) / RADIATION_GRACKLE_SOLAR_METAL_FRACTION;
  return D_relative_Z * (local_dust_to_gas_ratio /
                         (float)RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);
}

/**
 * @brief Band dust mass opacity in internal units (area/mass):
 * `sigma_d_band * D(Z) / (mu_H * m_H)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants (for the proton mass).
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust mass opacity, internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_mass_opacity(const struct unit_system *us,
                                const struct phys_const *phys_const, float Z,
                                float sigma_d_band_cgs,
                                float local_dust_to_gas_ratio) {

  const float D_relative = radiation_get_dust_to_gas_ratio_relative_to_MW(
      Z, local_dust_to_gas_ratio);
  return sigma_d_band_cgs * D_relative /
         (RADIATION_MU_H * phys_const->const_proton_mass *
          units_cgs_conversion_factor(us, UNIT_CONV_AREA));
}

/**
 * @brief Band dust extinction factor `exp(-kappa_eff * Sigma_gas_p)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us,
                                     const struct phys_const *phys_const,
                                     float Z, float sigma_d_band_cgs,
                                     float Sigma_gas_p,
                                     float local_dust_to_gas_ratio) {

  const float kappa_eff = radiation_get_dust_mass_opacity(
      us, phys_const, Z, sigma_d_band_cgs, local_dust_to_gas_ratio);
  return expf(-kappa_eff * Sigma_gas_p);
}

/**
 * @brief Local linear dust absorption rate `kappa_eff(Z) * rho` in 1/length.
 * It is local: no column density is involved.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param Z Gas metal mass fraction.
 * @param rho_p Physical gas density, internal units.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Local linear dust absorption rate, internal units (1/length).
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_linear_absorption_rate(const struct unit_system *us,
                                          const struct phys_const *phys_const,
                                          float Z, float rho_p,
                                          float sigma_d_band_cgs,
                                          float local_dust_to_gas_ratio) {
  return radiation_get_dust_mass_opacity(us, phys_const, Z, sigma_d_band_cgs,
                                         local_dust_to_gas_ratio) *
         rho_p;
}

/**
 * @brief Receiver-side LW/PE dust extinction factors for a gas particle.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param cooling The cooling function properties (for the resolved
 * chemistry_data.local_dust_to_gas_ratio).
 * @param extinction_path Comoving extinction path length, from
 * #radiation_get_comoving_extinction_path.
 * @param extinction (return) Extinction factor of each operator, indexed by
 * #radiation_isrf_operator.
 */
__attribute__((always_inline)) INLINE void
radiation_get_part_ISRF_extinction_factors(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct part *p, float Z,
    const struct cooling_function_data *cooling, const float extinction_path,
    float extinction[ISRF_OPERATOR_COUNT]) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p, extinction_path) *
      cosmo->a2_inv;
  /* Resolved value, not the raw `-1` sentinel. */
  const float local_dust_to_gas_ratio =
      (float)cooling->chemistry_data.local_dust_to_gas_ratio;

  extinction[ISRF_OPERATOR_PE] = radiation_get_dust_extinction_factor(
      us, phys_const, Z, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p,
      local_dust_to_gas_ratio);
  extinction[ISRF_OPERATOR_LW] = radiation_get_dust_extinction_factor(
      us, phys_const, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p,
      local_dust_to_gas_ratio);
}
