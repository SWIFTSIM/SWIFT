

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
#include "radiation_h2_shielding.h"
#include "radiation_isrf.h"
#include "radiation_isrf_extinction.h"
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
  fd->ISRF_H2_line_absorption = 0;
  /* An LW field of the ICs has no source: it dissociates at the constant. */
  fd->LW_dissociating_u = fd->isrf_moment[ISRF_MOMENT_LW].u;
  fd->H2_column_per_length = 0.f;
  fd->H2_shielding_exponent = 0.f;
  fd->H2_line_width_km_s = 0.f;
  fd->H2_velocity_limit = 0.f;
  fd->LW_line_H2_column = 0.f;
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
  fd->ISRF_H2_line_absorption =
      (e->feedback_props->ISRF_H2_line_absorption ? RADIATION_ISRF_H2_LINES
                                                  : 0) |
      (e->feedback_props->ISRF_packet_constants.LW_population_track
           ? RADIATION_ISRF_LW_POPULATION
           : 0);
}

/**
 * @brief Cache the H2 gas properties that the hop of a share and the merge of
 * a slot read, at the drifted state of the particle: the column per length of
 * path, and the fit if asked.
 *
 * The temperature and the mean molecular weight are read only for the fit.
 * Needs a cooling network that tracks H2. Without one the column is 0.
 *
 * @param p The #part.
 * @param xp The #xpart of the particle.
 * @param e The #engine.
 * @param set_fit Also set the fit (exponent, line width, velocity limit)?
 */
static void radiation_cache_H2_line_absorption(struct part *p,
                                               const struct xpart *xp,
                                               const struct engine *e,
                                               const int set_fit) {
  struct feedback_part_data *fd = &p->feedback_data;
  double X_H2 = 0., mu = 0., T = 0.;
#if defined(GEAR_COOLING) && COOLING_NETWORK > 1
  X_H2 = cooling_get_H2_mass_fraction_subgrid(xp);
  if (set_fit) {
    mu = cooling_get_mean_molecular_weight(
        e->physical_constants, e->internal_units, e->cosmology,
        e->hydro_properties, e->cooling_func, p, xp);
    T = cooling_get_temperature(e->physical_constants, e->hydro_properties,
                                e->internal_units, e->cosmology,
                                e->cooling_func, p, xp);
  }
#endif
  const struct unit_system *us = e->internal_units;
  const double rho_cgs = (double)hydro_get_physical_density(p, e->cosmology) *
                         units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  const double m_H_cgs = e->physical_constants->const_proton_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double k_B_cgs = e->physical_constants->const_boltzmann_k *
                         units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
                         units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  double column_per_length, exponent, b5, velocity_limit;
  radiation_h2_shielding_receiver_cache(
      X_H2, rho_cgs, mu, T, m_H_cgs, k_B_cgs,
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH),
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY), &column_per_length,
      &exponent, &b5, &velocity_limit);
  fd->H2_column_per_length = (float)column_per_length;
  if (!set_fit) return;
  fd->H2_shielding_exponent = (float)exponent;
  fd->H2_line_width_km_s = (float)b5;
  fd->H2_velocity_limit = (float)velocity_limit;
}

/**
 * @brief Once per step, before the density loop: cache the absorption rates,
 * the comoving density and `dt_prev`, and merge the packet.
 *
 * Runs at drift time, the last point where the metal fraction and `p->rho`
 * still hold the previous step's converged values.
 *
 * @param p The #part to reset.
 * @param xp The #xpart of the particle.
 * @param e The #engine.
 */
void radiation_snapshot_part_propagation(struct part *p, const struct xpart *xp,
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
  const float D_relative = cooling_get_dust_to_gas_ratio_relative_to_MW(
      e->cooling_func, Z, rho_phys);
  p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, D_relative, rho_phys,
          RADIATION_SIGMA_D_PE_CGS);
  p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, D_relative, rho_phys,
          RADIATION_SIGMA_D_LW_CGS);

  /* The fit is frozen over the particle's step (#H2_shielding_exponent). A
     particle with no fit yet takes it before its test beams and its merge. */
  const int with_lines = e->feedback_props->ISRF_H2_line_absorption;
  const int first_fit =
      with_lines && !(p->feedback_data.H2_line_width_km_s > 0.f);
  if (with_lines) radiation_cache_H2_line_absorption(p, xp, e, first_fit);

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

  const struct radiation_packet_constants *packet_constants =
      &e->feedback_props->ISRF_packet_constants;

  /* The drift runs for inactive particles too. Give them the step lattice
   * now, so that light arriving before their first merge is placed once, on
   * the final lattice. */
  if (!part_is_active(p, e) && p->time_bin > 0)
    radiation_packet_set_first_step(
        p, e->ti_current, ti_step,
        get_integer_time_end(e->ti_current, p->time_bin), packet_constants);

  if (part_is_active(p, e)) {
    if (!p->feedback_data.isrf_packet.has_merged_step)
      radiation_packet_inject_test_beams(
          p, e->feedback_props->ISRF_packet_test_beam,
          e->feedback_props->ISRF_packet_test_beam_origin_distance_h,
          e->ti_current, packet_constants);
    radiation_packet_merge(p, e->ti_current, ti_step, e->time_base,
                           with_cosmology ? e->cosmology : NULL,
                           packet_constants);
    radiation_packet_reset_field_mean(p);

    if (with_lines && !first_fit) {
      const float old_alpha = p->feedback_data.H2_shielding_exponent;
      const float old_b5 = p->feedback_data.H2_line_width_km_s;
      radiation_cache_H2_line_absorption(p, xp, e, /*set_fit=*/1);
      if (p->feedback_data.H2_shielding_exponent != old_alpha ||
          p->feedback_data.H2_line_width_km_s != old_b5)
        radiation_packet_line_refit(p, old_alpha, old_b5);
    }
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
  radiation_packet_couple_field(p, &e->feedback_props->ISRF_packet_constants);
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
      e->feedback_props->ISRF_packet_read_mean_time,
      &e->feedback_props->ISRF_packet_constants);
}
