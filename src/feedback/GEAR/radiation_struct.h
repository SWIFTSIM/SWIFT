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
#ifndef SWIFT_FEEDBACK_GEAR_RADIATION_STRUCT_H
#define SWIFT_FEEDBACK_GEAR_RADIATION_STRUCT_H

/**
 * @file src/feedback/GEAR/radiation_struct.h
 * @brief Types of the GEAR subgrid radiation shared by the GEAR feedback
 * modules: the radiation policy, the ISRF moments and the star-side radiation
 * state.
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "timeline.h"

/*! Maximum number of HEALPix pixels (12*nside_max^2) of a star's HII budget,
    set by ./configure --with-number-of-hii-angular-pixels. */
#ifndef HII_MAX_ANGULAR_PIXELS
#error "HII_MAX_ANGULAR_PIXELS should be defined by configure (config.h)"
#endif

/**
 * @brief The ISRF moments: PE (6-11.2 eV) and Lyman-Werner (11.2-13.6 eV)
 * specific energy, and the Lyman-Werner photon-number moment.
 *
 * Indexes #feedback_part_data.isrf_moment and
 * #feedback_spart_data.radiation.L_band.
 *
 * #ISRF_MOMENT_LW_PHOTON uses the operator of #ISRF_MOMENT_LW.
 */
enum radiation_isrf_moment {
  ISRF_MOMENT_PE = 0,
  ISRF_MOMENT_LW,
  ISRF_MOMENT_LW_PHOTON,
  ISRF_MOMENT_COUNT
};

/**
 * @brief The band operators (dust opacity, M1 closure, dissipation) a moment
 * is evaluated against. Indexes #feedback_part_data.isrf_operator.
 */
enum radiation_isrf_operator {
  ISRF_OPERATOR_PE = 0,
  ISRF_OPERATOR_LW,
  ISRF_OPERATOR_COUNT
};

/**
 * @brief Operator used by each #radiation_isrf_moment.
 *
 * Unsized, so that a missing initialiser fails the _Static_assert below.
 */
static const enum radiation_isrf_operator radiation_isrf_moment_to_operator[] =
    {ISRF_OPERATOR_PE, ISRF_OPERATOR_LW, ISRF_OPERATOR_LW};

_Static_assert(sizeof(radiation_isrf_moment_to_operator) /
                       sizeof(radiation_isrf_moment_to_operator[0]) ==
                   ISRF_MOMENT_COUNT,
               "radiation_isrf_moment_to_operator needs one entry per "
               "ISRF_MOMENT_COUNT.");

/**
 * @brief Owning moment of each #radiation_isrf_operator, the reverse of
 * #radiation_isrf_moment_to_operator.
 *
 * Writers of an #feedback_isrf_operator_data field loop over operators and
 * use this map, so a shared operator is written once. Each entry MUST be the
 * lowest-index moment mapping to that operator, see
 * #feedback_check_isrf_operator_owner_map(). Unsized, like the forward map.
 */
static const enum radiation_isrf_moment radiation_isrf_operator_owner[] = {
    ISRF_MOMENT_PE, ISRF_MOMENT_LW};

_Static_assert(sizeof(radiation_isrf_operator_owner) /
                       sizeof(radiation_isrf_operator_owner[0]) ==
                   ISRF_OPERATOR_COUNT,
               "radiation_isrf_operator_owner needs one entry per "
               "ISRF_OPERATOR_COUNT.");

/**
 * @brief The subgrid radiation feedback processes.
 */
enum radiation_policy {
  radiation_policy_none = 0,
  /*! Photoionization (Strömgren sphere). */
  radiation_policy_photoionization = (1 << 0),
  /*! Radiation pressure from the stars' bolometric luminosity */
  radiation_policy_radiation_pressure = (1 << 1),

  /*! Interstellar radiation field: photoelectric heating and H2
   * photodissociation. */
  radiation_policy_isrf = (1 << 2),
};

/**
 * @brief Per-moment ISRF transport state of a hydro particle.
 */
struct feedback_isrf_moment_data {

  /*! Local specific radiation field of the band (physical, per unit mass).
      Double: the per-step relaxation depth can be ~1e-8, too small for
      `1 - a` in float. */
  double u;

  /*! Snapshot of #u taken once per step. Double, like #u. */
  double u_prev;

  /*! Reduced specific flux `F_true/c_hyp`, in the units of #u. */
  float specific_flux[3];

  /*! `(1/rho) div(rho F)` accumulated by the force loop. */
  float div_specific_flux;

  /*! Artificial-dissipation source term from the force loop. Scratch. */
  float dissipation_u;

  /*! M1 pressure-tensor divergence `(1/rho) div(D rho u)` from the gradient
      loop, per physical length. Scratch. */
  float grad_u[3];

  /*! Mass-specific emission dose not yet injected into #u. */
  float u_dose_reservoir;

  /*! This step's source rate, drawn from #u_dose_reservoir. Scratch. */
  float u_source_rate;

  /*! Transport amount owed by finer neighbours, divided by c_hyp. Added and
      zeroed by the next own update. */
  float pending_transport_u;

  /*! Dissipation amount owed by finer neighbours, divided by c_hyp and
      weighted by their phi. Same life cycle as #pending_transport_u. */
  float pending_dissipation_u;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Most negative #u at the end of an update since the previous snapshot,
      0 if none. */
  float u_min_since_snapshot;

  /*! Cumulative dose handed to #u by the reservoir, rescaled by `c_hyp/c`
      but not relaxed. */
  float cumulative_injected;

  /*! Cumulative energy lost to relaxation (dust absorption and redshift),
      including transport and dissipation terms of optically thick steps. */
  float cumulative_absorbed;
#endif
};

/**
 * @brief Per-operator ISRF band coefficients of a hydro particle.
 */
struct feedback_isrf_operator_data {

  /*! Local linear dust absorption rate of the band, cached once per step. Not
      the injection-side kappa of #radiation_get_part_ISRF_extinction_factors.
   */
  float kappa;

  /*! M1 closure tensor `D(f)` of the owning moment. */
  float m1_closure_D[3][3];

  /*! Kernel mean of the neighbours' `|rho_prev*u_prev|`, the reference of the
      negativity trigger. Scratch. */
  float ngb_mean_abs_u_V;

  /*! Reactive part of the dissipation coefficient, raised by the negativity
      trigger and decayed otherwise. Dimensionless. */
  float dissipation_alpha_trigger;

  /*! Anticipatory part, the `h/lambda`-gated floor
      (#radiation_dissipation_alpha_floor_band). */
  float dissipation_alpha_floor;
};

/**
 * @brief Radiation pressure momentum received by a hydro particle.
 */
struct feedback_xpart_radiation_data {

  /*! Momentum received from radiation pressure */
  float delta_p[3];
};

/**
 * @brief HII ionization payload of a hydro particle, computed by the owner of
 * the particle.
 */
struct feedback_xpart_HII_region_data {

  /*! Mean photon energy of the tagging star above 13.6 eV, frozen at tag
      time, in erg. Only set with GEARFeedback:HII_couple_ionization_rate. */
  float excess_photon_energy_HI;

  /*! Photoionization rate coefficient Gamma_HI of the tagging star at this
      particle, frozen at tag time (internal 1/time). Only set with
      GEARFeedback:HII_couple_ionization_rate, 0 otherwise. */
  float photoionization_rate_HI;
};

/**
 * @brief Radiation state carried by each star particle.
 */
struct feedback_spart_radiation_data {

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Bolometric luminosity (physical units) from the stellar evolution */
  double L_bol;
#endif

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Ionizing photon rate of each active pixel (physical units): the star's
      total, split evenly. Double, since the value is huge. A pure rate,
      never debited. */
  double dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS];

  /*! dot_N_ion_pix at the previous HII rebuild pass, for the trapezoid rule
      of radiation_open_ionizing_photon_budget(). Negative until the first
      pass. */
  double dot_N_ion_pix_prev[HII_MAX_ANGULAR_PIXELS];

  /*! Photon count spendable per pixel this rebuild pass, `dot_N_ion_pix *
      dt_back`, debited by radiation_consume_ionizing_photons. */
  double N_ion_budget_pix[HII_MAX_ANGULAR_PIXELS];

  /*! Number of active angular pixels this star is currently using
      (1 = spherical/HEALPix disabled) */
  int n_HII_pixels;

  /*! Mass in the HII region generated by this star particle */
  float mass_HII_region;

  /*! Co-moving HII region radius when the star died or stopped being
      eligible. */
  float final_HII_radius;

  /*! Ionized gas mass of the final HII region. */
  float final_HII_mass;

  /*! Star age at this HII region's last rebuild pass. Double, since
      dt_back (age now minus this) scales the whole photon budget above. */
  double HII_region_last_rebuild;

  /*! Star age at the last rebuild attempt, whether or not gas was found.
      Anchors the interval dt_back. */
  double HII_region_last_attempt;

  /*! Mean photon energy above 13.6 eV, cached once per HII rebuild pass, in
      erg. Only set with GEARFeedback:HII_couple_ionization_rate. */
  float mean_excess_photon_energy_HI;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /*! Moment luminosity (physical units) by #radiation_isrf_moment, from the
      radiation table. 0 unless the interstellar radiation field is on. */
  double L_band[ISRF_MOMENT_COUNT];
#endif

#ifdef GEAR_SUBGRID_RADIATION
  /*! Photospheric effective temperature (internal units), a diagnostic no
      feedback channel uses. For a population particle, the hottest
      surviving star. 0 without a "Teff" table dataset. */
  float teff;

  /*! This star's feedback timestep (proper time, internal units), cached once
      per step by feedback_prepare_radiation_feedback. */
  float Delta_t;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Integer time at the start of the step #Delta_t was cached for. */
  integertime_t Delta_t_cached_ti_begin;
#endif
#endif
};

#endif /* SWIFT_FEEDBACK_GEAR_RADIATION_STRUCT_H */
