/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "../GEAR/feedback_tracers_struct.h"
#include "../GEAR/radiation_struct.h"
#include "chemistry_struct.h"
#include "timeline.h"

/**
 * @brief Feedback fields carried by each hydro particle.
 *
 * Lives in #part, not #xpart, so the HII tag follows MPI exchange and
 * restarts.
 */
struct feedback_part_data {

  /* Members are ordered from large to small alignment. The padding of struct
     part depends on whether the MPI HII floats are compiled in. */

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Id of the star that ionized this particle. */
  long long star_id;

  /*! Simulation time until which this particle stays flagged as ionized. */
  double end_time;
#endif

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  /*! Photons debited from the claiming star's budget for this claim. Read by
      the owner-side merge when this claim loses a collision, to count the
      forfeited budget. */
  double cost;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /*! Per-moment ISRF state, indexed by #radiation_isrf_moment. */
  struct feedback_isrf_moment_data isrf_moment[ISRF_MOMENT_COUNT];

  /*! Step (#engine.ti_current) at which #feedback_isrf_moment_data.u was
      last written. With ISRF_propagation off, a match with the current step
      means a further star sums into u, a mismatch means u is zeroed first.
      Set to -1 at first init. */
  integertime_t ISRF_last_touch_ti;

  /*! Integer time until which #is_illuminated_ISRF stays set, renewed on
      every injection touch with a margin of
      RADIATION_ISRF_TAG_LIFETIME_INTERVALS star steps. Set to -1 at first
      init. */
  integertime_t ISRF_illumination_end_ti;
#endif

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /* Trace the metals received from feedback events. This is similar to not
     diffusing metals */
  double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

#ifdef GEAR_COOLING
  /*! Neutral hydrogen mass fraction cached by every cooling call (grackle_0:
      1.0f). Read only for a foreign copy, which has no #xpart. Set to -1 at
      first init, so "not yet primed" is unambiguous. */
  float neutral_H_frac;
#endif

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  /*! HII eligibility temperature cached by the owner's cooling call, for a
      foreign copy with no #xpart. 0 before the first cooling call, which
      reads as cold: an unprimed copy of hot gas then reads as cold. */
  float T_eligibility;

  /*! Mean molecular weight cached by the same cooling call as
      #T_eligibility, so the two describe one internal energy. Set to -1 at
      first init. */
  float mu_eligibility;

  /*! Squared comoving distance to the star of the current claim, at claim
      time. The owner-side merge keeps the claim with the smallest
      (r2, star_id). */
  float r2;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /*! Per-operator ISRF band-physics coefficients, indexed by
      #radiation_isrf_operator; a moment reaches its own operator through
      #radiation_isrf_moment_to_operator. */
  struct feedback_isrf_operator_data isrf_operator[ISRF_OPERATOR_COUNT];

  /*! LW specific energy weighted by the H2 line transmission, the fraction of
      the dissociation rate that the lines let through; the kernel mean of the
      packet field. Written with #isrf_moment, read for the H2 dissociation
      rate with H2 line absorption on. */
  double LW_dissociating_u;

  /*! H2 column per unit internal length of path, in cm^-2, of the gas of this
      particle, cached at every drift. */
  float H2_column_per_length;

  /*! Exponent of the H2 self-shielding fit at the gas density and temperature
      of this particle. The fit (this, #H2_line_width_km_s and
      #H2_velocity_limit) is frozen over the particle's step: it is set after
      the merge of an active drift, so every share summed into a slot between
      two merges uses the fit that the next merge inverts with. */
  float H2_shielding_exponent;

  /*! Thermal line width of H2 in km/s, part of the fit. 0 until the first fit
      and if the temperature was not positive. */
  float H2_line_width_km_s;

  /*! Relative velocity in internal units beyond which light no longer sees
      the H2 column it carries, part of the fit. */
  float H2_velocity_limit;

  /*! Comoving density snapshot of the previous step, read by the packet
      transport. Seeded to 1 at first init, never 0. */
  float rho_prev;

  /*! This particle's own physical timestep, cached once per step at the
      drift. */
  float dt_prev;

  /*! Radiation packet of the packet transport. */
  struct feedback_isrf_packet_data isrf_packet;
#endif

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Tag to mark the particle as ionized. */
  char is_ionized;
#endif

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  /*! Set when this rank claims a foreign copy of this particle, so the claim
      is reported to the owner. Cleared when the report is packed. */
  char claimed_this_pass;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /*! Is the particle illuminated by a star's injection pass, with the
      episode still live? Gates a first-touch-only timestep_sync_part call and
      is cleared once #ISRF_illumination_end_ti lapses, like #is_ionized. */
  char is_illuminated_ISRF;

  /*! Copy of #feedback_props.ISRF_propagation, set at drift and in the extra
      ghost, so the neighbour loops read it without the engine. */
  char ISRF_packet_transport;

  /*! Bits of the LW dissociation transport, set at the drift:
      #RADIATION_ISRF_H2_LINES (#feedback_props.ISRF_H2_line_absorption) and
      #RADIATION_ISRF_LW_POPULATION
      (#radiation_packet_constants.LW_population_track).
      Any bit set turns on the dissociating energy bookkeeping. Read through
      #feedback_part_has_H2_lines and #feedback_part_tracks_LW_dissociating. */
  char ISRF_H2_line_absorption;

  /*! H2 column in cm^-2 of the own LW light, weighted by its dissociating
      energy, from the end of the gradient loop. 0 without dissociating light
      or without H2 line absorption. */
  float LW_line_H2_column;
#endif

  /* Save quantities computed in the #hydro density loop for feedback loop */
  struct {
    /*! Neighbour number count. */
    float wcount;
  } density;
};

/* The shared accessors read the fields above. */
#include "../GEAR/radiation_part.h"

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! Mass received from supernovae and stellar winds */
  float delta_mass;

  /*! Physical thermal energy received from supernovae and stellar winds */
  float delta_E_th;

  /*! Sum over the feedback events of the squared norm of the physical
      momentum received from the blastwave */
  float delta_p_norm_2_sum;

  /*! Sum over the feedback events of the physical kinetic energy that each
      event gave to the particle alone, in the peculiar frame. */
  float delta_E_kin_events;

  /*! Sum over the feedback events of the physical Hubble flow velocity
     (relative to the star) dotted with the feedback momentum of the event. */
  float delta_p_hubble_work;

  /*! Values of the events of this step for the tracers (empty without
      tracers) */
  struct feedback_tracers_pending tracers_pending;

  /*! Comoving feedback momentum received from the blastwave. */
  float delta_p[3];

  /*! Comoving momentum contribution brought solely by the transferred stellar
   * mass */
  float delta_p_ejecta[3];

  /*! Number of supernovae affecting this particle */
  unsigned int number_SN;

  /*! Number of stellar wind events affecting this particle */
  unsigned int number_winds;

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Radiation pressure momentum received this step */
  struct feedback_xpart_radiation_data radiation;
#endif

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! HII ionization payload computed by the owner, local to its rank. The tag
      itself lives in #feedback_part_data. */
  struct feedback_xpart_HII_region_data HII_region;
#endif

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Indicator if the particle receives momentum from radiation pressure */
  char hit_by_radiation;
#endif
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Is the star dead? */
  int is_dead;

  /*! Normalisation factor used for the enrichment. Corresponds to the
     denominator in eq (9) in https://arxiv.org/abs/1707.07010  */
  float enrichment_weight;

  /*! Comoving SPH gas density at the star position */
  float gas_density;

  /*! Does the particle needs the feedback loop? */
  char will_do_feedback;

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Does the particle needs the HII ionization loop? */
  char will_do_HII_ionization;
#endif

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Gas density gradient at the star location */
  float grad_rho_star[3];

  /*! Gas metallicity at the star location, weighted like #gas_density */
  float Z_star;
#endif

  /*! Parameters to be accumulated in the feedback loops. Used to compute the
     vector weights (isotropic distribution) */
  double f_sum_plus_term[3];
  double f_sum_minus_term[3];

  /*! Sums of the |w_j|-weighted comoving gas properties used to compute the
     terminal momentum. Normalize with feedback_get_weighted_gas_density() and
     feedback_get_weighted_gas_metallicity(). */
  float weighted_gas_density;
  double weighted_gas_metallicity;

  /*! Number of Ia supernovae */
  float number_snia;

  /*! Number of II supernovae */
  float number_snii;

  /* Supernovae data struct */
  struct {

    /*! Energy injected in the surrounding particles */
    float energy_ejected;

    /*! Total mass ejected by the supernovae */
    float mass_ejected;

  } supernovae;

  /*! Chemical composition of the mass ejected */
  double metal_mass_ejected[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Stellar winds data struct */
  struct {

    /*! Energy injected in the surrounding particles */
    float energy_ejected;

    /*! Mass injected in the surrounding particles */
    float mass_ejected;

  } winds;

#if FEEDBACK_GEAR_MECHANICAL_MODE == 2
  struct {
    /*! Accumulated value for the total energy available in the SN, taking into
       account gas-star motion. This is eq (A4) (lower formula) sum terms in
       https://arxiv.org/abs/2404.16987, without the 0.5*m_ej. */
    float E_total;

    /*! Parameters to determine the coupled energy, momentum and internal energy
       of the SN */
    float beta_1; /* Accumulated value for beta_1 */
    float beta_2; /* Accumulated value for beta_2 */
  } accumulator_sn;

  struct {
    /*! Accumulated value for the total energy available in the SN, taking into
       account gas-star motion. This is eq (A4) (lower formula) sum terms in
       https://arxiv.org/abs/2404.16987, without the 0.5*m_ej. */
    float E_total;

    /*! Parameters to determine the coupled energy, momentum and internal energy
       of the SN */
    float beta_1; /* Accumulated value for beta_1 */
    float beta_2; /* Accumulated value for beta_2 */
  } accumulator_winds;
#endif

  /*! Checks that the sum of the fluxes is 0. These ensures the weights are
     properly constructed. */
#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  struct {
    double delta_m;
    double delta_p_norm;
    double delta_p[3];
  } fluxes_conservation_check;
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */

  /*! Radiation data structs */
  struct feedback_spart_radiation_data radiation;
};

/**
 * @brief Feedback fields carried by each sink particle.
 */
struct feedback_sink_data {

#if defined(GEAR_SUBGRID_RADIATION_ISRF) && defined(SWIFT_DEBUG_CHECKS_ISRF)
  /*! Packet ledger of the gas the sink formed from and swallowed, and of the
      sinks it merged with. */
  struct feedback_isrf_packet_removed_data isrf_packet;
#endif
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H */
