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
 * @brief Feedback fields carried by each hydro particles
 *
 * The subgrid radiation fields have the names and the meaning of the GEAR
 * thermal module's #feedback_part_data: the shared GEAR radiation code reads
 * them by name. They live in #part, not #xpart, so the HII tag follows MPI
 * exchange and restarts.
 */
struct feedback_part_data {

  /*! Tag to mark the particle as ionized. */
  char is_ionized;

  /*! Largest #part.time_bin among this particle and its neighbours in the
      ISRF density loop of this h-iteration. Drives #c_hyp. Sits in the
      padding after #is_ionized, so #part does not grow. */
  timebin_t max_ngb_time_bin;

  /*! Id of the star that ionized this particle. */
  long long star_id;

  /*! Simulation time until which this particle stays flagged as ionized. */
  double end_time;

  /*! Neutral hydrogen mass fraction cached by the cooling step (grackle_0:
      1.0f). Not read by anything yet: an MPI-consistent consumer must read
      the PREVIOUS pass's value. 0 until first written, not "fully ionized". */
  float neutral_H_frac;

  /*! Per-moment ISRF transport state, indexed by #radiation_isrf_moment. */
  struct feedback_isrf_moment_data isrf_moment[ISRF_MOMENT_COUNT];

  /*! Per-operator ISRF band-physics coefficients, indexed by
      #radiation_isrf_operator; a moment reaches its own operator through
      #radiation_isrf_moment_to_operator. */
  struct feedback_isrf_operator_data isrf_operator[ISRF_OPERATOR_COUNT];

  /*! Comoving density snapshot of the previous step, used by the gradient
      and force loops so `grad` stays the adjoint of `div`. Seeded to 1 at
      first init, never 0. */
  float rho_prev;

  /*! Hyperbolic propagation speed, shared by both bands. Under
      #isrf_c_hyp_scheme_fixed_fraction it is `f*c`. Under
      #isrf_c_hyp_scheme_kernel_local_reduced_flux it is
      `min(C_hyp*h_i/dt_max(i), c)`, with `dt_max(i)` set by #max_ngb_time_bin.
      The CFL condition is not guaranteed for pairs with `H_i <= r < H_j`. */
  float c_hyp;

  /*! This particle's own physical timestep, cached once per step at the
      drift. */
  float dt_prev;

#ifdef SWIFT_DEBUG_CHECKS
  /*! #engine.snapshot_output_count at the last write of
      #feedback_isrf_moment_data.u_min_since_snapshot. */
  int u_min_snapshot_index;
#endif

  /*! Step (#engine.ti_current) at which #feedback_isrf_moment_data.u was
      last written. With ISRF_propagation off, a match with the current step
      means a further star sums into u, a mismatch means u is zeroed first.
      Set to -1 at first init. */
  integertime_t ISRF_last_touch_ti;

  /*! Is the particle illuminated by a star's injection pass, with the
      episode still live? Gates a first-touch-only timestep_sync_part call and
      is cleared once #ISRF_illumination_end_ti lapses, like #is_ionized. */
  char is_illuminated_ISRF;

  /*! Cross-bin pair booking state: < 0 inactive (may receive pending), > 0
      active (equals #dt_prev), 0 legacy (every pair booked as before). */
  float dt_active;

  /*! Integer time until which #is_illuminated_ISRF stays set, renewed on
      every injection touch with a margin of
      RADIATION_ISRF_TAG_LIFETIME_INTERVALS star steps. Set to -1 at first
      init. */
  integertime_t ISRF_illumination_end_ti;

  /*! Integer time by which every dose in
      #feedback_isrf_moment_data.u_dose_reservoir is fully drained, extended on
      every star touch. Set to -1 at first init. */
  integertime_t ISRF_reservoir_end_ti;

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /* Trace the metals received from feedback events. This is similar to not
     diffusing metals */
  double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

  /* Save quantities computed in the #hydro density loop for feedback loop */
  struct {
    /*! Neighbour number count. */
    float wcount;
  } density;
};

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

  /*! Radiation pressure momentum received this step */
  struct feedback_xpart_radiation_data radiation;

  /*! HII ionization payload computed by the owner, local to its rank. The tag
      itself lives in #feedback_part_data. */
  struct feedback_xpart_HII_region_data HII_region;

  /*! Indicator if the particle receives momentum from radiation pressure */
  char hit_by_radiation;
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

  /*! Does the particle needs the HII ionization loop? */
  char will_do_HII_ionization;

  /*! Gas density gradient at the star location */
  float grad_rho_star[3];

  /*! Gas metallicity at the star location, weighted like #gas_density */
  float Z_star;

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

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H */
