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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_H

#include "chemistry_struct.h"
#include "timeline.h"

/*! Maximum number of HEALPix pixels (12*nside_max^2) of a star's HII budget,
    set by ./configure --with-number-of-hii-angular-pixels. */
#ifndef HII_MAX_ANGULAR_PIXELS
#error "HII_MAX_ANGULAR_PIXELS should be defined by configure (config.h)"
#endif

/**
 * @brief The ISRF moments: PE (6-11.2 eV) and Lyman-Werner (11.2-13.6 eV)
 * specific energy, and the Lyman-Werner photon-number moment. Indexes
 * #feedback_part_data.isrf_moment and #feedback_spart_data.radiation.L_band.
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

#ifdef SWIFT_DEBUG_CHECKS
  /*! Most negative #u written since the previous snapshot, 0 if none. */
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
 * @brief Feedback fields carried by each hydro particle.
 *
 * Lives in #part, not #xpart, so the HII tag follows MPI exchange and
 * restarts.
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

  /*! Integer time until which #is_illuminated_ISRF stays set, renewed on
      every injection touch with a margin of
      RADIATION_ISRF_TAG_LIFETIME_INTERVALS star steps. Set to -1 at first
      init. */
  integertime_t ISRF_illumination_end_ti;

  /*! Integer time by which every dose in
      #feedback_isrf_moment_data.u_dose_reservoir is fully drained, extended on
      every star touch. Set to -1 at first init. */
  integertime_t ISRF_reservoir_end_ti;
};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;

  /*! Momentum received from a supernova */
  float delta_p[3];

  /*! Radiation struct */
  struct {

    /*! Momentum received from radiation pressure */
    float delta_p[3];
  } radiation;

  /*! HII ionization payload computed by the owner, local to its rank. The tag
      itself lives in #feedback_part_data. */
  struct {

    /*! Mean photon energy of the tagging star above 13.6 eV, frozen at tag
        time, in erg. Only set with GEARFeedback:HII_couple_ionization_rate. */
    float excess_photon_energy_HI;

    /*! Photoionization rate coefficient Gamma_HI of the tagging star at this
        particle, frozen at tag time (internal 1/time). Only set with
        GEARFeedback:HII_couple_ionization_rate, 0 otherwise. */
    float photoionization_rate_HI;

  } HII_region;

  /*! Indicator if the particule receive energy from SN specifically */
  char hit_by_SN;

  /*! Indicator if the particle receives energy from SW specifically */
  char hit_by_winds;

  /*! Indicator if the particle receives energy from radiation pressure
   * specifically */
  char hit_by_radiation;
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Is the star dead? */
  int is_dead;

  /*! Gas density at the star location. This is also the inverse of
   * normalisation factor used for the enrichment. */
  float enrichment_weight;

  /*! Does the particle needs to go through the feedback loops? */
  char will_do_feedback;

  /*! Does the particle needs to go through the HII ionization loop? */
  char will_do_HII_ionization;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Gas density gradient at the star location */
  float grad_rho_star[3];

  /*! Gas metallicity at the star location */
  float Z_star;

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

  /*! Radiation data structs */
  struct {

    /*! Bolometric luminosity (physical units) from the stellar evolution */
    double L_bol;

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

    /*! Moment luminosity (physical units) by #radiation_isrf_moment, from the
        radiation table. 0 unless the interstellar radiation field is on. */
    double L_band[ISRF_MOMENT_COUNT];

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

  } radiation;
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
