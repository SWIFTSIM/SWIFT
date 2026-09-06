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

/*! Maximum number of HEALPix angular pixels the HII ionization budget can
    be split across (12*nside_max^2). Every star carries a fixed-size
    array of this length regardless of the run's actual
    GEARFeedback:HII_angular_nside, so this is a memory/generality
    trade-off, not a physics one -- set via
    ./configure --with-number-of-hii-angular-pixels=N (default 12, i.e.
    nside<=1) rather than hardcoded here, so builds that only ever need
    nside<=1 don't pay for a finer split they'll never request.
    radiation_init() (src/feedback/GEAR/radiation.c) errors clearly at
    startup if GEARFeedback:HII_angular_nside implies more pixels than
    this build was configured for. */
#ifndef HII_MAX_ANGULAR_PIXELS
#error "HII_MAX_ANGULAR_PIXELS should be defined by configure (config.h)"
#endif

/**
 * @brief Feedback fields carried by each hydro particles
 *
 * Carries the HII ionization tag core (radiation.c's
 * radiation_tag_part_as_ionized() and friends). A struct part field, not
 * struct xpart, so it rides the particle's normal MPI exchange and restart
 * dump automatically, unlike xpart's owner-only payload
 * (tracers_xpart_data.HII_region: excess_photon_energy_HI,
 * photoionization_rate_HI), which stays local to the owning rank.
 */
struct feedback_part_data {

  /*! Tag to mark the particle as ionized. */
  char is_ionized;

  /*! Id of the star that ionized this particle. */
  long long star_id;

  /*! Simulation time until which this particle stays flagged as ionized. */
  double end_time;

  /*! Neutral hydrogen mass fraction, cached by the cooling step right after
      its species update (grackle_1+: HI_frac; grackle_0: 1.0f
      unconditionally, no species tracked). Not read by anything yet: a
      future MPI scheme's F3 latch must consume the PREVIOUS pass's value,
      not the current step's.

      The cache write is skipped on cooling_new_energy()'s early-return
      paths (subgrid-ionized floor, or pinned by
      IONIZATION_FEEDBACK_DEBUG_FIXED_*_TEMPERATURE_K), so on a particle's
      first such step this field is still its zero-init default, 0.0f, the
      OPPOSITE of the 1.0f "no data" sentinel above. A reader must not treat
      0.0f as "fully ionized" without checking the cache has actually been
      written for that particle. */
  float neutral_H_frac;

  /*! Local specific FUV-band (6-11.2 eV) radiation field, internal
      specific-energy units (per-unit-mass, like this codebase's own
      hydro `u`; NOT cgs, unlike mean_excess_photon_energy_HI above). An
      instantaneous field strength, not an accumulated dose: holds the
      illuminating star(s)' most recently computed contribution, summed
      across every star that touched this particle in the same step
      (#LW_FUV_last_touch_ti), then held unchanged until the next step any
      star touches it again. Never cleared by cooling: a reader gets
      whatever was last written, however long ago that was. */
  float u_FUV;

  /*! Local specific Lyman-Werner-band (11.2-13.6 eV) radiation field,
      internal specific-energy units. See #u_FUV; feeds Grackle's
     RT_H2_dissociation_rate (COOLING_GRACKLE_MODE > 1 only) separately from
     #u_FUV, since the two bands carry different dust opacities. */
  float u_LW;

  /*! Snapshot of #u_FUV/#u_LW taken once per step (feedback_reset_part,
      cell_drift.c), before the density loop's h-iterations begin. The
      propagation update reads and mixes these (not #u_FUV/#u_LW
      directly) so it stays correct no matter how many h-iterations a
      particle or its neighbours need: #u_FUV/#u_LW are its per-iteration
      output, safe to overwrite repeatedly since it is never read back as
      an input mid-step. */
  float u_FUV_prev;
  float u_LW_prev;

  /*! Band-specific local linear dust absorption rate (see
      #radiation_get_part_linear_absorption_rate), cached once per step
      (radiation_snapshot_part_propagation) so the propagation density loop
      does not recompute it, and the same unit conversion, per neighbour
      pair. Propagation-corrected (scaled by
      #feedback_props.LW_FUV_yukawa_lambda_correction): the value used by
      the Yukawa propagation (#radiation_end_density_propagation and
      radiation_propagation_iact.h) exclusively, distinct from and NOT
      interchangeable with the injection-side extinction's own,
      independently-computed, uncorrected kappa
      (#radiation_get_part_LW_FUV_extinction_factors). */
  float kappa_FUV;
  float kappa_LW;

  /*! Simulation step (#engine.ti_current) #u_FUV/#u_LW were last written
      at. radiation_iact_nonsym_feedback_apply compares this against the
      current step: a match means some star already wrote this step, so a
      further touch (a second illuminating star) sums into the existing
      value; a mismatch means this is the first touch this step, so
      #u_FUV/#u_LW are zeroed before summing. This is what makes the field
      an instantaneous strength rather than an ever-growing total, while
      still summing multiple simultaneously-illuminating stars correctly
      within one step. feedback_first_init_part sets this to -1 (never a
      valid step) so the very first touch of a particle's life also
      resets rather than summing onto uninitialized memory. */
  integertime_t LW_FUV_last_touch_ti;

  /*! Has this particle been illuminated (u_FUV or u_LW nonzero) by any
      star's injection pass, and is that illumination episode still live?
      Dedicated flag, not inferred from u_FUV/u_LW themselves, since those
      now reset every step a star touches this particle and so cannot
      signal "newly illuminated" via a zero-crossing. Mirrors #is_ionized's
      claimed/not-claimed cycle, including the reset half: gates a
      first-touch-only timestep_sync_part call in
      radiation_iact_nonsym_feedback_apply, mirroring
      feedback_hii_claim_part/feedback_iact_HII_maintain_ionized_part's own
      claim-vs-maintain split, and is cleared once #LW_FUV_illumination_end_ti
      lapses (radiation_gas.c:radiation_reset_part_LW_FUV_illumination_tag,
      called from feedback_reset_part), exactly as cooling clears #is_ionized
      once its own end_time lapses (cooling_gear_subgrid.h). A particle that
      leaves every illuminating star's kernel, then re-enters one later (a
      star re-approaches, a new star's kernel reaches it, or its own h
      changes), therefore gets a fresh sync on re-illumination instead of
      being silently skipped forever. */
  char is_illuminated_LW_FUV;

  /*! Absolute integer time (#engine.ti_current units) until which
      #is_illuminated_LW_FUV stays set. Renewed to
      `ti_current + RADIATION_LW_FUV_TAG_LIFETIME_INTERVALS * ti_step` on
      every injection touch (ti_step = the illuminating star's own
      integer timestep), whether or not this is the particle's first touch
      this episode -- mirrors #feedback_iact_HII_maintain_ionized_part's
      per-pass renewal of the HII tag's own end_time. The
      RADIATION_LW_FUV_TAG_LIFETIME_INTERVALS buffer (radiation.h) keeps
      the window from lapsing between two touches by a star on a coarser
      time bin than this particle's own once-per-step expiry check
      (feedback_reset_part). feedback_first_init_part sets this to -1 so a
      never-illuminated particle's garbage/zero-initialized state can never
      read as "still illuminated". */
  integertime_t LW_FUV_illumination_end_ti;

  /*! Yukawa propagation's per-step mixing accumulators, one
      denominator/numerator pair per band: sum_w_FUV = sum_j
      W(r_ij,h_i)*exp(-tau_ij), sum_wu_FUV = sum_j
      W(r_ij,h_i)*exp(-tau_ij)*u_j_FUV, accumulated pairwise over gas
      neighbours in the density loop (radiation_propagation_iact.h) and
      consumed once per particle by feedback_end_density(), which divides
      them down to the normalized mixing weight sum and resets them for
      the next step. No restart I/O: recomputed fresh every step, like
      density itself. */
  float isrf_prop_sum_w_FUV;
  float isrf_prop_sum_wu_FUV;
  float isrf_prop_sum_w_LW;
  float isrf_prop_sum_wu_LW;
};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;

  /*! Momemtum received from a supernovae */
  float delta_p[3];

  /*! Radiation struct */
  struct {

    /*! Momemtum received from a radiation_pressure */
    float delta_p[3];
  } radiation;

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

    /*! Number of ionizing photons per unit time, split evenly across the
        n_HII_pixels active angular pixels (this is a HUGE number, so must
        be a double) (physical units) from the stellar evolution. A pure
        rate: never debited, so it stays valid all pass for the rate-coupled
        flux calculation (GEARFeedback:HII_couple_ionization_rate). */
    double dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS];

    /*! Photon *count* spendable per pixel this HII rebuild pass:
        dot_N_ion_pix * dt_back (elapsed time since the last pass), debited
        by radiation_consume_ionizing_photons. Budgeting counts over the real
        elapsed interval, instead of comparing bare rates, is what makes the
        ionized extent independent of the rebuild cadence. */
    double N_ion_budget_pix[HII_MAX_ANGULAR_PIXELS];

    /*! Number of active angular pixels this star is currently using
        (1 = spherical/HEALPix disabled) */
    int n_HII_pixels;

    /*! Mass in the HII region generated by this star particle */
    float mass_HII_region;

    /*! Star age at this HII region's last rebuild pass. Double, since
        dt_back (age now minus this) scales the whole photon budget above. */
    double HII_region_last_rebuild;

    /*! Star age the last time this star was given a chance to rebuild its
        HII region, whether or not gas was found (unlike
        HII_region_last_rebuild, which only advances on an actual rebuild).
        Anchors the photon-budget interval dt_back instead, so a pass
        skipped by a gas-free working-level cell does not make the next
        real pass look like it covers the whole gap -- see
        runner_radiation_feedback.c. */
    double HII_region_last_attempt;

    /*! Mean photon energy above the 13.6 eV HI ionization threshold,
        cached once per HII rebuild pass (only computed when
        GEARFeedback:HII_couple_ionization_rate is on; 0 otherwise). Stored
        in cgs (erg), not internal units, since the absolute per-particle
        value underflows float precision in this project's internal unit
        system. */
    float mean_excess_photon_energy_HI;

    /*! Non-ionizing FUV band luminosity, 6-11.2 eV (physical units), split
        off L_bol via this star's own Teff (radiation_planck_band_fraction).
        Feeds the injection term together with #L_LW; only computed when
        GEARFeedback:with_photoelectric_heating is on, 0 otherwise. */
    double L_FUV;

    /*! Lyman-Werner band luminosity, 11.2-13.6 eV (physical units): H2
        photodissociating photons. See #L_FUV. */
    double L_LW;

  } radiation;
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
