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
 * radiation_tag_part_as_ionized() and friends): a struct part field, not
 * struct xpart, so it is included in the particle's normal MPI exchange and
 * restart dump automatically, unlike xpart's owner-only payload
 * (tracers_xpart_data.HII_region -- excess_photon_energy_HI,
 * photoionization_rate_HI) which stays local to whichever rank currently
 * owns the particle.
 */
struct feedback_part_data {

  /*! Tag to mark the particle as ionized. */
  char is_ionized;

  /*! Id of the star that ionized this particle. */
  long long star_id;

  /*! Simulation time until which this particle stays flagged as ionized. */
  double end_time;

  /*! Neutral hydrogen mass fraction, cached by the cooling step right after
      its species update (grackle_1+: HI_frac; grackle_0, which tracks no
      species: 1.0f unconditionally). Not read by anything yet -- this is
      the input a future MPI scheme's F3 latch will consume, and it MUST
      read the PREVIOUS pass's value rather than the current step's (see
      PLAN_radiation_mpi.md blocker F3: the owner's HII pass and a same-step
      cooling call race on this same-step value today, on one rank as much
      as across ranks).

      The cache write is skipped on cooling_new_energy()'s early-return
      paths (a particle held at the subgrid-ionized floor, or pinned by
      IONIZATION_FEEDBACK_DEBUG_FIXED_*_TEMPERATURE_K): on those steps this
      field keeps its previous value, which on a particle's very first such
      step is still the struct's zero-init default, 0.0f -- the OPPOSITE of
      the 1.0f "no data" sentinel above. A future reader must not treat 0.0f
      as "fully ionized" without first checking whether the cache has ever
      actually been written for that particle. */
  float neutral_H_frac;
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

  } radiation;
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
