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

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "../GEAR/feedback_tracers_struct.h"
#include "chemistry_struct.h"
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
 * @brief The band operators (dust opacity) a moment is evaluated against.
 * Indexes #feedback_part_data.isrf_operator.
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
/*! Equal solid-angle cones (HEALPix nside 1) a star splits its packet
    emission into. */
#define RADIATION_PACKET_EMISSION_CONES 12

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
 * @brief Per-moment ISRF state of a hydro particle.
 */
struct feedback_isrf_moment_data {

  /*! Local specific radiation field of the band (physical, per unit mass).
      Double, as the field it is set from. */
  double u;
};

/**
 * @brief Per-operator ISRF band coefficients of a hydro particle.
 */
struct feedback_isrf_operator_data {

  /*! Local linear dust absorption rate of the band, cached once per step. Not
      the injection-side kappa of #radiation_get_part_ISRF_extinction_factors.
   */
  float kappa;
};

/*! Number of test beams of the packet transport, see
    #radiation_packet_inject_test_beams. */
#define RADIATION_PACKET_TEST_BEAMS 2

/*! Numbers of a test beam: the lower and the upper corner of the box (3 + 3),
    the direction (3) and the energy per particle in the PE and LW band (2). */
#define RADIATION_PACKET_TEST_BEAM_FIELDS 11
#define RADIATION_PACKET_TEST_BEAM_LOWER 0
#define RADIATION_PACKET_TEST_BEAM_UPPER 3
#define RADIATION_PACKET_TEST_BEAM_DIRECTION 6
#define RADIATION_PACKET_TEST_BEAM_ENERGY_PE 9
#define RADIATION_PACKET_TEST_BEAM_ENERGY_LW 10

/*! Reception slots of the packet of a gas particle: one per cell of an
    equal-solid-angle partition of the directions, fixed in the frame
    #feedback_props.ISRF_packet_octant_frame (#radiation_packet_slot). Set by
    `--with-isrf-packet-slots`: 12 (the default),
    8 (the octants of the frame) or 16 (octants halved by the plane
    `|x| = |y|`). Each slot costs `sizeof(struct
    feedback_isrf_packet_slot_data)` bytes in every gas particle. */
#ifndef RADIATION_PACKET_SLOTS
#define RADIATION_PACKET_SLOTS 12
#endif

/*! Near arrival bins per slot, one step wide each, on the receiver's own step
    lattice. */
#ifndef RADIATION_PACKET_BINS
#define RADIATION_PACKET_BINS 4
#endif

/*! Far arrival cells per slot, for light that arrives beyond the near bins.
    Their width is a power of two times the step. */
#ifndef RADIATION_PACKET_FAR_CELLS
#define RADIATION_PACKET_FAR_CELLS 2
#endif

#if RADIATION_PACKET_SLOTS != 8 && RADIATION_PACKET_SLOTS != 12 && \
    RADIATION_PACKET_SLOTS != 16
#error "RADIATION_PACKET_SLOTS must be 8, 12 or 16."
#endif

#if RADIATION_PACKET_FAR_CELLS < 1
#error "RADIATION_PACKET_FAR_CELLS must be at least 1."
#endif

/*! Text of the snapshot descriptions: what a reception slot is. */
#if RADIATION_PACKET_SLOTS == 8
#define RADIATION_PACKET_SLOT_TEXT                                         \
  "8 reception slots at its last merge, away from its origin: slot 0 "     \
  "first, then x, y, z. Slot k is the octant of the fixed frame with bit " \
  "0 set for positive first frame coordinate, bit 1 for the second, bit "  \
  "2 for the third."
#elif RADIATION_PACKET_SLOTS == 12
#define RADIATION_PACKET_SLOT_TEXT                                          \
  "12 reception slots at its last merge, away from its origin: slot 0 "     \
  "first, then x, y, z. Slot k = 4 m + 2 a + b is the face of a rhombic "   \
  "dodecahedron in the fixed frame: m is the frame axis with the smallest " \
  "absolute coordinate (the lowest axis on a tie), a and b are 1 for a "    \
  "coordinate that is zero or positive on the two other axes, in axis "     \
  "order."
#else
#define RADIATION_PACKET_SLOT_TEXT                                          \
  "16 reception slots at its last merge, away from its origin: slot 0 "     \
  "first, then x, y, z. Slot k = o + 8 h: o is the octant of the fixed "    \
  "frame (bit 0 set for positive first frame coordinate, bit 1 for the "    \
  "second, bit 2 for the third) and h is 1 when the absolute second frame " \
  "coordinate exceeds the absolute first."
#endif

/**
 * @brief One band of the own field and, with SWIFT_DEBUG_CHECKS_ISRF, the
 * ledger of a gas particle.
 *
 * Photon numbers are stored as `N E_lo`, the number times the lower edge photon
 * energy of the band: a raw count overflows a float. The mean photon energy is
 * then `E_lo energy / photons`. Ledger terms are double.
 */
struct feedback_isrf_packet_band_data {

  /*! Own field: specific energy read off the held packet, summed over the
      slots, or its running mean over
      #feedback_props.ISRF_packet_read_mean_time. */
  float field_energy;

  /*! Own field of the photon number, as `N E_lo` per unit mass. */
  float field_photons;

  /*! Kernel mean of #field_energy over the neighbours and the particle,
      accumulated in the force loop and normalised at its end. */
  float ngb_mean_field_energy;

  /*! Kernel mean of #field_photons, as #ngb_mean_field_energy. */
  float ngb_mean_field_photons;

#ifdef SWIFT_DEBUG_CHECKS_ISRF
  /*! Energy absorbed by dust in this particle, cumulative. */
  double absorbed_energy;

  /*! Photon number absorbed by dust in this particle, as `N E_lo`. */
  double absorbed_photons;

  /*! Energy of packets with no downstream neighbour, cumulative. */
  double escaped_energy;

  /*! Photon number escaped, as `N E_lo`. */
  double escaped_photons;

  /*! Energy removed from this band by the expansion in the shares booked into
      this particle, cumulative: redshift and band-edge loss, minus the gain
      from the other band. Zero in a static box. */
  double redshifted_energy;

  /*! Photon number removed by the expansion, as `N E_lo`, net of the gain
      from the other band. */
  double redshifted_photons;
#endif
};

/**
 * @brief One band of one reception slot: the transported light.
 */
struct feedback_isrf_packet_slot_band_data {

  /*! Energy held, sent downstream at this particle's next active step. */
  float energy;

  /*! Photon number held, as `N E_lo`. */
  float photons;

  /*! Energy received and not yet released, per near arrival bin (see
      #RADIATION_PACKET_BINS). Without a speed limit everything is in bin 0. */
  float bin_energy[RADIATION_PACKET_BINS];

  /*! Photon number in the near arrival bins, as `N E_lo`. */
  float bin_photons[RADIATION_PACKET_BINS];

  /*! The speed of the light is limited (#radiation_packet_light_speed_limited)
      or not, for the whole run, so the two sets of fields below share their
      bytes. */
  union {
    struct {
      /*! Energy in the far arrival cells (see #RADIATION_PACKET_FAR_CELLS),
          with a speed limit. */
      float far_energy[RADIATION_PACKET_FAR_CELLS];

      /*! Photon number in the far arrival cells, as `N E_lo`. */
      float far_photons[RADIATION_PACKET_FAR_CELLS];
    };
    struct {
      /*! Energy received from coarser senders, released over their step,
          without a speed limit. */
      float reservoir_energy;

      /*! Photon number in the reservoir, as `N E_lo`. */
      float reservoir_photons;
    };
  };
};

/**
 * @brief One reception slot of a gas particle: the light that arrived from one
 * cell of the fixed frame, with one direction and origin for all bands.
 */
struct feedback_isrf_packet_slot_data {

  /*! Per-band contents, indexed by #radiation_isrf_operator. */
  struct feedback_isrf_packet_slot_band_data band[ISRF_OPERATOR_COUNT];

  /*! Unit direction of propagation, away from the origin. */
  float direction[3];

  /*! Comoving distance from the origin to this particle. */
  float origin_distance;

  /*! Sum over the inbox shares of their energy times the comoving offset
      from this particle to their origin. */
  double inbox_origin_offset[3];

  /*! Sum of the transmission weights over the downstream neighbours, from
      the gradient loop. */
  float weight_sum;

  /*! Weighted mean comoving advance along #direction per hop, from the
      gradient loop. */
  float forward_advance;

  /*! The speed of the light is limited (#radiation_packet_light_speed_limited)
      or not, for the whole run, so the two sets of fields below share their
      bytes. */
  union {
    struct {
      /*! Same sum as #inbox_origin_offset, over the reservoir, without a speed
          limit. */
      double reservoir_origin_offset[3];

      /*! Integer time by which the reservoir is fully released, without a
          speed limit. */
      integertime_t reservoir_end_ti;
    };
    struct {
      /*! Latest arrival time of the light of each near arrival bin, as the
          distance below the end of the bin. Rounded so that the time is
          later, never earlier. Valid for a bin that holds energy. */
      float bin_stamp_below_end[RADIATION_PACKET_BINS];

      /*! Distance from the end of each far cell down to the first arrival
          time of its light, when the cell holds the light of one share, which
          arrives evenly between this time and the one of #far_last_below. A
          negative value marks a cell that holds the light of several shares.
          Rounded so that the time is later. */
      float far_first_below[RADIATION_PACKET_FAR_CELLS];

      /*! Distance from the end of each far cell down to the latest arrival
          time of its light, rounded so that the time is later. Negative for
          light that arrives beyond the end of the last cell. Valid for a cell
          that holds energy. */
      float far_last_below[RADIATION_PACKET_FAR_CELLS];

      /*! Latest arrival time of the light in the held packet. A share booked
          from it arrives its path over the limit speed after this time, so
          that no part of the packet moves faster than the limit. */
      double arrival;
    };
  };
};

/**
 * @brief Radiation packet of a gas particle: the reception slots, the own field
 * and the ledger.
 */
struct feedback_isrf_packet_data {

  /*! Per-band field and ledger, indexed by #radiation_isrf_operator. */
  struct feedback_isrf_packet_band_data band[ISRF_OPERATOR_COUNT];

  /*! Reception slots, indexed by #radiation_packet_slot. */
  struct feedback_isrf_packet_slot_data slot[RADIATION_PACKET_SLOTS];

  /*! Sum of the kernel weights of the field mean. Scratch. */
  float ngb_field_weight_sum;

  /*! Copy of #feedback_props.ISRF_packet_cone_power, set with
      #feedback_part_data.ISRF_packet_transport. */
  float cone_power;

  /*! Copy of #feedback_props.ISRF_packet_octant_frame (row-major rotation
      matrix), set with #cone_power. */
  float octant_frame[9];

  /*! The width of the far arrival cells is the length of the step times two
      to this power (#radiation_packet_far_end). Placed before the 8-byte
      members, where it fills padding. */
  unsigned char far_scale;

  /*! End of the current step of the particle, set at its merge. A finer
      receiver releases its shares by then. */
  integertime_t step_end_ti;

  /*! Length of this particle's current step, set at its merge. */
  integertime_t step_ti;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Number of non-zero shares booked into this particle, reset by a test. */
  int shares_received;
#endif
};

#ifdef SWIFT_DEBUG_CHECKS_ISRF
/**
 * @brief Packet ledger of the gas that a star or a sink formed from or
 * swallowed, or that left a non-periodic box, per band
 * (#radiation_isrf_operator).
 *
 * The light in flight is booked as escaped, the gas's own ledgers as inherited.
 */
struct feedback_isrf_packet_removed_data {

  /*! Energy in flight in the gas when it was removed, lost with it. */
  double removed_energy[ISRF_OPERATOR_COUNT];

  /*! Photon number in flight when the gas was removed, as `N E_lo`. */
  double removed_photons[ISRF_OPERATOR_COUNT];

  /*! Absorbed, escaped and redshifted energy of the gas's own ledger, carried
      over. */
  double inherited_energy[ISRF_OPERATOR_COUNT];

  /*! Absorbed, escaped and redshifted photon number of the gas, as `N E_lo`. */
  double inherited_photons[ISRF_OPERATOR_COUNT];
};
#endif /* SWIFT_DEBUG_CHECKS_ISRF */

/**
 * @brief Feedback fields carried by each hydro particle.
 *
 * Lives in #part, not #xpart, so the HII tag follows MPI exchange and
 * restarts.
 */
struct feedback_part_data {

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Tag to mark the particle as ionized. */
  char is_ionized;
#endif

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

  /*! Set when this rank claims a foreign copy of this particle, so the claim
      is reported to the owner. Cleared when the report is packed. */
  char claimed_this_pass;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /*! Per-moment ISRF state, indexed by #radiation_isrf_moment. */
  struct feedback_isrf_moment_data isrf_moment[ISRF_MOMENT_COUNT];

  /*! Per-operator ISRF band-physics coefficients, indexed by
      #radiation_isrf_operator; a moment reaches its own operator through
      #radiation_isrf_moment_to_operator. */
  struct feedback_isrf_operator_data isrf_operator[ISRF_OPERATOR_COUNT];

  /*! Comoving density snapshot of the previous step, read by the packet
      transport. Seeded to 1 at first init, never 0. */
  float rho_prev;

  /*! This particle's own physical timestep, cached once per step at the
      drift. */
  float dt_prev;

  /*! Step (#engine.ti_current) at which #feedback_isrf_moment_data.u was
      last written. With ISRF_propagation off, a match with the current step
      means a further star sums into u, a mismatch means u is zeroed first.
      Set to -1 at first init. */
  integertime_t ISRF_last_touch_ti;

  /*! Is the particle illuminated by a star's injection pass, with the
      episode still live? Gates a first-touch-only timestep_sync_part call and
      is cleared once #ISRF_illumination_end_ti lapses, like #is_ionized. */
  char is_illuminated_ISRF;

  /*! Copy of #feedback_props.ISRF_propagation, set at drift and in the extra
      ghost, so the neighbour loops read it without the engine. */
  char ISRF_packet_transport;

  /*! Integer time until which #is_illuminated_ISRF stays set, renewed on
      every injection touch with a margin of
      RADIATION_ISRF_TAG_LIFETIME_INTERVALS star steps. Set to -1 at first
      init. */
  integertime_t ISRF_illumination_end_ti;

  /*! Radiation packet of the packet transport. */
  struct feedback_isrf_packet_data isrf_packet;
#endif

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /* Trace the metals received from feedback events. This is similar to not
     diffusing metals */
  double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif
};

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)

/**
 * @brief HII claim sent back by the rank that ran the star's search to the
 * rank that owns the gas (task_subtype_part_hii_tag).
 *
 * One entry per particle of the cell. The owner merges only the entries
 * stamped by #claimed_this_pass.
 */
struct hii_tag_report {

  /*! Id of the star that claimed the particle. */
  long long star_id;

  /*! Simulation time until which the particle stays tagged. */
  double end_time;

  /*! Photons charged for this claim, counted as forfeited if it loses. */
  double cost;

  /*! Squared distance to the claiming star at claim time. The owner may hold
      no copy of that star, so it is sent, not recomputed. */
  float r2;

  /*! Ionization tag of the claim. */
  char is_ionized;

  /*! Set if this entry is a claim of this pass; others are not merged. */
  char claimed_this_pass;
};

/**
 * @brief Post-cooling HII state sent by the gas owner to every rank that
 * holds a foreign copy (task_subtype_part_hii_state).
 *
 * The search on the other rank then sees the owner's tag and eligibility
 * caches of this step, not the pre-cooling snapshot of the xv/rho channel.
 * The unpack writes only these fields.
 */
struct hii_state_update {

  /*! Id of the star holding the tag. */
  long long star_id;

  /*! Simulation time until which the particle stays tagged. */
  double end_time;

  /*! Photons debited for the current claim. */
  double cost;

  /*! Eligibility temperature cache, see #feedback_part_data. */
  float T_eligibility;

  /*! Mean molecular weight cache, see #feedback_part_data. */
  float mu_eligibility;

  /*! Neutral hydrogen mass fraction cache, see #feedback_part_data. */
  float neutral_H_frac;

  /*! Squared distance to the star holding the tag, at claim time. */
  float r2;

  /*! Ionization tag. */
  char is_ionized;
};

#else

/*! Empty payloads: without HII under MPI these subtypes carry nothing. */
struct hii_tag_report {};
struct hii_state_update {};

#endif

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! Metal mass received from supernovae */
  double delta_metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Values of the events for the tracers, given at the update */
  struct feedback_tracers_pending tracers_pending;

  /*! Thermal energy (not specific) received from supernovae and winds */
  float delta_E_th;

  /*! Momentum received from a supernova */
  float delta_p[3];

  /*! Part of delta_p that the winds direct away from the stars, the part that
     the multiple-event correction rescales (comoving) */
  float delta_p_directed[3];

  /*! Sum over the events of the squared physical directed momentum */
  float delta_p_norm_2_sum;

  /*! Sum over the events of the physical kinetic energy that each event would
     give to the particle alone */
  float delta_E_kin_events;

  /*! Physical work of the directed momentum against the Hubble flow around
     each star */
  float delta_p_hubble_work;

  /*! Number of supernova events received in this timestep (counted only with
     the multiple-event correction) */
  unsigned int number_SN;

  /*! Number of stellar wind events received in this timestep (counted only
     with the multiple-event correction) */
  unsigned int number_winds;

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Radiation struct */
  struct {

    /*! Momentum received from radiation pressure */
    float delta_p[3];
  } radiation;
#endif

#ifdef GEAR_SUBGRID_RADIATION_HII
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
#endif

  /*! Indicator if the particule receive energy from SN specifically */
  char hit_by_SN;

  /*! Indicator if the particle receives energy from SW specifically */
  char hit_by_winds;

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Indicator if the particle receives energy from radiation pressure
   * specifically */
  char hit_by_radiation;
#endif
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

#ifdef GEAR_SUBGRID_RADIATION_HII
  /*! Does the particle needs to go through the HII ionization loop? */
  char will_do_HII_ionization;
#endif

  /*! Integer number of neighbours */
  int num_ngbs;

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
  /*! Gas density gradient at the star location */
  float grad_rho_star[3];

  /*! Gas metallicity at the star location */
  float Z_star;
#endif

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

    /*! Number of HII passes made by this star, the clock of the cone frame
        redraw. */
    int HII_pass_count;

    /*! Mean photon energy above 13.6 eV, cached once per HII rebuild pass, in
        erg. Only set with GEARFeedback:HII_couple_ionization_rate. */
    float mean_excess_photon_energy_HI;
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
    /*! Moment luminosity (physical units) by #radiation_isrf_moment, from the
        radiation table. 0 unless the interstellar radiation field is on. */
    double L_band[ISRF_MOMENT_COUNT];

    /*! Emission of the ISRF packet transport. */
    struct {
      /*! Sum of the receivers' weights in each emission cone, from the star
          density loop. */
      double cone_weight[RADIATION_PACKET_EMISSION_CONES];

      /*! `E_lo / E_bar` per band, the photon number (as `N E_lo`) emitted
          per unit energy. Indexed by #radiation_isrf_operator. */
      float photon_fraction[ISRF_OPERATOR_COUNT];

#ifdef SWIFT_DEBUG_CHECKS_ISRF
      /*! Energy emitted by the star, cumulative. */
      double emitted_energy[ISRF_OPERATOR_COUNT];

      /*! Photon number (as `N E_lo`) emitted by the star, cumulative. */
      double emitted_photons[ISRF_OPERATOR_COUNT];

      /*! Emitted energy that escapes the cones, cumulative. */
      double escaped_energy[ISRF_OPERATOR_COUNT];

      /*! Emitted photon number (as `N E_lo`) that escapes the cones. */
      double escaped_photons[ISRF_OPERATOR_COUNT];

      /*! Ledger of the gas particle the star was converted from. */
      struct feedback_isrf_packet_removed_data removed;
#endif

#ifdef SWIFT_DEBUG_CHECKS
      /*! Per cone, density loop: receivers, the sum of their IDs, and gas
          inside the kernel beyond the emission edge. */
      int cone_count[RADIATION_PACKET_EMISSION_CONES];
      long long cone_id_sum[RADIATION_PACKET_EMISSION_CONES];
      int cone_edge_count[RADIATION_PACKET_EMISSION_CONES];

      /*! The same over the feedback loop's deliveries, and their weights. */
      int cone_count_delivered[RADIATION_PACKET_EMISSION_CONES];
      long long cone_id_sum_delivered[RADIATION_PACKET_EMISSION_CONES];
      double cone_weight_delivered[RADIATION_PACKET_EMISSION_CONES];

      /*! Smoothing length of the density loop and of the feedback loop. */
      float h_density;
      float h_feedback;

      /*! Was this step's emission booked, so the deliveries can be checked? */
      char cones_booked;
#endif
    } packet;
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

  } radiation;
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

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
